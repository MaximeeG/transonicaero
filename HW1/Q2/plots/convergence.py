import csv
import math
import re
from itertools import islice
from pathlib import Path

import matplotlib.pyplot as plt


# data and figures folders are both inside this directory
PLOTS_DIR = Path(__file__).resolve().parent
# only affects display: log(0) is undefined
DISPLAY_FLOOR = 1e-16


def density_history(csv_path, nx):
    # one full spatial profile is written per iteration, including t = 0.
    # Read only one profile at a time instead of loading the entire CSV.
    residuals = []
    previous_rho = None
    grid = None
    invalid_iteration = None

    with csv_path.open(newline="") as file:
        reader = csv.reader(file)
        header = next(reader, [])
        if not {"t", "x", "rho"}.issubset(header):
            raise ValueError("CSV must contain t, x and rho columns")
        t_col, x_col, rho_col = (header.index(name) for name in ("t", "x", "rho"))
        snapshot = 0

        while True:
            rows = list(islice(reader, nx))
            if not rows:
                break
            if len(rows) != nx:
                raise ValueError(f"Incomplete snapshot {snapshot}: expected {nx} rows, got {len(rows)}")

            x = [float(row[x_col]) for row in rows]
            times = [float(row[t_col]) for row in rows]
            rho = [float(row[rho_col]) for row in rows]
            if not all(math.isfinite(v) for v in x):
                raise ValueError("Nonfinite grid coordinate")
            if any(b <= a for a, b in zip(x, x[1:])):
                raise ValueError("Grid points must increase within each snapshot; check Nx in filename")
            if grid is None:
                grid = x
            elif x != grid:
                raise ValueError("Grid changes between snapshots; check Nx in filename")

            # Invalid data must not appear as successful convergence.
            if not all(math.isfinite(v) for v in times + rho) or any(v <= 0 for v in rho):
                invalid_iteration = snapshot
                break
            if any(t != times[0] for t in times):
                raise ValueError("Mixed times within a snapshot; check Nx in filename")

            if previous_rho is not None:
                # lecture definition: L2 norm of the density change.
                # No division by dt, nx or the initial residual.
                residual = math.hypot(*(a - b for a, b in zip(rho, previous_rho)))
                if not math.isfinite(residual):
                    invalid_iteration = snapshot
                    break
                residuals.append(residual)
            previous_rho = rho
            snapshot += 1

    return residuals, invalid_iteration


def plot_convergence(csv_path, output_path):
    match = re.fullmatch(r"(.+)_(SUPER|SUB)_CFL_(\d+)_NX_(\d+)", csv_path.stem)
    if not match:
        raise ValueError("Expected ALGORITHM_SUPER/SUB_CFL_50_NX_201.csv")
    algorithm, outlet, cfl, nx = match.groups()
    algorithm = {"MACCORMACK": "MacCormack", "BEAM_WARMING": "Beam-Warming"}.get(
        algorithm, algorithm.replace("_", " ")
    )
    outlet = {"SUPER": "Supersonic", "SUB": "Subsonic"}[outlet]
    cfl, nx = int(cfl) / 100.0, int(nx)
    if nx < 2:
        raise ValueError("Nx must be at least 2")

    residuals, invalid_iteration = density_history(csv_path, nx)
    if not residuals and invalid_iteration is None:
        raise ValueError("At least two snapshots are needed to calculate convergence")

    fig, ax = plt.subplots(figsize=(9, 6))
    iterations = range(1, len(residuals) + 1)
    log_residual = [math.log10(max(value, DISPLAY_FLOOR)) for value in residuals]
    ax.plot(iterations, log_residual, linewidth=1.2, label=f"{algorithm}, CFL = {cfl:g}")
    ax.set_xlabel("Iterations")
    ax.set_ylabel(r"$\log_{10}(\|\rho^{n+1}-\rho^n\|_2)$")
    ax.set_title(f"Density convergence — {outlet.lower()} outlet\n$N_x$ = {nx}")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)

    if invalid_iteration is not None:
        ax.axvline(invalid_iteration, color="tab:red", linestyle="--",
                   label=f"Invalid density/time at iteration {invalid_iteration}")
        print(f"Warning: {csv_path.name}: stopped at invalid snapshot {invalid_iteration}")
    if any(value < DISPLAY_FLOOR for value in residuals):
        ax.text(0.02, 0.04, f"Values below {DISPLAY_FLOOR:g} shown at the display floor",
                transform=ax.transAxes, fontsize=9)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output_path} ({len(residuals)} iteration differences)")


def main():
    figures_dir = PLOTS_DIR / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    csv_files = sorted((PLOTS_DIR / "data").glob("*.csv"))
    if not csv_files:
        print(f"No CSV files found in {PLOTS_DIR / 'data'}")
    for csv_path in csv_files:
        output_path = figures_dir / f"{csv_path.stem}_convergence.png"
        # regenerate if the solver has written a newer CSV
        if output_path.exists() and output_path.stat().st_mtime_ns >= csv_path.stat().st_mtime_ns:
            print(f"Skipping {csv_path.name}: convergence plot is up to date")
            continue
        try:
            plot_convergence(csv_path, output_path)
        except (ValueError, IndexError) as error:
            print(f"Skipping {csv_path.name}: {error}")


if __name__ == "__main__":
    main()
