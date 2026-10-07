import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt


# data and figures folders are both inside this directory
PLOTS_DIR = Path(__file__).resolve().parent


def main():
    data_dir = PLOTS_DIR / "data"
    figures_dir = PLOTS_DIR / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    csv_files = sorted(data_dir.glob("*.csv"))

    if not csv_files:
        print(f"No CSV files found in {data_dir}")
        return

    for csv_path in csv_files:
        output_path = figures_dir / f"{csv_path.stem}_rho_p_u_M.png"
        # skip before reading the CSV, since these files can be quite large
        if output_path.exists():
            print(f"Skipping {csv_path.name}: plot already exists")
            continue
        plot_csv(csv_path, output_path)


def plot_csv(csv_path, output_path):
    # example: BEAM_WARMING_SUB_CFL_50_NX_201.csv
    # the algorithm can contain underscores; CFL is stored multiplied by 100
    match = re.fullmatch(r"(.+)_(SUPER|SUB)_CFL_(\d+)_NX_(\d+)", csv_path.stem)
    if not match:
        print(f"Skipping {csv_path.name}: expected ALGORITHM_SUPER/SUB_CFL_50_NX_201.csv")
        return

    algorithm_name, outlet_name, cfl_text, nx_text = match.groups()
    algorithm = {"MACCORMACK": "MacCormack", "BEAM_WARMING": "Beam-Warming"}.get(
        algorithm_name, algorithm_name.replace("_", " ")
    )
    cfl = int(cfl_text) / 100.0
    nx = int(nx_text)
    outlet = {"SUPER": "Supersonic", "SUB": "Subsonic"}[outlet_name]

    # the solver writes a full spatial profile at every time step.
    # Keep only the last profile, so a large CSV does not fill up memory.
    last_time = None
    last_x = None
    rows = []
    with csv_path.open(newline="") as file:
        reader = csv.DictReader(file)
        required_columns = {"t", "x", "rho", "p", "u", "M"}
        if not required_columns.issubset(reader.fieldnames or []):
            raise ValueError("CSV must contain the columns t, x, rho, p, u and M")

        for row in reader:
            time = float(row["t"])
            current_x = float(row["x"])
            # time may round to the same CSV value in consecutive snapshots
            if time != last_time or (last_x is not None and current_x <= last_x):
                rows = []
                last_time = time
            rows.append(row)
            last_x = current_x

    if not rows:
        raise ValueError(f"No flow data found in {csv_path}")

    # sort by x so the lines connect neighboring grid points
    rows.sort(key=lambda row: float(row["x"]))
    x = [float(row["x"]) for row in rows]

    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
    variables = [
        ("rho", r"Density $\rho$", "tab:blue"),
        ("p", r"Pressure $p$", "tab:orange"),
        ("u", r"Velocity $u$", "tab:green"),
        ("M", r"Mach number $M$", "tab:purple"),
    ]

    for ax, (name, label, color) in zip(axes, variables):
        values = [float(row[name]) for row in rows]
        ax.plot(x, values, color=color, linewidth=2)
        ax.set_ylabel(label)
        ax.grid(True, alpha=0.3)
        if name == "M":
            ax.axhline(1.0, color="black", linestyle="--", linewidth=1, label="Sonic: M = 1")
            ax.legend()

    axes[-1].set_xlabel("x")
    # run settings are read from the CSV filename
    fig.suptitle(
        f"{csv_path.stem} — final saved time: t = {last_time:.6g}\n"
        f"Algorithm: {algorithm} | Nx: {nx} | CFL: {cfl:g}\n"
        f"Outlet: {outlet}",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Plotted {len(rows)} grid points at t = {last_time:.6g}")
    print(f"Saved as {output_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
