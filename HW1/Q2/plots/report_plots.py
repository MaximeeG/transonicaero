"""Create P01-P20, C01-C20 and the eight paired comparison groups.

Run with the same Python environment as plotting.py. No C files are changed.
The solver must save one complete profile per iteration, starting at iteration 0.
"""
import csv
from itertools import islice
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PLOTS_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = PLOTS_DIR / "figures" / "report"
DISPLAY_FLOOR = 1e-16
# Used only to assess the reconstructed stopping residual, not the density norm.
STOPPING_TOLERANCE = 1e-8

# Exact run order from the report table: algorithm, outlet, CFL*100, Nx.
RUNS = [
    ("MACCORMACK", "SUPER", 10, 201),
    ("MACCORMACK", "SUPER", 50, 201),
    ("MACCORMACK", "SUPER", 110, 201),
    ("MACCORMACK", "SUB", 10, 201),
    ("MACCORMACK", "SUB", 50, 201),
    ("MACCORMACK", "SUB", 110, 201),
    ("BEAM_WARMING", "SUPER", 10, 201),
    ("BEAM_WARMING", "SUPER", 50, 201),
    ("BEAM_WARMING", "SUPER", 110, 201),
    ("BEAM_WARMING", "SUB", 10, 201),
    ("BEAM_WARMING", "SUB", 50, 201),
    ("BEAM_WARMING", "SUB", 110, 201),
    ("MACCORMACK", "SUPER", 50, 101),
    ("MACCORMACK", "SUPER", 50, 401),
    ("MACCORMACK", "SUB", 50, 101),
    ("MACCORMACK", "SUB", 50, 401),
    ("BEAM_WARMING", "SUPER", 50, 101),
    ("BEAM_WARMING", "SUPER", 50, 401),
    ("BEAM_WARMING", "SUB", 50, 101),
    ("BEAM_WARMING", "SUB", 50, 401),
]
GROUPS = [
    ("CFL", [1, 2, 3]), ("CFL", [4, 5, 6]),
    ("CFL", [7, 8, 9]), ("CFL", [10, 11, 12]),
    ("grid", [13, 2, 14]), ("grid", [15, 5, 16]),
    ("grid", [17, 8, 18]), ("grid", [19, 11, 20]),
]
VARIABLES = [("rho", r"Density $\rho$"), ("p", r"Pressure $p$"),
             ("u", r"Velocity $u$"), ("M", r"Mach number $M$")]


def run_name(settings):
    algorithm, outlet, cfl, nx = settings
    return f"{algorithm}_{outlet}_CFL_{cfl}_NX_{nx}"


def label(settings):
    algorithm, outlet, cfl, nx = settings
    algorithm = {"MACCORMACK": "MacCormack", "BEAM_WARMING": "Beam-Warming"}[algorithm]
    outlet = {"SUPER": "supersonic", "SUB": "subsonic"}[outlet]
    return f"{algorithm} | {outlet} outlet | CFL = {cfl / 100:g} | Nx = {nx}"


def read_run(path, nx):
    # Keep two snapshots and the scalar history: memory does not grow with CSV size.
    result = dict(history=[], final=None, failure=None, stopping=None,
                  time=None, iteration=0, status="insufficient data")
    previous = None
    with path.open() as file:
        header = next(csv.reader([file.readline()]))
        required = ["t", "x", "A", "rho", "u", "p", "M", "e", "mass_flow"]
        if not set(required).issubset(header):
            raise ValueError(f"Required columns: {', '.join(required)}")
        col = {name: header.index(name) for name in required}
        snapshot = 0
        while True:
            lines = list(islice(file, nx))
            if not lines:
                break
            if len(lines) != nx:
                raise ValueError(f"Incomplete snapshot {snapshot}: {len(lines)} of {nx} rows")
            values = np.loadtxt(lines, delimiter=",", ndmin=2)
            if values.shape != (nx, len(header)):
                raise ValueError(f"Incorrect row width in snapshot {snapshot}")
            state = {name: values[:, index] for name, index in col.items()}
            if not np.isfinite(values).all() or (state['rho'] <= 0).any() or (state['p'] <= 0).any():
                result.update(failure=snapshot, final=None,
                              status=f"invalid flow at iteration {snapshot}")
                break
            if not (np.diff(state['x']) > 0).all() or not (state['t'] == state['t'][0]).all():
                raise ValueError(f"Snapshot {snapshot} does not match Nx in filename")
            if previous is not None:
                if not np.array_equal(previous['x'], state['x']):
                    raise ValueError("Spatial grid changed between snapshots")
                result['history'].append(float(np.linalg.norm(state['rho'] - previous['rho'])))
                dt = state['t'][0] - previous['t'][0]
                # Reconstruct the C stopping measure, distinct from the plotted norm.
                result['stopping'] = None
                rounding_bound = 0.0
                if dt > 0:
                    old = np.column_stack((previous['rho'] * previous['A'], previous['mass_flow'], previous['e'] * previous['A']))
                    new = np.column_stack((state['rho'] * state['A'], state['mass_flow'], state['e'] * state['A']))
                    scales = np.abs(new[0])
                    if (scales > 0).all():
                        result['stopping'] = float(np.max(np.abs(new[1:-1] - old[1:-1]) / (dt * scales)))
                        # The writer saves 12 significant digits. Reconstructing
                        # products and subtracting snapshots loses precision near
                        # the stopping tolerance; do not claim failure from that.
                        rounding_bound = float(np.max(
                            1e-11 * (np.abs(new[1:-1]) + np.abs(old[1:-1])) / (dt * scales)))
                if result['stopping'] is None:
                    result['status'] = "convergence unverified (saved time resolution)"
                elif result['stopping'] + rounding_bound <= STOPPING_TOLERANCE:
                    result['status'] = "meets reconstructed stopping criterion"
                elif result['stopping'] - rounding_bound <= STOPPING_TOLERANCE:
                    result['status'] = "near stopping tolerance (CSV precision limited)"
                else:
                    result['status'] = "not converged at final saved iteration"
            previous = state
            result.update(final=state, time=float(state['t'][0]), iteration=snapshot)
            snapshot += 1
    if previous is None and result['failure'] is None:
        raise ValueError("No snapshots found")
    return result


def save(fig, path, footer):
    fig.text(0.5, 0.012, footer, ha="center", va="bottom", fontsize=8)
    footer_height = ((footer.count('\n') + 1) * 12 + 12) / 72
    bottom = 0.012 + footer_height / fig.get_figheight()
    fig.tight_layout(rect=(0, bottom, 1, 0.98))
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def profile_axes():
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
    for ax, (_, ylabel) in zip(axes, VARIABLES):
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.25)
    axes[-1].axhline(1, color="black", linestyle="--", linewidth=1, label="M = 1")
    axes[-1].set_xlabel("x")
    return fig, axes


def density_curve(ax, result, text, color):
    history = np.asarray(result['history'])
    ax.plot(np.arange(1, len(history) + 1), np.log10(np.maximum(history, DISPLAY_FLOOR)),
            color=color, linewidth=1, label=text)
    if result['failure'] is not None:
        ax.axvline(result['failure'], color=color, linestyle=":", linewidth=1,
                   label=f"Invalid flow: iteration {result['failure']}")
    ax.set_xlabel("Iterations")
    ax.set_ylabel(r"$\log_{10}(\|\rho^{n+1}-\rho^n\|_2)$")
    ax.grid(True, alpha=0.25)


def individual(number, settings, result):
    stem = run_name(settings)
    fig, axes = profile_axes()
    if result['final'] is not None:
        for ax, (key, _) in zip(axes, VARIABLES):
            ax.plot(result['final']['x'], result['final'][key], linewidth=1.4)
        subtitle = f"Final saved profile: iteration {result['iteration']}, t = {result['time']:.6g}"
    else:
        subtitle = "Invalid final profile omitted"
        for ax in axes:
            ax.text(.5, .5, result['status'], transform=ax.transAxes, ha="center")
    fig.suptitle(f"P{number:02d} — {label(settings)}\n{subtitle}", fontsize=11)
    axes[-1].legend()
    save(fig, OUTPUT_DIR / f"P{number:02d}_{stem}.png", result['status'])

    fig, ax = plt.subplots(figsize=(10, 6))
    density_curve(ax, result, label(settings), "tab:blue")
    ax.set_xlim(left=0)
    ax.legend(fontsize=8)
    fig.suptitle(f"C{number:02d} — Density convergence", fontsize=13)
    save(fig, OUTPUT_DIR / f"C{number:02d}_{stem}.png",
         f"{result['status']}\nUnnormalized density change; display floor {DISPLAY_FLOOR:g}; not the C stopping residual.")


def comparisons(results):
    colors = ["tab:blue", "tab:orange", "tab:green"]
    for group, (study, numbers) in enumerate(GROUPS, 1):
        first = RUNS[numbers[0] - 1]
        title = f"{first[0].replace('_', '-')} | {first[1]} outlet | {study} comparison"
        title += " | Nx = 201" if study == "CFL" else " | CFL = 0.5"
        fig_p, axes = profile_axes()
        fig_c, ax_c = plt.subplots(figsize=(10, 6))
        notes = []
        for number, color in zip(numbers, colors):
            settings, result = RUNS[number - 1], results[number]
            text = f"CFL = {settings[2]/100:g}" if study == "CFL" else f"Nx = {settings[3]}"
            notes.append(f"{text}: {result['status']}")
            density_curve(ax_c, result, text, color)
            if result['final'] is not None:
                for ax, (key, _) in zip(axes, VARIABLES):
                    ax.plot(result['final']['x'], result['final'][key], color=color, linewidth=1.2, label=text)
        for ax in axes:
            if ax.get_legend_handles_labels()[0]:
                ax.legend(fontsize=8)
        ax_c.legend(fontsize=8)
        # Set limits only after adding every curve, so shorter runs do not
        # disable autoscaling and cut off the longer histories.
        ax_c.set_xlim(left=0)
        fig_p.suptitle(f"G{group:02d} — Final saved profiles\n{title}", fontsize=11)
        fig_c.suptitle(f"G{group:02d} — Density convergence\n{title}", fontsize=11)
        save(fig_p, OUTPUT_DIR / f"G{group:02d}_{study}_profiles.png", "\n".join(notes))
        save(fig_c, OUTPUT_DIR / f"G{group:02d}_{study}_convergence.png",
             "\n".join(notes) + "\nUnnormalized density L2 norm (grid-size dependent); display floor 1e-16.")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {}
    for number, settings in enumerate(RUNS, 1):
        path = PLOTS_DIR / 'data' / f"{run_name(settings)}.csv"
        print(f"[{number}/20] Reading {path.name}", flush=True)
        try:
            result = read_run(path, settings[3])
        except (OSError, ValueError, IndexError) as error:
            result = dict(history=[], final=None, failure=None, stopping=None,
                          time=None, iteration=0, status=f"Could not process: {error}")
        results[number] = result
        individual(number, settings, result)
        print(f"  {result['status']}; {len(result['history'])} density differences", flush=True)
    comparisons(results)
    # Index gives the report-table numbering and flags unusable/nonconverged runs.
    lines = ["# Q2 report figures", "", "Each CSV is read once. Re-running regenerates these figures.", "",
             "Profiles are final saved states. Invalid profiles are replaced by explanatory panels.",
             "Convergence is log10 of the unnormalized L2 density change between consecutive saved iterations.",
             "Stopping status is reconstructed from rounded CSV values with tolerance 1e-8; it is not a solver-certified status.",
             "No gamma or dissipation settings are inferred from filenames.", "",
             "| Run | Settings | Outcome | Figures |", "|---|---|---|---|"]
    for number, settings in enumerate(RUNS, 1):
        stem = run_name(settings)
        lines.append(f"| {number} | {label(settings).replace('|', ';')} | {results[number]['status']} | "
                     f"[P{number:02d}](P{number:02d}_{stem}.png), [C{number:02d}](C{number:02d}_{stem}.png) |")
    lines += ["", "## Comparisons", ""]
    for group, (study, numbers) in enumerate(GROUPS, 1):
        lines.append(f"- G{group:02d}, {study}, runs {numbers}: [profiles](G{group:02d}_{study}_profiles.png), "
                     f"[convergence](G{group:02d}_{study}_convergence.png)")
    (OUTPUT_DIR / 'README.md').write_text('\n'.join(lines) + '\n')
    print(f"Saved 40 individual figure slots and 16 comparison figures in {OUTPUT_DIR}", flush=True)


if __name__ == '__main__':
    main()
