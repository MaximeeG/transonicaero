# this ideally needs to be run in an virtual environment. Reminder how to activate it:
# when in HW1 folder: source .venv/bin/activate
#
# if .venv for some reason doesnt exist:
# python3 -m venv .venv
# then install relevant libs. Dont ask copilot for this, it will fuck it up

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from pathlib import Path
import re

# font update
plt.rcParams.update({
    "font.family": "sans-serif",
    "mathtext.fontset": "dejavusans",
})

# --------------------------------------------------
# Load data
# --------------------------------------------------

data_dir = Path("HW1/plots/data")
figures_dir = Path("HW1/plots/figures")

csv_files = sorted(data_dir.glob("*.csv"))

if not csv_files:
    raise FileNotFoundError(f"No CSV files found in {data_dir}")

for filename_path in csv_files:

    filename = str(filename_path)

    # get algorithm name and CFL number from filename
    file_stem = Path(filename).stem

    match = re.fullmatch(
        r"WAVE_(.+)_CFL_(\d+)_NX_(\d+)",
        file_stem
    )

    if not match:
        print(f"Skipping file with unexpected name format: {filename}")
        continue

    algname = match.group(1)
    cfl = int(match.group(2)) / 100.0
    nx = int(match.group(3))

    wave_name = file_stem

    output_3d = Path(f"HW1/plots/figures/{wave_name}_3D.png")
    output_comparison = Path(f"HW1/plots/figures/{wave_name}_COMPARISON.png")

    if output_3d.exists() and output_comparison.exists():
        print(f"Skipping {wave_name}: plots already exist")
        continue

    times = []
    c_values = []
    u_vectors = []

    with open(filename, "r") as f:

        next(f)  # skip header: t,c,u

        for line in f:

            t_str, c_str, u_str = line.strip().split(",", 2)

            t = float(t_str)

            c_value = float(c_str)

            u = np.fromstring(u_str.strip().strip("[]"), sep=",")

            times.append(t)

            c_values.append(c_value)

            u_vectors.append(u)

    times = np.array(times)

    c_values = np.array(c_values)

    U = np.array(u_vectors)   # shape: (num_times, num_x_points)


    # check if numerical solution has blown up too much to plot
    max_u = np.nanmax(np.abs(U))

    if not np.all(np.isfinite(U)) or max_u > 1e100:
        print(
            f"Skipping {wave_name}: numerical solution blew up "
            f"(max |u| = {max_u:.3e})"
        )

        continue

    # c value taken directly from the csv file
    c = c_values[0]

    # c value taken directly from the csv file
    c = c_values[0]

    # make sure c is constant throughout the csv
    if not np.allclose(c_values, c):

        raise ValueError(f"c is not constant throughout the CSV file: {filename}")

    # x values over [0, pi]
    x = np.linspace(0, np.pi, U.shape[1])


    # --------------------------------------------------
    # Analytical solution
    # --------------------------------------------------

    def u_exact(x, t):

        x0 = x - c * t

        return np.where(

            (x0 >= 0.5) & (x0 <= 1.0),

            1.0,

            0.0

        )


    # --------------------------------------------------
    # Color map based on time
    # --------------------------------------------------

    norm = Normalize(vmin=times.min(), vmax=times.max())

    # cmap = cm.rainbow   # or cm.viridis, cm.plasma, cm.turbo
    # cmap = cm.plasma
    cmap = cm.viridis
    # cmap = cm.turbo


    # --------------------------------------------------
    # 3D plot
    # --------------------------------------------------

    fig = plt.figure(figsize=(12, 8))

    ax = fig.add_subplot(111, projection="3d")

    # maximum number of time slices shown in the 3D plot
    max_lines = 60

    plot_indices = np.linspace(
        0,
        len(times) - 1,
        min(max_lines, len(times)),
        dtype=int
    )

    for i in plot_indices:

        t = times[i]

        ax.plot(
            x,
            np.full_like(x, t),
            U[i],
            color=cmap(norm(t)),
            linewidth=1.2,
            alpha=0.85
        )

    ax.set_xlabel("x")

    ax.set_ylabel("Time t")

    ax.set_zlabel("u(x,t)")

    ax.set_title(

        f"Evolution of u(x,t) with Color Showing Time\n"

        f"Algorithm: {algname}, CFL = {cfl:g}, $N_x$ = {nx}"

    )

    # Add colorbar to show time-to-color mapping
    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)

    mappable.set_array([])

    cbar = fig.colorbar(mappable, ax=ax, pad=0.1, shrink=0.7)

    cbar.set_label("Time t")

    plt.tight_layout()

    plt.savefig(

        output_3d,

        dpi=300,

        bbox_inches="tight"

    )

    plt.close()

    print(f"Saved as {wave_name}_3D.png")


    # --------------------------------------------------
    # Analytical vs numerical comparison
    # --------------------------------------------------

    # use a denser x grid so the analytical discontinuities are sharp
    x_exact = np.linspace(0, np.pi, 2000)

    # middle and final numerical time steps
    middle_index = len(times) // 2

    final_index = len(times) - 1

    middle_time = times[middle_index]

    final_time = times[final_index]

    # numerical solution at t = 0 is the initial condition
    U_initial = np.where(

        (x >= 0.5) & (x <= 1.0),

        1.0,

        0.0

    )

    fig, axes = plt.subplots(

        3,

        1,

        figsize=(10, 10),

        sharex=True,

        sharey=True

    )


    # --------------------------------------------------
    # Top plot: t = 0
    # --------------------------------------------------

    axes[0].plot(

        x,

        U_initial,

        linewidth=2,

        label="Numerical"

    )

    axes[0].plot(

        x_exact,

        u_exact(x_exact, 0.0),

        color="black",

        linestyle=":",

        linewidth=2.5,

        label="Analytical"

    )

    axes[0].set_title("t = 0")

    axes[0].set_ylabel("u(x,t)")

    axes[0].grid(True, alpha=0.3)

    axes[0].legend()


    # --------------------------------------------------
    # Middle plot
    # --------------------------------------------------

    axes[1].plot(

        x,

        U[middle_index],

        linewidth=2,

        label="Numerical"

    )

    axes[1].plot(

        x_exact,

        u_exact(x_exact, middle_time),

        color="black",

        linestyle=":",

        linewidth=2.5,

        label="Analytical"

    )

    axes[1].set_title(f"t = {middle_time:.3f}")

    axes[1].set_ylabel("u(x,t)")

    axes[1].grid(True, alpha=0.3)

    axes[1].legend()


    # --------------------------------------------------
    # Bottom plot: final time step
    # --------------------------------------------------

    axes[2].plot(

        x,

        U[final_index],

        linewidth=2,

        label="Numerical"

    )

    axes[2].plot(

        x_exact,

        u_exact(x_exact, final_time),

        color="black",

        linestyle=":",

        linewidth=2.5,

        label="Analytical"

    )

    axes[2].set_title(f"t = {final_time:.3f}")

    axes[2].set_xlabel("x")

    axes[2].set_ylabel("u(x,t)")

    axes[2].grid(True, alpha=0.3)

    axes[2].legend()


    fig.suptitle(

        f"Numerical vs Analytical Solution\n"

        f"Algorithm: {algname}, CFL = {cfl:g}, $N_x$ = {nx}",

        fontsize=14

    )

    plt.tight_layout()

    plt.savefig(

        output_comparison,

        dpi=300,

        bbox_inches="tight"

    )

    plt.close()

    print(f"Saved as {wave_name}_COMPARISON.png")