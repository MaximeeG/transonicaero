# Transonic Aerodynamics — HW1

C solvers and Python plotting scripts for MEC6602 Homework 1.

- **Q1:** 1D advection of a square pulse, comparing numerical schemes with the analytical solution.
- **Q2:** Quasi-1D Euler flow through a nozzle, using MacCormack and Beam–Warming with supersonic or subsonic outlet conditions.

## Directory structure

```text
transonicaero/
├── lecture_slides/              # reference material for the Euler equations
└── HW1/
    ├── MEC6602E_HW1.pdf         # assignment statement
    ├── .venv/                  # local Python environment (if created)
    ├── Q1/
    │   ├── src/                # main.c, wavesolver.c, wavesolver.h
    │   ├── build/              # compiled executable: output
    │   ├── makefile            # compiles the wave solver
    │   ├── run.py              # one scheme and parameter combination
    │   ├── runall.py           # batch of schemes and CFL values
    │   └── plots/
    │       ├── plotting.py     # wave evolution and analytical comparisons
    │       ├── data/           # CSV solution histories
    │       └── figures/        # PNG plots
    └── Q2/
        ├── src/                # main.c, eulersolver.c, eulersolver.h
        ├── build/              # compiled executable: output
        ├── makefile            # compiles AND runs the Euler solver
        └── plots/
            ├── plotting.py    # final density, pressure, velocity, Mach
            ├── convergence.py # density convergence histories
            ├── report_plots.py # grouped figures for the report
            ├── data/          # CSV solution histories
            └── figures/
                └── report/    # report figures and their README index
```

In each `src/` directory, `main.c` sets up the run and writes the results; the solver `.c` file contains initialization, numerical methods, boundary treatment, and cleanup. The `.h` file defines the structs and function interfaces. For Q1, the file `run.py`can be used do run the C file, while also having direct access to simulation settings.

## Setup

Requires GCC, Make, and Python 3 with NumPy and Matplotlib. From the repository root:

```bash
python3 -m venv HW1/.venv
source HW1/.venv/bin/activate
python -m pip install numpy matplotlib
```

In later sessions, just activate the same environment. The commands below also start from the repository root.

## Q1 — wave equation

Set `algorithm`, `cfl`, and `nx` at the top of `HW1/Q1/run.py`, then run:

```bash
python HW1/Q1/run.py
```

This compiles the solver, runs the simulation, and plots the saved data. Available schemes are `BACKWARD`, `FORWARD`, `LAX`, `LAX_WENDROFF`, `LEAPFROG`, `THETA`, `OWN2SPACE4TIME`, and `OWN4SPACE2TIME`.

For `THETA`, the theta value is set in `HW1/Q1/src/main.c` (currently 1.0); it is not a Python argument. Theta is not included in the CSV filename, so preserve results separately when comparing theta values.

To run the batch defined in `runall.py`:

```bash
python HW1/Q1/runall.py
```

The batch uses 500 grid points and CFL values 0.8, 1.0, and 1.2, plus 1.8 for `OWN2SPACE4TIME`. It uses the configured theta value for `THETA`.

## Q2 — nozzle flow

Set the algorithm, outlet type, CFL, grid size, and stopping criteria in `HW1/Q2/src/main.c`. Keep the filename labels `algName` and `outletName` consistent with the selected algorithm and outlet (`SUPER` or `SUB`).

```bash
mkdir -p HW1/Q2/build HW1/Q2/plots/data
make -C HW1/Q2
python HW1/Q2/plots/plotting.py
python HW1/Q2/plots/convergence.py
```

The solver stops when its residual tolerance is met or the maximum iteration count is reached. The profile plots show the final saved state, which is not necessarily converged.

To regenerate the report comparisons from the predefined set of saved runs:

```bash
python HW1/Q2/plots/report_plots.py
```

See [the report figure index](HW1/Q2/plots/figures/report/README.md) for the run settings, figure links, and notes on convergence reconstructed from CSV data.

## Saved results

Each question keeps CSV data in `plots/data/` and figures in `plots/figures/`. Filenames identify the scheme, CFL multiplied by 100, and grid size; Q2 also includes the outlet type. Repeating a run with the same filename overwrites its CSV.

Both `plotting.py` scripts skip figures that already exist. Delete the corresponding PNGs before regenerating them from updated data. Q2's `convergence.py` checks file timestamps, while `report_plots.py` regenerates its figures each time.
