# transonicaero

Code for the Transonic Aerodynamics assignments. So far, HW1 contains a C solver for the 1D advection equation and a Python script to plot the results.

## Directory structure

```text
transonicaero/
├── run.py                  # compile, run the solver, and generate plots
├── makefile                # GCC build command
└── HW1/
    ├── MEC6602E_HW1.pdf     # assignment statement
    ├── src/
    │   ├── main.c          # arguments, simulation setup, and CSV output
    │   ├── wavesolver.c    # initial condition and numerical schemes
    │   └── wavesolver.h    # structs and function declarations
    ├── build/
    │   └── output          # compiled solver
    ├── plots/
    │   ├── plotting.py     # numerical vs analytical solution plots
    │   ├── data/           # solver output (.csv)
    │   └── figures/        # generated plots (.png)
    └── .venv/              # local Python environment, if created
```

## Running it

You need Python 3, GCC, and Make. The plotting script also needs NumPy and Matplotlib.
Run the commands below from the `transonicaero` directory; the build and solver use relative paths.

For the first run, set up the Python environment and output folders:

```bash
python3 -m venv HW1/.venv
source HW1/.venv/bin/activate
python -m pip install numpy matplotlib
mkdir -p HW1/build HW1/plots/data HW1/plots/figures
```

For later sessions, just activate the environment again:

```bash
source HW1/.venv/bin/activate
```

Change the settings at the top of `run.py` before running. The current defaults are:

```python
algorithm = "LAX_WENDROFF"
cfl = 0.5
nx = 10000
```

- `algorithm`: `BACKWARD`, `FORWARD`, `LAX`, or `LAX_WENDROFF`. `LEAPFROG` is recognised by the argument parser but is not implemented yet.
- `cfl`: CFL number, which sets the time step through `dt = cfl * dx / c`. Use a positive value.
- `nx`: number of grid points, at least 2. Try 100 or 500 for a quick run; 10000 produces a much larger CSV.

Then run:

```bash
python run.py
```

There are no command-line options for `run.py`; it reads the settings from the file. It calls `make`, runs `HW1/build/output` with those settings, then runs the plotting script using the same Python environment.

## Results

CSV files go into `HW1/plots/data/`. For example, `WAVE_LAX_WENDROFF_CFL_50_NX_10000.csv` is the output for Lax–Wendroff with CFL 0.5 and 10000 grid points. The CFL value in the filename is multiplied by 100 and converted to an integer.

The plotting script checks all CSV files in that folder and saves two figures per case in `HW1/plots/figures/`: a 3D view of the wave over time and a comparison with the analytical solution. Cases with non-finite or extremely large solution values are skipped.

Running the same settings again overwrites the CSV, but plotting is skipped if both figures already exist. Delete the corresponding PNGs if you want to regenerate them.

To run only the plotting step:

```bash
python HW1/plots/plotting.py
```
