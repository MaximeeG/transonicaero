import subprocess
import sys
from pathlib import Path

# SIMULATION SETTINGS
# algorithm = "LAX"
algorithm = "LAX_WENDROFF"
# algorithm = "FORWARD"
# algorithm = "BACKWARD"
cfl = 0.5
nx = 10000


# path of this directory
root = Path(__file__).resolve().parent

print("Compiling and running C solver...")

# compile
subprocess.run(
    ["make"]
)

# execute
subprocess.run(
    [
        "./HW1/build/output",
        algorithm,
        str(cfl),
        str(nx)
    ],
    check=True
)

print("Solver finished.")
print("Running plotting script...")

subprocess.run(
    [sys.executable, "HW1/plots/plotting.py"],
    cwd=root,
    check=True
)

print("Done.")