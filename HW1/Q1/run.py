import subprocess
import sys
from pathlib import Path

# SIMULATION SETTINGS
# algorithm = "LAX"
# algorithm = "LAX_WENDROFF"
# algorithm = "FORWARD"
# algorithm = "BACKWARD"
# algorithm = "LEAPFROG"
algorithm = "THETA"

cfl = 0.5
nx = 500


# path of this directory
root = Path(__file__).resolve().parent
(root / "build").mkdir(exist_ok=True)
(root / "plots" / "data").mkdir(parents=True, exist_ok=True)

print("Compiling and running C solver...")

# compile
subprocess.run(
    ["make"],
    cwd=root,
    check=True
)

# execute
subprocess.run(
    [
        "./build/output",
        algorithm,
        str(cfl),
        str(nx)
    ],
    cwd=root,
    check=True
)

print("Solver finished.")
print("Running plotting script...")

subprocess.run(
    [sys.executable, "plots/plotting.py"],
    cwd=root,
    check=True
)

print("Done.")
