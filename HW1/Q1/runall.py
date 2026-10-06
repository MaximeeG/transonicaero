import subprocess
import sys
from pathlib import Path

# 1. Setup paths and directories
root = Path(__file__).resolve().parent
(root / "build").mkdir(exist_ok=True)
(root / "plots" / "data").mkdir(parents=True, exist_ok=True)

# 2. Define Simulation Parameters
nx = 500
standard_cfls = [0.8, 1.0, 1.2]

# Dictionary defining each scheme and its specific CFL requirements.
runs = {
    "BACKWARD": standard_cfls,
    "FORWARD": standard_cfls,
    "LAX": standard_cfls,
    "LAX_WENDROFF": standard_cfls,
    "LEAPFROG": standard_cfls,
    "OWN2SPACE4TIME": [0.8, 1.0, 1.2, 1.8], 
    "OWN4SPACE2TIME": standard_cfls,
    "THETA": standard_cfls
}

# 3. Compile the C code exactly once
print("Compiling C solver...")
try:
    subprocess.run(["make"], cwd=root, check=True)
except subprocess.CalledProcessError as e:
    print(f"Compilation failed. Exiting. Error: {e}")
    sys.exit(1)

# 4. Execute all parameter combinations directly through the C binary
print("Starting simulations...")
for scheme, cfl_list in runs.items():
    for cfl in cfl_list:
        
        # Construct the command array strictly matching main.c (3 arguments)
        cmd = [
            "./build/output", 
            scheme, 
            str(cfl), 
            str(nx)
        ]
        
        print(f"Executing: {' '.join(cmd)}")
        
        try:
            subprocess.run(cmd, cwd=root, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Run failed for {scheme} at CFL {cfl}. Error: {e}")

# 5. Run the plotting script exactly once at the end
print("All simulations finished. Running plotting script...")
try:
    subprocess.run([sys.executable, "plots/plotting.py"], cwd=root, check=True)
    print("Done.")
except subprocess.CalledProcessError as e:
    print(f"Plotting script failed. Error: {e}")