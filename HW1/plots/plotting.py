import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize

# --------------------------------------------------
# Load data
# --------------------------------------------------
filename = "HW1/plots/WAVE_BACKWARD.csv"

times = []
u_vectors = []

with open(filename, "r") as f:
    next(f)  # skip header: t,u

    for line in f:
        t_str, u_str = line.strip().split(",", 1)
        t = float(t_str)
        u = np.fromstring(u_str.strip().strip("[]"), sep=",")
        times.append(t)
        u_vectors.append(u)

times = np.array(times)
U = np.array(u_vectors)   # shape: (num_times, num_x_points)

# x values over [0, pi]
x = np.linspace(0, np.pi, U.shape[1])

c = 1.0  # change if your simulation uses a different wave speed

def u_exact(x, t):
    x0 = x - c*t
    return np.where((x0 >= 0.5) & (x0 <= 1.0), 1.0, 0.0)

# --------------------------------------------------
# Color map based on time
# --------------------------------------------------
norm = Normalize(vmin=times.min(), vmax=times.max())
# cmap = cm.rainbow   # you can also try cm.viridis, cm.plasma, cm.turbo
# cmap = cm.plasma
cmap = cm.viridis
# cmap = cm.turbo

# --------------------------------------------------
# 3D plot
# --------------------------------------------------
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection="3d")

for i, t in enumerate(times):
    ax.plot(
        x,
        np.full_like(x, t),
        U[i],
        color=cmap(norm(t)),
        linewidth=1.5
    )

ax.set_xlabel("x")
ax.set_ylabel("Time t")
ax.set_zlabel("u(x,t)")
ax.set_title("Evolution of u(x,t) with Color Showing Time")

# Add colorbar to show time-to-color mapping
mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
mappable.set_array([])
cbar = fig.colorbar(mappable, ax=ax, pad=0.1, shrink=0.7)
cbar.set_label("Time t")

plt.tight_layout()
plt.savefig("HW1/plots/WAVE_BACKWARD_3D.png", dpi=300, bbox_inches="tight")
plt.close()

print("Saved as wave_3d_rainbow.png")