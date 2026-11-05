import numpy as np
import matplotlib.pyplot as plt

# ===== Settings =====
filename = "file_abundance.txt"

# ===== Load data =====
data = np.loadtxt(filename)

step = data[:, 0]
amino_acid = data[:, 1]
sugar = data[:, 2]
purine = data[:, 3]

# ===== Plot =====
fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(step, amino_acid, label="Amino acid",
        color="tab:blue", linewidth=2)

ax.plot(step, sugar, label="Sugar",
        color="tab:green", linewidth=2)

ax.plot(step, purine, label="Purine core",
        color="tab:red", linewidth=2)

# ===== Log scale =====
ax.set_yscale("symlog", linthresh=1e-4)

# ===== Style =====
ax.set_xlabel("Step", fontsize=14)
ax.set_ylabel("Abundance", fontsize=14)
ax.set_title("Molecular Abundance vs Step", fontsize=16)

ax.legend(fontsize=12)
ax.grid(True, which="both", alpha=0.3)
ax.tick_params(labelsize=12)

plt.tight_layout()

# ===== Save =====
plt.savefig("abundance_plot.png", dpi=300)

# ===== Display =====
plt.show()