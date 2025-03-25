import matplotlib.pyplot as plt
import numpy as np

# Hamming Loss values from logs
hamming_loss_values = {
    "Random Forest": 0.002367059176479412,
    "SVM": 0.0023999474042946003,
    "Gradient Boosting": 0.0020867211944592032
}

models = list(hamming_loss_values.keys())
values = list(hamming_loss_values.values())

# Create a clean and subtle line plot with useful axis scaling
fig, ax = plt.subplots(figsize=(6, 4))

ax.plot(models, values, marker="o", linestyle="-", color="black", linewidth=1.5, markerfacecolor="black")

# Labels and Title (subtle styling)
ax.set_title("Hamming Loss Comparison", fontsize=12, fontweight='bold', pad=10)
ax.set_ylabel("Hamming Loss", fontsize=11)
ax.set_xticks(models)
ax.set_xticklabels(models, fontsize=10, rotation=15)

# Set Y-axis ticks dynamically for better readability
ax.set_yticks(np.linspace(min(values), max(values), num=5))
ax.set_yticklabels([f"{y:.5f}" for y in np.linspace(min(values), max(values), num=5)], fontsize=10)

# Light grid for readability
ax.grid(axis='y', linestyle='--', alpha=0.3)

# Add subtle annotations for exact values
for i, v in enumerate(values):
    ax.text(i, v + 0.00002, f"{v:.5f}", ha='center', fontsize=9, color="black")

plt.show()
