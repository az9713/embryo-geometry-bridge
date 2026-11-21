"""Original mathematical schematic; no microscopy data are used."""
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent
x = np.linspace(-1, 1, 400)
basal = 1 + 0.1 * np.cos(np.pi * x)
apical = [-0.2 - 0.1 * np.cos(np.pi * x), -0.65 - 0.15 * np.cos(np.pi * x)]
plt.style.use("dark_background")
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
for axis, title, curve in zip(axes, ["Latent geometry A", "Latent geometry B"], apical):
    axis.axhspan(-1, 0, color="#313b4a", alpha=0.9, label="Outside observed crop")
    axis.axhline(0, color="#8b99aa", ls="--", label="Crop edge, not anatomy")
    axis.plot(x, curve, color="#f5bc66", label="Unobserved apical boundary")
    axis.plot(x, basal, color="#76d4e9", lw=2, label="Same observed basal boundary")
    axis.set_title(title)
    axis.set_xlabel("Arbitrary lateral coordinate")
    axis.set_ylim(1.4, -1)
axes[0].set_ylabel("Downward image coordinate (arbitrary units)")
axes[1].legend(loc="lower right", fontsize=7)
fig.suptitle("Boundary cropping can hide different heights\n"
             "Original schematic: no source image or biological measurement", fontsize=11)
fig.tight_layout()
fig.savefig(root / "censoring_nonidentifiability_schematic.png", dpi=160)
plt.close(fig)
