"""Original chart from retained six-packet enclosure-width comparison only."""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent
data = json.loads((OUT / "width_comparison.json").read_text(encoding="utf8"))
rows = data["rows"]
plt.style.use("dark_background")
fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
fig.patch.set_facecolor("#0e1622")
for ax, prefix, title in zip(axes, ("force", "Hessian"), ("Force", "Hessian")):
    ax.set_facecolor("#0e1622")
    x = np.arange(6)
    old = [float(row[f"max_{prefix}_interval_width_old_exact"].split("/")[0]) / float(row[f"max_{prefix}_interval_width_old_exact"].split("/")[1]) for row in rows]
    new = [float(row[f"max_{prefix}_interval_width_new_exact"].split("/")[0]) / float(row[f"max_{prefix}_interval_width_new_exact"].split("/")[1]) for row in rows]
    ax.bar(x - .19, old, .36, color="#edb563", label="Duffy representation")
    ax.bar(x + .19, new, .36, color="#6ed8c0", label="Physical polynomials")
    ax.set_yscale("log")
    ax.set_xticks(x, [str(row["element"]) for row in rows])
    ax.set_xlabel("Full-axis-edge element")
    ax.set_ylabel("Maximum component interval width\n(normalized FE coordinates)")
    ax.set_title(title, fontweight="bold")
    ax.grid(axis="y", alpha=.16)
    ax.set_axisbelow(True)
    ax.legend(loc="upper right", fontsize=9)
fig.suptitle("Local uniform enclosure widths by representation", fontsize=16, fontweight="bold")
fig.text(.5, .025, "Six identical-target contributions. Enclosure widths do not measure physical error or certify a root/PDE bound.", ha="center", fontsize=9)
fig.tight_layout(rect=(0, .07, 1, .92))
for suffix in ("png", "svg"):
    dest = OUT / ("widths." + suffix)
    if dest.exists():
        raise ValueError("preserve original visual")
    fig.savefig(dest, dpi=180, facecolor=fig.get_facecolor())
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
receipt = {"status": "ORIGINAL_SIX_LOCAL_ENCLOSURE_WIDTH_VISUAL", "source_sha256": sha(Path(__file__)), "comparison_sha256": sha(OUT / "width_comparison.json"), "figures": {"widths.png": sha(OUT / "widths.png"), "widths.svg": sha(OUT / "widths.svg")}, "scope": "Original plot of retained local interval widths only, no experimental images, global bound or root inference."}
dest = OUT / "visual_generation.json"
if dest.exists():
    raise ValueError("preserve receipt")
dest.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf8")
print(json.dumps(receipt))
