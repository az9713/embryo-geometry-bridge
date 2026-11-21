"""Original proof-margin/location figures from exact retained fractions.

Publication draft only. This does not rerun any scientific certificate.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[3]
BASE = ROOT / "runs/exact_integration_bridge_v1"
RESULT = BASE / "integral_root_B_safe_v2/results.json"
PROTOCOL = BASE / "integral_root_B_safe_v2/protocol.json"
VECTORS = BASE / "retube_B_preparation/candidate_vectors.npz"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write(p, obj):
    if p.exists():
        raise ValueError("preserve frozen publication preparation")
    p.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf8")


def main():
    if sha(RESULT) != "42d97691311c057d8152f54a6ea31a7218763e7b451df63358bdf423010e18d4":
        raise ValueError("changed intended scalar-composition result")
    result = json.loads(RESULT.read_text(encoding="utf8"))
    protocol = json.loads(PROTOCOL.read_text(encoding="utf8"))
    r = Q(result["r_exact"])
    eI = Q(result["root_error_from_translated_predictor_upper_exact"])
    eQ = Q(result["old_Q_root_error_upper_exact"])
    with np.load(VECTORS, allow_pickle=False) as z:
        decrement = [Q.from_float(float(a)) - Q.from_float(float(b)) for a, b in zip(z["u0"], z["u1"])]
    d = decrement[588]
    dnorm = max(map(abs, decrement))
    common = dnorm + r
    first, refined = result["first_representation_test"], result["physical_axis_test"]
    checks = {
        "exact_result_protocol_binding": result["protocol_sha256"] == sha(PROTOCOL),
        "source_alltheta_sameC_recordB": result["record_id"] == "optoGEF2_small_col_B" and result["parameter_cell_sha256"] == protocol["parameter_cell_sha256"],
        "first_z_below_one_but_selfmap_fails": Q(first["z"]) < 1 and Q(first["selfmap_lhs"]) > r and first["strict_sufficient_conditions"] is False,
        "physical_both_strict_conditions_pass": Q(refined["z"]) < 1 and Q(refined["selfmap_lhs"]) < r and refined["strict_sufficient_conditions"] is True,
        "normalized_selfmap_not_residual_alone": Q(refined["selfmap_lhs"]) == Q(refined["eta"]) + Q(refined["z"]) * r,
        "posterior_error_from_Banach": eI == Q(refined["eta"]) / (1 - Q(refined["z"])),
        "actual_binary_decrement_norm": dnorm == Q(result["d_norm_exact"]),
        "witness588_positive_decrement": d > 0,
        "I_root_interval_outside_old_cube": -d + eI < -r,
        "I_root_interval_inside_common": -common <= -d - eI <= -d + eI <= common,
        "Q_root_interval_inside_old_cube": eQ < r,
        "candidate_cube_inside_common": -common <= -d - r and -d + r <= common,
        "component_gap_matches_exact_result": d - eI - eQ == Q(result["I_Q_branch_component_separation_lower_exact"]),
        "exclusion_distance_matches_exact_result": d - eI == Q(result["integral_root_distance_from_old_predictor_lower_exact"]),
        "curvature_positive_normalized_only": Q(result["uniform_integral_energy_curvature_lower_exact"]) > 0,
        "not_all17_or_PDE_validation": result["all17_integral_or_PDE_or_biological_validation"] is False,
    }
    if not all(v is True for v in checks.values()):
        raise AssertionError(checks)
    plt.style.use("dark_background")
    bg, amber, green, blue = "#0e1622", "#efb869", "#72ddc6", "#8bbdff"
    fig, ax = plt.subplots(figsize=(11, 5.7))
    fig.patch.set_facecolor(bg)
    ax.set_facecolor(bg)
    x = np.arange(2)
    first_values = [Q(first["selfmap_lhs"]) / r, Q(first["z"])]
    refined_values = [Q(refined["selfmap_lhs"]) / r, Q(refined["z"])]
    a = ax.bar(x - .19, list(map(float, first_values)), .36, color=amber, label="First Duffy sufficient bounds")
    b = ax.bar(x + .19, list(map(float, refined_values)), .36, color=green, label="Physical-axis sufficient bounds")
    ax.axhline(1, color="#f28d8d", ls="--", lw=1.7)
    ax.text(1.05, 1.025, "Strict threshold = 1", color="#f0a0a0", fontsize=10)
    ax.set_ylim(0, 1.42)
    ax.set_xticks(x, [r"Self-map load $(\eta_I+zr)/r$", r"Contraction $z$"])
    ax.set_ylabel("Dimensionless proof quantities")
    ax.set_title("Both strict inequalities must pass", fontsize=17, fontweight="bold")
    ax.legend(loc="upper right", framealpha=.92)
    ax.grid(axis="y", alpha=.15)
    ax.set_axisbelow(True)
    for bars in (a, b):
        for rect in bars:
            value = rect.get_height()
            ax.text(rect.get_x() + rect.get_width() / 2, value + .025, f"{value:.6f}", ha="center", fontsize=10)
    fig.text(.5, .025, "Same energy, parameter cell and nodal domain. Display values are rounded from exact fractions; these are not observations.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .07, 1, 1))
    paths = []
    for ext in ("png", "svg"):
        path = OUT / ("root_margins." + ext)
        if path.exists():
            raise ValueError("preserve original figure")
        fig.savefig(path, dpi=180, facecolor=bg)
        paths.append(path)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 6.6))
    fig.patch.set_facecolor(bg)
    ax.set_facecolor(bg)
    data = [
        (4, -common / r, common / r, "#77889c", 10),
        (3, Q(-1), Q(1), amber, 7),
        (2, -eQ / r, eQ / r, blue, 12),
        (1, (-d - r) / r, (-d + r) / r, "#3f8e81", 7),
        (0, (-d - eI) / r, (-d + eI) / r, green, 12),
    ]
    for y, lo, hi, color, thick in data:
        ax.plot([float(lo), float(hi)], [y, y], color=color, lw=thick, solid_capstyle="butt")
        ax.plot([float(lo), float(lo)], [y - .12, y + .12], color=color, lw=2)
        ax.plot([float(hi), float(hi)], [y - .12, y + .12], color=color, lw=2)
    ax.scatter([0, -float(d / r)], [3, 1], marker="D", s=80, facecolors=bg, edgecolors="white", linewidths=1.8, zorder=5)
    ax.axvline(-1, color=amber, ls=":", alpha=.7)
    ax.axvline(1, color=amber, ls=":", alpha=.7)
    ax.axvline(0, color="#bdc7d2", lw=.8, alpha=.35)
    ax.set_xlim(-float(common / r) - .35, float(common / r) + .35)
    ax.set_ylim(-.55, 4.55)
    ax.set_yticks([4, 3, 2, 1, 0], ["Common cube projection", "Old predictor cube", "Cached-Q root enclosure", "Translated predictor cube", "Exact-integral root enclosure"])
    ax.set_xlabel(r"Coordinate $(u_{588}-p_{0,588}(\theta))/r$, for each fixed $\theta\in C$")
    ax.set_title("Proved root enclosures in coordinate 588", fontsize=17, fontweight="bold")
    ax.grid(axis="x", alpha=.13)
    ax.legend(handles=[Line2D([], [], marker="D", color="none", markerfacecolor=bg, markeredgecolor="white", label="Predictor marker (not an equilibrium)")], loc="upper right", fontsize=9)
    ax.text(-float(d / r), -.35, f"Integral enclosure ≈ [{float((-d-eI)/r):.3f}, {float((-d+eI)/r):.3f}]", ha="center", color=green, fontsize=9)
    fig.text(.5, .035, "Intervals are proof enclosures, not measurements or point estimates. Common-cube uniqueness is a separate theorem premise.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .09, 1, 1))
    for ext in ("png", "svg"):
        path = OUT / ("root_location." + ext)
        if path.exists():
            raise ValueError("preserve original figure")
        fig.savefig(path, dpi=180, facecolor=bg)
        paths.append(path)
    plt.close(fig)
    exact_data = {
        "r": str(r), "d588": str(d), "d_norm": str(dnorm), "common_radius": str(common),
        "eI": str(eI), "eQ": str(eQ),
        "first_selfmap_over_r": str(first_values[0]), "physical_selfmap_over_r": str(refined_values[0]),
        "first_z": str(first_values[1]), "physical_z": str(refined_values[1]),
        "I_coordinate_interval_over_r": [str((-d-eI)/r), str((-d+eI)/r)],
        "Q_coordinate_interval_over_r": [str(-eQ/r), str(eQ/r)],
        "I_Q_component_gap_lower": result["I_Q_branch_component_separation_lower_exact"],
    }
    write(OUT / "figure_data_exact.json", exact_data)
    receipt = {
        "status": "ORIGINAL_FIGURES_FROM_SUCCESSFUL_PRODUCER_PENDING_FINAL_INDEPENDENT_REVIEW",
        "checks": checks,
        "bindings": {path.relative_to(ROOT).as_posix(): sha(path) for path in (Path(__file__), RESULT, PROTOCOL, VECTORS, OUT / "figure_data_exact.json")},
        "figures": {path.name: sha(path) for path in paths},
        "scope": "Exact-fraction figure arithmetic and original plot generation only, no certificate replay or publication approval. The producer success awaits final independent review.",
    }
    write(OUT / "figure_generation.json", receipt)
    print(json.dumps({"status": receipt["status"], "checks": len(checks), "figures": receipt["figures"], "receipt_sha256": sha(OUT / "figure_generation.json")}))


if __name__ == "__main__":
    main()
