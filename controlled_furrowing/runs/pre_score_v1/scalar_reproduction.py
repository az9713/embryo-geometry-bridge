"""A deliberately limited, frozen scalar reproduction; not a 3D shape solver.

Run prepare -> predict -> score. Prediction never loads evaluation responses.
All paths are relative to this study; all receipts record exact file hashes.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = ROOT / path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def rows(path):
    with (ROOT / path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def read_json(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def response(intensity, height, nu, gain, bisections=60):
    """Broadcast inputs; unique x>=1. An analytic upper bound brackets root.

    Both nonnegative terms are increasing for x>=1. For nu<=1/2,
    the first term dominates (x^1.5-1)/2; for nu>=1/2 the second
    dominates (x-1)*x^1.5/2. Therefore x<=1+2*a*I is a common
    (loose) upper bound, since x^1.5-1 >= x-1.
    """
    intensity, gain = np.broadcast_arrays(np.asarray(intensity, float), np.asarray(gain, float))
    if height <= 0 or not 0 <= nu <= 1 or np.any(intensity < 0) or np.any(gain < 0):
        raise ValueError("Nonphysical null-model input")
    target = gain * intensity
    lo = np.ones_like(target)
    hi = 1 + 2 * target
    for _ in range(bisections):
        mid = (lo + hi) / 2
        lhs = (1 - nu) * (mid ** 1.5 - 1) + nu * (mid - 1) * mid ** 1.5
        below = lhs < target
        lo = np.where(below, mid, lo)
        hi = np.where(below, hi, mid)
    x = (lo + hi) / 2
    return height * (1 - 1 / x)


def provenance(paths):
    return {str(Path(p).as_posix()): digest(ROOT / p) for p in paths}


def prepare():
    # This operation accesses responses only to strip them; the restriction
    # concerns the subsequent fit/predict process, not alleged human blinding.
    source = rows("data/derived/evaluation_large_50um.csv")
    inputs = [{"record_id": r["record_id"],
               "activation_magnitude": float(r["activation_magnitude"]),
               "reported_diameter_um": float(r["reported_diameter_um"])} for r in source]
    write_json("data/derived/large_inputs_only.json", inputs)
    print(json.dumps({"prepared_input_records": len(inputs), "response_fields": False}))


def predict():
    protocol = read_json("PROTOCOL.json")
    train = rows(protocol["calibration"])
    inputs = read_json(protocol["evaluation_inputs"])
    intensity = np.array([float(r["activation_magnitude"]) for r in train])
    observed = np.array([float(r["maximum_shortening_um"]) for r in train])
    test_i = np.array([r["activation_magnitude"] for r in inputs])
    setting = protocol["mechanical_null"]
    grid = setting["a_grid"]
    gains = np.linspace(grid["minimum"], grid["maximum"],
                        round((grid["maximum"] - grid["minimum"]) / grid["step"]) + 1)
    fits = []
    for height in [setting["height_primary_um"], *setting["height_sensitivity_um"]]:
        for nu in setting["nu_grid"]:
            values = response(intensity[None, :], height, nu, gains[:, None], setting["root_bisections"])
            rss = np.sum((values - observed) ** 2, axis=1)
            best = int(np.argmin(rss))
            gain = float(gains[best])
            fits.append({"height_um": height, "nu": nu, "a": gain,
                         "calibration_RSS_um2": float(rss[best]),
                         "gain_at_boundary": best in (0, len(gains)-1),
                         "training_prediction_um": values[best].tolist(),
                         "evaluation_prediction_um": response(test_i, height, nu, gain).tolist()})
    primary_fits = [f for f in fits if f["height_um"] == setting["height_primary_um"]]
    selected = min(primary_fits, key=lambda f: f["calibration_RSS_um2"])
    design = np.column_stack([np.ones(len(intensity)), intensity])
    coefficients = np.linalg.lstsq(design, observed, rcond=None)[0]
    outputs = {
        "status": "FROZEN scalar predictions, not geometric predictions",
        "scope": protocol["scope"],
        "python": platform.python_version(), "numpy": np.__version__,
        "files_read": ["PROTOCOL.json", protocol["calibration"], protocol["evaluation_inputs"]],
        "hashes": provenance(["PROTOCOL.json", "scripts/scalar_reproduction.py", protocol["calibration"], protocol["evaluation_inputs"]]),
        "evaluation_record_ids": [r["record_id"] for r in inputs],
        "evaluation_intensity": test_i.tolist(),
        "finite_grid_fit_count": len(fits) * len(gains),
        "primary_mechanical": selected,
        "material_height_sensitivity": fits,
        "OLS": {"intercept_um": float(coefficients[0]), "slope_um": float(coefficients[1]),
                "training_prediction_um": (design @ coefficients).tolist(),
                "evaluation_prediction_um": (coefficients[0] + coefficients[1] * test_i).tolist()},
        "zero_shortening": {"evaluation_prediction_um": np.zeros(len(test_i)).tolist()},
        "geometric_prediction": None,
        "uncertainty": protocol["uncertainty"],
    }
    path = ROOT / "runs/frozen_scalar_predictions.json"
    if path.exists():
        raise RuntimeError("Frozen predictions already exist; do not silently overwrite a scored protocol")
    write_json("runs/frozen_scalar_predictions.json", outputs)
    write_json("evidence/scalar_prediction_lock.json", {
        "prediction_sha256": digest(path), "hashes": outputs["hashes"],
        "responses_loaded_by_prediction": False,
        "scope": "Computational data-flow separation; authors already know publication outcomes",
        "evaluation_response_path_exists": (ROOT / protocol["evaluation_responses"]).exists(),
    })
    print(json.dumps({"status": outputs["status"], "primary": {k: selected[k] for k in ("height_um", "nu", "a", "calibration_RSS_um2")}, "finite_grid_fit_count": outputs["finite_grid_fit_count"], "sha256": digest(path)}))


def score():
    frozen = read_json("runs/frozen_scalar_predictions.json")
    lock = read_json("evidence/scalar_prediction_lock.json")
    if digest(ROOT / "runs/frozen_scalar_predictions.json") != lock["prediction_sha256"]:
        raise RuntimeError("Frozen predictions changed")
    for path, expected in frozen["hashes"].items():
        if digest(ROOT / path) != expected:
            raise RuntimeError(f"Changed preregistered file: {path}")
    evaluation = rows("data/derived/evaluation_large_50um.csv")
    if [r["record_id"] for r in evaluation] != frozen["evaluation_record_ids"]:
        raise RuntimeError("Evaluation identity/order mismatch")
    observed = np.array([float(r["maximum_shortening_um"]) for r in evaluation])
    scores = {}
    for name in ["primary_mechanical", "OLS", "zero_shortening"]:
        estimate = np.array(frozen[name]["evaluation_prediction_um"])
        residual = estimate - observed
        scores[name] = {"RMSE_um": float(np.sqrt(np.mean(residual ** 2))),
                        "MAE_um": float(np.mean(np.abs(residual))),
                        "signed_bias_um": float(np.mean(residual)),
                        "residuals_um": residual.tolist()}
    result = {"status": "SCORED auxiliary scalar reproduction only", "n_calibration": 17,
              "n_evaluation": len(evaluation), "prediction_sha256": lock["prediction_sha256"],
              "evaluation_sha256": digest(ROOT / "data/derived/evaluation_large_50um.csv"),
              "scores": scores, "geometric_gate": "NOT PASSED; no coupled geometry prediction exists",
              "novelty": "None established: source Fig4h already compares this size dependence",
              "statistical_limit": "No independence/coverage guarantee, equivalence claim or biological mechanism identification"}
    write_json("evidence/scalar_scores.json", result)
    plot(frozen, evaluation, scores)
    print(json.dumps(result))


def plot(frozen, evaluation, scores):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("dark_background")
    train = rows("data/derived/calibration_small_18um.csv")
    i = np.linspace(0, 1.6, 200)
    f = frozen["primary_mechanical"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.9), constrained_layout=True)
    ax = axes[0]
    ax.scatter([float(r["activation_magnitude"]) for r in train], [float(r["maximum_shortening_um"]) for r in train], label="18 µm calibration (17)", color="#87cdea")
    ax.scatter([float(r["activation_magnitude"]) for r in evaluation], [float(r["maximum_shortening_um"]) for r in evaluation], label="50 µm evaluation (9)", marker="s", color="#f6bb70")
    ax.plot(i, response(i, f["height_um"], f["nu"], f["a"]), color="#91d4a8", label="Frozen independent-column null")
    ax.plot(i, frozen["OLS"]["intercept_um"] + frozen["OLS"]["slope_um"] * i, "--", color="#e2b5eb", label="Calibration-only linear comparator")
    ax.set(xlabel="Published normalized optical recruitment", ylabel="Maximum shortening (µm)", title="Scalar reproduction; not a shape prediction")
    ax.legend(fontsize=8)
    ax = axes[1]
    names = ["primary_mechanical", "OLS", "zero_shortening"]
    ax.bar(["Column null", "Linear", "Zero"], [scores[n]["RMSE_um"] for n in names], color=["#91d4a8", "#e2b5eb", "#999999"])
    ax.set(ylabel="50 µm evaluation RMSE (µm)", title="All declared comparators, no refitting")
    fig.suptitle("Known size-response comparison: a bounded reproduction, no novelty claim", fontsize=12)
    path = ROOT / "docs/assets/scalar_reproduction.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["prepare", "predict", "score"])
    args = parser.parse_args()
    {"prepare": prepare, "predict": predict, "score": score}[args.stage]()
