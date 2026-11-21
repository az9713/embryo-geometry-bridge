import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path

BASE = Path(__file__).resolve().parent
ROOT = next(p for p in BASE.parents if (p / "PLAN.md").is_file())
CERTS = (
    ("max_contrast_change_small_17", "cebd38784ea3976342517073fe4a62926b00fcf09cb63d23c5259a2989f1a45a", "ad1d74f6282f70842959ffe7c43fcb34181ab93e96ceba10566bb428448a165a"),
    ("max_contrast_change_large_03", "61674691643a7ed3e6ada51ddd6a5adc37db471187a5eff217d867ca9e1cd55f", "8c1c93adccd0f2d0e43e22d7a9e7020ff9d178c397bc87270bdbe3f451d2ecff"),
    ("max_negative_contrast_change_small_17", "665fab8f7b5c822eb0846c39d1f50f67aa151090dffc0ddf4be5abcf77b06b98", "149ed5718bb5563dd0b59358025044ec89bd8d43b0dc121a77c1b04e808f791d"),
    ("max_negative_contrast_change_large_03", "0871d6e6026c6dc2b6bc9837d2d7be75a7064e5694e7cdb1d1e59d7dd2c89558", "05fe679276d5e9daf6ec6d895ae2f17c638340dfd52f0f816d89b9b72cb2774f"),
)
hashed = {}


def verify(path, expected, parse=False):
    path = ROOT / path
    name = path.relative_to(ROOT).as_posix()
    payload = path.read_bytes() if parse or name not in hashed else None
    actual = hashlib.sha256(payload).hexdigest() if payload is not None else hashed[name]
    assert actual == expected, name
    hashed[name] = actual
    return json.loads(payload) if parse else name


results = []
for pilot, native_hash, os_hash in CERTS:
    folder = BASE / "actual_replay" / pilot
    v = verify(folder / "receipt.json", native_hash, parse=True)
    os = verify(folder / "os_exit_observation.json", os_hash, parse=True)
    end = verify(folder / "terminal.json", os["terminal_sha256"], parse=True)
    assert v["status"] == "PASS_INDEPENDENT_EXACT_SPATIAL_ROOT_AND_STRICT_LOCAL_MINIMUM"
    assert v["pilot_id"] == os["pilot_id"] == end["pilot_id"] == pilot
    assert v["pid"] == os["pid"] == end["pid"]
    assert os["status"] == "PASS_OBSERVED_INDEPENDENT_GLOBAL_REDUCTION_TRUE_OS_EXIT"
    assert type(os["observed_exit_code"]) is type(end["exit_code"]) is int and os["observed_exit_code"] == end["exit_code"] == 0
    assert os["receipt_sha256"] == end["receipt_sha256"] == native_hash
    assert v["root_admission_sha256"] == os["root_admission_sha256"] == end["root_admission_sha256"]
    for name, h in v["source_and_actual_sha256"].items():
        verify(name, h)
    r = Q(v["radius_exact"])
    b = v["conservative_compared_bounds"]
    eta, z, lhs, margin, error = map(Q, (b["eta_upper_exact"], b["z_upper_exact"], b["selfmap_lhs_upper_exact"],
                                           b["selfmap_margin_lower_exact"], b["posterior_error_upper_exact"]))
    m, norm, goal = map(Q, (v["whole_box_normalized_H_euclidean_lower_exact"], v["exact_R_infinity_norm"],
                              v["observer_state_design_goal_exact"]))
    assert r == Q(1, 2**34) and 0 <= z < 1 and 0 <= eta and lhs == eta+z*r < r
    assert margin == r-lhs > 0 and error == eta/(1-z) < min(r, goal) and m == (1-z)/norm > 0
    assert b["strict_selfmap"] is True and b["strict_exclusion_of_this_cube"] is False
    assert v["posterior_below_observer_state_design_goal"] is True and v["full_five_gate_goal_achieved"] is False
    item = {
        "pilot_id": pilot, "target_identity": v["target_identity"], "bits": v["precision_bits"],
        "radius_exact": str(r), "selfmap": b, "whole_box_normalized_H_min_exact": str(m),
        "state_observer_design_goal_exact": str(goal), "independent_proof_file": (folder / "receipt.json").relative_to(ROOT).as_posix(),
        "independent_proof_sha256": native_hash, "independent_os_observation_file": (folder / "os_exit_observation.json").relative_to(ROOT).as_posix(),
        "independent_os_observation_sha256": os_hash,
    }
    results.append(item)

dest = BASE / "all_four_readback.json"
with dest.open("x", encoding="utf-8", newline="\n") as f:
    f.write(json.dumps({
        "status": "PASS_FOUR_INDEPENDENT_FIXED_INPUT_EXACT_SPATIAL_ROOTS_AND_MINIMA",
        "scope": "One verified fixed-input finite-element root and normalized strict local minimum per admitted cube. This input does not evaluate experimental compatibility or complete calibration coverage.",
        "independent_exact_spatial_root_count": 4, "stronger_loaded_or_continuum_claim": False,
        "full_five_gate_goal_achieved": False, "cases": results,
        "current_file_binding_count": len(hashed), "current_file_sha256": dict(sorted(hashed.items())),
    }, indent=2)+"\n")
print(hashlib.sha256(dest.read_bytes()).hexdigest(), len(hashed))

