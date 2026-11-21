"""Record the sparse-Arb cleanup counterexample without rerunning science."""
from pathlib import Path
import datetime
import hashlib
import json
import sys

import flint
from flint import arb

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
import exact_fe_target as target


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    output = ROOT / "evidence/exact_integration_final_delivery_independent_review.json"
    if output.exists():
        raise RuntimeError("write-once output already exists")
    x = arb(0, 1)
    old = target.clean({(0, 0): x})
    checks = {
        "nonexact_Arb_zero_crossing_coefficient": not x.is_zero(),
        "Arb_ne_zero_is_false_for_overlapping_ball": (x != 0) is False,
        "Arb_eq_zero_is_also_false": (x == 0) is False,
        "actual_frozen_clean_drops_nonexact_ball": old == {},
        "exact_zero_test_distinguishes_exact_zero": arb(0).is_zero(),
    }
    files = [
        "scripts/exact_fe_target.py",
        "runs/exact_integration_bridge_v1/bridge_design_probe/polynomial_jet.py",
        "runs/exact_integration_bridge_v1/surface_pilot_independent_replay/replay.py",
        "evidence/exact_integration_bridge_final_delivery.json",
        "evidence/exact_complete_point_B_independent_review.json",
        "evidence/exact_integration_documentation_review.json",
    ]
    receipt = {
        "status": "WITHDRAWN_SCIENTIFIC_APPROVAL_PENDING_CORRECTED_ARITHMETIC_AND_REPLAY",
        "recorded_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "checks": checks,
        "finding": {
            "source": "scripts/exact_fe_target.py:34",
            "expression": "if v != 0",
            "counterexample": str(x),
            "consequence": "An interval containing zero but not identically zero is silently removed from polynomial coefficients. This operation need not enclose the original polynomial.",
            "required_fix": "Discard only proven exact-zero coefficients; retain nonexact zero-crossing balls. Version and replay all affected polynomial-enclosure computations and dependent conclusions.",
        },
        "bound_files": {name: sha(ROOT / name) for name in files},
        "review_source_sha256": sha(Path(__file__)),
        "python": sys.executable,
        "python_flint_version": flint.__version__,
        "affected_claims": [
            "Point B exact-integral/cache force and Hessian discrepancy enclosures derived through this cleanup",
            "Dependent negative direct-self-map numerical certificate and delivered scientific approval",
            "Our six-edge replay is affected by the shared target arithmetic; agreement is not an independent safeguard against this common defect",
            "Our mathematical documentation PASS is a historical conditional content review, not current approval of its affected numerical premises",
        ],
        "not_inferred": [
            "No conclusion that every affected numerical value is numerically false",
            "No invalidation of unrelated exact-Fraction geometry/moment identities or old cached-Q dense root certificates solely from this defect",
            "No completion of original five gates, no biological refutation, no new scientific result",
        ],
        "delivery_binding_audit_completed": False,
        "scientific_approval": False,
        "goal_achieved": False,
        "scope": "Counterexample and dependency inspection only; no corrected integration, full-proof replay, production document edit, browser action or Space edit.",
    }
    if not all(checks.values()):
        raise RuntimeError("counterexample checks did not reproduce")
    output.write_text(json.dumps(receipt, indent=2, ensure_ascii=True) + "\n", encoding="utf8")
    print(json.dumps({"status": receipt["status"], "checks": len(checks), "receipt_sha256": sha(output), "bindings": receipt["bound_files"]}, ensure_ascii=True))


if __name__ == "__main__":
    main()
