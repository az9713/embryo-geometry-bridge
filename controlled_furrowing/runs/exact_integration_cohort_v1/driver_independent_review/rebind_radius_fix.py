"""Rebind the L-only prelaunch gate after the audited radius enclosure correction."""

import argparse
import difflib
import json
from fractions import Fraction as Q
from pathlib import Path

import check_driver_gate as gate_check
import review_repaired as prior_review
import review_prelaunch as exact_check

P = prior_review.P
HERE = Path(__file__).resolve().parent
OUT = prior_review.OUT
ARCHIVE = exact_check.ROOT / "runs/exact_integration_cohort_v1/driver_preparation/prelaunch_revision_2"
RADIUS_REVIEW = HERE / "cohort_radius_ball_containment.json"
RADIUS_SHA = "6a55576790a4ae39a499f3802c76215827099d82c8834da4437541d08fc101f0"
PRIOR_GATE_SHA = "92d8c5317f718aa4d74c5939df8a96d4a02d1f607ccc79bce8d615380fd3de9e"
PRIOR_BOX_SHA = "161ba86e1004626417b627deb213c72bc1f5f75b6e587fe0a89b7dcbda152761"


def driver_diff():
    old = (ARCHIVE / "cohort_exact_pilot.py").read_text().splitlines()
    new = P.SOURCE.read_text().splitlines()
    diff = list(difflib.unified_diff(old, new, n=0))
    removed = [line[1:].strip() for line in diff[2:] if line.startswith("-")]
    added = [line[1:].strip() for line in diff[2:] if line.startswith("+")]
    assert len(added) == len(removed) == 1
    assert removed[0] == "if F(bounds['radius']['upper_exact'])!=context.old_radius:raise ValueError('audited radius must equal selected original radius')"
    assert added[0].startswith("if not F(bounds['radius']['lower_exact'])<=context.old_radius<=F(bounds['radius']['upper_exact']):raise ValueError(")
    return diff


def main(source_sha, protocol_sha, preflight_sha):
    gate_file = OUT / "domains_prelaunch_independent_review.json"
    assert not gate_file.exists() and not (OUT / "domains_started.json").exists()
    current_files = [P.SOURCE, OUT / "protocol.json", OUT / "preflight.json"]
    assert [P.sha(f) for f in current_files] == [source_sha, protocol_sha, preflight_sha]
    old_gate_file = HERE / "repaired_L_gate_before_radius_refreeze.json"
    assert P.sha(old_gate_file) == PRIOR_GATE_SHA
    old_gate = P.read(old_gate_file)
    assert old_gate["status"] == "PASS" and P.alltrue(old_gate["checks"])
    assert P.sha(P.ROOT / old_gate["support_file"]) == old_gate["support_sha256"]
    old_source = ARCHIVE / "cohort_exact_pilot.py"
    assert P.sha(old_source) == prior_review.SOURCE
    assert P.sha(ARCHIVE / "protocol.json") == prior_review.PROTOCOL
    assert P.sha(ARCHIVE / "preflight.json") == prior_review.PREFLIGHT
    diff = driver_diff()
    old = P.read(ARCHIVE / "protocol.json")
    current = P.read(OUT / "protocol.json")
    assert P.sha(ARCHIVE / "coordinate_box.json") == P.sha(OUT / "coordinate_box.json") == PRIOR_BOX_SHA
    old_hashes = old.pop("source_and_input_hashes")
    new_hashes = current.pop("source_and_input_hashes")
    assert old_hashes.pop("scripts/cohort_exact_pilot.py") == prior_review.SOURCE
    assert new_hashes.pop("scripts/cohort_exact_pilot.py") == source_sha
    assert old_hashes == new_hashes
    assert old.pop("utc") != current.pop("utc")
    assert old == current
    preflight = P.read(OUT / "preflight.json")
    old_preflight = P.read(ARCHIVE / "preflight.json")
    for v, s, p in ((preflight, source_sha, protocol_sha),
                    (old_preflight, prior_review.SOURCE, prior_review.PROTOCOL)):
        assert v.pop("source_sha256") == s and v.pop("protocol_sha256") == p
    assert preflight == old_preflight
    assert P.sha(RADIUS_REVIEW) == RADIUS_SHA
    radius_receipt = P.read(RADIUS_REVIEW)
    assert radius_receipt["status"] == "PASS" and radius_receipt["record_count"] == 17
    assert len({r["record_id"] for r in radius_receipt["records"]}) == 17
    assert all(Q(r["lower_exact"]) <= Q(r["selected_radius_exact"]) <= Q(r["upper_exact"])
               for r in radius_receipt["records"])
    checks = {
        **old_gate["checks"],
        "archived_repaired_L_source_protocol_and_gate_remain_exact": True,
        "only_radius_enclosure_guard_changed_since_repaired_surface_smoke": True,
        "current_preparation_keeps_exact_L_math_box_and_input_hashes": True,
        "all_17_exact_selected_radii_lie_in_independent_radius_enclosures": True,
    }
    guards = gate_check.verify(OUT)
    assert guards["status"] == "PASS" and all(v is True for v in guards["checks"].values())
    support = {"driver_radius_guard_diff": diff, "previous_gate_sha256": PRIOR_GATE_SHA,
               "cohort_radius_enclosure_receipt_sha256": RADIUS_SHA, "current_gate_checks": guards}
    support_file = HERE / "final_radius_refreeze_checks.json"
    prior_review.write_new(support_file, support)
    gate = {
        "status": "PASS", "checks": checks, "source_sha256": source_sha,
        "protocol_sha256": protocol_sha, "preflight_sha256": preflight_sha,
        "record_id": exact_check.RECORD,
        "support_file": support_file.relative_to(P.ROOT).as_posix(),
        "support_sha256": P.sha(support_file),
        "prior_reviewed_surface_smoke_sha256": old_gate["support_sha256"],
        "review_source_sha256": P.sha(__file__),
        "scope": "One scoped prelaunch only, for L's original exact domain and selector job. "
                 "No full job was launched by the reviewer.",
    }
    prior_review.write_new(gate_file, gate)
    P.gate(OUT, "domains")
    print(json.dumps({"status": "PASS", "checks": len(checks), "gate_sha256": P.sha(gate_file),
                      "support_sha256": P.sha(support_file)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("source_sha")
    parser.add_argument("protocol_sha")
    parser.add_argument("preflight_sha")
    args = parser.parse_args()
    main(args.source_sha, args.protocol_sha, args.preflight_sha)
