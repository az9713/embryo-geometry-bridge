"""Admit the complete scaled L manifest, then reduce exact signed intersections."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement
import json
import os
from pathlib import Path
import sys
from time import monotonic


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[4]
sys.path.insert(0, str(ROOT / "scripts"))
import review_cohort_exact_record as D


STAGE = OUT.parent.relative_to(ROOT).as_posix()
RECORD = "optoGEF2_small_col_L"
ORIGINAL = D._record_dir(RECORD)
OLD = D.BASE + "/reduction_preparation/" + RECORD
AXIS = [*range(6), *range(72, 78)]
OFFAXIS = [i for i in range(144) if i not in AXIS]
PINS = {
    f"{STAGE}/preparation/protocol.json": "1c07ccbebcbed226c3f8eeb29c707ea0d8524c516de74f86f22dd8267db5db61",
    f"{STAGE}/preparation/preflight.json": "54b8776d29337d872287f35aa73f57ae6ccd40505cab848ebf9c04fb79341be3",
    f"{STAGE}/preparation/selectors.json": "f95c210036c3fbc13e1acb25477c58fa1b785f77e90ec2b93f3b6fee3e60a2be",
    f"{STAGE}/preparation/prelaunch_independent_review.json": "29f1c54c70eda884b3aa150b79cad708d304bac70049649939753139d291d18c",
    OLD + "_input_admission.json": "13fcd8e1a469d20a49ec050d985c7ffed8af4989af6229d7e45d9b15f13b9735",
    OLD + "_reduction_review.json": "5a55e0d139edb39206ab3999681d9858d6b0162f73c7e5103d7a5bed3ebc4514",
    "scripts/review_cohort_exact_record.py": "67db748a5d40f31471ddd3d6c84fae280ff7ead650e6eb4e950a4766e54e0285",
}


def _relative(path):
    return Path(path).relative_to(ROOT).as_posix()


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _packed(lo, hi):
    return {"lower_exact": str(lo), "upper_exact": str(hi)}


def intersection(a, b):
    """Intersection of two enclosing intervals in exact signed endpoint arithmetic."""
    (al, ah), (bl, bh) = D._packed(a), D._packed(b)
    lower, upper = max(al, bl), min(ah, bh)
    D.require(lower <= upper, "scaled and original certified interval intersection is empty")
    return _packed(lower, upper)


def _validate_manifest(manifest, ids, expected_names, fields):
    D.require(all(manifest.get(key) == value for key, value in fields.items()),
              "scaled actual manifest does not bind the reviewed preparation")
    rows = manifest["files"]
    D.require(len(rows) == len(ids) == len(expected_names)
              and [p["element"] for p in rows] == list(ids)
              and all(type(p["element"]) is int for p in rows)
              and [p["file"] for p in rows] == list(expected_names)
              and len({p["file"] for p in rows}) == len(rows),
              "exactly one canonical packet for every off-axis element required")
    for p in rows:
        D._sha(p["sha256"])
    return rows


def _bound_preparation():
    hashes = {}
    docs = {p: D._bound(p, h, hashes, parse=p.endswith(".json")) for p, h in PINS.items()}
    scale = docs[f"{STAGE}/preparation/protocol.json"]
    selectors = docs[f"{STAGE}/preparation/selectors.json"]
    old_admission = docs[OLD + "_input_admission.json"]
    old_review = docs[OLD + "_reduction_review.json"]
    prelaunch = docs[f"{STAGE}/preparation/prelaunch_independent_review.json"]
    p = scale["source_file"]
    D._bound(p, scale["source_sha256"], hashes, parse=False)
    D._bound(prelaunch["support_file"], prelaunch["support_sha256"], hashes, parse=False)
    expected = {
        "source_sha256": hashes[p], "protocol_sha256": PINS[f"{STAGE}/preparation/protocol.json"],
        "preflight_sha256": PINS[f"{STAGE}/preparation/preflight.json"],
        "selectors_sha256": PINS[f"{STAGE}/preparation/selectors.json"],
        "prelaunch_sha256": PINS[f"{STAGE}/preparation/prelaunch_independent_review.json"],
    }
    D.require(prelaunch["status"] == "PASS" and D._yes(prelaunch["checks"])
              and all(prelaunch[k] == v for k, v in expected.items() if k != "prelaunch_sha256")
              and prelaunch["record_id"] == scale["record_id"] == old_admission["record_id"] == old_review["record_id"] == RECORD
              and scale["offaxis_element_indices"] == [v["element"] for v in selectors] == OFFAXIS
              and scale["original_axis_element_indices"] == AXIS, "same reviewed complete off-axis L preparation required")
    D.require(old_review["status"] == "PASS_REDUCTION_UNRESOLVED_SUFFICIENT_ROOT_TEST"
              and old_review["component_count"] == 168 and old_review["force_entries"] == 611
              and old_review["upper_Hessian_entries"] == 6676
              and old_review["input_admission_sha256"] == PINS[OLD + "_input_admission.json"],
              "same previously reviewed original endpoint reduction required")
    old_bulk = old_admission["expected_packets"]["bulk"]
    old_surface = old_admission["expected_packets"]["surface"]
    D._manifest(f"{ORIGINAL}/bulk", old_bulk, [f"element_{i:03d}.json" for i in range(144)])
    D._manifest(f"{ORIGINAL}/surface", old_surface, [f"face{f}_edge{i:02d}.json" for f in (0, 1) for i in range(12)])
    D.require({x["file"]: x["sha256"] for x in old_bulk} == scale["original_bulk_packet_hashes"]
              and all(old_review["source_input_and_component_hashes"][p["file"]] == p["sha256"]
                      for p in old_bulk + old_surface), "matching original 168 packet hashes required")
    for path, digest in scale["input_hashes"].items():
        D._bound(path, digest, hashes, parse=False)
    return hashes, scale, selectors, old_admission, old_review, expected


def _matches_prelaunch(review_sha, hashes):
    path = _relative(OUT / "prelaunch_independent_review.json")
    gate = D._bound(path, review_sha, hashes)
    source = _sha(__file__)
    protocol, preflight = (_sha(OUT / f"{name}.json") for name in ("protocol", "preflight"))
    D.require(gate["status"] == "PASS" and D._yes(gate["checks"])
              and (gate["source_sha256"], gate["protocol_sha256"], gate["preflight_sha256"])
                  == (source, protocol, preflight), "matching independent prelaunch review of this reduction required")
    local_protocol = json.loads((OUT / "protocol.json").read_bytes())
    pre = json.loads((OUT / "preflight.json").read_bytes())
    D.require(local_protocol["source_sha256"] == pre["source_sha256"] == source
              and pre["protocol_sha256"] == protocol and D._yes(pre["checks"]),
              "frozen source and preparation prerequisites required")
    hashes.update({_relative(__file__): source, _relative(OUT / "protocol.json"): protocol,
                   _relative(OUT / "preflight.json"): preflight})


def admit(review_sha, terminal_sha, manifest_sha):
    hashes, scale, _, _, _, expected = _bound_preparation()
    _matches_prelaunch(review_sha, hashes)
    term_path = f"{STAGE}/actual/terminal.json"
    manifest_path = f"{STAGE}/actual/results.json"
    terminal = D._bound(term_path, terminal_sha, hashes)
    D._terminal(terminal, manifest_sha)
    manifest = D._bound(manifest_path, manifest_sha, hashes)
    fields = {**expected, "record_id": RECORD, "status": "ALL_132_SCALED_OFFAXIS_PACKETS_UNREVIEWED",
              "unchanged_axis_elements": AXIS, "original_surface_packets_retained": 24, "global_reduction_run": False}
    packets = _validate_manifest(manifest, OFFAXIS, scale["actual_expected_files"], fields)
    receipt = {"status": "BOUND_COMPLETE_SCALED_MANIFEST_BEFORE_LOCAL_READ", "record_id": RECORD,
               "review_sha256": review_sha, "scaled_terminal_sha256": terminal_sha, "scaled_manifest_sha256": manifest_sha,
               "inputs_and_manifests": hashes, "expected_scaled_packets": packets,
               "local_scaled_payloads_opened": False}
    digest = D._write(OUT / "input_admission.json", receipt)
    print(json.dumps({"admission_sha256": digest, "scaled_packets": len(packets)}), flush=True)


def _entry_intersections(old, new, ids):
    n = len(ids)
    f0, f1 = (p["force_integral_minus_cache"] for p in (old, new))
    H0, H1 = (D._entry_map(p["H_integral_minus_cache_upper"], n) for p in (old, new))
    exact_keys = set(combinations_with_replacement(range(n), 2))
    D.require(len(f0) == len(f1) == n and H0.keys() == H1.keys() == exact_keys,
              "full original and scaled local F/H coverage required")
    f = [intersection(a, b) for a, b in zip(f0, f1)]
    H = [{"i": i, "j": j, "enclosure": intersection(H0[i, j]["enclosure"], H1[i, j]["enclosure"])}
         for i, j in sorted(exact_keys)]
    D.require(len(new["force_intersection_with_original"]) == n
              and list(map(D._packed, new["force_intersection_with_original"])) == list(map(D._packed, f)),
              "reported scaled force intersection differs from signed endpoints")
    declared = D._entry_map(new["H_intersection_with_original"], n)
    D.require(declared.keys() == exact_keys and all(
        D._packed(declared[v["i"], v["j"]]["enclosure"]) == D._packed(v["enclosure"]) for v in H),
        "reported scaled Hessian intersection differs from signed endpoints")
    return f, H


def _scaled_packet(old, new, proof, proof_pins, old_sha, ids):
    same = ("record_id", "dose_exact", "parameter_cell_sha256", "model_observation_contract_sha256",
            "point_cache_sha256", "predictor_sha256", "precision_bits", "element", "chart", "global_indices",
            "free_local", "active", "nodal_hulls_sha256", "auxiliary_element_rule_sha256")
    D.require(all(new.get(key) == old[key] == proof[key] for key in same)
              and new["status"] == "SCALED_OFFAXIS_COMPONENT_UNREVIEWED"
              and new["chart"] == "off_axis" and proof["global_indices"] == ids
              and all(new.get(k) == v for k, v in proof_pins.items() if k != "prelaunch_sha256")
              and new["original_packet_sha256"] == old_sha and new["selector_proof_sha256"] == D._digest(proof)
              and new["degree"] == old["degree"] == proof["unchanged_degree"]
              and new["rho_exact"] == old["rho_exact"] == proof["unchanged_rho_exact"],
              "scaled packet provenance, record, free projection or original domain mismatch")
    force_tail = [_nonnegative_tail(v, proof["force_tail_budget_exact"]) for v in new["force_two_mass_tail_upper_exact"]]
    H_tail = [_nonnegative_tail(v["tail_upper_exact"], proof["H_entry_tail_budget_exact"])
              for v in new["H_integral_minus_cache_upper"]]
    D.require(len(force_tail) == len(ids) and len(H_tail) == len(ids) * (len(ids) + 1) // 2
              and max(force_tail) == Q(proof["force_tail_max_exact"]) and max(H_tail) == Q(proof["H_entry_tail_max_exact"]),
              "full exact entrywise tails must match unchanged frozen local budgets")
    cache = new["same_Y_cache_consistency"]
    D.require(cache["all64_convex_strain_lines_positive"] is True
              and D._packed(cache["minimum_cache_convex_line_J_lower"])[0] > Q(1, 10**6)
              and D._packed(cache["max_coefficientwise_delta_q_upper"])[0] >= 0,
              "same original 64-site analytic cache domains required")
    return _entry_intersections(old, new, ids)


def _nonnegative_tail(value, budget):
    exact = D._fraction(value)
    D.require(0 <= exact <= D._fraction(budget), "nonnegative exact local tail must meet unchanged budget")
    return exact


def _reduce_admitted(admission, hashes, scale, selectors, old_admission, old_review, expected):
    frozen = old_admission["source_and_manifest_hashes"]
    protocol, domains, body = [
        D._bound(f"{ORIGINAL}/{name}", frozen[f"{ORIGINAL}/{name}"], hashes)
        for name in ("protocol.json", "domains.json", "body_geometry.json")
    ]
    context = D.C.load_context(RECORD)
    box = D._bound(protocol["coordinate_box_file"], protocol["coordinate_box_sha256"], hashes)
    exact_box = D.C.exact_node_box(context, context.old_radius)
    D.require(protocol["row"] == context.row and context.row["parameter_cell_sha256"] == scale["parameter_cell_sha256"]
              and all(D._digest(box[k]) == D._digest(v) for k, v in exact_box.items())
              and domains["body_sha256"] == frozen[f"{ORIGINAL}/body_geometry.json"]
              and body["status"] == "CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING"
              and protocol["precision_bits"] == scale["precision_bits"] == 640,
              "same independently reviewed original box, body, predictor, precision and parameter cell required")
    eta0, z0, radius, bound_provider = D._old_premises(context, protocol, hashes)
    scaled_files = {p["element"]: p for p in admission["expected_scaled_packets"]}
    proofs = {p["element"]: p for p in selectors}
    f, H, count, hcount = [D.ZERO] * 611, {}, [0] * 611, [0] * 611
    replaced, kept_axis, kept_surface = [], [], []
    for part, cases in (("bulk", domains["bulk_selectors"]), ("surface", domains["surface_cases"])):
        lookup = {p["file"]: p["sha256"] for p in old_admission["expected_packets"][part]}
        for i, case in enumerate(cases):
            face, edge = divmod(i, 12)
            name = f"element_{i:03d}.json" if part == "bulk" else f"face{face}_edge{edge:02d}.json"
            path = f"{ORIGINAL}/{part}/{name}"
            old = D._bound(path, lookup[path], hashes)
            ids = D._packet_meta(context, part, i, case, old, body, protocol)
            gp, Hp = old["force_integral_minus_cache"], old["H_integral_minus_cache_upper"]
            if part == "bulk" and i in scaled_files:
                value = scaled_files[i]
                new = D._bound(value["file"], value["sha256"], hashes)
                gp, Hp = _scaled_packet(old, new, proofs[i], expected, lookup[path], ids)
                replaced.append(i)
            elif part == "bulk":
                kept_axis.append(i)
            else:
                kept_surface.append((face, edge))
            D._scatter(f, H, ids, gp, Hp)
            for j in ids:
                count[j] += 1
                hcount[j] += len(ids)
    D.require(replaced == OFFAXIS and kept_axis == AXIS and len(kept_surface) == 24
              and len(f) == 611 and len(H) == 6676 and min(count) > 0
              and count == protocol["force_incidence"] and hcount == protocol["H_row_incidence"],
              "full exact 132+12+24 signed assembly and incidence required")
    H_rows = D._h_rows(611, H)
    rf, rows, norm = D._dyadic_reduce(context.R, f, H_rows)
    D.require(norm == Q(protocol["R_infinity_norm_exact"]) == Q(old_review["R_infinity_norm_exact"]),
              "unmodified exact binary preconditioner norm required")
    RF, b = max(map(D._abs_upper, rf)), max(rows)
    D.require(RF <= Q(old_review["RF_uniform_upper_exact"]) and b <= Q(old_review["b_upper_exact"]),
              "exact intersections cannot widen these endpoint norm bounds")
    eta, z, lhs, ok, error = D._strict_test(eta0, z0, RF, b, radius)
    e = D._minimum(context, hashes, z0, z, radius, norm, bound_provider) if ok else None
    return {
        "status": "CONDITIONAL_UNREVIEWED_SCALED_ORIGINAL_CUBE_SUFFICIENT_TEST_PASSES" if ok
                  else "UNREVIEWED_SCALED_ORIGINAL_CUBE_SUFFICIENT_TEST_UNRESOLVED",
        "record_id": RECORD, "scale_correction_elements": replaced, "original_axis_elements": kept_axis,
        "original_surface_packet_keys": kept_surface, "component_count": 168, "force_entries": len(f),
        "upper_Hessian_entries": len(H), "F_endpoint_sha256": D._digest(f),
        "H_endpoint_sha256": D._digest([(i, j, H[i, j]) for i, j in sorted(H)]),
        "RF_endpoint_sha256": D._digest(rf), "abs_R_H_row_sha256": D._digest(rows),
        "R_infinity_norm_exact": norm, "original_R_force_upper_exact": old_review["RF_uniform_upper_exact"],
        "original_b_upper_exact": old_review["b_upper_exact"],
        "RF_uniform_upper_exact": RF, "b_upper_exact": b, "cached_eta_upper_exact": eta0, "cached_Z_upper_exact": z0,
        "eta_upper_exact": eta, "z_upper_exact": z, "selfmap_lhs_exact": lhs, "normal_radius_exact": radius,
        "sufficient_conditions_pass": ok, "conditional_root_error_upper_exact": error, "conditional_curvature": e,
        "same_point_cell_and_cube": True, "root_assertion_independently_reviewed": False,
        "arithmetic": "Exact signed interval intersection and local entry identity checks here. Original independently "
                      "reviewed endpoint scatter, |H| row sums, integer-dyadic R reductions, and old cache proof reused.",
    }


def reduce(admission_sha):
    hashes, scale, selectors, old_admission, old_review, expected = _bound_preparation()
    admission = D._bound(_relative(OUT / "input_admission.json"), admission_sha, hashes)
    D.require(admission["status"] == "BOUND_COMPLETE_SCALED_MANIFEST_BEFORE_LOCAL_READ"
              and admission["record_id"] == RECORD and admission["local_scaled_payloads_opened"] is False,
              "complete digest-bound admission required before local scaled reads")
    _matches_prelaunch(admission["review_sha256"], hashes)
    for n, h in admission["inputs_and_manifests"].items():
        D._bound(n, h, hashes, parse=False)
    terminal = D._bound(f"{STAGE}/actual/terminal.json", admission["scaled_terminal_sha256"], hashes)
    D._terminal(terminal, admission["scaled_manifest_sha256"])
    manifest = D._bound(f"{STAGE}/actual/results.json", admission["scaled_manifest_sha256"], hashes)
    fields = {**expected, "record_id": RECORD, "status": "ALL_132_SCALED_OFFAXIS_PACKETS_UNREVIEWED",
              "unchanged_axis_elements": AXIS, "original_surface_packets_retained": 24, "global_reduction_run": False}
    D.require(_validate_manifest(manifest, OFFAXIS, scale["actual_expected_files"], fields)
              == admission["expected_scaled_packets"], "changed admission or terminal manifest packet expectations")
    stamp = lambda: datetime.now(timezone.utc).isoformat()
    started = monotonic()
    D._write(OUT / "started.json", {"pid": os.getpid(), "utc": stamp(), "admission_sha256": admission_sha})
    try:
        result = _reduce_admitted(admission, hashes, scale, selectors, old_admission, old_review, expected)
        result.update({"admission_sha256": admission_sha, "source_inputs_and_packets": hashes,
                       "seconds": monotonic() - started})
        digest = D._write(OUT / "results.json", result)
        D._write(OUT / "terminal.json", {"status": "FINISHED", "exit_code": 0, "pid": os.getpid(), "utc": stamp(),
                                        "results_sha256": digest})
        print(json.dumps({"status": result["status"], "results_sha256": digest}), flush=True)
    except BaseException as error:
        D._write(OUT / "failure.json", {"message": str(error), "admission_sha256": admission_sha})
        D._write(OUT / "terminal.json", {"status": "FAILED", "exit_code": 1, "pid": os.getpid(), "utc": stamp()})
        raise


def tests():
    checks = {}

    def check(name, condition):
        D.require(condition, name)
        checks[name] = True

    def reject(name, action):
        try:
            action()
        except ValueError:
            checks[name] = True
        else:
            raise AssertionError(name)

    check("signed_negative_overlap", intersection(_packed(-4, -1), _packed(-3, 2)) == _packed(-3, -1))
    check("signed_zero_crossing_overlap", intersection(_packed(-4, 3), _packed(-1, 2)) == _packed(-1, 2))
    check("closed_contact_is_retained", intersection(_packed(-1, 0), _packed(0, 1)) == _packed(0, 0))
    reject("empty_overlap_rejected", lambda: intersection(_packed(-2, -1), _packed(0, 1)))
    a, b = _packed(-1, 2), _packed(-2, 3)
    hp = lambda v: [{"i": 0, "j": 0, "enclosure": v}]
    old = {"force_integral_minus_cache": [a], "H_integral_minus_cache_upper": hp(a)}
    new = {"force_integral_minus_cache": [b], "H_integral_minus_cache_upper": hp(b),
           "force_intersection_with_original": [a], "H_intersection_with_original": hp(a)}
    check("full_1x1_replaces_by_entry_identity", _entry_intersections(old, new, [0]) == ([a], hp(a)))
    reject("Hessian_partial_coverage_rejected", lambda: _entry_intersections(old, {**new, "H_integral_minus_cache_upper": []}, [0]))
    reject("Hessian_duplicate_entry_rejected",
           lambda: _entry_intersections(old, {**new, "H_integral_minus_cache_upper": hp(b) * 2}, [0]))
    reject("reported_intersection_not_equal_exact_rejected",
           lambda: _entry_intersections(old, {**new, "force_intersection_with_original": [_packed(-1, 1)]}, [0]))
    manifest = {"state": "ready", "files": [{"element": 1, "file": "part/a.json", "sha256": "a" * 64},
                                            {"element": 2, "file": "part/b.json", "sha256": "b" * 64}]}
    verify = lambda p: _validate_manifest(p, [1, 2], ["part/a.json", "part/b.json"], {"state": "ready"})
    check("complete_manifest_accepted", verify(manifest) == manifest["files"])
    reject("partial_manifest_rejected", lambda: verify({**manifest, "files": manifest["files"][:1]}))
    reject("duplicate_scoped_element_rejected", lambda: verify({**manifest, "files": manifest["files"][:1] * 2}))
    reject("path_substitution_rejected", lambda: verify({**manifest, "files": [
        {**manifest["files"][0], "file": "old/a.json"}, manifest["files"][1]]}))
    reject("wrong_prelaunch_context_rejected", lambda: verify({**manifest, "state": "partial"}))
    reject("failed_terminal_cannot_admit", lambda: D._terminal(
        {"status": "FAILED", "exit_code": 1, "pid": 123, "results_sha256": "b" * 64}, "b" * 64))
    reject("nonmatching_terminal_digest_cannot_admit", lambda: D._terminal(
        {"status": "FINISHED", "exit_code": 0, "pid": 123, "results_sha256": "b" * 64}, "a" * 64))
    check("strict_gate_and_error", D._strict_test(Q(1, 8), Q(1, 4), Q(1, 8), Q(1, 8), Q(1))[-2:]
          == (True, Q(2, 5)))
    check("strict_selfmap_boundary_unresolved", D._strict_test(Q(1, 4), Q(1, 4), Q(0), Q(1, 4), Q(1, 2))[-2:]
          == (False, None))
    reject("overflowed_local_tail_rejected", lambda: _nonnegative_tail("1/4", "1/8"))
    return checks


def prepare():
    hashes, _, _, _, _, _ = _bound_preparation()
    checks = tests()
    checks["frozen_original_168_and_scaled_132_preparation_match"] = True
    source = Path(__file__).resolve()
    snapshot = OUT / "reduce_scaled_L_frozen.py"
    protocol = OUT / "protocol.json"
    preflight = OUT / "preflight.json"
    D.require(not any(p.exists() for p in (snapshot, protocol, preflight)), "preserve write-once reduction preparation")
    with snapshot.open("xb") as fp:
        fp.write(source.read_bytes())
    p = {
        "status": "FROZEN_SCALED_L_REDUCTION_PREPARATION_REVIEW_PENDING", "record_id": RECORD,
        "source_file": _relative(source), "source_sha256": _sha(source), "input_hashes": hashes,
        "offaxis_intersections_required": OFFAXIS, "original_axis_packets_retained": AXIS,
        "original_surface_packets_retained": 24, "expected_force_entries": 611, "expected_upper_Hessian_entries": 6676,
        "admission": "External SHA-256 for the true terminal exit-zero receipt and matching full 132-packet results "
                     "manifest, plus bound independent reduction prelaunch review. Freeze expected SHA for every packet "
                     "before opening any of the scaled payloads.",
        "math": "Exact signed entrywise intersections before scattering. Reuse independently reviewed endpoint sums, "
                "integer dyadic R*F and |R| times Hessian absolute-row upper bounds, old cached eta/Z, same original radius, "
                "strict inequality, and conditional error/curvature.",
        "actual_scaled_packet_payloads_opened": False, "actual_scaled_reduction_run": False,
    }
    p_sha = D._write(protocol, p)
    receipt = {"status": "PASS_SCALED_REDUCTION_PREPARATION_ONLY", "checks": checks, "source_sha256": p["source_sha256"],
               "protocol_sha256": p_sha, "check_count": len(checks), "actual_scaled_reduction_run": False}
    pre_sha = D._write(preflight, receipt)
    print(json.dumps({"status": receipt["status"], "source_sha256": p["source_sha256"], "protocol_sha256": p_sha,
                      "preflight_sha256": pre_sha, "checks": len(checks)}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_subparsers(dest="action", required=True)
    action.add_parser("prepare")
    admit_args = action.add_parser("admit")
    for key in ("review", "terminal", "manifest"):
        admit_args.add_argument(f"--{key}-sha256", required=True)
    reduce_args = action.add_parser("reduce")
    reduce_args.add_argument("--admission-sha256", required=True)
    args = parser.parse_args()
    if args.action == "prepare":
        prepare()
    elif args.action == "admit":
        admit(args.review_sha256, args.terminal_sha256, args.manifest_sha256)
    else:
        reduce(args.admission_sha256)
