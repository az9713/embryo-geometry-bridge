"""Scaled off-axis component adapter and frozen 132-element L preparation."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
from time import monotonic
from unittest.mock import patch

from flint import arb, ctx as precision


ROOT = Path(__file__).resolve().parents[1]
BASE = "runs/exact_integration_cohort_v1/L_defect_diagnostic"
PILOT = BASE + "/scaled_chart_v1"
STAGE = ROOT / BASE / "scaled_offaxis_L_v1"
PREP = STAGE / "preparation"
RECORD = "optoGEF2_small_col_L"
ORIGINAL = f"runs/exact_integration_cohort_v1/records/{RECORD}/original_cube_v1"
PILOT_HASHES = {
    PILOT + "/scaled_element21.py": "a173b3759fb66fe66e5f3e442c1bcce9a5a185e67515649a94499f82d1d67f60",
    PILOT + "/protocol.json": "76e5086c374fe86f3659ef1d1b099c4a733211245fbb1fa70251b01462c3fc5c",
    PILOT + "/preflight.json": "a35814e39106e14514f2afa6134996781d0e08f5036bcd6e22b1aafd5587a505",
    PILOT + "/actual_independent_review.json": "f87cc324c6469f0717cf399ebd4f695302ca8495698296ea886987b7cf45e9d4",
    PILOT + "/results.json": "f537a5f635d6e1f0d124bb9e52e135ecf3d88cc854559e9a27725e80cfe9683a",
    PILOT + "/terminal.json": "3d99e5a292b47d919e00e815b41d596fcb5ae5dcb7c245704ea98b7813c47ab9",
    ORIGINAL + "/bulk_terminal.json": "68837ae925b7ab406221e5b678fc2c950a3f8d6abc6ab601b57f56c704734b8a",
    ORIGINAL + "/surface_terminal.json": "84b5c6510762fc37e72c7a94133b251d668ca49d4014537ce2d10bf1bdc0b571",
    ORIGINAL + "/assembled.json": "13e06c3f23c0e840a9b3043f8b09d40757cff5b84831f728c4b413ae9a8bc438",
}
sys.path.insert(0, str(ROOT / PILOT))
import scaled_element21 as S


D, E, C = S.D, S.E, S.F.C
U, M, T, JET, core, A = S.U, S.M, S.T, S.JET, S.core, S.A


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _part(path):
    return Path(path).relative_to(ROOT).as_posix()


def _endpoint_intersection(left, right):
    a, b = D._packed(left), D._packed(right)
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    D.require(lo <= hi, "disjoint scaled and original chart enclosures")
    return {"lower_exact": str(lo), "upper_exact": str(hi)}


def _geometry(context, nodes, rule, body_element, selector, element):
    D.require(type(element) is int and 0 <= element < len(context.target.a["triangles"]),
              "explicit valid integer bulk element required")
    D.require(selector.get("element") == element and selector.get("chart") == "off_axis",
              "same immutable off-axis selector required")
    setup = E._bulk_setup(context, nodes, rule, element)
    E._validate_selector(selector, setup)
    identity, fields, drive, _, mass, _, _ = setup
    lo, hi, raw, rounded = E._body_rho(body_element)
    D.require(fields["chart"]["kind"] == "off_axis"
              and selector["J_bounds_exact"] == [str(lo), str(hi)]
              and Q(selector["raw_rho_exact"]) == raw and Q(selector["rho_exact"]) == rounded,
              "same independently reviewed element geometry and convergent rho required")
    ids = tuple(map(int, context.target.a["triangles"][element]))
    geometry = S.chart(fields["corners"])
    D.require(geometry.Rhat_min > 0 and geometry.det_physical > 0, "strict positive off-axis material chart required")
    return identity, fields, drive, mass, ids, geometry


def _tail_bounds(case, mass, pointwise_gradient_norms, pointwise_H_norms):
    # Either coordinate representation bounds the same physical derivative.
    # Keep the original certified maxima available as per-entry fallbacks.
    g0, h0 = (Q(case[k]) for k in ("max_DJ_norm_upper_exact", "max_D2J_norm_upper_exact"))
    first, second = (Q(case[k]) for k in ("phi_prime_tail_upper_exact", "phi_second_tail_upper_exact"))
    g = [min(g0, norm) for norm in pointwise_gradient_norms]
    h = {key: min(h0, norm) for key, norm in pointwise_H_norms.items()}
    f = [mass * first * norm for norm in g]
    H = {(i, j): mass * (second * g[i] * g[j] + first * norm) for (i, j), norm in h.items()}
    D.require(max(f) <= Q(case["force_tail_budget_exact"]) and max(H.values()) <= Q(case["hessian_tail_budget_exact"]),
              "retain original local force and Hessian tail budgets")
    return f, H


def prove_selector(context, nodes, rule, body_element, selector, element, center=None):
    """Pure exact identities and interval chart derivative majorants; no integration."""
    identity, fields, _, mass, ids, geometry = _geometry(context, nodes, rule, body_element, selector, element)
    center = context.target.nodes() if center is None else center
    free = fields["free_local"]
    rational_checks = S.rational_identities(geometry, [center[j] for j in ids], free)
    fixed = set(range(12)) - set(free)
    for index in fixed:
        i, axis = divmod(index, 2)
        c = center[ids[i]][axis]
        D.require(E.SUR.rational_hull(nodes[ids[i]][axis]) == (c, c), "fixed local coordinate lost its exact value")
    # Equality of q at the exact origin and every free linear column B
    # establishes equality for every supplied fixed-coordinate nodal field.
    locations = [S.chart_point(geometry, v) for v in rule["volume_barycentric"]]
    D.require(len(locations) == 64 and all(abs(x) <= 1 and abs(y) <= 1 and r > 0 for x, y, r in locations),
              "all original exact auxiliary sites must remain in the scaled rectangle")
    q, B = S.numerators(geometry, [nodes[j] for j in ids], free)
    jet = JET.bulk_jacobian_jet(q, B)
    g = [S._norm(p, geometry) for p in jet["gradient"]]
    H = {key: S._norm(p, geometry) for key, p in jet["hessian_upper"].items()}
    ft, ht = _tail_bounds(selector, mass, g, H)
    return {
        **identity, "physical_selector_sha256": D._digest(selector),
        "body_element_sha256": D._digest(body_element), "geometry": S._metadata(geometry),
        "local_nodes": ids, "fixed_local_positions": sorted(fixed), "polynomial_identity_count": rational_checks,
        "invariant_Qstar_node_count": len(locations), "scaled_Qstar_coordinates_sha256": D._digest(locations),
        "unchanged_degree": selector["degree"], "unchanged_rho_exact": selector["rho_exact"],
        "force_tail_budget_exact": selector["force_tail_budget_exact"], "H_entry_tail_budget_exact": selector["hessian_tail_budget_exact"],
        "scaled_max_DJ_majorant_exact": max(g), "scaled_max_D2J_majorant_exact": max(H.values()),
        "force_tail_max_exact": max(ft), "H_entry_tail_max_exact": max(ht.values()),
    }


def scaled_bulk_component(context, nodes, rule, body_element, selector, proof, old_packet):
    """Compute one explicit off-axis component for an admitted frozen chart/box."""
    element = proof["element"]
    fresh = prove_selector(context, nodes, rule, body_element, selector, element)
    D.require(D._digest(fresh) == D._digest(proof), "scaled selector geometry or polynomial identity proof changed")
    identity, fields, drive, mass, ids, geometry = _geometry(context, nodes, rule, body_element, selector, element)
    D.require(old_packet["record_id"] == context.record_id and old_packet["element"] == element
              and old_packet["global_indices"] == fields["global_indices"] and old_packet["chart"] == "off_axis"
              and old_packet["nodal_hulls_sha256"] == identity["nodal_hulls_sha256"]
              and old_packet["auxiliary_element_rule_sha256"] == identity["auxiliary_element_rule_sha256"]
              and old_packet["degree"] == selector["degree"], "original off-axis packet binding mismatch")
    started = monotonic()
    fun = S.Functional(geometry, context.target.energy_scale, rule["volume_barycentric"], rule["volume_weights"][element])
    try:
        analytic_mass = geometry.det_physical * sum(value[0] for value in geometry.corners) / 6
        D.require(fun.integral(0, 0, 0).overlaps(context.target.normalization() * A(analytic_mass))
                  and fun.weight_mass > 0, "same positive original material and auxiliary masses required")
        F, H, g, h = S._eval_parts(geometry, [nodes[j] for j in ids], fields["free_local"], fun, drive, selector["degree"])
        force_tails, H_tails = _tail_bounds(selector, mass, g, h)
        cf, ch, cache = U.uniform_cache(U.RuleTarget(context.target, rule), fields, element, nodes, drive)
        D.require(cache["all64_convex_strain_lines_positive"] is True, "same-node cache lines must stay analytic")

        def enclosure(value, tail):
            return core.pack(value + arb(0, core.upper(A(tail))))

        force = [enclosure(v + cf[i], force_tails[i]) for i, v in enumerate(F)]
        Hscaled = [{"i": i, "j": j, "enclosure": enclosure(value + ch[i, j], H_tails[i, j]),
                    "tail_upper_exact": H_tails[i, j]} for (i, j), value in H.items()]
        old_H = D._entry_map(old_packet["H_integral_minus_cache_upper"], len(force))
        D.require(old_H.keys() == H.keys() and len(old_packet["force_integral_minus_cache"]) == len(force),
                  "matching complete original force/Hessian packet required")
        f_intersection = [_endpoint_intersection(a, b)
                          for a, b in zip(force, old_packet["force_integral_minus_cache"])]
        H_intersection = [{"i": p["i"], "j": p["j"],
                           "enclosure": _endpoint_intersection(p["enclosure"], old_H[p["i"], p["j"]]["enclosure"])}
                          for p in Hscaled]
        return {
            **identity, "status": "SCALED_OFFAXIS_COMPONENT_UNREVIEWED",
            "degree": selector["degree"], "rho_exact": selector["rho_exact"],
            "selector_proof_sha256": D._digest(proof),
            "force_integral_minus_cache": force, "H_integral_minus_cache_upper": Hscaled,
            "force_two_mass_tail_upper_exact": force_tails, "same_Y_cache_consistency": cache,
            "force_intersection_with_original": f_intersection, "H_intersection_with_original": H_intersection,
            "signed_moment_defects": fun.defect.cache_info().misses, "seconds": monotonic() - started,
        }
    finally:
        fun.integral.cache_clear()
        fun.defect.cache_clear()
        fun.site_power.cache_clear()
        M._triangle_moment.cache_clear()
        core.A.cache_clear()


def _pinned_L_inputs():
    precision.prec = 640
    core.A.cache_clear()
    hashes = {}
    pinned = {name: D._bound(name, digest, hashes, parse=name.endswith(".json"))
              for name, digest in PILOT_HASHES.items()}
    checked = pinned[PILOT + "/actual_independent_review.json"]
    D.require(checked["status"] == "PASS_ACTUAL_SCALED_E21_ENDPOINT_REVIEW" and D._yes(checked["checks"])
              and checked["results_sha256"] == PILOT_HASHES[PILOT + "/results.json"]
              and _sha(S.__file__) == PILOT_HASHES[PILOT + "/scaled_element21.py"], "reviewed scaled pilot kernel required")
    D._terminal(pinned[PILOT + "/terminal.json"], PILOT_HASHES[PILOT + "/results.json"])
    D._bound(checked["support_file"], checked["support_sha256"], hashes, parse=False)
    source = S.F.INPUTS
    D._bound(f"{BASE}/element21_force.py", S.PREMISES[f"{BASE}/element21_force.py"], hashes, parse=False)
    for name, digest in source.items():
        D._bound(name, digest, hashes, parse=False)
    protocol = D._bound(ORIGINAL + "/protocol.json", source[ORIGINAL + "/protocol.json"], hashes)
    for name, digest in protocol["source_and_input_hashes"].items():
        D._bound(name, digest, hashes, parse=False)
    domains = D._bound(ORIGINAL + "/domains.json", source[ORIGINAL + "/domains.json"], hashes)
    box = D._bound(protocol["coordinate_box_file"], protocol["coordinate_box_sha256"], hashes)
    body = D._bound(ORIGINAL + "/body_geometry.json", domains["body_sha256"], hashes)
    admitted_path = S.F.REDUCTION + RECORD + "_input_admission.json"
    admitted = D._bound(admitted_path, source[admitted_path], hashes)
    review_path = S.F.REDUCTION + RECORD + "_reduction_review.json"
    reviewed = D._bound(review_path, source[review_path], hashes)
    D.require(reviewed["status"] == "PASS_REDUCTION_UNRESOLVED_SUFFICIENT_ROOT_TEST"
              and reviewed["component_count"] == 168 and reviewed["input_admission_sha256"] == source[admitted_path]
              and admitted["external_pins"]["assembled"] == PILOT_HASHES[ORIGINAL + "/assembled.json"]
              and admitted["external_pins"]["bulk_terminal"] == PILOT_HASHES[ORIGINAL + "/bulk_terminal.json"]
              and admitted["external_pins"]["surface_terminal"] == PILOT_HASHES[ORIGINAL + "/surface_terminal.json"],
              "exact reviewed original 168-packet context required")
    for part in ("bulk", "surface"):
        D._terminal(pinned[f"{ORIGINAL}/{part}_terminal.json"])
    ctx = C.load_context(RECORD)
    D.require(ctx.row == protocol["row"] and domains["coordinate_box_sha256"] == protocol["coordinate_box_sha256"]
              and body["status"] == "CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING"
              and domains["protocol_sha256"] == hashes[ORIGINAL + "/protocol.json"], "same original cube required")
    for name, digest in ctx.input_hashes.items():
        D._bound(name, digest, hashes, parse=False)
    packets = admitted["expected_packets"]["bulk"]
    D._manifest(ORIGINAL + "/bulk", packets, [f"element_{i:03d}.json" for i in range(144)])
    expected = {item["file"]: item["sha256"] for item in packets}
    D.require(all(reviewed["source_input_and_component_hashes"][p] == h for p, h in expected.items()),
              "frozen original bulk packets differ from reviewed reduction")
    nodes = E.SUR.common_nodes(ctx.target, box)
    old_axis = [i for i, tri in enumerate(ctx.target.a["triangles"])
                if any(T.rational(ctx.target.a["reference"][int(j), 0]) == 0 for j in tri[:3])]
    offaxis = [i for i in range(144) if i not in old_axis]
    D.require(len(offaxis) == 132 and old_axis == [*range(6), *range(72, 78)]
              and len(domains["bulk_selectors"]) == len(body["triangles"]) == 144,
              "scope must be all 132 off-axis L elements and original 12 axis elements")
    return hashes, ctx, nodes, ctx.target.auxiliary_rule(), body, domains, expected, offaxis, old_axis


def prepare():
    forbidden = AssertionError("no actual cohort integration or local cache sums during preparation")
    with patch.object(S.Functional, "__init__", side_effect=forbidden), \
            patch.object(U, "uniform_cache", side_effect=forbidden), \
            patch.object(U.C, "exact_bulk_maps", side_effect=forbidden), \
            patch.object(M, "triangle_moment", side_effect=forbidden):
        hashes, ctx, nodes, rule, body, domains, originals, elements, original_axis = _pinned_L_inputs()
        center = ctx.target.nodes()
        proofs = [
            prove_selector(ctx, nodes, rule, body["triangles"][i], domains["bulk_selectors"][i], i, center)
            for i in elements
        ]
        checks = {
            "all_132_explicit_offaxis_selectors": [p["element"] for p in proofs] == elements,
            "all_frozen_Qstar_sites_exact_and_unchanged": sum(p["invariant_Qstar_node_count"] for p in proofs) == 132 * 64,
            "all_free_columns_and_exact_Q_J_DJ_D2J_polynomial_identities": min(p["polynomial_identity_count"] for p in proofs) > 0,
            "all_scaled_material_triangles_strict_positive_R": all(p["geometry"]["Rhat_min"] > 0 for p in proofs),
            "original_per_element_force_and_H_tail_budgets": all(
                p["force_tail_max_exact"] <= Q(p["force_tail_budget_exact"])
                and p["H_entry_tail_max_exact"] <= Q(p["H_entry_tail_budget_exact"]) for p in proofs),
            "original_12_axis_and_24_surface_packets_retained": original_axis == [*range(6), *range(72, 78)],
            "no_actual_132_M_integrals_or_cache_sums_called": True,
        }
        D.require(D._yes(checks), "all off-axis chart proofs required")
    files = ("exact_scaled_offaxis_cohort_frozen.py", "selectors.json", "protocol.json", "preflight.json")
    D.require(not any((PREP / name).exists() for name in files), "write-once off-axis preparation exists")
    PREP.mkdir(parents=True, exist_ok=True)
    src = Path(__file__).resolve()
    with (PREP / files[0]).open("xb") as stream:
        stream.write(src.read_bytes())
    selectors_sha = D._write(PREP / files[1], proofs)
    protocol = {
        "status": "FROZEN_132_OFFAXIS_L_PREPARATION_REVIEW_PENDING", "record_id": RECORD,
        "source_file": _part(src), "source_sha256": _sha(src), "input_hashes": hashes,
        "selectors_file": _part(PREP / files[1]), "selectors_sha256": selectors_sha,
        "precision_bits": 640, "parameter_cell_sha256": ctx.row["parameter_cell_sha256"],
        "original_cube_radius_exact": ctx.old_radius, "dose_exact": ctx.dose,
        "offaxis_element_indices": elements, "original_axis_element_indices": original_axis,
        "original_bulk_packet_hashes": originals,
        "actual_expected_files": [_part(STAGE / "actual" / "elements" / f"element_{i:03d}.json") for i in elements],
        "actual_job_output": _part(STAGE / "actual"),
        "algebra": "Reviewed scaled element-21 primitives generalized by exact off-axis chart geometry and explicit record context. "
                   "Exact q at base and all exact B columns identify the affine strain field for the full unchanged nodal box. "
                   "Exact J, DJ, D2J numerator identities are tested on each element. The normalized Rhat denominator is retained.",
        "tails": "Reuse each original degree, dyadic rho, scalar tails and positive mass. A per-entry minimum of original global "
                 "derivative majorants and new scaled derivative majorants preserves both frozen budgets.",
        "reduction": "After actual terminal exit zero, independently review charts, 132 packet identities and endpoint arithmetic. "
                     "Admit only matching original/validated scaled intersections by element ID; retain original 12 axis and "
                     "24 surface packets. Recompute 611 force and sparse symmetric Hessian sums and strict root/curvature tests "
                     "on the SAME original cube. No replacement or root assertion during preparation.",
        "actual_force_or_Hessian_run": False,
    }
    protocol_sha = D._write(PREP / files[2], protocol)
    report = {"status": "PASS_132_OFFAXIS_PREPARATION_ONLY", "source_sha256": protocol["source_sha256"],
              "protocol_sha256": protocol_sha, "selectors_sha256": selectors_sha, "checks": checks,
              "elements": len(proofs), "exact_polynomial_identities": sum(p["polynomial_identity_count"] for p in proofs),
              "actual_force_or_Hessian_run": False}
    preflight_sha = D._write(PREP / files[3], report)
    print(json.dumps({"status": report["status"], "source_sha256": report["source_sha256"],
                      "protocol_sha256": protocol_sha, "preflight_sha256": preflight_sha,
                      "selectors_sha256": selectors_sha, "elements": len(proofs),
                      "identities": report["exact_polynomial_identities"]}), flush=True)


def run(review_sha256):
    review = D._bound(_part(PREP / "prelaunch_independent_review.json"), review_sha256, {})
    expected = {"source_sha256": _sha(__file__), "protocol_sha256": _sha(PREP / "protocol.json"),
                "preflight_sha256": _sha(PREP / "preflight.json"), "selectors_sha256": _sha(PREP / "selectors.json")}
    D.require(review["status"] == "PASS" and D._yes(review["checks"])
              and all(review[k] == v for k, v in expected.items()), "bound independent scaled off-axis prelaunch review required")
    p = json.loads((PREP / "protocol.json").read_bytes())
    hashes, ctx, nodes, rule, body, domains, originals, elements, original_axis = _pinned_L_inputs()
    proofs = json.loads((PREP / "selectors.json").read_bytes())
    D.require(hashes == p["input_hashes"] and originals == p["original_bulk_packet_hashes"]
              and [k["element"] for k in proofs] == elements and original_axis == p["original_axis_element_indices"]
              and expected["selectors_sha256"] == p["selectors_sha256"], "complete matching frozen scaled proof set required")
    out = STAGE / "actual"
    (out / "elements").mkdir(parents=True, exist_ok=True)
    stamp = lambda: datetime.now(timezone.utc).isoformat()
    started = monotonic()
    D._write(out / "started.json", {**expected, "prelaunch_sha256": review_sha256, "pid": os.getpid(), "utc": stamp()})
    packets = []
    try:
        for i, proof in zip(elements, proofs):
            name = f"element_{i:03d}.json"
            original_file = f"{ORIGINAL}/bulk/{name}"
            old = D._bound(original_file, originals[original_file], {})
            value = scaled_bulk_component(ctx, nodes, rule, body["triangles"][i], domains["bulk_selectors"][i], proof, old)
            path = out / "elements" / name
            digest = D._write(path, {**value, "original_packet_sha256": originals[original_file], **expected})
            packets.append({"element": i, "file": _part(path), "sha256": digest})
            print(json.dumps({"element": i, "complete": len(packets), "total": len(elements), "seconds": monotonic() - started}), flush=True)
        D.require([f["file"] for f in packets] == p["actual_expected_files"], "exact complete scaled packet manifest required")
        digest = D._write(out / "results.json", {"status": "ALL_132_SCALED_OFFAXIS_PACKETS_UNREVIEWED",
                                                "record_id": RECORD, **expected, "prelaunch_sha256": review_sha256,
                                                "files": packets, "unchanged_axis_elements": original_axis,
                                                "original_surface_packets_retained": 24,
                                                "global_reduction_run": False, "seconds": monotonic() - started})
        D._write(out / "terminal.json", {"status": "FINISHED", "exit_code": 0, "pid": os.getpid(),
                                        "results_sha256": digest, "utc": stamp()})
    except BaseException as error:
        D._write(out / "failure.json", {"error": str(error), "completed_packets": packets, **expected})
        D._write(out / "terminal.json", {"status": "FAILED", "exit_code": 1, "pid": os.getpid(), "utc": stamp()})
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest="command", required=True)
    command.add_parser("prepare")
    work = command.add_parser("run")
    work.add_argument("--review-sha256", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare()
    else:
        run(args.review_sha256)
