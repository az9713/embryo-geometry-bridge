"""Prepare or run one bounded record-L element-21 force diagnostic."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
from time import monotonic
from types import SimpleNamespace
from unittest.mock import patch

from flint import arb, arb_mat, ctx as precision


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import exact_cohort_context as C
import exact_cohort_integrator as E
import review_cohort_exact_record as D


RECORD = "optoGEF2_small_col_L"
ELEMENT = 21
FOLDER = f"{D.BASE}/records/{RECORD}/original_cube_v1"
REDUCTION = D.BASE + "/reduction_preparation/"
INPUTS = {
    FOLDER + "/protocol.json": "66bfdb7dbe49d3ece787e7ecc9f16e004f7611a8b946be0e078f8a6806a4c937",
    FOLDER + "/domains.json": "7568831341a79caec28cb12c8a5921f40845a982cef4971b5c064c07bb8d1ac4",
    FOLDER + "/coordinate_box.json": "161ba86e1004626417b627deb213c72bc1f5f75b6e587fe0a89b7dcbda152761",
    FOLDER + "/body_geometry.json": "a26897e1f45513f6231c4ac81b58c8550b677a13b7398da1ac1003a029511f8e",
    FOLDER + "/bulk/element_021.json": "8d4fa6f4a8c93cfaffb089bd5ab22efa9ebce70221d262324de5077b1545654d",
    REDUCTION + RECORD + "_input_admission.json": "13fcd8e1a469d20a49ec050d985c7ffed8af4989af6229d7e45d9b15f13b9735",
    REDUCTION + RECORD + "_reduction_review.json": "5a55e0d139edb39206ab3999681d9858d6b0162f73c7e5103d7a5bed3ebc4514",
    REDUCTION + RECORD + "_signed_endpoint_diagnostics.json": "d251d4b797877e2a5d588b33f33483a7fe6bc8f3af2c1b27d4a7358590ebc619",
    "scripts/review_cohort_exact_record.py": "67db748a5d40f31471ddd3d6c84fae280ff7ead650e6eb4e950a4766e54e0285",
}
U, M, T, JET, core = E.U, E.M, E.T, E.JET, E.core
A = U.A


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _asset(relative, bindings, parse=True):
    return D._bound(relative, INPUTS[relative], bindings, parse=parse)


def _geometry():
    precision.prec = 640
    core.A.cache_clear()
    hashes = {}
    for path in INPUTS:
        _asset(path, hashes, parse=False)
    p = _asset(FOLDER + "/protocol.json", hashes)
    for name, value in p["source_and_input_hashes"].items():
        D._bound(name, value, hashes, parse=False)
    context = C.load_context(RECORD)
    for name, value in context.input_hashes.items():
        D._bound(name, value, hashes, parse=False)
    admission = _asset(REDUCTION + RECORD + "_input_admission.json", hashes)
    reduction = _asset(REDUCTION + RECORD + "_reduction_review.json", hashes)
    diagnostic = _asset(REDUCTION + RECORD + "_signed_endpoint_diagnostics.json", hashes)
    packet_name = FOLDER + "/bulk/element_021.json"
    if (p["row"] != context.row or reduction["record_id"] != RECORD
            or reduction["status"] != "PASS_REDUCTION_UNRESOLVED_SUFFICIENT_ROOT_TEST"
            or reduction["component_count"] != 168
            or reduction["input_admission_sha256"] != INPUTS[REDUCTION + RECORD + "_input_admission.json"]
            or reduction["source_input_and_component_hashes"][packet_name] != INPUTS[packet_name]
            or diagnostic["top_five_local_width_contributions"][0]["file"] != packet_name
            or {item["file"]: item["sha256"] for item in admission["expected_packets"]["bulk"]}[packet_name] != INPUTS[packet_name]):
        raise ValueError("mismatched frozen reduction or diagnostic context")
    domains = _asset(FOLDER + "/domains.json", hashes)
    box = _asset(FOLDER + "/coordinate_box.json", hashes)
    body = _asset(FOLDER + "/body_geometry.json", hashes)
    old = _asset(packet_name, hashes)
    case = domains["bulk_selectors"][ELEMENT]
    if (domains["body_sha256"] != INPUTS[FOLDER + "/body_geometry.json"]
            or domains["coordinate_box_sha256"] != INPUTS[FOLDER + "/coordinate_box.json"]
            or domains["protocol_sha256"] != INPUTS[FOLDER + "/protocol.json"]
            or case["J_bounds_exact"] != [body["triangles"][ELEMENT]["J_bound"][k]
                                          for k in ("lower_exact", "upper_exact")]
            or old["status"] != "UNIFORM_COHORT_BULK_COMPONENT_UNREVIEWED"
            or old["degree"] != case["degree"] or old["rho_exact"] != case["rho_exact"]):
        raise ValueError("unchanged whole-box analytic domain and order required")
    nodes = E.SUR.common_nodes(context.target, box)
    center = context.target.nodes()
    if D._digest(center) != D._digest(box["current_centers_exact"]):
        raise ValueError("explicit original unshifted predictor base required")
    rule = context.target.auxiliary_rule()
    setup = E._bulk_setup(context, nodes, rule, ELEMENT)
    E._validate_selector(case, setup)
    point_fields, ids = E._bulk_fields(context.target, center, ELEMENT)
    point_map = point_fields["global_indices"]
    if (point_fields["chart"]["kind"] != "off_axis" or point_map != case["global_indices"]
            or point_map != old["global_indices"] or context.base["common_gain"] <= 0
            or not context.cell["common_gain"][0] <= context.base["common_gain"] <= context.cell["common_gain"][1]):
        raise ValueError("same off-axis chart, projection and in-cell exact gain required")
    differences = []
    for node in ids:
        for axis in (0, 1):
            center_exact = center[node][axis]
            radius = Q(box["coordinate_radii_exact"][node][axis])
            if not D._contains(E.SUR.rational_hull(nodes[node][axis]), (center_exact - radius, center_exact + radius)):
                raise ValueError("outward ball must contain every exact declared local node value")
            differences.append(A(-radius).union(A(radius)))
    hessians = D._entry_map(old["H_integral_minus_cache_upper"], len(point_map))
    if len(hessians) != len(point_map) * (len(point_map) + 1) // 2:
        raise ValueError("complete frozen local interval Hessian required")
    return hashes, context, case, old, rule, ids, nodes, center, setup, differences, hessians


def _series_first(poly, order):
    coefficients, _ = M.volumetric_derivative_coefficients(order, 1, 5000)
    power, out = T.ONE, {}
    poly = {k: A(v) for k, v in poly.items()}
    for n, coefficient in enumerate(coefficients):
        out = T.add(out, T.scale(power, A(coefficient)))
        if n != order:
            power = T.mul(power, poly)
    return out


def _force_cache_site(star, cache, nodes, select, drive):
    """Signed auxiliary-minus-cache force for one site, same full Y on both maps."""
    delta = [[s - c for s, c in zip(a, b)] for a, b in zip(star, cache)]
    q = [sum((A(v) * y for v, y in zip(row, nodes)), arb(0)) for row in cache]
    dq = [sum((A(v) * y for v, y in zip(row, nodes)), arb(0)) for row in delta]
    line = [v + arb(0).union(w) for v, w in zip(q, dq)]
    mu, kappa = A(1), A(5000)
    g, _, _ = core.bulk_local(q, mu, kappa, drive)
    _, H, J = core.bulk_local(line, mu, kappa, drive)
    dg = arb_mat(5, 5, [v for row in H for v in row]) * arb_mat(5, 1, dq)
    return [sum((A(delta[k][i]) * g[k] + A(star[k][i]) * dg[k, 0] for k in range(5)), arb(0))
            for i in select], J


def _force_parts(fields, jet, fun, box_y, drive, maps, case):
    """Only the 12 force columns. Use the same signed moment functional."""
    degree = case["degree"]
    first = _series_first(T.sub(jet["value"], T.ONE), degree)
    gradient = jet["gradient"]
    support = {k for v in gradient for k in v}
    dual = JET.moment_dual(first, support, fun.defect)
    volume = [JET.contract_factor(v, dual) for v in gradient]
    diagonal = [A(1), A(1) + drive, A(1), A(1) + drive, A(1)]
    shear = []
    for column in zip(*fields["B"]):
        p = {}
        for diag, q, B in zip(diagonal, fields["q"], column):
            p = T.add(p, T.scale(T.mul(q, B), diag))
        shear.append(fun.apply(p))
    cache = [arb(0) for _ in fields["free_local"]]
    minimum_J = None
    max_weight = Q(0)
    for site in maps:
        if site["global_indices"] != fields["global_indices"] or site["selected_local"] != fields["free_local"]:
            raise ValueError("same auxiliary projection required")
        local, J = _force_cache_site(site["star_full"], site["cache_full"], box_y, fields["free_local"], drive)
        weight = site["auxiliary_normalized_weight_exact"]
        if weight <= 0:
            raise ValueError("positive weight required")
        max_weight = max(max_weight, weight)
        cache = [total + A(weight) * value for total, value in zip(cache, local)]
        minimum_J = J.lower() if minimum_J is None else min(minimum_J, J.lower())
    norms = [E._norm(v, fields) for v in gradient]
    mass = Q(case["geometry_two_mass_upper_exact"])
    first_tail = Q(case["phi_prime_tail_upper_exact"])
    tails = [mass * first_tail * v for v in norms]
    if max(tails) > Q(case["force_tail_budget_exact"]):
        raise ValueError("same scalar force tail budget required")
    total = [a + b + c + arb(0, core.upper(A(d))) for a, b, c, d in zip(volume, shear, cache, tails)]
    return {"volumetric_truncation": volume, "quadratic_strain_defect": shear,
            "finite_same_Y_cache_correction": cache, "normalized_force_total": total,
            "two_mass_tail_upper_exact": tails, "cache_min_J_bound": core.pack(minimum_J),
            "largest_cached_weight_exact": max_weight}


def _gain_derivative(context, fields, fun, maps, whole_y):
    """Exact affine axial stiffness difference applied to whole-box Y."""
    phi = T.material_basis(fields["corners"])["phi"]
    z_grad = [T.diff(p, 1) for p in phi]
    outputs = []
    for col in fields["free_local"]:
        val = arb(0)
        axis = col % 2
        k = 1 if axis == 0 else 3
        for j in range(axis, 12, 2):
            integrand = T.mul(z_grad[col // 2], z_grad[j // 2])
            finite = sum((m["auxiliary_normalized_weight_exact"] *
                          (m["star_full"][k][col] * m["star_full"][k][j]
                           - m["cache_full"][k][col] * m["cache_full"][k][j]) for m in maps), Q(0))
            val += (fun.apply(integrand) + A(finite)) * whole_y[j]
        outputs.append(A(context.dose * fields["active"]) * val)
    return outputs


def _state_increment(displacements, hessian):
    n = len(displacements)
    return [sum((A(lo).union(A(hi)) * d for j, d in enumerate(displacements)
                 for lo, hi in (hessian[min(i, j), max(i, j)],)), arb(0)) for i in range(n)]


def _bounds(value):
    return E.SUR.rational_hull(value)


def _intersect(a, b):
    x, y = _bounds(a), _bounds(b)
    lo, hi = max(x[0], y[0]), min(x[1], y[1])
    if lo > hi:
        raise ValueError("independent valid local enclosure paths must overlap")
    return lo, hi


def _pack_profile(profile):
    return {key: [core.pack(v) for v in values] if key in (
        "volumetric_truncation", "quadratic_strain_defect", "finite_same_Y_cache_correction", "normalized_force_total")
        else values for key, values in profile.items()}


def tests():
    checks = {}
    eye = [[Q(int(i == j)) for j in range(5)] for i in range(5)]
    q = list(map(A, (1, Q(1, 8), 0, 1, 1)))
    result, J = _force_cache_site(eye, eye, q, list(range(5)), A(2))
    checks["single_synthetic_matched_cache_site_exact_zero"] = all(v.is_zero() for v in result) and J > 0
    changed = [[Q(v) for v in row] for row in eye]
    changed[1][1] += Q(1, 32)
    delta, _ = _force_cache_site(changed, eye, q, list(range(5)), A(2))
    star_g, _, _ = core.bulk_local([q[0], q[1] * A(Q(33, 32)), *q[2:]], A(1), A(5000), A(2))
    cache_g, _, _ = core.bulk_local(q, A(1), A(5000), A(2))
    exact = [a - b for a, b in zip(star_g, cache_g)]
    exact[1] = A(Q(33, 32)) * star_g[1] - cache_g[1]
    checks["single_synthetic_changed_cache_site_contains_direct_pullback"] = all(x.overlaps(y) for x, y in zip(delta, exact))
    e = {(1, 0): Q(1, 10)}
    a = _series_first(e, 3)
    b, _ = U.series(e, 3)
    checks["force_only_series_matches_reviewed_first_series"] = a.keys() == b.keys() and all(_bounds(a[k]) == _bounds(b[k]) for k in a)
    corners = ((Q(1), Q(0)), (Q(2), Q(0)), (Q(2), Q(1)))
    evaluation_point = (Q(7, 4), Q(0))
    # A synthetic signed moment test with geometry terms only.
    functional = SimpleNamespace(apply=lambda p: A(T.evaluate(p, *evaluation_point)))
    context = SimpleNamespace(dose=Q(3, 2))
    field = {"corners": corners, "free_local": [2, 3], "active": Q(1)}
    Y = [A(0)] * 12
    Y[2] = A(2)
    derivative = _gain_derivative(context, field, functional, [], Y)
    p = T.material_basis(corners)["phi"][1]
    expected = 3 * T.evaluate(T.diff(p, 1), *evaluation_point)**2
    checks["affine_gain_derivative_dose_and_component_parity"] = (
        expected > 0 and derivative[0].contains(A(expected)) and derivative[1].is_zero())
    inc = _state_increment([A(-2), A(1)], {(0, 0): (Q(1), Q(1)), (0, 1): (Q(-2), Q(-2)), (1, 1): (Q(3), Q(3))})
    checks["signed_upper_triangle_state_term"] = [_bounds(v) for v in inc] == [(Q(-4), Q(-4)), (Q(7), Q(7))]
    checks["exact_intersection_accepts_signed_overlap"] = _intersect(A(-2).union(A(0)), A(-1).union(A(1)))[0] <= -1
    try:
        _intersect(A(-1), A(1))
        checks["disjoint_enclosures_rejected"] = False
    except ValueError:
        checks["disjoint_enclosures_rejected"] = True
    if not D._yes(checks):
        raise AssertionError(checks)
    return checks


def prepare():
    forbidden = AssertionError("preparation attempted actual moments or a cached quadrature site")
    with patch.object(E.P.Functional, "__init__", side_effect=forbidden), \
            patch.object(U.C, "exact_bulk_maps", side_effect=forbidden):
        hashes, context, case, old, _, ids, _, _, _, _, _ = _geometry()
        checks = tests()
    width = [hi - lo for lo, hi in map(D._packed, old["force_integral_minus_cache"])]
    tail = [D._packed(v)[1] for v in old["force_two_mass_tail_upper"]]
    checks["dominant_existing_local_width_many_orders_above_certified_tail"] = max(width) > 10**10 * max(tail)
    checks["same_frozen_chart_degree_and_12_local_forces"] = case["degree"] == 12 and len(width) == 12 and len(ids) == 6
    if not D._yes(checks):
        raise AssertionError(checks)
    source = Path(__file__).resolve()
    source_sha = _sha(source)
    protocol = {
        "status": "FROZEN_LOCAL_FORCE_DIAGNOSTIC_REVIEW_PENDING",
        "source_file": source.relative_to(ROOT).as_posix(), "source_sha256": source_sha,
        "input_hashes": hashes, "record_id": RECORD, "element": ELEMENT, "precision_bits": 640,
        "element_nodes": ids, "chart": case["chart"], "order": case["degree"], "rho_exact": case["rho_exact"],
        "same_C_sha256": context.row["parameter_cell_sha256"], "same_dose_exact": context.dose,
        "origin": "Frozen cache proof center Y0=exact reference+u, at exact base gain times dose.",
        "whole_box": "The original exact interval nodal hulls times the same closed C.",
        "direct_diagnostic": "Return force-only signed contributions for phi-prime truncation, quadratic "
                             "strains, finite cache correction and scalar tail on point and whole box.",
        "centered_enclosure": "F defect at (Y0,g0) plus whole-box Hessian defect times (Y-Y0), "
                              "plus (g-g0) times exact affine gain derivative at Y. "
                              "Gain derivative uses signed axial-basis stiffness differences on identical geometry/rule.",
        "run_budget": "One off-axis element, two force-only evaluations with shared signed moments at order 12, "
                      "plus small affine gain calculation and frozen local H multiplication. Exactly 64 sites "
                      "per force profile, with no fresh Hessian moments or global assembly.",
        "decision": "Compare component widths at fixed precision and scope. An unresolved centered bound motivates "
                    "rebuilding exact basis functions in local scaled material coordinates before interval insertion, "
                    "while preserving Laurent denominator, material measure and auxiliary rule.",
        "existing_max_local_width_exact": max(width), "existing_max_local_force_tail_exact": max(tail),
        "actual_record_force_integration_run": False,
    }
    names = ("element21_force_frozen.py", "protocol.json", "preflight.json")
    if any((OUT / name).exists() for name in names):
        raise ValueError("write-once preparation already exists")
    with (OUT / names[0]).open("xb") as stream:
        stream.write(source.read_bytes())
    p_sha = D._write(OUT / names[1], protocol)
    pre = {"status": "PASS_PREPARATION_ONLY", "source_sha256": source_sha, "protocol_sha256": p_sha,
           "checks": checks, "actual_record_force_integration_run": False}
    pre_sha = D._write(OUT / names[2], pre)
    print(json.dumps({"status": pre["status"], "checks": len(checks), "source_sha256": source_sha,
                      "protocol_sha256": p_sha, "preflight_sha256": pre_sha}), flush=True)


def run(expected_review_sha256):
    review_path = OUT / "prelaunch_independent_review.json"
    gate = D._bound(review_path.relative_to(ROOT).as_posix(), expected_review_sha256, {})
    protocol_path = OUT / "protocol.json"
    preflight_path = OUT / "preflight.json"
    expected = {"source_sha256": _sha(__file__), "protocol_sha256": _sha(protocol_path),
                "preflight_sha256": _sha(preflight_path)}
    if gate["status"] != "PASS" or not D._yes(gate["checks"]) or any(gate[k] != v for k, v in expected.items()):
        raise ValueError("matching independently reviewed element-21 preparation required")
    p = json.loads(protocol_path.read_bytes())
    hashes, context, case, old, rule, ids, nodes, center, setup, differences, hessians = _geometry()
    if p["input_hashes"] != hashes or p["source_sha256"] != expected["source_sha256"]:
        raise ValueError("frozen element-21 source or input set changed")
    started = monotonic()
    D._write(OUT / "started.json", {"record_id": RECORD, "element": ELEMENT, **expected,
                                    "review_sha256": expected_review_sha256, "pid": os.getpid(),
                                    "utc": datetime.now(timezone.utc).isoformat()})
    try:
        _, whole, wide_drive, whole_jet, _, _, _ = setup
        local_box = [A(nodes[k][axis]) for k in ids for axis in (0, 1)]
        local_center = [A(center[k][axis]) for k in ids for axis in (0, 1)]
        fun = E.P.Functional(context.target, whole, ELEMENT, rule)
        target = U.RuleTarget(context.target, rule)
        maps = [U.C.exact_bulk_maps(target, ELEMENT, iq) for iq in range(64)]
        point, _ = E._bulk_fields(context.target, center, ELEMENT)
        point_jet = JET.bulk_jacobian_jet(point["q"], point["B"])
        point_profile = _force_parts(point, point_jet, fun, local_center, A(context.drive), maps, case)
        print(json.dumps({"phase": "point_force_complete", "seconds": monotonic() - started}), flush=True)
        direct_profile = _force_parts(whole, whole_jet, fun, local_box, wide_drive, maps, case)
        print(json.dumps({"phase": "whole_box_direct_force_complete", "seconds": monotonic() - started}), flush=True)
        lo, hi = context.cell["common_gain"]
        gain_offset = A(lo - context.base["common_gain"]).union(A(hi - context.base["common_gain"]))
        gain = _gain_derivative(context, whole, fun, maps, local_box)
        gain_effect = [gain_offset * value for value in gain]
        state = _state_increment([differences[j] for j in whole["free_local"]],
                                 {key: D._packed(value["enclosure"]) for key, value in hessians.items()})
        centered = [a + b + c for a, b, c in zip(point_profile["normalized_force_total"], state, gain_effect)]
        direct = direct_profile["normalized_force_total"]
        overlap = [_intersect(a, b) for a, b in zip(direct, centered)]
        for a, b in zip(direct, old["force_integral_minus_cache"]):
            lo, hi = D._packed(b)
            if not D._contains((lo, hi), _bounds(a)):
                # Different associations can widen either enclosure. It is the
                # intersection (not total containment) which the two imply.
                _intersect(a, A(lo).union(A(hi)))
        result = {
            "status": "LOCAL_FORCE_DIAGNOSTIC_COMPLETE_UNREVIEWED",
            "record_id": RECORD, "element": ELEMENT, "precision_bits": 640,
            "global_indices": whole["global_indices"], "domain_hashes": hashes, **expected,
            "review_sha256": expected_review_sha256, "point": _pack_profile(point_profile),
            "whole_box_direct": _pack_profile(direct_profile),
            "whole_box_centered_state_increment": [core.pack(v) for v in state],
            "whole_box_gain_derivative": [core.pack(v) for v in gain],
            "whole_box_gain_increment": [core.pack(v) for v in gain_effect],
            "whole_box_force_centered": [core.pack(v) for v in centered],
            "whole_box_intersection_exact": overlap,
            "moment_defects_needed_for_force_and_affine_gain_only": len(fun.requested_defect_keys),
            "fresh_Hessian_moments_or_assembly": False, "new_root_test": False,
            "seconds": monotonic() - started,
        }
        digest = D._write(OUT / "results.json", result)
        D._write(OUT / "terminal.json", {"status": "FINISHED", "exit_code": 0, "results_sha256": digest,
                                        "pid": os.getpid(), "utc": datetime.now(timezone.utc).isoformat()})
        print(json.dumps({"status": result["status"], "sha256": digest}), flush=True)
    except BaseException as error:
        D._write(OUT / "failure.json", {"error": type(error).__name__, "message": str(error), **expected})
        D._write(OUT / "terminal.json", {"status": "FAILED", "exit_code": 1, "pid": os.getpid(),
                                        "utc": datetime.now(timezone.utc).isoformat()})
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_subparsers(dest="mode", required=True)
    mode.add_parser("prepare")
    run_args = mode.add_parser("run")
    run_args.add_argument("--review-sha256", required=True)
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare()
    else:
        run(args.review_sha256)
