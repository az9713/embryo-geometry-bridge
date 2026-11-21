"""Prepare and evaluate centered bulk force only on L's reviewed original cube."""
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

from flint import arb, ctx as precision


ROOT = Path(__file__).resolve().parents[1]
BASE = "runs/exact_integration_cohort_v1/L_defect_diagnostic"
SCALED = BASE + "/scaled_offaxis_L_v1"
STAGE = ROOT / BASE / "centered_bulk_L_v1"
PREP = STAGE / "preparation"
RECORD = "optoGEF2_small_col_L"
sys.path.insert(0, str(ROOT / BASE))
sys.path.insert(0, str(ROOT / SCALED / "reduction_preparation"))
import element21_force as F
import reduce_scaled_L as B


D, E, C, U, M, T, JET, core, A = F.D, F.E, F.C, F.U, F.M, F.T, F.JET, F.core, F.A
ORIGINAL = B.ORIGINAL
PIN = {
    f"{BASE}/element21_force.py": "9a5f93bb3d4a39020e28bee165ba606d1d17f63e44a17ee758959aa05109d329",
    f"{BASE}/actual_independent_review.json": "148903b4f29c87d27ce88c8a4c81e3bbbedcf2f55a5ad37efd1ecf0dadd32a4d",
    f"{SCALED}/reduction_preparation/reduce_scaled_L.py": "0ebaf0c1dfd65a3962841be51b714f3540b6a7f07e226dab130591bd4b7f70e9",
    f"{SCALED}/reduction_preparation/input_admission.json": "1b6d28f43a732b64c0a716f007cc9db9dfb898a3e3ac0691b25f23ed460a86de",
    f"{SCALED}/reduction_preparation/results.json": "51d12c1ee61901f9164dbd3ab4e366b4887dc115398446d379835bdeb4df094a",
    f"{SCALED}/reduction_preparation/terminal.json": "89c43fbd476c6fc7453432b642faf8c14fe1c564d3ce219d34ab14ea5aad5b4f",
    f"{SCALED}/reduction_preparation/actual_independent_review.json": "54dee8659601c387712c16732ab6aed90fbf313c1e3b0f77930b3e12b19220e6",
    f"{SCALED}/actual_independent_review.json": "bcd6fd98baf2d611ce224054e51ce5cdff366e41547ec09a5f500e9b4cd74bb2",
}


def _relative(path):
    return Path(path).relative_to(ROOT).as_posix()


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _source_inputs():
    precision.prec = 640
    core.A.cache_clear()
    hashes, scale, selectors, old, _, expected = B._bound_preparation()
    documents = {n: D._bound(n, h, hashes, parse=n.endswith(".json")) for n, h in PIN.items()}
    D.require(_sha(F.__file__) == PIN[f"{BASE}/element21_force.py"]
              and _sha(B.__file__) == PIN[f"{SCALED}/reduction_preparation/reduce_scaled_L.py"],
              "same reviewed force-only kernel and scaled intersection implementation required")
    audit = documents[f"{SCALED}/reduction_preparation/actual_independent_review.json"]
    local_review = documents[f"{SCALED}/actual_independent_review.json"]
    pilot = documents[f"{BASE}/actual_independent_review.json"]
    D._bound(local_review["support_file"], local_review["support_sha256"], hashes, parse=False)
    D._bound(pilot["support_file"], pilot["support_sha256"], hashes, parse=False)
    result_name = f"{SCALED}/reduction_preparation/results.json"
    admission_name = f"{SCALED}/reduction_preparation/input_admission.json"
    terminal_name = f"{SCALED}/reduction_preparation/terminal.json"
    scaled_result, scaled = documents[result_name], documents[admission_name]
    D._terminal(documents[terminal_name], PIN[result_name])
    D.require(audit["status"] == "PASS_ACTUAL_SCALED_ORIGINAL_CUBE_REDUCTION_REVIEW"
              and D._yes(audit["checks"]) and audit["result_sha256"] == PIN[result_name]
              and audit["admission_sha256"] == PIN[admission_name]
              and audit["terminal_sha256"] == PIN[terminal_name]
              and audit["component_replay_sha256"] == local_review["support_sha256"]
              and D._yes(local_review["checks"]) and D._yes(pilot["checks"])
              and scaled_result["root_assertion_independently_reviewed"] is False
              and scaled_result["status"] == "UNREVIEWED_SCALED_ORIGINAL_CUBE_SUFFICIENT_TEST_UNRESOLVED"
              and scaled_result["scale_correction_elements"] == B.OFFAXIS
              and scaled_result["original_axis_elements"] == B.AXIS
              and scaled_result["force_entries"] == 611 and scaled_result["upper_Hessian_entries"] == 6676
              and scaled_result["component_count"] == 168 and scaled_result["same_point_cell_and_cube"] is True,
              "same independently reviewed combined L derivative domain required")
    scaled_files = scaled["expected_scaled_packets"]
    D.require([p["element"] for p in scaled_files] == B.OFFAXIS
              and [p["file"] for p in scaled_files] == scale["actual_expected_files"],
              "complete previously admitted 132 scaled packet manifest required")
    packet_hashes = scaled_result["source_inputs_and_packets"]
    old_files = old["expected_packets"]["bulk"]
    D.require(all(packet_hashes[p["file"]] == p["sha256"] for p in old_files + scaled_files),
              "the same actual packet set must feed the centered enclosure")
    protocol, domains, body = [
        D._bound(f"{ORIGINAL}/{name}", old["source_and_manifest_hashes"][f"{ORIGINAL}/{name}"], hashes)
        for name in ("protocol.json", "domains.json", "body_geometry.json")
    ]
    context = C.load_context(RECORD)
    box = D._bound(protocol["coordinate_box_file"], protocol["coordinate_box_sha256"], hashes)
    computed_box = C.exact_node_box(context, context.old_radius)
    D.require(all(D._digest(box[key]) == D._digest(value) for key, value in computed_box.items())
              and protocol["row"] == context.row and context.row["parameter_cell_sha256"] == scale["parameter_cell_sha256"],
              "same exact original unshifted nodal box, cell and record required")
    for path, digest in context.input_hashes.items():
        D._bound(path, digest, hashes, parse=False)
    center = context.target.nodes()
    D.require(D._digest(center) == D._digest(box["current_centers_exact"]),
              "stored exact Y0 must equal the original box center")
    g0 = context.base["common_gain"]
    low, high = context.cell["common_gain"]
    D.require(0 < low <= g0 <= high and context.drive == g0 * context.dose and context.dose > 0,
              "exact base drive and finite in-cell gain segment required")
    # Use the original nodal-ball constructor already reviewed for this box.
    wide = E.SUR.common_nodes(context.target, box)
    rule = context.target.auxiliary_rule()
    return SimpleNamespace(
        hashes=hashes, context=context, box=box, center=center, wide=wide, rule=rule, protocol=protocol,
        domains=domains, body=body, old_files={p["file"]: p["sha256"] for p in old_files},
        scaled_files={p["element"]: p for p in scaled_files}, proofs={s["element"]: s for s in selectors},
        expected=expected, gain_bounds=(low - g0, high - g0),
    )


def axial_gradient(fields):
    """Exact physical Z derivative of every full nodal basis, in the chosen chart."""
    base = T.material_basis(fields["corners"])
    gradients = [T.diff(p, 1) for p in base["phi"]]
    chart = fields["chart"]
    kind = chart["kind"]
    if kind == "axis_vertex":
        zeta = T.sub(chart["Z"], T.constant(base["zcenter"]))
        return [T.compose(p, chart["R"], zeta) for p in gradients]
    D.require(kind in ("off_axis", "axis_edge_physical"), "explicit supported material chart required")
    return gradients


def affine_gain_derivative(dose, active, fields, functional, maps, full_y):
    """At fixed Y, d_g(F_I-F_Q) equals dose times an axial stiffness difference."""
    free = fields["free_local"]
    if active == 0:
        return [arb(0) for _ in free]
    D.require(active == 1 and dose > 0 and len(full_y) == 12,
              "positive fixed dose, elementconstant activation and full same-Y coordinates required")
    z = axial_gradient(fields)
    output = []
    for i in free:
        parity = i % 2
        k = 1 if parity == 0 else 3
        row = arb(0)
        for j in range(parity, 12, 2):
            difference = sum((m["auxiliary_normalized_weight_exact"] *
                              (m["star_full"][k][i] * m["star_full"][k][j]
                               - m["cache_full"][k][i] * m["cache_full"][k][j]) for m in maps), Q(0))
            row += (functional.apply(T.mul(z[i // 2], z[j // 2])) + A(difference)) * full_y[j]
        output.append(A(dose) * row)
    return output


def centered_force(point, displacements, H, gain, gain_interval):
    state = F._state_increment(displacements, H)
    gain_effect = [gain_interval * derivative for derivative in gain]
    return [a + b + c for a, b, c in zip(point, state, gain_effect)], state, gain_effect


def _polynomials_digest(fields):
    def packed(p):
        return [[i, j, str(v)] for (i, j), v in sorted(p.items())]
    return D._digest({"q": [packed(p) for p in fields["q"]], "B": [[packed(p) for p in row] for row in fields["B"]]})


def _element(data, index):
    D.require(type(index) is int and 0 <= index < 144, "explicit bulk element index required")
    old_name = f"{ORIGINAL}/bulk/element_{index:03d}.json"
    old_hash = data.old_files[old_name]
    old = D._bound(old_name, old_hash, {})
    case = data.domains["bulk_selectors"][index]
    body = data.body["triangles"][index]
    context = data.context
    ids = D._packet_meta(context, "bulk", index, case, old, data.body, data.protocol)
    fields, nodes = E._bulk_fields(context.target, data.center, index)
    D.require(fields["global_indices"] == ids and fields["free_local"] == old["free_local"],
              "frozen exact center and uniform H must share the same local coordinates")
    hrows, prior_force = old["H_integral_minus_cache_upper"], old["force_integral_minus_cache"]
    scaled_hash = None
    if fields["chart"]["kind"] == "off_axis":
        binding = data.scaled_files[index]
        scaled_hash = binding["sha256"]
        scaled = D._bound(binding["file"], scaled_hash, {})
        prior_force, hrows = B._scaled_packet(old, scaled, data.proofs[index], data.expected, old_hash, ids)
    else:
        D.require(index not in data.scaled_files and index in B.AXIS, "retain exactly the reviewed twelve axis packets")
    H = {key: D._packed(row["enclosure"]) for key, row in D._entry_map(hrows, len(ids)).items()}
    D.require(len(H) == len(ids) * (len(ids) + 1) // 2, "complete local upper H defect required")
    rule_hash = E._bulk_rule(context.target, data.rule, index)
    if rule_hash != case["auxiliary_element_rule_sha256"] or not body["positive_J_for_entire_box"]:
        raise ValueError("exact same strictly positive auxiliary rule and fresh J domain required")
    local_radii = [
        Q(data.box["coordinate_radii_exact"][nodes[j // 2]][j % 2])
        for j in fields["free_local"]
    ]
    for k, node in enumerate(nodes):
        for axis in (0, 1):
            y = data.center[node][axis]
            r = Q(data.box["coordinate_radii_exact"][node][axis])
            D.require(r >= 0 and D._contains(E.SUR.rational_hull(data.wide[node][axis]), (y - r, y + r)),
                      "center and whole canonical local nodal segment must stay in original box")
            if 2 * k + axis not in fields["free_local"]:
                D.require(r == 0, "constrained coordinate must stay fixed along path")
    # Verify the exact all-column map for active-gain derivative in this
    # chart. No auxiliary map evaluations or quadrature sums are needed.
    dz = axial_gradient(fields)
    for col, free in enumerate(fields["free_local"]):
        k = 1 if free % 2 == 0 else 3
        D.require(fields["B"][k][col] == dz[free // 2],
                  "same axial P2 derivative in selected chart required")
    jet = JET.bulk_jacobian_jet(fields["q"], fields["B"])
    norms = [E._norm(p, fields) for p in jet["gradient"]]
    first, mass = (Q(case[k]) for k in ("phi_prime_tail_upper_exact", "geometry_two_mass_upper_exact"))
    tails = [first * mass * norm for norm in norms]
    D.require(max(tails) <= Q(case["force_tail_budget_exact"]), "original point force tail budget exceeded")
    proof = {
        "element": index, "chart": fields["chart"]["kind"], "active": str(fields["active"]),
        "global_indices": ids, "free_local": fields["free_local"], "local_nodes": nodes,
        "point_field_and_B_sha256": _polynomials_digest(fields),
        "uniform_hessian_endpoint_sha256": D._digest([(i, j, H[i, j]) for i, j in sorted(H)]),
        "prior_force_endpoint_sha256": D._digest(list(map(D._packed, prior_force))),
        "original_packet_sha256": old_hash, "reviewed_scaled_packet_sha256": scaled_hash,
        "original_selector_sha256": D._digest(case), "original_element_geometry_sha256": D._digest(body),
        "normalized_Qstar_rule_sha256": rule_hash,
        "order": case["degree"], "rho_exact": case["rho_exact"], "point_force_tail_max_exact": max(tails),
        "original_force_tail_budget_exact": case["force_tail_budget_exact"], "original_exact_radii": local_radii,
    }
    return proof, fields, jet, nodes, case, H, prior_force, local_radii


def _compute(data, proof):
    fresh, fields, jet, nodes, case, H, prior_force, radii = _element(data, proof["element"])
    D.require(D._digest(fresh) == D._digest(proof), "frozen centered panel proof changed")
    target = data.context.target
    element = proof["element"]
    functional = E._functional_type(fields)(target, fields, element, data.rule)
    local_center = [A(data.center[node][axis]) for node in nodes for axis in (0, 1)]
    local_wide = [A(data.wide[node][axis]) for node in nodes for axis in (0, 1)]
    finite_target = U.RuleTarget(target, data.rule)
    maps = [U.C.exact_bulk_maps(finite_target, element, iq) for iq in range(64)]
    try:
        point = F._force_parts(fields, jet, functional, local_center, A(fields["active"] * data.context.drive), maps, case)
        gain = affine_gain_derivative(data.context.dose, fields["active"], fields, functional, maps, local_wide)
        gain_offset = A(data.gain_bounds[0]).union(A(data.gain_bounds[1]))
        displacements = [A(-r).union(A(r)) for r in radii]
        combined, state, gain_effect = centered_force(point["normalized_force_total"], displacements, H, gain, gain_offset)
        intervals = [core.pack(v) for v in combined]
        intersected = [B.intersection(a, b) for a, b in zip(intervals, prior_force)]
        return {
            "status": "CENTERED_BULK_FORCE_COMPONENT_UNREVIEWED", "record_id": RECORD,
            "dose_exact": data.context.dose, "parameter_cell_sha256": data.context.row["parameter_cell_sha256"],
            "precision_bits": 640, **proof, "local_proof_sha256": D._digest(proof),
            "point_force": F._pack_profile(point),
            "same_original_box_Hessian_times_displacement": [core.pack(v) for v in state],
            "affine_gain_derivative_at_Y": [core.pack(v) for v in gain],
            "affine_gain_increment": [core.pack(v) for v in gain_effect],
            "centered_force_enclosure": intervals, "force_intersection_with_prior": intersected,
            "moments": len(functional.requested_defect_keys), "new_Hessian_integrations": False,
        }
    finally:
        for method in (E.P.Functional.integral, E.P.Functional.defect, E.P.Functional.node_power,
                       E.AXbridge.PhysicalFunctional.integral, M._triangle_moment, core.A):
            method.cache_clear()


def tests():
    checks = {}
    base = [[Q(0), Q(0)], [Q(2), Q(0)], [Q(2), Q(1)]]
    shifted = [[r + Q(1), z] for r, z in base]
    for kind, corners in (("axis_vertex", base), ("off_axis", shifted),
                          ("axis_edge_physical", [[Q(0), Q(0)], [Q(2), Q(0)], [Q(0), Q(1)]])):
        nodes = corners + [[(corners[i][k] + corners[j][k]) / 2 for k in (0, 1)] for i, j in ((0, 1), (1, 2), (2, 0))]
        free = [i for i in range(12) if i % 2 or nodes[i // 2][0] > 0]
        fields = E.AX.physical_fields(corners, nodes, free) if kind == "axis_edge_physical" else T.element_fields(corners, nodes, free)
        basis = axial_gradient(fields)
        checks[kind + ": exact axial basis equals map"] = all(
            fields["B"][1 if q % 2 == 0 else 3][i] == basis[q // 2] for i, q in enumerate(free))
        y = [A(v) for node in nodes for v in node]
        killed = SimpleNamespace(apply=lambda p: (_ for _ in ()).throw(AssertionError("inactive gain must be exact zero")))
        checks[kind + ": inactive gain zero before moments"] = all(
            v.is_zero() for v in affine_gain_derivative(Q(2), 0, fields, killed, [], y))
    force, state, gain = centered_force([A(2), A(-3)], [A(-2), A(1)],
                                        {(0, 0): (Q(1), Q(1)), (0, 1): (Q(-2), Q(-2)), (1, 1): (Q(3), Q(3))},
                                        [A(4), A(-2)], A(Q(1, 4)))
    checks["signed_H_and_gain_paths"] = [E.SUR.rational_hull(v) for v in force] == [(Q(-1), Q(-1)), (Q(7, 2), Q(7, 2))]
    crossing = T.constant(arb(0, A(Q(1, 4))))
    checks["signed_uncertain_coefficient_kept"] = (0, 0) in crossing and not crossing[0, 0].is_zero()
    checks.update({f"reviewed_force_only_{label}": result for label, result in F.tests().items()})
    D.require(D._yes(checks), "all centered method synthetic guards required")
    return checks


def prepare():
    forbidden = AssertionError("preparation called actual moments, cache sites or force evaluation")
    with patch.object(E.P.Functional, "__init__", side_effect=forbidden), \
            patch.object(M, "triangle_moment", side_effect=forbidden), \
            patch.object(U.C, "exact_bulk_maps", side_effect=forbidden), \
            patch.object(F, "_force_parts", side_effect=forbidden):
        data = _source_inputs()
        proofs = [_element(data, i)[0] for i in range(144)]
        checks = tests()
    counts = {kind: sum(v["chart"] == kind for v in proofs)
              for kind in ("off_axis", "axis_vertex", "axis_edge_physical")}
    checks.update({
        "all144_point_maps_local_H_measures_and_tails": [v["element"] for v in proofs] == list(range(144)),
        "expected_three_chart_counts": counts == {"off_axis": 132, "axis_vertex": 6, "axis_edge_physical": 6},
        "active_mask_24_elements_and_all_remaining_zero_gain": sum(v["active"] == "1" for v in proofs) == 24,
        "same_exact_center_and_in_cell_finite_parameter_path": data.gain_bounds[0] <= 0 <= data.gain_bounds[1],
        "no_actual_moments_cache_sites_or_point_force_sums": True,
    })
    D.require(D._yes(checks), "complete 144-map centered preparation required")
    names = ("exact_centered_bulk_L_frozen.py", "selectors.json", "protocol.json", "preflight.json")
    D.require(not any((PREP / name).exists() for name in names), "preserve write-once centered preparation")
    PREP.mkdir(parents=True, exist_ok=True)
    with (PREP / names[0]).open("xb") as fp:
        fp.write(Path(__file__).read_bytes())
    selectors_sha = D._write(PREP / names[1], proofs)
    protocol = {
        "status": "FROZEN_CENTERED_BULK_L_PREPARATION_REVIEW_PENDING", "record_id": RECORD,
        "source_file": _relative(__file__), "source_sha256": _sha(__file__), "input_hashes": data.hashes,
        "selectors_file": _relative(PREP / names[1]), "selectors_sha256": selectors_sha,
        "point_center_sha256": D._digest(data.center), "exact_base_gain": data.context.base["common_gain"],
        "exact_source_dose": data.context.dose, "closed_gain_shift_exact": data.gain_bounds,
        "same_parameter_cell_sha256": data.context.row["parameter_cell_sha256"],
        "original_coordinate_box_sha256": data.protocol["coordinate_box_sha256"],
        "chart_counts": counts, "precision_bits": 640, "elements": 144,
        "point_policy": "At exact stored Y0 and exact base g0*dose*active, compute local integral-cache force only. "
                        "Use original order/rho, analytic force tail, same-Y same-rule finite cache correction, "
                        "physical axis-edge or regularized axis-vertex chart as applicable. No H moment integration.",
        "centered_policy": "First vary Y from Y0 at g0. Eliminate this segment via local uniform H from reviewed original "
                           "axis packets or reviewed original/scaled intersection on off-axis cells. Then vary gain "
                           "at fixed final Y. The activation defect is affine in gain, with a signed axial stiffness "
                           "difference using exact physical P2 dZ bases and same frozen normalized cache. Inactive gain term is zero.",
        "composition": "Intersect centered force with the previously reviewed force on the same element, full box and C. "
                       "Keep the corresponding original/combined Hessian entry unchanged. A later independent global "
                       "reduction can combine these packets with reviewed surface components. No root assertion here.",
        "actual_expected_files": [_relative(STAGE / "actual" / "elements" / f"element_{i:03d}.json") for i in range(144)],
        "actual_record_force_integrals_run": False,
        "failed_attempts": {_relative(p): _sha(p) for p in sorted(PREP.glob("development_revision_*/*")) if p.is_file()},
    }
    psha = D._write(PREP / names[2], protocol)
    pre = {"status": "PASS_CENTERED_BULK_L_PREPARATION_ONLY", "source_sha256": protocol["source_sha256"],
           "protocol_sha256": psha, "selectors_sha256": selectors_sha, "checks": checks,
           "elements": 144, "actual_record_force_integrals_run": False, "actual_cache_site_sums_run": False}
    presha = D._write(PREP / names[3], pre)
    print(json.dumps({"status": pre["status"], "source_sha256": pre["source_sha256"], "protocol_sha256": psha,
                      "preflight_sha256": presha, "selectors_sha256": selectors_sha, "checks": len(checks)}), flush=True)


def run(review_sha):
    review = D._bound(_relative(PREP / "prelaunch_independent_review.json"), review_sha, {})
    hashes = {"source_sha256": _sha(__file__), "protocol_sha256": _sha(PREP / "protocol.json"),
              "preflight_sha256": _sha(PREP / "preflight.json"), "selectors_sha256": _sha(PREP / "selectors.json")}
    D.require(review["status"] == "PASS" and D._yes(review["checks"])
              and all(review[k] == v for k, v in hashes.items()), "bound independent centered-L bulk prelaunch required")
    p = json.loads((PREP / "protocol.json").read_bytes())
    proofs = json.loads((PREP / "selectors.json").read_bytes())
    data = _source_inputs()
    D.require(p["input_hashes"] == data.hashes and p["selectors_sha256"] == hashes["selectors_sha256"]
              and [v["element"] for v in proofs] == list(range(144)), "complete frozen 144-element centered domain required")
    out = STAGE / "actual"
    (out / "elements").mkdir(parents=True, exist_ok=True)
    utc = lambda: datetime.now(timezone.utc).isoformat()
    D._write(out / "started.json", {**hashes, "prelaunch_sha256": review_sha, "pid": os.getpid(), "utc": utc()})
    start = monotonic()
    packets = []
    try:
        for proof in proofs:
            i = proof["element"]
            value = _compute(data, proof)
            path = out / "elements" / f"element_{i:03d}.json"
            digest = D._write(path, {**value, **hashes, "prelaunch_sha256": review_sha})
            packets.append({"element": i, "file": _relative(path), "sha256": digest})
            print(json.dumps({"element": i, "complete": len(packets), "total": 144,
                              "elapsed_seconds": monotonic() - start}), flush=True)
        D.require([v["file"] for v in packets] == p["actual_expected_files"], "complete exact centered 144-packet manifest required")
        digest = D._write(out / "results.json", {"status": "ALL_144_CENTERED_BULK_FORCE_PACKETS_UNREVIEWED",
                                                "record_id": RECORD, "files": packets, **hashes,
                                                "prelaunch_sha256": review_sha, "new_Hessian_integrations": False,
                                                "root_asserted": False, "elapsed_seconds": monotonic() - start})
        D._write(out / "terminal.json", {"status": "FINISHED", "exit_code": 0, "pid": os.getpid(),
                                        "utc": utc(), "results_sha256": digest})
    except BaseException as error:
        D._write(out / "failure.json", {"error": str(error), "completed_packets": packets, **hashes})
        D._write(out / "terminal.json", {"status": "FAILED", "exit_code": 1, "pid": os.getpid(), "utc": utc()})
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
