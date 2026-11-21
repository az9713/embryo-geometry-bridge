"""Explicit-record centered force-only bulk component on a reviewed original cube."""
from dataclasses import dataclass
from fractions import Fraction as Q
from hashlib import sha256
import importlib
import json
from pathlib import Path
from time import monotonic
from unittest.mock import patch

from flint import arb, ctx as precision


ROOT = Path(__file__).resolve().parents[1]
BASE = "runs/exact_integration_cohort_v1"
OUT = ROOT / BASE / "remaining15_centered_preparation"
DIRECT = BASE + "/remaining15_reduction_preparation"
L_CENTRED = "scripts/exact_centered_bulk_L.py"
REDUCER = "scripts/exact_cohort_record_reduction.py"
L_REDUCER = BASE + "/L_defect_diagnostic/scaled_offaxis_L_v1/reduction_preparation/reduce_scaled_L.py"
CD = {
    "optoGEF2_small_col_C": "b66b821666c74ee68d2c0105df3087c4bec9e3a408d92f14cfa88231c9fc7d08",
    "optoGEF2_small_col_D": "03017ddf88803327d2bf24cb4d22c8b4957ae4969d999ed062c66924d1e5b411",
}
PINS = {
    L_CENTRED: "a83c89ffc302be341df220fcb97c5be82f191a70e317fea9af839d2917a87492",
    REDUCER: "8e3177434299fbffeee2fa5161dd4025fc0e9cd0b7f46be85fceaeccae204ed2",
    L_REDUCER: "0ebaf0c1dfd65a3962841be51b714f3540b6a7f07e226dab130591bd4b7f70e9",
    DIRECT + "/independent_actual/CD_width_results.json": "faf52e086282f86ad2632b6f932d9a438442ef86916aaf890188d8f764e386f6",
}


@dataclass(frozen=True)
class Domain:
    context: object
    nodes: list
    box: dict
    rule: dict
    body: dict
    protocol: dict


def require(condition, why):
    if not condition:
        raise ValueError(why)


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _relative(path):
    return Path(path).relative_to(ROOT).as_posix()


def _boot():
    for name, digest in PINS.items():
        require(_sha(ROOT / name) == digest, f"frozen centered prerequisite differs: {name}")
    R = importlib.import_module("exact_cohort_record_reduction")
    require(Path(R.__file__).resolve() == (ROOT / REDUCER).resolve(), "wrong generic reviewed reducer module")
    op, plan, jobs, expect, hashes = R._boot()
    L = importlib.import_module("exact_centered_bulk_L")
    require(Path(L.__file__).resolve() == (ROOT / L_CENTRED).resolve()
            and Path(L.B.__file__).resolve() == (ROOT / L_REDUCER).resolve()
            and L.F is op.F and L.D is op.D and L.C is op.C and L.E is op.E
            and L.U is op.E.U and L.M is op.E.M and L.JET is op.E.JET and L.T is op.E.T,
            "reviewed force-only pure helper aliases differ")
    for n, h in PINS.items():
        op.D._bound(n, h, hashes, parse=False)
    return R, op, L, plan, jobs, expect, hashes


def _intersect(D, a, b):
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    require(lo <= hi, "disjoint prior-uniform and centered force enclosures")
    return lo, hi


def _record_gain(domain):
    ctx = domain.context
    g0 = ctx.base["common_gain"]
    lo, hi = ctx.cell["common_gain"]
    require(0 < lo <= g0 <= hi and ctx.dose > 0 and ctx.drive == g0 * ctx.dose, "same in-cell exact base gain/dose required")
    return lo - g0, hi - g0


def _local(runtime, domain, selector, route, uniform, uniform_sha, index):
    R, op, L = runtime
    D, E, T = op.D, op.E, op.E.T
    ctx = domain.context
    require(type(index) is int and index == selector["element"] == route["element"] and 0 <= index < 144
            and uniform["record_component_proof_sha256"] == D._digest(route)
            and route["body_element_sha256"] == D._digest(domain.body["triangles"][index])
            and route["original_selector_sha256"] == D._digest(selector),
            "same record-specific terminal operator route proof required")
    if route["integration_chart"] == "scaled_positive_R_affine_material":
        ids = R._scaled_packet(D, ctx, index, selector, domain.body["triangles"][index], route, uniform, domain.protocol)
    else:
        ids = D._packet_meta(ctx, "bulk", index, selector, uniform, domain.body, domain.protocol)
        require(uniform["integration_chart"] == route["integration_chart"] == selector["chart"]
                and selector["chart"] in ("axis_vertex", "axis_edge_physical"), "same original axis kernel required")
    center = ctx.target.nodes()
    fields, nodes = E._bulk_fields(ctx.target, center, index)
    require(fields["global_indices"] == ids == route["global_indices"]
            and fields["free_local"] == route["free_local"] and fields["chart"]["kind"] == selector["chart"]
            and fields["active"] in (0, 1), "exact center, old uniform box and route must use one free projection")
    radii = []
    for k, node in enumerate(nodes):
        for a in (0, 1):
            center_value = center[node][a]
            radius = Q(domain.box["coordinate_radii_exact"][node][a])
            require(Q(domain.box["current_centers_exact"][node][a]) == center_value and radius >= 0
                    and D._contains(E.SUR.rational_hull(domain.nodes[node][a]), (center_value - radius, center_value + radius)),
                    "whole line from the exact selected center must lie in the admitted original nodal box")
            if 2 * k + a not in fields["free_local"]:
                require(radius == 0, "fixed coordinate lies outside its declared subspace")
            else:
                radii.append(radius)
    h = D._entry_map(uniform["H_integral_minus_cache_upper"], len(ids))
    H = {k: D._packed(v["enclosure"]) for k, v in h.items()}
    require(len(H) == len(ids) * (len(ids) + 1) // 2, "complete independently reviewed local H enclosure required")
    force = list(map(D._packed, uniform["force_integral_minus_cache"]))
    require(len(force) == len(ids), "complete independently reviewed local uniform force required")
    dz = L.axial_gradient(fields)
    require(all(fields["B"][1 if q % 2 == 0 else 3][i] == dz[q // 2] for i, q in enumerate(fields["free_local"])),
            "physical axial derivative must agree with the original finite-element projection in this chart")
    rule_digest = E._bulk_rule(ctx.target, domain.rule, index)
    require(rule_digest == selector["auxiliary_element_rule_sha256"] == route["auxiliary_element_rule_sha256"]
            and selector["nodal_hulls_sha256"] == route["nodal_hulls_sha256"],
            "same original Qstar/cache rule and uniform nodal domain required")
    jet = E.JET.bulk_jacobian_jet(fields["q"], fields["B"])
    mass, first = (Q(selector[k]) for k in ("geometry_two_mass_upper_exact", "phi_prime_tail_upper_exact"))
    tails = [mass * first * E._norm(p, fields) for p in jet["gradient"]]
    require(max(tails) <= Q(selector["force_tail_budget_exact"]), "point force tail exceeds original per-element allowance")
    proof = {
        "record_id": ctx.record_id, "element": index, "chart": fields["chart"]["kind"],
        "uniform_source_integration_chart": route["integration_chart"],
        "global_indices": ids, "free_local": fields["free_local"], "local_nodes": nodes,
        "dose_exact": str(ctx.dose), "base_drive_exact": str(ctx.drive), "active": str(fields["active"]),
        "common_gain_shift_exact": list(map(str, _record_gain(domain))),
        "same_C_sha256": ctx.row["parameter_cell_sha256"], "same_cache_sha256": ctx.row["point_cache_sha256"],
        "same_predictor_sha256": ctx.row["predictor_sha256"], "point_polynomials_sha256": L._polynomials_digest(fields),
        "uniform_packet_sha256": D._sha(uniform_sha), "original_route_proof_sha256": D._digest(route),
        "uniform_force_bounds_sha256": D._digest(force),
        "uniform_H_bounds_sha256": D._digest([(i, j, H[i, j]) for i, j in sorted(H)]),
        "original_local_radii_exact": list(map(str, radii)),
        "auxiliary_rule_sha256": rule_digest, "order": selector["degree"], "rho_exact": selector["rho_exact"],
        "point_force_tail_max_exact": str(max(tails)), "force_tail_budget_exact": selector["force_tail_budget_exact"],
    }
    return proof, fields, nodes, center, jet, H, force, radii


def prepare_component(runtime, domain, selector, route, uniform, uniform_sha, element):
    """Freeze a record-specific point/segment proof without numerical component work."""
    return _local(runtime, domain, selector, route, uniform, uniform_sha, element)[0]


def bulk_component(runtime, domain, selector, route, uniform, uniform_sha, element, proof):
    """Evaluate one reviewed centered bulk force component after an external launch gate."""
    before = monotonic()
    _, op, L = runtime
    D, E, F = op.D, op.E, op.F
    A, core, U = E.U.A, E.core, E.U
    expected, fields, nodes, center, jet, H, prior, radii = _local(runtime, domain, selector, route, uniform, uniform_sha, element)
    require(D._digest(expected) == D._digest(proof), "bound record-specific centered launch proof changed")
    ctx = domain.context
    fun = E._functional_type(fields)(ctx.target, fields, element, domain.rule)
    try:
        base = [A(center[n][i]) for n in nodes for i in (0, 1)]
        full_Y = [A(domain.nodes[n][i]) for n in nodes for i in (0, 1)]
        as_rule = U.RuleTarget(ctx.target, domain.rule)
        maps = [U.C.exact_bulk_maps(as_rule, element, iq) for iq in range(64)]
        point = F._force_parts(fields, jet, fun, base, A(fields["active"] * ctx.drive), maps, selector)
        gain = L.affine_gain_derivative(ctx.dose, fields["active"], fields, fun, maps, full_Y)
        lo, hi = _record_gain(domain)
        intervals, state, effect = L.centered_force(point["normalized_force_total"], [A(-r).union(A(r)) for r in radii],
                                                   H, gain, A(lo).union(A(hi)))
        packed = [core.pack(v) for v in intervals]
        intersection = [_intersect(D, D._packed(p), q) for p, q in zip(packed, prior)]
        return {
            **proof, "status": "CENTERED_RECORD_BULK_FORCE_UNREVIEWED", "precision_bits": 640,
            "centered_proof_sha256": D._digest(proof), "point_force": F._pack_profile(point),
            "uniform_H_times_nodal_displacement": [core.pack(v) for v in state],
            "affine_gain_derivative_at_Y": [core.pack(v) for v in gain], "affine_gain_increment": [core.pack(v) for v in effect],
            "centered_force_enclosure": packed, "force_intersection_with_prior_exact": intersection,
            "new_Hessian_integrations": False, "unique_point_and_gain_moments": len(fun.requested_defect_keys),
            "seconds": monotonic() - before,
        }
    finally:
        for method in (E.P.Functional.integral, E.P.Functional.defect, E.P.Functional.node_power,
                       E.AXbridge.PhysicalFunctional.integral, E.M._triangle_moment, core.A):
            method.cache_clear()


def bind_reviewed_record(runtime, plan, jobs, expect, record, review_sha, hashes):
    """Open only a completed record whose full actual component/reduction review is pinned."""
    R, op, _ = runtime
    D = op.D
    approved_file = f"{DIRECT}/independent_actual/{record}.json"
    audit = D._bound(approved_file, review_sha, hashes)
    require(audit["status"].startswith("PASS_DIRECT_REDUCTION_") and audit["record_id"] == record
            and audit["F_entries"] == 611 and audit["H_upper_entries"] == 6676
            and audit["packet_route_counts"] == {"bulk": 144, "surface": 24, "axis_vertex": 6,
                                               "axis_edge_physical": 6, "scaled": 132},
            "completed independent actual whole direct operator/reduction required")
    pins = audit["producer_pins"]
    D._bound(pins["observation_file"], pins["observation_sha256"], hashes, parse=False)
    original_result = D._bound(f"{DIRECT}/{record}/results.json", pins["results_sha256"], hashes)
    terminal = D._bound(f"{DIRECT}/{record}/terminal.json", pins["terminal_sha256"], hashes)
    D._terminal(terminal, pins["results_sha256"])
    admission = D._bound(f"{DIRECT}/{record}_input_admission.json", pins["admission_sha256"], hashes)
    require(original_result["record_id"] == admission["record_id"] == record
            and original_result["admission_sha256"] == pins["admission_sha256"]
            and audit["same_closed_C_sha256"] == original_result["same_closed_C_sha256"]
            and original_result["component_count"] == 168
            and admission["status"] == "ALL_168_EXPECTED_PACKETS_BOUND_BEFORE_LOCAL_READ",
            "same actual scoped proof and full 168 terminal input admission required")
    ctx, p, domains, body, routes, links, stage = R._load_record(op, plan, jobs, expect, record, hashes)
    packets = {part: R._terminal_part(D, stage, links, record, part, admission["external_part_pins"][part], hashes)
               for part in R.PARTS}
    require(packets == admission["expected_packets"] and all(original_result["input_hashes"][row["file"]] == row["sha256"]
                                                            for seq in packets.values() for row in seq),
            "prior independently reviewed raw packet identities changed")
    box = D._bound(p["coordinate_box_file"], p["coordinate_box_sha256"], hashes)
    domain = Domain(ctx, op.E.SUR.common_nodes(ctx.target, box), box, ctx.target.auxiliary_rule(), body, p)
    _record_gain(domain)
    return domain, domains, routes, packets, original_result


def _static_profiles(op):
    result = []
    for record in list("CDEFGHIJKMNOPQR"):
        ctx = op.C.load_context("optoGEF2_small_col_" + record)
        a = ctx.target.a
        k = [sum(a["reference"][int(node), 0] == 0 for node in tri[:3]) for tri in a["triangles"]]
        require(len(k) == 144 and k.count(0) == 132 and k.count(1) == k.count(2) == 6
                and ctx.cell["common_gain"][0] <= ctx.base["common_gain"] <= ctx.cell["common_gain"][1]
                and ctx.drive == ctx.base["common_gain"] * ctx.dose,
                "all selected exact base gains and chart routes required")
        result.append({"record_id": ctx.record_id, "context_input_hashes": ctx.input_hashes,
                       "exact_dose": ctx.dose, "base_gain_exact": ctx.base["common_gain"],
                       "same_C_sha256": ctx.row["parameter_cell_sha256"],
                       "routes": {"off_axis": 132, "axis_vertex": 6, "axis_edge_physical": 6},
                       "exact_original_radius": ctx.old_radius, "terminal_actual_packets_opened": ctx.record_id in CD})
    return result


def prepare():
    R, op, L, plan, jobs, expect, hashes = _boot()
    runtime = R, op, L
    D, E, F = op.D, op.E, op.F
    precision.prec = 640
    E.core.A.cache_clear()
    no = AssertionError("attempted actual force, moment or cache evaluation during centered preparation")
    with (
        patch.object(F, "_force_parts", side_effect=no),
        patch.object(E.P.Functional, "__init__", side_effect=no),
        patch.object(E.U.C, "exact_bulk_maps", side_effect=no),
        patch.object(E.U, "uniform_cache", side_effect=no),
        patch.object(E.M, "triangle_moment", side_effect=no),
    ):
        checks = L.tests()
        profiles = _static_profiles(op)
        proofs = {}
        for record, sha in CD.items():
            domain, actual, routes, packets, original = bind_reviewed_record(runtime, plan, jobs, expect, record, sha, hashes)
            payloads = [D._bound(row["file"], row["sha256"], hashes) for row in packets["bulk"]]
            proofs[record] = {
                "record_id": record, "actual_direct_review_sha256": sha, "closed_cell_sha256": domain.context.row["parameter_cell_sha256"],
                "bulk_proofs": [prepare_component(runtime, domain, actual["bulk_selectors"][i], routes["bulk"][i],
                                                  p, packets["bulk"][i]["sha256"], i) for i, p in enumerate(payloads)],
                "retained_surface_manifest": packets["surface"], "retained_upper_Hessian_entries": original["upper_Hessian_entries"],
                "retained_fullbox_b_exact": original["b_upper_exact"], "new_force_or_cache_sums_started": False,
            }
    checks.update({"generic_static_full15_contexts": len(profiles) == 15,
                   "288_real_exact_point_projections_and_H_packet_bindings": all(
                       [r["element"] for r in p["bulk_proofs"]] == list(range(144)) for p in proofs.values()),
                   "only_completed_and_reviewed_CD_terminal_payloads_opened": set(proofs) == set(CD),
                   "no_actual_new_force_H_moments_or_cache_work": True})
    require(D._yes(checks), "all generic centered force prelaunch guards required")
    OUT.mkdir(parents=True, exist_ok=True)
    names = ("exact_centered_cohort_bulk_frozen.py", "protocol.json", "preflight.json")
    require(not any((OUT / k).exists() for k in names), "preserve write-once centered preparation")
    with (OUT / names[0]).open("xb") as stream:
        stream.write(Path(__file__).read_bytes())
    records = {}
    for rid, p in proofs.items():
        path = OUT / f"{rid}_proofs.json"
        records[_relative(path)] = D._write(path, p)
    protocol = {
        "status": "FROZEN_REMAINING15_CENTERED_HELPER_PREPARATION_REVIEW_PENDING",
        "source_file": _relative(__file__), "source_sha256": _sha(__file__), "input_hashes": hashes,
        "static_profiles": profiles, "completed_record_proof_hashes": records, "precision_bits": 640,
        "scope": "All 144 bulk charts in each record's own original physical box and C. Only C/D have bound actual prerequisites now. "
                 "Others require completed bulk/surface terminals and actual independent review before binding their proofs.",
        "math": "Signed point force at exact stored Y0 and base gain. Add independently reviewed whole-box H times the local "
                "nodal increment at fixed g0 and the affine gain derivative times exact gain shift at fixed Y. Intersect with "
                "prior signed force. Full future 168 assembly retains the existing 24 surface force/H and full original b/z.",
        "prelaunch_scope": "Independent operator and record-specific centered proof prelaunch before calling bulk_component. "
                           "No actual call in this preparation.",
        "actual_new_force_H_or_cache_sum_started": False,
        "failed_attempts": {_relative(p): _sha(p) for p in sorted(OUT.glob("development_revision_*/*")) if p.is_file()},
    }
    psha = D._write(OUT / names[1], protocol)
    fsha = D._write(OUT / names[2], {
        "status": "PASS_REMAINING15_CENTERED_HELPER_PREPARATION_ONLY", "source_sha256": protocol["source_sha256"],
        "protocol_sha256": psha, "checks": checks, "records_with_bound_144_proofs": list(proofs),
        "new_force_H_or_cache_sum_started": False,
    })
    print(json.dumps({"source_sha256": protocol["source_sha256"], "protocol_sha256": psha, "preflight_sha256": fsha,
                      "checks": len(checks), "bound_record_proofs": records}), flush=True)


if __name__ == "__main__":
    prepare()
