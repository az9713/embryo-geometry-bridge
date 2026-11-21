"""Review actual L domains before launching either integral part."""

import hashlib
import json
import math
import time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
from unittest.mock import patch

import rebind_radius_fix as prior

P, OUT = prior.P, prior.OUT
ROOT = P.ROOT
HERE = Path(__file__).resolve().parent
SOURCE = "a21b680117f236cc38efc46e086be1d9aff80cfb09ae613691b1642db8bed25e"
PROTOCOL = "66bfdb7dbe49d3ece787e7ecc9f16e004f7611a8b946be0e078f8a6806a4c937"
PREFLIGHT = "ad5e6dc6ca0d1fc0f8eafd58cd4691bb19d360f18a03dd4dd1e3dde63d09cb41"
DOMAINS = "7568831341a79caec28cb12c8a5921f40845a982cef4971b5c064c07bb8d1ac4"
BODY = "a26897e1f45513f6231c4ac81b58c8550b677a13b7398da1ac1003a029511f8e"
PREVIOUS_GATE = "29c15e6d0cae1428d2860950a66b27db957f24fdc8c5dea44bcb14d3205636e9"
GENERIC = {
    "evidence/exact_cohort_context_independent_review.json":
        "898fb6828c16fc7f0db410c9d677373cf1696db5d0d418c08ef0004a30e215c5",
    "evidence/exact_cohort_integrator_independent_review.json":
        "1d233c61e5b7283ddfeb856e206a064d28cb2731246f13bff716e9a4063dbfa8",
}
RECORD = "optoGEF2_small_col_L"


def q(value):
    if hasattr(value, "dtype") and value.dtype.kind == "f":
        return Q.from_float(float(value))
    return Q(value)


def digest(values):
    return hashlib.sha256(json.dumps(values, separators=(",", ":"), sort_keys=True).encode()).hexdigest()


def projection(context, nodes):
    indices = {int(v): i for i, v in enumerate(context.free)}
    dofs = [2 * int(node) + k for node in nodes for k in (0, 1)]
    free = [i for i, dof in enumerate(dofs) if dof in indices]
    return free, [indices[dofs[i]] for i in free]


def ceil_dyadic(value):
    return Q(math.ceil(value * 2**48), 2**48)


def upper(value):
    return q(str(value.abs_upper().fmpq()))


def bulk_tail(order, rho):
    def harmonic(n):
        return sum((Q(1, i) for i in range(1, n + 1)), Q(0))
    m = order + 1
    gap = 1 - rho
    term0, term1, term2 = (1/gap, rho/gap**2, rho*(1+rho)/gap**3)
    first = (1 + 5000*harmonic(m))*term0 + Q(5000, m+1)*term1
    A = 1 + 5000*harmonic(m+1)
    B = Q(5000, m+2)
    second = (m+1)*A*term0 + (A+(m+1)*B)*term1 + B*term2
    return rho**m*first, rho**m*second


def bulk_check(context, nodes, rule, case, body_elem, i, ft, ht):
    a = context.target.a
    ids = list(map(int, a["triangles"][i]))
    free, projected = projection(context, ids)
    corners = [[q(v) for v in a["reference"][node]] for node in ids[:3]]
    axis = sum(r == 0 for r, _ in corners)
    chart = ("off_axis", "axis_vertex", "axis_edge_physical")[axis]
    assert case["element"] == i and case["active"] == str(q(a["active"][i, 0]))
    assert case["chart"] == chart and case["global_indices"] == projected and case["free_local"] == free
    local_hulls = [[list(map(str, P.S.rational_hull(v))) for v in nodes[node]] for node in ids]
    assert case["nodal_hulls_sha256"] == digest(local_hulls)
    weights = [q(v) / context.target.energy_scale for v in a["volume_weights"][i]]
    sites = []
    for row in a["shape"]:
        n = list(map(q, row))
        xi, eta = n[1] + (n[3]+n[4])/2, n[2] + (n[4]+n[5])/2
        sites.append((1-xi-eta, xi, eta))
    assert len(weights) == len(sites) == 64 and min(weights) > 0 and all(min(v)>0 for v in sites)
    assert list(map(tuple, rule["volume_barycentric"])) == sites and rule["volume_weights"][i] == weights
    assert case["auxiliary_element_rule_sha256"] == digest([[[str(x) for x in s] for s in sites], list(map(str, weights))])
    norms = P.E._bulk_setup(context, nodes, rule, i)
    identity, _, _, _, shared_mass, grad, H = norms
    assert all(case[k] == v for k, v in identity.items())
    g, h = max(grad), max(H.values())
    assert g == q(case["max_DJ_norm_upper_exact"]) and h == q(case["max_D2J_norm_upper_exact"])
    (r0, z0), (r1, z1), (r2, z2) = corners
    det = (r1-r0)*(z2-z0)-(r2-r0)*(z1-z0)
    material_mass = context.target.normalization()*P.core.A(det*(r0+r1+r2)/6)
    mass = upper(material_mass+P.core.A(sum(weights)))
    assert q(case["geometry_two_mass_upper_exact"]) == shared_mass == mass > 0
    assert body_elem["positive_J_for_entire_box"] is True
    jlo, jhi = [q(body_elem["J_bound"][k]) for k in ("lower_exact", "upper_exact")]
    assert q(case["force_tail_budget_exact"]) == ft and q(case["hessian_tail_budget_exact"]) == ht
    raw = max(abs(jlo-1), abs(jhi-1))
    rho = ceil_dyadic(raw)
    assert q("0.000001") < jlo <= jhi and raw <= rho < 1
    assert list(map(q, case["J_bounds_exact"])) == [jlo, jhi]
    assert q(case["raw_rho_exact"]) == raw and q(case["rho_exact"]) == rho
    attempts = case["order_attempts"]
    assert 2 <= case["degree"] <= case["max_order"] == 32
    assert [row["degree"] for row in attempts] == list(range(2, case["degree"]+1))
    for j, row in enumerate(attempts):
        B, A = bulk_tail(row["degree"], rho)
        f_value, h_value = mass*B*g, mass*(A*g*g+B*h)
        assert q(row["phi_prime_tail_upper_exact"]) == B and q(row["phi_second_tail_upper_exact"]) == A
        assert q(row["force_tail_upper_exact"]) == f_value and q(row["H_entry_tail_upper_exact"]) == h_value
        assert (f_value <= ft and h_value <= ht) == (j == len(attempts)-1)
    assert all(case[k] == v for k, v in attempts[-1].items())
    return len(attempts), chart


def hull(intervals):
    return min(v[0] for v in intervals), max(v[1] for v in intervals)


def square(interval):
    lo, hi = interval
    return (Q(0) if lo <= 0 <= hi else min(lo*lo, hi*hi)), max(lo*lo, hi*hi)


def linear_hull(coefficients, vals):
    terms = [sorted((c*l, c*h)) for c, (l, h) in zip(coefficients, vals)]
    return sum(b[0] for b in terms), sum(b[1] for b in terms)


def cache_site(context, box, rule, ids, selected, global_indices, site, face, edge, iq):
    a = context.target.a
    n = list(map(q, a[f"surface_{face}_N"][iq]))
    d = list(map(q, a[f"surface_{face}_dN"][iq]))
    s = n[2]+n[1]/2
    lo, hi = (q(a["reference"][i, 0]) for i in (ids[0], ids[2]))
    L = hi-lo
    R = lo+L*s
    l0 = q(a[f"surface_{face}_l0"][edge])
    r0 = q(a[f"surface_{face}_r0"][edge, iq])
    assert 0 < s < 1 and L > 0 and R > 0 and r0 > 0 and l0 > 0
    assert s == rule["surface_s"][face][iq]
    D = (L/(r0*l0), 1/L, 1/L)
    star_n = (1-3*s+2*s*s, 4*s-4*s*s, 2*s*s-s)
    star_d = (-3+4*s, 4-8*s, 4*s-1)
    z = Q(0)
    old, canonical = [[z]*6 for _ in range(3)], [[z]*6 for _ in range(3)]
    for j in range(3):
        old[0][2*j] = n[j]*D[0]
        old[1][2*j] = d[j]/L
        old[2][2*j+1] = d[j]/L
        canonical[0][2*j] = star_n[j]/R
        canonical[1][2*j] = star_d[j]/L
        canonical[2][2*j+1] = star_d[j]/L
    Y = []
    for i in ids:
        for k in (0, 1):
            c = q(box["current_centers_exact"][i][k])
            r = q(box["coordinate_radii_exact"][i][k])
            Y.append((c-r, c+r))
    cached = [linear_hull(row, Y) for row in old]
    delta = [linear_hull([ss-cc for ss, cc in zip(A, B)], Y) for A, B in zip(canonical, old)]
    line = [(lo+min(0, dl), hi+max(0, dh)) for (lo, hi), (dl, dh) in zip(cached, delta)]
    tangent = square(line[1])[0] + square(line[2])[0]
    w = q(a[f"surface_{face}_weights"][edge, iq])/context.target.energy_scale
    assert site["quadrature_index"] == iq and site["nodes"] == ids
    assert site["selected_local"] == selected and site["global_indices"] == global_indices
    assert tuple(map(q, site["normalization_D_exact"])) == D
    assert [tuple(map(q, v)) for v in site["cache_q_exact"]] == cached
    assert [tuple(map(q, v)) for v in site["delta_q_exact"]] == delta
    assert [tuple(map(q, v)) for v in site["convex_line_exact"]] == line
    assert q(site["tangent_norm_squared_lower_exact"]) == tangent > 0
    assert q(site["auxiliary_normalized_weight_exact"]) == rule["surface_weights"][face][edge][iq] == w > 0
    assert line[0][0] > 0 and site["status"] == "CERTIFIED_LOCAL_FULLY_COEFFICIENTWISE_LINE_DOMAIN_NO_CACHE_SUM"


def half_tail(order, rho, three_halves=False):
    m = order+1
    c = Q(math.comb(2*m, m), 4**m)
    if not three_halves:
        return c*rho**m/(1-rho)
    ratio = rho*Q(2*m+3, 2*m+2)
    bounds = [rho**m*(Q(2*m+1)/(1-rho)+2*rho/(1-rho)**2)]
    if ratio < 1:
        bounds.append(c*(2*m+1)*rho**m/(1-ratio))
    return min(bounds)


def surface_check(context, nodes, box, rule, case, face, edge, ft, ht):
    a = context.target.a
    ids = list(map(int, a[f"surface_{face}_edges"][edge]))
    free, global_indices = projection(context, ids)
    fields = context.target.surface(face, edge, nodes=nodes)
    assert case["record_id"] == RECORD and (case["face"], case["edge"]) == (face, edge)
    assert case["global_indices"] == fields["global_indices"] == global_indices
    actual = P.S.uniform_G_bounds(fields, max_depth=8)
    old = case["domain"]
    assert all(old[k] == v for k, v in actual.items() if k != "rho_exact")
    assert old["status"] == "CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN"
    width = 2**old["depth"]
    assert 0 <= old["depth"] <= 8 and len(old["leaves"]) == width
    G = []
    hoop = []
    for i, leaf in enumerate(old["leaves"]):
        assert list(map(q, leaf["s_exact"])) == [Q(i, width), Q(i+1, width)]
        A, B, C = (hull([tuple(map(q, v)) for v in part])
                   for part in leaf["q_numerator_Bernstein_coefficients_exact"])
        den = hull([tuple(map(q, v)) for v in leaf["hoop_denominator_Bernstein_coefficients_exact"]])
        assert den[0] > 0
        ratios = [x/y for x in A for y in den]
        A = min(ratios), max(ratios)
        assert [tuple(map(q, v)) for v in leaf["q_bounds_exact"]] == [A, B, C]
        aa, bb, cc = map(square, (A, B, C))
        predicted = aa[0]*(bb[0]+cc[0]), aa[1]*(bb[1]+cc[1])
        assert tuple(map(q, leaf["G_bounds_exact"])) == predicted
        G.append(predicted)
        hoop.append(A[0])
    lo, hi = hull(G)
    assert 0 < lo <= hi and min(hoop) == q(old["hoop_lower_exact"]) > 0
    assert tuple(map(q, old["G_bounds_exact"])) == (lo, hi)
    q0 = (lo+hi)/2
    raw = (hi-lo)/(hi+lo)
    rho = ceil_dyadic(raw)
    assert q0 == q(old["q0_exact"]) and q(actual["rho_exact"]) == q(old["raw_rho_exact"]) == raw
    assert raw <= rho == q(old["rho_exact"]) < 1
    material = context.target.normalization()*P.core.A((fields["Rhi"]**2-fields["Rlo"]**2)/2)
    mass = upper(material + P.core.A(sum(rule["surface_weights"][face][edge])))
    sqrtlo = q(str(P.core.A(q0).sqrt().lower().fmpq()))
    assert 0 < sqrtlo and sqrtlo*sqrtlo <= q0
    prefix = ("basal", "apical")[face]
    K, T = (context.cell[prefix + suffix] for suffix in ("_area_modulus", "_surface_tension"))
    dmax = max(abs(t-k) for t in T for k in K)
    jet = P.E.JET.surface_squared_jacobian_jet(fields["q_physical_R"], fields["B_physical_R"])
    norm = lambda poly: P.S.poly_norm_upper(poly, fields["Rlo"], fields["Rhi"])
    g = max(map(norm, jet["gradient"]))
    h = max(map(norm, jet["hessian_upper"].values()))
    records = case["order_attempts"]
    order = case["order"]
    assert 2 <= order <= 64 and [row["degree"] for row in records] == list(range(2, order+1))
    for i, row in enumerate(records):
        A = dmax*half_tail(row["degree"], rho)/(2*sqrtlo)
        B = dmax*half_tail(row["degree"], rho, three_halves=True)/(4*q0*sqrtlo)
        f_value, h_value = mass*A*g, mass*(A*h+B*g*g)
        assert q(row["force_tail_upper_exact"]) == f_value and q(row["H_entry_tail_upper_exact"]) == h_value
        assert (f_value <= ft and h_value <= ht) == (i == len(records)-1)
    assert len(case["cache_site_domains"]) == 12
    for iq, site in enumerate(case["cache_site_domains"]):
        cache_site(context, box, rule, ids, free, global_indices, site, face, edge, iq)
    return len(records), width


def main():
    elapsed = time.perf_counter()
    checks = {}
    def check(key, value):
        assert value is True, key
        checks[key] = True
    files = {P.SOURCE: SOURCE, OUT/"protocol.json": PROTOCOL, OUT/"preflight.json": PREFLIGHT,
             OUT/"domains.json": DOMAINS, OUT/"body_geometry.json": BODY,
             OUT/"domains_prelaunch_independent_review.json": PREVIOUS_GATE}
    check("source_protocol_preflight_domain_body_and_prior_gate_hashes", all(P.sha(f)==v for f,v in files.items()))
    P.gate(OUT, "domains")
    protocol, domain, body, started, terminal = [P.read(OUT/f) for f in (
        "protocol.json", "domains.json", "body_geometry.json", "domains_started.json", "domains_terminal.json")]
    P.assert_inputs(protocol)
    check("both_previous_generic_reviews_current_bound_and_pass",
          all(protocol["source_and_input_hashes"][p] == h == P.sha(ROOT/p) and
              P.read(ROOT/p)["status"] == "PASS" for p, h in GENERIC.items()))
    check("actual_domains_pid_exit_and_current_result_binding",
          terminal["status"] == "FINISHED" and started["pid"] == terminal["pid"] == 35096 and
          terminal["exit_code"] == 0 and terminal["results_sha256"] == DOMAINS and
          domain["status"] == "PASS_FRESH_DOMAINS_AND_SELECTORS_UNREVIEWED" and
          domain["actual_force_H_integrals_run"] is False and domain["protocol_sha256"] == PROTOCOL and
          domain["record_id"] == protocol["record_id"] == RECORD)
    check("no_integration_started_before_independent_gate",
          not any((OUT/f).exists() for f in ("bulk_started.json", "surface_started.json",
                                            "integration_prelaunch_independent_review.json")))
    box = P.read(ROOT / protocol["coordinate_box_file"])
    check("unchanged_original_nodal_box_bound_to_actual_domain_and_geometry",
          P.sha(ROOT / protocol["coordinate_box_file"]) == domain["coordinate_box_sha256"] ==
          protocol["coordinate_box_sha256"] == prior.PRIOR_BOX_SHA and domain["body_sha256"] == BODY and
          box["root_membership_assumed"] is False)
    P.precision.prec = 640
    P.core.A.cache_clear()
    context = P.C.load_context(RECORD)
    nodes = P.S.common_nodes(context.target, box)
    rule = context.target.auxiliary_rule()
    check("exact_context_row_and_unchanged_original_corrector_radius",
          context.row == protocol["row"] and context.old_radius == q(protocol["normal_radius_exact"]))
    geometry = P.G.certify_box(box["reference_nodes_exact"], box["current_centers_exact"],
                               box["triangles"], box["coordinate_radii_exact"])
    check("fresh_reviewed_body_geometry_kernel_complete_exact_replay",
          body == geometry and body["status"] == "CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING" and
          body["axis_preserved_for_entire_box"] is True and body["global_injectivity_for_entire_box"] is True and
          all(e["positive_J_for_entire_box"] is True for e in body["triangles"]) and len(body["triangles"]) == 144)
    print("Reviewed body replay matches", flush=True)
    ft, ht = q(protocol["force_tail_budget_exact"]), q(protocol["H_entry_tail_budget_exact"])
    count = Counter()
    attempts_bulk = 0
    assert len(domain["bulk_selectors"]) == 144
    for i, b in enumerate(domain["bulk_selectors"]):
        attempts, chart = bulk_check(context, nodes, rule, b, body["triangles"][i], i, ft, ht)
        attempts_bulk += attempts
        count[chart] += 1
        if (i+1) % 24 == 0:
            print(f"Reviewed {i+1}/144 bulk selectors", flush=True)
    check("all_144_bulk_element_context_topology_qstar_positive_weights_and_derivative_norms",
          dict(count) == {"off_axis":132, "axis_vertex":6, "axis_edge_physical":6})
    check("all_bulk_J_to_rho_and_exact_two_mass_force_Hessian_order_attempts", attempts_bulk > 144)
    attempts_surface = leaves = 0
    assert len(domain["surface_cases"]) == 24
    for i, case in enumerate(domain["surface_cases"]):
        attempts, width = surface_check(context, nodes, box, rule, case, i//12, i%12, ft, ht)
        attempts_surface += attempts
        leaves += width
    check("all_24_closed_G_domains_and_bins_raw_to_outward_dyadic_rho", leaves >= 24)
    check("all_surface_two_mass_force_Hessian_order_attempts", attempts_surface >= 24)
    check("all_288_same_nodal_box_unweighted_cache_strain_line_domains_and_provenance",
          sum(len(row["cache_site_domains"]) for row in domain["surface_cases"]) == 288)
    print("Reviewed 24 surface selectors and 288 cache sites", flush=True)
    fake = {"status":"PASS", "checks":{"fixture":True}, "source_sha256":SOURCE,
            "protocol_sha256":PROTOCOL, "preflight_sha256":PREFLIGHT, "domains_result_sha256":"0"*64}
    with patch.object(P, "read", return_value=fake):
        try:
            P.gate(OUT, "integration")
        except ValueError as error:
            assert "geometry/selectors" in str(error)
        else:
            raise AssertionError("Incorrect domains hash was accepted")
    check("integration_guard_rejects_wrong_actual_domains_hash", True)
    check("after_review_current_digests_still_match", all(P.sha(f)==v for f,v in files.items()))
    support = {"record_id": RECORD, "hashes": {str(f.relative_to(ROOT)): h for f,h in files.items()},
               "phase_results": "Only geometry, local field polynomials and selector scalar bounds evaluated.",
               "actual_domain_run_pid": terminal["pid"], "actual_domain_elapsed_seconds": domain["seconds"],
               "review_elapsed_seconds": time.perf_counter()-elapsed, "body_J": body["J_bound"],
               "bulk_charts": dict(count), "bulk_selectors": 144,
               "bulk_exact_scalar_order_attempts": attempts_bulk, "surface_edges": 24,
               "surface_closed_interval_leaves": leaves, "surface_exact_scalar_order_attempts": attempts_surface,
               "cache_same_box_sites": 288, "engine_reuse": "Reviewed body certification, context/Arb hull construction, "
               "geometry Bernstein conversion, low-degree polynomial jets and norm engines reused. Cache maps and "
               "interval arithmetic, affine geometry masses, finite-rule bindings, dyadic radii and tail arithmetic "
               "reconstructed in this review."}
    support_file = HERE / "actual_L_integration_prelaunch_checks.json"
    prior.prior_review.write_new(support_file, support)
    gate = {"status":"PASS", "checks":checks, "source_sha256":SOURCE, "protocol_sha256":PROTOCOL,
            "preflight_sha256":PREFLIGHT, "domains_result_sha256":DOMAINS, "body_sha256":BODY,
            "record_id":RECORD, "support_file":support_file.relative_to(ROOT).as_posix(),
            "support_sha256":P.sha(support_file), "review_source_sha256":P.sha(__file__),
            "scope":"Prelaunch of L bulk and surface components on this original box only. Review performed no force "
                     "or Hessian integration, moment defects, cache derivative sums, global assembly or root proof."}
    gate_file = OUT / "integration_prelaunch_independent_review.json"
    prior.prior_review.write_new(gate_file, gate)
    P.gate(OUT, "integration")
    print(json.dumps({"status":"PASS", "checks":len(checks), "gate_sha256":P.sha(gate_file),
                      "support_sha256":P.sha(support_file), "bulk_attempts":attempts_bulk,
                      "surface_attempts":attempts_surface, "closed_surface_leaves":leaves}, indent=2), flush=True)


if __name__ == "__main__":
    main()
