"""Record-neutral adapter for scaled off-axis and reviewed axis bulk components."""
from fractions import Fraction as Q
from pathlib import Path
import sys
from time import monotonic

from flint import arb

import exact_cohort_integrator as E
import review_cohort_exact_record as D


ROOT = Path(__file__).resolve().parents[1]
SCALED_KERNEL = (
    "runs/exact_integration_cohort_v1/L_defect_diagnostic/scaled_chart_v1/scaled_element21.py"
)
sys.path.insert(0, str((ROOT / SCALED_KERNEL).parent))
import scaled_element21 as S


U, M, T, JET, core, A = E.U, E.M, E.T, E.JET, E.core, E.U.A


def _setup(context, nodes, rule, body_element, selector, element):
    D.require(type(element) is int and 0 <= element < len(context.target.a["triangles"])
              and selector["element"] == element, "exact explicit element selector required")
    setup = E._bulk_setup(context, nodes, rule, element)
    E._validate_selector(selector, setup)
    identity, fields, drive, _, mass, _, _ = setup
    lo, hi, raw, rounded = E._body_rho(body_element)
    D.require(list(map(Q, selector["J_bounds_exact"])) == [lo, hi]
              and Q(selector["raw_rho_exact"]) == raw and Q(selector["rho_exact"]) == rounded,
              "same independently reviewed per-element J bounds and outward rho required")
    ids = list(map(int, context.target.a["triangles"][element]))
    return identity, fields, drive, mass, ids


def _tail_bounds(selector, mass, g_scaled, h_scaled):
    g0, h0 = (Q(selector[key]) for key in ("max_DJ_norm_upper_exact", "max_D2J_norm_upper_exact"))
    first, second = (Q(selector[key]) for key in ("phi_prime_tail_upper_exact", "phi_second_tail_upper_exact"))
    g = [min(g0, value) for value in g_scaled]
    h = {key: min(h0, value) for key, value in h_scaled.items()}
    f = [mass * first * value for value in g]
    H = {(i, j): mass * (second * g[i] * g[j] + first * value) for (i, j), value in h.items()}
    D.require(max(f) <= Q(selector["force_tail_budget_exact"])
              and max(H.values()) <= Q(selector["hessian_tail_budget_exact"]),
              "scaled element must retain original force and Hessian tail budgets")
    return f, H


def _proof(context, nodes, rule, body_element, selector, setup):
    identity, fields, _, mass, ids = setup
    route = fields["chart"]["kind"]
    proof = {
        **identity, "body_element_sha256": D._digest(body_element),
        "original_selector_sha256": D._digest(selector), "local_nodes": ids,
        "unchanged_degree": selector["degree"], "unchanged_rho_exact": selector["rho_exact"],
        "original_force_tail_budget_exact": selector["force_tail_budget_exact"],
        "original_H_entry_tail_budget_exact": selector["hessian_tail_budget_exact"],
    }
    if route != "off_axis":
        D.require(route in ("axis_vertex", "axis_edge_physical"), "explicit reviewed axis route required")
        return {**proof, "integration_chart": route, "scaled_polynomial_identity_count": 0}, None

    geometry = S.chart(fields["corners"])
    center = context.target.nodes()
    free = fields["free_local"]
    identities = S.rational_identities(geometry, [center[j] for j in ids], free)
    for k in set(range(12)) - set(free):
        node, axis = ids[k // 2], k % 2
        c = center[node][axis]
        D.require(E.SUR.rational_hull(nodes[node][axis]) == (c, c),
                  "fixed coordinate must remain at the stored exact center value")
    locations = [S.chart_point(geometry, p) for p in rule["volume_barycentric"]]
    D.require(len(locations) == 64 and all(abs(x) <= 1 and abs(y) <= 1 and r > 0 for x, y, r in locations),
              "same strictly off-axis 64-node auxiliary rule required")
    q, B = S.numerators(geometry, [nodes[i] for i in ids], free)
    jet = JET.bulk_jacobian_jet(q, B)
    g = [S._norm(p, geometry) for p in jet["gradient"]]
    h = {key: S._norm(p, geometry) for key, p in jet["hessian_upper"].items()}
    f, H = _tail_bounds(selector, mass, g, h)
    return {
        **proof, "integration_chart": "scaled_positive_R_affine_material",
        "scaled_chart_geometry": S._metadata(geometry), "scaled_polynomial_identity_count": identities,
        "scaled_exact_Qstar_node_count": len(locations), "scaled_Qstar_node_sha256": D._digest(locations),
        "scaled_force_tail_max_exact": max(f), "scaled_H_entry_tail_max_exact": max(H.values()),
    }, geometry


def prepare_component(context, nodes, rule, body_element, selector, element):
    """Select and prove one route on an already reviewed original domain without moments."""
    setup = _setup(context, nodes, rule, body_element, selector, element)
    return _proof(context, nodes, rule, body_element, selector, setup)[0]


def _scaled_component(context, nodes, rule, element, selector, setup, geometry):
    identity, fields, drive, mass, ids = setup
    fun = S.Functional(geometry, context.target.energy_scale, rule["volume_barycentric"], rule["volume_weights"][element])
    try:
        exact_mass = geometry.det_physical * sum(r for r, _ in geometry.corners) / 6
        D.require(fun.integral(0, 0, 0).overlaps(context.target.normalization() * A(exact_mass))
                  and fun.weight_mass > 0, "same original positive material and auxiliary measures required")
        f, H, g, h = S._eval_parts(geometry, [nodes[i] for i in ids], fields["free_local"], fun, drive, selector["degree"])
        Ft, Ht = _tail_bounds(selector, mass, g, h)
        cf, ch, cache = U.uniform_cache(U.RuleTarget(context.target, rule), fields, element, nodes, drive)
        D.require(cache["all64_convex_strain_lines_positive"] is True,
                  "all same-Y frozen cache paths must stay in the analytic domain")

        def inflate(value, tail):
            return core.pack(value + arb(0, core.upper(A(tail))))

        return {
            **identity, "status": "UNIFORM_COHORT_BULK_COMPONENT_UNREVIEWED",
            "degree": selector["degree"], "rho_exact": selector["rho_exact"],
            "force_integral_minus_cache": [inflate(value + cf[i], Ft[i]) for i, value in enumerate(f)],
            "H_integral_minus_cache_upper": [
                {"i": i, "j": j, "enclosure": inflate(value + ch[i, j], Ht[i, j]),
                 "two_mass_tail_upper": core.pack(A(Ht[i, j])), "two_mass_tail_upper_exact": str(Ht[i, j])}
                for (i, j), value in H.items()
            ],
            "force_two_mass_tail_upper": [core.pack(A(value)) for value in Ft],
            "force_two_mass_tail_upper_exact": list(map(str, Ft)),
            "cache_consistency": cache, "unique_moment_defects": fun.defect.cache_info().misses,
            "scope": "One bulk element on this explicit original nodal box and closed cell.",
        }
    finally:
        fun.integral.cache_clear()
        fun.site_power.cache_clear()
        fun.defect.cache_clear()
        M._triangle_moment.cache_clear()
        core.A.cache_clear()


def bulk_component(context, nodes, rule, body_element, selector, element, proof):
    """Compute one component after a caller's independent domain and proof review."""
    start = monotonic()
    setup = _setup(context, nodes, rule, body_element, selector, element)
    current, geometry = _proof(context, nodes, rule, body_element, selector, setup)
    D.require(D._digest(current) == D._digest(proof), "changed record-specific route preflight or nodal proof")
    if geometry is None:
        result = E.bulk_component(context, nodes, rule, selector)
    else:
        result = _scaled_component(context, nodes, rule, element, selector, setup, geometry)
    identity = setup[0]
    D.require(all(result[k] == value for k, value in identity.items()),
              "computed packet differs from current record, box, rule or map")
    return {
        **result, "integration_chart": proof["integration_chart"],
        "record_component_proof_sha256": D._digest(proof), "seconds": monotonic() - start,
    }
