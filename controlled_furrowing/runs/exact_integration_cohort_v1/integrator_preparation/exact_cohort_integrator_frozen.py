"""Explicit-record uniform exact-integral/cache components for frozen cohort domains.

The caller binds the common nodal box to a reviewed body/trace certificate and
sets Arb precision to 640, clearing core.A's precision-dependent cache first.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import ceil
from time import monotonic

from flint import arb, ctx as arb_context

import exact_axis_edge_physical as AX
import exact_uniform_axis_edge_bridge_B_safe as AXbridge
import exact_uniform_bulk_defect_B_safe as U
import exact_uniform_surface_bridge_B_safe as SUR


P, T, M, JET, core = U.P, U.T, U.M, U.JET, U.core
PRECISION_BITS = 640
RHO_DENOMINATOR = 2**48
MAX_BULK_ORDER = 32


class TailBudgetUnresolved(ValueError):
    """Preserve every attempted order when the declared budget is unresolved."""

    def __init__(self, element, attempts):
        super().__init__(f"bulk element {element} exceeded both-budget order policy through {attempts[-1]['degree']}")
        self.attempts = attempts


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _rational(value):
    try:
        return T.rational(value)
    except (TypeError, ValueError, ZeroDivisionError, OverflowError) as error:
        raise ValueError("finite non-Boolean rational required") from error


def _digest(value):
    return sha256(json.dumps(value, separators=(",", ":"), sort_keys=True).encode()).hexdigest()


def _context_identity(ctx):
    cell = {key: [str(_rational(value)) for value in ctx.cell[key]] for key in T.ORDER}
    row_cell = {key: [str(_rational(value)) for value in ctx.row["parameter_cell"][key]] for key in T.ORDER}
    _require(arb_context.prec == PRECISION_BITS, "cohort components require caller-set 640-bit precision")
    _require(ctx.record_id == ctx.row["record_id"] and ctx.dose == _rational(ctx.row["dose_exact"]) > 0
             and cell == row_cell and _digest(cell) == ctx.row["parameter_cell_sha256"],
             "explicit record, exact dose and unchanged cell required")
    return {
        "record_id": ctx.record_id,
        "dose_exact": str(ctx.dose),
        "parameter_cell_sha256": ctx.row["parameter_cell_sha256"],
        "model_observation_contract_sha256": ctx.row["model_observation_contract_sha256"],
        "point_cache_sha256": ctx.row["point_cache_sha256"],
        "predictor_sha256": ctx.row["predictor_sha256"],
        "precision_bits": PRECISION_BITS,
    }


def _bulk_fields(target, nodes, index):
    _require(type(index) is int and 0 <= index < len(target.a["triangles"]),
             "valid explicit integer bulk element required")
    ids = list(map(int, target.a["triangles"][index]))
    corners = target.a["reference"][ids[:3]]
    axis_count = sum(T.rational(row[0]) == 0 for row in corners)
    chart = {0: "off_axis", 1: "axis_vertex", 2: "axis_edge_physical"}.get(axis_count)
    _require(chart is not None, "triangle needs positive material radial extent")
    fields = AX.physical_element(target, index, nodes=nodes) if axis_count == 2 else target.element(index, nodes=nodes)
    _require(fields["chart"]["kind"] == chart and fields["global_indices"]
             and len(set(fields["global_indices"])) == len(fields["global_indices"])
             and len(fields["global_indices"]) == len(fields["free_local"]),
             "fresh element chart and unique nonempty free projection required")
    return fields, ids


def _norm(poly, fields):
    chart = fields["chart"]["kind"]
    _require(all(type(i) is int and type(j) is int and j >= 0 and (i >= 0 or chart == "off_axis")
                 for i, j in poly), "axis and axial polynomial exponents must be nonnegative integers")
    return AXbridge.pnorm(poly, fields) if chart == "axis_edge_physical" else U.pnorm(poly, fields)


def _functional_type(fields):
    return AXbridge.PhysicalFunctional if fields["chart"]["kind"] == "axis_edge_physical" else P.Functional


def _bulk_rule(target, rule, index):
    """Bind the explicit local rule to the authoritative normalized cache weights and sites."""
    barycentric = rule["volume_barycentric"]
    weights = rule["volume_weights"][index]
    _require(len(barycentric) == len(weights) == 64, "all 64 positive frozen bulk sites required")
    coordinates = []
    for cached, site in zip(target.a["shape"], barycentric):
        N = list(map(_rational, cached))
        xi, eta = N[1] + (N[3] + N[4]) / 2, N[2] + (N[4] + N[5]) / 2
        exact = (1 - xi - eta, xi, eta)
        _require(tuple(site) == exact and min(exact) > 0, "changed or noninterior auxiliary bulk site")
        coordinates.append(list(map(str, exact)))
    exact_weights = [_rational(w) / target.energy_scale for w in target.a["volume_weights"][index]]
    _require(list(weights) == exact_weights and min(exact_weights) > 0,
             "changed or nonpositive normalized auxiliary bulk weight")
    return _digest([coordinates, list(map(str, exact_weights))])


def _geometry_mass_upper(target, fields, index, rule):
    integral_r = fields["det"] * sum(row[0] for row in fields["corners"]) / 6
    auxiliary = sum(rule["volume_weights"][index], Fraction(0))
    _require(integral_r > 0 and auxiliary > 0, "positive material and auxiliary masses required")
    mass = core.upper(target.normalization() * U.A(integral_r) + U.A(auxiliary))
    return Fraction(str(mass.fmpq()))


def _local_node_digest(nodes, ids):
    return _digest([[list(map(str, SUR.rational_hull(value))) for value in nodes[node]] for node in ids])


def _body_rho(body_element):
    _require(body_element.get("positive_J_for_entire_box") is True
             and isinstance(body_element.get("J_bound"), dict),
             "fresh whole-box positive-J element certificate required")
    lower, upper = (_rational(body_element["J_bound"][key]) for key in ("lower_exact", "upper_exact"))
    _require(Fraction(1, 10**6) < lower <= upper, "bulk J bounds must stay above the analytic cutoff")
    raw = max(abs(lower - 1), abs(upper - 1))
    rho = Fraction(ceil(raw * RHO_DENOMINATOR), RHO_DENOMINATOR)
    _require(raw <= rho < 1, "bulk scalar expansion must converge after outward rho rounding")
    return lower, upper, raw, rho


def _bulk_materials(ctx, fields):
    parameters = ctx.target.packet["parameters"]
    _require(_rational(parameters["mu"]) == 1 and _rational(parameters["bulk"]) == 5000,
             "reviewed uniform bulk engine requires frozen mu=1, kappa=5000")
    active = fields["active"]
    lo, hi = map(_rational, ctx.cell["common_gain"])
    _require(active in (0, 1) and 0 < lo <= hi, "frozen element mask and ordered positive gain required")
    return U.A(active * ctx.dose * lo).union(U.A(active * ctx.dose * hi))


def _bulk_setup(ctx, nodes, rule, index):
    identity = _context_identity(ctx)
    fields, ids = _bulk_fields(ctx.target, nodes, index)
    rule_hash = _bulk_rule(ctx.target, rule, index)
    drive = _bulk_materials(ctx, fields)
    jet = JET.bulk_jacobian_jet(fields["q"], fields["B"])
    mass = _geometry_mass_upper(ctx.target, fields, index, rule)
    gradient = [_norm(poly, fields) for poly in jet["gradient"]]
    hessian = {key: _norm(poly, fields) for key, poly in jet["hessian_upper"].items()}
    identity.update({
        "element": index,
        "chart": fields["chart"]["kind"],
        "global_indices": fields["global_indices"],
        "free_local": fields["free_local"],
        "active": str(fields["active"]),
        "nodal_hulls_sha256": _local_node_digest(nodes, ids),
        "auxiliary_element_rule_sha256": rule_hash,
    })
    return identity, fields, drive, jet, mass, gradient, hessian


def _select(identity, body_element, mass, gradient, hessian, force_budget, hessian_budget, max_order):
    _require(type(max_order) is int and 2 <= max_order <= MAX_BULK_ORDER,
             "maximum bulk order must be an integer in 2..32")
    force_budget, hessian_budget = map(_rational, (force_budget, hessian_budget))
    _require(force_budget > 0 and hessian_budget > 0, "positive force and Hessian tail budgets required")
    lower, upper, raw, rho = _body_rho(body_element)
    g, h = max(gradient), max(hessian.values())
    attempts = []
    for order in range(2, max_order + 1):
        first, second = M.volumetric_derivative_tail_bounds(order, rho, 1, 5000)
        force = mass * first * g
        H = mass * (second * g**2 + first * h)
        attempts.append({
            "degree": order,
            "phi_prime_tail_upper_exact": str(first),
            "phi_second_tail_upper_exact": str(second),
            "force_tail_upper_exact": str(force),
            "H_entry_tail_upper_exact": str(H),
        })
        if force <= force_budget and H <= hessian_budget:
            return {
                **identity,
                **attempts[-1],
                "status": "SELECTED_UNIFORM_BULK_ORDER",
                "J_bounds_exact": [str(lower), str(upper)],
                "raw_rho_exact": str(raw),
                "rho_exact": str(rho),
                "force_tail_budget_exact": str(force_budget),
                "hessian_tail_budget_exact": str(hessian_budget),
                "geometry_two_mass_upper_exact": str(mass),
                "max_DJ_norm_upper_exact": str(g),
                "max_D2J_norm_upper_exact": str(h),
                "max_order": max_order,
                "order_attempts": attempts,
            }
    raise TailBudgetUnresolved(identity["element"], attempts)


def prepare_bulk_selector(ctx, nodes, rule, body_element, index, force_tail_budget, hessian_tail_budget, max_order=32):
    """Choose a bounded order using positive geometry masses and exact scalar-tail majorants."""
    identity, _, _, _, mass, gradient, hessian = _bulk_setup(ctx, nodes, rule, index)
    return _select(identity, body_element, mass, gradient, hessian,
                   force_tail_budget, hessian_tail_budget, max_order)


def _validate_selector(selector, setup):
    identity, _, _, _, mass, gradient, hessian = setup
    _require(type(selector.get("degree")) is int and type(selector.get("element")) is int
             and selector.get("status") == "SELECTED_UNIFORM_BULK_ORDER"
             and all(selector.get(key) == value for key, value in identity.items()),
             "selected record, element, chart, free projection, nodes or rule changed")
    lower, upper = selector["J_bounds_exact"]
    body_element = {"positive_J_for_entire_box": True, "J_bound": {"lower_exact": lower, "upper_exact": upper}}
    expected = _select(identity, body_element, mass, gradient, hessian,
                       selector["force_tail_budget_exact"], selector["hessian_tail_budget_exact"],
                       selector["max_order"])
    _require(selector == expected, "bulk selector tail or selector metadata changed")


def bulk_component(ctx, nodes, rule, selector):
    """Compute one reviewed-domain component with the shared safe uniform kernels."""
    started = monotonic()
    index = selector["element"]
    setup = _bulk_setup(ctx, nodes, rule, index)
    _validate_selector(selector, setup)
    identity, fields, drive, jet, mass, norms, hessian_norms = setup
    degree = selector["degree"]
    first = Fraction(selector["phi_prime_tail_upper_exact"])
    second = Fraction(selector["phi_second_tail_upper_exact"])
    force_tails = [mass * first * value for value in norms]
    H_tails = {(i, j): mass * (second * norms[i] * norms[j] + first * value)
               for (i, j), value in hessian_norms.items()}
    mu = U.A(1)
    diagonal = [mu, mu + drive, mu, mu + drive, mu]

    def shear(left, right):
        polynomial = {}
        for d, a, b in zip(diagonal, left, right):
            polynomial = T.add(polynomial, T.scale(T.mul(a, b), d))
        return polynomial

    B = list(zip(*fields["B"]))
    try:
        functional = _functional_type(fields)(ctx.target, fields, index, rule)
        poly_a, poly_b = U.series(T.sub(jet["value"], T.ONE), degree)
        defect = JET.contract_volumetric_defect(jet, poly_a, poly_b, functional.defect)
        cache_g, cache_H, cache_info = U.uniform_cache(U.RuleTarget(ctx.target, rule), fields, index, nodes, drive)
        _require(cache_info["all64_convex_strain_lines_positive"] is True,
                 "all same-node cache strain segments must be positive")

        def enclose(value, tail):
            return core.pack(value + arb(0, core.upper(U.A(tail))))

        force = [
            enclose(defect["gradient"][i] + functional.apply(shear(fields["q"], column)) + cache_g[i],
                    force_tails[i])
            for i, column in enumerate(B)
        ]
        H = [
            {
                "i": i, "j": j,
                "enclosure": enclose(value + functional.apply(shear(B[i], B[j])) + cache_H[i, j], H_tails[i, j]),
                "two_mass_tail_upper": core.pack(U.A(H_tails[i, j])),
            }
            for (i, j), value in defect["hessian_upper"].items()
        ]
        return {
            **identity,
            "status": "UNIFORM_COHORT_BULK_COMPONENT_UNREVIEWED",
            "degree": degree,
            "rho_exact": selector["rho_exact"],
            "force_integral_minus_cache": force,
            "H_integral_minus_cache_upper": H,
            "force_two_mass_tail_upper": [core.pack(U.A(value)) for value in force_tails],
            "cache_consistency": cache_info,
            "unique_moment_defects": len(functional.requested_defect_keys),
            "seconds": monotonic() - started,
            "scope": "One bulk element on the reviewed explicit whole nodal box and unchanged physical cell.",
        }
    finally:
        # The existing class-level caches otherwise retain each processed functional.
        for method in (P.Functional.integral, P.Functional.defect, P.Functional.node_power,
                       AXbridge.PhysicalFunctional.integral):
            method.cache_clear()


def _surface_rule(target, rule, face, edge):
    a = target.a
    cache_sites = a[f"surface_{face}_N"]
    sites = [_rational(row[2]) + _rational(row[1]) / 2 for row in cache_sites]
    weights = [_rational(w) / target.energy_scale for w in a[f"surface_{face}_weights"][edge]]
    _require(len(sites) == len(weights) == 12 and all(0 < s < 1 for s in sites)
             and min(weights) > 0 and list(rule["surface_s"][face]) == sites
             and list(rule["surface_weights"][face][edge]) == weights,
             "all twelve unchanged interior surface sites and positive weights required")


def surface_component(ctx, nodes, rule, case):
    """Use the existing uniform trace kernel on the explicitly supplied record and case."""
    identity = _context_identity(ctx)
    face, edge, order = (case[key] for key in ("face", "edge", "order"))
    _require(type(face) is int and face in (0, 1) and type(edge) is int and 0 <= edge < 12
             and type(order) is int and order >= 2, "valid integer surface face, edge and order required")
    if "record_id" in case:
        _require(case["record_id"] == ctx.record_id, "surface case belongs to a different record")
    _surface_rule(ctx.target, rule, face, edge)
    try:
        result = SUR.run_edge(ctx.target, nodes, ctx.row, case, rule, order)
        return {
            **result,
            **identity,
            "status": "UNIFORM_COHORT_SURFACE_COMPONENT_UNREVIEWED",
            "scope": "One surface edge on the reviewed explicit whole nodal box and unchanged physical cell.",
        }
    finally:
        SUR.S.SurfaceFunctional.integral.cache_clear()
        SUR.S.SurfaceFunctional.defect.cache_clear()
