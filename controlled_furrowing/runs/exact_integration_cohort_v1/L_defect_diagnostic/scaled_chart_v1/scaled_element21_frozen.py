"""Exact scaled material coordinates for one frozen off-axis L force/H component."""
import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
from functools import lru_cache
from fractions import Fraction as Q
import json
from math import comb
import os
from pathlib import Path
import sys
from time import monotonic
from types import SimpleNamespace
from unittest.mock import patch

from flint import arb


OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parent))
import element21_force as F


D, E, U, M, T, JET, core, A = F.D, F.E, F.U, F.M, F.T, F.JET, F.core, F.A
ROOT = F.ROOT
PARENT = OUT.parent.relative_to(ROOT).as_posix()
PREMISES = {
    f"{PARENT}/element21_force.py": "9a5f93bb3d4a39020e28bee165ba606d1d17f63e44a17ee758959aa05109d329",
    f"{PARENT}/actual_independent_review.json": "148903b4f29c87d27ce88c8a4c81e3bbbedcf2f55a5ad37efd1ecf0dadd32a4d",
    f"{PARENT}/results.json": "e1f0a45dc3e3a7554e5fc84b60a2e670a21169d282da8e6de5f21d52e9eecdf0",
    f"{PARENT}/terminal.json": "2101a60986db57577906fa60edd6f4ed62df8261bb6307f3b3ac2b5c117692e2",
    f"{PARENT}/conditioning_witness.json": "6d93f21f25096c6da41f8db8d03f75ca1ff7f506decbb365cabf225ed27e9dd8",
}


@dataclass(frozen=True)
class Chart:
    corners: tuple
    R_center: Q
    R_scale: Q
    Z_center: Q
    Z_scale: Q
    Rhat: dict
    Rhat_min: Q
    det_physical: Q
    det_scaled: Q
    local_corners: tuple
    shape: tuple
    dR: tuple
    dZ: tuple

    def to_local(self, R, Z):
        return (R - self.R_center) / self.R_scale, (Z - self.Z_center) / self.Z_scale


def chart(corners):
    material = T.material_basis(corners)
    corners = tuple(map(tuple, material["corners"]))
    low, high = min(v[0] for v in corners), max(v[0] for v in corners)
    zc = material["zcenter"]
    rc, h = (low + high) / 2, (high - low) / 2
    k = max(abs(v[1] - zc) for v in corners)
    if not 0 < low < high or k <= 0:
        raise ValueError("strictly off-axis nondegenerate scaled material chart required")
    physical_R = T.add(T.constant(rc), T.scale(T.X, h))
    physical_zeta = T.scale(T.Y, k)
    shape = tuple(T.compose(p, physical_R, physical_zeta) for p in material["phi"])
    dR = tuple(T.scale(T.diff(p, 0), 1 / h) for p in shape)
    dZ = tuple(T.scale(T.diff(p, 1), 1 / k) for p in shape)
    local = tuple(((R - rc) / h, (Z - zc) / k) for R, Z in corners)
    a, b = (local[1][i] - local[0][i] for i in (0, 1))
    c, d = (local[2][i] - local[0][i] for i in (0, 1))
    det_scaled = a * d - b * c
    if h * k * det_scaled != material["det"] or any(abs(v) > 1 for row in local for v in row):
        raise ValueError("exact scaled-chart measure or containing rectangle mismatch")
    return Chart(corners, rc, h, zc, k, T.scale(physical_R, 1 / rc), low / rc,
                 material["det"], det_scaled, local, shape, dR, dZ)


def numerators(geometry, current, free_local):
    """The first four strain rows are polynomial; hoop rows have common Rhat denominator."""
    if len(current) != 6 or any(len(row) != 2 for row in current):
        raise ValueError("six P2 nodes with two coordinates required")
    free = list(free_local)
    if not free or len(set(free)) != len(free) or any(type(i) is not int or not 0 <= i < 12 for i in free):
        raise ValueError("unique nonempty exact local free map required")
    radial, axial = zip(*current)
    hoop_basis = tuple(T.scale(p, 1 / geometry.R_center) for p in geometry.shape)
    q = [
        T.linear_combination(geometry.dR, radial), T.linear_combination(geometry.dZ, radial),
        T.linear_combination(geometry.dR, axial), T.linear_combination(geometry.dZ, axial),
        T.linear_combination(hoop_basis, radial),
    ]
    columns = []
    for i in free:
        j, component = divmod(i, 2)
        R, Z, hoop = geometry.dR[j], geometry.dZ[j], hoop_basis[j]
        columns.append((R, Z, {}, {}, hoop) if component == 0 else ({}, {}, R, Z, {}))
    return q, [list(row) for row in zip(*columns)]


def chart_point(geometry, barycentric):
    if len(barycentric) != 3 or sum(barycentric) != 1 or min(barycentric) <= 0:
        raise ValueError("strictly interior exact three-component barycentric point required")
    R, Z = (sum(lam * v[k] for lam, v in zip(barycentric, geometry.corners)) for k in (0, 1))
    x, y = geometry.to_local(R, Z)
    Rhat = R / geometry.R_center
    if T.evaluate(geometry.Rhat, x, y) != Rhat or Rhat < geometry.Rhat_min:
        raise ValueError("unchanged physical auxiliary radius required")
    return x, y, Rhat


def material_moment(geometry, a, b, d):
    """Unnormalized integral of R*x^a*y^b/Rhat^d in the original physical triangle."""
    if any(type(i) is not int or i < 0 for i in (a, b, d)):
        raise ValueError("ordinary scaled powers and nonnegative denominator degree required")
    # The integration primitive is in physical dR dZ, so there is no
    # additional h*k factor. Its R exponent includes the radial measure once.
    expr = M.LogMoment()
    for j in range(a + 1):
        term = M.triangle_moment(geometry.corners, 1 + j - d, b, geometry.Z_center)
        expr += term * (comb(a, j) * (-geometry.R_center)**(a - j))
    return expr * (geometry.R_center**d / (geometry.R_scale**a * geometry.Z_scale**b))


class Functional:
    """Signed moment defect with a denominator power in the cache key."""

    def __init__(self, geometry, sigma, barycentric, weights):
        if len(weights) != len(barycentric) or not weights or min(weights) <= 0 or sigma <= 0:
            raise ValueError("positive frozen normalization and matching auxiliary quadrature required")
        self.geometry = geometry
        self.normalization = 2 * arb.pi() / A(sigma)
        self.sites = [(tuple(map(A, chart_point(geometry, p))), A(w))
                      for p, w in zip(barycentric, weights)]
        self.weight_mass = sum(weights, Q(0))
        # Instance-owned caches are released with this one element.
        self.integral = lru_cache(maxsize=None)(self._integral)
        self.site_power = lru_cache(maxsize=None)(self._site_power)
        self.defect = lru_cache(maxsize=None)(self._defect)

    def _site_power(self, index, axis, exponent):
        return self.sites[index][0][axis]**exponent

    def _integral(self, a, b, d):
        return self.normalization * material_moment(self.geometry, a, b, d).ball(640)

    def _defect(self, a, b, d):
        material = self.integral(a, b, d)
        finite = sum((w * self.site_power(i, 0, a) * self.site_power(i, 1, b) * self.site_power(i, 2, -d)
                      for i, (_, w) in enumerate(self.sites)), arb(0))
        return material - finite

    def apply(self, poly, d):
        return sum((T.cmul(value, self.defect(a, b, d)) for (a, b), value in poly.items()), arb(0))


def series_numerators(e_numerator, Rhat, order):
    """Horner numerators for phi' and phi'' over the shared exact Rhat^order."""
    first, second = M.volumetric_derivative_coefficients(order, 1, 5000)
    A_poly, B_poly, R_power = T.constant(A(first[-1])), T.constant(A(second[-1])), T.ONE
    for degree in range(order - 1, -1, -1):
        R_power = T.mul(R_power, Rhat)
        A_poly = T.add(T.mul(e_numerator, A_poly), T.scale(R_power, A(first[degree])))
        B_poly = T.add(T.mul(e_numerator, B_poly), T.scale(R_power, A(second[degree])))
    return A_poly, B_poly


def _norm(poly, geometry):
    return sum((max(map(abs, E.SUR.rational_hull(value))) for value in poly.values()), Q(0)) / geometry.Rhat_min


def _eval_parts(geometry, local, free, fun, drive, order):
    q, B = numerators(geometry, local, free)
    jet = JET.bulk_jacobian_jet(q, B)
    # JET is applied to the fifth strain numerator h_num. J, J_u, and
    # J_uu are linear in that strain or its derivatives, hence divide by Rhat.
    e_numerator = T.sub(jet["value"], geometry.Rhat)
    a, b = series_numerators(e_numerator, geometry.Rhat, order)
    grad, hess = jet["gradient"], jet["hessian_upper"]
    products = {(i, j): T.mul(grad[i], grad[j]) for i, j in hess}
    support_a = {k for poly in [*grad, *hess.values()] for k in poly}
    support_b = {k for poly in products.values() for k in poly}
    dual_a = JET.moment_dual(a, support_a, lambda i, j: fun.defect(i, j, order + 1))
    dual_b = JET.moment_dual(b, support_b, lambda i, j: fun.defect(i, j, order + 2))
    diagonal = [A(1), A(1) + drive, A(1), A(1) + drive]

    def shear(left, right):
        numerator = {}
        for coefficient, p, q in zip(diagonal, left, right):
            numerator = T.add(numerator, T.scale(T.mul(p, q), coefficient))
        return fun.apply(numerator, 0) + fun.apply(T.mul(left[4], right[4]), 2)

    columns = list(zip(*B))
    force = [JET.contract_factor(p, dual_a) + shear(q, column) for p, column in zip(grad, columns)]
    H = {
        (i, j): JET.contract_factor(hess[i, j], dual_a) + JET.contract_factor(products[i, j], dual_b)
                + shear(columns[i], columns[j])
        for i, j in hess
    }
    return force, H, [_norm(p, geometry) for p in grad], {key: _norm(p, geometry) for key, p in hess.items()}


def _certificate():
    hashes, context, case, old, rule, ids, nodes, center, setup, _, _ = F._geometry()
    support = {p: D._bound(p, digest, hashes) for p, digest in PREMISES.items() if p.endswith(".json")}
    D._bound(f"{PARENT}/element21_force.py", PREMISES[f"{PARENT}/element21_force.py"], hashes, parse=False)
    review = support[f"{PARENT}/actual_independent_review.json"]
    if (review["status"] != "PASS_ACTUAL_LOCAL_FORCE_DIAGNOSTIC_ENDPOINT_REVIEW"
            or not D._yes(review["checks"]) or review["results_sha256"] != PREMISES[f"{PARENT}/results.json"]
            or F._sha(F.__file__) != PREMISES[f"{PARENT}/element21_force.py"]):
        raise ValueError("same reviewed element-21 force diagnostic required")
    D._terminal(support[f"{PARENT}/terminal.json"], PREMISES[f"{PARENT}/results.json"])
    D._bound(review["support_file"], review["support_sha256"], hashes, parse=False)
    _, physical, drive, _, mass, _, _ = setup
    geometry = chart(physical["corners"])
    local = [nodes[j] for j in ids]
    q, B = numerators(geometry, local, physical["free_local"])
    jet = JET.bulk_jacobian_jet(q, B)
    g = [_norm(p, geometry) for p in jet["gradient"]]
    h = {key: _norm(p, geometry) for key, p in jet["hessian_upper"].items()}
    # The selector's already audited geometry mass is for the same exact
    # material/rule measures. The only changed bounds are the chart majorants.
    first, second = (Q(case[k]) for k in ("phi_prime_tail_upper_exact", "phi_second_tail_upper_exact"))
    Ft = [mass * first * norm for norm in g]
    Ht = {(i, j): mass * (second * g[i] * g[j] + first * bound) for (i, j), bound in h.items()}
    if max(Ft) > Q(case["force_tail_budget_exact"]) or max(Ht.values()) > Q(case["hessian_tail_budget_exact"]):
        raise ValueError("scaled chart fails original individual tail budgets")
    return hashes, context, case, old, rule, ids, nodes, center, physical, drive, geometry, Ft, Ht


def _to_physical(geometry, p, d):
    """Exact numerator over physical R**d for symbolic chart-identity tests."""
    x = T.scale(T.sub(T.X, T.constant(geometry.R_center)), 1 / geometry.R_scale)
    y = T.scale(T.Y, 1 / geometry.Z_scale)
    return T.mul(T.compose(p, x, y), {(-d, 0): geometry.R_center**d})


def rational_identities(geometry, point, free):
    physical = T.element_fields(geometry.corners, point, free)
    q, B = numerators(geometry, point, free)
    for i, (value, columns) in enumerate(zip(q, B)):
        d = int(i == 4)
        if _to_physical(geometry, value, d) != physical["q"][i]:
            raise AssertionError(f"exact strain chart identity {i}")
        for j, col in enumerate(columns):
            if _to_physical(geometry, col, d) != physical["B"][i][j]:
                raise AssertionError(f"exact derivative chart identity {i},{j}")
    scaled = JET.bulk_jacobian_jet(q, B)
    original = JET.bulk_jacobian_jet(physical["q"], physical["B"])
    for x, y in [(scaled["value"], original["value"]),
                 *zip(scaled["gradient"], original["gradient"]),
                 *((v, original["hessian_upper"][key]) for key, v in scaled["hessian_upper"].items())]:
        if _to_physical(geometry, x, 1) != y:
            raise AssertionError("exact J, J_u or J_uu chart identity")
    return 5 * (1 + len(free)) + 1 + len(free) + len(scaled["hessian_upper"])


def tests(geometry, point, free, barycentric):
    checks = {}
    checks["actual_frozen_exact_q_B_J_Ju_Juu_polynomial_identities"] = rational_identities(geometry, point, free) > 100
    reference = ((Q(1), Q(0)), (Q(2), Q(0)), (Q(2), Q(1)))
    toy = chart(reference)
    midpoint = tuple(tuple((reference[i][j] + reference[k][j]) / 2 for j in (0, 1))
                     for i, k in ((0, 1), (1, 2), (2, 0)))
    synthetic = [(r + z / 7, z - r / 11) for r, z in reference + midpoint]
    checks["nonidentity_synthetic_cross_chart_exact_polynomial_equality"] = rational_identities(
        toy, synthetic, list(range(12))) > 100
    checks["determinant_chart_measure_reproduces_original"] = toy.R_scale * toy.Z_scale * toy.det_scaled == toy.det_physical
    # Independent analytic triangle slices R in [1,2], Z in [0,R-1].
    checks["symbolic_zero_power_material_mass"] = material_moment(toy, 0, 0, 0) == M.LogMoment(Q(5, 6))
    analytic_log = M.LogMoment(Q(1), ((Q(2), Q(-1)),))
    checks["symbolic_inverse_square_denominator_with_log_measure"] = material_moment(toy, 0, 0, 2) == analytic_log * toy.R_center**2
    signed = (M.LogMoment(Q(1, 2)) + analytic_log * (-toy.R_center)) * (toy.R_center**2 / toy.R_scale)
    checks["symbolic_signed_centered_radial_log_moment"] = material_moment(toy, 1, 0, 2) == signed
    sigma, weights = Q(3), [Q(1, 7), Q(1, 13)]
    sites = [(Q(1, 3), Q(1, 3), Q(1, 3)), (Q(1, 2), Q(1, 4), Q(1, 4))]
    f = Functional(toy, sigma, sites, weights)
    old = E.P.Functional(SimpleNamespace(energy_scale=sigma), {**T.material_basis(reference),
                                  "chart": {"kind": "off_axis"}}, 0,
                         {"volume_barycentric": sites, "volume_weights": [weights]})
    checks["synthetic_defect_zeroth_mass_matches_physical_engine"] = f.defect(0, 0, 0).overlaps(old.defect(0, 0))
    for a, b, d in ((1, 0, 2), (2, 1, 3), (1, 2, 0), (0, 2, 1)):
        ambient = _to_physical(toy, {(a, b): Q(1)}, d)
        checks[f"synthetic_signed_defect_invariance_{a}_{b}_{d}"] = f.defect(a, b, d).overlaps(old.apply(ambient))
    for l in barycentric:
        x, y, r = chart_point(geometry, l)
        R, Z = (sum(v * c[k] for v, c in zip(l, geometry.corners)) for k in (0, 1))
        if (geometry.R_center + geometry.R_scale * x, geometry.Z_center + geometry.Z_scale * y) != (R, Z):
            raise AssertionError("actual exact auxiliary site altered")
        if geometry.R_center * r != R:
            raise AssertionError("actual exact auxiliary radius altered")
    checks["all_64_actual_sites_mapped_exactly_without_sums"] = len(barycentric) == 64
    p = {(1, 0): arb(0, A(Q(1, 8)))}
    a, b = series_numerators(p, T.ONE, 3)
    old_a, old_b = U.series(p, 3)
    traces = ((a, old_a), (b, old_b))
    checks["signed_uncertain_series_coefficient_kept_at_unit_denominator"] = all(
        (1, 0) in v and not v[1, 0].is_zero() and all(w[k].overlaps(z) for k, z in v.items())
        for v, w in traces)
    e = {(0, 0): Q(1, 7), (1, 0): Q(-2, 9)}
    a, b = series_numerators(e, toy.Rhat, 4)
    x, y = Q(2, 7), Q(0)
    denominator = T.evaluate(toy.Rhat, x, y)
    r = T.evaluate(e, x, y) / denominator
    coeffs = M.volumetric_derivative_coefficients(4, 1, 5000)
    checks["Horner_two_constitutive_derivatives_keep_exact_Rhat_powers"] = all(
        T.evaluate(v, x, y).contains(A(sum(c * r**j for j, c in enumerate(coef)) * denominator**4))
        for v, coef in zip((a, b), coeffs))
    try:
        material_moment(toy, 0, 0, -1)
        checks["invalid_denominator_degree_rejected"] = False
    except ValueError:
        checks["invalid_denominator_degree_rejected"] = True
    if not D._yes(checks):
        raise AssertionError(checks)
    return checks


def _metadata(geometry):
    return {k: getattr(geometry, k) for k in (
        "corners", "R_center", "R_scale", "Z_center", "Z_scale", "Rhat_min",
        "det_physical", "det_scaled", "local_corners")}


def prepare():
    forbidden = AssertionError("no actual cache sums during scaled-chart preparation")
    with patch.object(U, "uniform_cache", side_effect=forbidden), patch.object(U.C, "exact_bulk_maps", side_effect=forbidden):
        hashes, context, case, _, rule, ids, _, center, physical, _, geometry, Ft, Ht = _certificate()
        checks = tests(geometry, [center[i] for i in ids], physical["free_local"], rule["volume_barycentric"])
    checks["same_old_force_and_Hessian_tail_budgets"] = max(Ft) <= Q(case["force_tail_budget_exact"]) and max(
        Ht.values()) <= Q(case["hessian_tail_budget_exact"])
    if not D._yes(checks):
        raise AssertionError(checks)
    source = Path(__file__).resolve()
    names = ("scaled_element21_frozen.py", "protocol.json", "preflight.json")
    if any((OUT / x).exists() for x in names):
        raise ValueError("preserve write-once scaled-chart preparation")
    protocol = {
        "status": "FROZEN_SCALED_CHART_E21_PREPARATION_REVIEW_PENDING", "source_file": source.relative_to(ROOT).as_posix(),
        "source_sha256": F._sha(source), "input_hashes": hashes, "record_id": F.RECORD, "element": F.ELEMENT,
        "chart_geometry": _metadata(geometry), "precision_bits": 640, "global_indices": physical["global_indices"],
        "same_parameter_cell_sha256": context.row["parameter_cell_sha256"],
        "unchanged_order": case["degree"], "unchanged_rho_exact": case["rho_exact"],
        "unchanged_force_tail_budget_exact": case["force_tail_budget_exact"],
        "unchanged_H_entry_tail_budget_exact": case["hessian_tail_budget_exact"],
        "scaled_max_force_tail_exact": max(Ft), "scaled_max_H_entry_tail_exact": max(Ht.values()),
        "representation": "R=Rc+h*x and Z=Zc+k*y exactly. Construct rational scaled P2 polynomials before inserting interval Y. "
                          "Use physical derivatives by chain rule, Rhat=R/Rc and exact denominator powers.",
        "jet_and_series": "J, DJ, D2J are numerator polynomials over Rhat. Write J-1=(J_num-Rhat)/Rhat. "
                          "The two constitutive Horner numerators share Rhat^N. Force and H terms carry Rhat^(N+1) or Rhat^(N+2).",
        "measure": "Symbolically expand x^a about physical R, multiply Rc^d, integrate physical Laurent monomials "
                   "with the original R measure and 2*pi/sigma. Qstar uses exact same sites and old normalized binary weights.",
        "cache": "Retain exact same Y in original shared strain-cache correction, which still uses unchanged physical maps.",
        "run_budget": "One element, full 12 force and 78 upper Hessian entries on original interval box and C. "
                      "One fresh same-Y 64-site local cache correction. No global assembly or root assertion.",
        "actual_record_force_H_integration_run": False,
    }
    with (OUT / names[0]).open("xb") as stream:
        stream.write(source.read_bytes())
    protocol_sha = D._write(OUT / names[1], protocol)
    pre = {"status": "PASS_SCALED_CHART_PREPARATION_ONLY", "checks": checks, "source_sha256": protocol["source_sha256"],
           "protocol_sha256": protocol_sha, "actual_record_force_H_integration_run": False}
    pre_sha = D._write(OUT / names[2], pre)
    print(json.dumps({"status": pre["status"], "source_sha256": pre["source_sha256"], "protocol_sha256": protocol_sha,
                      "preflight_sha256": pre_sha, "checks": len(checks)}), flush=True)


def run(review_sha):
    relative = (OUT / "prelaunch_independent_review.json").relative_to(ROOT).as_posix()
    accepted = D._bound(relative, review_sha, {})
    expected = {"source_sha256": F._sha(__file__), "protocol_sha256": F._sha(OUT / "protocol.json"),
                "preflight_sha256": F._sha(OUT / "preflight.json")}
    if accepted["status"] != "PASS" or not D._yes(accepted["checks"]) or any(accepted[k] != v for k, v in expected.items()):
        raise ValueError("bound independent prelaunch review required")
    protocol = json.loads((OUT / "protocol.json").read_bytes())
    hashes, context, case, old, rule, ids, nodes, _, physical, drive, geometry, Ft, Ht = _certificate()
    if hashes != protocol["input_hashes"] or D._digest(_metadata(geometry)) != D._digest(protocol["chart_geometry"]):
        raise ValueError("scaled chart source inputs or coordinate map changed")
    started = monotonic()
    stamp = lambda: datetime.now(timezone.utc).isoformat()
    D._write(OUT / "started.json", {"pid": os.getpid(), "utc": stamp(), **expected, "prelaunch_sha256": review_sha})
    try:
        fun = Functional(geometry, context.target.energy_scale, rule["volume_barycentric"], rule["volume_weights"][F.ELEMENT])
        material_upper = core.upper(fun.integral(0, 0, 0) + A(fun.weight_mass))
        # The analytic geometry mass used in the force/H tails was frozen in
        # the original selector. Here overlap checks the measure implementation.
        reference_mass = geometry.det_physical * sum(r for r, _ in geometry.corners) / 6
        if (not fun.integral(0, 0, 0).overlaps(context.target.normalization() * A(reference_mass))
                or fun.weight_mass <= 0):
            raise ValueError("scaled signed moment uses different material or auxiliary mass")
        f, H, _, _ = _eval_parts(geometry, [nodes[j] for j in ids], physical["free_local"], fun, drive, case["degree"])
        cf, ch, cache = U.uniform_cache(U.RuleTarget(context.target, rule), physical, F.ELEMENT, nodes, drive)
        if cache["all64_convex_strain_lines_positive"] is not True:
            raise ValueError("these exact same-Y cache lines must be positive")
        def enclosure(value, tail):
            return value + arb(0, core.upper(A(tail)))
        force = [core.pack(enclosure(v + cf[i], Ft[i])) for i, v in enumerate(f)]
        hessian = [{"i": i, "j": j, "enclosure": core.pack(enclosure(v + ch[i, j], Ht[i, j])),
                    "tail_upper_exact": Ht[i, j]} for (i, j), v in H.items()]
        def intersection(a, b):
            x, y = D._packed(a), D._packed(b)
            lo, hi = max(x[0], y[0]), min(x[1], y[1])
            if lo > hi:
                raise ValueError("validated physical and scaled chart enclosures must intersect")
            return lo, hi
        previous = D._entry_map(old["H_integral_minus_cache_upper"], len(f))
        force_intersect = [intersection(v, p) for v, p in zip(force, old["force_integral_minus_cache"])]
        H_intersect = [intersection(v["enclosure"], previous[v["i"], v["j"]]["enclosure"]) for v in hessian]
        report = {
            "status": "SCALED_CHART_LOCAL_FORCE_H_COMPONENT_UNREVIEWED", **expected, "prelaunch_sha256": review_sha,
            "record_id": F.RECORD, "element": F.ELEMENT, "global_indices": physical["global_indices"],
            "coordinate_geometry": _metadata(geometry), "precision_bits": 640, "degree": case["degree"],
            "rho_exact": case["rho_exact"], "input_hashes": hashes,
            "force_integral_minus_cache": force, "H_integral_minus_cache_upper": hessian, "force_tail_upper_exact": Ft,
            "force_intersections_with_original_exact": force_intersect,
            "H_intersections_in_packet_order_exact": H_intersect, "positive_same_Y_cache": cache,
            "material_plus_rule_mass_upper": core.pack(material_upper), "signed_moment_defects": fun.defect.cache_info().misses,
            "new_global_assembly_or_root_test": False, "seconds": monotonic() - started,
        }
        digest = D._write(OUT / "results.json", report)
        D._write(OUT / "terminal.json", {"status": "FINISHED", "exit_code": 0, "pid": os.getpid(), "utc": stamp(),
                                        "results_sha256": digest})
        print(json.dumps({"status": report["status"], "sha256": digest}), flush=True)
    except BaseException as error:
        D._write(OUT / "failure.json", {"error": type(error).__name__, "message": str(error), **expected})
        D._write(OUT / "terminal.json", {"status": "FAILED", "exit_code": 1, "pid": os.getpid(), "utc": stamp()})
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest="command", required=True)
    command.add_parser("prepare")
    run_args = command.add_parser("run")
    run_args.add_argument("--review-sha256", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare()
    else:
        run(args.review_sha256)
