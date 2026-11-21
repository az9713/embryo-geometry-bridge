"""SAFE uniform surface bridge on the fixed common B nodal box and SAME C.

prepare: selector/domain proofs and synthetic tests only; no actual integrals
or cache sums. run: requires independent source/protocol/preflight approval.
This is a surface component, not a root or continuum certificate.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import argparse
import datetime
import hashlib
import json
import time

from flint import arb, arb_mat, ctx

import exact_surface_bridge_safe as S
import exact_surface_cache_consistency as L

T, JET, M, core = S.T, S.JET, S.M, S.core
ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(__file__)
OUT = ROOT / "runs/exact_integration_bridge_v1/uniform_surface_B_safe_v2"
TARGET_PREPARATION = ROOT / S.PREP
BOXFILE = ROOT / "runs/exact_integration_bridge_v1/retube_B_geometry/coordinate_boxes.json"
GEOMETRY_PROTOCOL = BOXFILE.with_name("protocol.json")
PRECISION = 640
ORDER = 28
MAX_DEPTH = 8
RHO_GOAL = F(1, 8)
SCOPE = "Uniform surface-only integral/cache component over SAME explicit common B nodal box and physical parameter cell C. No equilibrium/root membership, complete bulk+surface bound, PDE accuracy, biological validation or five-gate completion."


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding="utf8"))


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise ValueError("preserve write-once output: " + str(path))
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf8")


def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def A(value):
    if isinstance(value, arb):
        if not value.is_finite():
            raise ValueError("finite coefficient required")
        return value
    return core.A(T.rational(value))


def set_precision(bits):
    ctx.prec = bits
    # core.A caches conversions at the active precision.
    core.A.cache_clear()


def rational_hull(value):
    if isinstance(value, arb):
        if not value.is_finite():
            raise ValueError("finite interval coefficient required")
        return F(str(value.lower().fmpq())), F(str(value.upper().fmpq()))
    q = T.rational(value)
    return q, q


def iadd(x, y):
    return x[0] + y[0], x[1] + y[1]


def iscale(x, scalar):
    vals = x[0] * scalar, x[1] * scalar
    return min(vals), max(vals)


def ihull(values):
    return min(x[0] for x in values), max(x[1] for x in values)


def iball(bounds):
    if bounds[0] > bounds[1]:
        raise ValueError("ordered coefficient bounds required")
    return A(bounds[0]).union(A(bounds[1]))


def pack_bounds(bounds):
    return list(map(str, bounds))


def coefficient_hulls(poly):
    if any(not isinstance(i, int) or j != 0 or i < 0 for i, j in poly):
        raise ValueError("nonnegative univariate chart polynomial required")
    return {i: rational_hull(value) for (i, _), value in poly.items()}


def bernstein_bounds(poly, lo, hi):
    """Exact rational hull of every interval Bernstein coefficient on [lo,hi].

    Convert p(lo+(hi-lo)t) to power coefficients d_j, then to degree-n
    Bernstein coefficients b_i=sum_{j<=i} C(i,j)/C(n,j)*d_j. Every conversion
    acts on exact rational endpoint intervals; no Arb overlap comparison is
    used to delete or select coefficients.
    """
    lo, hi = F(lo), F(hi)
    if not 0 <= lo < hi <= 1:
        raise ValueError("closed chart subinterval in [0,1] required")
    coefficients = coefficient_hulls(poly)
    degree = max(coefficients, default=0)
    d = []
    for j in range(degree + 1):
        value = (F(0), F(0))
        for k, coefficient in coefficients.items():
            if k >= j:
                value = iadd(value, iscale(coefficient, comb(k, j) * lo**(k - j) * (hi - lo)**j))
        d.append(value)
    bernstein = []
    for i in range(degree + 1):
        value = (F(0), F(0))
        for j in range(i + 1):
            value = iadd(value, iscale(d[j], F(comb(i, j), comb(degree, j))))
        bernstein.append(value)
    return ihull(bernstein), bernstein


def uniform_G_bounds(fields, max_depth=MAX_DEPTH, rho_goal=RHO_GOAL):
    if isinstance(max_depth, bool) or not isinstance(max_depth, int) or not 0 <= max_depth <= 12 or not 0 < rho_goal < 1:
        raise ValueError("finite subdivision and convergent rho policy required")
    polys = fields["q_numerators"]
    denominator = fields["hoop_denominator"]
    history = []
    selected = None
    for depth in range(max_depth + 1):
        leaves = []
        for i in range(2**depth):
            lo, hi = F(i, 2**depth), F(i + 1, 2**depth)
            qbounds, bernstein = zip(*(bernstein_bounds(p, lo, hi) for p in polys))
            den, den_coefficients = bernstein_bounds(denominator, lo, hi)
            if den[0] <= 0:
                raise ValueError("regularized hoop denominator not positive")
            qbounds = [S.idiv(qbounds[0], den), qbounds[1], qbounds[2]]
            aa, bb, cc = map(S.isquare, qbounds)
            G = aa[0] * (bb[0] + cc[0]), aa[1] * (bb[1] + cc[1])
            leaves.append({
                "s_exact": pack_bounds((lo, hi)),
                "q_bounds_exact": [pack_bounds(q) for q in qbounds],
                "G_bounds_exact": pack_bounds(G),
                "q_numerator_Bernstein_coefficients_exact": [[pack_bounds(v) for v in row] for row in bernstein],
                "hoop_denominator_Bernstein_coefficients_exact": [pack_bounds(v) for v in den_coefficients],
            })
        lower = min(F(x["G_bounds_exact"][0]) for x in leaves)
        upper = max(F(x["G_bounds_exact"][1]) for x in leaves)
        hoop = min(F(x["q_bounds_exact"][0][0]) for x in leaves)
        q0 = (lower + upper) / 2
        rho = (upper - lower) / (upper + lower) if lower > 0 else None
        history.append({"depth": depth, "leaves": len(leaves), "G_bounds_exact": pack_bounds((lower, upper)), "hoop_lower_exact": str(hoop), "rho_exact": str(rho) if rho is not None else None})
        selected = {
            "status": "UNRESOLVED_UNIFORM_DOMAIN",
            "G_bounds_exact": pack_bounds((lower, upper)), "hoop_lower_exact": str(hoop),
            "q0_exact": str(q0), "rho_exact": str(rho) if rho is not None else None,
            "depth": depth, "leaves": leaves, "history": history,
            "coefficient_hulls_exact": [{str(k): pack_bounds(v) for k, v in coefficient_hulls(p).items()} for p in polys],
            "denominator_coefficient_hulls_exact": {str(k): pack_bounds(v) for k, v in coefficient_hulls(denominator).items()},
            "rho_goal_exact": str(rho_goal),
            "scope": "Whole closed edge and every nodal value in supplied coefficient box; not point-only positivity.",
        }
        if lower > 0 and hoop > 0 and rho is not None and rho < 1:
            selected["status"] = "CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN"
            if rho <= rho_goal:
                break
    if selected is None:
        raise AssertionError("domain policy produced no attempt")
    return selected


def verify_frozen_domain(fields, domain):
    if domain.get("status") != "CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN":
        raise ValueError("certified whole-box domain required")
    if F(domain["hoop_lower_exact"]) <= 0 or F(domain["G_bounds_exact"][0]) <= 0 or not 0 <= F(domain["rho_exact"]) < 1:
        raise ValueError("positive uniform domain invariants required")
    for p, frozen in zip(fields["q_numerators"], domain["coefficient_hulls_exact"]):
        actual = coefficient_hulls(p)
        if {str(k) for k in actual} != set(frozen):
            raise ValueError("changed interval polynomial support")
        for k, bounds in actual.items():
            lo, hi = map(F, frozen[str(k)])
            if not lo <= bounds[0] <= bounds[1] <= hi:
                raise ValueError("fresh polynomial coefficient leaves frozen uniform domain")
    frozen = domain["denominator_coefficient_hulls_exact"]
    actual = coefficient_hulls(fields["hoop_denominator"])
    if {str(k): pack_bounds(v) for k, v in actual.items()} != frozen:
        raise ValueError("changed exact hoop denominator")


def poly_norm_upper(poly, lo, hi):
    lo, hi = F(lo), F(hi)
    if not 0 <= lo < hi or any(j != 0 or (i < 0 and lo == 0) for i, j in poly):
        raise ValueError("regular physical-radius Laurent polynomial required")
    def magnitude(value):
        return max(map(abs, rational_hull(value)))
    raw = sum((magnitude(v) * (hi**i if i >= 0 else lo**i) for (i, _), v in poly.items()), F(0))
    denominator = max(0, -min((i for i, _ in poly), default=0))
    numerator = T.mul(poly, T.power(T.X, denominator))
    shifted = T.compose(numerator, {(0, 0): (lo + hi) / 2, (1, 0): (hi - lo) / 2}, {})
    centered = sum((magnitude(v) for v in shifted.values()), F(0)) / (lo**denominator if denominator else 1)
    return min(raw, centered)


def uniform_series(e, q0, order, K, tension):
    if F(q0) <= 0:
        raise ValueError("positive expansion center required")
    K, tension = A(K), A(tension)
    root = A(q0).sqrt()
    c1 = S.binomial_coefficients(F(-1, 2), order)
    c3 = S.binomial_coefficients(F(-3, 2), order)
    prefA = (tension - K) / (2 * root)
    prefB = -(tension - K) / (4 * A(q0) * root)
    p = T.ONE
    outA, outB = T.constant(K / 2), {}
    for n in range(order + 1):
        outA = T.add(outA, T.scale(p, prefA * A(c1[n])))
        outB = T.add(outB, T.scale(p, prefB * A(c3[n])))
        if n < order:
            p = T.mul(p, e)
    return outA, outB


def material_boxes(row, face):
    prefix = "basal" if face == 0 else "apical"
    K = tuple(map(F, row["parameter_cell"][prefix + "_area_modulus"]))
    tension = tuple(map(F, row["parameter_cell"][prefix + "_surface_tension"]))
    if not 0 < K[0] <= K[1] or not 0 <= tension[0] <= tension[1]:
        raise ValueError("same physical cell with positive K/nonnegative tension required")
    return iball(K), iball(tension), K, tension


def common_nodes(target, box):
    if box["root_membership_assumed"] is not False:
        raise ValueError("no assumed root membership permitted")
    reference = [[str(T.rational(v)) for v in row] for row in target.a["reference"]]
    if reference != box["reference_nodes_exact"] or target.a["triangles"].tolist() != box["triangles"]:
        raise ValueError("same reference geometry/connectivity required")
    centers, radii = box["current_centers_exact"], box["coordinate_radii_exact"]
    if len(centers) != 325 or len(radii) != 325 or any(len(x) != 2 for x in centers + radii):
        raise ValueError("325 two-coordinate common nodes required")
    nodes = []
    for center, radius in zip(centers, radii):
        row = []
        for c, r in zip(center, radius):
            c, r = F(c), F(r)
            if r < 0:
                raise ValueError("nonnegative exact coordinate radii required")
            row.append(iball((c - r, c + r)))
        nodes.append(row)
    for node in range(325):
        if target.a["reference"][node, 0] == 0 and not nodes[node][0].is_zero():
            raise ValueError("axis radial trace must remain EXACT zero")
    return nodes


def local_site_domain(target, box, face, edge, iq):
    maps = L.actual_maps(target, face, edge, iq)
    centers = [F(box["current_centers_exact"][node][c]) for node in maps["nodes"] for c in (0, 1)]
    radii = [F(box["coordinate_radii_exact"][node][c]) for node in maps["nodes"] for c in (0, 1)]
    Y = [(c - r, c + r) for c, r in zip(centers, radii)]
    Bs, Bc = maps["star_full"], maps["cache_normalized_full"]
    delta = [[s - c for s, c in zip(sr, cr)] for sr, cr in zip(Bs, Bc)]
    def apply(B):
        return [sum_intervals([iscale(y, b) for b, y in zip(row, Y)]) for row in B]
    qc, dq = apply(Bc), apply(delta)
    line = [(v[0] + min(F(0), d[0]), v[1] + max(F(0), d[1])) for v, d in zip(qc, dq)]
    tangent2 = iadd(S.isquare(line[1]), S.isquare(line[2]))
    if line[0][0] <= 0 or tangent2[0] <= 0 or maps["auxiliary_normalized_weight_exact"] <= 0:
        raise ValueError("entire local common-box line must have positive hoop/tangent and weight")
    return {
        "quadrature_index": iq, "nodes": maps["nodes"], "global_indices": maps["global_indices"],
        "selected_local": maps["selected_local"],
        "normalization_D_exact": list(map(str, maps["normalization_D_exact"])),
        "cache_q_exact": [pack_bounds(v) for v in qc], "delta_q_exact": [pack_bounds(v) for v in dq],
        "convex_line_exact": [pack_bounds(v) for v in line],
        "tangent_norm_squared_lower_exact": str(tangent2[0]),
        "auxiliary_normalized_weight_exact": str(maps["auxiliary_normalized_weight_exact"]),
        "status": "CERTIFIED_LOCAL_FULLY_COEFFICIENTWISE_LINE_DOMAIN_NO_CACHE_SUM",
    }


def sum_intervals(values):
    value = F(0), F(0)
    for x in values:
        value = iadd(value, x)
    return value


def cache_edge_correction(target, nodes, face, edge, K, tension, frozen_sites):
    fields = target.surface(face, edge, nodes=nodes)
    n = len(fields["free_local"])
    g = [arb(0) for _ in range(n)]
    H = arb_mat(n, n)
    summaries = []
    if len(frozen_sites) != 12:
        raise ValueError("all twelve frozen auxiliary sites required")
    for iq in range(12):
        maps = L.actual_maps(target, face, edge, iq)
        frozen = frozen_sites[iq]
        if iq != frozen["quadrature_index"] or maps["global_indices"] != fields["global_indices"] or maps["global_indices"] != frozen["global_indices"] or list(map(str, maps["normalization_D_exact"])) != frozen["normalization_D_exact"]:
            raise ValueError("changed local selector/normalization")
        Y = [nodes[node][c] for node in maps["nodes"] for c in (0, 1)]
        local = L.consistency_bound(maps["star_full"], maps["cache_normalized_full"], Y, maps["selected_local"], K, tension)
        weight = maps["auxiliary_normalized_weight_exact"]
        if str(weight) != frozen["auxiliary_normalized_weight_exact"] or weight <= 0:
            raise ValueError("same strictly positive old-normalized weight required")
        for i in range(n):
            g[i] += A(weight) * local["force_difference"][i, 0]
            for j in range(n):
                H[i, j] += A(weight) * local["H_difference"][i, j]
        summaries.append({"quadrature_index": iq, "weight_applied_exactly_once": str(weight), "unweighted_F_max": core.pack(L.maxabs(local["force_difference"].entries())), "unweighted_H_max": core.pack(L.maxabs(local["H_difference"].entries()))})
    return g, H, summaries


def run_edge(target, nodes, row, case, rule, order):
    face, edge = case["face"], case["edge"]
    fields = target.surface(face, edge, nodes=nodes)
    if fields["global_indices"] != case["global_indices"]:
        raise ValueError("same edge free-column selector required")
    verify_frozen_domain(fields, case["domain"])
    jet = JET.surface_squared_jacobian_jet(fields["q_physical_R"], fields["B_physical_R"])
    K, tension, _, _ = material_boxes(row, face)
    q0, rho = F(case["domain"]["q0_exact"]), F(case["domain"]["rho_exact"])
    e = T.scale(T.sub(jet["value"], T.constant(q0)), 1 / q0)
    polyA, polyB = uniform_series(e, q0, order, K, tension)
    functional = S.SurfaceFunctional(target, fields, face, edge, rule)
    truncation = JET.contract_volumetric_defect(jet, polyA, polyB, functional.defect)
    mass = core.upper(functional.mass_integral + functional.mass_rule)
    dmax = max(map(abs, rational_hull(tension - K)))
    root = A(q0).sqrt()
    tailA = A(dmax * S.binomial_tail(F(-1, 2), order, rho)) / (2 * root)
    tailB = A(dmax * S.binomial_tail(F(-3, 2), order, rho)) / (4 * A(q0) * root)
    norms = [poly_norm_upper(p, fields["Rlo"], fields["Rhi"]) for p in jet["gradient"]]
    Ft = [mass * tailA * A(v) for v in norms]
    Ht = {(i, j): mass * (tailA * A(poly_norm_upper(p, fields["Rlo"], fields["Rhi"])) + tailB * A(norms[i] * norms[j])) for (i, j), p in jet["hessian_upper"].items()}
    cacheg, cacheH, sites = cache_edge_correction(target, nodes, face, edge, K, tension, case["cache_site_domains"])
    def inflate(value, radius):
        return A(value) + arb(0, core.upper(radius))
    Fout = [inflate(value, Ft[i]) + cacheg[i] for i, value in enumerate(truncation["gradient"])]
    Hout = {key: inflate(value, Ht[key]) + cacheH[key[0], key[1]] for key, value in truncation["hessian_upper"].items()}
    return {
        "status": "UNIFORM_SURFACE_EDGE_COMPONENT_UNREVIEWED", "face": face, "edge": edge,
        "global_indices": fields["global_indices"], "q0_exact": str(q0), "rho_exact": str(rho),
        "precision_bits": ctx.prec, "order": order,
        "force_integral_minus_cache": [core.pack(v) for v in Fout],
        "force_auxiliary_minus_cache": [core.pack(v) for v in cacheg],
        "force_two_mass_tail_upper": [core.pack(v) for v in Ft],
        "H_integral_minus_cache_upper": [{"i": i, "j": j, "enclosure": core.pack(v), "auxiliary_minus_cache": core.pack(cacheH[i, j]), "two_mass_tail_upper": core.pack(Ht[i, j])} for (i, j), v in Hout.items()],
        "actual_cache_sites": sites,
        "material_mass": core.pack(functional.mass_integral), "auxiliary_mass": core.pack(functional.mass_rule),
        "unique_moment_defects": len(functional.requested_defect_keys),
        "uniform_quantifier": "Every supplied common-box nodal field and SAME physical K/T parameter cell; no root membership assumed.",
        "scope": SCOPE,
    }


def synthetic_tests():
    checks = {}
    def ck(name, value):
        if value is not True:
            raise AssertionError(name)
        checks[name] = True
    for n in range(7):
        p = {(n, 0): F(1)}
        for lo, hi in ((F(0), F(1)), (F(1, 3), F(3, 4))):
            hull, bernstein = bernstein_bounds(p, lo, hi)
            ck(f"monomial_Bernstein_endpoints_{n}_{lo}", hull[0] == lo**n and hull[1] == hi**n and len(bernstein) == n + 1)
    # Exact quadratic min is zero, despite positive endpoint samples.
    p = {(0, 0): F(1), (1, 0): F(-4), (2, 0): F(4)}
    hull, _ = bernstein_bounds(p, F(0), F(1))
    ck("Bernstein_does_not_assume_positive_endpoints_sufficient", hull[0] <= 0 <= hull[1])
    toy = {"q_numerators": [T.ONE, p, {}], "hoop_denominator": T.ONE}
    ck("interior_collapse_remains_unresolved", uniform_G_bounds(toy, max_depth=3)["status"] == "UNRESOLVED_UNIFORM_DOMAIN")
    uncertain = {(0, 0): arb(1, A(2))}
    toy = {"q_numerators": [uncertain, T.ONE, {}], "hoop_denominator": T.ONE}
    ck("positive_midpoint_does_not_certify_wide_hoop_box", uniform_G_bounds(toy, max_depth=2)["status"] == "UNRESOLVED_UNIFORM_DOMAIN")
    for radius in (F(1, 100), F(1, 10**80)):
        p = {(0, 0): arb(0, A(radius)), (1, 0): F(1)}
        hull, b = bernstein_bounds(p, F(0), F(1))
        ck("zero_crossing_coefficient_retained_Bernstein_" + str(radius), hull[0] < 0 and hull[1] > 1)
        ck("safe_polynomial_zero_crossing_retained_" + str(radius), (0, 0) in T.clean(p))
    for lo in (F(0), F(1)):
        ref = [[lo, F(0)], [lo + F(1, 2), F(0)], [lo + F(1), F(0)]]
        field = T.surface_fields(ref, ref)
        domain = uniform_G_bounds(field, max_depth=3)
        ck("identity_uniform_domain_" + str(lo), domain["status"] == "CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN" and F(domain["G_bounds_exact"][0]) <= 1 <= F(domain["G_bounds_exact"][1]))
    for rho in (F(0), F(1, 100), F(1, 8)):
        for alpha in (F(-1, 2), F(-3, 2)):
            tail = S.binomial_tail(alpha, 8, rho)
            for x in (-rho, rho):
                exact = 1 / (1 + A(x)).sqrt()
                if alpha == F(-3, 2):
                    exact /= 1 + A(x)
                trunc = sum((a * x**n for n, a in enumerate(S.binomial_coefficients(alpha, 8))), F(0))
                ck(f"reviewed_binomial_endpoint_tail_{rho}_{alpha}_{x}", core.upper(exact - A(trunc)) <= A(tail))
    interval_poly = {(0, 0): iball((F(-1, 100), F(1, 100))), (1, 0): iball((F(9, 10), F(11, 10)))}
    bound = poly_norm_upper(interval_poly, F(0), F(1))
    ck("uniform_interval_polynomial_norm_covers_endpoint", bound >= F(111, 100))
    polyA, polyB = uniform_series({(1, 0): F(1, 100)}, F(1), 4, iball((F(1), F(2))), iball((F(2), F(3))))
    ck("uniform_series_keeps_KT_zero_crossing_terms", bool(polyA) and bool(polyB) and any(not v.is_zero() and v.contains(0) for v in polyB.values()))
    # The nilpotent helper's direct cleanup must not become the former !=0 bug.
    ck("local_cache_jet_zero_crossing_retained", L.ZERO in L.jclean({L.ZERO: arb(0, 1)}))
    # Domain guard explicitly rejects caller-supplied point-only proofs.
    try:
        verify_frozen_domain({"q_numerators": [T.ONE] * 3, "hoop_denominator": T.ONE}, {"status": "CERTIFIED_POSITIVE_POINT_G_ON_FULL_EDGE"})
        ck("point_only_domain_rejected", False)
    except ValueError:
        ck("point_only_domain_rejected", True)
    valid_review = {"status": "PASS", "checks": {"one": True}, "source_sha256": "s", "protocol_sha256": "p", "preflight_sha256": "f"}
    validate_review(valid_review, "s", "p", "f")
    ck("valid_prelaunch_guard_accepted", True)
    for label, changed in (("empty", {}), ("false", {"one": False}), ("integer_truthy", {"one": 1}), ("list", [True])):
        try:
            validate_review({**valid_review, "checks": changed}, "s", "p", "f")
            ck("invalid_prelaunch_checks_rejected_" + label, False)
        except ValueError:
            ck("invalid_prelaunch_checks_rejected_" + label, True)
    try:
        validate_review(valid_review, "changed", "p", "f")
        ck("wrong_source_prelaunch_rejected", False)
    except ValueError:
        ck("wrong_source_prelaunch_rejected", True)
    for label, function in (("nonfinite_Bernstein", lambda: bernstein_bounds({(0, 0): arb("nan")}, F(0), F(1))), ("axis_Laurent_norm", lambda: poly_norm_upper({(-1, 0): F(1)}, F(0), F(1)))):
        try:
            function()
            ck("invalid_input_rejected_" + label, False)
        except ValueError:
            ck("invalid_input_rejected_" + label, True)
    return checks


def bind(bindings, path, digest=None):
    path = Path(path).resolve()
    key = path.relative_to(ROOT.resolve()).as_posix()
    current = sha(path)
    if digest is not None and current != digest:
        raise ValueError("changed frozen dependency: " + key)
    if key in bindings and bindings[key] != current:
        raise ValueError("conflicting dependency hash")
    bindings[key] = current


def validate_review(review, source_hash, protocol_hash, preflight_hash):
    S.validate_prelaunch(review, source_hash, protocol_hash, preflight_hash)


def prepare():
    set_precision(PRECISION)
    if (OUT / "protocol.json").exists():
        raise ValueError("preserve frozen preparation")
    checks = synthetic_tests()
    preparation, geometry, box = read(TARGET_PREPARATION), read(GEOMETRY_PROTOCOL), read(BOXFILE)
    row = next(x for x in preparation["records"] if x["record_id"] == "optoGEF2_small_col_B")
    if row["parameter_cell_sha256"] != box["parameter_cell_sha256"] or row["parameter_cell_sha256"] != geometry["parameter_cell_sha256"] or sha(BOXFILE) != geometry["prepared_box_sha256"]:
        raise ValueError("SAME unchanged C and common geometry required")
    if sha(ROOT / row["point_cache_file"]) != row["point_cache_sha256"]:
        raise ValueError("same old Q cache required")
    target = T.Target(read(ROOT / row["point_cache_file"]))
    nodes = common_nodes(target, box)
    bindings = {}
    for group in (preparation["source_hashes"], preparation["input_hashes"]):
        for name, digest in group.items():
            bind(bindings, ROOT / name, digest)
    for file in (SOURCE, Path(S.__file__), Path(T.__file__), Path(JET.__file__), Path(M.__file__), Path(core.__file__), Path(L.__file__), Path(L.target.__file__), TARGET_PREPARATION, GEOMETRY_PROTOCOL, BOXFILE, ROOT / "evidence/exact_surface_cache_consistency_independent_review.json", ROOT / "evidence/exact_safe_polynomial_kernel_independent_review.json", ROOT / "evidence/exact_fe_moments_independent_review.json"):
        bind(bindings, file)
    local_review = read(ROOT / "evidence/exact_surface_cache_consistency_independent_review.json")
    if local_review["status"] != "PASS_INDEPENDENT_LOCAL_SURFACE_CACHE_CONSISTENCY_REVIEW" or not local_review["checks"] or not all(v is True for v in local_review["checks"].values()) or local_review["source_sha256"] != sha(Path(L.__file__)):
        raise ValueError("exact-source local cache review required")
    safe_review = read(ROOT / "evidence/exact_safe_polynomial_kernel_independent_review.json")
    if safe_review["status"] != "PASS" or not safe_review["checks"] or not all(v is True for v in safe_review["checks"].values()) or safe_review["source_sha256"] != sha(Path(T.__file__)) or safe_review["helper_source_sha256"] != sha(Path(JET.__file__)):
        raise ValueError("safe arithmetic independent review required")
    cases = []
    for face in (0, 1):
        material_boxes(row, face)
        for edge in range(12):
            ids = list(map(int, target.a[f"surface_{face}_edges"][edge]))
            if len({T.rational(target.a["reference"][node, 1]) for node in ids}) != 1:
                raise ValueError("horizontal exact material reference required")
            fields = target.surface(face, edge, nodes=nodes)
            domain = uniform_G_bounds(fields)
            sites = [local_site_domain(target, box, face, edge, iq) for iq in range(12)]
            cases.append({"face": face, "edge": edge, "global_indices": fields["global_indices"], "domain": domain, "cache_site_domains": sites})
            print(json.dumps({"prepared_domain": [face, edge], "status": domain["status"], "rho_display_only": float(F(domain["rho_exact"])) if domain["rho_exact"] is not None else None, "depth": domain["depth"]}), flush=True)
    checks["all24_uniform_domains_certified"] = all(x["domain"]["status"] == "CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN" for x in cases)
    checks["all288_exact_fullY_local_line_domains_certified"] = all(len(x["cache_site_domains"]) == 12 for x in cases)
    checks["same_full24_predeclared_selection"] = [(x["face"], x["edge"]) for x in cases] == [(f, e) for f in (0, 1) for e in range(12)]
    protocol = {
        "version": 2, "utc": stamp(), "source_and_input_hashes": bindings,
        "selected_record": row, "point_cache_file": row["point_cache_file"],
        "common_geometry_box_file": BOXFILE.relative_to(ROOT).as_posix(), "common_geometry_box_sha256": sha(BOXFILE),
        "parameter_cell_sha256": row["parameter_cell_sha256"], "common_geometry_domain_id": box["domain_id"],
        "precision_bits": PRECISION, "truncation_order": ORDER,
        "domain_policy": {"max_depth": MAX_DEPTH, "rho_goal_exact": str(RHO_GOAL), "method": "Exact rational interval coefficient hulls, power-to-Bernstein convex hull on closed chart leaves, regularized hoop denominator, positivity for whole nodal box."},
        "cases": cases,
        "coefficient_policy": "SAFE exact-zero-only sparse polynomial cleanup. All nonexact zero-crossing interval coefficients retained. Shared nodal/material dependence may be overestimated, never replaced by midpoint values.",
        "cache_policy": "At each of 12 auxiliary sites: fullY deltaq before intervalization, reviewed third-derivative FTC bound and all3 signed Hessian terms; multiply actual old-normalized positive weight once, then assemble selected free indices.",
        "tail_policy": "Both positive material and auxiliary masses multiply uniform scalar tails and interval-coefficient derivative polynomial sup norms; no Gaussian exactness assumed.",
        "commonbox_ancestry": "Explicit supplied set only. Old unsafe point conclusions in historical geometry provenance are not premises for root membership or this uniform component.",
        "scope": SCOPE, "root_membership_assumed": False,
    }
    write(OUT / "protocol.json", protocol)
    write(OUT / "synthetic_tests.json", {"status": "PASS" if all(checks.values()) else "UNRESOLVED_UNIFORM_DOMAINS", "checks": checks, "source_sha256": sha(SOURCE), "protocol_sha256": sha(OUT / "protocol.json")})
    preflight = {"status": "PASS_UNIFORM_SURFACE_PREPARATION_ONLY" if all(checks.values()) else "UNRESOLVED_UNIFORM_SURFACE_PREPARATION_ONLY", "source_sha256": sha(SOURCE), "protocol_sha256": sha(OUT / "protocol.json"), "synthetic_tests_sha256": sha(OUT / "synthetic_tests.json"), "checks": checks, "actual_edge_integrals_run": False, "actual_cache_sums_run": False, "uniform_complete_integral_root_proved": False, "scope": SCOPE}
    write(OUT / "preflight.json", preflight)
    print(json.dumps({"status": preflight["status"], "checks": len(checks), "source_sha256": sha(SOURCE), "protocol_sha256": sha(OUT / "protocol.json"), "preflight_sha256": sha(OUT / "preflight.json"), "bindings": len(bindings)}), flush=True)


def run():
    p = read(OUT / "protocol.json")
    preflight = read(OUT / "preflight.json")
    if preflight["status"] != "PASS_UNIFORM_SURFACE_PREPARATION_ONLY":
        raise ValueError("all uniform domains must pass before actual integration")
    validate_review(read(OUT / "prelaunch_independent_review.json"), sha(SOURCE), sha(OUT / "protocol.json"), sha(OUT / "preflight.json"))
    for name, digest in p["source_and_input_hashes"].items():
        if sha(ROOT / name) != digest:
            raise ValueError("changed frozen source/input: " + name)
    set_precision(p["precision_bits"])
    target = T.Target(read(ROOT / p["point_cache_file"]))
    box = read(ROOT / p["common_geometry_box_file"])
    nodes = common_nodes(target, box)
    rule = target.auxiliary_rule()
    if p["parameter_cell_sha256"] != box["parameter_cell_sha256"]:
        raise ValueError("same shared parameter cell required")
    write(OUT / "started.json", {"status": "RUNNING", "utc": stamp(), "source_sha256": sha(SOURCE), "protocol_sha256": sha(OUT / "protocol.json")})
    start = time.monotonic()
    g, H, files = [arb(0) for _ in range(611)], {}, []
    for case in p["cases"]:
        result = run_edge(target, nodes, p["selected_record"], case, rule, p["truncation_order"])
        indices = result["global_indices"]
        for i, v in zip(indices, result["force_integral_minus_cache"]):
            g[i] += iball(tuple(map(F, (v["lower_exact"], v["upper_exact"]))))
        for item in result["H_integral_minus_cache_upper"]:
            key = tuple(sorted((indices[item["i"]], indices[item["j"]])))
            v = item["enclosure"]
            H[key] = H.get(key, arb(0)) + iball((F(v["lower_exact"]), F(v["upper_exact"])))
        file = OUT / "edges" / f"face{case['face']}_edge{case['edge']:02d}.json"
        write(file, result)
        files.append({"face": case["face"], "edge": case["edge"], "file": file.relative_to(ROOT).as_posix(), "sha256": sha(file)})
        print(json.dumps({"uniform_edge_complete": [case["face"], case["edge"]]}), flush=True)
    for name, digest in p["source_and_input_hashes"].items():
        if sha(ROOT / name) != digest:
            raise ValueError("source/input changed during uniform run: " + name)
    result = {"status": "UNIFORM_ALL24_SURFACE_COMPONENT_ASSEMBLED_UNREVIEWED", "source_sha256": sha(SOURCE), "protocol_sha256": sha(OUT / "protocol.json"), "edges": files, "force_integral_minus_cache": [core.pack(v) for v in g], "H_upper": [{"i": i, "j": j, "enclosure": core.pack(v)} for (i, j), v in sorted(H.items())], "parameter_cell_sha256": p["parameter_cell_sha256"], "common_geometry_domain_id": p["common_geometry_domain_id"], "all_frozen_bindings_current": True, "seconds": time.monotonic() - start, "root_transfer_proved": False, "scope": SCOPE}
    write(OUT / "results.json", result)
    print(json.dumps({"status": result["status"], "seconds": result["seconds"]}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("prepare", "run", "test"))
    args = parser.parse_args()
    if args.command == "prepare":
        prepare()
    elif args.command == "run":
        run()
    else:
        set_precision(PRECISION)
        checks = synthetic_tests()
        print(json.dumps({"status": "PASS", "checks": len(checks)}))
