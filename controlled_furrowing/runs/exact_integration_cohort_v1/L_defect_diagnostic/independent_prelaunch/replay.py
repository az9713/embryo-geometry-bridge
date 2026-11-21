"""Independent review of the fixed L element-21 force diagnostic, without integration."""
import ast
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch

from flint import arb

OUT = Path(__file__).resolve().parent
JOB = OUT.parent
ROOT = JOB.parents[2]
FROZEN = {
    "element21_force.py": "9a5f93bb3d4a39020e28bee165ba606d1d17f63e44a17ee758959aa05109d329",
    "protocol.json": "4b2f6063653ed89c1e2c0c5d4ba6e74077fae92c62ef09c38f35101f4559ef8b",
    "preflight.json": "2b9ce8dfa112198c549e513e19ae1c34084e1a9d8bce9e784d689a50db2f4909",
}


def err(action):
    try:
        action()
    except ValueError:
        return
    raise AssertionError("invalid fixture accepted")


def cache_site(X):
    unit = [[Q(i == j) for j in range(5)] for i in range(5)]
    star = [row.copy() for row in unit]
    star[1][1] = Q(7, 8)
    columns = [1, 0, 3]
    y = list(map(X.A, [1, Q(3, 16), 0, 1, 1]))
    drive = Q(5, 4)
    actual, J = X._force_cache_site(star, unit, y, columns, X.A(drive))
    # Both point strains have J=1 and differ only in b. The gradient component
    # is (1+drive)*b and the star pullback supplies the second factor 7/8.
    expected = (Q(7, 8)**2 - 1) * (1 + drive) * Q(3, 16)
    assert actual[0].contains(X.A(expected))
    assert all(v.contains(0) for v in actual[1:]) and J > 0
    err(lambda: X._force_cache_site(star, unit, [arb(0), *y[1:]], columns, X.A(drive)))


def gain_oracle(X):
    # Analytic P2 derivative values from barycentric lambda=(5/8,1/4,1/8)
    # at physical (R,Z)=(3/2,1/4). Input to polynomial evaluation is zeta=Z-2/3.
    corners = ((Q(1), Q(0)), (Q(3), Q(0)), (Q(1), Q(2)))
    sample = (Q(3, 2), Q(1, 4)-Q(2, 3))
    lam = (Q(5, 8), Q(1, 4), Q(1, 8))
    d = (Q(-1, 2), Q(0), Q(1, 2))
    grad = [(4*l-1)*v for l, v in zip(lam, d)] + [
        4 * (d[i]*lam[j] + lam[i]*d[j]) for i, j in [(0, 1), (1, 2), (2, 0)]]
    exact_diff_eval = lambda p: X.A(sum(v * sample[0]**i * sample[1]**j for (i,j),v in p.items()))
    fun = SimpleNamespace(apply=exact_diff_eval)
    star, frozen = [[[Q(0)]*12 for _ in range(5)] for _ in range(2)]
    for slot in (1, 3):
        axis = (slot - 1)//2
        for j in range(6):
            star[slot][2*j+axis] = grad[j]
            frozen[slot][2*j+axis] = grad[j] * Q(3, 4 if axis else 5)
    weight = Q(3, 16)
    rule = [{"star_full": star, "cache_full": frozen, "auxiliary_normalized_weight_exact": weight}]
    ys = [Q(i - 4, 16) for i in range(12)]
    delta = Q(1, 64)
    balls = [X.A(y-delta).union(X.A(y+delta)) for y in ys]
    cols = [0, 3, 4, 11]
    fields = {"corners": corners, "free_local": cols, "active": Q(1)}
    dose = Q(9, 8)
    result = X._gain_derivative(SimpleNamespace(dose=dose), fields, fun, rule, balls)
    for ys_end in ([v-delta for v in ys], [v+delta for v in ys]):
        for k, actual in zip(cols, result):
            axis = k % 2
            s = 1 if axis == 0 else 3
            material = grad[k//2] * sum(grad[j] * ys_end[2*j+axis] for j in range(6))
            aux_s = star[s][k] * sum(v*y for v, y in zip(star[s], ys_end))
            aux_c = frozen[s][k] * sum(v*y for v, y in zip(frozen[s], ys_end))
            expected = dose * (material + weight * (aux_s - aux_c))
            assert actual.contains(X.A(expected))
    inactive = X._gain_derivative(SimpleNamespace(dose=dose), {**fields, "active": 0}, fun, rule, balls)
    assert all(v.is_zero() for v in inactive)


def centered_oracle(X):
    # Toy symmetric Hessian and exact affine g dependence. Compare an
    # independent mixed quadratic force law at every state/gain corner.
    f = lambda x, y, g: ((2+g)*x + 3*y + x*x, 3*x + (5+2*g)*y)
    x0, y0, g0 = Q(1, 2), Q(-1, 4), Q(7, 8)
    xran, yran, gran = (Q(1, 4), Q(3, 4)), (Q(-3, 4), Q(1, 4)), (Q(1, 2), Q(1))
    hess = {(0,0):(2+gran[0]+2*xran[0],2+gran[1]+2*xran[1]), (0,1):(Q(3),Q(3)),
            (1,1):(5+2*gran[0],5+2*gran[1])}
    ball = lambda a, b: X.A(a).union(X.A(b))
    increment = X._state_increment([ball(xran[0]-x0,xran[1]-x0),ball(yran[0]-y0,yran[1]-y0)], hess)
    dg = ball(gran[0]-g0, gran[1]-g0)
    a, b = f(x0,y0,g0)
    centered = (X.A(a) + increment[0] + dg*ball(*xran), X.A(b) + increment[1] + 2*dg*ball(*yran))
    for x, y, g in product(xran, yran, gran):
        assert all(ball.contains(X.A(exact)) for ball, exact in zip(centered, f(x,y,g)))
    assert X._intersect(X.A(-1).union(X.A(0)), X.A(0).union(X.A(1)))[0] <= 0
    err(lambda: X._intersect(X.A(-1), X.A(1)))


def run():
    if not __debug__:
        raise RuntimeError("Replay requires enabled assertions")
    bindings = {}

    def raw(path, expected=None):
        filename = ROOT / path
        assert filename.resolve().is_relative_to(ROOT)
        data = filename.read_bytes()
        h = sha256(data).hexdigest()
        assert expected is None or h == expected, path
        assert bindings.setdefault(path, h) == h
        return data

    rel = JOB.relative_to(ROOT).as_posix()
    source = raw(f"{rel}/element21_force.py", FROZEN["element21_force.py"])
    assert source == raw(f"{rel}/element21_force_frozen.py", FROZEN["element21_force.py"])
    protocol = json.loads(raw(f"{rel}/protocol.json", FROZEN["protocol.json"]))
    pre = json.loads(raw(f"{rel}/preflight.json", FROZEN["preflight.json"]))
    assert pre["status"] == "PASS_PREPARATION_ONLY" and all(v is True for v in pre["checks"].values())
    assert pre["source_sha256"] == protocol["source_sha256"] == FROZEN["element21_force.py"]
    assert pre["protocol_sha256"] == FROZEN["protocol.json"]
    assert protocol["actual_record_force_integration_run"] is pre["actual_record_force_integration_run"] is False
    for name, expected in protocol["input_hashes"].items():
        raw(name, expected)
    sys.path.insert(0, str(JOB))
    import element21_force as X
    forbidden = AssertionError("review attempted a real integral, force profile or all-site correction")
    with patch.object(X.E.P.Functional, "__init__", side_effect=forbidden), \
            patch.object(X.U.C, "exact_bulk_maps", side_effect=forbidden), \
            patch.object(X, "_force_parts", side_effect=forbidden), \
            patch.object(X.U, "uniform_cache", side_effect=forbidden):
        assets, context, case, old, _, ids, nodes, zero, setup, differences, H = X._geometry()
    assert assets == protocol["input_hashes"] and protocol["record_id"] == context.record_id == X.RECORD
    assert case["element"] == protocol["element"] == X.ELEMENT == 21
    assert case["degree"] == old["degree"] == protocol["order"] == 12
    assert case["rho_exact"] == old["rho_exact"] == protocol["rho_exact"]
    assert case["chart"] == old["chart"] == protocol["chart"] == "off_axis"
    assert protocol["same_C_sha256"] == context.row["parameter_cell_sha256"]
    assert Q(protocol["same_dose_exact"]) == context.dose
    point, point_ids = X.E._bulk_fields(context.target, zero, X.ELEMENT)
    _, box, wide_drive, _, mass, _, _ = setup
    assert int(context.target.a["active"][X.ELEMENT,0]) == point["active"] == box["active"] == 1
    assert case["active"] == old["active"] == "1"
    assert context.drive == context.dose * context.base["common_gain"]
    for endpoint in context.cell["common_gain"]:
        assert wide_drive.contains(X.A(context.dose * endpoint))
    assert ids == point_ids == protocol["element_nodes"] and len(ids) == 6
    assert set(box["free_local"]) == set(range(12)) and len(differences) == 12 and len(H) == 78
    for i, (node, axis) in enumerate(product(ids, range(2))):
        l, r = X.E.SUR.rational_hull(nodes[node][axis])
        assert l <= zero[node][axis] <= r
        a, b = X.E.SUR.rational_hull(differences[i])
        assert a <= 0 <= b
    assert mass == Q(case["geometry_two_mass_upper_exact"]) > 0
    actual_tail, _ = X.M.volumetric_derivative_tail_bounds(12, Q(case["rho_exact"]), 1, 5000)
    assert actual_tail == Q(case["phi_prime_tail_upper_exact"]) > 0
    assert int(protocol["precision_bits"]) == X.precision.prec == 640
    diagnostic = json.loads(raw(X.REDUCTION + X.RECORD + "_signed_endpoint_diagnostics.json"))
    previous = json.loads(raw(X.REDUCTION + X.RECORD + "_reduction_review.json"))
    assert previous["status"] == "PASS_REDUCTION_UNRESOLVED_SUFFICIENT_ROOT_TEST"
    assert diagnostic["coordinates_strictly_excluding_zero"] == 0
    assert Q(diagnostic["certified_RF_infinity_norm_lower_exact"]) == 0
    eta, z, b, r = [Q(previous[k]) for k in ("cached_eta_upper_exact", "cached_Z_upper_exact",
                                           "b_upper_exact", "normal_radius_exact")]
    force_zero_lhs = eta + (z+b)*r
    assert force_zero_lhs >= r

    definitions = {f.name: f for f in ast.parse(source).body if isinstance(f, ast.FunctionDef)}
    force_calls = {ast.unparse(n.func) for n in ast.walk(definitions["_force_parts"]) if isinstance(n, ast.Call)}
    assert {"_series_first", "JET.moment_dual", "JET.contract_factor", "fun.apply", "_force_cache_site", "E._norm"} <= force_calls
    assert "JET.contract_volumetric_defect" not in force_calls
    synthetics = X.tests()
    assert synthetics and all(v is True for v in synthetics.values())
    uncertainty = arb(0, X.A(Q(1,64)))
    a = X._series_first({(1,0):uncertainty}, 5)
    first, _ = X.U.series({(1,0):uncertainty}, 5)
    assert {(i,0) for i in range(6)} == a.keys() == first.keys()
    assert all(X._bounds(v) == X._bounds(first[k]) for k, v in a.items())
    cache_site(X)
    gain_oracle(X)
    centered_oracle(X)
    replay = Path(__file__).relative_to(ROOT).as_posix()
    raw(replay)
    for name, expected in tuple(bindings.items()):
        raw(name, expected)
    report = {
        "status": "PASS",
        "scope": "One already-bound active off-axis L element, force-only diagnostic, on the original exact box/C at 640 bits. No root test, fresh Hessian moment or global assembly.",
        "record_id": X.RECORD, "element": X.ELEMENT, "source_sha256": FROZEN["element21_force.py"],
        "protocol_sha256": FROZEN["protocol.json"], "preflight_sha256": FROZEN["preflight.json"],
        "checks": {
            "frozen_source_protocol_preflight_and_all_inputs_match": True,
            "exact_same_box_C_radius_projection_active_mask_order_and_rho": True,
            "point_and_direct_signed_force_terms_and_two_positive_mass_tail_correct": True,
            "zero_crossing_series_coefficients_retained": True,
            "single_toy_cache_site_direct_pullback_and_nonpositive_J_guard": True,
            "independent_synthetic_axial_stiffness_oracle_with_cache_correction_and_column_parity": True,
            "centered_MVT_matches_independent_symmetric_toy_force_law_at_all_corners": True,
            "force_zero_still_fails_original_cube_at_retained_b": True,
            "actual_force_profile_and_site_calls_blocked_during_review": True,
        },
        "input_hashes": dict(sorted(bindings.items())),
        "review_source_file": replay, "review_source_sha256": bindings[replay],
        "actual_record_force_integration_run": False,
        "force_zero_selfmap_lhs_over_radius_exact": str(force_zero_lhs/r),
        "specialization": "The point call supplies full exact base drive. Element 21 has active mask exactly one in the authoritative cache, selector, old packet and both field contexts. This review covers that frozen element only.",
        "centered_path": "Hold g at g0 on the straight state segment from Y0 to Y inside the convex original box, using the retained same-box H21. Then hold Y fixed and use the exact affine gain derivative to reach g. Fixed coordinates do not move.",
    }
    if "--freeze" in sys.argv:
        receipt = OUT / "replay.json"
        gate = JOB / "prelaunch_independent_review.json"
        assert not receipt.exists() and not gate.exists()
        with receipt.open("x", encoding="utf8", newline="\n") as stream:
            stream.write(json.dumps(report, indent=2) + "\n")
        report["review_receipt_file"] = receipt.relative_to(ROOT).as_posix()
        report["review_receipt_sha256"] = sha256(receipt.read_bytes()).hexdigest()
        with gate.open("x", encoding="utf8", newline="\n") as stream:
            stream.write(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "checks": len(report["checks"]),
                      "force_zero_selfmap_lhs_over_radius": float(force_zero_lhs/r)}), flush=True)


if __name__ == "__main__":
    run()
