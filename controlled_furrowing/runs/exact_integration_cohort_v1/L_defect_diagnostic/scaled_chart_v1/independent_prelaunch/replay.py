"""Independent scaled-chart prelaunch replay without real integral/cache work."""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from math import comb
from pathlib import Path
import sys
from unittest.mock import patch

from flint import arb

OUT = Path(__file__).resolve().parent
JOB = OUT.parent
ROOT = JOB.parents[3]
SOURCE = "scaled_element21.py"
EXPECTED = {
    "source_sha256": "a173b3759fb66fe66e5f3e442c1bcce9a5a185e67515649a94499f82d1d67f60",
    "protocol_sha256": "76e5086c374fe86f3659ef1d1b099c4a733211245fbb1fa70251b01462c3fc5c",
    "preflight_sha256": "a35814e39106e14514f2afa6134996781d0e08f5036bcd6e22b1aafd5587a505",
}


def exact_cover(S, ball, value):
    left, right = S.E.SUR.rational_hull(ball)
    assert left <= value <= right, (left, value, right)


def closed_form_moments(S):
    # For material R in [1,2], Z in [0,R-1], x=2R-3 and y=3Z/2-1/2.
    # Integrate y^b directly, expand (R-1) afterward, and integrate the
    # remaining signed radial powers with a separate exact primitive.
    geom = S.chart(((1,0),(2,0),(2,1)))
    for a,b,d in ((0,0,0),(0,0,2),(4,3,14),(7,1,13),(3,2,0),(2,0,4)):
        constant, log_coefficient = Q(0), Q(0)
        for axial in range(b+1):
            z_integral = Q(comb(b,axial), axial+1) * Q(3,2)**axial * Q(-1,2)**(b-axial)
            for axial_power in range(axial+2):
                z_term = z_integral * comb(axial+1,axial_power) * (-1)**(axial+1-axial_power)
                for x_power in range(a+1):
                    x_term = comb(a,x_power) * 2**x_power * (-3)**(a-x_power)
                    term = z_term * x_term * Q(3,2)**d
                    p = 1 - d + axial_power + x_power
                    if p == -1:
                        log_coefficient += term
                    else:
                        constant += term * (Q(2)**(p+1)-1)/(p+1)
        expected = S.M.LogMoment(constant, ((Q(2),log_coefficient),))
        assert S.material_moment(geom,a,b,d) == expected
    try:
        S.chart(((0,0),(2,0),(2,1)))
    except ValueError:
        pass
    else:
        raise AssertionError("axis triangle admitted by a positive-radius Laurent chart")
    return geom


def synthetic_force_hessian(S, geom):
    ref = list(geom.corners) + [tuple((a+b)/2 for a,b in zip(geom.corners[i],geom.corners[j]))
                                for i,j in ((0,1),(1,2),(2,0))]
    current = [(R + R*Z/19 - Z/13, Z - R*R/31 + 1/Q(7)) for R,Z in ref]
    free = [11,4,1,2,6]
    chart_q, chart_B = S.numerators(geom,current,free)
    original = S.T.element_fields(geom.corners,current,free)
    points = [(Q(1,2), Q(1,4), Q(1,4)), (Q(1,7),Q(1,7),Q(5,7))]
    coordinates = [S.chart_point(geom,p) for p in points]
    weights = [Q(7,5), Q(-5,3)]
    drive = Q(5,7)
    order = 5
    classes = set()

    class PointDifference:
        def defect(self, i,j,d):
            classes.add(d)
            return S.A(sum(w*x**i*y**j/r**d for w,(x,y,r) in zip(weights,coordinates)))

        def apply(self, poly,d):
            return sum((S.A(c) * self.defect(i,j,d) for (i,j),c in poly.items()),arb(0))

    F,H,g_norm,H_norm = S._eval_parts(geom,current,free,PointDifference(),S.A(drive),order)
    assert classes == {0,2,order+1,order+2}
    n = len(free)
    pairs = list(combinations_with_replacement(range(n),2))
    assert set(H) == set(pairs) == set(H_norm)
    direct_f = [Q(0)]*n
    direct_h = dict.fromkeys(pairs,Q(0))

    for w, barycentric, site in zip(weights,points,coordinates):
        R,Z = [sum(p*v[k] for p,v in zip(barycentric, geom.corners)) for k in (0,1)]
        x,y,d = site
        assert (R,Z) == (geom.R_center+geom.R_scale*x, geom.Z_center+geom.Z_scale*y)
        value = lambda p: S.T.evaluate(p,R,Z-geom.Z_center)
        q = [value(p) for p in original["q"]]
        B = [[value(p) for p in row] for row in original["B"]]
        assert q == [S.T.evaluate(p,x,y)/(d if k == 4 else 1) for k,p in enumerate(chart_q)]
        assert B == [[S.T.evaluate(p,x,y)/(d if k == 4 else 1) for p in row] for k,row in enumerate(chart_B)]

        def J(shifts):
            a,b,c,d,h = [v+sum(row[i]*s for i,s in shifts.items()) for v,row in zip(q,B)]
            return (a*d-b*c)*h

        base = J({})
        grad = [(8*(J({i:1})-J({i:-1})) - (J({i:2})-J({i:-2})))/12 for i in range(n)]
        hess = {}
        for i,j in pairs:
            if i == j:
                hess[i,j] = J({i:1}) - 2*base + J({i:-1})
            else:
                hess[i,j] = (J({i:1,j:1})-J({i:1,j:-1})-J({i:-1,j:1})+J({i:-1,j:-1}))/4
        e = base-1
        f1 = sum((-1)**(k+1)*(1+5000*sum(Q(1,j) for j in range(1,k+1)))*e**k for k in range(order+1))
        f2 = sum((-1)**k*(k+1)*(1+5000*sum(Q(1,j) for j in range(1,k+2)))*e**k for k in range(order+1))
        diagonal = (1,1+drive,1,1+drive,1)
        for i in range(n):
            assert abs(grad[i]) <= g_norm[i]
            direct_f[i] += w*(f1*grad[i] + sum(diagonal[k]*q[k]*B[k][i] for k in range(5)))
        for i,j in pairs:
            assert abs(hess[i,j]) <= H_norm[i,j]
            direct_h[i,j] += w*(f1*hess[i,j] + f2*grad[i]*grad[j] + sum(diagonal[k]*B[k][i]*B[k][j] for k in range(5)))
    assert any(i != j and value for (i,j),value in direct_h.items())
    for v,expect in zip(F,direct_f):
        exact_cover(S,v,expect)
    for key,value in H.items():
        exact_cover(S,value,direct_h[key])


def run():
    if not __debug__:
        raise RuntimeError("Review requires enabled assertions")
    bindings = {}

    def raw(name, expected=None):
        f = ROOT/name
        assert f.resolve().is_relative_to(ROOT)
        data = f.read_bytes()
        h = sha256(data).hexdigest()
        assert expected is None or expected == h, name
        assert bindings.setdefault(name,h) == h, name
        return data

    rel = JOB.relative_to(ROOT).as_posix()
    before = raw(f"{rel}/{SOURCE}",EXPECTED["source_sha256"])
    assert before == raw(f"{rel}/scaled_element21_frozen.py",EXPECTED["source_sha256"])
    p = json.loads(raw(f"{rel}/protocol.json",EXPECTED["protocol_sha256"]))
    pre = json.loads(raw(f"{rel}/preflight.json",EXPECTED["preflight_sha256"]))
    assert pre["status"] == "PASS_SCALED_CHART_PREPARATION_ONLY"
    assert len(pre["checks"]) == 16 and all(v is True for v in pre["checks"].values())
    assert pre["source_sha256"] == p["source_sha256"] == EXPECTED["source_sha256"]
    assert pre["protocol_sha256"] == EXPECTED["protocol_sha256"]
    assert pre["actual_record_force_H_integration_run"] is p["actual_record_force_H_integration_run"] is False
    for path,h in p["input_hashes"].items():
        raw(path,h)

    sys.path.insert(0,str(JOB))
    import scaled_element21 as S
    forbidden = AssertionError("actual moments, contraction, or cache sites entered during preparation replay")
    with patch.object(S,"material_moment",side_effect=forbidden), \
            patch.object(S,"_eval_parts",side_effect=forbidden), \
            patch.object(S.U,"uniform_cache",side_effect=forbidden), \
            patch.object(S.U.C,"exact_bulk_maps",side_effect=forbidden):
        hashes,c,case,old,rule,ids,nodes,center,physical,drive,g,Ft,Ht = S._certificate()
    assert hashes == p["input_hashes"] and c.record_id == p["record_id"] == "optoGEF2_small_col_L"
    assert p["element"] == case["element"] == old["element"] == 21
    assert p["unchanged_order"] == case["degree"] == old["degree"] == 12
    assert p["unchanged_rho_exact"] == case["rho_exact"] == old["rho_exact"]
    assert c.row["parameter_cell_sha256"] == p["same_parameter_cell_sha256"]
    assert p["global_indices"] == physical["global_indices"] == old["global_indices"]
    assert S.D._digest(S._metadata(g)) == S.D._digest(p["chart_geometry"])
    Rmin = min(r for r,z in g.corners)
    assert g.Rhat_min == Rmin/g.R_center > 0 and g.det_scaled*g.R_scale*g.Z_scale == g.det_physical > 0
    assert all(abs(x) <= 1 and abs(y) <= 1 for x,y in g.local_corners)
    assert len(ids) == 6 and len(Ft) == 12 and len(Ht) == 78
    assert 0 <= max(Ft) == Q(p["scaled_max_force_tail_exact"]) <= Q(p["unchanged_force_tail_budget_exact"])
    assert 0 <= max(Ht.values()) == Q(p["scaled_max_H_entry_tail_exact"]) <= Q(p["unchanged_H_entry_tail_budget_exact"])
    assert len(rule["volume_barycentric"]) == 64
    weights = rule["volume_weights"][case["element"]]
    assert len(weights) == 64 and min(weights) > 0
    assert weights == [Q.from_float(float(w))/c.target.energy_scale for w in c.target.a["volume_weights"][21]]
    for barycentric in rule["volume_barycentric"]:
        x,y,rhat = S.chart_point(g,barycentric)
        R,Z = [sum(v[k]*l for v,l in zip(g.corners,barycentric)) for k in (0,1)]
        assert g.R_center+g.R_scale*x == rhat*g.R_center == R
        assert g.Z_center+g.Z_scale*y == Z and rhat >= g.Rhat_min
    synth = closed_form_moments(S)
    synthetic_force_hessian(S,synth)
    uncertainty = arb(0,S.A(Q(1,16)))
    A,B = S.series_numerators({(1,0):uncertainty},S.T.ONE,8)
    assert A.keys() == B.keys() == {(j,0) for j in range(9)}
    assert not A[1,0].is_zero() and not B[1,0].is_zero()
    own = Path(__file__).relative_to(ROOT).as_posix()
    raw(own)
    for name,h in tuple(bindings.items()):
        raw(name,h)
    report = {
        **EXPECTED, "status": "PASS",
        "scope": "One active off-axis L element on the frozen original 640-bit box/C at order 12. Prelaunch review covers changes of basis and local force/Hessian enclosure only.",
        "checks": {
            "frozen_source_protocol_preflight_and_all_bound_assets_match": True,
            "same_L_original_box_C_free_projection_order12_Rho_and_normalization": True,
            "positive_radial_denominator_and_exact_physical_measure_with_no_extra_chart_factor": True,
            "all64_exact_same_sites_and_old_positive_binary_weights_without_site_sums": True,
            "signed_closed_form_log_moments_match_independent_slice_primitive": True,
            "independent_cubic_derivative_and_signed_point_oracle_for_all_synthetic_F_H_entries": True,
            "nonzero_mixed_H_and_distinct_denominator_powers_verified": True,
            "uncertain_crossing_zero_coefficients_retained": True,
            "same_positive_two_measure_bounds_and_DJ_D2J_tails_fit_original_budgets": True,
            "unchanged_physical_same_Y_cache_correction_reviewed": True,
        },
        "signed_log_fixtures": 6, "all64_sites_checked_without_sums": True,
        "actual_record_force_H_integration_run": False, "actual_cache_correction_run": False,
        "new_global_assembly_or_root_test": False, "input_hashes": dict(sorted(bindings.items())),
        "review_source_file": own, "review_source_sha256": bindings[own],
    }
    if "--freeze" in sys.argv:
        receipt = OUT/"replay.json"
        gate = JOB/"prelaunch_independent_review.json"
        assert not receipt.exists() and not gate.exists()
        with receipt.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(report,indent=2)+"\n")
        report["review_receipt_file"] = receipt.relative_to(ROOT).as_posix()
        report["review_receipt_sha256"] = sha256(receipt.read_bytes()).hexdigest()
        with gate.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":report["status"],"checks":len(report["checks"])}),flush=True)


if __name__ == "__main__":
    run()
