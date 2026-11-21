"""Cold independent synthetic checks for frozen transverse surface-stress scan."""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
P=ROOT/"runs/nonaxisymmetric_adequacy_audit_v1"
PINS={
    "surface_stress.py":"f33aaf277eea10fd364a529f95d369551c99b73747f7856cee3af4318b60b3a3",
    "protocol.json":"a09fb767d9990fce5816203d45579a5b400a1271248788b685a91aebda355e32",
    "preflight.json":"7eb3c81037d1ecac47f5e5ca45122f8abd0a316f5022314c14dd6e1470220206",
}


def h(path):
    return sha256(path.read_bytes()).hexdigest()


def fixed(x):
    return Q(x),Q(x)


def main():
    assert all(h(P/k)==v for k,v in PINS.items())
    spec=importlib.util.spec_from_file_location("_cold_surface_stress",P/"surface_stress.py")
    M=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(M)
    q=json.loads((P/"protocol.json").read_bytes())
    f=json.loads((P/"preflight.json").read_bytes())
    assert q["source_sha256"]==f["source_sha256"]==PINS["surface_stress.py"]
    assert h(P/"READBACK.md")==q["derivation_sha256"] and f["protocol_sha256"]==PINS["protocol.json"]
    assert q["actual_stress_evaluations"]==f["actual_stress_evaluations"]==0
    assert len(q["premise_sha256"])==26 and all(h(ROOT/k)==v for k,v in q["premise_sha256"].items())
    assert all(v is True for v in f["checks"].values())
    windows=[tuple(map(Q,w)) for w in q["fixed_windows_on_each_ordered_edge"]]
    assert windows==[(Q(i,16),Q(i+1,16)) for i in range(1,15)]
    assert all(0<x<y<1 for x,y in windows) and q["total_windows"]==2*2*12*len(windows)==672
    with patch.object(M,"evaluate",side_effect=AssertionError("no actual stress in cold review")):
        cases,notes,pins=M.inputs()
        assert M.protocol(notes,pins)==q
        assert [a["case_id"].rsplit("_",1)[-1] for a in cases]==["16","17"]
    # R(t)=1+t, even though the archived material middle still lies strictly
    # between the endpoints without being their arithmetic midpoint.
    ref=[Q(1),Q(71,50),Q(2)]
    radius=[Q(1),Q(3,2),Q(2)]
    trace=[[fixed(R/2),fixed(3*R/8)] for R in radius]
    win=(Q(1,4),Q(3,4))
    c,=M.screen_edge(ref,trace,Q(8),Q(2),(win,))
    assert c["standing"]=="STRICT_NEGATIVE_PATCH"
    assert Q(-7,2)<=Q(c["tau_negative_upper_exact"])<0
    assert 0<Q(c["minus_tau_r_R_squared_over_Js_lower_exact"])<=Q(14,5)
    assert tuple(map(Q,c["ref_radial_interval_exact"]))==(Q(5,4),Q(7,4))
    shape=[[fixed(R*R/2),fixed(0)] for R in radius]
    cross,=M.screen_edge(ref,shape,Q(4),Q(0),((Q(1,4),Q(1,2)),))
    # The true D changes sign across these closed endpoints.
    assert (Q(5,4)**2/2)<1<(Q(3,2)**2/2) and cross["standing"]=="UNRESOLVED"
    flat=[[fixed(R),fixed(0)] for R in radius]
    zero,=M.screen_edge(ref,flat,Q(8),Q(0),(win,))
    assert zero["standing"]=="NONNEGATIVE_ON_WINDOW" and zero["tau_negative_upper_exact"] is None
    large,=M.screen_edge(ref,trace,Q(8),Q(9),(win,))
    assert large["standing"]=="NONNEGATIVE_ON_WINDOW" and large["D_R2_times_Js2_minus_alpha2_exact"] is None
    d={"status":"PASS_COLD_INDEPENDENT_STRESS_AUTHORITY_AND_SYNTHETICS",
       "verified_current_pins_sha256":PINS,"input_closure_files":len(pins),"review_fixture_source_sha256":h(Path(__file__)),
       "checks":{"same_K_L_independent_posteriors_and_matched_whole_spatial_geometry":True,
                 "all_672_ordered_predeclared_nondegenerate_interior_windows":True,
                 "hand_known_negative_with_nonzero_slope_and_nonzero_axial_tangent":True,
                 "affine_reference_uses_corner_radii":True,
                 "true_sign_crossing_unresolved_and_exact_zero_or_Nonpositive_alpha_never_negative":True},
       "actual_stress_AD_global_matrix_or_PDE_jobs":0}
    out=HERE/"cold_source_and_fixture_readback.json"
    with out.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(d,f,indent=2)
        f.write("\n")
    print(h(out),flush=True)


if __name__=="__main__":
    main()
