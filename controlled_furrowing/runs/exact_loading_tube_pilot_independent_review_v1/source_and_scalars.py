"""Bind retained source and exact scalar endpoint diagnostics. No spatial or capsule evaluation."""
import base64
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
EX=ROOT/"runs/transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1"
ACTIVE=EX/"negative_KL_active_energy_v1"
D=ROOT/"runs/exact_loading_tube_design_v1"
PATHS={
"D/METHOD.md":D/"METHOD.md","D/crosswalk":D/"receipt.json","D/root_review":D/"root_readback.json",
"active memo":EX/"loading_affine_identity_v1/READBACK.md",
"active review":EX/"loading_affine_identity_v1/independent_readback/review_receipt.json",
"active proof":ACTIVE/"independent_actual/receipt.json","active true exit":ACTIVE/"independent_actual/os_exit_observation.json",
"active phase":ACTIVE/"protocol.json",
"pure body":"scripts/root_body_geometry.py","target":"scripts/exact_fe_target_safe.py",
"bulk surface":"scripts/exact_cohort_integrator.py","surface bounds":"scripts/exact_uniform_surface_bridge_B_safe.py",
"lower bulk":"scripts/exact_uniform_bulk_defect_B_safe.py","scaled bulk":"scripts/exact_scaled_cohort_bulk.py",
"explicit bridge":"scripts/transfer_ambiguity_exact_bridge.py","remaining verified context":"scripts/transfer_ambiguity_remaining32_exact_bridge.py",
"frozen point and cached variation":"scripts/algebraic_root_feasibility.py",
"prior tapered path":"scripts/energy_gain_cell.py","prior proof scope":"docs/ENERGY-GAIN-GEOMETRY.md",
"exact active derivative example":"scripts/exact_centered_bulk_L.py",
"existing independent local hull primitive":EX/"bridge_v1/independent_actual_preparation/replay_math.py",
}
PINS={"active proof":"84bbe79fac39fbf0ce72b45b9e0c18f76310cad981da004c5f0b5e2a7f6d1abf",
"active true exit":"060f31afa87c5571760ba0d57c4c51052bafbd2994d76a3328951c5dcef89316",
"active phase":"06dd35434358ad7dcb993349281ed6bf8a84fd0a34f7ca2d4932ec0792842361",
"D/crosswalk":"dbfb9871da2b384fcd57e44015a35bfeebbfccf25e23a6f45178f3aef488ff3c"}
LINKS={}


def check(p,expect=None):
    p=ROOT/p
    content=p.read_bytes()
    digest=sha256(content).hexdigest()
    assert expect is None or digest==expect
    LINKS[p.relative_to(ROOT).as_posix()]=digest
    return json.loads(content) if p.suffix==".json" else None


def polynomial_bounds(q,a,b,h):
    def positions(v,w):
        x=[Q(0),h]
        if w:
            rv=-v/(2*w)
            if 0<rv<h:
                x.append(rv)
        return x
    return (min(q[0]+a[0]*t+b[0]*t*t for t in positions(a[0],b[0])),
            max(q[1]+a[1]*t+b[1]*t*t for t in positions(a[1],b[1])))


def main():
    d={k:check(p,PINS.get(k)) for k,p in PATHS.items()}
    assert d["D/root_review"]["receipt_sha256"]==PINS["D/crosswalk"]
    assert d["active true exit"]["receipt_sha256"]==PINS["active proof"]
    assert d["active true exit"]["observed_exit_code"]==0 and d["active proof"]["certified_loading_path"] is False
    assert Q(d["active proof"]["difference"]["independent_interval_exact"][1])<0
    p=d["active phase"]
    hh={h:ROOT/k for k,h in p["inputs_sha256"].items()}
    k,l=p["cases"]
    K,L=(check(hh[x["exact_spatial_proof_sha256"]],x["exact_spatial_proof_sha256"]) for x in (k,l))
    u=[]
    for case in (k,l):
        cache=check(hh[case["cache_sha256"]],case["cache_sha256"])["frozen_binary_cache"]
        data=base64.b64decode(cache["data"],validate=True)
        assert sha256(data).hexdigest()==cache["sha256"]
        with np.load(BytesIO(data),allow_pickle=False) as arrays:
            u.append(list(map(Q,arrays["u"].tolist())))
    I,Ld=(Q(x["dose_exact"]) for x in (k,l))
    r=Q(K["radius_exact"])
    eta,z,old_margin=map(Q,(K["arithmetic"]["eta_upper_exact"],K["arithmetic"]["z_upper_exact"],
                           K["arithmetic"]["selfmap_margin_lower_exact"]))
    ek=eta/(1-z)
    el=Q(L["independent_bounds"]["posterior_error_upper_exact"])
    assert ek==Q(k["posterior_exact"])<r and el==Q(l["posterior_exact"])<r
    assert eta+z*r<r and old_margin==r-eta-z*r
    span=max(abs(q-p) for p,q in zip(*u,strict=True))
    h=Ld-I
    assert h==Q(22159027,10**8) and r==Q(1,2**34)
    candidate=Q(1,2**12)
    box=candidate+span/2
    # Synthetic signed quadratic endpoint/stationary-point only.
    assert polynomial_bounds((0,0),(-1,2),(1,0),Q(1))==(-Q(1,4),Q(2))
    assert polynomial_bounds((0,0),(0,2),(0,-2),Q(1))==(Q(0),Q(1,2))
    out={"status":"PASS_COLD_SCALARS_AND_CURRENT_SOURCE_BINDINGS","source_sha256":sha256(Path(__file__).read_bytes()).hexdigest(),
         "file_sha256":LINKS,
         "same_fixed_reference_K_L":True,"positive_normalization_exact":k["normalization_exact"],
         "dose_span_exact":str(h),"exact_max_trial_center_motion":str(span),"exact_max_secant_slope":str(span/h),
         "prospective_radius_exact":str(candidate),"maximum_rectangular_coordinate_radius_exact":str(box),
         "old_K_cube_radius_exact":str(r),"old_K_eta_exact":str(eta),"old_K_z_exact":str(z),
         "old_K_spatial_selfmap_headroom_exact":str(old_margin),"K_posterior_exact":str(ek),"L_posterior_exact":str(el),
         "display":{"max_center_span":float(span),"max_secant_slope":float(span/h),"max_rectangular_coordinate_radius":float(box),
                    "old_K_selfmap_headroom":float(old_margin),"K_posterior":float(ek),"L_posterior":float(el),"old_K_z":float(z)},
         "new_spatial_geometry_quadrature_matrix_AD_or_root_jobs":0,"new_loading_interval_certified":False}
    p=HERE/"source_and_scalar_checkpoint.json"
    with p.open("x",encoding="utf-8",newline="\n") as stream:
        json.dump(out,stream,indent=2)
        stream.write("\n")
    print(sha256(p.read_bytes()).hexdigest())


if __name__=="__main__":
    main()
