"""Review transfer-ambiguity prelaunch using frozen sources and bounded fixtures."""
import csv
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
PRE=HERE.parent
ROOT=HERE.parents[3]
SOURCE=ROOT/"scripts/transfer_ambiguity_linear.py"
PINS={
    SOURCE:"4e89db4c6e46f219f790408a45a71deb7fba4a641e7b359d7734e09b5082cb8b",
    PRE/"protocol.json":"ece0d458af81906c95faa231633254c3c69d55b3f0026074ccb8e9a0adeed314",
    PRE/"preflight.json":"23e1812067bf2a643ae133a1f7a05bfe1a17a93458f8cf10676ec0ae9ce29bce",
}
FIELD="runs/transfer_pair_pilot_v1/jobs/baseline_large/checkpoint_03/field.json"
SMALL="runs/material_cohort_descent_v1/protocol.json"
TANGENT="runs/material_cohort_descent_v1/results.json"
PILOT="runs/transfer_pair_pilot_v1/preparation/protocol.json"
SMALL_REVIEW="evidence/material_cohort_descent_independent_review.json"
LARGE_REVIEW="runs/transfer_pair_pilot_v1/independent_actual/actual_review.json"
SMALL_DATA="data/derived/calibration_small_18um.csv"
LARGE_INPUT="data/derived/spatial_inputs_only.json"
COLUMNS=["log_apical_area","log_basal_area","apical_tension_additive","log_basal_tension","log_common_gain"]
OBS=["maximum_shortening_um","region_averaged_shortening_um","apical_max_shortening_um","basal_max_shortening_um"]


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def write(path,payload):
    with path.open("x",encoding="utf-8") as f:
        json.dump(payload,f,indent=2)
        f.write("\n")


def exact_binary(x):
    return Q.from_float(x) if isinstance(x,float) else Q(x)


def direct_basis(t):
    return (2*(t-Q(1,2))*(t-1),4*t*(1-t),2*t*(t-Q(1,2)))


def slope(edge,t,k):
    ds=(4*t-3,4-8*t,4*t-1)
    return sum(n*pt[k] for n,pt in zip(ds,edge))


def projected_value(edge,t,k):
    return sum(n*pt[k] for n,pt in zip(direct_basis(t),edge))


def contains(ball,x):
    return ball.lo<=x<=ball.hi


def guard_on_binding_zero(m):
    # Rounded floating proposal violates an exact saturated inequality by a
    # small positive amount. The exact radial repair can safely return zero.
    fake=SimpleNamespace(success=True,message="",x=[0.5,2**-54],
                         ineqlin=SimpleNamespace(marginals=[0,0]))
    with patch.object(m,"linprog",return_value=fake):
        a=m.max_affine([1,0],[[0,1]],[1],1,[(-1,1),(0,1)])
    assert Q(a["repair_alpha_exact"])==Q(a["primal_lower_exact"])==0
    assert Q(a["dual_upper_exact"])==1


def guard_with_exact_rational_optimum(m):
    # Independently enumerable corners in box intersected with d0+d1<=1/2.
    vertices=[(-1,0),(-1,1),(Q(-1,2),1),(Q(1,2),0)]
    exact=max(2*x+y for x,y in vertices)
    a=m.max_affine([2,1],[[1,1]],[1],Q(3,2),[(-1,1),(0,1)])
    assert Q(a["primal_lower_exact"])==exact==Q(a["dual_upper_exact"])==1


def independent_actual_graph(state,review,protocol,m):
    # Extract by an edge-incidence map independent of the producer observer.
    ref=state["reference_nodes"]
    coords=state["current_nodes"]
    counts={}
    for triangle in state["triangles"]:
        a,b,c,ab,bc,ca=triangle
        for u,v,mid in ((a,b,ab),(b,c,bc),(c,a,ca)):
            key=tuple(sorted((u,v)))
            counts.setdefault(key,[]).append(mid)
    edge_ids=[]
    for (a,b),mids in counts.items():
        if len(mids)==1 and ref[a][1]==ref[b][1]==state["parameters"]["height"]:
            l,r=sorted((a,b),key=lambda i: ref[i][0])
            edge_ids.append((l,mids[0],r))
    edge_ids.sort(key=lambda ids:ref[ids[0]][0])
    L=Q("33.5")
    edges=[[(L*exact_binary(coords[n][0]),L*exact_binary(coords[n][1])) for n in ids] for ids in edge_ids]
    assert len(edges)==12 and edges[0][0][0]==0
    assert all(slope(edge,t,0)>0 for edge,t in product(edges,(0,1)))
    a=Q(protocol["shoulder_um_exact"])
    ((edge_index,edge),)=[(i,e) for i,e in enumerate(edges) if e[0][0]<a<e[2][0]]
    got=protocol["preparation_shoulder_details"]
    assert edge_index==got["edge_index"]==3
    lo,hi=map(Q,got["t_exact"])
    assert 0<lo<=hi<1 and projected_value(edge,lo,0)<=a<=projected_value(edge,hi,0)
    assert Q(got["endpoint_margin_um_exact"])==min(a-edge[0][0],edge[2][0]-a)
    assert got["Q_um_exact"]==m.state_gradient(state,a)[1]["Q_um_exact"]
    v=[tuple(map(Q,got["Q_um_exact"])),tuple(map(Q,review["Q_um"]["baseline_large"]))]
    assert max(x[0] for x in v)<=min(x[1] for x in v)
    gradient,detail=m.state_gradient(state,a)
    assert detail==got
    # Independent analytic derivative on a second, nonuniform-radius curve.
    u=Q(1,3)
    curve=[(t*t+2*t,t-3*t*t) for t in (Q(0),Q(1,2),Q(1))]
    radial=u*u+2*u
    local,_=m.contrast_curve_gradient([curve],radial)
    ratio=(1-6*u)/(2+2*u)
    for j,n in enumerate(direct_basis(u)):
        assert contains(local[0][j][0],ratio*n)
        assert contains(local[0][j][1],(1 if j==0 else 0)-n)
    # Check mapped support and outward derivative enclosures at the bracket.
    used=set(edge_ids[edge_index])|{edge_ids[0][0]}
    assert all(x.lo==x.hi==0 for i,row in enumerate(gradient) if i not in used for x in row)
    assert gradient[edge_ids[0][0]][0].lo==gradient[edge_ids[0][0]][0].hi==0
    assert gradient[edge_ids[0][0]][1].lo==gradient[edge_ids[0][0]][1].hi==1
    return {"edge_index":edge_index,"endpoint_margin_um":float(Q(got["endpoint_margin_um_exact"])),
            "top_boundary_edges":len(edges),"selected_edge_node_ids":list(edge_ids[edge_index]),
            "axis_top_node":edge_ids[0][0]}


def main():
    assert all(sha(p)==h for p,h in PINS.items())
    p,pre=[read(PRE/f) for f in ("protocol.json","preflight.json")]
    assert pre["status"]=="PASS_PREPARATION_ONLY" and pre["checks"] and all(v is True for v in pre["checks"].values())
    assert pre["protocol_sha256"]==PINS[PRE/"protocol.json"]
    assert p["source_sha256"]==pre["source_sha256"]==PINS[SOURCE]
    assert p["heldout_response_tables_read"] is p["new_equilibrium_solves"] is False
    assert not (PRE.parent/"actual").exists()
    pins=p["input_hashes"]
    assert all(sha(ROOT/k)==h for k,h in pins.items())
    assert {k for k in pins if not k.startswith("scripts/")}=={FIELD,SMALL,TANGENT,PILOT,SMALL_REVIEW,LARGE_REVIEW,"PLAN.md"}
    small,arch,pilot,f,sr,lr=[read(ROOT/k) for k in (SMALL,TANGENT,PILOT,FIELD,SMALL_REVIEW,LARGE_REVIEW)]
    assert sr["status"]=="PASS_INDEPENDENT_NUMERIC_COHORT_TANGENT_REVIEW"
    assert sr["reviewed_metadata_sha256"][SMALL]==pins[SMALL]==arch["protocol_sha256"]
    assert sr["reviewed_metadata_sha256"][TANGENT]==pins[TANGENT]
    assert lr["status"]=="PASS_ACTUAL_FOUR_SAVED_FIELD_REPLAY" and lr["all_replayed_hashes"][FIELD]==pins[FIELD]
    assert small["active_LP_columns"]==p["columns"]==sr["active_coordinates"]==COLUMNS
    assert small["active_column_indices"]==[1,2,3,4,5]
    ((a),)=[v for v in arch["candidates"] if v["variant"]=="current_disk__spatial_max"]
    lp=a["affine_LP"]
    labels=a["all68_constraint_order"]
    assert labels==[{"record_id":r["record_id"],"observable":v} for r,v in product(small["records"],OBS)]
    assert len(lp["affine_residual_intercepts_exact"])==len(lp["affine_derivative_coefficients_exact"])==len(labels)==68
    assert all(len(row)==5 for row in lp["affine_derivative_coefficients_exact"])
    assert p["trust_box_exact"]==lp["trust_box_exact"]
    assert [list(map(Q,x)) for x in p["trust_box_exact"]]==[list(map(exact_binary,x)) for x in small["trust_box"]]
    assert p["one_sided_apical_tension"] is True and Q(p["trust_box_exact"][2][0])==0
    assert sr["source_csv_sha256"]==sha(ROOT/SMALL_DATA)
    with (ROOT/SMALL_DATA).open(encoding="utf-8-sig",newline="") as stream:
        csvrows=list(csv.DictReader(stream))
    assert len(csvrows)==17 and all(row["record_id"]==rec["record_id"] for row,rec in zip(csvrows,small["records"],strict=True))
    sources={SMALL_DATA:sha(ROOT/SMALL_DATA),LARGE_INPUT:sha(ROOT/LARGE_INPUT)}
    for i,(row,record,archive) in enumerate(zip(csvrows,small["records"],arch["records"],strict=True)):
        assert archive["record_id"]==record["record_id"]==row["record_id"]
        assert all(Q(row[k])==Q(record["observed_exact"][k]) for k in OBS)
        fpath=ROOT/archive["file"]
        assert archive["sha256"]==sha(fpath)
        sources[archive["file"]]=archive["sha256"]
        packet=read(fpath)
        for j,key in enumerate(OBS):
            lo,hi=map(Q,packet["baseline_endpoint_intervals"][a["variant"]][key])
            expected=(lo+hi)/2-Q(row[key])
            assert Q(lp["affine_residual_intercepts_exact"][4*i+j])==expected
            ds=[d["two_steps"][1]["central_difference_intervals"][a["variant"]][key] for d in packet["coordinate_diagnostics"]]
            assert [Q(x) for x in lp["affine_derivative_coefficients_exact"][4*i+j]]==[Q.from_float(float(sum(map(Q,d))/2)) for d in ds]
    case=next(x for x in pilot["cases"] if x["case_id"]=="baseline_large")
    assert f["parameters"]==case["parameters"] and f["reached_target"] is True
    assert Q(case["gain_exact"])==Q(small["baseline_gain_exact"])==Q(sr["baseline_common_gain_exact"])
    assert [k for k in small["parameters"] if small["parameters"][k]!=case["parameters"][k]]==["patch_radius"]
    assert f["parameters"]["apical_surface_tension"]==0 and f["parameters"]["mu"]==1 and f["parameters"]["bulk"]==5000
    assert f["integration_policy"]["volume_duffy_gauss_order"]==8
    assert f["integration_policy"]["surface_gauss_order"]==12
    assert len(p["source_discrepancy_offsets_um_exact"])+len(p["small_prediction_change_budgets_um_exact"])==10
    assert min(Q(x) for x in p["source_discrepancy_offsets_um_exact"])>=0
    assert min(Q(x) for x in p["small_prediction_change_budgets_um_exact"])>0
    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_linear as m
    from sparse_quadratic import AxisymmetricModel
    from inspect import signature
    policy=signature(AxisymmetricModel.__init__).parameters["integration_policy"].default
    assert (policy.volume_duffy_gauss_order,policy.surface_gauss_order)==(8,12)
    graph=independent_actual_graph(f,lr,p,m)
    guard_on_binding_zero(m)
    guard_with_exact_rational_optimum(m)
    assert all(sha(ROOT/k)==h for k,h in pins.items())
    support={"status":"PASS_BOUNDED_PRELAUNCH_AUDIT", "initial_frozen_inputs":pins,
             "independent_fixture_and_row_inputs":sources,"graph":graph,
             "scope":"Source review, exact archived-row reconstruction, independent graph geometry/analytic fixture and small LP fixtures only. No new tangent, Hessian recovery, equilibrium or actual 20-LP diagnostic.",
             "rounding_scope":"Proven stored-graph derivative interval composed with floating nodal tangent, followed by an explicit binary midpoint objective and rounded small Jacobian. Exact affine proof does not certify nonlinear extrema.",
             "LP_math":"A=[J;-J], b=[delta-r;delta+r]. Feasible zero permits radial shrink; if any saturated face forces alpha=0, bound stays valid and gap reveals lost quality. For arbitrary nonnegative lambda, c·x <= lambda·b + h_box(c-A^T lambda). Independent rational corner fixture and a mock binding-zero solver exercise exact upper/lower checks.",
             "LP_scope":"The residual grid starts at the archived zero-step baseline discrepancy. The small-change grid sets the surrogate residual intercept to zero. Neither is an experimental tolerance. All 20 proposed optimizations retain exact widths and gaps.",
             "shared_premises":"Existing reviewed sparse force/tangent helpers and exact P2 graph kernels; no claim of a new verification of their mechanics.",
             "source_sha256":PINS[SOURCE],"protocol_sha256":PINS[PRE/"protocol.json"],
             "preflight_sha256":PINS[PRE/"preflight.json"],"audit_source_sha256":sha(Path(__file__))}
    out=HERE/"support.json"
    write(out,support)
    gate={"status":"PASS","checks":{
        "frozen_nonempty_preflight_and_allowlisted_inputs_rehashed":True,
        "all_68_archived_residual_and_five_coordinate_J_rows_rebuilt":True,
        "same_baseline_gain_materials_dose_case_and_default_8_12_policy":True,
        "fixed_axis_strict_interior_edge_and_moving_radius_derivative_verified":True,
        "one_sided_Ta_exact_primal_repair_and_weak_dual_bound_verified":True,
        "maxima_nonsmooth_and_rounded_affine_scope_explicit":True,
        "actual_diagnostic_not_run_during_prelaunch":True},
        "source_sha256":PINS[SOURCE],"protocol_sha256":PINS[PRE/"protocol.json"],
        "preflight_sha256":PINS[PRE/"preflight.json"],"support_file":out.relative_to(ROOT).as_posix(),
        "support_sha256":sha(out),"scope":support["scope"]}
    path=PRE/"prelaunch_independent_review.json"
    write(path,gate)
    print(sha(path),sha(out))


if __name__=="__main__":
    main()

