"""Frozen independent replay of saved nonlinear transfer fields."""
import argparse
from collections import Counter
from dataclasses import asdict
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
import math
import os
from pathlib import Path
import sys
from time import monotonic

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
ROOT=HERE.parents[2]
PRE=BASE/"preparation"
WORK=BASE/"actual"
PINS={
 ROOT/"scripts/transfer_ambiguity_nonlinear_probe.py":"87ca09b46f9569903276923d937bef63d00bba9ad1b99d78b7e1708168f7e631",
 PRE/"protocol.json":"9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94",
 PRE/"preflight.json":"f225616c4274a856576f587d40a93c89c64f88a9e4b3f1522d9f026434477211",
 PRE/"prelaunch_independent_review.json":"df6daf32ca582e03ce7d4ab2ecc97e4ce44ca6420fe85134f352b5a64aa21b26",
 PRE/"independent/support.json":"4e3794e035b5b49f7978c1df51fcbe3ffd41713886d852c6205633b67351a0cb",
 ROOT/"runs/transfer_pair_pilot_v1/independent_actual/replay.py":"6bd07976cb43528bdc148dd3a0c8d5e3eca57ff30b4c234f227e2a0180a34828"
}
ENDPOINTS=("max_contrast_change","max_negative_contrast_change")
OBS=("maximum_shortening_um","region_averaged_shortening_um","apical_max_shortening_um","basal_max_shortening_um")
SCOPE=("Verified stored-coordinate P2 fields and numerical Q8/Q12 stationarity samples. Shared reviewed force, "
       "whole-body geometry and exact graph supremum/thickness kernels; independent checkpoint controls, "
       "boundary incidence, disk means, component peaks and residual/comparison arithmetic. "
       "This does not enclose positions of true equilibria, integration error, a continuum PDE, continuous "
       "stability, source-compatible fitting or experimental predictions.")
L=Q(67,2)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def write(path,payload):
    with path.open("x",encoding="utf-8") as stream:
        json.dump(payload,stream,indent=2,allow_nan=False)
        stream.write("\n")


def rel(path):
    return path.relative_to(ROOT).as_posix()


def bound(path,h,mapping):
    k=rel(path)
    assert sha(path)==h==mapping.setdefault(k,h),k
    return read(path)


def yes(x):
    return isinstance(x,dict) and bool(x) and all(v is True for v in x.values())


def add(a,b):
    return a[0]+b[0],a[1]+b[1]


def minus(a,b):
    return a[0]-b[1],a[1]-b[0]


def mul(a,b):
    vals=[x*y for x,y in product(a,b)]
    return min(vals),max(vals)


def scale(a,x):
    return tuple(sorted((a[0]*x,a[1]*x)))


def overlap(a,b):
    return max(a[0],b[0])<=min(a[1],b[1])


def pack(x):
    return list(map(str,x))


def conv(a,b):
    out=[Q(0)]*(len(a)+len(b)-1)
    for i,j in product(range(len(a)),range(len(b))):
        out[i+j]+=a[i]*b[j]
    return out


def poly(p,t):
    return sum(v*t**i for i,v in enumerate(p))


def interval_horner(p,ab):
    out=(Q(0),Q(0))
    for v in reversed(p):
        out=add(mul(out,ab),(v,v))
    return out


def power(edge,k):
    a,m,b=(node[k] for node in edge)
    return a,4*m-3*a-b,2*(a-2*m+b)


def bracket(p,r,steps=80):
    if r==p[0]:
        return Q(0),Q(0)
    if r==sum(p):
        return Q(1),Q(1)
    assert p[0]<r<sum(p)
    lo,hi=Q(0),Q(1)
    for _ in range(steps):
        mid=(lo+hi)/2
        x=poly(p,mid)
        if x==r:
            return mid,mid
        if x<r:
            lo=mid
        else:
            hi=mid
    return lo,hi


def disk_and_extrema(edges,radius):
    """Independent exact antiderivative and single-surface extrema."""
    integral=(Q(0),Q(0))
    minima,maxima=[],[]
    assert edges[0][0][0]==0 and radius<=edges[-1][-1][0]
    for e in edges:
        rad,z=power(e,0),power(e,1)
        assert min(rad[1],rad[1]+2*rad[2])>0
        if e[0][0]>=radius:
            break
        cut=bracket(rad,min(radius,e[-1][0]))
        integrand=conv(conv(z,rad),[rad[1],2*rad[2]])
        primitive=[Q(0)]+[v/Q(i+1) for i,v in enumerate(integrand)]
        integral=add(integral,interval_horner(primitive,cut))
        va,vb=poly(z,cut[0]),poly(z,cut[1])
        certain=[z[0],min(va,vb),max(va,vb)] if cut[0]==cut[1] else [z[0],va]
        possible=[z[0],va,vb]
        if z[2]:
            critical=-z[1]/(2*z[2])
            if 0<critical<=cut[1]:
                possible.append(poly(z,critical))
                if critical<=cut[0]:
                    certain.append(poly(z,critical))
        minima.append((min(possible),min(certain)))
        maxima.append((max(certain),max(possible)))
    mass=radius**2/2
    return scale(integral,1/mass),tuple(map(min,zip(*minima))),tuple(map(max,zip(*maxima)))


def residual(y,x,delta):
    lo,hi=x
    assert lo<=hi
    d0=max(lo-y,y-hi,Q(0))
    d1=max(abs(lo-y),abs(hi-y))
    return d0,d1,("OUT" if d0>delta else ("IN" if d1<=delta else "UNRESOLVED"))


def gap(x,y):
    return x[0]-y[1],x[1]-y[0]


def prematches():
    hashes={}
    for path,h in PINS.items():
        assert sha(path)==h
        hashes[rel(path)]=h
    p,pre,g=[read(PRE/f) for f in ("protocol.json","preflight.json","prelaunch_independent_review.json")]
    assert pre["status"]==g["status"]=="PASS" and yes(pre["checks"]) and yes(g["checks"])
    assert p["source_sha256"]==g["source_sha256"]==PINS[ROOT/"scripts/transfer_ambiguity_nonlinear_probe.py"]
    assert g["protocol_sha256"]==PINS[PRE/"protocol.json"] and g["preflight_sha256"]==PINS[PRE/"preflight.json"]
    assert g["support_sha256"]==PINS[PRE/"independent/support.json"]
    assert pre["nonlinear_jobs_launched"]==0
    assert [x["endpoint_id"] for x in p["endpoints"]]==list(ENDPOINTS)
    for name,digest in {**p["source_hashes"],**p["input_hashes"]}.items():
        assert sha(ROOT/name)==digest==hashes.setdefault(name,digest)
    return p,hashes


def prepare():
    p,hashes=prematches()
    own=Path(__file__).resolve()
    frozen=HERE/"replay_frozen.py"
    with frozen.open("xb") as file:
        file.write(own.read_bytes())
    # Flat Z=2, R=t, and a variable-height case Z=R with disk mean 2a/3.
    for make,expected in ((lambda t:(t,Q(2)),Q(2)),(lambda t:(t,t),Q(1,2))):
        edge=[make(t) for t in (Q(0),Q(1,2),Q(1))]
        a,b,c=disk_and_extrema([edge],Q(3,4))
        assert a==(expected,expected)
        assert b[0]<=b[1]<=c[0]<=c[1]
    assert residual(Q(3),(Q(2),Q(4)),Q(1))==(0,1,"IN")
    assert residual(Q(3),(Q(3),Q(5)),Q(1))[2]=="UNRESOLVED"
    assert residual(Q(3),(Q(5),Q(6)),Q(1))[2]=="OUT"
    assert gap((Q(2),Q(3)),(Q(0),Q(1)))==(1,3)
    protocol={"status":"FROZEN_WAITING_FOR_EXTERNAL_TERMINAL_ADMISSION","source_sha256":sha(own),
              "prerequisite_sha256":hashes,"expected_endpoint_ids":list(ENDPOINTS),"fields_per_endpoint":22,
              "supervisor_pid":35436,"supervisor_unified_exec_session":90479,
              "required_terminal_shape":"Exact external root observation and 13 hashes. No field read before all terminal/results hashes and exits pass.",
              "scope":SCOPE}
    write(HERE/"protocol.json",protocol)
    preflight={"status":"PASS_INDEPENDENT_ACTUAL_REPLAY_PREPARATION_ONLY","source_sha256":sha(own),
               "protocol_sha256":sha(HERE/"protocol.json"),"checks":{
                  "current_frozen_producer_inputs_and_prelaunch_bindings":True,
                  "independent_polynomial_mean_and_scalar_residual_targeted_checks":True,
                  "actual_results_and_fields_unread":True,"no_numerical_jobs_run":True}}
    write(HERE/"preflight.json",preflight)
    print(json.dumps({"source":sha(own),"protocol":sha(HERE/"protocol.json"),"preflight":sha(HERE/"preflight.json")}),flush=True)


def admit(path,h,protocol,hashes):
    """All external/OS/worker receipts and both results before any field read."""
    admission=bound(path,h,hashes)
    assert admission["status"]=="PASS_ROOT_OBSERVED_NONLINEAR_TERMINAL"
    assert admission["supervisor_pid"]==35436 and admission["unified_exec_session"]==90479
    assert type(admission["observed_exit_code"]) is int and admission["observed_exit_code"]==0
    names=[WORK/f for f in ("started.json","terminal.json","results.json")]
    child=("supervisor_child_started.json","worker_started.json","supervisor_os_exit.json","terminal.json","results.json")
    names.extend(WORK/k/f for k,f in product(ENDPOINTS,child))
    assert set(admission["sha256"])=={rel(x) for x in names}
    frozen={x:bound(x,admission["sha256"][rel(x)],hashes) for x in names}
    start,t,result=[frozen[WORK/f] for f in ("started.json","terminal.json","results.json")]
    assert start["pid"]==t["pid"]==35436 and t["status"]=="FINISHED"
    assert start["bindings"]==t["bindings"]==result["bindings"]
    assert start["serial_endpoint_order"]==list(ENDPOINTS)
    assert t["results_sha256"]==admission["sha256"][rel(WORK/"results.json")]
    assert result["status"]=="COMPLETE_CONDITIONAL_STORED_FIELD_COMPARISON"
    assert result["scope"]==protocol["scope"]
    for k in ("source_witness_certified","exact_integral_FE_root_certified","continuum_accuracy_certified","heldout_response_tables_read"):
        assert result[k] is False
    paths={k:PRE/f for k,f in (("protocol_sha256","protocol.json"),("preflight_sha256","preflight.json"),("prelaunch_review_sha256","prelaunch_independent_review.json"))}
    binds={k:sha(f) for k,f in paths.items()}
    binds["source_sha256"]=protocol["source_sha256"]
    binds["linear_actual_review_sha256"]=protocol["authority"]["independent_linear_actual_review_sha256"]
    assert binds==start["bindings"]
    assert len(t["workers"])==len(result["endpoints"])==2
    assert [x["endpoint_id"] for x in t["workers"]]==list(ENDPOINTS)
    assert [x["endpoint_id"] for x in result["endpoints"]]==list(ENDPOINTS)
    prior=start["utc"]
    for k,linked,component in zip(ENDPOINTS,t["workers"],result["endpoints"],strict=True):
        folder=WORK/k
        s,worker,os,done,res=[frozen[folder/f] for f in child]
        assert s["pid"]==worker["pid"]==os["pid"]==done["pid"]==linked["pid"]
        assert s["supervisor_pid"]==worker["supervisor_pid"]==os["supervisor_pid"]==done["supervisor_pid"]==35436
        assert done["binds"]==worker["bindings"]==binds
        assert done["endpoint_id"]==worker["endpoint_id"]==res["endpoint_id"]==k
        assert done["status"]=="FINISHED" and type(os["observed_exit_code"]) is int and os["observed_exit_code"]==0
        assert os["terminal_sha256"]==admission["sha256"][rel(folder/"terminal.json")]
        assert done["results_sha256"]==admission["sha256"][rel(folder/"results.json")]
        assert linked=={**os,"os_exit_receipt_sha256":admission["sha256"][rel(folder/"supervisor_os_exit.json")]}
        assert res==component and res["status"]=="COMPLETE_CONDITIONAL_STORED_FIELDS"
        assert prior<=s["utc"]<=os["utc"]<=t["utc"] and worker["utc"]<=done["utc"]<=os["utc"]
        prior=os["utc"]
        assert len(res["checkpoints"])==22 and res["verified_final_outcome_fields"]==18
    return result


def calc_small(state,controls,graph):
    from quadratic_endpoint_observer import Curve, Field, field_supremum, positive_thickness, shared_radius_difference
    a,b=controls["apical"],controls["basal"]
    top,bottom=Curve(a,80),Curve(b,80)
    roi=(Q(0),L*Q(state["parameters"]["patch_radius"]))
    H=L*Q(state["parameters"]["height"])
    assert positive_thickness(top,bottom,roi,Q(1,1000),8192)>0
    u,al,_=disk_and_extrema(a,roi[1])
    d,_,bu=disk_and_extrema(b,roi[1])
    f=shared_radius_difference(top,bottom,H)
    if f is not None:
        # Same radial parametrization. Independent extremum on the exact
        # shared-height-difference polynomial as a scalar surface.
        paired=[]
        for u0,d0 in zip(a,b,strict=True):
            assert [p[0] for p in u0]==[p[0] for p in d0]
            paired.append([(x[0],H-x[1]+y[1]) for x,y in zip(u0,d0,strict=True)])
        _,_,mx=disk_and_extrema(paired,roi[1])
    else:
        # The reviewed two-graph joint bound preserves the sum at common R.
        joined=Field(H,[(-1,top),(1,bottom)])
        supremum=field_supremum(joined,roi,Q(1,1000),8192)
        assert supremum["requested_width_achieved"] is True
        mx=(supremum["interval"].lo,supremum["interval"].hi)
    raw={OBS[0]:mx,OBS[1]:add((H,H),minus(d,u)),OBS[2]:minus((H,H),al),OBS[3]:bu}
    assert all(lo<=hi for lo,hi in raw.values())
    return raw


def load_graph_review():
    import importlib.util
    path=ROOT/"runs/transfer_pair_pilot_v1/independent_actual/replay.py"
    spec=importlib.util.spec_from_file_location("reviewed_P2_readback",path)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def condition_book(endpoint,actual,original,p,hashes,graph):
    from sparse_quadratic import AxisymmetricModel as MainModel, IntegrationPolicy, SparsePolicy
    from robust_quadratic import Parameters, SolvePolicy
    from quadratic_quadrature_study import AxisymmetricModel as CrossModel, exact_geometry, evaluate
    import numpy as np
    import torch

    name=endpoint["endpoint_id"]
    folder=WORK/name
    policy=p["numerical_policy"]
    branch={}
    for level in ("small","large"):
        params=Parameters(**endpoint["materials"][level]["parameters"])
        solve=SolvePolicy(**policy["solve"])
        main=MainModel(params,solve,IntegrationPolicy(**policy["solve_quadrature"]),SparsePolicy(**policy["sparse"]))
        cross=CrossModel(params,solve,IntegrationPolicy(**policy["cross_force_quadrature"]))
        for key in ("reference","triangles","fixed","free","activation_elements"):
            assert np.array_equal(getattr(main,key),getattr(cross,key))
        branch[level]=(main,cross,params)
    gain=Q(endpoint["materials"]["small"]["gain_exact"])
    assert gain==Q(endpoint["materials"]["large"]["gain_exact"])
    plan=[("small_zero",endpoint["small_zero"],"small_zero")]
    plan.extend((f"small_{i:02d}",row,"small_record") for i,row in enumerate(endpoint["small_records"],1))
    plan.extend((f"large_{i:02d}",row,("large_zero","large_checkpoint","large_checkpoint","large_final")[i])
                for i,row in enumerate(endpoint["large_checkpoints"]))
    previous={}
    baseline={}
    delta=Q(p["selection"]["conditional_source_allowance_um_exact"])
    qvalue=None
    records,predictions,samples=[],{},[]
    for (label,step,role),row in zip(plan,actual["checkpoints"],strict=True):
        stage=folder/label
        level="small" if role.startswith("small") else "large"
        main,cross,params=branch[level]
        assert row["checkpoint"]==label and row["step"]==step and row["role"]==role
        request=bound(stage/"requested.json",sha(stage/"requested.json"),hashes)
        expect_previous=previous.get(level)
        prev=None if expect_previous is None else {k:expect_previous[k] for k in ("field","field_sha256","verification","verification_sha256")}
        assert request=={"level":level,"role":role,"step":step,"previous":prev}
        assert bound(stage/"passed.json",sha(stage/"passed.json"),hashes)==row
        assert row["field"]==rel(stage/"field.json") and row["verification"]==rel(stage/"verification.json")
        state=bound(stage/"field.json",row["field_sha256"],hashes)
        v=bound(stage/"verification.json",row["verification_sha256"],hashes)
        req=Q(step["dose_exact"])*gain
        reached=Q.from_float(float(req))
        assert req==Q(step["requested_strength_exact"]) and reached==Q(step["solver_strength_binary64_exact"])
        assert reached-req==Q(step["binary64_rounding_difference_exact"])
        assert state["parameters"]==asdict(params) and state["integration_policy"]==policy["solve_quadrature"]
        assert Q(state["activation_strength"])==Q(row["actual_strength_binary64_exact"])==reached
        assert reached/gain==Q(row["actual_dose_implied_by_binary_strength_exact"])
        assert state["reached_target"] is True and state["converged"] is True
        for k in ("reference","triangles","activation_elements"):
            assert np.array_equal(state[{"reference":"reference_nodes"}.get(k,k)],getattr(main,k))
        vector=np.asarray(state["free_displacement"],dtype="<f8")
        assert vector.shape==(611,) and np.isfinite(vector).all()
        coords=main.coordinates(torch.tensor(vector,dtype=torch.float64)).detach().numpy()
        assert np.array_equal(coords,state["current_nodes"])
        if role.endswith("zero"):
            assert expect_previous is None and reached==0
            assert state["zero_initialization"]=="Independent zero vector in this endpoint and geometry."
            assert row["attempts"]==row["rejected_attempts"]==0
            assert not row["accepted_path_implied_optical_doses"]
        else:
            assert expect_previous is not None
            assert bound(stage/"solver/final.json",sha(stage/"solver/final.json"),hashes)==state
            audit=state["robust_continuation"]
            assert audit["failed_state_reuse"] is False
            ledger=bound(stage/"solver/attempt_ledger.json",sha(stage/"solver/attempt_ledger.json"),hashes)
            assert ledger["attempts"]==audit["attempts"] and ledger["accepted"]==audit["accepted_increments"]
            assert Q(ledger["requested_strength"])==Q(ledger["last_verified_strength"])==reached
            present=Q(expect_previous["actual_strength_binary64_exact"])
            saved_hash=baseline[level]
            assert row["attempts"]==len(ledger["attempts"])
            assert row["rejected_attempts"]==sum(not x["accepted"] for x in ledger["attempts"])
            for attempt in ledger["attempts"]:
                assert Q(attempt["from_strength"])==present and attempt["initial_verified_vector_sha256"]==saved_hash
                target=Q(attempt["to_strength"])
                assert present<target<=reached
                if "result_path" in attempt:
                    log_path=Path(attempt["result_path"])
                    if not log_path.is_absolute():
                        log_path=ROOT/log_path
                    assert log_path.resolve().parent==stage/"solver"
                    packet=bound(log_path,sha(log_path),hashes)
                    assert packet["parameters"]==asdict(params)
                    assert Q(packet["activation_strength"])==target
                    ph=sha256(np.asarray(packet["free_displacement"],dtype="<f8").tobytes()).hexdigest()
                    assert ph==attempt["final_vector_sha256"]
                if attempt["accepted"]:
                    assert "result_path" in attempt
                    present,saved_hash=target,attempt["final_vector_sha256"]
            assert present==reached and saved_hash==sha256(vector.tobytes()).hexdigest()
            assert [(Q(x["from_strength"])/gain,Q(x["to_strength"])/gain) for x in ledger["accepted"]]==[
                (Q(y["from_dose_exact"]),Q(y["to_dose_exact"])) for y in row["accepted_path_implied_optical_doses"]]
        baseline[level]=sha256(vector.tobytes()).hexdigest()
        previous[level]=row
        u,c=evaluate(main,state),evaluate(cross,state)
        assert all(math.isfinite(x["scaled_force_linf"]) and 0<=x["scaled_force_linf"]<1e-7 for x in (u,c))
        certificate=exact_geometry(state)
        assert len(certificate["triangle_certificates"])==144
        assert certificate["all_exact_J_positive"] is True and certificate["global_embedding_certified"] is True
        assert Q(certificate["J_lower_exact"])>0
        if level=="small":
            ve=v["verification"]
            assert certificate==ve["geometry"] and v["variant"]=="current_disk__spatial_max"
            assert math.isclose(u["scaled_force_linf"],ve["Q8_force_linf"],rel_tol=5e-12,abs_tol=3e-14)
            graph.assert_mechanical_close(c,ve["Q12_cross_evaluation"])
            controls=graph.geometric_controls(state,L)
            values=calc_small(state,controls,graph)
            archived={k:tuple(map(Q,vs)) for k,vs in v["observer_intervals_um_exact"].items()}
            assert set(archived)==set(values)==set(OBS)
            assert all(overlap(values[k],archived[k]) for k in OBS)
            if role=="small_record":
                records.append(step)
                predictions[step["record_id"]]=archived
        else:
            assert v["status"]=="PASS" and yes(v["checks"]) and v["root_validated"] is False
            assert certificate==v["geometry"]
            graph.assert_mechanical_close(u,v["Q8_mechanical_evaluation"])
            graph.assert_mechanical_close(c,v["Q12_cross_evaluation"])
            shoulder=Q(p["large_observable"]["shoulder_radius_um_exact"])
            q,_,_,_=graph.replay_surface(state,v["surface"],L,shoulder,p["large_observable"]["root_bisections"],role=="large_final")
            if role=="large_final":
                qvalue=q
        samples.append({"checkpoint":label,"Q8":u["scaled_force_linf"],"Q12":c["scaled_force_linf"],
                        "J_lower_exact":certificate["J_lower_exact"],"attempts":row["attempts"],
                        "rejected_attempts":row["rejected_attempts"],"field_sha256":row["field_sha256"]})
        print(name,label,"PASS",flush=True)
    assert tuple(map(Q,actual["large_Q_um_exact"]))==qvalue
    summary=actual["small_residual"]
    all68=summary["all68_residuals"]
    assert all68["records"]==len(records)==17 and all68["constraints"]==68
    intervals=[]
    for (record,k),obs in zip(product(records,OBS),all68["constraint_rows"],strict=True):
        rid,y=record["record_id"],Q(record["observed_exact"][k])
        xi=predictions[rid][k]
        lower,upper,label=residual(y,xi,delta)
        assert obs["record_id"]==rid and obs["observable"]==k and Q(obs["observed_exact"])==y
        assert tuple(map(Q,obs["prediction_interval_exact"]))==xi
        assert lower==Q(obs["minimum_interval_residual_exact"]) and upper==Q(obs["maximum_interval_residual_exact"])
        assert obs["status"]==label and Q(obs["discrepancy_allowance_exact"])==delta
        intervals.append((k,lower,upper,label))
    maxima=(max(x[1] for x in intervals),max(x[2] for x in intervals))
    assert maxima==tuple(map(Q,summary["required_allowance_um_exact"]))
    for k in OBS:
        assert max(x[1] for x in intervals if x[0]==k)==Q(all68["discrepancy_lower_requirements_exact"][k])
        assert max(x[2] for x in intervals if x[0]==k)==Q(all68["discrepancy_sufficient_requirements_exact"][k])
    scalar="OUT" if maxima[0]>delta else ("IN" if maxima[1]<=delta else "UNRESOLVED")
    assert scalar==all68["status"]
    assert summary["status"]=={"OUT":"ALLOWANCE_EXCEEDED","IN":"ALLOWANCE_SATISFIED","UNRESOLVED":"UNRESOLVED"}[scalar]
    return {"endpoint_id":name,"checkpoints":samples,"large_Q_um_exact":pack(qvalue),
            "required_allowance_um_exact":pack(maxima),"conditional_stored_allowance_classification":scalar}


def replay(args):
    p,hashes=prematches()
    protocol=bound(HERE/"protocol.json",args.protocol_sha256,hashes)
    pre=bound(HERE/"preflight.json",args.preflight_sha256,hashes)
    assert pre["status"]=="PASS_INDEPENDENT_ACTUAL_REPLAY_PREPARATION_ONLY" and yes(pre["checks"])
    assert pre["protocol_sha256"]==args.protocol_sha256
    assert sha(Path(__file__))==sha(HERE/"replay_frozen.py")==protocol["source_sha256"]==pre["source_sha256"]
    assert protocol["prerequisite_sha256"]=={rel(f):sha(f) for f in PINS}|p["source_hashes"]|p["input_hashes"]
    hashes[rel(Path(__file__).resolve())]=sha(Path(__file__))
    hashes[rel(HERE/"replay_frozen.py")]=protocol["source_sha256"]
    result=admit(args.admission.resolve(),args.admission_sha256,p,hashes)
    start=monotonic()
    output=HERE/"actual_replay"
    output.mkdir(exist_ok=False)
    write(output/"started.json",{"pid":os.getpid(),"admission_sha256":args.admission_sha256})
    try:
        sys.path.insert(0,str(ROOT/"scripts"))
        import torch
        torch.set_num_threads(1)
        graph=load_graph_review()
        summary=[condition_book(x,a,protocol,p,hashes,graph) for x,a in zip(p["endpoints"],result["endpoints"],strict=True)]
        assert sum(len(x["checkpoints"]) for x in summary)==44
        raw=[tuple(map(Q,x["large_Q_um_exact"])) for x in summary]
        observed=gap(*raw)
        assert observed==tuple(map(Q,result["large_Q_positive_direction_minus_negative_direction_um_exact"]))
        assert result["large_Q_stored_intervals_disjoint"] is (observed[0]>0 or observed[1]<0)
        allowances=[tuple(map(Q,x["required_allowance_um_exact"])) for x in summary]
        assert tuple(map(Q,result["both_endpoints_required_allowance_um_exact"]))==tuple(map(max,zip(*allowances)))
        assert result["both_endpoints_satisfy_declared_allowance"] is all(x["conditional_stored_allowance_classification"]=="IN" for x in summary)
        assert result["verified_outcome_fields"]==36 and result["small_observer_intervals"]==136
        affine=Q(p["selection"]["half_pair_affine_large_Q_range_um_exact"])
        assert all(sha(ROOT/k)==h for k,h in hashes.items())
        result_out={"status":"PASS_INDEPENDENT_ACTUAL_NONLINEAR_STORED_FIELD_REPLAY","scope":SCOPE,
                    "actual_results_sha256":hashes[rel(WORK/"results.json")],
                    "control_scope":"Endpoint small records are loaded in common ascending-dose order from a fresh material-specific small zero. The large cases start at separate zero and use the frozen 1/3,2/3,1 dose checkpoints.",
                    "both_endpoints_satisfy_conditional_stored_allowance":result["both_endpoints_satisfy_declared_allowance"],
                    "endpoints":summary,"stored_large_Q_signed_difference_um_exact":pack(observed),
                    "half_step_affine_Q_difference_um_exact":str(affine),
                    "nonlinear_minus_affine_signed_difference_um_exact":pack(minus(observed,(affine,affine))),
                    "checkpoints":44,"numerical_force_evaluations":88,"stored_exact_element_checks":6336,
                    "independently_recomputed_small_residuals":136,"all_hashes":hashes,
                    "new_equilibrium_solves":False,"heldout_response_tables_read":False,"source_compatible_witness":False,
                    "exact_equilibrium_position_bounded":False,"seconds":monotonic()-start}
        write(output/"receipt.json",result_out)
        write(output/"terminal.json",{"status":"FINISHED","pid":os.getpid(),"exit_code":0,"receipt_sha256":sha(output/"receipt.json")})
        print(json.dumps({"status":result_out["status"],"receipt_sha256":sha(output/"receipt.json"),"elapsed_seconds":result_out["seconds"]}),flush=True)
    except BaseException as err:
        write(output/"terminal.json",{"status":"FAILED","pid":os.getpid(),"exit_code":1,"type":type(err).__name__,"error":str(err)})
        raise


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    sub=parser.add_subparsers(dest="stage",required=True)
    sub.add_parser("prepare")
    r=sub.add_parser("replay")
    r.add_argument("--admission",type=Path,required=True)
    r.add_argument("--admission-sha256",required=True)
    r.add_argument("--protocol-sha256",required=True)
    r.add_argument("--preflight-sha256",required=True)
    args=parser.parse_args()
    if args.stage=="prepare":
        prepare()
    else:
        replay(args)

