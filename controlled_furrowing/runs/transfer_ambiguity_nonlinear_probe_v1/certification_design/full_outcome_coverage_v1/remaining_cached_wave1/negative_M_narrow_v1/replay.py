"""Independent −M replay on a separately declared smaller cached cube."""
import argparse
import base64
from fractions import Fraction as Q
from hashlib import sha256
import io
import json
import os
from pathlib import Path
import sys
import time
import traceback

SOURCE=Path(__file__).resolve()
HERE=SOURCE.parent
WAVE=HERE.parent
DESIGN=WAVE.parents[1]
ROOT=DESIGN.parents[2]
OLD=WAVE/"independent_actual_preparation"
LAST=OLD/"actual_replay"
RUN=HERE/"actual"
P="max_negative_contrast_change_small_09"
R=Q(1,2**34)
N=611
COUNTS={"bulk":9216,"surface":288}
AD=ROOT/"scripts/review_algebraic_root.py"
PINS={
    AD:"6af2311d0667a28af263b6fb6237063481aa8e65209ae9722c2b0246a8cff27a",
    OLD/"replay.py":"1895a08aa52b3d007e871e8e4d0f15d032fcb9abf209228fd4e37137b5283a15",
    OLD/"protocol.json":"b04ccb595239fc9715f8e67b819725ede554cf7d602188bc025f2c95ed63315b",
    OLD/"preflight.json":"49e4a5ec1ec7202db731191285d1c08beb37c507627c250603520ae1ff192331",
    LAST/"receipt.json":"9854b5033a055ea33d7af818cf8da27b8eee35b3f4bb57606fcbf9bca96e3b56",
    LAST/"terminal.json":"745374287982959495edb20217fbf8d2a6828acaceae9b87b51d4f7f7b917fde",
    LAST/"os_exit_observation.json":"42e11a4b7065b07eca318e5127951e884e8bd7703773b37df184bc9b747706ef",
    LAST/f"{P}.json":"17f68e45c56df1a19bfa8233520e23a5ea6453128d3095780d017deb7b61214a",
    WAVE/"actual/root_terminal_admission.json":"c7b71addf9a035babe6f6bf5d91dae750dbc187fda83328a3f14265ed42b051f",
    WAVE/"actual/results.json":"f647b57696f69fc6260134a954f2880061803508c8eb3d2ed96dab46c3577284",
    WAVE/"actual/terminal.json":"381c2c20ff5670f2319d0b3b564c143432f2ca7c6bdf27aca0ca953727b342fa",
    WAVE/"preparation/protocol.json":"c1543cbbd2520664313b56114b66a885faaf4d75869f56dbb1aa83a43a78fe3a",
}
SCOPE=("Independent smaller-cube local root test for the same −M exact-real binary Q8/S12 potential. "
       "The initial r=2^-32 replay remains unresolved. No new model, R, nonlinear solve, exact spatial "
       "integration or stability, full-body, continuum or source conclusion.")


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def name(p):
    return p.relative_to(ROOT).as_posix()


def read(p,h):
    b=p.read_bytes()
    assert sha256(b).hexdigest()==h,p
    return json.loads(b)


def save(p,o):
    with p.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(o,f,indent=2,allow_nan=False)
        f.write("\n")
    return sha(p)


def phase(radius):
    assert type(radius) is Q and radius==R


def bindings():
    assert all(sha(p)==h for p,h in PINS.items())
    a,t,o,point=[read(p,PINS[p]) for p in (LAST/"receipt.json",LAST/"terminal.json",
                                          LAST/"os_exit_observation.json",LAST/f"{P}.json")]
    assert a["status"]=="INDEPENDENT_TWO_M_CACHED_REPLAY_WITH_UNRESOLVED_CASE"
    assert a["cached_local_roots"]==1 and a["spatial_root_count"]==0
    assert t["status"]=="FINISHED" and type(t["exit_code"]) is int and t["exit_code"]==0
    assert o["status"]=="PASS_OBSERVED_INDEPENDENT_M_REPLAY_TRUE_OS_EXIT"
    assert o["unified_exec_session"]==29833 and type(o["observed_exit_code"]) is int and o["observed_exit_code"]==0
    assert t["pid"]==o["pid"]==a["source_and_phase"]["pid"]==34020
    assert o["receipt_sha256"]==t["receipt_sha256"]==PINS[LAST/"receipt.json"]
    assert o["terminal_sha256"]==PINS[LAST/"terminal.json"]
    x,=[v for v in a["results"] if v["case_id"]==P]
    assert x["native_file"]==name(LAST/f"{P}.json") and x["native_sha256"]==PINS[LAST/f"{P}.json"]
    assert x["status"]==point["status"]=="UNRESOLVED_INDEPENDENT_M_CACHED_ROOT"
    assert point["independent_posterior_error_upper_exact"] is None
    assert point["producer_candidate_at_same_radius"] is True
    assert point["record_id"]=="optoGEF2_small_col_M" and point["case_id"]==P
    assert Q(point["radius_exact"])==Q(1,2**32)==4*R
    assert Q(point["independent_ZQ_upper_exact"])>1
    records=read(WAVE/"preparation/protocol.json",PINS[WAVE/"preparation/protocol.json"])
    c,=[v for v in records["cases"] if v["case_id"]==P]
    assert c["field"]==point["field"] and c["field_sha256"]==point["field_sha256"]
    goal=Q(c["e_goal_dimensionless_exact"])
    drive=Q(c["drive"]["common_gain_exact"])*Q(78208913,100000000)
    assert drive==Q(c["drive"]["activation_strength_exact"])==Q(point["exact_strength"])
    assert drive!=Q(c["drive"]["activation_strength_binary_exact"]) and point["rounded_strength_used_in_proof"] is False
    assert Q(point["observer_goal_exact"])==goal and 0<R<Q(1,1000)<goal
    packet=WAVE/"actual"/P/"cache_packet.json"
    assert a["producer_result_bytes_and_payloads"][name(packet)]==point["cache_packet_sha256"]
    assert a["producer_result_bytes_and_payloads"][point["field"]]==point["field_sha256"]
    assert sha(ROOT/point["field"])==point["field_sha256"]
    eta,z0=map(Q,(point["independent_eta_upper_exact"],point["independent_Z0_upper_exact"]))
    assert 0<=eta<R and 0<=z0<1
    rows=point["signed_center_force_intervals"]
    assert len(rows)==N
    assert max(max(abs(Q(x["lower_exact"])),abs(Q(x["upper_exact"]))) for x in rows)<=eta
    return c,point,packet


def plan():
    c,p,packet=bindings()
    return {"status":"FROZEN_INDEPENDENT_NEGATIVE_M_RADIUS_FALLBACK",
        "source_sha256":sha(SOURCE),"input_sha256":{name(f):h for f,h in PINS.items()},
        "phase":{"pilot_id":P,"radius_exact":str(R),"source_unresolved_radius_exact":p["radius_exact"],
           "site_counts":COUNTS,"precision_bits":256,"cached_packet":name(packet),
           "cache_packet_sha256":p["cache_packet_sha256"],"field":p["field"],"field_sha256":p["field_sha256"],
           "exact_strength":p["exact_strength"],"observer_goal_exact":p["observer_goal_exact"],
           "prior_eta_upper_exact":p["independent_eta_upper_exact"],"prior_center_Z0_upper_exact":p["independent_Z0_upper_exact"],
           "embedded_R_content_sha256":p["cache_tensor_content_sha256"]["preconditioner"]},
        "scope":SCOPE}


def prepare():
    p=plan()
    assert not RUN.exists()
    for wrong in (R/2,R*4,float(R),True):
        try:
            phase(wrong)
        except AssertionError:
            continue
        raise AssertionError("different radius or type admitted")
    phase(R)
    c,x,_=bindings()
    drive=Q(x["exact_strength"])
    rounded=Q(c["drive"]["activation_strength_binary_exact"])
    assert drive!=rounded
    g=save(HERE/"protocol.json",p)
    h=save(HERE/"preflight.json",{"status":"PASS","source_sha256":sha(SOURCE),"protocol_sha256":g,
            "reject_changed_r_or_nondeterministic_types":True,
            "exact_requested_drive_is_distinct_from_rounded_solver_drive":True,
            "input_cache_data_decoded":0,"new_model_or_AD_calculations":0})
    print(g,h,sha(SOURCE),flush=True)


def upper(ball):
    assert ball.is_finite()
    return Q(str(ball.upper().fmpq()))


def compute(c,proof,cache_path):
    import numpy as np
    sys.path.insert(0,str(ROOT/"scripts"))
    import review_algebraic_root as ad
    assert Path(ad.__file__).resolve()==AD
    ad.ctx.prec=256
    ad.A.cache_clear()
    packet=read(cache_path,proof["cache_packet_sha256"])
    assert packet["case_id"]==P and packet["record_id"]==proof["record_id"]
    assert packet["input_state_sha256"]==proof["field_sha256"] and packet["input_state"]==proof["field"]
    assert packet["closed_parameter_cell_exact"]==c["closed_parameter_cell_exact"]
    assert packet["parameters"]==c["material"]["parameters"]
    force_drive=Q(packet["activation_strength_exact"])
    assert force_drive==Q(proof["exact_strength"])==Q(c["step"]["requested_strength_exact"])
    assert force_drive!=Q(packet["activation_strength_binary_exact"])
    blob=packet["frozen_binary_cache"]
    data=base64.b64decode(blob["data"],validate=True)
    assert len(data)==blob["bytes"] and sha256(data).hexdigest()==blob["sha256"]==proof["cache_array_contents_sha256"]
    with np.load(io.BytesIO(data),allow_pickle=False) as arch:
        arrays={k:arch[k].copy() for k in arch.files}
    assert {k:sha256(v.tobytes(order="C")).hexdigest() for k,v in arrays.items()}==proof["cache_tensor_content_sha256"]
    assert all(v.dtype.kind in ("i","u","f") and np.isfinite(v).all() for v in arrays.values())
    assert {k:list(v.shape) for k,v in arrays.items()}==blob["array_shapes"]
    assert arrays["reference"].shape==(325,2) and arrays["triangles"].shape==(144,6)
    assert arrays["free"].tolist()==c["free_map"] and arrays["u"].shape==(N,)
    assert arrays["preconditioner"].dtype==np.dtype("<f8") and arrays["preconditioner"].shape==(N,N)
    eta,z0=map(Q,(proof["independent_eta_upper_exact"],proof["independent_Z0_upper_exact"]))
    assert eta<R
    centers=[]
    counts=dict.fromkeys(COUNTS,0)
    for term in ad.terms_from_arrays(arrays,packet["parameters"],force_drive):
        assert term["weight"]>0 and term["indices"]
        fn=ad.bulk_AD if term["kind"]=="bulk" else ad.surface_AD
        term["Hcenter"]=fn(term["q"],*term["parameters"])[1]
        term["absB"]=[[ad.UP(x) for x in row] for row in term["B"]]
        term["rowsumB"]=[sum(row,ad.arb(0)).upper() for row in term["absB"]]
        centers.append(term)
        counts[term["kind"]]+=1
    assert counts==COUNTS and len(centers)==9504
    rows,domain=ad.variation(centers,N,ad.A(R))
    assert domain["counts"]==COUNTS
    assert all(Q(v["lower_exact"])>0 for v in domain["strict_lower_margins"].values())
    assert Q(domain["strict_lower_margins"]["bulk_J_min"]["lower_exact"])>Q(1,10**6)
    qrows=[upper(rows[j,0]) for j in range(N)]
    assert min(qrows)>=0
    M=[[abs(Q(x)) for x in row] for row in arrays["preconditioner"].tolist()]
    norm=max(sum(row) for row in M)
    assert norm==Q(proof["preconditioner_infinity_norm_exact"])>0
    out=[sum(a*b for a,b in zip(row,qrows,strict=True)) for row in M]
    dz=max(out)
    z=z0+dz
    lhs=eta+z*R
    assert z<1 and lhs<R<Q(1,1000)<Q(proof["observer_goal_exact"])
    posterior=eta/(1-z)
    assert posterior<R
    return {"case_id":P,"record_id":proof["record_id"],"scope":SCOPE,"precision_bits":256,
        "radius_exact":str(R),"observer_goal_exact":proof["observer_goal_exact"],"exact_strength":str(force_drive),
        "field":proof["field"],"field_sha256":proof["field_sha256"],"cache_packet_sha256":proof["cache_packet_sha256"],
        "cache_array_contents_sha256":proof["cache_array_contents_sha256"],
        "preconditioner_binary_sha256":proof["cache_tensor_content_sha256"]["preconditioner"],
        "preconditioner_infinity_norm_exact":str(norm),"independent_whole_cube_domains":domain,
        "reused_signed_cached_anchor_file":name(LAST/f"{P}.json"),"reused_signed_cached_anchor_sha256":PINS[LAST/f"{P}.json"],
        "eta_upper_exact":str(eta),"center_Z0_upper_exact":str(z0),"independent_variation_upper_exact":str(dz),
        "independent_ZQ_upper_exact":str(z),"strict_selfmap_lhs_upper_exact":str(lhs),
        "posterior_error_upper_exact":str(posterior),"original_unresolved_radius_exact":proof["radius_exact"],
        "original_unresolved_test_preserved":True,"checks":{
            "same_center_exact_drive_binary_cache_and_R":True,"fresh_9504_site_whole_box_variation":True,
            "all_cached_analytic_domains_positive":True,"strict_selfmap_and_observer_goal":True},
        "display":{"eta":float(eta),"z":float(z),"lhs":float(lhs),"posterior":float(posterior)},
        "whole_rows_exact":list(map(str,qrows)),"mapped_rows_exact":list(map(str,out)),
        "new_models":0,"new_preconditioners":0,"new_equilibrium_solves":0,"spatial_integration_or_root_claim":False}


def run(ph,fh):
    p=read(HERE/"protocol.json",ph)
    pre=read(HERE/"preflight.json",fh)
    assert p==plan() and pre["status"]=="PASS" and pre["source_sha256"]==sha(SOURCE)
    assert pre["protocol_sha256"]==ph and pre["input_cache_data_decoded"]==pre["new_model_or_AD_calculations"]==0
    c,proof,path=bindings()
    RUN.mkdir(exist_ok=False)
    stamp={"pid":os.getpid(),"source_sha256":sha(SOURCE),"protocol_sha256":ph,"preflight_sha256":fh}
    save(RUN/"started.json",stamp|{"status":"RUNNING"})
    end=stamp|{"status":"FAILED","exit_code":1}
    start=time.perf_counter()
    try:
        x=compute(c,proof,path)
        assert p==plan() and sha(path)==proof["cache_packet_sha256"]
        h=save(RUN/"receipt.json",x|stamp|{"status":"PASS_INDEPENDENT_NEGATIVE_M_SMALLER_CACHED_ROOT",
                "seconds":time.perf_counter()-start})
        end|={"status":"FINISHED","exit_code":0,"receipt_sha256":h}
        print(json.dumps(x["display"]))
    except BaseException as e:
        save(RUN/"failure.json",{"error":str(e),"traceback":traceback.format_exc()})
        raise
    finally:
        save(RUN/"terminal.json",end)
        print(json.dumps(end),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    sub=parser.add_subparsers(dest="step",required=True)
    sub.add_parser("prepare")
    a=sub.add_parser("run")
    a.add_argument("--protocol-sha256",required=True)
    a.add_argument("--preflight-sha256",required=True)
    args=parser.parse_args()
    prepare() if args.step=="prepare" else run(args.protocol_sha256,args.preflight_sha256)
