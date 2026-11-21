"""Cold prepare and separately gated exact symmetric-part proofs for four cached R matrices."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback
from unittest.mock import patch

import numpy as np
from flint import ctx as arb_context

import cohort_energy_stability as producer
import review_cohort_energy_stability as replay
import transfer_ambiguity_cached_derivative_extension as DX
import transfer_ambiguity_exact_bridge as bridge
import transfer_ambiguity_exact_targets as EX

SOURCE = Path(__file__).resolve()
ROOT = SOURCE.parents[1]
OUT = EX.OUT / "preconditioner_spd_v1"
PREP = OUT / "preparation"
ACTUAL = OUT / "actual"
DESIGN = EX.OUT / "minimum_design_v1"
METHOD = OUT / "METHOD.md"
DESIGN_SHA = "fb7f36d8f88fad7781baff8d1d11c8aee481c50f4eb6a386444a347ecc1551b0"
DESIGN_MEMO_SHA = "d78c9fa8082c4ded27de1b9d85e420dc5f01c7801e815d99e5339db463155c80"
GEOMETRY_SHA = "884e0bc76a27465c86b62a202cab491d764aedf030623592238e09adee86d77d"
ROOT_REVIEW_SHA = "4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17"
MANIFEST_SHA = "cb25fbeaef753cea700c187cd8b377b7dcf4aa1fe2b0b246ed17134a07a89bc3"
REVIEW_STATUS = "PASS_INDEPENDENT_SINGLETON_PRECONDITIONER_SPD_PRELAUNCH"
PRECISION = 640
NFREE = 611
PILOTS = tuple(EX.PILOTS)
SCOPE = ("Exact SPD of the symmetric part of the frozen binary preconditioner, one matrix per approved "
         "fixed transfer target. An inertia premise only. No exact-spatial H calculation, root, minimum, "
         "stability, continuum, source fit or loading-branch claim.")


def rel(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def require(value,message):
    if not value:
        raise ValueError(message)


def yes(checks):
    return type(checks) is dict and bool(checks) and all(v is True for v in checks.values())


def true_zero(value):
    return type(value) is int and value == 0


def now():
    return datetime.now(timezone.utc).isoformat()


def save(path,value):
    with Path(path).open("x",encoding="utf8",newline="\n") as f:
        json.dump(value,f,indent=2,allow_nan=False)
        f.write("\n")
    return sha(path)


def pinned(path,digest):
    require(type(digest) is str and len(digest)==64 and all(c in "0123456789abcdef" for c in digest),
            "explicit lowercase SHA-256 required")
    raw = Path(path).read_bytes()
    require(sha256(raw).hexdigest()==digest,f"changed frozen input: {path}")
    return json.loads(raw)


def link(links,name,h):
    path = (ROOT/name).resolve()
    require(type(name) is str and path.is_relative_to(ROOT) and rel(path)==name,
            "canonical in-study source or artifact required")
    require(sha(path)==h and links.setdefault(name,h)==h,f"changed or conflicting input: {name}")


def boot():
    EX.bootstrap()
    require(arb_context.prec==PRECISION and DX.EX is bridge.EX is EX and DX is bridge.DX
            and EX.T.core.ctx is arb_context and producer.ctx is replay.ctx is arb_context,
            "same scoped exact target and shared Arb640 context required")


def selected(pilot,cases):
    require(type(pilot) is str and pilot in PILOTS,"unsupported target identity")
    matched=[case for case in cases if case["pilot_id"]==pilot]
    require(len(matched)==1,"one explicit current pilot record required")
    return matched[0]


def find_matrix(a,rank=NFREE):
    require(type(a) is np.ndarray and type(rank) is int and rank>0 and a.shape==(rank,rank)
            and a.dtype.str=="<f8" and np.isfinite(a).all(),
            "finite square little-endian binary64 matrix with declared rank required")
    return sha256(a.tobytes(order="C")).hexdigest()


def matched_target(expected,observed):
    require(observed==expected,"different pilot, cache, exact drive or binary R identity")


def authority():
    boot()
    links,cases,manifest,data=DX.authority()
    links=dict(links)
    design=pinned(DESIGN/"receipt.json",DESIGN_SHA)
    require(design["status"]=="READ_ONLY_SAME_TARGET_MINIMUM_DESIGN_COMPLETE"
            and design["memo_sha256"]==DESIGN_MEMO_SHA and design["actual_force_Hessian_AD_matrix_or_solve_jobs"]==0
            and design["no_stage_gate_or_root_or_minimum_claimed"] is True,
            "read-only minimum design and precise scope required")
    for path,h in design["source_and_metadata_sha256"].items():
        link(links,path,h)
    for path,h in bridge.PINS.items():
        link(links,path,h)
    for path,h in (
        (DESIGN/"receipt.json",DESIGN_SHA),(DESIGN/"README.md",DESIGN_MEMO_SHA),
        (EX.OUT/"preparation/target_manifest.json",MANIFEST_SHA),
        (ROOT/DX.GEOMETRY_REVIEW,GEOMETRY_SHA),(ROOT/EX.CACHE_REVIEW,ROOT_REVIEW_SHA),
    ):
        link(links,rel(path),h)
    links[rel(SOURCE)]=sha(SOURCE)
    links[rel(METHOD)]=sha(METHOD)
    for mod in (EX,DX,bridge,producer,replay):
        path=Path(mod.__file__).resolve()
        require(links.get(rel(path))==sha(path),"current runtime source alias differs")
    require(tuple(case["pilot_id"] for case in cases)==tuple(manifest["pilot_order"])==PILOTS,
            "all four approved singleton targets required")
    derivatives=bridge.derivative_inputs(bridge.DERIV_ADMISSION,bridge.DERIV_ADMISSION_SHA,cases)
    rows=[]
    for case in cases:
        pid=case["pilot_id"]
        entry=selected(pid,manifest["targets"])
        contract=pinned(ROOT/entry["model_contract"],entry["model_contract_sha256"])
        pt=pinned(ROOT/case["verified_point"],case["verified_point_sha256"])
        _,_,norm=DX.validate_center(pt,case)
        require(case["cache_packet_sha256"]==contract["point_cache_sha256"]
                and case["embedded_cache_sha256"]==contract["frozen_binary_cache_sha256"]
                and case["preconditioner_binary_sha256"]==contract["binary_preconditioner_sha256"]
                and case["verified_point_sha256"]==contract["cached_root_independent_receipt_sha256"]
                and case["model_contract_sha256"]==entry["model_contract_sha256"]
                and case["box_sha256"]==entry["box_1_sha256"],
                "same verified derivative and frozen preconditioner authority required")
        params=contract["exact_requested_drive"]
        dose,gain,exact=(Q(params[k]) for k in ("optical_dose_exact","common_gain_exact","activation_strength_exact"))
        require(0<exact==dose*gain==Q(case["actual_exact_drive"])==Q(pt["exact_strength"])
                and exact!=Q(params["activation_strength_binary_exact"])
                and Q(contract["energy_scale_binary_exact"])>0
                and contract["parameter_cell_sha256"]==EX.E._digest(contract["parameter_cell"]),
                "same exact singleton drive, cell and energy normalization required")
        identity={"pilot_id":pid,"record_id":case["record_id"],
                  "field_sha256":case["field_sha256"],"cache_packet_sha256":case["cache_packet_sha256"],
                  "embedded_cache_sha256":case["embedded_cache_sha256"],
                  "R_binary_sha256":case["preconditioner_binary_sha256"],
                  "model_contract_sha256":entry["model_contract_sha256"],
                  "parameter_cell_sha256":contract["parameter_cell_sha256"],
                  "dose_exact":str(dose),"gain_exact":str(gain),"drive_exact":str(exact),
                  "energy_scale_binary_exact":str(Q(contract["energy_scale_binary_exact"]))}
        rows.append({"pilot_id":pid,"identity":identity,"point_proof_sha256":case["verified_point_sha256"],
                     "cached_reference_R_norm_upper_exact":str(norm),
                     "larger_radius_box_sha256":case["box_sha256"],
                     "native_cached_derivative_case_sha256":derivatives[pid]["native_case_sha256"]})
    return dict(sorted(links.items())),rows,manifest,data


def protocol(state):
    links,rows,_,_=state
    return {"status":"FROZEN_FOUR_SINGLETON_SYMMETRIC_PRECONDITIONER_PART_PROOFS",
            "source_sha256":sha(SOURCE),"scope":SCOPE,"precision_bits":PRECISION,"rank":NFREE,
            "source_and_input_sha256":links,"ordered_targets":rows,
            "target_manifest_sha256":MANIFEST_SHA,"minimum_design_sha256":DESIGN_SHA,
            "cached_root_review_sha256":ROOT_REVIEW_SHA,"geometry_review_sha256":GEOMETRY_SHA,
            "larger_derivative_admission_sha256":bridge.DERIV_ADMISSION_SHA,
            "versions":{x:importlib.metadata.version(x) for x in ("numpy","scipy","python-flint")},
            "python_version":sys.version.split()[0],
            "exact_matrix":"P=(R+R^T)/2 using exact binary64 entries and outward Arb640 products.",
            "proposal":"K is an untrusted binary64 inverse Cholesky proposal; preserve its NPY and entry-byte hashes.",
            "acceptance":"Both exact midpoint-free outward congruence bounds on the same R,K must be strictly <1.",
            "runtime":"Four serial subprocesses. Match parent PID, source gate, child terminal and OS exits.",
            "actual_production_matrix_jobs_before_launch":0,
            "matrix_phase_cannot_promote_root_or_minimum":True}


def gate(review_path,review_sha):
    current={"source_sha256":sha(SOURCE),"protocol_sha256":sha(PREP/"protocol.json"),
             "preflight_sha256":sha(PREP/"preflight.json"),"minimum_design_sha256":DESIGN_SHA,
             "cached_root_review_sha256":ROOT_REVIEW_SHA,"target_manifest_sha256":MANIFEST_SHA}
    review=pinned(review_path,review_sha)
    require(review["status"]==REVIEW_STATUS and yes(review["checks"])
            and all(review[k]==v for k,v in current.items()),
            "current independently reviewed Sym(R) launch approval required")
    pf=pinned(PREP/"preflight.json",current["preflight_sha256"])
    require(pf["status"]=="PASS" and yes(pf["checks"])
            and pf["source_sha256"]==current["source_sha256"] and pf["protocol_sha256"]==current["protocol_sha256"]
            and true_zero(pf["actual_production_matrix_jobs"])
            and sha(PREP/"transfer_ambiguity_preconditioner_spd_frozen.py")==current["source_sha256"],
            "current cold matrix preflight required")
    snapshot=authority()
    require(protocol(snapshot)==pinned(PREP/"protocol.json",current["protocol_sha256"]),
            "current source, authority or runtime drift")
    return snapshot,{**current,"prelaunch_review_sha256":review_sha}


def actual_matrix(pilot,snapshot):
    _,cases,manifest,data=snapshot
    row=selected(pilot,cases)
    expected=row["identity"]
    # Reset shared precision and A cache before each factory, even in isolated workers.
    boot()
    ctx,_,entry=EX.open_context(pilot,manifest,data)
    require(arb_context.prec==PRECISION,"factory altered Arb precision")
    R=ctx.R
    require(R.flags.writeable is False and ctx.u.shape==ctx.free.shape==(NFREE,),
            "immutable selected exact matrix and matched 611-coordinate center required")
    got={"pilot_id":ctx.pilot_id,"record_id":ctx.record_id,
         "field_sha256":ctx.target.packet["input_state_sha256"],"cache_packet_sha256":ctx.row["point_cache_sha256"],
         "embedded_cache_sha256":ctx.row["frozen_binary_cache_sha256"],"R_binary_sha256":find_matrix(R),
         "model_contract_sha256":entry["model_contract_sha256"],
         "parameter_cell_sha256":ctx.row["parameter_cell_sha256"],
         "dose_exact":str(ctx.dose),"gain_exact":str(ctx.base["common_gain"]),"drive_exact":str(ctx.drive),
         "energy_scale_binary_exact":str(ctx.target.energy_scale)}
    matched_target(expected,got)
    return R,row


def check_child(observed,terminal,identity,manifest_sha):
    require(type(observed["pid"]) is int and observed["pid"]>0
            and type(observed["supervisor_pid"]) is int and observed["supervisor_pid"]>0
            and all(terminal[k]==v for k,v in identity.items())
            and all(observed[k]==identity[k] for k in ("pid","supervisor_pid","pilot_id"))
            and true_zero(observed["observed_exit_code"]) and true_zero(terminal["exit_code"])
            and terminal["status"]=="FINISHED" and terminal["result_sha256"]==manifest_sha,
            "real successful OS completion with matched child/phase required")


def worker(args):
    state,bindings=gate(args.review,args.review_sha256)
    parent=json.loads((ACTUAL/"started.json").read_bytes())
    require(parent=={"status":"RUNNING","pid":args.supervisor_pid,"bindings":bindings,
                     "ordered_pilots":list(PILOTS)} and os.getppid()==args.supervisor_pid,
            "live approved parent process required")
    folder=ACTUAL/args.pilot_id
    ident={"pilot_id":args.pilot_id,"pid":os.getpid(),"supervisor_pid":args.supervisor_pid,"bindings":bindings}
    save(folder/"started.json",{"status":"RUNNING",**ident})
    terminal={"status":"FAILED","exit_code":1,**ident}
    try:
        R,row=actual_matrix(args.pilot_id,state)
        r_hash=row["identity"]["R_binary_sha256"]
        p,K=producer.certify_sympart(R,precision_bits=PRECISION)
        require(arb_context.prec==PRECISION and p.get("R_binary_sha256",r_hash)==r_hash
                and find_matrix(R)==r_hash,"changed matrix/precision during proposal or verified product")
        replay_result=None
        k_links={}
        if K is not None:
            k_hash=find_matrix(K)
            file=folder/"congruence_K.npy"
            with file.open("xb") as f:
                np.save(f,K,allow_pickle=False)
            k_links={"K_file":rel(file),"K_file_sha256":sha(file),"K_binary_sha256":k_hash}
            frozen=np.load(file,allow_pickle=False)
            require(find_matrix(frozen)==k_hash,"K binary proposal serialization changed")
            q,rows=replay.congruence_bounds(R,frozen,PRECISION)
            norm=replay.norm_R_integer(R)
            require(arb_context.prec==PRECISION and len(rows)==NFREE and min(rows)>=0 and q==max(rows)
                    and 0<norm<=Q(row["cached_reference_R_norm_upper_exact"])
                    and find_matrix(R)==r_hash and sha(file)==k_links["K_file_sha256"],
                    "invalid second association, exact R norm or unchanged proposal")
            qr=Q(p["q_infinity_upper_exact"])
            require(len(p["row_absolute_sum_upper_exact"])==NFREE
                    and qr==max(map(Q,p["row_absolute_sum_upper_exact"])) and qr>=0,
                    "complete producer congruence row bounds required")
            replay_result={"q_infinity_upper_exact":str(q),"row_absolute_sum_upper_exact":list(map(str,rows)),
                           "exact_R_infinity_norm":str(norm),"independent_association":"K*(P*K^T)",
                           "source_sha256":sha(replay.__file__)}
        accepted=(replay_result is not None and p["status"]=="CERTIFIED_EXACT_SYMMETRIC_PART_SPD"
                  and Q(replay_result["q_infinity_upper_exact"])<1)
        result={"status":"CANDIDATE_EXACT_SYMMETRIC_PART_SPD_AWAITING_REVIEW" if accepted else
                         "UNRESOLVED_SYMMETRIC_PART_SUFFICIENT_TEST",
                **ident,"scope":SCOPE,"target":row,"precision_bits":PRECISION,
                "producer_congruence":p,"second_association_readback":replay_result,**k_links,
                "external_independent_actual_review_required":True,
                "preconditioner_changed":False,"actual_spatial_H_or_force_jobs":0,
                "actual_AD_or_equilibrium_solves":0,"root_or_minimum_claimed":False}
        h=save(folder/"result.json",result)
        terminal.update(status="FINISHED",exit_code=0,result_sha256=h)
    except BaseException as error:
        save(folder/"failure.json",{"error":str(error),"traceback":traceback.format_exc(),**ident})
        raise
    finally:
        save(folder/"terminal.json",terminal)


def run(args):
    _,bindings=gate(args.review,args.review_sha256)
    ACTUAL.mkdir(exist_ok=False)
    pid=os.getpid()
    active={"status":"RUNNING","pid":pid,"bindings":bindings,"ordered_pilots":list(PILOTS)}
    save(ACTUAL/"started.json",active)
    final={"status":"FAILED","exit_code":1,"pid":pid,"bindings":bindings,"workers":[]}
    env={**os.environ,"OPENBLAS_NUM_THREADS":"1","MKL_NUM_THREADS":"1","OMP_NUM_THREADS":"1"}
    try:
        rows=[]
        for pilot in PILOTS:
            folder=ACTUAL/pilot
            folder.mkdir(exist_ok=False)
            cmd=[sys.executable,"-B",str(SOURCE),"worker","--pilot-id",pilot,"--supervisor-pid",str(pid),
                 "--review",str(args.review),"--review-sha256",args.review_sha256]
            with (folder/"stdout.log").open("x",encoding="utf8") as stdout:
                child=subprocess.Popen(cmd,cwd=ROOT,env=env,stdout=stdout,stderr=subprocess.STDOUT)
                save(folder/"supervisor_child_started.json",{"pilot_id":pilot,"pid":child.pid,"supervisor_pid":pid,"command":cmd})
                code=child.wait()
            t=folder/"terminal.json"
            observed={"pilot_id":pilot,"pid":child.pid,"supervisor_pid":pid,"observed_exit_code":code,
                      "terminal_sha256":sha(t) if t.exists() else None}
            oh=save(folder/"observed_os_exit.json",observed)
            final["workers"].append({**observed,"os_exit_sha256":oh})
            require(observed["terminal_sha256"] is not None,"missing child terminal")
            terminal=pinned(t,observed["terminal_sha256"])
            ident={"pilot_id":pilot,"pid":child.pid,"supervisor_pid":pid,"bindings":bindings}
            rfile=folder/"result.json"
            require(rfile.exists(),"missing child matrix outcome")
            rh=sha(rfile)
            check_child(observed,terminal,ident,rh)
            result=pinned(rfile,rh)
            require(all(result[k]==v for k,v in ident.items()) and result["root_or_minimum_claimed"] is False,
                    "mixed worker matrix result")
            rows.append({"pilot_id":pilot,"status":result["status"],"result_sha256":rh,
                         "terminal_sha256":observed["terminal_sha256"],"observed_os_exit_sha256":oh})
            print(json.dumps(rows[-1]),flush=True)
        rh=save(ACTUAL/"results.json",{"status":"COMPLETE_FOUR_PRECONDITIONER_SPD_ATTEMPTS_AWAITING_REVIEW",
                    "scope":SCOPE,"cases":rows,"bindings":bindings,"no_F_H_AD_solve_or_minimum":True})
        final.update(status="FINISHED",exit_code=0,results_sha256=rh)
    except BaseException as e:
        final["error"]=str(e)
        raise
    finally:
        save(ACTUAL/"terminal.json",final)


def rejects(f):
    try:
        f()
    except ValueError:
        return True
    return False


def fixtures():
    ident=np.eye(2,dtype="<f8")
    nonsym=np.array([[2.0,0.75],[-0.25,2.0]],dtype="<f8")
    proposal=ident/np.sqrt(2)
    p,K=producer.certify_sympart(nonsym,proposal,PRECISION)
    qi,rows=replay.congruence_bounds(nonsym,K,PRECISION)
    bad=np.diag(np.array([1.0,-1.0],dtype="<f8"))
    no,_=producer.certify_sympart(bad,ident,PRECISION)
    boundary=np.diag(np.array([1.0,0.0],dtype="<f8"))
    tie,_=producer.certify_sympart(ident,boundary,PRECISION)
    check={"pilot_id":PILOTS[0],"cache":"same","drive":"1/3","R_binary_sha256":find_matrix(nonsym,2)}
    match=lambda **kw: matched_target(check,{**check,**kw})
    assert match() is None
    who={"pilot_id":PILOTS[0],"pid":7,"supervisor_pid":4,"bindings":{}}
    term={**who,"status":"FINISHED","exit_code":0,"result_sha256":"a"*64}
    ob={"pilot_id":PILOTS[0],"pid":7,"supervisor_pid":4,"observed_exit_code":0}
    assert check_child(ob,term,who,"a"*64) is None
    return {
        "nonsymmetric_binary_R_exact_positive_symmetric_part_two_associations":
            p["status"]=="CERTIFIED_EXACT_SYMMETRIC_PART_SPD" and Q(p["q_infinity_upper_exact"])<1 and qi==max(rows)<1,
        "reject_exact_symmetric_indefinite_or_equality":Q(no["q_infinity_upper_exact"])>=1 and Q(tie["q_infinity_upper_exact"])==1,
        "typed_finite_binary_matrix_and_rank_guards":all(rejects(fn) for fn in (
            lambda:find_matrix(np.eye(2,dtype=int),2),lambda:find_matrix(np.eye(2,dtype=bool),2),
            lambda:find_matrix(np.array([[float("inf")]],dtype="<f8"),1),lambda:find_matrix(nonsym,3))),
        "reject_wrong_pilot_cache_drive_and_R_identity":all(rejects(fn) for fn in (
            lambda:match(pilot_id=PILOTS[1]),lambda:match(cache="foreign"),
            lambda:match(drive="0.3333333333333333"),lambda:match(R_binary_sha256=find_matrix(ident,2)),
            lambda:selected("foreign",[]))),
        "reject_nonliteral_or_mixed_worker_exit":(rejects(lambda:check_child({**ob,"observed_exit_code":False},term,who,"a"*64))
            and rejects(lambda:check_child(ob,{**term,"pid":8},who,"a"*64))),
    }


def prepare():
    PREP.mkdir(exist_ok=False)
    with (PREP/"transfer_ambiguity_preconditioner_spd_frozen.py").open("xb") as f:
        f.write(SOURCE.read_bytes())
    try:
        stop=AssertionError("cold preflight attempted a production matrix factory or mechanics job")
        with patch.object(EX,"open_context",side_effect=stop),patch.object(EX.T.core,"assemble_point",side_effect=stop), \
                patch.object(EX.T.core,"local_terms",side_effect=stop):
            snap=authority()
            checks=fixtures()
        checks["four_current_point_authorities_and_no_production_matrices"]=len(snap[1])==4 and not ACTUAL.exists()
        checks["outside_study_source_rejected"]=rejects(lambda:link({},"../PLAN.md","a"*64))
        require(yes(checks),"cold matrix preparation fixture failed")
        ph=save(PREP/"protocol.json",protocol(snap))
        pf=save(PREP/"preflight.json",{"status":"PASS","source_sha256":sha(SOURCE),"protocol_sha256":ph,
              "checks":checks,"actual_production_matrix_jobs":0,"new_F_H_AD_or_solve_jobs":0,
              "matrix_fixtures_scope":"Synthetic 2-by-2 only; no actual field/matrix factory.",
              "root_or_minimum_claimed":False})
        print(json.dumps({"source_sha256":sha(SOURCE),"protocol_sha256":ph,"preflight_sha256":pf}))
    except BaseException as e:
        save(PREP/"failure.json",{"error":str(e),"traceback":traceback.format_exc()})
        raise


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    sub=cli.add_subparsers(dest="stage",required=True)
    sub.add_parser("prepare")
    for mode in ("run","worker"):
        p=sub.add_parser(mode)
        p.add_argument("--review",type=lambda s:Path(s).resolve(),required=True)
        p.add_argument("--review-sha256",required=True)
        if mode=="worker":
            p.add_argument("--pilot-id",choices=PILOTS,required=True)
            p.add_argument("--supervisor-pid",type=int,required=True)
    args=cli.parse_args()
    if args.stage=="prepare":
        prepare()
    elif args.stage=="run":
        run(args)
    else:
        worker(args)
