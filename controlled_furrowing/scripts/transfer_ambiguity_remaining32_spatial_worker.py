"""Cold first-wave raw exact-spatial worker with independent selector admission."""
import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from contextlib import contextmanager
from itertools import combinations_with_replacement
import os
from pathlib import Path
import signal
from threading import Event
import sys
import traceback
from types import SimpleNamespace
from unittest.mock import patch

import transfer_ambiguity_remaining32_exact_bridge as S

SOURCE=Path(__file__).resolve()
ROOT=S.ROOT
V=S.V
B=S.B
PREP=S.OUT/"worker_preparation"
ACTUAL=S.OUT/"actual_first_wave"
INDEPENDENT=PREP/"independent_prelaunch/review.json"
WAVE_GATE=S.OUT/"admission/first_wave_worker_admission.json"
BRIDGE_SOURCE="4d52e8fb8dbc01ba58c38052bbff537443693234499a9a3b9b3e138c82d9b51b"
BRIDGE_PROTOCOL="ade6158b96bec63665880beaae1447057b36f162863211acf1f8d24f95519664"
BRIDGE_MANIFEST="d42b7260283ba02ccc34b9e9ba69ce76d58507115c22615c90d28fa38020cb41"
BRIDGE_PREFLIGHT="a6741fb642a4006b126ac71db4fdb964761f5245e2563cd70a0175db6ee18bda"
BRIDGE_REVIEW="89f9b5cacf89b4055fc404c00d3e3c42cf10fed20337fe455e996c2944ef73b2"
LAUNCHER_SOURCE="68bca869ee486d7dc25b34cae3ceee1ddc2adc53500a20301c0038c259a2b26a"
PINS={S.SOURCE:BRIDGE_SOURCE,S.PREP/"protocol.json":BRIDGE_PROTOCOL,S.PREP/"manifest.json":BRIDGE_MANIFEST,
      S.PREP/"preflight.json":BRIDGE_PREFLIGHT,S.COLD_REVIEW_PATH:BRIDGE_REVIEW,
      Path(V.__file__).resolve():LAUNCHER_SOURCE}
FIRST=("max_contrast_change_small_01","max_contrast_change_small_13",
       "max_negative_contrast_change_small_01","max_negative_contrast_change_small_16")
COUNTS={"bulk":144,"surface":24}
MAX_WORKERS=4
JOBS=tuple((case,part) for case in FIRST for part in COUNTS)
STAGE="remaining32_spatial_first_wave"
HOLD=S.OUT/"worker_development_revision_0/independent_prelaunch/hold_300923.json"
HOLD_SHA="1b2b81878bfa1af5663461ac4a51be0f20c02976151a21357379d605b2658993"
HANDOFF_HOLD=S.OUT/"worker_development_revision_1/independent_prelaunch/hold_301fc5.json"
HANDOFF_HOLD_SHA="c8c9237b02a1dba3ed2226e2cbb7eec385ab5fcddc691a13ee084ffa4eeddc5e"
REVIEW="PASS_INDEPENDENT_REMAINING32_FIRST_WAVE_SPATIAL_WORKER_COLD"
ADMISSION="PASS_INDEPENDENT_REMAINING32_FIRST_WAVE_SPATIAL_ADMISSION"
STATUS="COMPLETE_RAW_REMAINING32_SPATIAL_PART_AWAITS_INDEPENDENT_REPLAY"
SCOPE="First four explicit targets only. Raw separate exact center F and whole admitted r34 H defects with same-deformation auxiliary/cache correction."
ZERO={"new_spatial_integrations":0,"new_derivative_jobs":0,"new_solves":0}
read,save,sha,require=S.read,S.save,S.sha,S.require
canonical,rel=S.canonical,S.rel


def bridge_args(case,admission):
    return SimpleNamespace(case_id=case,protocol_sha256=BRIDGE_PROTOCOL,manifest_sha256=BRIDGE_MANIFEST,
                           preflight_sha256=BRIDGE_PREFLIGHT,review=S.COLD_REVIEW_PATH,review_sha256=BRIDGE_REVIEW,
                           admission=S.JOINED,admission_sha256=admission)


def bridge_manifest():
    require(all(sha(p)==h for p,h in PINS.items()),"unchanged reviewed bridge and supervisor lifecycle sources")
    m,b=S.gate(bridge_args(FIRST[0],"unused-in-cold-stage"))
    require(len(m["cases"])==32 and tuple(x["case_id"] for x in m["cases"])==V.A.IDS
            and m["components_catalogued"]==5376 and m["no_spatial_workers_authorized"] is True,
            "full frozen catalogue, separate first-wave authorization")
    cases={x["case_id"]:x for x in m["cases"]}
    require(tuple(case for case in V.A.IDS if case in FIRST)==FIRST
            and all(cases[a]["cached_input_kind"]=="FRESH_R34" for a in FIRST)
            and tuple(x for x in FIRST if S.Q(cases[x]["trial"]["archived_center_shift_exact"])>0)
                ==(FIRST[1],FIRST[3]),"exact declared ordinary and recentered wave cases")
    return cases,b


def case_payload(entry):
    t=entry["trial"]
    keys=("trial_box_sha256","model_contract_sha256","point_cache_sha256","R_binary_sha256",
          "center_binary_sha256","archived_center_shift_exact","exact_drive","energy_scale_binary_exact")
    return {"case_id":entry["case_id"],"record_id":entry["record_id"],
            "original_observer_field_sha256":entry["original_observer_field_sha256"],
            "original_addition_gap_exact":entry["original_addition_gap_exact"],
            "observer_goal_exact":entry["observer_goal_exact"],**{k:t[k] for k in keys}}


def plan_docs(a,case):
    row=V.A.pick(a["plans"],case)
    target=a["targets"][case]
    require(row["target"]==target and row["case_id"]==case,"same original observer and exact proof identity")
    folder=S.PLAN/case
    outer=read(folder/"root_os_exit_observation.json",row["observed_os_exit_sha256"])
    require(row["plan_sha256"]==outer["plan_sha256"] and row["terminal_sha256"]==outer["terminal_sha256"],
            "exact separate root observed selector and terminal hashes")
    review_file=canonical(row["independent_plan_review_file"])
    require(review_file.is_relative_to(S.OUT) and review_file.name.endswith(".json"),"local actual plan review required")
    review=read(review_file,row["independent_plan_review_sha256"])
    p=read(folder/"plan.json",row["plan_sha256"])
    t=read(folder/"terminal.json",row["terminal_sha256"])
    require(review["observed_os_exit_sha256"]==row["observed_os_exit_sha256"]
            and p["bindings"]["bridge_input_admission_sha256"]==a["joined_bridge_input_admission_sha256"]
            and p["trial_box_sha256"]==target["trial_box_sha256"]
            and p["record_id"]==t["record_id"]==target["record_id"],
            "same joined admission, independent reviewed plan and trial")
    S.validate_plan(outer,t,p,review,row["plan_sha256"],row["terminal_sha256"],p["bindings"],case,S.OPS)
    return row,review_file


def worker_gate_data(review,admission,expected,targets):
    require(review["status"]==REVIEW and V.yes(review["checks"])
            and all(review[k]==v for k,v in expected.items())
            and review["permitted_first_wave_case_ids"]==list(FIRST)
            and review["maximum_owned_workers"]==MAX_WORKERS
            and review["no_later_waves_authorized"] is True,
            "external independent cold source review of first wave only")
    require(admission["status"]==ADMISSION and V.yes(admission["checks"])
            and all(admission[k]==v for k,v in expected.items())
            and tuple(r["case_id"] for r in admission["plans"])==FIRST
            and admission["targets"]==targets
            and admission["no_later_waves_authorized"] is True
            and admission["independent_packet_replay_required"] is True, "external first-wave selected plans admission")


def gate(args):
    cases,b=bridge_manifest()
    p=read(PREP/"protocol.json",args.protocol_sha256)
    f=read(PREP/"preflight.json",args.preflight_sha256)
    require(args.review.resolve()==INDEPENDENT.resolve() and args.admission.resolve()==WAVE_GATE.resolve(),
            "canonical separate first-wave external admissions")
    r=read(args.review,args.review_sha256)
    a=read(args.admission,args.admission_sha256)
    wanted={"source_sha256":sha(SOURCE),"protocol_sha256":args.protocol_sha256,"preflight_sha256":args.preflight_sha256,
            "bridge_manifest_sha256":BRIDGE_MANIFEST,"bridge_cold_review_sha256":BRIDGE_REVIEW}
    require(p["status"]=="FROZEN_RAW_EXACT_SPATIAL_FIRST_WAVE_WORKER"
            and p["source_sha256"]==f["source_sha256"]==wanted["source_sha256"]
            and f["protocol_sha256"]==args.protocol_sha256 and f["status"]=="PASS" and V.yes(f["checks"])
            and all(V.zero(f[k]) for k in ZERO) and p["catalogue_binding"]==b
            and p["lifecycle_source_sha256"]==LAUNCHER_SOURCE and sha(PREP/"worker_frozen.py")==wanted["source_sha256"]
            and p["first_wave"]==list(FIRST) and p["max_owned_workers"]==MAX_WORKERS,
            "current frozen first-wave worker preparation")
    targets={case:case_payload(cases[case]) for case in FIRST}
    worker_gate_data(r,a,wanted,targets)
    require(a["worker_source_review_sha256"]==args.review_sha256,"same external reviewed worker source")
    joined=read(S.JOINED,a["joined_bridge_input_admission_sha256"])
    require(joined["status"]==S.ADMISSION and V.yes(joined["checks"]) and joined["spatial_workers_authorized"] is False
            and all(joined[k]==v for k,v in b.items()) and joined["same_full_32_geometry_10_fresh_derivative_required"] is True,
            "same joined proof, with independent separate worker approval")
    return a,{**wanted,"worker_source_review_sha256":args.review_sha256,"worker_wave_admission_sha256":args.admission_sha256,
              "joined_bridge_input_admission_sha256":a["joined_bridge_input_admission_sha256"]}


def reviewed_context(a,case):
    row,path=plan_docs(a,case)
    d,tails,plans=S.open_reviewed_plan(bridge_args(case,a["joined_bridge_input_admission_sha256"]),
                          row["observed_os_exit_sha256"],path,row["independent_plan_review_sha256"])
    require(case_payload(d["entry"])==row["target"] and tuple(plans)==S.OPS and len(tails)==2 and min(tails)>0,
            "same explicit trial and complete independent selected components")
    return row,d,tails,plans


def path_for(case,part):
    return ACTUAL/"cases"/case/part


def rational_interval(x):
    require(type(x) is list and len(x)==2,"two exact raw endpoints")
    a,b=map(S.Q,x)
    require(a<=b,"ordered signed raw endpoints")


def validate_packet(packet,proof,case,slot,part,tail):
    n=len(proof["global_indices"])
    require(type(packet["slot"]) is int and packet["slot"]==slot and packet["kind"]==part
            and packet["record_id"]==case["record_id"] and packet["precision_bits"]==640
            and packet["global_indices"]==proof["global_indices"]
            and packet["proof_sha256"]==B.digest(proof) and packet["root_asserted"] is False,
            "same typed reviewed raw proof and projected component only")
    f=packet["point_force"]
    h=packet["whole_box_hessian"]
    F=f["F_I_minus_F_cache_intervals_exact"]
    correction=f["Qstar_minus_cache_correction_exact"]
    t=f["force_tail_upper_exact"]
    U=h["H_I_minus_H_cache_upper_exact"]
    require(f["scope"]==S.POINT and h["scope"]==S.BOX and S.Q(h["radius_exact"])==S.RADIUS
            and len(F)==len(correction)==len(t)==n and all(0<=S.Q(x)<=tail[0] for x in t),
            "point force/whole-H roles and exact tail and cache corrections")
    for v in (*F,*correction):
        rational_interval(v)
    require(all(type(x) is list and len(x)==4 and type(x[0]) is int and type(x[1]) is int for x in U)
            and [(i,j) for i,j,*_ in U]==list(combinations_with_replacement(range(n),2)),
            "full ordered typed local symmetric upper Hessian")
    for _,_,lo,hi in U:
        rational_interval([lo,hi])
    return n,len(U)


def worker(args):
    a,links=gate(args)
    start=read(ACTUAL/"started.json",args.supervisor_started_sha256)
    job=(args.case_id,args.part)
    require(job in JOBS and V.positive_int(start["pid"]) and V.positive_int(args.supervisor_pid)
            and start["pid"]==args.supervisor_pid==os.getppid() and start["bindings"]==links
            and start["stage"]==STAGE and start["jobs"]==[list(j) for j in JOBS]
            and start["status"]=="RUNNING","only authorized live current first-wave supervisor")
    folder=path_for(*job)
    mark={"pid":os.getpid(),"supervisor_pid":args.supervisor_pid,"stage":STAGE,"case_id":job[0],"part":job[1],
          "bindings":links}
    save(folder/"started.json",{"status":"RUNNING",**mark})
    terminal={"status":"FAILED","exit_code":1,**mark}
    done=[]
    try:
        row,d,tail,proofs=reviewed_context(a,job[0])
        force,scaled=S.kernels()
        nforce=0
        nH=0
        for slot in range(COUNTS[job[1]]):
            proof=proofs[(job[1],slot)]
            entry=B.evaluate(d,proof,tail,force,scaled)
            nf,nupper=validate_packet(entry,proof,row["target"],slot,job[1],tail)
            nforce+=nf
            nH+=nupper
            name=f"{job[1]}_{slot:03d}.json"
            h=save(folder/name,V.A.doc({**entry,"worker_identity":mark,
                              "trial_box_sha256":row["target"]["trial_box_sha256"],"selector_plan_sha256":row["plan_sha256"]}))
            done.append({"kind":job[1],"slot":slot,"file":name,"sha256":h,"force_entries":nf,"upper_H_entries":nupper})
        manifest={"status":STATUS,**mark,"target":row["target"],"plan":row,"components":done,
                  "counts":{"components":len(done),"force_entries":nforce,"upper_H_entries":nH},
                  "point_force_scope":S.POINT,"whole_H_scope":S.BOX,"radius_exact":str(S.RADIUS),
                  "root_asserted":False,"packet_independent_replay_required":True}
        h=save(folder/"manifest.json",manifest)
        terminal.update(status="FINISHED",exit_code=0,results_sha256=h)
    except BaseException as error:
        save(folder/"failure.json",{**mark,"error":str(error),"traceback":traceback.format_exc(),"completed_packets":done})
        raise
    finally:
        save(folder/"terminal.json",terminal)

def observed_part(a,job,pid,links,started_hash):
    case,part=job
    folder=path_for(case,part)
    folder.mkdir(parents=True,exist_ok=False)
    mark={"supervisor_pid":pid,"stage":STAGE,"case_id":case,"part":part,"bindings":links}
    cmd=[sys.executable,"-B",str(SOURCE),"worker","--case-id",case,"--part",part,
         "--supervisor-pid",str(pid),"--supervisor-started-sha256",started_hash,
         "--protocol-sha256",links["protocol_sha256"],"--preflight-sha256",links["preflight_sha256"],
         "--review",str(INDEPENDENT),"--review-sha256",links["worker_source_review_sha256"],
         "--admission",str(WAVE_GATE),"--admission-sha256",links["worker_wave_admission_sha256"]]
    owned=[]
    identity,code=V.run_owned_child(cmd,folder,mark,owned)
    t=folder/"terminal.json"
    native_outcome={**identity,"observed_exit_code":code,"terminal_sha256":sha(t) if t.is_file() else None}
    oh=save(folder/"observed_os_exit.json",native_outcome)
    require(native_outcome["terminal_sha256"] is not None and V.zero(code),"literal child process failure retained")
    term=read(t,native_outcome["terminal_sha256"])
    h=term["results_sha256"]
    m=read(folder/"manifest.json",h)
    V.child_admit(code,m,term,identity,h)
    row=V.A.pick(a["plans"],case)
    require(m["status"]==STATUS and m["target"]==row["target"] and m["plan"]==row
            and m["radius_exact"]==str(S.RADIUS) and m["packet_independent_replay_required"] is True
            and (m["point_force_scope"],m["whole_H_scope"])==(S.POINT,S.BOX),
            "same canonical part and independently reviewed plan")
    require(m["counts"]["components"]==len(m["components"])==COUNTS[part] and
            all(type(x["slot"]) is int and x["slot"]==i and x["kind"]==part
                and x["file"]==f"{part}_{i:03d}.json"
                and V.positive_int(x["force_entries"]) and V.positive_int(x["upper_H_entries"])
                and x["upper_H_entries"]==x["force_entries"]*(x["force_entries"]+1)//2
                for i,x in enumerate(m["components"])), "complete ordered scoped raw part packet metadata")
    require(all(sha(folder/x["file"])==x["sha256"] for x in m["components"])
            and all(sum(x[k] for x in m["components"])==m["counts"][k] for k in ("force_entries","upper_H_entries")),
            "every raw part packet digest and entry count")
    return {**native_outcome,"manifest_sha256":h,"observed_os_exit_sha256":oh,
            "manifest_file":rel(folder/"manifest.json"),"counts":m["counts"]}


@contextmanager
def deferred_interrupts():
    """Record console interrupts without asynchronous exceptions in the owned-work region."""
    signs=[signal.SIGINT]
    if hasattr(signal,"SIGBREAK"):
        signs.append(signal.SIGBREAK)
    previous={}
    requested=[False]

    def on_interrupt(signum,frame):
        requested[0]=True

    try:
        for sig in signs:
            previous[sig]=signal.getsignal(sig)
            signal.signal(sig,on_interrupt)
        yield lambda:requested[0]
    finally:
        for sig,handler in previous.items():
            signal.signal(sig,handler)


def drain_submitted(futures):
    """Only Future completion proves all owned child functions really finished."""
    while True:
        try:
            pending=[f for f in futures if not f.done()]
            if not pending:
                return
            wait(pending,timeout=.2,return_when=FIRST_COMPLETED)
        except KeyboardInterrupt:
            continue


def bounded(items,run_job,limit,interrupted=lambda:False):
    require(type(limit) is int and 1<=limit<=MAX_WORKERS,"bounded owned worker limit")
    waiting=iter(items)
    complete={}
    submitted=[]
    pending={}
    pool=ThreadPoolExecutor(max_workers=limit)

    def check_stop():
        if interrupted():
            raise KeyboardInterrupt("console interrupt deferred until every owned worker exits")

    def add():
        check_stop()
        for item in waiting:
            ready=Event()
            armed=False

            def owned_task():
                ready.wait()
                if armed:
                    return run_job(item)
                return None

            try:
                # CPython can queue the callable before submit() raises while
                # creating a thread. Such an unreturned item stays disarmed.
                f=pool.submit(owned_task)
                submitted.append(f)
                pending[f]=item
                check_stop()
                armed=True
            finally:
                ready.set()
            break

    try:
        for _ in range(min(limit,len(items))):
            add()
        while pending:
            check_stop()
            batch,_=wait(tuple(pending),timeout=.2,return_when=FIRST_COMPLETED)
            check_stop()
            for future in batch:
                item=pending.pop(future)
                complete[item]=future.result()
            for _ in batch:
                add()
        check_stop()
    finally:
        # A second or later raised interrupt cannot skip this observed drain.
        # The console handlers remain installed around the entire owner region.
        drain_submitted(submitted)
        # No native child can be live now. Do not rely on Thread.join's
        # interrupt-sensitive bookkeeping to establish completion.
        pool.shutdown(wait=False,cancel_futures=True)
    return [complete[item] for item in items]


def run(args):
    a,links=gate(args)
    for case in FIRST:
        reviewed_context(a,case)
    with V.stage_lock(ACTUAL),deferred_interrupts() as interrupted:
        require(not (ACTUAL/"started.json").exists() and not (ACTUAL/"cases").exists(),
                "exclusive raw first wave. Preserve prior aborted or completed reservations.")
        pid=os.getpid()
        ident={"pid":pid,"stage":STAGE,"jobs":[list(j) for j in JOBS],"bindings":links}
        h=save(ACTUAL/"started.json",{"status":"RUNNING",**ident})
        end={"status":"FAILED","exit_code":1,**ident}
        try:
            parts=bounded(JOBS,lambda job:observed_part(a,job,pid,links,h),MAX_WORKERS,interrupted)
            if interrupted():
                raise KeyboardInterrupt("console interrupt observed after owned workers were reaped")
            require([(p["case_id"],p["part"]) for p in parts]==list(JOBS)
                    and len({r["pid"] for r in parts})==len(JOBS),"all eight distinct true OS parts required")
            target_rows=[]
            for i in range(0,len(parts),2):
                bulk,surface=parts[i:i+2]
                case=bulk["case_id"]
                require(case==surface["case_id"] and (bulk["part"],surface["part"])==("bulk","surface")
                        and bulk["counts"]["components"]+surface["counts"]["components"]==168,
                        "one complete pair of actual parts per same record")
                target_rows.append({"case_id":case,"target":a["targets"][case],"parts":[bulk,surface]})
            complete={"status":"COMPLETE_FOUR_RAW_SPATIAL_BRIDGES_AWAITING_INDEPENDENT_PACKET_REPLAY",
                      **ident,"cases":target_rows,"raw_components":672,"later_waves_launched":False,
                      "independent_packet_replay_required":True,"root_asserted":False}
            mh=save(ACTUAL/"manifest.json",complete)
            end.update(status="FINISHED",exit_code=0,manifest_sha256=mh)
            print({"manifest_sha256":mh,"native_supervisor_pid":pid,"completed_parts":len(parts)},flush=True)
        except BaseException as error:
            files=(path_for(c,p)/f for c,p in JOBS
                   for f in ("observed_os_exit.json","abnormal_child_reap.json","terminal.json","failure.json"))
            # bounded() has drained all futures before this failure receipt.
            save(ACTUAL/"failure.json",{**ident,"error":str(error),"traceback":traceback.format_exc(),
                 "closed_worker_attempt_sha256":{rel(p):sha(p) for p in files if p.is_file()}})
            raise
        finally:
            save(ACTUAL/"terminal.json",end)


def fixture_gate():
    wanted={"source_sha256":"a"*64,"bridge_manifest_sha256":BRIDGE_MANIFEST}
    targets={case:{"case_id":case} for case in FIRST}
    review={"status":REVIEW,**wanted,"checks":{"synthetic":True},
            "permitted_first_wave_case_ids":list(FIRST),"maximum_owned_workers":MAX_WORKERS,
            "no_later_waves_authorized":True}
    admission={"status":ADMISSION,**wanted,"checks":{"synthetic":True},
               "plans":[{"case_id":case} for case in FIRST],"targets":targets,
               "no_later_waves_authorized":True,"independent_packet_replay_required":True}
    worker_gate_data(review,admission,wanted,targets)
    return all(V.A.reject(lambda:worker_gate_data(review,x,wanted,targets)) for x in (
        {**admission,"checks":{"synthetic":1}},
        {**admission,"targets":{**targets,FIRST[0]:{"case_id":FIRST[1]}}},
        {**admission,"plans":[*admission["plans"],{"case_id":"later"}]},
    ))


def fixture_raw():
    proof={"kind":"bulk","slot":0,"global_indices":[1,0]}
    case={"record_id":"synthetic"}
    packet={"kind":"bulk","slot":0,"record_id":"synthetic","precision_bits":640,
        "global_indices":[1,0],"proof_sha256":B.digest(proof),"root_asserted":False,
        "point_force":{"scope":S.POINT,"F_I_minus_F_cache_intervals_exact":[["-2","1"],["0","0"]],
             "Qstar_minus_cache_correction_exact":[["0","0"],["1/8","1/4"]],"force_tail_upper_exact":["1/8","0"]},
        "whole_box_hessian":{"scope":S.BOX,"radius_exact":str(S.RADIUS),
              "H_I_minus_H_cache_upper_exact":[[0,0,"1","2"],[0,1,"-1","0"],[1,1,"0","1"]]}}
    require(validate_packet(packet,proof,case,0,"bulk",(S.Q("1/8"),S.Q(1)))==(2,3),"synthetic valid raw map")
    for bad in (
        {**packet,"record_id":"sibling"},
        {**packet,"slot":False},
        {**packet,"point_force":{**packet["point_force"],"scope":S.BOX}},
        {**packet,"point_force":{**packet["point_force"],"Qstar_minus_cache_correction_exact":[["0","0"]]}},
        {**packet,"whole_box_hessian":{**packet["whole_box_hessian"],"radius_exact":"0"}},
        {**packet,"whole_box_hessian":{**packet["whole_box_hessian"],"H_I_minus_H_cache_upper_exact":[[0,0,"1","2"]]}},
    ):
        require(V.A.reject(lambda:validate_packet(bad,proof,case,0,"bulk",(S.Q("1/8"),S.Q(1)))),
                "reject mixed scope, wrong case, incomplete finite correction or H, bool slot")
    return True


def fixture_parallel():
    from threading import Barrier, Lock
    barrier=Barrier(MAX_WORKERS)
    lock=Lock()
    active=peak=0
    def task(i):
        nonlocal active,peak
        with lock:
            active+=1
            peak=max(peak,active)
        if i<MAX_WORKERS:
            barrier.wait(timeout=10)
        with lock:
            active-=1
        return i*i
    items=list(range(8))
    values=bounded(items,task,MAX_WORKERS)
    return peak==MAX_WORKERS and active==0 and values==[x*x for x in items] and V.A.reject(
        lambda:bounded(items,task,MAX_WORKERS+1))


def fixture_repeated_interrupts():
    started=Event()
    release=Event()
    finished=Event()
    events=[]
    real_wait=wait

    def owned(name):
        require(name=="first","queued task scheduled after interrupted wait")
        started.set()
        require(release.wait(10),"cleanup fixture release timeout")
        finished.set()

    def corrupt_wait(*a,**kw):
        require(started.wait(10),"future did not actually enter its task")
        events.append("wait")
        if len(events)<=3:
            raise KeyboardInterrupt("repeat interrupt while draining")
        release.set()
        return real_wait(*a,**kw)

    with patch.object(sys.modules[__name__],"wait",side_effect=corrupt_wait):
        try:
            bounded(("first","queued"),owned,1)
        except KeyboardInterrupt as error:
            require(str(error)=="repeat interrupt while draining","original interrupt lost")
        else:
            raise AssertionError("interrupted main future wait should be retained")
    require(finished.is_set() and len(events)>=4,"pending future left live after repeated interrupts")
    prior=signal.getsignal(signal.SIGINT)
    with deferred_interrupts() as requested:
        require(requested() is False,"unexpected initial interrupt")
        signal.raise_signal(signal.SIGINT)
        signal.raise_signal(signal.SIGINT)
        require(requested() is True,"real repeated console interrupt was lost")
    return signal.getsignal(signal.SIGINT)==prior


def fixture_submit_handoff():
    consumed=Event()
    first_finished=Event()
    unblock_first=Event()
    blocked_returned=Event()
    attempted=[]
    executor=ThreadPoolExecutor

    class FailAfterQueue(executor):
        calls=0

        def submit(self,fn,*args,**kwargs):
            self.calls+=1
            if self.calls==2:
                # The existing worker genuinely reaches the second callable
                # while the main submit cannot have returned its Future.
                def queued():
                    consumed.set()
                    try:
                        return fn(*args,**kwargs)
                    finally:
                        blocked_returned.set()
                return super().submit(queued)
            return super().submit(fn,*args,**kwargs)

        def _adjust_thread_count(self):
            if self.calls==2:
                unblock_first.set()
                require(consumed.wait(10),"native queue fixture was not consumed")
                raise RuntimeError("injected thread creation failure after queue")
            return super()._adjust_thread_count()

    def actual_job(name):
        attempted.append(name)
        if name=="tracked":
            require(unblock_first.wait(10),"first owned task did not drain")
            first_finished.set()
        else:
            V.subprocess.Popen(["must-not-start"])

    with patch.object(sys.modules[__name__],"ThreadPoolExecutor",FailAfterQueue), \
         patch.object(V.subprocess,"Popen",side_effect=AssertionError("unreturned task started a child")) as child:
        try:
            bounded(("tracked","unreturned"),actual_job,2)
        except RuntimeError as error:
            require(str(error)=="injected thread creation failure after queue","original submit failure lost")
        else:
            raise AssertionError("expected live queue-before-return error")
        require(blocked_returned.wait(10) and first_finished.is_set() and not child.called,
                "queued unreturned wrapper did not quiesce harmlessly")
    return consumed.is_set() and attempted==["tracked"]


def prepare():
    cases,binding=bridge_manifest()
    require(sha(HOLD)==HOLD_SHA and sha(HANDOFF_HOLD)==HANDOFF_HOLD_SHA,"retain independent predecessor HOLD receipts")
    checks={"record_target_scope_first_wave_only":fixture_gate(),
            "repeated_future_wait_interrupts_drained_before_shutdown":fixture_repeated_interrupts(),
            "queued_unreturned_submit_cannot_start_unowned_child":fixture_submit_handoff(),
            "literal_maps_scopes_complete_H_and_same_Y_correction":fixture_raw(),
            "bounded_scheduler_real_thread_order":fixture_parallel(),
            "cold_no_live_selector_or_spatial_payload_opened":True,
            "previous_reviewed_owned_child_exit_lifecycle_pinned":sha(Path(V.__file__).resolve())==LAUNCHER_SOURCE}
    require(V.yes(checks) and not ACTUAL.exists(),"cold checked first wave, no worker stage")
    PREP.mkdir(parents=True,exist_ok=False)
    p={"status":"FROZEN_RAW_EXACT_SPATIAL_FIRST_WAVE_WORKER","scope":SCOPE,"source_sha256":sha(SOURCE),
       **ZERO,"first_wave":list(FIRST),"jobs":[list(x) for x in JOBS],"max_owned_workers":MAX_WORKERS,
       "targets":{case:case_payload(cases[case]) for case in FIRST},
       "catalogue_binding":binding,"lifecycle_source_sha256":LAUNCHER_SOURCE,
       "source_and_protocol_sha256":{rel(p):h for p,h in PINS.items()},
       "worker_start_requires_external_same_case_plan_and_joined_proof":True,
       "native_child_and_external_outer_exits_are_separate":True,"later_waves_authorized":False,
       "superseded_cold_holds":[{"file":rel(path),"sha256":digest} for path,digest in
                               ((HOLD,HOLD_SHA),(HANDOFF_HOLD,HANDOFF_HOLD_SHA))],
       "root_asserted":False}
    ph=save(PREP/"protocol.json",p)
    pf=save(PREP/"preflight.json",{"status":"PASS","source_sha256":sha(SOURCE),"protocol_sha256":ph,
                                **ZERO,"checks":checks})
    with (PREP/"worker_frozen.py").open("xb") as f:
        f.write(SOURCE.read_bytes())
    print({"source_sha256":sha(SOURCE),"protocol_sha256":ph,"preflight_sha256":pf},flush=True)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    sub=cli.add_subparsers(dest="phase",required=True)
    sub.add_parser("prepare")
    for name in ("run","worker"):
        target=sub.add_parser(name)
        for key in ("protocol","preflight","review","admission"):
            target.add_argument("--"+key+"-sha256",required=True)
            if key in ("review","admission"):
                target.add_argument("--"+key,type=lambda x:Path(x).resolve(),required=True)
        if name=="worker":
            target.add_argument("--case-id",required=True,choices=FIRST)
            target.add_argument("--part",required=True,choices=tuple(COUNTS))
            target.add_argument("--supervisor-pid",required=True,type=int)
            target.add_argument("--supervisor-started-sha256",required=True)
    args=cli.parse_args()
    if args.phase=="prepare":
        prepare()
    elif args.phase=="run":
        run(args)
    else:
        worker(args)
