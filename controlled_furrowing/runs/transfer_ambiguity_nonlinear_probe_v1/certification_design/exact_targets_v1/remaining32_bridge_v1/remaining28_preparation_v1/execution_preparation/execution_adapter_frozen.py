"""Cold wave-aware selector command and raw spatial worker adapter for remaining28."""
import argparse
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import os
import sys
import traceback

import profiles as P

SOURCE=Path(__file__).resolve()
HERE=SOURCE.parent
ROOT=P.ROOT
sys.path.insert(0,str(ROOT/"scripts"))
import transfer_ambiguity_remaining32_spatial_worker as W

S,V=W.S,W.V
PREP=HERE/"execution_preparation"
RUNS=HERE/"waves"
LOCK=HERE/"scientific_process_budget"
REVIEW_PATH=PREP/"independent_prelaunch/review.json"
RELEASE_PATH=HERE/"admission/later_waves_release.json"
IND_FIRST=S.OUT/"first_wave_reduction_v1/independent_actual_preparation"
NATIVE_FIRST=IND_FIRST/"actual_replay"
JOINED_SHA="2e2e388e798ed509bb52bc72a6169b8b6ee5db55173b3529ac73c3812b23a6cb"
PROFILE_SHA="31c680040387762f835a235e678435400ac163081509abf9b012494d20de04f0"
PROFILE_CODE_SHA="56ec56ea88d5880bea9f93ebf448165f391635d157345a64d119fe2ec9d00a42"
SOURCE_REVIEW="PASS_INDEPENDENT_REMAINING28_WAVE_EXECUTION_COLD"
RELEASE="PASS_ROOT_RELEASE_REMAINING28_AFTER_FIRST_WAVE_GLOBAL_REPLAY"
FIRST_OS="PASS_ROOT_OBSERVED_FIRST_WAVE_INDEPENDENT_GLOBAL_EXIT"
WAVE_REVIEW="PASS_INDEPENDENT_LATER_WAVE_SPATIAL_ADMISSION"
RAW_STATUS="COMPLETE_RAW_LATER_WAVE_SPATIAL_PART_AWAITS_INDEPENDENT_REPLAY"
ZERO={"actual_selector_jobs":0,"actual_raw_spatial_jobs":0,"actual_hull_replays":0,"actual_global_reductions":0}
COUNTS=W.COUNTS
WAVES=tuple(f"wave_{i:02}" for i in range(2,9))
read,save,sha,require=S.read,S.save,S.sha,S.require
rel,canonical=S.rel,S.canonical


def cold():
    require(sha(HERE/"profiles.py")==PROFILE_CODE_SHA and sha(HERE/"profiles.json")==PROFILE_SHA,
            "same frozen inherited profiles")
    for p,h in P.PINS.items():
        require(sha(p)==h,"unchanged reviewed scientific or lifecycle source")
    original,b=W.bridge_manifest()
    result={}
    for wave in WAVES:
        profile=P.load_profile(HERE/"profiles.json",PROFILE_SHA,wave)
        ids=profile["case_order"]
        identity={case:W.case_payload(original[case]) for case in ids}
        for name,row in identity.items():
            saved=profile["targets"][name]
            require(all(saved[k]==v for k,v in row.items()),"same inherited original archive and proof authority")
        result[wave]={**profile,"raw_worker_targets":identity}
    return result,b


def source_gate(review,plan,pre,expected):
    require(review["status"]==SOURCE_REVIEW and V.yes(review["checks"]) and all(review[k]==v for k,v in expected.items())
            and review["permitted_waves"]==list(WAVES) and review["maximum_owned_scientific_processes"]==4
            and review["actual_release_required"] is True
            and plan["status"]=="FROZEN_REMAINING28_EXECUTION_ADAPTER"
            and plan["source_sha256"]==pre["source_sha256"]==expected["source_sha256"]
            and plan["profiles_sha256"]==expected["profiles_sha256"]==PROFILE_SHA
            and pre["protocol_sha256"]==expected["protocol_sha256"] and pre["status"]=="PASS"
            and all(V.zero(pre[k]) for k in ZERO) and V.yes(pre["checks"]),
            "same independently reviewed execution adapter and cap")


def release_summary(x,observed,end,native,bind):
    require(x["status"]==RELEASE and V.yes(x["checks"]) and all(x[k]==v for k,v in bind.items())
            and x["case_order"]==list(W.FIRST) and x["following_waves"]==list(WAVES)
            and x["selector_invocations_only_under_this_release"] is True
            and x["per_wave_plan_review_and_raw_admission_required"] is True,
            "root release must bind complete first wave and all later profiles")
    require(observed["status"]==FIRST_OS and V.positive_int(observed["pid"]) and V.zero(observed["observed_exit_code"])
            and V.positive_int(observed["unified_exec_session"])
            and all(V.positive_int(z["pid"]) and z["pid"]==observed["pid"] for z in (end,native))
            and end["status"]=="FINISHED" and V.zero(end["exit_code"])
            and observed["receipt_sha256"]==end["receipt_sha256"]==x["first_wave_independent_receipt_sha256"]
            and native["status"]=="COMPLETE_INDEPENDENT_FIRST_WAVE_SIGNED_ENDPOINT_REPLAY"
            and tuple(x["case_id"] for x in native["case_results"])==W.FIRST
            and native["full_five_gate_complete"] is False,"literal entire first-wave independent actual closure")


def release_gate(args):
    allprofiles,base=cold()
    p=read(PREP/"protocol.json",args.protocol_sha256)
    f=read(PREP/"preflight.json",args.preflight_sha256)
    require(args.review.resolve()==REVIEW_PATH.resolve() and args.release.resolve()==RELEASE_PATH.resolve()
            and sha(PREP/"execution_adapter_frozen.py")==sha(SOURCE),
            "exact canonical frozen source review and first-wave release")
    check=read(args.review,args.review_sha256)
    x=read(args.release,args.release_sha256)
    expected={"source_sha256":sha(SOURCE),"protocol_sha256":args.protocol_sha256,"preflight_sha256":args.preflight_sha256,
              "profiles_sha256":PROFILE_SHA,"joined_bridge_input_admission_sha256":JOINED_SHA}
    source_gate(check,p,f,expected)
    want={**expected,"source_review_sha256":args.review_sha256}
    current=read(S.JOINED,JOINED_SHA)
    require(current["status"]==S.ADMISSION and current["spatial_workers_authorized"] is False
            and all(current[k]==v for k,v in base.items()),"same actual geometry and cached derivative joined proof")
    obs=read(NATIVE_FIRST/"os_exit_observation.json",x["first_wave_independent_observed_exit_sha256"])
    end=read(NATIVE_FIRST/"terminal.json",obs["terminal_sha256"])
    native=read(NATIVE_FIRST/"receipt.json",obs["receipt_sha256"])
    release_summary(x,obs,end,native,want)
    actual_src=IND_FIRST/"replay_global.py"
    require(native["source_sha256"]==end["source_sha256"]==sha(actual_src)==P.PINS[actual_src]
            and all(native[k]==end[k] for k in ("protocol_sha256","preflight_sha256",
                      "producer_root_terminal_admission_sha256","producer_reduction_input_admission_sha256"))
            and all(sha(IND_FIRST/(k+".json"))==native[k+"_sha256"] for k in ("protocol","preflight")),
            "same observed independent reviewed global replay phase")
    return allprofiles[args.wave],{**want,"later_waves_release_sha256":args.release_sha256,
                "first_wave_independent_receipt_sha256":obs["receipt_sha256"],
                "first_wave_independent_observed_exit_sha256":x["first_wave_independent_observed_exit_sha256"]}


def wave_folder(wave):
    require(wave in WAVES,"explicit canonical later wave")
    return RUNS/wave


def raw_gate(args):
    profile,bind=release_gate(args)
    path=wave_folder(args.wave)/"admission/raw_wave_admission.json"
    a=read(path,args.wave_admission_sha256)
    require(a["status"]==WAVE_REVIEW and V.yes(a["checks"]) and all(a[k]==v for k,v in bind.items())
            and a["wave_id"]==args.wave and a["targets"]==profile["raw_worker_targets"]
            and tuple(x["case_id"] for x in a["plans"])==tuple(profile["case_order"])
            and a["independent_component_hull_required"] is True and a["other_waves_authorized"] is False,
            "separate case-scoped worker authorization after independent selectors")
    return profile,a,{**bind,"wave_id":args.wave,"raw_wave_admission_sha256":args.wave_admission_sha256}


def selector_args(case):
    return W.bridge_args(case,JOINED_SHA)


def independent_plan_exit(row,case):
    folder=S.OUT/"independent_plan_review"/case
    h=row["independent_plan_replay_observed_os_exit_sha256"]
    doc=read(folder/"os_exit_observation.json",h)
    term=read(folder/"terminal.json",doc["terminal_sha256"])
    require(doc["status"]=="PASS_ROOT_OBSERVED_REMAINING32_INDEPENDENT_SELECTOR_REPLAY_EXIT"
            and doc["case_id"]==term["case_id"]==case and V.positive_int(doc["pid"]) and V.positive_int(term["pid"])
            and doc["pid"]==term["pid"] and V.positive_int(doc["unified_exec_session"])
            and V.zero(doc["observed_exit_code"]) and V.zero(term["exit_code"]) and term["status"]=="FINISHED"
            and doc["review_sha256"]==term["review_sha256"]==row["independent_plan_review_sha256"]
            and doc["selector_observed_os_sha256"]==row["observed_os_exit_sha256"]
            and row["independent_plan_review_file"]==rel(folder/"review.json")
            and all(doc[k]==term[k] for k in ("source_sha256","protocol_sha256","preflight_sha256")),
            "true independent selector plan replay before payload read")


def exact_context(a,case):
    row=V.A.pick(a["plans"],case)
    independent_plan_exit(row,case)
    checked,path=W.plan_docs(a,case)
    require(checked==row,"same wave-specific selector receipt")
    d,tails,proofs=S.open_reviewed_plan(selector_args(case),row["observed_os_exit_sha256"],path,row["independent_plan_review_sha256"])
    require(W.case_payload(d["entry"])==row["target"] and tuple(proofs)==S.OPS and len(tails)==2 and min(tails)>0,
            "exact later target and complete 168 same-scope selectors")
    return row,d,tails,proofs


def raw_folder(wave):
    return wave_folder(wave)/"raw"


def part_folder(wave,case,part):
    return raw_folder(wave)/"cases"/case/part


def native_raw(args):
    profile,a,bind=raw_gate(args)
    root=raw_folder(args.wave)
    start=read(root/"started.json",args.supervisor_started_sha256)
    job=[args.case_id,args.part]
    require(job in profile["part_jobs"] and start["jobs"]==profile["part_jobs"] and start["bindings"]==bind
            and V.positive_int(start["pid"]) and V.positive_int(args.supervisor_pid)
            and start["pid"]==args.supervisor_pid==os.getppid() and start["status"]=="RUNNING"
            and start["wave_id"]==args.wave,"only genuinely owned same-wave process")
    leaf=part_folder(args.wave,*job)
    mark={"pid":os.getpid(),"supervisor_pid":args.supervisor_pid,"wave_id":args.wave,"case_id":job[0],"part":job[1],"bindings":bind}
    save(leaf/"started.json",{"status":"RUNNING",**mark})
    term={"status":"FAILED","exit_code":1,**mark}
    components=[]
    try:
        selected,d,tails,proofs=exact_context(a,job[0])
        force,scaled=S.kernels()
        for i in range(COUNTS[job[1]]):
            proof=proofs[job[1],i]
            item=W.B.evaluate(d,proof,tails,force,scaled)
            nf,nh=W.validate_packet(item,proof,selected["target"],i,job[1],tails)
            file=f"{job[1]}_{i:03d}.json"
            digest=save(leaf/file,V.A.doc({**item,"worker_identity":mark,
                    "trial_box_sha256":selected["target"]["trial_box_sha256"],"selector_plan_sha256":selected["plan_sha256"]}))
            components.append({"kind":job[1],"slot":i,"file":file,"sha256":digest,"force_entries":nf,"upper_H_entries":nh})
        total={"components":len(components),**{k:sum(x[k] for x in components) for k in ("force_entries","upper_H_entries")}}
        h=save(leaf/"manifest.json",{"status":RAW_STATUS,**mark,"target":selected["target"],"plan":selected,
             "components":components,"counts":total,"point_force_scope":S.POINT,"whole_H_scope":S.BOX,
             "radius_exact":str(S.RADIUS),"root_asserted":False,"packet_independent_replay_required":True})
        term.update(status="FINISHED",exit_code=0,results_sha256=h)
    except BaseException as err:
        save(leaf/"failure.json",{**mark,"traceback":traceback.format_exc(),"error":str(err),"completed":components})
        raise
    finally:
        save(leaf/"terminal.json",term)

def part_identity(doc,base):
    return all(doc[k]==v for k,v in base.items()) and V.positive_int(doc["pid"]) and V.positive_int(doc["supervisor_pid"])


def check_part(pm,terminal,obs,identity,manifest_hash,selected):
    require(all(part_identity(x,identity) for x in (pm,terminal,obs)) and V.zero(obs["observed_exit_code"]),
            "true native owned child identity and exit")
    V.child_admit(obs["observed_exit_code"],pm,terminal,identity,manifest_hash)
    part=identity["part"]
    components=pm["components"]
    require(pm["status"]==RAW_STATUS and pm["target"]==selected["target"] and pm["plan"]==selected
            and pm["packet_independent_replay_required"] is True
            and pm["point_force_scope"]==S.POINT and pm["whole_H_scope"]==S.BOX
            and pm["radius_exact"]==str(S.RADIUS)
            and len(components)==pm["counts"]["components"]==COUNTS[part]
            and all(type(x["slot"]) is int and x["slot"]==i and x["kind"]==part and x["file"]==f"{part}_{i:03d}.json"
                and V.positive_int(x["force_entries"]) and V.positive_int(x["upper_H_entries"])
                and x["upper_H_entries"]==x["force_entries"]*(x["force_entries"]+1)//2
                for i,x in enumerate(components))
            and all(sum(x[k] for x in components)==pm["counts"][k] for k in ("force_entries","upper_H_entries")),
            "complete typed raw part component references and entry counts")


def command_flags(bind,args):
    return ["--protocol-sha256",bind["protocol_sha256"],"--preflight-sha256",bind["preflight_sha256"],
            "--review",str(args.review),"--review-sha256",bind["source_review_sha256"],
            "--release",str(args.release),"--release-sha256",bind["later_waves_release_sha256"],
            "--wave",args.wave]


def launch_part(args,a,job,pid,bind,started_sha):
    case,part=job
    path=part_folder(args.wave,case,part)
    path.mkdir(parents=True,exist_ok=False)
    base={"supervisor_pid":pid,"wave_id":args.wave,"case_id":case,"part":part,"bindings":bind}
    cmd=[sys.executable,"-B",str(SOURCE),"native-raw",*command_flags(bind,args),
         "--wave-admission-sha256",bind["raw_wave_admission_sha256"],"--case-id",case,"--part",part,
         "--supervisor-pid",str(pid),"--supervisor-started-sha256",started_sha]
    ident,code=V.run_owned_child(cmd,path,base,[])
    terminal=path/"terminal.json"
    t=sha(terminal) if terminal.is_file() else None
    observed={**ident,"observed_exit_code":code,"terminal_sha256":t}
    exit_hash=save(path/"observed_os_exit.json",observed)
    require(t is not None and V.zero(code),"failed native worker retained")
    term=read(terminal,t)
    m=term["results_sha256"]
    pm=read(path/"manifest.json",m)
    check_part(pm,term,observed,ident,m,V.A.pick(a["plans"],case))
    require(all(sha(path/x["file"])==x["sha256"] for x in pm["components"]),"all raw local signed packet bytes present")
    return {**observed,"observed_os_exit_sha256":exit_hash,"manifest_sha256":m,
            "manifest_file":rel(path/"manifest.json"),"counts":pm["counts"]}


def check_closed_prior_waves():
    """A hard interruption without the supervisor terminal blocks all later raw groups."""
    for wave in WAVES:
        here=raw_folder(wave)
        start=here/"started.json"
        end=here/"terminal.json"
        if start.is_file():
            require(end.is_file(),"unclosed earlier raw supervisor requires separately observed recovery")
            x,y=read(start,sha(start)),read(end,sha(end))
            require(x["status"]=="RUNNING" and y["status"] in ("FINISHED","FAILED")
                    and V.positive_int(x["pid"]) and V.positive_int(y["pid"])
                    and x["pid"]==y["pid"] and x["wave_id"]==y["wave_id"]==wave
                    and type(y["exit_code"]) is int,"prior ownership terminal does not match live-start identity")


def raw_supervisor(args):
    profile,a,bind=raw_gate(args)
    for c in profile["case_order"]:
        exact_context(a,c)
    with V.stage_lock(LOCK),W.deferred_interrupts() as interrupted:
        check_closed_prior_waves()
        here=raw_folder(args.wave)
        here.mkdir(parents=True,exist_ok=False)
        pid=os.getpid()
        mark={"pid":pid,"wave_id":args.wave,"case_order":profile["case_order"],"jobs":profile["part_jobs"],"bindings":bind}
        h=save(here/"started.json",{"status":"RUNNING",**mark})
        end={"status":"FAILED","exit_code":1,**mark}
        jobs=[tuple(x) for x in profile["part_jobs"]]
        try:
            parts=W.bounded(jobs,lambda item:launch_part(args,a,item,pid,bind,h),4,interrupted)
            require([(p["case_id"],p["part"]) for p in parts]==jobs
                    and len({p["pid"] for p in parts})==8,"complete eight distinct native exits")
            if interrupted():
                raise KeyboardInterrupt("later wave stopped after submitted owned tasks truly exited")
            rows=[{"case_id":case,"target":a["targets"][case],"parts":[p for p in parts if p["case_id"]==case]}
                  for case in profile["case_order"]]
            result={"status":"COMPLETE_RAW_LATER_WAVE_SPATIAL_BRIDGES_AWAITING_INDEPENDENT_REPLAY",
                    **mark,"cases":rows,"raw_components":672,"root_asserted":False,"independent_packet_replay_required":True}
            rh=save(here/"manifest.json",result)
            end.update(status="FINISHED",exit_code=0,manifest_sha256=rh)
            print({"raw_manifest_sha256":rh,"wave_id":args.wave,"native_supervisor_pid":pid},flush=True)
        except BaseException as err:
            receipts=(part_folder(args.wave,c,p)/name for c,p in jobs
                      for name in ("observed_os_exit.json","abnormal_child_reap.json","terminal.json","failure.json"))
            save(here/"failure.json",{**mark,"error":str(err),"traceback":traceback.format_exc(),
                  "closed_attempt_sha256":{rel(p):sha(p) for p in receipts if p.is_file()}})
            raise
        finally:
            save(here/"terminal.json",end)


def open_raw_phase(args,outer_sha):
    """Read-only typed raw metadata adapter for a separately approved independent replay."""
    p,a,bind=raw_gate(args)
    root=raw_folder(args.wave)
    d=read(root/"root_terminal_admission.json",outer_sha)
    require(d["status"]=="PASS_ROOT_OBSERVED_LATER_WAVE_RAW_SPATIAL_TERMINALS"
            and d["wave_id"]==args.wave and V.positive_int(d["pid"]) and V.positive_int(d["unified_exec_session"])
            and V.zero(d["observed_exit_code"]) and d["raw_wave_admission_sha256"]==args.wave_admission_sha256,
            "genuine root observed native raw phase exit")
    files=d["sha256"]
    expected={rel(root/(k+".json")) for k in ("started","terminal","manifest")}
    for c,k in p["part_jobs"]:
        folder=part_folder(args.wave,c,k)
        expected|={rel(folder/(x+".json")) for x in ("launched","started","terminal","observed_os_exit","manifest")}
        expected|={rel(folder/f"{k}_{i:03d}.json") for i in range(COUNTS[k])}
    require(set(files)==expected and len(files)==715,"complete externally hashed raw wave, no alien file")
    def load(path):
        return read(path,files[rel(path)])
    started,terminal,manifest=(load(root/(k+".json")) for k in ("started","terminal","manifest"))
    expected_parent={"pid":d["pid"],"wave_id":args.wave,"bindings":bind,"jobs":p["part_jobs"],"case_order":p["case_order"]}
    require(all(V.positive_int(x["pid"]) and all(x[k]==v for k,v in expected_parent.items())
           for x in (started,terminal,manifest))
           and started["status"]=="RUNNING" and terminal["status"]=="FINISHED" and V.zero(terminal["exit_code"])
           and terminal["manifest_sha256"]==files[rel(root/"manifest.json")]
           and manifest["status"]=="COMPLETE_RAW_LATER_WAVE_SPATIAL_BRIDGES_AWAITING_INDEPENDENT_REPLAY"
           and manifest["raw_components"]==672 and manifest["root_asserted"] is False
           and manifest["independent_packet_replay_required"] is True
           and [x["case_id"] for x in manifest["cases"]]==p["case_order"],"whole current raw phase is terminal")
    parts={}
    seen=set()
    for c in manifest["cases"]:
        case=c["case_id"]
        selected=V.A.pick(a["plans"],case)
        require(c["target"]==selected["target"] and tuple(x["part"] for x in c["parts"])==tuple(COUNTS),
                "same per-case two part mapping")
        for row in c["parts"]:
            kind=row["part"]
            folder=part_folder(args.wave,case,kind)
            launch,begin,term,osr,rawm=(load(folder/(k+".json")) for k in
                 ("launched","started","terminal","observed_os_exit","manifest"))
            ident={"pid":row["pid"],"supervisor_pid":d["pid"],"wave_id":args.wave,"case_id":case,"part":kind,"bindings":bind}
            check_part(rawm,term,osr,ident,row["manifest_sha256"],selected)
            require(row["pid"] not in seen and all(part_identity(x,ident) for x in (row,launch,begin))
                    and begin["status"]=="RUNNING" and V.zero(row["observed_exit_code"])
                    and osr["terminal_sha256"]==row["terminal_sha256"]==files[rel(folder/"terminal.json")]
                    and row["manifest_sha256"]==files[rel(folder/"manifest.json")]
                    and row["observed_os_exit_sha256"]==files[rel(folder/"observed_os_exit.json")]
                    and all(files[rel(folder/x["file"])]==x["sha256"] for x in rawm["components"]),
                    "replayed parent/child metadata and every raw digest reference")
            seen.add(row["pid"])
            parts[(case,kind)]=(folder,rawm)
    return p,a,parts,d


def selector_commands(args):
    p,_=release_gate(args)
    for c in p["case_order"]:
        a=selector_args(c)
        print(V.A.doc([sys.executable,"-B",str(S.SOURCE),"plan","--case-id",c,
             "--protocol-sha256",a.protocol_sha256,"--preflight-sha256",a.preflight_sha256,
             "--manifest-sha256",a.manifest_sha256,"--review",str(a.review),"--review-sha256",a.review_sha256,
             "--admission",str(a.admission),"--admission-sha256",a.admission_sha256]))

def fixtures(profiles):
    p=profiles[WAVES[0]]
    require(len(p["raw_worker_targets"])==4 and all(len(x["part_jobs"])==8 for x in profiles.values()),
            "all seven inherited scopes set before launch")
    expected={"source_sha256":"f"*64}
    x={"status":RELEASE,**expected,"checks":{"fixture":True},"case_order":list(W.FIRST),"following_waves":list(WAVES),
       "selector_invocations_only_under_this_release":True,"per_wave_plan_review_and_raw_admission_required":True,
       "first_wave_independent_receipt_sha256":"a"*64}
    osr={"status":FIRST_OS,"pid":13,"observed_exit_code":0,"unified_exec_session":19,"receipt_sha256":"a"*64}
    term={"status":"FINISHED","pid":13,"exit_code":0,"receipt_sha256":"a"*64}
    native={"status":"COMPLETE_INDEPENDENT_FIRST_WAVE_SIGNED_ENDPOINT_REPLAY","pid":13,
            "case_results":[{"case_id":c} for c in W.FIRST],"full_five_gate_complete":False}
    release_summary(x,osr,term,native,expected)
    require(all(V.A.reject(lambda:release_summary(x,bad,term,native,expected)) for bad in
               ({**osr,"observed_exit_code":False},{**osr,"pid":13.0},{**osr,"receipt_sha256":"b"*64}))
            and V.A.reject(lambda:release_summary(x,osr,term,{**native,"case_results":native["case_results"][:-1]},expected)),
            "incomplete or type-spoofed first-wave release cannot pass")
    # Reuse completed independent selector metadata validation without opening
    # their old leaf plans or any spatial outputs.
    first=read(W.WAVE_GATE,"8468da6cd555797298df182d6213a93a4d054029bdfca9b332aea80ca39f6855")
    for r in first["plans"]:
        independent_plan_exit(r,r["case_id"])
    require(V.A.reject(lambda:independent_plan_exit(first["plans"][0],first["plans"][1]["case_id"])),
            "can't borrow sibling plan proof")
    with TemporaryDirectory(dir=HERE) as scratch:
        base=Path(scratch)
        def inert(job):
            folder=base/job
            folder.mkdir()
            ident,code=V.run_owned_child([sys.executable,"-c","pass"],folder,
                             {"wave_id":"inert","part":"fixture","case_id":job,"supervisor_pid":os.getpid(),"bindings":{}},[])
            require(V.positive_int(ident["pid"]) and V.zero(code),"real harmless owned subprocess exit")
            return ident["pid"]
        with V.stage_lock(base):
            handles=W.bounded(["zero","one","two","three"],inert,4)
            require(len(set(handles))==4,"each inert subprocess genuinely individually owned")
            def locked():
                with V.stage_lock(base):
                    raise AssertionError("owner lost stage lock")
            require(V.A.reject(locked),"global scientific ownership lock remains held through native completion")
        with V.stage_lock(base):
            pass
    return {"all_seven_scope_profiles_same_28_targets_and_job_order":True,
            "reject_mixed_partial_or_alias_first_wave_release":True,
            "real_independent_selector_metadata_sibling_swap_rejected":True,
            "four_inert_owned_OS_process_exits_and_stage_lock":True,
            "no_new_scientific_kernels_invoked":True}


def prepare():
    profiles,base=cold()
    checks=fixtures(profiles)
    require(V.yes(checks) and not RUNS.exists() and not LOCK.exists(),"cold only, no actual wave output")
    PREP.mkdir(exist_ok=False)
    source=sha(SOURCE)
    protocol={"status":"FROZEN_REMAINING28_EXECUTION_ADAPTER","source_sha256":source,"profiles_sha256":PROFILE_SHA,
              "inherited_profiles_source_sha256":PROFILE_CODE_SHA,"joined_bridge_input_admission_sha256":JOINED_SHA,
              "root_independent_first_wave_release_required":True,"raw_per_wave_external_authority_required":True,
              "wave_ids":list(WAVES),"max_owned_native_scientific_processes":4,"bridge_binding":base,
              "downstream":("open_raw_phase supplies complete typed raw authority for existing local hull math; "
                            "future independent replay and signed reduction have their separate proof gates."),
              "per_wave_jobs":{k:v["part_jobs"] for k,v in profiles.items()},**ZERO}
    ph=save(PREP/"protocol.json",protocol)
    fh=save(PREP/"preflight.json",{"status":"PASS","source_sha256":source,"protocol_sha256":ph,
                                   "profiles_sha256":PROFILE_SHA,**ZERO,"checks":checks})
    with (PREP/"execution_adapter_frozen.py").open("xb") as f:
        f.write(SOURCE.read_bytes())
    print({"source_sha256":source,"protocol_sha256":ph,"preflight_sha256":fh},flush=True)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    subs=cli.add_subparsers(dest="stage",required=True)
    subs.add_parser("prepare")
    for name in ("selector-commands","raw","native-raw"):
        config=subs.add_parser(name)
        config.add_argument("--wave",choices=WAVES,required=True)
        for key in ("protocol","preflight","review","release"):
            config.add_argument("--"+key+"-sha256",required=True)
            if key in ("review","release"):
                config.add_argument("--"+key,required=True,type=lambda v:Path(v).resolve())
        if name!="selector-commands":
            config.add_argument("--wave-admission-sha256",required=True)
        if name=="native-raw":
            config.add_argument("--case-id",required=True)
            config.add_argument("--part",required=True,choices=tuple(COUNTS))
            config.add_argument("--supervisor-pid",required=True,type=int)
            config.add_argument("--supervisor-started-sha256",required=True)
    args=cli.parse_args()
    {"prepare":lambda _:prepare(),"selector-commands":selector_commands,"raw":raw_supervisor,"native-raw":native_raw}[args.stage](args)
