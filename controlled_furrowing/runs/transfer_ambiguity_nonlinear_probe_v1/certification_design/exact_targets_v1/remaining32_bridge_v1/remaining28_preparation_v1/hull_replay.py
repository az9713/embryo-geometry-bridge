"""Independent local component replay for one exact case from an admitted later wave."""
import argparse
import importlib.util
from itertools import combinations_with_replacement
from pathlib import Path
import os
import sys
from time import monotonic
import traceback
from types import SimpleNamespace

import wave_execution as X

S,W,V=X.S,X.W,X.V
SOURCE=Path(__file__).resolve()
HERE=SOURCE.parent
PREP=HERE/"hull_replay_preparation"
REVIEW=PREP/"independent_prelaunch/review.json"
HULL_SOURCE=S.OUT/"independent_component_replay_preparation/hull_math.py"
LOWER=S.OUT.parent/"bridge_v1/independent_actual_preparation/replay_math.py"
PINS={
 X.SOURCE:"44c4fea65168a77926f93a08c3817edc0645c94ec7ef625fa34b4656beade8b7",
 X.PREP/"protocol.json":"5057fecf76205bd0eee84cde3aff730b0011299d87c36fdebaad1be53cf10c72",
 X.PREP/"preflight.json":"aa03be5fd4d2341326a830bea1081e421b3739711d51b50d9398dbe6c764e214",
 X.REVIEW_PATH:"bf1af4575dffa36a7ddade7df0137156ec7a9ae610265e4f4d17386b9956149a",
 X.RELEASE_PATH:"002b471dceed8666aa0f70e1e9673744f49d0bf3a70da65a194c3faf34b0adbe",
 HULL_SOURCE:"f8b4e1e8b4fb9069d235303f3fd0a45ad629ccd9d0bf37bf30ceca509c96b38b",
 LOWER:"e7bb66066d54e131fd009d2bf48e94f3f1fb7bc1955023f66bc96cbe3c6d78e9",
 S.OUT/"independent_component_replay_preparation/root_source_review.json":
     "987577c3b4eaa0c7dce0655a04a88ece28bf0e891fbc84fb573acd26db95edec",
 S.OUT/"admission/root_hull_design_review.json":
     "6442253147ebc1a73dff1358f8a44ebba8425ddb3e0af18ba01fb3b6a8fefc03",
}
STATUS="FROZEN_LATER_WAVE_INDEPENDENT_COMPONENT_HULL_WRAPPER"
SOURCE_STATUS="PASS_INDEPENDENT_LATER_WAVE_COMPONENT_HULL_ADAPTER"
RAW_PASS="PASS_INDEPENDENT_ACTUAL_LATER_WAVE_COMPONENT_HULLS"
RAW_HOLD="HOLD_TRUE_DISJOINT_LATER_WAVE_COMPONENT_DIAGNOSTIC"
SCOPE="Own independently valid signed spatial/cached component intervals on the same trial, expanded to exact full union hulls. Point F and entire-cube H separate. Disjoint anywhere is HOLD."
CHANNELS=("point_force","whole_box_H","point_same_Y_cache_correction")
ZERO={"actual_packets_read":0,"new_moments_or_spatial_jobs":0,"new_AD_or_solves":0}
read,sha,save,require=X.read,X.sha,X.save,X.require
rel=X.rel


def import_pinned(path,name):
    require(sha(path)==PINS[path],"same frozen local mathematics")
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def math(args=None):
    local=import_pinned(HULL_SOURCE,"_later_wave_frozen_hull")
    if args is None:
        return local
    lower=import_pinned(LOWER,"_later_wave_frozen_primitives")
    require(lower.E is S.E and lower.T is S.E.T and lower.RADIUS==S.RADIUS
            and (lower.POINT,lower.BOX)==(S.POINT,S.BOX),"same explicit lower mathematical primitives")
    return lower,local


def release_args(wave):
    return SimpleNamespace(wave=wave,protocol_sha256=PINS[X.PREP/"protocol.json"],
         preflight_sha256=PINS[X.PREP/"preflight.json"],review=X.REVIEW_PATH,review_sha256=PINS[X.REVIEW_PATH],
         release=X.RELEASE_PATH,release_sha256=PINS[X.RELEASE_PATH])


def static():
    require(all(sha(p)==h for p,h in PINS.items()),"reviewed independent local source and current release")
    profiles,b=X.cold()
    # Read the immutable first-wave release and native closure through its existing
    # authenticated reader. This does not dispatch a wave or open a later profile payload.
    p,release=X.release_gate(release_args(X.WAVES[0]))
    require(p==profiles[X.WAVES[0]] and all(len(p["case_order"])==4 for p in profiles.values())
            and len({c for p in profiles.values() for c in p["case_order"]})==28,
            "all seven exact later profiles and first-wave reference closure")
    prior=read(S.OUT/"independent_component_replay_preparation/root_source_review.json",
               PINS[S.OUT/"independent_component_replay_preparation/root_source_review.json"])
    memo=read(S.OUT/"admission/root_hull_design_review.json",PINS[S.OUT/"admission/root_hull_design_review.json"])
    require(prior["status"]=="PASS" and V.yes(prior["checks"])
            and prior["local_bound_and_hull_source_sha256"]==PINS[HULL_SOURCE]
            and memo["status"]=="PASS_ROOT_CONDITIONAL_INDEPENDENT_ENCLOSURE_HULL_DESIGN"
            and V.yes(memo["checks"]) and memo["actual_component_jobs_authorized"] is False,
            "same independently approved mathematical own-interval/union proof")
    return profiles,release


def wave_scope(profile,case):
    require(type(case) is str and profile["case_order"].count(case)==1 and
            [p for p in profile["part_jobs"] if p[0]==case]==[[case,"bulk"],[case,"surface"]],
            "known same-wave target and both exact part roles")
    return profile["raw_worker_targets"][case]


def review_fields(p,check,pre,b):
    require(check["status"]==SOURCE_STATUS and V.yes(check["checks"]) and all(check[k]==v for k,v in b.items())
            and check["permitted_waves"]==list(X.WAVES)
            and check["requires_external_raw_observed_exit_for_each_wave"] is True
            and check["new_raw_spatial_workers_authorized"] is False
            and pre["status"]=="PASS" and V.yes(pre["checks"]) and all(V.zero(pre[k]) for k in ZERO)
            and pre["protocol_sha256"]==b["protocol_sha256"]
            and p["status"]==STATUS and p["source_sha256"]==pre["source_sha256"]==b["source_sha256"],
            "separate current independent hull replay source review only")


def gate(args):
    profiles,release=static()
    p=read(PREP/"protocol.json",args.protocol_sha256)
    f=read(PREP/"preflight.json",args.preflight_sha256)
    require(args.review.resolve()==REVIEW.resolve() and sha(PREP/"hull_replay_frozen.py")==sha(SOURCE),
            "canonical independently reviewed later wave hull source")
    r=read(args.review,args.review_sha256)
    b={"source_sha256":sha(SOURCE),"protocol_sha256":args.protocol_sha256,"preflight_sha256":args.preflight_sha256,
       "later_waves_release_sha256":release["later_waves_release_sha256"],
       "wave_execution_source_sha256":PINS[X.SOURCE],"local_bound_and_hull_source_sha256":PINS[HULL_SOURCE]}
    review_fields(p,r,f,b)
    require(p["wave_case_order"]=={n:profiles[n]["case_order"] for n in X.WAVES}
            and p["static_source_sha256"]=={rel(path):h for path,h in PINS.items()},
            "current cold scope and local math bytes")
    profile=profiles[args.wave]
    wave_scope(profile,args.case_id)
    return profile,{**b,"hull_adapter_review_sha256":args.review_sha256}


def exact_worker_fields(wave,case,kind,packet,part,row,proof):
    fields=("pid","supervisor_pid","wave_id","case_id","part","bindings")
    own=packet["worker_identity"]
    require(type(packet["slot"]) is type(proof["slot"]) is int and packet["slot"]==proof["slot"]
            and packet["kind"]==proof["kind"]==kind and packet["global_indices"]==proof["global_indices"]
            and packet["record_id"]==row["target"]["record_id"]
            and part["case_id"]==own["case_id"]==case and part["part"]==own["part"]==kind
            and part["wave_id"]==own["wave_id"]==wave
            and all(own[k]==part[k] for k in fields) and V.positive_int(own["pid"]) and V.positive_int(own["supervisor_pid"])
            and packet["trial_box_sha256"]==row["target"]["trial_box_sha256"] and packet["selector_plan_sha256"]==row["plan_sha256"],
            "same raw native, wave/case/part, independent selector and target")


def require_hull_shape(proof,own,local,kind,slot):
    """Mechanical validation only. Interval validity is from the pinned independent local computation."""
    need=set(own["signed_channels"])==set(CHANNELS) and own["kind"]==kind and own["slot"]==slot
    require(need and type(own["slot"]) is int and own["global_indices"]==proof["global_indices"]
            and own["checked_cache_sites"]==(64 if kind=="bulk" else 12)
            and own["root_asserted"] is False and own["only_independently_valid_enclosures_used_for_hull_truth"] is True
            and own["overlap_not_used_as_proof"] is True and own["point_force_scope"]==S.POINT
            and own["whole_H_scope"]==S.BOX and S.Q(own["radius_exact"])==S.RADIUS,"preserve exact signed hull scopes")
    found=False
    for k in CHANNELS:
        ch=own["signed_channels"][k]
        expected=len(proof["global_indices"])
        if k=="whole_box_H":
            expected=expected*(expected+1)//2
        require(len(ch["rows"])==expected,"complete three signed channels and upper symmetry")
        expected_misses=[]
        expected_disjoint=[]
        n=len(proof["global_indices"])
        index=(list(map(list,combinations_with_replacement(range(n),2))) if k=="whole_box_H" else list(range(n)))
        for row,i in zip(ch["rows"],index,strict=True):
            got=row["index"]
            require((type(got) is int if type(i) is int else type(got) is list and all(type(j) is int for j in got))
                    and got==i,"typed ordered signed channel index")
            v,nested,separate=local.hull_and_compare(row["independent_exact"],row["producer_exact"])
            require(tuple(map(S.Q,row["hull_exact"]))==v and row["producer_contains_independent"] is nested
                    and row["disjoint"] is separate,"exact full union and typed diagnostics; no endpoint narrowing")
            if separate:
                expected_disjoint.append(row["index"])
            if not nested:
                expected_misses.append(row["index"])
        require(ch["disjoint_indices"]==expected_disjoint and ch["producer_containment_misses"]==expected_misses,
                "all forced disjoint/containment flags agree with exact endpoints")
        found|=bool(expected_disjoint)
    require(own["status"]==("HOLD_DISJOINT_SIGNED_ENCLOSURES" if found else "PASS_INDEPENDENT_COMPONENT_WITH_EXACT_HULLS"),
            "disjointness independently forces HOLD")

def context_args(wave,raw_admission_sha):
    a=release_args(wave)
    a.wave_admission_sha256=raw_admission_sha
    return a


def run(args):
    preprofile,binding=gate(args)
    xa=context_args(args.wave,args.raw_wave_admission_sha256)
    profile,a,parts,root_raw=X.open_raw_phase(xa,args.raw_root_observed_exit_sha256)
    require(profile==preprofile,"native raw wave must match independently reviewed profile")
    case=args.case_id
    target=wave_scope(profile,case)
    row,d,tails,selectors=X.exact_context(a,case)
    require(row["target"]==target and tuple(selectors)==S.OPS and len(tails)==2 and min(tails)>0,
            "same independently reviewed all168 selectors, geometry box, active material and tails")
    # Whole externally observed wave is authenticated before any component math.
    raw_files=root_raw["sha256"]
    for name,h in raw_files.items():
        require(sha(X.canonical(name))==h,"exact root-observed complete later-wave raw bytes")
    lower,local=math(args)
    _,scaled=S.kernels()
    require(S.E.core.ctx.prec==640,"same context precision before actual coefficient or moment calculation")
    lower.ad.ctx.prec=640
    folder=X.wave_folder(args.wave)/"independent_hulls"/case
    folder.mkdir(parents=True,exist_ok=False)
    pid=os.getpid()
    phase={**binding,"pid":pid,"case_id":case,"wave_id":args.wave,
           "raw_root_observed_exit_sha256":args.raw_root_observed_exit_sha256,
           "raw_wave_admission_sha256":args.raw_wave_admission_sha256}
    end={"status":"FAILED","exit_code":1,**phase}
    save(folder/"started.json",{"status":"RUNNING",**phase})
    t0=monotonic()
    completed=[]
    misses=dict.fromkeys(CHANNELS,0)
    holds=dict.fromkeys(CHANNELS,0)
    try:
        for kind,slot in S.OPS:
            rawbase,part=parts[(case,kind)]
            meta=part["components"][slot]
            rawpath=rawbase/meta["file"]
            s=selectors[(kind,slot)]
            packet=read(rawpath,meta["sha256"])
            W.validate_packet(packet,s,target,slot,kind,tails)
            exact_worker_fields(args.wave,case,kind,packet,part,row,s)
            S.E.core.A.cache_clear()
            lower.ad.A.cache_clear()
            own=local.local(lower,d,kind,slot,s,packet,tails,scaled.S)
            require_hull_shape(s,own,local,kind,slot)
            savepath=folder/meta["file"]
            h=save(savepath,{**own,"phase":phase,"packet_sha256":meta["sha256"],
                          "selector_leaf_content_sha256":W.B.digest(s)})
            completed.append({"kind":kind,"slot":slot,"native_file":rel(savepath),"native_sha256":h,
                "packet_file":rel(rawpath),"packet_sha256":meta["sha256"],"status":own["status"],
                "global_indices":own["global_indices"]})
            for k,v in own["signed_channels"].items():
                misses[k]+=len(v["producer_containment_misses"])
                holds[k]+=len(v["disjoint_indices"])
            print({"wave_id":args.wave,"case_id":case,"kind":kind,"slot":slot,"status":own["status"]},flush=True)
        require(tuple((x["kind"],x["slot"]) for x in completed)==S.OPS,
                "all 168 whole point/H component files in model order")
        held=any(holds.values())
        record={"status":RAW_HOLD if held else RAW_PASS,**phase,"scope":SCOPE,"target":target,
            "record_id":target["record_id"],"case_components":completed,
            "producer_observed_pid":root_raw["pid"],"producer_observed_unified_exec_session":root_raw["unified_exec_session"],
            "raw_root_observed_mapping_sha256":S.E._digest(raw_files),
            "raw_root_observed_component_files":len(raw_files),"selector_plan_sha256":row["plan_sha256"],
            "plan_review_sha256":row["independent_plan_review_sha256"],"trial_box_sha256":target["trial_box_sha256"],
            "model_contract_sha256":target["model_contract_sha256"],"exact_drive":target["exact_drive"],
            "cache_packet_sha256":target["point_cache_sha256"],"R_binary_sha256":target["R_binary_sha256"],
            "center_binary_sha256":target["center_binary_sha256"],"original_archived_shift_exact":target["archived_center_shift_exact"],
            "original_addition_gap_exact":target["original_addition_gap_exact"],"observer_goal_exact":target["observer_goal_exact"],
            "original_observer_field_sha256":target["original_observer_field_sha256"],"radius_exact":str(S.RADIUS),
            "producer_containment_miss_counts":misses,"disjoint_counts":holds,"full_independent_hulls_computed":True,
            "admit_hull_only_global_reduction_after_separate_review":not held,
            "overlap_or_producer_containment_used_as_proof":False,"root_asserted":False,"seconds":monotonic()-t0}
        h=save(folder/"receipt.json",record)
        end.update(status="FINISHED_WITH_HOLD" if held else "FINISHED",exit_code=2 if held else 0,receipt_sha256=h)
        print({"receipt_sha256":h,"disjoint_counts":holds,"containment_misses":misses},flush=True)
        return end["exit_code"]
    except BaseException as err:
        save(folder/"failure.json",{"error":str(err),"traceback":traceback.format_exc(),"completed":completed,**phase})
        raise
    finally:
        save(folder/"terminal.json",end)


def fixtures(profiles):
    hm=math()
    zero=(0,1)
    overlap=(1,2)
    disjoint=(2,3)
    require(hm.hull_and_compare(zero,overlap)==((0,2),False,False)
            and hm.hull_and_compare(zero,disjoint)==((0,3),False,True),"whole overlap and genuine disjoint treated differently")
    p={"kind":"bulk","slot":0,"global_indices":[5]}
    def channel(raw):
        return {"point_force":hm.labelled_channel({0:zero},{0:raw}),
                "whole_box_H":hm.labelled_channel({(0,0):zero},{(0,0):overlap}),
                "point_same_Y_cache_correction":hm.labelled_channel({0:zero},{0:overlap})}
    a={**p,"checked_cache_sites":64,"root_asserted":False,"only_independently_valid_enclosures_used_for_hull_truth":True,
       "overlap_not_used_as_proof":True,"point_force_scope":S.POINT,"whole_H_scope":S.BOX,"radius_exact":str(S.RADIUS)}
    good={**a,"signed_channels":channel(overlap),"status":"PASS_INDEPENDENT_COMPONENT_WITH_EXACT_HULLS"}
    require_hull_shape(p,good,hm,"bulk",0)
    held={**a,"signed_channels":channel(disjoint),"status":"HOLD_DISJOINT_SIGNED_ENCLOSURES"}
    require_hull_shape(p,held,hm,"bulk",0)
    shrunk={**a,"signed_channels":channel(overlap),"status":good["status"]}
    shrunk["signed_channels"]["point_force"]["rows"][0]["hull_exact"]=["1","1"]
    require(V.A.reject(lambda:require_hull_shape(p,shrunk,hm,"bulk",0))
            and V.A.reject(lambda:require_hull_shape(p,{**held,"status":good["status"]},hm,"bulk",0))
            and V.A.reject(lambda:require_hull_shape(p,{**good,"whole_H_scope":S.POINT},hm,"bulk",0)),
            "narrowing substitution, disjoint PASS and scope swap all forbidden")
    for wave,profile in profiles.items():
        for c in profile["case_order"]:
            wave_scope(profile,c)
        require(V.A.reject(lambda:wave_scope(profile,W.FIRST[0])),"old case foreign to every later wave")
    ident={"pid":67,"supervisor_pid":23,"wave_id":"wave_02","case_id":"synthetic","part":"bulk","bindings":{"scope":"same"}}
    row={"target":{"record_id":"synthetic-r","trial_box_sha256":"c"*64},"plan_sha256":"d"*64}
    raw={**p,"record_id":"synthetic-r","worker_identity":ident,"trial_box_sha256":"c"*64,"selector_plan_sha256":"d"*64}
    exact_worker_fields("wave_02","synthetic","bulk",raw,ident,row,p)
    require(all(V.A.reject(f) for f in (
        lambda:exact_worker_fields("wave_03","synthetic","bulk",raw,ident,row,p),
        lambda:exact_worker_fields("wave_02","synthetic","surface",raw,ident,row,p),
        lambda:exact_worker_fields("wave_02","synthetic","bulk",{**raw,"worker_identity":{**ident,"pid":67.0}},ident,row,p),
        lambda:exact_worker_fields("wave_02","synthetic","bulk",{**raw,"trial_box_sha256":"b"*64},ident,row,p),
    )), "raw worker identity and physical scope guards")
    return {"all_seven_explicit_profiles_and_28_case_parts":True,
            "all_three_union_diagnostics_and_exact_hold_transition":True,
            "reject_narrowed_endpoint_wrong_H_scope_and_fake_disjoint_pass":True,
            "reject_foreign_wave_part_PID_alias_and_box_digest":True,"cold_no_actual_raw_or_moments_read":True}


def prepare():
    profiles,b=static()
    checks=fixtures(profiles)
    require(V.yes(checks),"cold later wave hull wrapper")
    PREP.mkdir(parents=True,exist_ok=False)
    phase={"source_sha256":sha(SOURCE),"later_waves_release_sha256":b["later_waves_release_sha256"],
           "wave_execution_source_sha256":PINS[X.SOURCE],"local_bound_and_hull_source_sha256":PINS[HULL_SOURCE]}
    p=save(PREP/"protocol.json",{"status":STATUS,**phase,"scope":SCOPE,"static_source_sha256":{rel(k):h for k,h in PINS.items()},
                  "wave_case_order":{k:profiles[k]["case_order"] for k in X.WAVES},"case_components":168,"raw_wave_files":715,
                  "requires_separate_root_external_observation_for_each_raw_wave":True,
                  "outer_observation_of_this_replay_still_required":True,"raw_workers_or_later_global_jobs_authorized":False,**ZERO})
    f=save(PREP/"preflight.json",{"status":"PASS",**phase,"protocol_sha256":p,"checks":checks,**ZERO})
    with (PREP/"hull_replay_frozen.py").open("xb") as stream:
        stream.write(SOURCE.read_bytes())
    print({**phase,"protocol_sha256":p,"preflight_sha256":f},flush=True)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    sub=cli.add_subparsers(dest="stage",required=True)
    sub.add_parser("prepare")
    actual=sub.add_parser("run")
    actual.add_argument("--wave",required=True,choices=X.WAVES)
    actual.add_argument("--case-id",required=True)
    for k in ("protocol","preflight","review","raw-wave-admission","raw-root-observed-exit"):
        actual.add_argument("--"+k+"-sha256",required=True)
        if k=="review":
            actual.add_argument("--review",type=lambda v:Path(v).resolve(),required=True)
    args=cli.parse_args()
    if args.stage=="prepare":
        prepare()
    else:
        sys.exit(run(args))
