"""Read back all observed selectors for one wave. No selector or worker launch."""
import argparse
import json
from pathlib import Path
import sys
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import wave_execution as X

S,W,V=X.S,X.W,X.V
REPLAY=HERE.parent/"selector_replay_preparation"
EXEC={
    "protocol_sha256":"5057fecf76205bd0eee84cde3aff730b0011299d87c36fdebaad1be53cf10c72",
    "preflight_sha256":"aa03be5fd4d2341326a830bea1081e421b3739711d51b50d9398dbe6c764e214",
    "review":X.REVIEW_PATH,"review_sha256":"bf1af4575dffa36a7ddade7df0137156ec7a9ae610265e4f4d17386b9956149a",
    "release":X.RELEASE_PATH,"release_sha256":"002b471dceed8666aa0f70e1e9673744f49d0bf3a70da65a194c3faf34b0adbe",
}
NATIVE={
    "source_sha256":"b8cf7a1ec27f4a342473b40bdd8ec894f1183ab26c4a4bf2976115bbf52b6c63",
    "protocol_sha256":"bf9f215e8adbfbd4cb01444177d38e9ff6e877e77bbcf0f4268782b743adba04",
    "preflight_sha256":"e0c01b71bb7355e50b940d84d3babe32ef943ee227eee0fa3fb82096b6aba433",
    "root_source_review_sha256":"add734c73252220d81a2f6ee25780c475189a6972a1e986153553820a56d6364",
    "later_release_sha256":EXEC["release_sha256"],
}
SOURCE_PINS={
    X.SOURCE:"44c4fea65168a77926f93a08c3817edc0645c94ec7ef625fa34b4656beade8b7",
    REPLAY/"replay.py":NATIVE["source_sha256"],
    REPLAY/"protocol.json":NATIVE["protocol_sha256"],
    REPLAY/"preflight.json":NATIVE["preflight_sha256"],
    REPLAY/"root_source_review.json":NATIVE["root_source_review_sha256"],
    REPLAY/"independent_source_review/receipt.json":"2997c283627146b5447321dd3490c5bb0ab7e1950155a0c63cf0a7a661606dfc",
}
read,require,save,sha,rel=X.read,X.require,X.save,X.sha,X.rel


def cold(wave):
    a=SimpleNamespace(wave=wave,**EXEC)
    profile,bind=X.release_gate(a)
    require(all(sha(p)==h for p,h in SOURCE_PINS.items()),"unchanged independent and selector phase sources")
    phase=read(REPLAY/"root_source_review.json",NATIVE["root_source_review_sha256"])
    require(phase["status"]=="PASS_ROOT_REMAINING28_SELECTOR_REPLAY_SOURCE" and V.yes(phase["checks"])
            and all(phase[k]==NATIVE[k] for k in ("source_sha256","protocol_sha256","preflight_sha256","later_release_sha256"))
            and phase["independent_source_review_sha256"]==SOURCE_PINS[REPLAY/"independent_source_review/receipt.json"]
            and phase["new_spatial_integrations_authorized"] is False,"same separately reviewed later selector source")
    return profile,bind


def index(x,expected):
    return type(x["slot"]) is int and (x["kind"],x["slot"])==expected


def orders(leaf,kind):
    return {scope:(leaf[scope]["selector"]["degree"] if kind=="bulk" else leaf[scope]["order"]) for scope in (S.POINT,S.BOX)}


def pair(manifest,review,leaf,route,expected):
    require(all(index(x,expected) for x in (manifest,review,leaf,route))
            and manifest["file"]==f"{expected[0]}_{expected[1]:03d}.json"
            and manifest["sha256"]==review["checked_leaf_sha256"]
            and leaf["global_indices"]==route["global_indices"] and all(type(i) is int for i in leaf["global_indices"])
            and all(V.positive_int(v) for v in orders(leaf,expected[0]).values())
            and review["point_and_box_selected_orders"]==orders(leaf,expected[0]),
            "typed same projection, order and hash for each independently reviewed selector leaf")
    return S.E._digest(review)


def inspect(wave,case,ref_hash,targets,entry,bind):
    f=S.OUT/"independent_plan_review"/case
    out=read(f/"os_exit_observation.json",ref_hash)
    require(out["status"]=="PASS_ROOT_OBSERVED_REMAINING32_INDEPENDENT_SELECTOR_REPLAY_EXIT"
            and out["case_id"]==case and V.positive_int(out["pid"]) and V.positive_int(out["unified_exec_session"])
            and V.zero(out["observed_exit_code"]) and all(out[k]==NATIVE[k] for k in ("source_sha256","protocol_sha256","preflight_sha256")),
            "genuine literal same-case independent selector replay exit")
    end=read(f/"terminal.json",out["terminal_sha256"])
    review=read(f/"review.json",out["review_sha256"])
    start_review=f/"started.json"
    start_review_hash=sha(start_review)
    started=read(start_review,start_review_hash)
    common={**NATIVE,"case_id":case,"wave_id":wave,"pid":out["pid"],"observed_os_exit_sha256":out["selector_observed_os_sha256"]}
    require(all(all(z[k]==v for k,v in common.items()) and V.positive_int(z["pid"]) for z in (started,end,review))
            and started["status"]=="RUNNING" and end["status"]=="FINISHED" and V.zero(end["exit_code"])
            and end["review_sha256"]==out["review_sha256"]
            and review["status"]==S.PLAN_REVIEW and V.yes(review["checks"])
            and V.zero(review["force_H_AD_moment_or_solve_jobs"]) and review["spatial_worker_authorized_by_this_receipt"] is False,
            "complete native independent selector replay, unchanged wave and reviewed phase")
    observed=out["selector_observed_os_sha256"]
    p=S.PLAN/case
    raw_os=read(p/"root_os_exit_observation.json",observed)
    start=read(p/"started.json",raw_os["started_sha256"])
    term=read(p/"terminal.json",raw_os["terminal_sha256"])
    raw=read(p/"plan.json",raw_os["plan_sha256"])
    require(raw_os["started_sha256"]==review["selector_started_sha256"]
            and start["status"]=="RUNNING_SELECTOR_PLAN" and V.positive_int(start["pid"])
            and start["pid"]==raw_os["pid"]==review["observed_selector_pid"]
            and start["case_id"]==case and start["bindings"]==term["bindings"]==raw["bindings"]
            and raw_os["unified_exec_session"]==review["observed_selector_unified_exec_session"]
            and raw_os["plan_sha256"]==review["plan_sha256"] and raw_os["terminal_sha256"]==review["terminal_sha256"]
            and targets[case]==W.case_payload(entry),"same seeded target and source selector process")
    row={"case_id":case,"target":targets[case],"plan_sha256":raw_os["plan_sha256"],
         "terminal_sha256":raw_os["terminal_sha256"],"observed_os_exit_sha256":observed,
         "independent_plan_review_file":rel(f/"review.json"),"independent_plan_review_sha256":out["review_sha256"],
         "independent_plan_replay_observed_os_exit_sha256":ref_hash}
    checked,path=W.plan_docs({**bind,"targets":targets,"plans":[row]},case)
    X.independent_plan_exit(row,case)
    require(checked==row and path==f/"review.json" and review["independently_verified_component_count"]==168,
            "actual current adapter accepts this independent full plan")
    src=review["checked_input_sha256"]
    pins={rel(p/name):raw_os[k] for name,k in (("started.json","started_sha256"),("terminal.json","terminal_sha256"),("plan.json","plan_sha256"))}
    pins[rel(p/"root_os_exit_observation.json")]=observed
    require(all(src[name]==h for name,h in pins.items()),"independent replay bound every selector process artifact")
    leaf_review={}
    for desc,certificate,route,key in zip(raw["components"],review["case_component_review"],entry["component_catalogue"],S.OPS,strict=True):
        file=p/desc["file"]
        require(src[rel(file)]==desc["sha256"],"actual local packet in independent hash closure")
        leaf=read(file,desc["sha256"])
        pins[rel(file)]=desc["sha256"]
        leaf_review[desc["file"]]=pair(desc,certificate,leaf,route,key)
    require(len(leaf_review)==168 and len(pins)==172,"complete 168 exact independent component reviews")
    return row, {"case_id":case,"selector_file_sha256":pins,"embedded_independent_component_review_sha256":leaf_review,
                 "native_start_sha256":start_review_hash,"native_receipt_sha256":out["review_sha256"],"native_terminal_sha256":out["terminal_sha256"],
                 "observed_independent_exit_sha256":ref_hash,"selector_pid":raw_os["pid"],"independent_pid":out["pid"]}


def run(wave,obs):
    profile,bind=cold(wave)
    require(len(obs)==len(set(obs))==4,"four unique external observed replay SHA pins in wave order")
    model,_=W.bridge_manifest()
    targets=profile["raw_worker_targets"]
    rows,links=zip(*(inspect(wave,c,h,targets,model[c],bind) for c,h in zip(profile["case_order"],obs,strict=True)),strict=True)
    require(len({r["selector_pid"] for r in links})==len({r["independent_pid"] for r in links})==4,
            "four different selector and review process identities")
    folder=X.wave_folder(wave)/"admission"
    folder.mkdir(parents=True,exist_ok=True)
    doc={"status":X.WAVE_REVIEW,**bind,"wave_id":wave,"targets":targets,"plans":list(rows),
         "independent_component_hull_required":True,"other_waves_authorized":False,
         "independent_component_count":672,"review_source_sha256":sha(Path(__file__)),
         "source_files_sha256":{rel(k):v for k,v in SOURCE_PINS.items()},"readback_cases":list(links),
         "checks":{"same_frozen_later_profiles_and_exact_original_target_identities":True,
                   "all_four_genuine_independent_and_producer_selector_exits_verified":True,
                   "all_672_hashes_and_independent_component_scopes_orders_maps_current":True,
                   "new_raw_only_in_this_wave_and_separate_independent_components_required":True},
         "selector_spatial_or_mechanical_jobs_by_admission":0}
    h=save(folder/"raw_wave_admission.json",doc)
    print(json.dumps({"wave_id":wave,"raw_wave_admission_sha256":h,"components":672}),flush=True)


def prepare():
    waves={w:cold(w)[0]["case_order"] for w in X.WAVES}
    require(len({c for group in waves.values() for c in group})==28,"all 28 separate targets")
    desc={"kind":"bulk","slot":4,"file":"bulk_004.json","sha256":"a"*64}
    leaf={**desc,"global_indices":[11,5],S.POINT:{"selector":{"degree":12}},S.BOX:{"selector":{"degree":24}}}
    review={**desc,"checked_leaf_sha256":"a"*64,"point_and_box_selected_orders":{S.POINT:12,S.BOX:24}}
    assert pair(desc,review,leaf,leaf,("bulk",4))
    reject=V.A.reject
    require(reject(lambda:pair(desc,review,{**leaf,"slot":4.0},leaf,("bulk",4)))
            and reject(lambda:pair(desc,{**review,"checked_leaf_sha256":"b"*64},leaf,leaf,("bulk",4)))
            and reject(lambda:pair(desc,review,leaf,{**leaf,"global_indices":[5,11]},("bulk",4))),
            "typed component alias, changed selector hash and changed projection rejected")
    out={"status":"PASS_COLD_LATER_WAVE_RAW_ADMISSION_HELPER","source_sha256":sha(Path(__file__)),
         "source_files_sha256":{rel(k):v for k,v in SOURCE_PINS.items()},"wave_case_order":waves,
         "checks":{"all_seven_cold_authority_and_latest_release":True,"typed_local_manifest_order_hash_and_map_guards":True},
         "actual_plan_payloads_opened":0,"actual_new_scientific_jobs":0}
    print(save(HERE/"preflight.json",out),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    sub=parser.add_subparsers(dest="action",required=True)
    sub.add_parser("prepare")
    actual=sub.add_parser("run")
    actual.add_argument("--wave",required=True,choices=X.WAVES)
    actual.add_argument("--observed-independent-sha256",required=True,nargs=4)
    args=parser.parse_args()
    prepare() if args.action=="prepare" else run(args.wave,args.observed_independent_sha256)
