"""Admit first-wave raw work after four independently replayed selector plans."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
BR=HERE.parent
ROOT=next(p for p in BR.parents if (p/"PLAN.md").is_file())
GP=BR/"worker_preparation"
GS=ROOT/"scripts/transfer_ambiguity_remaining32_spatial_worker.py"
IR=BR/"independent_plan_review"
FIRST=("max_contrast_change_small_01","max_contrast_change_small_13",
       "max_negative_contrast_change_small_01","max_negative_contrast_change_small_16")
ORDER=tuple((kind,i) for kind,n in (("bulk",144),("surface",24)) for i in range(n))
R=Q(1,2**34)
EXPECTED={
    "bridge_manifest_sha256":"d42b7260283ba02ccc34b9e9ba69ce76d58507115c22615c90d28fa38020cb41",
    "bridge_cold_review_sha256":"89f9b5cacf89b4055fc404c00d3e3c42cf10fed20337fe455e996c2944ef73b2",
}
JOIN="2e2e388e798ed509bb52bc72a6169b8b6ee5db55173b3529ac73c3812b23a6cb"
SELECTOR_PHASE={
    "source_sha256":"c4b096ba3d93ffa17abe5eb5a0b20bec2bd41a0c97ddb0fe0b13ed9c99616e9b",
    "protocol_sha256":"6e5b866df43c30147a97d4878c13f8680e01f954c41be6a29f532948de129263",
    "preflight_sha256":"957647263377fea5be85333f4a2f6cf54faaecc7888e125c0e9e7d1143ba86fb",
}
ROOT_SELECTOR_APPROVAL="da4939cbad4c64c7ad9a0b3d7dd909dae2bb87e9f76e53d665122bb33194dfc5"
RAW=(
    "369561fd90514458ea4d1cac8e0f2a1a42d371923f6baa7b5485f05195969f31",
    "39c5926add26212a638568c486f42ff659ee3d8dcec31440d1d98502a976d71c",
    "d241dc50639c877e0f7d81586ebc4551555b4e347292c2836992c3aef5299176",
    "374be4e910dab74c406efc78c100abf0684b8bd99154ce86d2610eb18781cc0a",
)
BASIS={
    BR/"preparation/manifest.json":EXPECTED["bridge_manifest_sha256"],
    BR/"independent_prelaunch/review.json":EXPECTED["bridge_cold_review_sha256"],
    HERE/"bridge_input_admission.json":JOIN,
    IR/"replay.py":SELECTOR_PHASE["source_sha256"],
    IR/"protocol.json":SELECTOR_PHASE["protocol_sha256"],
    IR/"preflight.json":SELECTOR_PHASE["preflight_sha256"],
    IR/"root_source_review.json":ROOT_SELECTOR_APPROVAL,
}
HASHES={}


def demand(ok,why):
    if not ok:
        raise ValueError(why)


def sha(v):
    return sha256(v).hexdigest()


def rel(path):
    return path.resolve().relative_to(ROOT).as_posix()


def read(path,digest,parse=True):
    path=(ROOT/path).resolve()
    key=rel(path)
    raw=path.read_bytes()
    h=sha(raw)
    demand(h==digest==HASHES.get(key,h),f"changed {key}")
    HASHES[key]=h
    return json.loads(raw) if parse else raw


def yes(d):
    return type(d) is dict and bool(d) and all(v is True for v in d.values())


def zero(v):
    return type(v) is int and v==0


def pid(v):
    return type(v) is int and v>0


def one(rows,k,v):
    item,=[a for a in rows if a[k]==v]
    return item


def dimensions(entry):
    tr=entry["trial"]
    fixed=("trial_box_sha256","model_contract_sha256","point_cache_sha256","R_binary_sha256",
           "center_binary_sha256","archived_center_shift_exact","exact_drive","energy_scale_binary_exact")
    arch=("case_id","record_id","original_observer_field_sha256","original_addition_gap_exact","observer_goal_exact")
    return {**{k:tr[k] for k in fixed},**{k:entry[k] for k in arch}}


def actual_selector(case,root_hash,own_hash,entry,joined):
    p=BR/"plans"/case
    q=IR/case
    outer=read(p/"root_os_exit_observation.json",root_hash)
    demand(outer["status"]=="PASS_OBSERVED_REMAINING32_SELECTOR_TRUE_OS_EXIT"
           and outer["case_id"]==case and zero(outer["observed_exit_code"])
           and pid(outer["pid"]) and pid(outer["unified_exec_session"]),"no root true OS selector admission")
    plan=read(p/"plan.json",outer["plan_sha256"])
    term=read(p/"terminal.json",outer["terminal_sha256"])
    demand(plan["status"]=="COMPLETE_REMAINING32_CASE_SELECTOR_PLAN_AWAITS_INDEPENDENT_REVIEW"
           and term["status"]=="FINISHED" and zero(term["exit_code"])
           and term["plan_sha256"]==outer["plan_sha256"] and zero(plan["new_spatial_integrations"])
           and plan["spatial_workers_authorized"] is False
           and all(x["case_id"]==case and x["pid"]==outer["pid"] and pid(x["pid"]) for x in (plan,term))
           and plan["trial_box_sha256"]==joined["trial_box_sha256"]==entry["trial"]["trial_box_sha256"]
           and plan["record_id"]==term["record_id"]==entry["record_id"]
           and plan["whole_box_radius_exact"]==str(R) and plan["bindings"]==term["bindings"]
           and plan["bindings"]["bridge_input_admission_sha256"]==JOIN,
           "wrong selector terminal, same box or joined identity")
    independent=read(q/"os_exit_observation.json",own_hash)
    demand(independent["status"]=="PASS_ROOT_OBSERVED_REMAINING32_INDEPENDENT_SELECTOR_REPLAY_EXIT"
           and independent["case_id"]==case and zero(independent["observed_exit_code"])
           and pid(independent["pid"]) and pid(independent["unified_exec_session"]),
           "independently replayed selector must have genuinely observed exit0")
    end=read(q/"terminal.json",independent["terminal_sha256"])
    result=read(q/"review.json",independent["review_sha256"])
    demand(result["status"]=="PASS_INDEPENDENT_REMAINING32_ACTUAL_SELECTOR_PLAN" and yes(result["checks"])
           and end["status"]=="FINISHED" and zero(end["exit_code"]) and end["review_sha256"]==independent["review_sha256"]
           and all(x["pid"]==independent["pid"] and pid(x["pid"]) and x["case_id"]==case for x in (end,result))
           and all(end[k]==v==result[k] for k,v in SELECTOR_PHASE.items())
           and result["observed_os_exit_sha256"]==end["observed_os_exit_sha256"]==root_hash
           and result["bindings"]==plan["bindings"] and result["record_id"]==plan["record_id"]
           and result["plan_sha256"]==outer["plan_sha256"] and result["terminal_sha256"]==outer["terminal_sha256"]
           and type(result["independently_verified_component_count"]) is int
           and result["independently_verified_component_count"]==168
           and zero(result["force_H_AD_moment_or_solve_jobs"]) and result["spatial_worker_authorized_by_this_receipt"] is False,
           "missing same true independent scoped selectors")
    items=plan["components"]
    reviewed=result["case_component_review"]
    for seq in (items,reviewed):
        demand(tuple((x["kind"],x["slot"]) for x in seq)==ORDER and all(type(x["slot"]) is int for x in seq),
               "partial or mixed selector leaf listing")
    for raw,leaf,route in zip(items,reviewed,entry["component_catalogue"],strict=True):
        k,i=raw["kind"],raw["slot"]
        f=p/f"{k}_{i:03d}.json"
        digest=raw["sha256"]
        proof=read(f,digest)
        demand(raw["file"]==f.name and digest==leaf["checked_leaf_sha256"]==result["checked_input_sha256"][rel(f)]
               and proof["global_indices"]==route["global_indices"] and type(proof["slot"]) is int
               and proof["slot"]==i and proof["kind"]==k and leaf["coefficient_hulls_checked"]>0,
               "independent selector bound does not match actual routed leaf")
        bounds=leaf["point_and_box_selected_orders"].values()
        demand(len(bounds)==2 and all(type(v) is int and 2<=v<=(32 if k=="bulk" else 64) for v in bounds),
               "selected order outside declared actual policy")
    return {"case_id":case,"target":dimensions(entry),"plan_sha256":outer["plan_sha256"],
            "terminal_sha256":outer["terminal_sha256"],"observed_os_exit_sha256":root_hash,
            "independent_plan_review_file":rel(q/"review.json"),"independent_plan_review_sha256":independent["review_sha256"],
            "independent_plan_replay_observed_os_exit_sha256":own_hash}


def run(worker_source_review_sha,independent_sha):
    gate=read(GP/"independent_prelaunch/review.json",worker_source_review_sha)
    fields={k:gate[k] for k in ("source_sha256","protocol_sha256","preflight_sha256")}
    read(GS,fields["source_sha256"],False)
    read(GP/"worker_frozen.py",fields["source_sha256"],False)
    p=read(GP/"protocol.json",fields["protocol_sha256"])
    pf=read(GP/"preflight.json",fields["preflight_sha256"])
    required={**EXPECTED,**fields}
    demand(gate["status"]=="PASS_INDEPENDENT_REMAINING32_FIRST_WAVE_SPATIAL_WORKER_COLD"
           and all(gate[k]==v for k,v in required.items()) and yes(gate["checks"])
           and gate["permitted_first_wave_case_ids"]==list(FIRST)
           and type(gate["maximum_owned_workers"]) is int and gate["maximum_owned_workers"]==4
           and gate["no_later_waves_authorized"] is True
           and p["status"]=="FROZEN_RAW_EXACT_SPATIAL_FIRST_WAVE_WORKER" and p["first_wave"]==list(FIRST)
           and p["later_waves_authorized"] is p["root_asserted"] is False
           and p["source_sha256"]==pf["source_sha256"]==fields["source_sha256"]
           and pf["protocol_sha256"]==fields["protocol_sha256"] and pf["status"]=="PASS" and yes(pf["checks"])
           and all(zero(pf[k]) for k in ("new_spatial_integrations","new_derivative_jobs","new_solves")),
           "current frozen worker prelaunch incomplete")
    docs={loc:read(loc,h,loc.suffix==".json") for loc,h in BASIS.items()}
    j=docs[HERE/"bridge_input_admission.json"]
    m=docs[BR/"preparation/manifest.json"]
    root=docs[IR/"root_source_review.json"]
    demand(j["status"]=="PASS_INDEPENDENT_REMAINING32_SAME_TARGET_BRIDGE_INPUTS" and yes(j["checks"])
           and j["spatial_workers_authorized"] is False and len(j["cases"])==32
           and root["status"]=="PASS_ROOT_FIRST_WAVE_SELECTOR_REPLAY_SOURCE"
           and yes(root["checks"]) and all(root[k]==v for k,v in SELECTOR_PHASE.items()),
           "same independent joined/context prerequisites and prior cold selector source review required")
    plans=[]
    for case,raw,own in zip(FIRST,RAW,independent_sha,strict=True):
        entry=one(m["cases"],"case_id",case)
        joined=one(j["cases"],"case_id",case)
        plans.append(actual_selector(case,raw,own,entry,joined))
    targets={v["case_id"]:v["target"] for v in plans}
    demand(targets==p["targets"],"actual admitted trial is not the frozen worker payload")
    out={"status":"PASS_INDEPENDENT_REMAINING32_FIRST_WAVE_SPATIAL_ADMISSION",**required,
         "worker_source_review_sha256":worker_source_review_sha,"joined_bridge_input_admission_sha256":JOIN,
         "plans":plans,"targets":targets,"checks":{"current_cold_worker_scope_first_four_only":True,
              "both_true_observed_selector_processes_per_case":True,
              "all672_sealed_typed_point_and_whole_selector_leaf_plans":True,
              "joined_case_identity_signed_scope_and_independent_geometry_derivative_current":True},
         "no_later_waves_authorized":True,"independent_packet_replay_required":True,
         "scope":"Admit only the first four separately planned raw component runs. A later independent packet review is required before proof reduction.",
         "all_current_sha256":HASHES,"admission_source_sha256":sha(Path(__file__).read_bytes()),
         "new_spatial_jobs_in_this_admission":0}
    dst=HERE/"first_wave_worker_admission.json"
    with dst.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(out,f,indent=2)
        f.write("\n")
    print(json.dumps({"admission_sha256":sha(dst.read_bytes()),"rehash_count":len(HASHES)}))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--worker-source-review-sha256",required=True)
    p.add_argument("--independent-selector-os-sha256",action="append",required=True)
    a=p.parse_args()
    demand(len(a.independent_selector_os_sha256)==4,"need all four orders")
    run(a.worker_source_review_sha256,a.independent_selector_os_sha256)

