"""Independent exact observer readback for each of the seven later waves."""
import argparse
from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
APP = HERE.parent
LATER = APP.parent
BR = LATER.parent
KERNEL = BR/"first_wave_observer_v1/independent_readback/readback.py"
KERNEL_SHA = "9781008eb4c384ae10dfc205fdb977d64d6964168c7ffb1e6e7731d7cf5524e3"
spec = spec_from_file_location("prior_exact_pure_archived_geometry", KERNEL)
H = module_from_spec(spec)
spec.loader.exec_module(H)
ROOT = H.ROOT
q, pair, doc, rel, sha, must, save, pick = H.q,H.pair,H.doc,H.rel,H.sha,H.must,H.save,H.pick
R, L, VAR, OBS = H.R,H.L,H.VAR,H.OBS
WAVES = tuple(f"wave_{i:02}" for i in range(2,9))
REVIEW_SHA = "a3e26e0bfcb02afbbbfbabfd14476cd830d081484302074ab8c832d0a38b150b"
PHASE = {
    "source_sha256":"df0ae901ad4c5bc9c976002c4c62e11e46a5a41b33c747779dfca4fa22871e6b",
    "protocol_sha256":"d20799dc6c69d9f8ab4f36e4d9bd085b2f48b5abf68d821296c3a9ddca74231e",
    "preflight_sha256":"e4978572f081a662adfbd77789171fe66da40365df47c66054f7d90dd45ff357",
    "later_waves_release_sha256":"002b471dceed8666aa0f70e1e9673744f49d0bf3a70da65a194c3faf34b0adbe",
    "profiles_sha256":"31c680040387762f835a235e678435400ac163081509abf9b012494d20de04f0",
    "global_source_sha256":"3051ac6d35ee18407e29e8b82e96d789e685073dd3f067b01b50b459c0b970e5",
    "global_source_review_sha256":"f39bb78a6ee4ee66cfdf25dc287a4856a1f7b5f94298c612fbbc8c17794d171f",
}
ROOT_PASS = "PASS_INDEPENDENT_LATER_WAVE_EXACT_SPATIAL_ROOT"
EXCLUDE = "SUFFICIENT_INDEPENDENT_EXCLUSION_OF_ONLY_THIS_CUBE"
UNRESOLVED = "UNRESOLVED_INDEPENDENT_LATER_WAVE_EXACT_SPATIAL_TEST"
PHASE_KEYS = ("source_sha256","protocol_sha256","preflight_sha256","source_review_sha256",
              "later_waves_release_sha256","profiles_sha256","hull_admission_sha256",
              "raw_wave_admission_sha256","raw_root_observed_exit_sha256")


def hash_only(path,h,links):
    must(sha(path)==h and links.setdefault(rel(path),h)==h,"immutable bytes "+str(path))


def same_shape(value, expected):
    return all(value.get(k)==v for k,v in expected.items())


def phase_zero(os_row, wave, stamp=None):
    return (os_row["wave_id"]==wave and type(os_row["pid"]) is int and os_row["pid"]>0
            and type(os_row["observed_exit_code"]) is int and os_row["observed_exit_code"]==0
            and (type(os_row.get("unified_exec_session")) is int and os_row["unified_exec_session"]>0
                 or type(os_row.get("unified_exec_chunk_id")) is str and bool(os_row["unified_exec_chunk_id"].strip()))
            and (stamp is None or same_shape(os_row,stamp)))


def closure(dir, o, stamp, links, result_status):
    m=[doc(dir/(k+".json"),o[k+"_sha256"],links) for k in ("started","terminal","results")]
    s,t,r=m
    must(all(type(d["pid"]) is int and d["pid"]==o["pid"] and same_shape(d,stamp) for d in m)
         and s["status"]=="RUNNING" and t["status"]=="FINISHED"
         and type(t["exit_code"]) is int and t["exit_code"]==0
         and t["results_sha256"]==o["results_sha256"] and r["status"]==result_status,
         "same observed actual terminal, source and phase")
    return r


def statics(links):
    files={KERNEL:KERNEL_SHA,APP/"apply.py":PHASE["source_sha256"],
           APP/"protocol.json":PHASE["protocol_sha256"],APP/"preflight.json":PHASE["preflight_sha256"],
           APP/"independent_source_review/receipt.json":REVIEW_SHA,
           LATER/"profiles.json":PHASE["profiles_sha256"],
           LATER/"admission/later_waves_release.json":PHASE["later_waves_release_sha256"],
           LATER/"global_reduction_v1/reduce.py":PHASE["global_source_sha256"],
           LATER/"global_reduction_v1/independent_source_review/receipt.json":PHASE["global_source_review_sha256"]}
    for p,h in files.items():
        hash_only(p,h,links)
    p=doc(APP/"protocol.json",PHASE["protocol_sha256"],links)
    pf=doc(APP/"preflight.json",PHASE["preflight_sha256"],links)
    review=doc(APP/"independent_source_review/receipt.json",REVIEW_SHA,links)
    profiles=doc(LATER/"profiles.json",PHASE["profiles_sha256"],links)
    global_phase=p["global_phase"]
    must(same_shape(review,PHASE) and review["approved_waves"]==list(WAVES) and review["checks"]
         and review["status"]=="PASS_INDEPENDENT_LATER_WAVE_OBSERVER_APPLICATION_SOURCE"
         and all(x is True for x in review["checks"].values()) and pf["status"]=="PASS"
         and global_phase["source_sha256"]==PHASE["global_source_sha256"]
         and global_phase["source_review_sha256"]==PHASE["global_source_review_sha256"],"reviewed observer/source phase")
    for n in ("protocol","preflight"):
        hash_only(LATER/"global_reduction_v1"/(n+".json"),global_phase[n+"_sha256"],links)
    must(tuple(z["wave_id"] for z in profiles["waves"])==WAVES and p["ordered_cases"]==
         {row["wave_id"]:row["case_order"] for row in profiles["waves"]},"exact later profiles")
    return p,{z["wave_id"]:z for z in profiles["waves"]}


def read_observer(wave, folders, outer, native_outer, links):
    o=doc(folders/"os_exit_observation.json",outer,links)
    stamp={**PHASE,"wave_id":wave,"review_sha256":REVIEW_SHA,"native_global_observed_exit_sha256":native_outer}
    must(o["status"]=="PASS_ROOT_OBSERVED_LATER_WAVE_OBSERVER_APPLICATION_EXIT" and phase_zero(o,wave,stamp)
         and o["direct_root_tool_observation"] is True and o["independent_actual_scalar_readback_required"] is True,
         "genuinely observed same later observer")
    return closure(folders,o,stamp,links,"LATER_WAVE_EXACT_SPATIAL_OBSERVERS_AWAIT_INDEPENDENT_READBACK")


def native(wave, root, h, profile, p, links):
    n=root/"independent"
    o=doc(n/"root_observed_exit.json",h,links)
    order=profile["case_order"]
    must(o["status"]=="PASS_ROOT_OBSERVED_LATER_WAVE_INDEPENDENT_GLOBAL_EXIT"
         and phase_zero(o,wave) and o["mode"]=="independent","observed later independent global")
    expected={rel(n/(k+".json")) for k in ("started","terminal","results")}
    expected|={rel(n/c/(k+".json")) for c in order for k in ("operands","results")}
    must(len(expected)==11 and set(o["sha256"])==expected,"exact eleven native externally pinned files")
    for path,h in o["sha256"].items():
        hash_only(ROOT/path,h,links)
    for k in ("started","terminal","results"):
        must(o[k+"_sha256"]==o["sha256"][rel(n/(k+".json"))],"same metadata before native reads")
    idx=closure(n,o,{**p["global_phase"],"mode":"independent","wave_id":wave},links,
                "COMPLETE_INDEPENDENT_LATER_WAVE_SIGNED_ENDPOINT_REPLAY")
    must([r["case_id"] for r in idx["cases"]]==order and idx["signed_hull_components"]==672,
         "exact four global closures")
    s=doc(n/"started.json",o["started_sha256"],links)
    t=doc(n/"terminal.json",o["terminal_sha256"],links)
    stamp={k:idx[k] for k in (*PHASE_KEYS,"producer_observed_os_sha256")}
    must(all(same_shape(v,stamp) for v in (s,t)),"same independent native input phase")
    pr=root/"producer"
    po=doc(pr/"root_observed_exit.json",idx["producer_observed_os_sha256"],links)
    must(po["status"]=="PASS_ROOT_OBSERVED_LATER_WAVE_GLOBAL_PRODUCER_EXIT"
         and po["mode"]=="producer" and phase_zero(po,wave),"separately observed producer native phase")
    pi=closure(pr,po,{k:idx[k] for k in PHASE_KEYS}|{"mode":"producer","wave_id":wave},
               links,"COMPLETE_LATER_WAVE_GLOBAL_PRODUCER_AWAITS_INDEPENDENT_ENDPOINT_REPLAY")
    must([r["case_id"] for r in pi["cases"]]==order,"separate producer exact ordered cases")
    for prr in pi["cases"]:
        folder=pr/prr["case_id"]
        for key,k in (("case_results_sha256","results"),("operands_sha256","operands")):
            hash_only(folder/(k+".json"),prr[key],links)
    outputs=[]
    for entry in idx["cases"]:
        c=entry["case_id"]
        require_index={a:entry[b] for a,b in (("operands","operands_sha256"),("results","case_results_sha256"))}
        must(all(o["sha256"][rel(n/c/(k+".json"))]==h for k,h in require_index.items()),"matching exact index pins")
        leaf=doc(n/c/"results.json",entry["case_results_sha256"],links)
        must(same_shape(leaf,stamp) and leaf["wave_id"]==wave and leaf["mode"]=="independent"
             and type(leaf["pid"]) is int and leaf["pid"]==o["pid"]
             and leaf["case_id"]==c and leaf["status"]==entry["status"] and leaf["archived_observer_field_replaced"] is False
             and leaf["operands_sha256"]==entry["operands_sha256"]
             and same_shape(profile["targets"][c],leaf["target"]),"same physical target and phase")
        outputs.append((entry,leaf))
    return outputs



def read_case(entry,leaf,sources,reported,links):
    name=entry["case_id"]
    endpoint,number=name.rsplit("_small_",1)
    checkpoint="small_"+number
    table,audit,old,old_a,secondary,secondary_a,primary,geom,allow=sources
    tgt=leaf["target"]
    trial=pick(geom["trials"],case_id=name)
    ctr=doc(trial["model_contract"],trial["model_contract_sha256"],links)
    bdoc=pick(table["field_budgets"],endpoint_id=endpoint,checkpoint=checkpoint)
    bd=pick(audit["fields"],endpoint_id=endpoint,checkpoint=checkpoint)
    prior=pick(old["cases"],case_id=name)
    file=rel(H.SEC/"actual"/endpoint/(checkpoint+".json"))
    sp=pick(secondary["field_observer_receipts"],path=file)["sha256"]
    saved=doc(file,sp,links)
    must(leaf["case_id"]==tgt["case_id"]==name and tgt["record_id"]==trial["record_id"]==ctr["record_id"]
         and all(tgt[k]==trial[k] for k in ("trial_box_sha256","model_contract_sha256","point_cache_sha256",
                          "R_binary_sha256","center_binary_sha256","archived_center_shift_exact","exact_drive","energy_scale_binary_exact"))
         and q(leaf["radius_exact"])==q(trial["trial_radius_exact"])==R,
         "same later root, proof center, material and R")
    box=doc(trial["trial_box"],trial["trial_box_sha256"],links)
    must(box["record_id"]==tgt["record_id"] and box["point_cache_sha256"]==tgt["point_cache_sha256"]
         and box["model_observation_contract_sha256"]==tgt["model_contract_sha256"] and q(box["corrector_radius_exact"])==R
         and bdoc["record_id"]==prior["record_id"]==bd["record_id"]==ctr["physical_source_record_id"]
         and saved["step"]==bdoc["step"] and saved["role"]==bdoc["role"]=="small_record"
         and secondary_a["input_and_output_sha256"][file]==sp
         and saved["endpoint_id"]==endpoint and saved["checkpoint"]==checkpoint
         and saved["field_sha256"]==prior["field_sha256"]==bd["field_sha256"]==ctr["archived_field_sha256"]
                                  ==tgt["original_observer_field_sha256"]
         and saved["field"]==bd["field"]==ctr["archived_field"]
         and saved["verification"]==ctr["archived_verification"] and saved["verification_sha256"]==ctr["archived_verification_sha256"]
         and all(bdoc[k]==saved[k] for k in ("field","field_sha256","verification","verification_sha256","record_id"))
         and q(tgt["observer_goal_exact"])==q(bd["e_goal_dimensionless_exact"])==q(bdoc["e_goal_dimensionless_exact"])
               ==q(ctr["observer_goal_exact"])
         and q(tgt["original_addition_gap_exact"])==q(bd["binary_gap_exact"])
               ==q(bdoc["center_and_model"]["node_addition_gap_dimensionless_exact"])==q(ctr["archived_addition_gap_exact"])
         and q(tgt["archived_center_shift_exact"])==q(ctr["archived_shift_exact"])==q(prior["root"]["sigma"])
         and ctr["archived_observer_variant"]==VAR[0] and q(ctr["archived_observer_allowance_um_exact"])==allow,
         "same source-linked original archive, verification, original observer budget and recenter provenance")
    requested=ctr["exact_requested_drive"]
    must(q(tgt["exact_drive"])==q(requested["activation_strength_exact"])==q(bdoc["center_and_model"]["exact_target_strength"])
         and q(requested["common_gain_exact"])==q(bdoc["center_and_model"]["gain_exact"])
         and q(requested["optical_dose_exact"])==q(box["dose_exact"])==q(saved["step"]["dose_exact"])
         and q(requested["activation_strength_exact"])==q(requested["common_gain_exact"])*q(requested["optical_dose_exact"]),
         "exact gain-dose product, original fixed input")
    for k in ("field","verification"):
        path=saved[k]
        expected=saved[k+"_sha256"]
        must(primary["all_hashes"][path]==old_a["current_input_sha256"][path]
            ==secondary_a["input_and_output_sha256"][path]==expected,"all three independent recorded archived input digests agree")
    field=doc(saved["field"],saved["field_sha256"],links)
    verify=doc(saved["verification"],saved["verification_sha256"],links)

    e,sigma,xi,eps=root_bounds(entry,leaf)
    roi,archived=H.archived_geometry(field,box,bdoc,prior,saved,q(0),sigma,xi)
    block="ROOT_NOT_CERTIFIED_IN_THIS_CUBE" if e is None else None
    graph=None
    if eps is not None:
        if any(q(v["stored_radial_min_um_exact"])<=8*eps for v in archived.values()):
            block="STRICT_RADIAL_GRAPH_GUARD_NOT_ESTABLISHED"
        else:
            graph={}
            for k,v in archived.items():
                m,s=(q(v[n]) for n in ("stored_radial_min_um_exact","slope_cap_exact"))
                graph[k]={**v, "strict_radial_margin_um_exact":str(m-8*eps),
                          "graph_error_um_exact":str(q("5/4")*(1+s)*eps)}
    same={"case_id":name,"root_proof_sha256":entry["case_results_sha256"],"record_id":bdoc["record_id"],
          "native_root_status":leaf["status"],"original_field_sha256":saved["field_sha256"],
          "posterior_exact":None if e is None else str(e),
          "archived_sigma_exact":str(sigma),"original_xi_exact":str(xi),
          "physical_nodal_error_um_exact":None if eps is None else str(eps),
          "ROI_um_exact":list(map(str,roi)),"graph_premises":graph}
    must(all(reported[k]==value for k,value in same.items()) and len(reported["comparisons"])==16
         and verify["variant"]==VAR[0] and set(saved["variant_observables_um_exact"])==set(VAR),
         "observer metadata exactly reproduces independent root and archive geometry")
    comps=[]
    for var in VAR:
        stored=saved["variant_observables_um_exact"][var]
        must(set(stored)==set(OBS),"four frozen archived observables")
        for i,obs in enumerate(OBS):
            actual=pick(reported["comparisons"],variant=var,observable=obs)
            sample=q(saved["step"]["observed_exact"][obs])
            if var==VAR[0]:
                prev=pick(prior["comparisons"],observable=obs)
                must(pair(stored[obs])==pair(verify["observer_intervals_um_exact"][obs])==pair(prev["stored_interval_um_exact"])
                     and sample==q(prev["observed_um_exact"]),"stored primary comparator identity")
            if block:
                error=roots=res=None
                standing=block
            else:
                ea,eb=(q(graph[k]["graph_error_um_exact"]) for k in ("apical","basal"))
                error=(ea+eb,ea+eb,ea,eb)[i]
                roots,res,standing=H.interval_check(stored[obs],sample,error,allow)
            cell={"variant":var,"observable":obs,"stored_interval_um_exact":list(map(str,pair(stored[obs]))),
                  "observed_um_exact":str(sample),"extra_observer_error_um_exact":None if error is None else str(error),
                  "exact_spatial_root_observer_interval_um_exact":None if roots is None else list(map(str,roots)),
                  "absolute_residual_um_exact":None if res is None else list(map(str,res)),"standing":standing}
            must(cell==actual,"exact observer interval, residual or standing differs from producer")
            comps.append(cell)
    return {**same,"comparisons":comps}


def scalars(rows, hrows, zq, radius):
    must(rows and len(rows)==len(hrows) and min(hrows)>=0 and zq>=0 and radius>0,"strict scalar domain")
    lo=max(max(0,l,-u) for l,u in rows)
    eta=max(max(abs(l),abs(u)) for l,u in rows)
    bh=max(hrows)
    z=zq+bh
    lhs=eta+z*radius
    outside=lo>(1+z)*radius
    inside=z<1 and lhs<radius
    must(not (inside and outside),"disjoint root decisions")
    e=eta/(1-z) if inside else None
    return lo,eta,bh,z,lhs,inside,outside,e


def root_bounds(entry,leaf):
    a,t,v=(leaf[k] for k in ("arithmetic","target","observer_and_optional_minimum"))
    rf=list(map(pair,a["signed_R_F_I_rows_exact"]))
    rh=list(map(q,a["RH_absolute_row_bounds_exact"]))
    must(len(rf)==len(rh)==611,"all signed rows")
    lo,eta,bh,z,lhs,inside,outside,e=scalars(rf,rh,q(a["ZQ_upper_exact"]),q(leaf["radius_exact"]))
    status=ROOT_PASS if inside else EXCLUDE if outside else UNRESOLVED
    must(all(q(a[k])==v for k,v in {
              "RF_lower_exact":lo,"eta_upper_exact":eta,"b_H_upper_exact":bh,"z_upper_exact":z,
              "selfmap_lhs_upper_exact":lhs,"selfmap_margin_lower_exact":R-lhs,
              "exclusion_margin_lower_exact":lo-(1+z)*R}.items())
         and a["strict_selfmap"] is inside and a["strict_cube_exclusion"] is outside
         and entry["status"]==leaf["status"]==status
         and entry["exact_spatial_root_certified"] is inside and leaf["root_asserted"] is inside,
         "independently rebuilt eta/z/selfmap/exclusion from signed rows")
    sigma,xi,goal=(q(t[k]) for k in ("archived_center_shift_exact","original_addition_gap_exact","observer_goal_exact"))
    must(0<=sigma<goal and xi>=0,"original shift bound")
    eps=None if e is None else L*(e+sigma+xi)
    budget=e is not None and sigma+e<goal
    expected={"archived_shift_plus_posterior":None if e is None else str(e+sigma),
              "physical_nodal_error_um_upper_exact":None if eps is None else str(eps)}
    must(a["posterior_error_upper_exact"]==(None if e is None else str(e))
         and same_shape(v,expected) and v["archived_state_budget_proved"] is budget
         and entry["shift_inclusive_observer_budget_proved"] is budget
         and (e is None or 0<=e<R),"exact e, sigma and xi including unavailable roots")
    return e,sigma,xi,eps


def fixture():
    checks=H.analytical_fixtures()
    must(scalars([(q("1/4"),q("1/4"))],[q("1/4")],q(0),q(1))[-3:]==(True,False,q("1/3")),
         "synthetic strict root")
    must(scalars([(-q(2),q(1))],[q(0)],q(0),q(1))[-3:]==(False,False,None),"zero-crossing unresolved")
    must(scalars([(q(2),q(3))],[q(0)],q(0),q(1))[-3:]==(False,True,None),"strict excluded cube")
    must(scalars([(q(1),q(1))],[q(0)],q(0),q(1))[-3:]==(False,False,None),"no strict root or exclusion at equality")
    return {**checks,"signed_global_rows_all_three_cases_and_equality":True}


def prepare():
    links={}
    p,profiles=statics(links)
    H.original_sources(links)
    protocol={"status":"FROZEN_LATER_WAVE_INDEPENDENT_ARCHIVE_SCALAR_READBACK",
              "source_sha256":sha(Path(__file__)), "pins":links, "profile_order":list(profiles),
              "later_wave_cases":p["ordered_cases"],"checks":fixture(),"opened_actual_later_numeric_payloads":0,
              "requires_external_hash_of_both_observed_phases":True,"producer_evaluator_imported":False}
    print(save(HERE/"protocol.json",protocol),flush=True)


def run(args):
    links={}
    protocol=doc(HERE/"protocol.json",args.protocol_sha256,links)
    p,profiles=statics(links)
    sources=H.original_sources(links)
    must(protocol["source_sha256"]==sha(Path(__file__)) and protocol["checks"]==fixture()
         and protocol["later_wave_cases"]==p["ordered_cases"]
         and protocol["pins"]=={k:v for k,v in links.items() if k!=rel(HERE/"protocol.json")},
         "frozen independent cold phase")
    wave=args.wave
    folder=LATER/"waves"/wave
    profile=profiles[wave]
    rows=native(wave,folder/"global_reduction_v1",args.native_os_sha256,profile,p,links)
    actual=read_observer(wave,folder/"observer_application_v1/actual",
                         args.producer_os_sha256,args.native_os_sha256,links)
    must([x["case_id"] for x in actual["cases"]]==profile["case_order"]
         and q(actual["allowance_um_exact"])==sources[-1] and type(actual["new_spatial_jobs"]) is int
         and actual["new_spatial_jobs"]==0 and actual["full_calibration_or_source_fit_claimed"] is False,
         "same four observer cases and exact common archived allowance")
    out=folder/"observer_application_v1/independent_actual_readback"
    out.mkdir(exist_ok=False)
    stamp={"source_sha256":sha(Path(__file__)),"protocol_sha256":args.protocol_sha256,"pid":os.getpid(),"wave_id":wave,
           "producer_os_sha256":args.producer_os_sha256,"native_global_os_sha256":args.native_os_sha256}
    save(out/"started.json",{"status":"RUNNING",**stamp})
    terminal={"status":"FAILED","exit_code":1,**stamp}
    try:
        cases=[read_case(e,l,sources,r,links) for (e,l),r in zip(rows,actual["cases"],strict=True)]
        counts={v:dict(Counter(c["standing"] for a in cases for c in a["comparisons"] if c["variant"]==v)) for v in VAR}
        must(counts==actual["counts"] and all(sum(v.values())==16 for v in counts.values())
             and actual["total_primary_comparisons"]==16 and actual["total_alternative_comparisons"]==48,
             "exact 64 retained comparison rows")
        outs=[{"case_id":a["case_id"],"variant":c["variant"],"observable":c["observable"],
               "strict_OUT_margin_um_exact":str(pair(c["absolute_residual_um_exact"])[0]-sources[-1]),
               "absolute_residual_um_exact":c["absolute_residual_um_exact"]}
              for a in cases for c in a["comparisons"] if c["standing"]=="OUT"]
        must(all(sha(ROOT/k)==v for k,v in links.items()),"source bytes at closure")
        h=save(out/"receipt.json",{"status":"PASS_INDEPENDENT_ACTUAL_LATER_WAVE_SPATIAL_OBSERVER_READBACK",
               **stamp,"observer_producer_results_sha256":links[rel(folder/"observer_application_v1/actual/results.json")],
               "allowance_um_exact":str(sources[-1]),"cases":cases,"counts":counts,
               "strict_out":outs,"all_64_exact_comparisons_recreated":True,
               "source_and_data_sha256":links,"shared_pure_archived_geometry_source_sha256":KERNEL_SHA,
               "producer_apply_imported":False,"new_force_FE_AD_spatial_integral_or_solve_jobs":0,
               "complete_calibration_or_source_fit_claimed":False})
        terminal.update(status="FINISHED",exit_code=0,receipt_sha256=h)
        print({"receipt_sha256":h,"counts":counts},flush=True)
    finally:
        save(out/"terminal.json",terminal)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    sub=ap.add_subparsers(dest="action",required=True)
    sub.add_parser("prepare")
    a=sub.add_parser("run")
    a.add_argument("--wave",required=True,choices=WAVES)
    for key in ("protocol","producer-os","native-os"):
        a.add_argument("--"+key+"-sha256",required=True)
    n=ap.parse_args()
    prepare() if n.action=="prepare" else run(n)

