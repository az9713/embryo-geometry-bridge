"""Explicit remaining32 exact-component cold preparation and guarded kernels."""
import argparse
from fractions import Fraction
import importlib.util
import os
from pathlib import Path
import traceback

import transfer_ambiguity_exact_bridge as B
import transfer_ambiguity_remaining32_geometry_variation as V
from transfer_ambiguity_remaining32_exact_targets import canonical, rel

SOURCE=Path(__file__).resolve()
ROOT=V.ROOT
OUT=V.A.X.OUT/"remaining32_bridge_v1"
PREP=OUT/"preparation"
PLAN=OUT/"plans"
JOINED=OUT/"admission/bridge_input_admission.json"
COLD_REVIEW_PATH=OUT/"independent_prelaunch/review.json"
RADIUS=Fraction(1,2**34)
E=V.E
POINT,BOX=B.POINT,B.BOX
N=611
Q=V.Q
read,save,require,sha=V.read,V.save,V.require,V.sha
OPS=tuple((kind,j) for kind,n in (("bulk",144),("surface",24)) for j in range(n))
PINS={
 Path(B.__file__).resolve():"58ae92f849414f185c4e28bbb2ec3b9d05361391c4d6e3fa71f519162e30e5e7",
 Path(V.__file__).resolve():"68bca869ee486d7dc25b34cae3ceee1ddc2adc53500a20301c0038c259a2b26a",
 V.GOUT/"preparation/protocol.json":"754da64bb96cd8a19707181315e6c3e3f014e8d37a898b83c1693342c82c9850",
 V.GOUT/"preparation/preflight.json":"d70a7217c4dbbf68f789bb59193321252ea4192576c3f41790cc7daea38fb15a",
 V.DOUT/"preparation/protocol.json":"6d04fc82779430e66c2d27b5ab7aa508f522f7ffac03b87f8f2e8e2587e1d128",
 V.DOUT/"preparation/preflight.json":"abbbeb7c0c303db2d2675398f8ea09155a19fd969b121e7936d7f4a1a468e264",
 V.TARGET_REVIEW:"67d33bbc3171839b1c578b02b5cfbf9ea6ba5d6cb351d0655dd14aaf99df75f9",
}
FORCE=ROOT/"runs/exact_integration_cohort_v1/L_defect_diagnostic/element21_force.py"
SCALED=ROOT/"scripts/exact_scaled_cohort_bulk.py"
SCALED_SERIES=ROOT/"runs/exact_integration_cohort_v1/L_defect_diagnostic/scaled_chart_v1/scaled_element21.py"
PINS[FORCE]="9a5f93bb3d4a39020e28bee165ba606d1d17f63e44a17ee758959aa05109d329"
PINS[SCALED]="6639bb2fae3c3892ec44115ff90e27242586744a2900daa16177cc5a091f1526"
PINS[SCALED_SERIES]="a173b3759fb66fe66e5f3e442c1bcce9a5a185e67515649a94499f82d1d67f60"
GEOM_REC=V.GOUT/"independent_actual_preparation/actual_replay"
DER_REC=V.DOUT/"independent_actual_preparation/actual_replay"
NATIVE_CODE={
 "variation":("replay.py",
      "e33c067b49c428d36babb2e75858e4b18fde42164af5dfc2f9244d500d368028",
      "3baf882ddab602fe19437a2af2e78cafd6df4b0b0182d883a72847ea3b0e855d",
      "12a44c45e576a5e3dcc2e04b3a34276bc3e342a5173e50aa02ceb09a00abd2cd"),
}
code,sh,ph,fh=NATIVE_CODE["variation"]
folder=DER_REC.parent
PINS.update({folder/code:sh,folder/"protocol.json":ph,folder/"preflight.json":fh})
PINS[ROOT/"scripts/transfer_ambiguity_cached_derivative_extension.py"]="7ad0b7d4d42efe634123989c6092ccf768cedb92a79244cc7aaa5c10a62f370c"
REVIEW="PASS_INDEPENDENT_REMAINING32_EXACT_BRIDGE_COLD"
ADMISSION="PASS_INDEPENDENT_REMAINING32_SAME_TARGET_BRIDGE_INPUTS"
PLAN_REVIEW="PASS_INDEPENDENT_REMAINING32_ACTUAL_SELECTOR_PLAN"
SCOPE="Fixed-target signed exact-spatial FE point force minus Qcache, and H defect over whole same-center r34 box, including finite same-Y corrections."
ZERO={"selector_components_prepared":0,"new_spatial_integrations":0,"cached_AD_jobs":0,"solves":0}


def kernels():
    for p in (FORCE,SCALED,SCALED_SERIES):
        require(sha(p)==PINS[p],"verified explicit reusable spatial kernel")
    spec=importlib.util.spec_from_file_location("_remaining32_force_parts",FORCE)
    force=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(force)
    import exact_scaled_cohort_bulk as scaled
    require(Path(scaled.__file__).resolve()==SCALED and Path(scaled.S.__file__).resolve()==SCALED_SERIES
            and force.E is scaled.E is E and force.T is scaled.T is E.T
            and force.core is scaled.core is E.core and B.RADIUS==RADIUS,
            "reuse of pure same-radius kernels only")
    return force,scaled


def routes_order(rows,n=N):
    expected=OPS if n==N else tuple((x["kind"],x["slot"]) for x in rows)
    require(tuple((x["kind"],x["slot"]) for x in rows)==expected
        and len(set(expected))==len(rows)>0 and all(type(x["slot"]) is int for x in rows),
        "typed ordered unique complete component catalogue")
    require(all(isinstance(x["global_indices"],list) and len(x["global_indices"])>0
            and len(set(x["global_indices"]))==len(x["global_indices"])
            and all(type(i) is int and 0<=i<n for i in x["global_indices"]) for x in rows)
            and {j for x in rows for j in x["global_indices"]}==set(range(n)),
            "authentic full free projection required")
    return expected


def sizes(catalogue,z,rnorm,n=N):
    routes_order(catalogue,n)
    z,rnorm=Q(z),Q(rnorm)
    require(0<=z and rnorm>0,"actual nonnegative ZQ and exact R infinity norm required")
    f=[0]*n
    H=[0]*n
    for entry in catalogue:
        ids=entry["global_indices"]
        for i in ids:
            f[i]+=1
            H[i]+=len(ids)
    h=1-z
    base={"ZQ_exact":str(z),"R_inf_exact":str(rnorm),"positive_headroom":h>0,
          "headroom_exact":str(h),"c_F":max(f),"c_H":max(H),"strictly_tail_only":True,
          "quarter_ZQ_admission_rule":False}
    if h<=0:
        return {**base,"caps_exact":None}
    cap_H=h/2
    cap_F=(h-cap_H)*RADIUS/2
    return {**base,"conditional_total_force_target_exact":str(cap_F),
            "conditional_H_defect_target_exact":str(cap_H),
            "caps_exact":list(map(str,(cap_F/(10000*rnorm*max(f)),cap_H/(10000*rnorm*max(H)))))}


def catalogue(routes):
    out=[]
    for row in routes["bulk"]:
        out.append({"kind":"bulk","slot":row["element"],"nodes":row["nodes"],
                    "global_indices":row["global_indices"],"chart":row["chart"],
                    "normalized_rule_sha256":row["auxiliary_element_rule_sha256"]})
    for row in routes["surface"]:
        out.append({"kind":"surface","slot":12*row["face"]+row["edge"],"nodes":row["nodes"],
                    "global_indices":row["global_indices"],"chart":"surface_physical_R",
                    "normalized_rule_sha256":row["trace_rule_sha256"]})
    routes_order(out)
    return out


def premises():
    m,phase,_=V.bindings()
    for p,h in PINS.items():
        require(sha(p)==h,"frozen neutral source or stage prep changed")
    V.target_gate(V.TARGET_REVIEW,PINS[V.TARGET_REVIEW],m)
    target,w,mp,rec,meta_links=V.A.load_authority()
    require(tuple(target)==V.A.IDS and tuple(phase["geometry"][0]["case_order"])==V.A.IDS,
            "current independent context authority and trial population")
    new={t["case_id"]:t for t in phase["variation"][0]["trials"]}
    out=[]
    for trial in phase["geometry"][0]["trials"]:
        case=trial["case_id"]
        c,row,_,_,_,_=V.A.case_authority(case,(target,w,mp,rec,meta_links))
        contract=read(canonical(trial["model_contract"]),trial["model_contract_sha256"])
        require(Q(trial["trial_radius_exact"])==RADIUS and (case in new)==(not trial["native_cached_bound_usable_at_trial_radius"])
                and (case not in new or trial==new[case]),"same trial for all geometry and variation")
        route=read(canonical(trial["routes"]),trial["routes_sha256"])
        bounds=read(canonical(trial["trial_box"]),trial["trial_box_sha256"])
        require(Q(bounds["corrector_radius_exact"])==RADIUS and bounds["point_cache_sha256"]==trial["point_cache_sha256"]
                and bounds["record_id"]==trial["record_id"] and route["identity"]["record_id"]==trial["record_id"]
                and Q(trial["exact_drive"])==Q(row["requested_drive_exact"])
                and row["binary_R_sha256"]==trial["R_binary_sha256"]==contract["binary_R_sha256"]
                and row["unsigned_binary_center_sha256"]==trial["center_binary_sha256"]==contract["center_binary_sha256"]
                and Q(trial["archived_center_shift_exact"])==Q(contract["archived_shift_exact"])
                    ==Q(row["archived_field_center_shift_sigma_exact"]),
                "new authority cannot borrow old pilot context or trial box")
        out.append({"case_id":case,"record_id":trial["record_id"],"trial":trial,
            "source_record_id":c["record_id"],"cached_input_kind":"FRESH_R34" if case in new else "NATIVE_R34",
            "native_case":row["native_case"],"native_case_sha256":row["native_case_sha256"],
            "native_cached_phase_sha256":row["native_phase_sha256"],"native_radius_exact":row["cached_proof_radius_exact"],
            "native_ZQ_exact":row["cached_ZQ_upper_exact"],"observer_goal_exact":row["observer_goal_exact"],
            "original_observer_field":row["field"],"original_observer_field_sha256":row["field_sha256"],
            "original_verification_sha256":row["verification_sha256"],"original_addition_gap_exact":row["node_addition_gap_exact"],
            "original_allowance_um_exact":c["frozen_conditional_allowance_um_exact"],"observer_variant":c["observer_variant"],
            "component_catalogue":catalogue(route),"point_force_scope":POINT,"whole_H_scope":BOX})
    require(len(out)==32 and sum(x["cached_input_kind"]=="FRESH_R34" for x in out)==10,"full 22+10 identity split")
    return out,meta_links

RAW={"geometry":(V.GOUT/"actual/root_terminal_admission.json",
                  "50ce32f50b5e330b9c78fc1863527c49eb6c5e4394bdb1c4e4313cea8628efe0"),
     "variation":(V.DOUT/"actual/root_terminal_admission.json",
                  "91b1db64f27c89f50643fdde4f16c70a0dcac9ee51566fcf5f352b5b1f663f18")}
NATIVE_REVIEW={
 "geometry":(GEOM_REC,"PASS_ROOT_OBSERVED_REMAINING32_GEOMETRY_REPLAY_EXIT",
             "PASS_INDEPENDENT_ACTUAL_REMAINING32_SINGLETON_GEOMETRY","geometry_records","pilot_id"),
 "variation":(DER_REC,"PASS_OBSERVED_INDEPENDENT_REMAINING32_VARIATION_TRUE_OS_EXIT",
              "PASS_INDEPENDENT_ACTUAL_TEN_FRESH_CACHED_DERIVATIVES","case_results","case_id"),
}


def replay_phase(stage,result):
    if stage=="geometry":
        return result["source_and_phase"]
    keys=("pid","source_sha256","protocol_sha256","preflight_sha256",
          "producer_admission_sha256","root_source_review_sha256")
    return {k:result[k] for k in keys}


def native_stage(stage,outer_sha,expected_order,approved_geometry=None):
    folder,observed_status,result_status,population,key=NATIVE_REVIEW[stage]
    obs=read(folder/"os_exit_observation.json",outer_sha)
    require(obs["status"]==observed_status and V.positive_int(obs["pid"])
            and V.positive_int(obs["unified_exec_session"]) and V.zero(obs["observed_exit_code"]),
            "independent actual replay requires externally observed literal OS success")
    terminal=read(folder/"terminal.json",obs["terminal_sha256"])
    root_raw,raw_sha=RAW[stage]
    require(V.positive_int(terminal["pid"]) and terminal["pid"]==obs["pid"] and V.zero(terminal["exit_code"])
            and terminal["status"]=="FINISHED" and terminal["receipt_sha256"]==obs["receipt_sha256"],
            "same independent actual native replay terminal")
    result=read(folder/"receipt.json",obs["receipt_sha256"])
    phase=replay_phase(stage,result)
    require(result["status"]==result_status and V.positive_int(phase["pid"]) and phase["pid"]==obs["pid"]
            and all(terminal.get(k,terminal.get("source_and_phase",{}).get(k))==phase[k]
                    for k in ("protocol_sha256","preflight_sha256","source_sha256")),
            "admit only the exact independent actual native replay output")
    if stage=="geometry":
        static=approved_geometry
        require(type(static) is dict and set(static)=={"source_sha256","protocol_sha256","preflight_sha256"},
                "joined external gate must identify the approved native geometry source and phase")
        paths={"source_sha256":"replay_geometry.py","protocol_sha256":"protocol.json","preflight_sha256":"preflight.json"}
        require(all(sha(folder.parent/name)==static[k] for k,name in paths.items()),
                "same currently approved geometry verifier on disk")
    else:
        _,src,ph,fh=NATIVE_CODE[stage]
        static={"source_sha256":src,"protocol_sha256":ph,"preflight_sha256":fh}
    require(all(phase[k]==terminal[k]==v for k,v in static.items()),"same frozen independent native verifier")
    raw_key="producer_root_admission_sha256" if stage=="geometry" else "producer_admission_sha256"
    source_key="source_review_sha256" if stage=="geometry" else "root_source_review_sha256"
    require(phase[raw_key]==terminal[raw_key]==raw_sha and terminal[source_key]==phase[source_key],
            "independent actual proof must bind original observed producer and source approval")
    review=read(folder.parent/"root_source_review.json",phase[source_key])
    require(review["status"]==("PASS_ROOT_REMAINING32_GEOMETRY_REPLAY_SOURCE" if stage=="geometry" else "PASS")
            and V.yes(review["checks"]) and all(review[k]==v for k,v in static.items()),
            "reviewed native verifier execution phase required")
    raw=read(root_raw,raw_sha)
    require(raw["stage"]==stage and V.zero(raw["observed_exit_code"]) and V.positive_int(raw["supervisor_pid"])
            and V.positive_int(raw["pid"]) and raw["supervisor_pid"]==raw["pid"]
            and V.positive_int(raw["unified_exec_session"]) and raw["launch_group"]==E._digest(list(expected_order)),
            "same externally observed stage batch, no partial or foreign authority")
    rows=result[population]
    require(tuple(x[key] for x in rows)==expected_order,"full independent stage population required")
    raw_map=raw["sha256"]
    linked=result["checked_input_sha256" if stage=="geometry" else "input_sha256"]
    require(len(raw_map)==(195 if stage=="geometry" else 53) and all(linked.get(k)==v for k,v in raw_map.items()),
            "native reviewer did not authenticate every external producer artifact")
    return {x[key]:x for x in rows},raw_map,phase


def packed_signed(rows):
    require(type(rows) is list and len(rows)==N,"complete signed cached center anchor")
    checked=[]
    for x in rows:
        raw=(x["lower_exact"],x["upper_exact"]) if isinstance(x,dict) else x
        lo,hi=map(Q,raw)
        require(lo<=hi,"ordered independently enclosed signed point RF")
        checked.append([str(lo),str(hi)])
    return checked


def match_derivative(entry,leaf):
    tr=entry["trial"]
    want={"case_id":entry["case_id"],"record_id":entry["record_id"],
          "native_cached_anchor_sha256":entry["native_case_sha256"],"trial_box_sha256":tr["trial_box_sha256"],
          "model_contract_sha256":tr["model_contract_sha256"],"cache_packet_sha256":tr["point_cache_sha256"],
          "R_binary_sha256":tr["R_binary_sha256"],"center_binary_sha256":tr["center_binary_sha256"],
          "proof_drive_exact":tr["exact_drive"],"radius_exact":str(RADIUS),
          "archived_center_shift_exact":tr["archived_center_shift_exact"],"archived_xi_exact":entry["original_addition_gap_exact"]}
    require(leaf["status"]=="PASS_REVIEWED_FRESH_R34_CACHED_DERIVATIVE"
            and all(leaf[k]==v for k,v in want.items()) and leaf["site_counts"]==V.COUNT,
            "fresh replayed same-r34 cached derivative and signed anchor required")
    z=Q(leaf["reviewed_ZQ_upper_exact"])
    require(z==max(Q(leaf["independent_ZQ_upper_exact"]),Q(leaf["producer_ZQ_upper_exact"]))>=0,
            "admitted ZQ must conservatively cover both fresh independent bounds")
    return z,Q(leaf["R_inf_exact"]),packed_signed(leaf["signed_native_cached_center_rows_exact"])


def native_derivative(entry):
    tr=entry["trial"]
    require(entry["cached_input_kind"]=="NATIVE_R34" and Q(entry["native_radius_exact"])==RADIUS
            and tr["native_cached_bound_usable_at_trial_radius"] is True,"same original native r34 only")
    p=read(canonical(entry["native_case"]),entry["native_case_sha256"])
    require(p["cache_packet_sha256"]==tr["point_cache_sha256"] and p["case_id"]==entry["case_id"]
            and p["field_sha256"]==entry["original_observer_field_sha256"] and p["record_id"]==entry["source_record_id"]
            and Q(p.get("exact_strength",p.get("exact_requested_strength")))==Q(tr["exact_drive"]),
            "original same binary and drive native r34 proof")
    if p["status"]=="PASS_INDEPENDENT_CACHED_LOCAL_ROOT":
        at,=[a for a in p["own_radius_attempts"] if Q(a["radius_exact"])==RADIUS]
        require(at["status"]=="PASS_INDEPENDENT_CACHED_LOCAL_ROOT" and at["whole_cube_AD_recomputed_for_this_radius"] is True,
                "matched positive original whole-cube native replay")
        require(p["per_tensor_sha256"]["preconditioner"]==tr["R_binary_sha256"]
                and p["per_tensor_sha256"]["u"]==tr["center_binary_sha256"],"original exact R and center")
        z=Q(at["ZQ_upper_exact"])
        norm=Q(at["preconditioner_infinity_norm_exact"])
        signed=packed_signed(at["signed_center_force_intervals"])
    else:
        require(p["status"]=="PASS_INDEPENDENT_NEGATIVE_M_SMALLER_CACHED_ROOT"
                and Q(p["radius_exact"])==RADIUS and p["preconditioner_binary_sha256"]==tr["R_binary_sha256"],
                "exception only for native negative M narrowed admissible radius")
        witness=read(canonical(p["reused_signed_cached_anchor_file"]),p["reused_signed_cached_anchor_sha256"])
        require(witness["cache_packet_sha256"]==tr["point_cache_sha256"]
                and witness["cache_tensor_content_sha256"]["u"]==tr["center_binary_sha256"]
                and witness["cache_tensor_content_sha256"]["preconditioner"]==tr["R_binary_sha256"],
                "same original negative M signed cached center")
        z=Q(p["independent_ZQ_upper_exact"])
        norm=Q(p["preconditioner_infinity_norm_exact"])
        signed=packed_signed(witness["signed_center_force_intervals"])
    require(0<=z==Q(entry["native_ZQ_exact"])<1,"same independent native derivative bound")
    return z,norm,signed


def match_geometry(entry,leaf,summary):
    tr=entry["trial"]
    expect={k:tr[k] for k in ("record_id","model_contract_sha256","point_cache_sha256","R_binary_sha256",
                              "center_binary_sha256","exact_drive","energy_scale_binary_exact")}
    expect.update(pilot_id=entry["case_id"],radius_exact=str(RADIUS),box_sha256=tr["trial_box_sha256"],
            archived_field_sha256=entry["original_observer_field_sha256"],
            archived_addition_gap_exact=entry["original_addition_gap_exact"],
            archived_shift_exact=tr["archived_center_shift_exact"],
            archived_observer_goal_exact=entry["observer_goal_exact"])
    require(leaf["status"]==summary["status"]=="PASS_INDEPENDENT_REMAINING32_WHOLE_BOX_GEOMETRY_ONLY"
            and all(leaf[k]==v for k,v in expect.items())
            and all(summary[k]==leaf[k] for k in ("pilot_id","record_id","radius_exact","box_sha256","body_geometry_sha256"))
            and leaf["all_144_body_and_24_trace_domains"] is True and leaf["closed_full_nodal_coordinates_rebuilt"]==650
            and leaf["raw_surface_cache_geometry_lines"]==288 and leaf["root_asserted"] is False,
            "admit complete independently reviewed same trial geometry")


def identity_key(p):
    return {k:p[k] for k in ("source_sha256","protocol_sha256","preflight_sha256","manifest_sha256","target_manifest_sha256")}


def gate(args):
    p=read(PREP/"protocol.json",args.protocol_sha256)
    f=read(PREP/"preflight.json",args.preflight_sha256)
    m=read(PREP/"manifest.json",args.manifest_sha256)
    require(args.review.resolve()==COLD_REVIEW_PATH.resolve(),"canonical independent cold bridge review")
    review=read(args.review,args.review_sha256)
    ids={"source_sha256":sha(SOURCE),"protocol_sha256":args.protocol_sha256,"preflight_sha256":args.preflight_sha256,
         "manifest_sha256":args.manifest_sha256,"target_manifest_sha256":V.T_SHA}
    require(p["source_sha256"]==m["source_sha256"]==f["source_sha256"]==ids["source_sha256"]
            and all(x["target_manifest_sha256"]==V.T_SHA for x in (p,f,m))
            and f["manifest_sha256"]==ids["manifest_sha256"] and f["protocol_sha256"]==ids["protocol_sha256"]
            and m["protocol_sha256"]==ids["protocol_sha256"] and f["status"]=="PASS" and V.yes(f["checks"])
            and all(V.zero(f[k]) for k in ZERO) and sha(PREP/"bridge_frozen.py")==ids["source_sha256"]
            and tuple(x["case_id"] for x in m["cases"])==V.A.IDS
            and review["status"]==REVIEW and V.yes(review["checks"]) and identity_key(review)==ids
            and review["selector_only_authorized_in_this_gate"] is True
            and review["spatial_workers_authorized_in_this_gate"] is False,
            "exact frozen preparer independently reviewed for selector stage only")
    for name,h in p["static_input_sha256"].items():
        require(sha(canonical(name))==h,"current frozen authority drift")
    return m,{**ids,"cold_review_sha256":args.review_sha256}


def admitted_domain(args):
    m,bindings=gate(args)
    entry=V.A.pick(m["cases"],args.case_id)
    require(args.admission.resolve()==JOINED.resolve(),"canonical joined external bridge input gate")
    phase=read(args.admission,args.admission_sha256)
    require(phase["status"]==ADMISSION and V.yes(phase["checks"]) and identity_key(phase)==identity_key(bindings)
            and phase["cold_review_sha256"]==args.review_sha256
            and phase["same_full_32_geometry_10_fresh_derivative_required"] is True
            and phase["spatial_workers_authorized"] is False,"reviewed joined input admission is not interchangeable")
    cases=tuple(e["case_id"] for e in m["cases"] if e["cached_input_kind"]=="FRESH_R34")
    geom,raw_geom,geo_phase=native_stage("geometry",phase["geometry_observed_exit_sha256"],V.A.IDS,phase["geometry_replay_phase"])
    fresh,raw_fresh,der_phase=native_stage("variation",phase["variation_observed_exit_sha256"],cases)
    case=args.case_id
    tr=entry["trial"]
    geo=geom[case]
    geopath=GEOM_REC/(case+".json")
    require(geo["receipt"]==rel(geopath),"same canonical independent geometry leaf")
    g=read(geopath,geo["sha256"])
    match_geometry(entry,g,geo)
    gfolder=V.GOUT/"actual/cases"/case
    gh=raw_geom[rel(gfolder/"results.json")]
    bodyh=raw_geom[rel(gfolder/"body_geometry.json")]
    require(g["body_geometry_sha256"]==bodyh and g["producer_case_results_sha256"]==gh
            and g["producer_case_results_file"]==rel(gfolder/"results.json")
            and g["producer_case_body_file"]==rel(gfolder/"body_geometry.json")
            and g["producer_worker_terminal_sha256"]==raw_geom[rel(gfolder/"terminal.json")]
            and g["producer_worker_observed_os_exit_sha256"]==raw_geom[rel(gfolder/"observed_os_exit.json")]
            and g["replay_phase"]==geo_phase,
            "same root admitted geometry producer bytes")
    gp=read(gfolder/"results.json",gh)
    body=read(gfolder/"body_geometry.json",bodyh)
    require(gp["trial"]==tr and gp["whole_body_sha256"]==bodyh and gp["case_id"]==case
            and gp["status"]=="PASS_CANDIDATE_WHOLE_BOX_GEOMETRY_AWAITS_INDEPENDENT_REVIEW"
            and body["status"]=="CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING"
            and len(gp["trace_geometry_only"])==24 and len(body["triangles"])==len(gp["bulk_geometry_only"])==144,
            "whole volume and trace domains at same original trial")
    if entry["cached_input_kind"]=="FRESH_R34":
        fs=fresh[case]
        leafpath=DER_REC/(case+".json")
        require(fs["native_file"]==rel(leafpath),"canonical fresh variation leaf")
        leaf=read(leafpath,fs["native_sha256"])
        z,norm,signed=match_derivative(entry,leaf)
        vp=V.DOUT/"actual/cases"/case/"results.json"
        require(fs["producer_result_sha256"]==leaf["producer_result_sha256"]==raw_fresh[rel(vp)]
                and fs["status"]==leaf["status"] and fs["reviewed_ZQ_upper_exact"]==leaf["reviewed_ZQ_upper_exact"]
                and leaf["phase"]==der_phase,
                "same native fresh replay and original candidate")
        derivative_hash=fs["native_sha256"]
    else:
        z,norm,signed=native_derivative(entry)
        derivative_hash=entry["native_case_sha256"]
    ctx,rule,box,route=V.current_context(tr)
    require(catalogue(route)==entry["component_catalogue"] and ctx.record_id==entry["record_id"],
            "authenticated reference route and target")
    rnorm=max(sum(abs(V.A.T.rational(x)) for x in row) for row in ctx.R)
    require(rnorm==norm,"independently supplied R norm must match all exact binary entries")
    size=sizes(entry["component_catalogue"],z,rnorm)
    point=ctx.target.nodes()
    d={"ctx":ctx,"rule":rule,"entry":entry,"box":box,"point":point,"point_box":B.exact_point_inside(point,box),
       "nodes":E.SUR.common_nodes(ctx.target,box),"routes":route,"geometry":gp,"body":body}
    inputs={**bindings,"bridge_input_admission_sha256":args.admission_sha256,
            "geometry_observed_exit_sha256":phase["geometry_observed_exit_sha256"],
            "variation_observed_exit_sha256":phase["variation_observed_exit_sha256"],
            "geometry_replay_case_sha256":geo["sha256"],"cached_derivative_case_sha256":derivative_hash}
    return d,size,signed,inputs

def plan_case(args):
    d,sizing,signed,links=admitted_domain(args)
    folder=PLAN/args.case_id
    folder.mkdir(parents=True,exist_ok=False)
    pid=os.getpid()
    mark={"pid":pid,"case_id":args.case_id,"record_id":d["entry"]["record_id"],"bindings":links}
    save(folder/"started.json",{"status":"RUNNING_SELECTOR_PLAN",**mark})
    end={"status":"FAILED","exit_code":1,**mark}
    components=[]
    try:
        if sizing["positive_headroom"] is not True:
            h=save(folder/"unresolved.json",{"status":"UNRESOLVED_CACHED_DERIVATIVE_HEADROOM",**mark,"sizing":sizing})
            end.update(status="FINISHED_UNRESOLVED",exit_code=0,unresolved_sha256=h)
            return
        tail=tuple(map(Q,sizing["caps_exact"]))
        force,scaled=kernels()
        with B.planning_only(force,scaled):
            for expected in d["entry"]["component_catalogue"]:
                kind,slot=expected["kind"],expected["slot"]
                try:
                    val=B.build_bulk(d,slot,tail,scaled) if kind=="bulk" else B.build_surface(d,slot,tail)
                except (E.TailBudgetUnresolved,B.OrderUnresolved) as error:
                    h=save(folder/"unresolved.json",{"status":"UNRESOLVED_SELECTOR_ORDER",**mark,
                           "kind":kind,"slot":slot,"attempts":error.attempts,"component_scope":error.component_scope})
                    end.update(status="FINISHED_UNRESOLVED",exit_code=0,unresolved_sha256=h)
                    return
                require((val["kind"],val["slot"],val["global_indices"])==(kind,slot,expected["global_indices"]),
                        "exact projected component order and map unchanged")
                name=f"{kind}_{slot:03d}.json"
                h=save(folder/name,V.A.doc(val))
                components.append({"kind":kind,"slot":slot,"file":name,"sha256":h})
        require(tuple((r["kind"],r["slot"]) for r in components)==OPS,"complete 168 selectors required")
        plan={"status":"COMPLETE_REMAINING32_CASE_SELECTOR_PLAN_AWAITS_INDEPENDENT_REVIEW",**mark,
              "whole_box_radius_exact":str(RADIUS),"trial_box_sha256":d["entry"]["trial"]["trial_box_sha256"],
              "point_box_sha256":B.digest(d["point_box"]),"sizing":sizing,"signed_cached_RF_center_rows_exact":signed,
              "point_force_scope":POINT,"whole_H_scope":BOX,"components":components,
              "new_spatial_integrations":0,"spatial_workers_authorized":False}
        h=save(folder/"plan.json",plan)
        end.update(status="FINISHED",exit_code=0,plan_sha256=h)
        print({"case_id":args.case_id,"plan_sha256":h,"components":len(components)},flush=True)
    except BaseException as e:
        save(folder/"failure.json",{**mark,"error":str(e),"traceback":traceback.format_exc(),
                                     "completed_components":components})
        raise
    finally:
        save(folder/"terminal.json",end)


def validate_plan(outer,term,p,review,plan_sha,term_sha,bindings,case,ids):
    require(outer["status"]=="PASS_OBSERVED_REMAINING32_SELECTOR_TRUE_OS_EXIT"
            and V.zero(outer["observed_exit_code"]) and V.zero(term["exit_code"])
            and all(V.positive_int(x["pid"]) for x in (outer,term,p))
            and outer["pid"]==term["pid"]==p["pid"] and V.positive_int(outer["unified_exec_session"])
            and outer["terminal_sha256"]==term_sha
            and outer["plan_sha256"]==term["plan_sha256"]==plan_sha
            and all(x["case_id"]==case for x in (outer,term,p,review))
            and p["status"]=="COMPLETE_REMAINING32_CASE_SELECTOR_PLAN_AWAITS_INDEPENDENT_REVIEW"
            and term["status"]=="FINISHED" and p["bindings"]==term["bindings"]==bindings
            and review["status"]==PLAN_REVIEW and V.yes(review["checks"])
            and review["independently_verified_component_count"]==len(ids)
            and review["plan_sha256"]==plan_sha and review["terminal_sha256"]==term_sha
            and review["bindings"]==bindings,
            "actual selector status, OS and independent full plan review required")
    require(len(p["components"])==len(ids) and all(type(x["slot"]) is int for x in p["components"])
            and tuple((x["kind"],x["slot"]) for x in p["components"])==ids,
            "complete typed ordered current plan")


def open_reviewed_plan(args,outer_sha,review_path,review_sha):
    """Future worker input path. Integrated numbers need separate packet replay/reduction."""
    d,sizing,signed,links=admitted_domain(args)
    folder=PLAN/args.case_id
    os_path=folder/"root_os_exit_observation.json"
    outer=read(os_path,outer_sha)
    p=read(folder/"plan.json",outer["plan_sha256"])
    term=read(folder/"terminal.json",outer["terminal_sha256"])
    review=read(review_path,review_sha)
    require(review["observed_os_exit_sha256"]==outer_sha,"same independently observed real selector exit")
    validate_plan(outer,term,p,review,outer["plan_sha256"],outer["terminal_sha256"],links,args.case_id,OPS)
    require(p["whole_box_radius_exact"]==str(RADIUS) and p["trial_box_sha256"]==d["entry"]["trial"]["trial_box_sha256"]
            and p["point_box_sha256"]==B.digest(d["point_box"]) and signed==p["signed_cached_RF_center_rows_exact"]
            and p["sizing"]==sizing and p["spatial_workers_authorized"] is False
            and (p["point_force_scope"],p["whole_H_scope"])==(POINT,BOX),
            "same center/cache/geometry and exact whole-box tail planning")
    packets={}
    for chosen,cold in zip(p["components"],d["entry"]["component_catalogue"],strict=True):
        kind,slot=chosen["kind"],chosen["slot"]
        require(chosen["file"]==f"{kind}_{slot:03d}.json","canonical independent plan leaf path")
        val=read(folder/chosen["file"],chosen["sha256"])
        require((val["kind"],val["slot"],val["global_indices"])==(kind,slot,cold["global_indices"]),
                "same complete projected component scope")
        packets[(kind,slot)]=val
    return d,tuple(map(Q,sizing["caps_exact"])),packets
        

def evaluate_reviewed_component(args,outer_sha,review_path,review_sha,kind,slot):
    """No CLI dispatch. A separately reviewed worker can use this exact math path."""
    require(type(slot) is int and (kind,slot) in OPS,"explicit canonical component required")
    d,tail,plans=open_reviewed_plan(args,outer_sha,review_path,review_sha)
    force,scaled=kernels()
    return B.evaluate(d,plans[(kind,slot)],tail,force,scaled)


def synthetic():
    rows=[{"kind":"bulk","slot":0,"global_indices":[0,1]},
          {"kind":"surface","slot":0,"global_indices":[1]}]
    out=sizes(rows,Fraction(3,8),Fraction(2),2)
    require(out["c_F"]==2 and out["c_H"]==3 and out["positive_headroom"] is True
            and out["quarter_ZQ_admission_rule"] is False
            and tuple(map(Q,out["caps_exact"]))==(Fraction(5,32)*RADIUS/40000,Fraction(5,16)/60000),
            "actual norm/incident tail accounting synthetic")
    blocked=sizes(rows,1,2,2)
    altered=[{**rows[0],"global_indices":[False,1]},rows[1]]
    case="fixture"
    h="c"*64
    t="d"*64
    b={"source_sha256":"1"*64}
    p={"status":"COMPLETE_REMAINING32_CASE_SELECTOR_PLAN_AWAITS_INDEPENDENT_REVIEW","pid":5,
       "case_id":case,"bindings":b,"components":[{"kind":"bulk","slot":0}]}
    term={"status":"FINISHED","pid":5,"case_id":case,"bindings":b,"exit_code":0,"plan_sha256":h}
    outer={"status":"PASS_OBSERVED_REMAINING32_SELECTOR_TRUE_OS_EXIT","unified_exec_session":6,
           "case_id":case,"pid":5,"observed_exit_code":0,"terminal_sha256":t,"plan_sha256":h}
    review={"status":PLAN_REVIEW,"case_id":case,"checks":{"synthetic":True},
            "independently_verified_component_count":1,"plan_sha256":h,"terminal_sha256":t,"bindings":b}
    validate_plan(outer,term,p,review,h,t,b,case,(("bulk",0),))
    def reject_alt(**new):
        return V.A.reject(lambda:validate_plan({**outer,**new},term,p,review,h,t,b,case,(("bulk",0),)))
    return {"signed_scopes_and_component_projection_type_guards":V.A.reject(lambda:routes_order(altered,2)),
        "above_quarter_headroom_and_exact_incidence_caps_checked":True,
        "no_optimization_budgets_for_nonpositive_headroom":blocked["caps_exact"] is None,
        "reject_wrong_plan_Flag_PID_or_case":all((reject_alt(observed_exit_code=False),reject_alt(pid=5.0),
                                                 reject_alt(case_id="sibling"),reject_alt(plan_sha256="0"*64)))}


def prepare():
    cases,links=premises()
    force,scaled=kernels()
    checks=synthetic()
    with B.planning_only(force,scaled):
        try:
            force._force_parts()
        except AssertionError:
            checks["cold_cannot_invoke_expensive_spatial_kernels"]=True
    require(V.yes(checks) and len(cases)==32 and sum(len(z["component_catalogue"]) for z in cases)==5376,
            "complete explicit cold population and checks required")
    PREP.mkdir(parents=True,exist_ok=False)
    binding={"source_sha256":sha(SOURCE),"target_manifest_sha256":V.T_SHA}
    protocol={"status":"FROZEN_REMAINING32_EXACT_BRIDGE_COLD_PROTOCOL",**binding,**ZERO,"scope":SCOPE,
              "static_input_sha256":{**links,**{rel(p):h for p,h in PINS.items()}},
              "closed_radius_exact":str(RADIUS),"case_order":list(V.A.IDS),
              "cached_derivative_coverage":{"native_equal_r34":22,"fresh_independent_r34":10},
              "joined_admission_canonical_file":rel(JOINED),"cold_review_canonical_file":rel(COLD_REVIEW_PATH),
              "required_independent_actual_evidence":{"all_geometry_case_count":32,"all_fresh_derivative_case_count":10,
                   "geometry_raw_external_admission_sha256":RAW["geometry"][1],
                   "geometry_replay_phase_keys":["source_sha256","protocol_sha256","preflight_sha256"],
                   "derivative_raw_external_admission_sha256":RAW["variation"][1],
                   "joined_status":ADMISSION,"cold_selector_review_status":REVIEW,"plan_review_status":PLAN_REVIEW},
              "spatial_worker_cli":False,"root_asserted":False,"full_five_gate_complete":False}
    ph=save(PREP/"protocol.json",protocol)
    m={"status":"COLD_CANONICAL_ALL_32_COMPONENT_CATALOGUES",**binding,**ZERO,"protocol_sha256":ph,
       "case_count":32,"components_per_case":168,"components_catalogued":5376,"cases":cases,
       "no_spatial_workers_authorized":True}
    mh=save(PREP/"manifest.json",m)
    pf=save(PREP/"preflight.json",{"status":"PASS",**binding,**ZERO,"protocol_sha256":ph,
                                  "manifest_sha256":mh,"checks":checks})
    with (PREP/"bridge_frozen.py").open("xb") as f:
        f.write(SOURCE.read_bytes())
    print({**binding,"protocol_sha256":ph,"manifest_sha256":mh,"preflight_sha256":pf},flush=True)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    sub=cli.add_subparsers(dest="phase",required=True)
    sub.add_parser("prepare")
    plan=sub.add_parser("plan")
    plan.add_argument("--case-id",required=True,choices=V.A.IDS)
    for item in ("protocol","preflight","manifest","review","admission"):
        plan.add_argument("--"+item+"-sha256",required=True)
        if item in ("review","admission"):
            plan.add_argument("--"+item,required=True,type=lambda v:Path(v).resolve())
    args=cli.parse_args()
    prepare() if args.phase=="prepare" else plan_case(args)
