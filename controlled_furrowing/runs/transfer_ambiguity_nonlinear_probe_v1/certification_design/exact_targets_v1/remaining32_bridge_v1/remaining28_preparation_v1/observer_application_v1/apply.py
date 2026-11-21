"""Transfer later-wave exact spatial root bounds to the frozen archived observers."""
import argparse
from collections import Counter
from fractions import Fraction as Q
import importlib.util
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
LATER=HERE.parent
BR=LATER.parent
sys.path[:0]=[str(LATER/"global_reduction_v1"),str(LATER)]
import reduce as G

X=G.X
S=G.S
V=G.V
FIRST=BR/"first_wave_observer_v1/apply.py"
FIRST_SHA="6db0999d7a11c6bcc7fbd4030a890c7ea80dcca0a41264f40329fb7070264b37"
FIRST_REVIEW="9876af4b080d59c78045c0d1b815b6dc7374438efc751d8352991cd40d7125cc"
GLOBAL_PHASE={
    "source_sha256":"3051ac6d35ee18407e29e8b82e96d789e685073dd3f067b01b50b459c0b970e5",
    "protocol_sha256":"1198124fe22265e81113a3f383a82c6b816334e5d3a54df837bd16eb491aadea",
    "preflight_sha256":"b3310545e2ab4e15e5988c90d4011565eebd5797f0c9449c75bcfad52a5f8b88",
    "source_review_sha256":"f39bb78a6ee4ee66cfdf25dc287a4856a1f7b5f94298c612fbbc8c17794d171f",
    "later_waves_release_sha256":G.RELEASE,"profiles_sha256":X.PROFILE_SHA,
}
STATUS="FROZEN_LATER_WAVE_EXACT_SPATIAL_OBSERVER_TRANSFER"
REVIEW="PASS_INDEPENDENT_LATER_WAVE_OBSERVER_APPLICATION_SOURCE"
OBSERVED="PASS_ROOT_OBSERVED_LATER_WAVE_INDEPENDENT_GLOBAL_EXIT"
ROOT_PASS="PASS_INDEPENDENT_LATER_WAVE_EXACT_SPATIAL_ROOT"
EXCLUDE="SUFFICIENT_INDEPENDENT_EXCLUSION_OF_ONLY_THIS_CUBE"
UNRESOLVED="UNRESOLVED_INDEPENDENT_LATER_WAVE_EXACT_SPATIAL_TEST"
require,sha,rel,save=G.require,G.sha,G.rel,G.save
RADIUS=G.RADIUS
SCOPE="Frozen primary and alternative model observers of fixed spatial finite-element roots. No new fit, loading path, continuum, or source compatibility result."


def helpers():
    prior=G.read(FIRST.with_name("independent_source_review")/"receipt.json",FIRST_REVIEW)
    require(sha(FIRST)==FIRST_SHA==prior["source_sha256"] and V.yes(prior["checks"]),
            "previously reviewed pure observer geometry and interval helpers")
    spec=importlib.util.spec_from_file_location("_prior_pure_observer_functions",FIRST)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    require(m.RADIUS==RADIUS and m.LENGTH==Q(67,2), "same dimensionless radius and physical length")
    return m


def archive(M,links):
    # Only archive, budget, and observer theorem artifacts are inherited.
    bnames=(M.PRIMARY,M.CACHED/"results.json",M.CACHED/"independent_readback/receipt.json",
            M.BUDGET/"actual/results.json",M.BUDGET/"independent_actual/receipt.json",
            M.SECONDARY/"actual/results.json",M.SECONDARY/"independent_actual/receipt.json",M.GEOMETRY)
    docs={p:M.read(p,M.PINS[p],links) for p in bnames}
    for p in (M.ROOT/"docs/NODAL-OBSERVER-ERROR-TRANSFER.md",M.ROOT/"evidence/nodal_observer_error_independent_review.json"):
        M.read(p,M.PINS[p],links,False)
    ca,ci=(docs[M.CACHED/x] for x in ("results.json","independent_readback/receipt.json"))
    ba,bi=(docs[M.BUDGET/x] for x in ("actual/results.json","independent_actual/receipt.json"))
    sa,si=(docs[M.SECONDARY/x] for x in ("actual/results.json","independent_actual/receipt.json"))
    a=M.exact(ca["frozen_allowance_um_exact"])
    require(ci["status"]=="PASS_INDEPENDENT_ACTUAL_36_CACHED_ROOT_OBSERVERS"
            and ci["current_input_sha256"][rel(M.CACHED/"results.json")]==M.PINS[M.CACHED/"results.json"]
            and bi["status"]=="PASS_INDEPENDENT_ACTUAL_CONDITIONAL_OBSERVER_BUDGETS"
            and bi["source_and_input_hashes"][rel(M.BUDGET/"actual/results.json")]==M.PINS[M.BUDGET/"actual/results.json"]
            and si["status"]=="PASS_INDEPENDENT_SECONDARY_STORED_OBSERVER_SWEEP"
            and si["input_and_output_sha256"][rel(M.SECONDARY/"actual/results.json")]==M.PINS[M.SECONDARY/"actual/results.json"]
            and 0<a==M.exact(ci["frozen_model_discrepancy_um_exact"])==M.exact(ba["conditional_source_allowance_um_exact"])
            ==M.exact(si["source_configured_delta_um_exact"])==M.exact(sa["frozen_primary_source_allowance_um_exact"]),
            "same source-conditioned allowance and independently observed archive")
    return docs,a


def source_inputs():
    M=helpers()
    p,profiles=G.static()
    names={"source_sha256":G.SOURCE,"protocol_sha256":G.HERE/"protocol.json",
           "preflight_sha256":G.HERE/"preflight.json","source_review_sha256":G.AUTH}
    hashes={}
    for k,path in names.items():
        M.read(path,GLOBAL_PHASE[k],hashes,False)
    G.approved(G.read(names["protocol_sha256"],GLOBAL_PHASE["protocol_sha256"]),
               G.read(names["preflight_sha256"],GLOBAL_PHASE["preflight_sha256"]),
               G.read(G.AUTH,GLOBAL_PHASE["source_review_sha256"]),{k:v for k,v in GLOBAL_PHASE.items() if k!="source_review_sha256"})
    require(G.read(names["protocol_sha256"],GLOBAL_PHASE["protocol_sha256"])==p,
            "same frozen generic global source phase")
    M.read(FIRST,FIRST_SHA,hashes,False)
    M.read(FIRST.with_name("independent_source_review")/"receipt.json",FIRST_REVIEW,hashes,False)
    docs,allowance=archive(M,hashes)
    return M,profiles,docs,allowance,hashes


def bind_archive(M,docs,hashes,name,target):
    endpoint,number=name.rsplit("_small_",1)
    checkpoint="small_"+number
    trial=M.pick(docs[M.GEOMETRY]["trials"],case_id=name)
    contract=M.read(trial["model_contract"],trial["model_contract_sha256"],hashes)
    keys=("trial_box_sha256","model_contract_sha256","point_cache_sha256","R_binary_sha256",
          "center_binary_sha256","archived_center_shift_exact","exact_drive","energy_scale_binary_exact")
    require(name==target["case_id"]==contract["pilot_id"]
            and target["record_id"]==trial["record_id"]==contract["record_id"]
            and all(target[k]==trial[k] for k in keys)
            and M.exact(trial["trial_radius_exact"])==RADIUS,"same later singleton and exact proof cube")
    box=M.read(trial["trial_box"],trial["trial_box_sha256"],hashes)
    require(box["record_id"]==target["record_id"] and box["point_cache_sha256"]==target["point_cache_sha256"]
            and box["model_observation_contract_sha256"]==target["model_contract_sha256"]
            and M.exact(box["corrector_radius_exact"])==RADIUS,"same case-specific exact nodal cube")
    prior=M.pick(docs[M.CACHED/"results.json"]["cases"],case_id=name)
    old=docs[M.CACHED/"independent_readback/receipt.json"]
    old_terms=M.pick(old["observer_transfer_terms"],case_id=name)
    primary=docs[M.PRIMARY]
    pilot=docs[M.BUDGET/"actual/results.json"]
    b=M.pick(pilot["field_budgets"],endpoint_id=endpoint,checkpoint=checkpoint)
    audit=M.pick(docs[M.BUDGET/"independent_actual/receipt.json"]["fields"],endpoint_id=endpoint,checkpoint=checkpoint)
    sweep=M.SECONDARY/"actual"/endpoint/(checkpoint+".json")
    sp=M.pick(docs[M.SECONDARY/"actual/results.json"]["field_observer_receipts"],path=rel(sweep))["sha256"]
    secondary=docs[M.SECONDARY/"independent_actual/receipt.json"]["input_and_output_sha256"]
    saved=M.read(sweep,sp,hashes)
    model=b["center_and_model"]
    require(saved["role"]==b["role"]=="small_record" and saved["endpoint_id"]==endpoint and saved["checkpoint"]==checkpoint
            and saved["step"]==b["step"] and secondary[rel(sweep)]==sp
            and all(saved[k]==b[k] for k in ("field","field_sha256","verification","verification_sha256","record_id"))
            and prior["record_id"]==audit["record_id"]==b["record_id"]==contract["physical_source_record_id"]
            and prior["field_sha256"]==audit["field_sha256"]==saved["field_sha256"]
              ==contract["archived_field_sha256"]==target["original_observer_field_sha256"]
            and contract["archived_field"]==audit["field"]==saved["field"]
            and contract["archived_verification"]==saved["verification"]
            and contract["archived_verification_sha256"]==saved["verification_sha256"]
            and target["observer_goal_exact"]==audit["e_goal_dimensionless_exact"]==b["e_goal_dimensionless_exact"]==contract["observer_goal_exact"]
            and target["archived_center_shift_exact"]==contract["archived_shift_exact"]==prior["root"]["sigma"]==old_terms["sigma_exact"]
            and target["original_addition_gap_exact"]==contract["archived_addition_gap_exact"]==audit["binary_gap_exact"]
              ==model["node_addition_gap_dimensionless_exact"]==old_terms["archived_xi_exact"]
            and M.exact(contract["archived_observer_allowance_um_exact"])==M.exact(pilot["conditional_source_allowance_um_exact"])
            and contract["archived_observer_variant"]==M.VARIANTS[0],"exact archived source field and observer contract")
    dt=contract["exact_requested_drive"]
    dose,gain,strength=(M.exact(dt[k]) for k in ("optical_dose_exact","common_gain_exact","activation_strength_exact"))
    require(dose==M.exact(saved["step"]["dose_exact"])==M.exact(box["dose_exact"])
            and gain==M.exact(model["gain_exact"]) and 0<strength==dose*gain==M.exact(model["exact_target_strength"])
            ==M.exact(target["exact_drive"]),"same exact dose and gain")
    opened=[]
    for name,h in ((saved["field"],saved["field_sha256"]),(saved["verification"],saved["verification_sha256"])):
        require(primary["all_hashes"][name]==old["current_input_sha256"][name]==secondary[name]==h,"three current sources agree")
        opened.append(M.read(name,h,hashes))
    state,verification=opened
    require(verification["variant"]==M.VARIANTS[0] and set(saved["variant_observables_um_exact"])==set(M.VARIANTS)
            and all(set(saved["variant_observables_um_exact"][v])==set(M.OBS) for v in M.VARIANTS),"all frozen sixteen original outputs")
    for v in M.OBS:
        primary_value=M.pick(prior["comparisons"],observable=v)
        require(M.pair(saved["variant_observables_um_exact"][M.VARIANTS[0]][v])==M.pair(primary_value["stored_interval_um_exact"])
                ==M.pair(verification["observer_intervals_um_exact"][v])
                and M.exact(b["step"]["observed_exact"][v])==M.exact(primary_value["observed_um_exact"]),"same exact primary observation")
    return b,prior,saved,state,box


def native_gate(M,wave,profile,outer_sha,links):
    root=X.wave_folder(wave)/"global_reduction_v1/independent"
    obs=M.read(root/"root_observed_exit.json",outer_sha,links)
    cases=profile["case_order"]
    require(obs["status"]==OBSERVED and obs["wave_id"]==wave and obs["mode"]=="independent"
            and G.known_exit(obs),"genuine externally observed independent generic global terminal")
    pins=obs["sha256"]
    expected={rel(root/(k+".json")) for k in ("started","terminal","results")}
    expected|={rel(root/c/(k+".json")) for c in cases for k in ("operands","results")}
    require(set(pins)==expected and len(expected)==11,"all signed global actual artifacts externally pinned")
    start,end,final=(M.read(root/(k+".json"),obs[k+"_sha256"],links) for k in ("started","terminal","results"))
    phase={**GLOBAL_PHASE,"mode":"independent","wave_id":wave,"pid":obs["pid"]}
    for k in ("hull_admission_sha256","raw_wave_admission_sha256","raw_root_observed_exit_sha256","producer_observed_os_sha256"):
        phase[k]=final[k]
    require(all(V.positive_int(x["pid"]) and all(x[k]==v for k,v in phase.items()) for x in (start,end,final))
            and start["status"]=="RUNNING" and end["status"]=="FINISHED" and V.zero(end["exit_code"])
            and end["results_sha256"]==obs["results_sha256"] and final["status"]==G.INDEP_STATUS
            and final["signed_hull_components"]==672 and V.zero(final["new_force_AD_or_spatial_integrals"])
            and final["full_five_gate_complete"] is False
            and tuple(x["case_id"] for x in final["cases"])==tuple(cases),"complete later independent global phase")
    for k in ("started","terminal","results"):
        require(obs[k+"_sha256"]==pins[rel(root/(k+".json"))],"externally pinned closure metadata")
    n=SimpleNamespace(wave=wave,protocol_sha256=phase["protocol_sha256"],preflight_sha256=phase["preflight_sha256"],
          review_sha256=phase["source_review_sha256"],hull_admission_sha256=phase["hull_admission_sha256"],
          raw_wave_admission_sha256=phase["raw_wave_admission_sha256"],raw_root_observed_exit_sha256=phase["raw_root_observed_exit_sha256"],
          producer_os_sha256=phase["producer_observed_os_sha256"])
    G.producer_observed(n,root.parent,cases)
    rows={}
    for e in final["cases"]:
        case=e["case_id"]
        f=root/case
        require(e["case_results_sha256"]==pins[rel(f/"results.json")] and e["operands_sha256"]==pins[rel(f/"operands.json")],
                "native per-case results and signed operands exact external digests")
        leaf=M.read(f/"results.json",e["case_results_sha256"],links)
        M.read(f/"operands.json",e["operands_sha256"],links,False)
        require(all(leaf[k]==v for k,v in phase.items()) and V.positive_int(leaf["pid"])
                and leaf["case_id"]==case and leaf["status"]==e["status"]
                and leaf["operands_sha256"]==e["operands_sha256"] and leaf["target"]==profile["raw_worker_targets"][case]
                and M.exact(leaf["radius_exact"])==RADIUS and leaf["full_five_gate_complete"] is False
                and leaf["archived_observer_field_replaced"] is False,"independent generic leaf same target and exact phase")
        rows[case]=(e,leaf)
    return rows


def locate(M,row,leaf):
    t=leaf["target"]
    a=leaf["arithmetic"]
    v=leaf["observer_and_optional_minimum"]
    eta,z,lo=(M.exact(a[k]) for k in ("eta_upper_exact","z_upper_exact","RF_lower_exact"))
    r=M.exact(leaf["radius_exact"])
    lhs=eta+z*r
    exists=z<1 and lhs<r
    separate=lo>(1+z)*r
    status=ROOT_PASS if exists else EXCLUDE if separate else UNRESOLVED
    require(0<=lo<=eta and 0<=z and r==RADIUS and not (exists and separate)
            and a["strict_selfmap"] is exists and a["strict_cube_exclusion"] is separate
            and M.exact(a["selfmap_lhs_upper_exact"])==lhs
            and leaf["status"]==row["status"]==status and leaf["root_asserted"] is exists
            and row["exact_spatial_root_certified"] is exists,"strict exact same-cube scalar root decision")
    if not exists:
        require(a["posterior_error_upper_exact"] is None and row["shift_inclusive_observer_budget_proved"] is False,
                "missing proof retains unresolved observer application")
        return None
    e=eta/(1-z)
    sigma,xi,goal=(M.exact(t[k]) for k in ("archived_center_shift_exact","original_addition_gap_exact","observer_goal_exact"))
    epsilon=M.LENGTH*(sigma+e+xi)
    require(0<=e<r and sigma>=0 and xi>=0 and M.exact(a["posterior_error_upper_exact"])==e
            and v["archived_state_budget_proved"] is (sigma+e<goal)
            and row["shift_inclusive_observer_budget_proved"] is (sigma+e<goal)
            and M.exact(v["archived_shift_plus_posterior"])==sigma+e
            and M.exact(v["physical_nodal_error_um_upper_exact"])==epsilon,"preserved archive sigma, added xi, and spatial posterior")
    return e,sigma,xi,epsilon


def apply_case(M,case,target,row,leaf,archived,allowance):
    b,prior,saved,state,box=archived
    p=locate(M,row,leaf)
    roi=M.pair(prior["ROI_um_exact"])
    graph=None
    reason="ROOT_NOT_CERTIFIED_IN_THIS_CUBE" if p is None else None
    if p is not None:
        e,sigma,xi,epsilon=p
        if any(M.exact(b["stored_face_derivatives"][s]["radial_min_um_exact"])<=8*epsilon for s in ("top","bottom")):
            reason="STRICT_RADIAL_GRAPH_GUARD_NOT_ESTABLISHED"
        else:
            roi,graph=M.geometry(b,prior,saved,state,box,sigma,xi,epsilon)
    errors={}
    if graph is not None:
        ea,eb=(M.exact(graph[s]["graph_error_um_exact"]) for s in ("apical","basal"))
        errors=dict(zip(M.OBS,(ea+eb,ea+eb,ea,eb),strict=True))
    out=[]
    for v in M.VARIANTS:
        for name in M.OBS:
            stored=M.pair(saved["variant_observables_um_exact"][v][name])
            sample=M.exact(b["step"]["observed_exact"][name])
            if reason:
                eps=wide=err=None
                standing=reason
            else:
                eps=str(errors[name])
                wide,err,standing=M.compare(stored,sample,errors[name],allowance)
                wide,err=list(map(str,wide)),list(map(str,err))
            out.append({"variant":v,"observable":name,"stored_interval_um_exact":list(map(str,stored)),
                        "observed_um_exact":str(sample),"extra_observer_error_um_exact":eps,
                        "exact_spatial_root_observer_interval_um_exact":wide,"absolute_residual_um_exact":err,"standing":standing})
    return {"case_id":case,"record_id":b["record_id"],"root_proof_sha256":row["case_results_sha256"],
            "native_root_status":leaf["status"],"original_field_sha256":target["original_observer_field_sha256"],
            "posterior_exact":None if p is None else str(p[0]),"archived_sigma_exact":target["archived_center_shift_exact"],
            "original_xi_exact":target["original_addition_gap_exact"],"physical_nodal_error_um_exact":None if p is None else str(p[3]),
            "ROI_um_exact":list(map(str,roi)),"graph_premises":graph,"comparisons":out}


def binding(args):
    return {"source_sha256":sha(Path(__file__)),"protocol_sha256":args.protocol_sha256,"preflight_sha256":args.preflight_sha256,
            "later_waves_release_sha256":G.RELEASE,"profiles_sha256":X.PROFILE_SHA,
            "global_source_sha256":GLOBAL_PHASE["source_sha256"],"global_source_review_sha256":GLOBAL_PHASE["source_review_sha256"]}


def prepare():
    M,profiles,docs,a,source=source_inputs()
    cases={}
    links=dict(source)
    for wave,p in profiles.items():
        cases[wave]=p["case_order"]
        for case in p["case_order"]:
            bind_archive(M,docs,links,case,p["raw_worker_targets"][case])
    require(len({x for c in cases.values() for x in c})==28,"cold exact archive and later profile scope")
    checks=M.fixtures()
    checks.update({"all28_same_archived_field_dose_box_ROI_source_metadata":True,
                   "generic_later_schema_expected_no_first_wave_native_authority":True,
                   "cold_not_opened_later_global_or_generated_new_observer":True})
    p=save(HERE/"protocol.json",{"status":STATUS,"source_sha256":sha(Path(__file__)),"source_and_archive_sha256":links,
        "global_phase":GLOBAL_PHASE,"ordered_cases":cases,"radius_exact":str(RADIUS),"allowance_um_exact":str(a),
        "variants":list(M.VARIANTS),"observables":list(M.OBS),"actual_observer_jobs":0})
    f=save(HERE/"preflight.json",{"status":"PASS","source_sha256":sha(Path(__file__)),"protocol_sha256":p,"checks":checks,
         "actual_observer_jobs":0,"opened_later_global_native_payloads":0,"actual_force_AD_or_spatial_jobs":0})
    print(json.dumps({"source_sha256":sha(Path(__file__)),"protocol_sha256":p,"preflight_sha256":f}),flush=True)


def run(args):
    expected=binding(args)
    p=G.read(HERE/"protocol.json",args.protocol_sha256)
    f=G.read(HERE/"preflight.json",args.preflight_sha256)
    g=G.read(HERE/"independent_source_review/receipt.json",args.review_sha256)
    require(g["status"]==REVIEW and V.yes(g["checks"]) and all(g[k]==v for k,v in expected.items())
            and g["approved_waves"]==list(X.WAVES) and g["new_spatial_integrations_authorized"] is False
            and p["status"]==STATUS and p["global_phase"]==GLOBAL_PHASE
            and p["source_sha256"]==f["source_sha256"]==expected["source_sha256"]
            and f["status"]=="PASS" and V.yes(f["checks"]) and f["protocol_sha256"]==args.protocol_sha256
            and all(V.zero(f[k]) for k in ("actual_observer_jobs","opened_later_global_native_payloads","actual_force_AD_or_spatial_jobs")),
            "current reviewed scalar observer source")
    M,profiles,docs,allowance,links=source_inputs()
    profile=profiles[args.wave]
    cases=profile["case_order"]
    require(p["ordered_cases"]=={k:v["case_order"] for k,v in profiles.items()} and Q(p["allowance_um_exact"])==allowance,
            "same later target catalogue and frozen allowance")
    rows=native_gate(M,args.wave,profile,args.observed_global_sha256,links)
    archived={k:bind_archive(M,docs,links,k,profile["raw_worker_targets"][k]) for k in cases}
    require(all(links[k]==p["source_and_archive_sha256"][k] for k in links if k in p["source_and_archive_sha256"]),
            "frozen archive readback")
    dest=X.wave_folder(args.wave)/"observer_application_v1/actual"
    dest.mkdir(parents=True,exist_ok=False)
    stamp={**expected,"pid":os.getpid(),"wave_id":args.wave,"review_sha256":args.review_sha256,
           "native_global_observed_exit_sha256":args.observed_global_sha256}
    save(dest/"started.json",{"status":"RUNNING",**stamp})
    term={"status":"FAILED","exit_code":1,**stamp}
    try:
        results=[apply_case(M,k,profile["raw_worker_targets"][k],*rows[k],archived[k],allowance) for k in cases]
        counts={v:dict(Counter(x["standing"] for r in results for x in r["comparisons"] if x["variant"]==v)) for v in M.VARIANTS}
        require(len(results)==4 and all(sum(d.values())==16 for d in counts.values()),"all 16 primary and 48 alternative comparisons")
        h=save(dest/"results.json",{"status":"LATER_WAVE_EXACT_SPATIAL_OBSERVERS_AWAIT_INDEPENDENT_READBACK",
              **stamp,"scope":SCOPE,"allowance_um_exact":str(allowance),"cases":results,"counts":counts,
              "total_primary_comparisons":16,"total_alternative_comparisons":48,"source_and_actual_sha256":links,
              "new_spatial_jobs":0,"full_calibration_or_source_fit_claimed":False,"full_five_gate_complete":False})
        term.update(status="FINISHED",exit_code=0,results_sha256=h)
        print(json.dumps({"results_sha256":h,"counts":counts}),flush=True)
    except Exception as err:
        save(dest/"failure.json",{"error":str(err),**stamp})
        raise
    finally:
        save(dest/"terminal.json",term)


if __name__=="__main__":
    cl=argparse.ArgumentParser()
    sub=cl.add_subparsers(dest="action",required=True)
    sub.add_parser("prepare")
    p=sub.add_parser("run")
    p.add_argument("--wave",choices=X.WAVES,required=True)
    for k in ("protocol","preflight","review","observed-global"):
        p.add_argument("--"+k+"-sha256",required=True)
    args=cl.parse_args()
    prepare() if args.action=="prepare" else run(args)
