"""Independent cold review of the remaining32 bridge adapter, with no selector work."""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace
import sys

HERE=Path(__file__).resolve().parent
BR=HERE.parent
EX=BR.parent
ROOT=next(p for p in EX.parents if (p/"PLAN.md").is_file())
P=BR/"preparation"
SOURCE=ROOT/"scripts/transfer_ambiguity_remaining32_exact_bridge.py"
TARGET=EX/"remaining32_targets_v1"
D=EX/"remaining32_derivative_extension_v1"
G=EX/"remaining32_geometry_v1"
DR=D/"independent_actual_preparation/actual_replay"
DP=D/"preparation/protocol.json"
GP=G/"preparation/protocol.json"
R=Q(1,2**34)
N=611
ORDER=tuple((k,i) for k,m in (("bulk",144),("surface",24)) for i in range(m))
CASES=tuple(f"{name}_small_{i:02d}" for name in ("max_contrast_change","max_negative_contrast_change") for i in range(1,17))
PHASE={
    "source_sha256":"4d52e8fb8dbc01ba58c38052bbff537443693234499a9a3b9b3e138c82d9b51b",
    "protocol_sha256":"ade6158b96bec63665880beaae1447057b36f162863211acf1f8d24f95519664",
    "manifest_sha256":"d42b7260283ba02ccc34b9e9ba69ce76d58507115c22615c90d28fa38020cb41",
    "preflight_sha256":"a6741fb642a4006b126ac71db4fdb964761f5245e2563cd70a0175db6ee18bda",
    "target_manifest_sha256":"dfea134640f99b45e9ca75d1ce718043d67f2917375506997d2465abda8c1242",
}
PINS={
    SOURCE:PHASE["source_sha256"], P/"bridge_frozen.py":PHASE["source_sha256"],
    P/"protocol.json":PHASE["protocol_sha256"], P/"manifest.json":PHASE["manifest_sha256"],
    P/"preflight.json":PHASE["preflight_sha256"], TARGET/"preparation/manifest.json":PHASE["target_manifest_sha256"],
    DP:"6d04fc82779430e66c2d27b5ab7aa508f522f7ffac03b87f8f2e8e2587e1d128",
    GP:"754da64bb96cd8a19707181315e6c3e3f014e8d37a898b83c1693342c82c9850",
    TARGET/"independent_prelaunch/review.json":"67d33bbc3171839b1c578b02b5cfbf9ea6ba5d6cb351d0655dd14aaf99df75f9",
    DR/"receipt.json":"2a1be56b3050215e05da8e11b4785a65cdc8a229210fa82ebee37f1aa8236d05",
    DR/"terminal.json":"7664d8609c005c262e5f18a4dc84e251d67bf315440c7c0307e4ea67b056f8eb",
    DR/"os_exit_observation.json":"69ededa4a104e07ac85e251e6a019b93fc11016cdd83954aaee9f7a987f83325",
    D/"actual/root_terminal_admission.json":"91b1db64f27c89f50643fdde4f16c70a0dcac9ee51566fcf5f352b5b1f663f18",
}
LINKS={}


def sha(b):
    return sha256(b).hexdigest()


def verify(ok,why):
    if not ok:
        raise ValueError(why)


def relative(path):
    return path.resolve().relative_to(ROOT).as_posix()


def read(path,expected,parse=True):
    p=(ROOT/path).resolve()
    name=relative(p)
    verify(not isinstance(path,str) or name==path,"noncanonical study-relative context path")
    raw=p.read_bytes()
    got=sha(raw)
    verify(got==expected==LINKS.get(name,got),f"changed {name}")
    LINKS[name]=got
    return json.loads(raw) if parse else raw


def truth(c):
    return isinstance(c,dict) and bool(c) and all(v is True for v in c.values())


def zero(v):
    return type(v) is int and v==0


def rational(v):
    verify(type(v) in (str,int,Q),"exact rational required")
    return Q(v)


def single(rows,key,value):
    r,=[x for x in rows if x[key]==value]
    return r


def decode_rf(rows):
    verify(type(rows) is list and len(rows)==N,"complete signed point RF required")
    values=[]
    for v in rows:
        lo,hi=(rational(v[k]) for k in ("lower_exact","upper_exact")) if type(v) is dict else map(rational,v)
        verify(lo<=hi,"reversed signed cached point")
        values.append((lo,hi))
    return values


def exact_budget(catalogue,z,rnorm):
    inc=Counter(i for e in catalogue for i in e["global_indices"])
    hinc=Counter()
    for c in catalogue:
        hinc.update({i:len(c["global_indices"]) for i in c["global_indices"]})
    verify(sorted(inc)==sorted(hinc)==list(range(N)),"complete 611 incidence")
    f=max(inc.values())
    h=max(hinc.values())
    reserve=1-z
    verify(f==6 and h==72 and 0<=z<1 and rnorm>0,"wrong physical incidence or headroom")
    cap_h=reserve/2
    cap_f=reserve*R/4
    local=(cap_f/(10000*rnorm*f),cap_h/(10000*rnorm*h))
    verify(cap_f+(z+cap_h)*R==R-cap_f<R and all(x>0 for x in local),
           "planned full derivative and force reserves do not leave strict margin")
    return f,h,local


def catalogue_match(e):
    tr=e["trial"]
    route=read(tr["routes"],tr["routes_sha256"])
    ret=[]
    for row in route["bulk"]:
        i=row["element"]
        verify(type(i) is int and row["chart"] in {"off_axis","axis_vertex","axis_edge_physical"},"malformed chart")
        ret.append({"kind":"bulk","slot":i,"nodes":row["nodes"],"global_indices":row["global_indices"],
                    "chart":row["chart"],"normalized_rule_sha256":row["auxiliary_element_rule_sha256"]})
    for row in route["surface"]:
        f,i=row["face"],row["edge"]
        verify(type(f) is type(i) is int and f in (0,1) and 0<=i<12,"malformed surface ordering")
        ret.append({"kind":"surface","slot":12*f+i,"nodes":row["nodes"],"global_indices":row["global_indices"],
                    "chart":"surface_physical_R","normalized_rule_sha256":row["trace_rule_sha256"]})
    verify(ret==e["component_catalogue"] and tuple((k["kind"],k["slot"]) for k in ret)==ORDER
           and all(all(type(i) is int and 0<=i<N for i in row["global_indices"])
                   and len(row["global_indices"])==len(set(row["global_indices"]))>0 for row in ret)
           and (e["point_force_scope"],e["whole_H_scope"])==("EXACT_CENTER_POINT_ONLY","ENTIRE_CLOSED_FIXED_PARAMETER_CUBE"),
           "current all168 signed scopes and projections do not reproduce")
    return ret


def independently_read_derivative(e,dr,raw):
    case=e["case_id"]
    tr=e["trial"]
    anchor=read(e["native_case"],e["native_case_sha256"])
    if e["cached_input_kind"]=="FRESH_R34":
        a=single(dr["case_results"],"case_id",case)
        loc=DR/(case+".json")
        verify(a["native_file"]==relative(loc) and a["status"]=="PASS_REVIEWED_FRESH_R34_CACHED_DERIVATIVE","bad actual fresh derivative leaf")
        leaf=read(loc,a["native_sha256"])
        producer=relative(D/"actual/cases"/case/"results.json")
        verify(a["producer_result_sha256"]==raw["sha256"][producer]==leaf["producer_result_sha256"],
               "fresh actual native and raw producer differ")
        matches={"trial_box_sha256":tr["trial_box_sha256"],"model_contract_sha256":tr["model_contract_sha256"],
                 "cache_packet_sha256":tr["point_cache_sha256"],"R_binary_sha256":tr["R_binary_sha256"],
                 "center_binary_sha256":tr["center_binary_sha256"],"proof_drive_exact":tr["exact_drive"],
                 "archived_center_shift_exact":tr["archived_center_shift_exact"],
                 "archived_xi_exact":e["original_addition_gap_exact"],"native_cached_anchor_sha256":e["native_case_sha256"],
                 "case_id":case,"record_id":e["record_id"]}
        verify(all(leaf[k]==v for k,v in matches.items()) and rational(leaf["radius_exact"])==R
               and leaf["site_counts"]=={"bulk":9216,"surface":288}
               and leaf["status"]==a["status"] and leaf["root_asserted"] is leaf["exact_spatial_derivative_claimed"] is False,
               "changed fresh cache, R, norm scope, point, radius or proof")
        z0,dz=map(rational,(leaf["native_center_Z0_upper_exact"],leaf["independent_variation_upper_exact"]))
        own=rational(leaf["independent_ZQ_upper_exact"])
        z=rational(leaf["reviewed_ZQ_upper_exact"])
        verify(own==z0+dz and z==max(own,rational(leaf["producer_ZQ_upper_exact"]))
               ==rational(a["reviewed_ZQ_upper_exact"]) and 0<=z0<=z<1,"invalid independent derivative arithmetic")
        source_r=rational(e["native_radius_exact"])
        if "own_radius_attempts" in anchor:
            old=single(anchor["own_radius_attempts"],"radius_exact",str(source_r))
        else:
            old=anchor
        prior=decode_rf(old["signed_center_force_intervals"])
        verify(prior==decode_rf(leaf["signed_native_cached_center_rows_exact"]),"signed native cached point anchor was changed")
        return z,rational(leaf["R_inf_exact"]),prior
    verify(e["cached_input_kind"]=="NATIVE_R34" and rational(e["native_radius_exact"])==R
           and tr["native_cached_bound_usable_at_trial_radius"] is True and anchor["cache_packet_sha256"]==tr["point_cache_sha256"]
           and anchor["field_sha256"]==e["original_observer_field_sha256"],
           "untouched native r34 premise required")
    if anchor["status"]=="PASS_INDEPENDENT_CACHED_LOCAL_ROOT":
        leaf=single(anchor["own_radius_attempts"],"radius_exact",str(R))
        tensors=anchor["per_tensor_sha256"]
        rows=decode_rf(leaf["signed_center_force_intervals"])
        value=rational(leaf["ZQ_upper_exact"])
        verify(leaf["whole_cube_AD_recomputed_for_this_radius"] is True,"old cached side is not closed-cube")
    else:
        verify(anchor["status"]=="PASS_INDEPENDENT_NEGATIVE_M_SMALLER_CACHED_ROOT","foreign cached premise")
        leaf=anchor
        old=read(anchor["reused_signed_cached_anchor_file"],anchor["reused_signed_cached_anchor_sha256"])
        tensors=old["cache_tensor_content_sha256"]
        rows=decode_rf(old["signed_center_force_intervals"])
        value=rational(leaf["independent_ZQ_upper_exact"])
    verify(tensors["u"]==tr["center_binary_sha256"] and tensors["preconditioner"]==tr["R_binary_sha256"]
           and value==rational(e["native_ZQ_exact"])<1,"mismatched native source center or derivative")
    return value,rational(leaf["preconditioner_infinity_norm_exact"]),rows


def reject(fn):
    try:
        fn()
    except (ValueError,KeyError):
        return True
    return False


def run():
    static={f:read(f,h,f.suffix==".json") for f,h in PINS.items()}
    p,m,pf=(static[P/(x+".json")] for x in ("protocol","manifest","preflight"))
    verify(p["status"]=="FROZEN_REMAINING32_EXACT_BRIDGE_COLD_PROTOCOL"
           and m["status"]=="COLD_CANONICAL_ALL_32_COMPONENT_CATALOGUES" and pf["status"]=="PASS"
           and truth(pf["checks"]) and p["spatial_worker_cli"] is False and m["no_spatial_workers_authorized"] is True
           and p["root_asserted"] is False and rational(p["closed_radius_exact"])==R
           and all(zero(o[k]) for o in (p,m,pf) for k in ("selector_components_prepared","new_spatial_integrations","cached_AD_jobs","solves"))
           and all(pf[k]==v for k,v in PHASE.items() if k!="preflight_sha256") and tuple(p["case_order"])==CASES,
           "bridge source cold binding or no-work statement mismatch")
    for name,h in p["static_input_sha256"].items():
        read(name,h,False)
    geo_pre=static[GP]
    var_pre=static[DP]
    context_review=static[TARGET/"independent_prelaunch/review.json"]
    derm,dert,dero=(static[DR/(x+".json")] for x in ("receipt","terminal","os_exit_observation"))
    raw=static[D/"actual/root_terminal_admission.json"]
    verify(dero["status"]=="PASS_OBSERVED_INDEPENDENT_REMAINING32_VARIATION_TRUE_OS_EXIT"
           and derm["status"]=="PASS_INDEPENDENT_ACTUAL_TEN_FRESH_CACHED_DERIVATIVES"
           and dert["status"]=="FINISHED" and zero(dero["observed_exit_code"]) and zero(dert["exit_code"])
           and type(dero["pid"]) is int and dero["pid"]==derm["pid"]==dert["pid"]
           and dero["receipt_sha256"]==dert["receipt_sha256"]==PINS[DR/"receipt.json"]
           and dero["terminal_sha256"]==PINS[DR/"terminal.json"]
           and derm["producer_admission_sha256"]==PINS[D/"actual/root_terminal_admission.json"]
           and all(derm["input_sha256"].get(k)==v for k,v in raw["sha256"].items())
           and len(raw["sha256"])==53 and derm["fresh_valid_derivative_count"]==10,
           "first native stage supply is not the complete reviewed variation10 phase")
    for name,h in raw["sha256"].items():
        read(name,h,False)
    ordered=tuple(x["case_id"] for x in m["cases"])
    verify(ordered==CASES and m["components_catalogued"]==5376 and m["case_count"]==32
           and [x["case_id"] for x in derm["case_results"]]==var_pre["case_order"]
           and len(set(var_pre["case_order"]))==10
           and tuple(geo_pre["case_order"])==CASES,"same complete canonical 32/10 case layout")
    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_remaining32_exact_bridge as api
    verify(Path(api.__file__).resolve()==SOURCE,"wrong scope adapter")
    cases=[]
    for e,g in zip(m["cases"],geo_pre["trials"],strict=True):
        case=e["case_id"]
        tr=e["trial"]
        assert g==tr and rational(tr["trial_radius_exact"])==R
        underlying=single(context_review["per_case"],"case_id",case)
        for k in ("box_sha256","routes_sha256","rule_sha256","model_contract_sha256","record_id"):
            verify(tr[k]==underlying[k],"previous independent cold context changed")
        if e["cached_input_kind"]=="FRESH_R34":
            verify(tr==single(var_pre["trials"],"case_id",case),"different derivative trial and geometry trial")
        else:
            verify(case not in var_pre["case_order"],"retained native scope mixed with fresh10")
        cat=catalogue_match(e)
        z,rnorm,sg=independently_read_derivative(e,derm,raw)
        iF,iH,caps=exact_budget(cat,z,rnorm)
        sizes=api.sizes(cat,z,rnorm)
        verify((sizes["c_F"],sizes["c_H"],tuple(map(rational,sizes["caps_exact"])))==(iF,iH,caps)
               and sizes["strictly_tail_only"] is True and sizes["quarter_ZQ_admission_rule"] is False,
               "independent incidence and reserve algebra differs from API")
        cases.append({"case_id":case,"record_id":e["record_id"],"artifact_identity":tr["model_contract_sha256"],
                      "trial_box_sha256":tr["trial_box_sha256"],"cached_input_kind":e["cached_input_kind"],
                      "observed_native_cached_signed_rows_count":len(sg),"ZQ_upper_exact":str(z),"R_inf_exact":str(rnorm),
                      "c_F":iF,"c_H":iH,"per_component_force_and_hessian_tail_caps_exact":list(map(str,caps))})
    checked=api.synthetic()
    verify(truth(checked),"producer-approved synthetic scope guards fail current replay")
    broken=[dict(row) for row in m["cases"][0]["component_catalogue"]]
    broken[0]["global_indices"]=[False,1]
    verify(reject(lambda:api.sizes(broken,Q(1,2),1)) and
           reject(lambda:api.sizes(m["cases"][0]["component_catalogue"][:-1],Q(1,2),1)),
           "bool alias or incomplete production catalogue accepted")
    blocked=api.sizes(m["cases"][0]["component_catalogue"],Q(1),Q(2))
    verify(blocked["positive_headroom"] is False and blocked["caps_exact"] is None,"insufficient contraction headroom still issued a budget")
    args=SimpleNamespace(review=HERE/"not_a_review.json",protocol_sha256=PHASE["protocol_sha256"],
         manifest_sha256=PHASE["manifest_sha256"],preflight_sha256=PHASE["preflight_sha256"],review_sha256="0"*64)
    verify(reject(lambda:api.gate(args)),"off-path independent gate accepted")
    payload={"status":"PASS_INDEPENDENT_REMAINING32_EXACT_BRIDGE_COLD",**PHASE,
             "checks":{"all_5376_signed_component_routes_equal_prior_independent_context":True,
                       "exact_32_trial_radii_and_no_case_switch":True,"fresh10_and_native22_signed_anchor_scopes_current":True,
                       "strict_tail_only_incidence_and_budget_arithmetic_for_all32":True,
                       "invalid_pid_scope_type_foreign_gate_and_partial_catalogue_rejected":True,
                       "same_Y_split_reference_and_pure_moment_firewall_source_read":True,
                       "separate_future_full_geometry_actual_input_and_per_case_selector_gates":True},
             "input_sha256":LINKS,"cold_review_source_sha256":sha(Path(__file__).read_bytes()),"case_rows":cases,
             "selector_only_authorized_in_this_gate":True,"spatial_workers_authorized_in_this_gate":False,
             "joined_full_geometry_and_derivative_authorized":False,
             "new_moment_cache_AD_or_solve_jobs":0,"new_actual_selector_components":0,
             "scope":"Cold bridge adapter, signed point-vs-whole scopes, projection and tail arithmetic. Actual joined 32 geometry authority still requires a separate receipt."}
    dst=HERE/"review.json"
    with dst.open("x",encoding="utf-8",newline="\n") as out:
        json.dump(payload,out,indent=2)
        out.write("\n")
    print(json.dumps({"review_sha256":sha(dst.read_bytes()),"case_count":len(cases),"input_count":len(LINKS)}))

if __name__=="__main__":
    run()

