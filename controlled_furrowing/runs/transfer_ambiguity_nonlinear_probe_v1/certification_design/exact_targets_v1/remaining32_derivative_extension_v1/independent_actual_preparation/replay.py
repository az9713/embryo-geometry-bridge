"""Independent Jet replay of ten cached derivative extensions using prior signed point anchors."""
import argparse
from base64 import b64decode
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import io
import json
import os
from pathlib import Path
import sys
import time
import traceback

import numpy as np

HERE=Path(__file__).resolve().parent
D=HERE.parent
EX=D.parent
ROOT=next(p for p in EX.parents if (p/"PLAN.md").is_file())
SOURCE=Path(__file__).resolve()
ACTUAL=D/"actual"
T=EX/"remaining32_targets_v1"
G=EX/"remaining32_geometry_v1"
AD=ROOT/"scripts/review_algebraic_root.py"
LAUNCHER=ROOT/"scripts/transfer_ambiguity_remaining32_geometry_variation.py"
PINS={
    LAUNCHER:"68bca869ee486d7dc25b34cae3ceee1ddc2adc53500a20301c0038c259a2b26a",
    D/"preparation/protocol.json":"6d04fc82779430e66c2d27b5ab7aa508f522f7ffac03b87f8f2e8e2587e1d128",
    D/"preparation/preflight.json":"abbbeb7c0c303db2d2675398f8ea09155a19fd969b121e7936d7f4a1a468e264",
    G/"preparation/protocol.json":"754da64bb96cd8a19707181315e6c3e3f014e8d37a898b83c1693342c82c9850",
    G/"preparation/preflight.json":"d70a7217c4dbbf68f789bb59193321252ea4192576c3f41790cc7daea38fb15a",
    G/"independent_prelaunch/review.json":"0c30f413a4796703f6f4a8248924587fcb2152a2dfe7449464fb69c99d763511",
    T/"preparation/manifest.json":"dfea134640f99b45e9ca75d1ce718043d67f2917375506997d2465abda8c1242",
    T/"independent_prelaunch/review.json":"67d33bbc3171839b1c578b02b5cfbf9ea6ba5d6cb351d0655dd14aaf99df75f9",
    AD:"6af2311d0667a28af263b6fb6237063481aa8e65209ae9722c2b0246a8cff27a",
}
R36,R34,R32=(Q(1,2**n) for n in (36,34,32))
N=611
BITS=320
SITE={"bulk":9216,"surface":288}
MARGINS={"bulk_det_min","bulk_hoop_min","bulk_J_min","surface_radius_min","surface_tangent_squared_min",
         "surface_J_min","normalized_weight_min"}
RECENTERED={"max_contrast_change_small_13","max_negative_contrast_change_small_16"}
FRESH="FRESH_VARIATION_REQUIRED_BEFORE_SPATIAL_LAUNCH"
SMALL="EXISTING_MATCHED_SMALL_CUBE_LARGER_TRIAL_REQUIRES_NEW_AD"
STAGE="variation"
SCOPE=("Ten fixed-input equivalent exact-real Q8/S12 cached derivatives at whole r=2^-34, "
       "with the identical independently certified center/point force/center derivative and frozen binary R. "
       "No spatial force or derivative, body embedding, solve, stability or source result.")


def need(ok,why):
    if not ok:
        raise ValueError(why)


def canonical(path):
    p=(ROOT/path).resolve()
    name=p.relative_to(ROOT).as_posix()
    if isinstance(path,str):
        need(path==name,"noncanonical study-relative artifact path")
    return p,name


def sha(raw):
    return sha256(raw).hexdigest()


def read(path,digest,links,parse=True):
    p,name=canonical(path)
    raw=p.read_bytes()
    h=sha(raw)
    need(h==digest==links.get(name,h),f"changed bytes {name}")
    links[name]=h
    return json.loads(raw) if parse else raw


def save(path,value):
    with path.open("x",encoding="utf-8",newline="\n") as out:
        json.dump(value,out,indent=2)
        out.write("\n")
    return sha(path.read_bytes())


def q(value):
    need(type(value) in (str,int,Q),"exact rational input only")
    return Q(value)


def interval(row):
    a,b=(q(row[k]) for k in ("lower_exact","upper_exact"))
    need(a<=b,"invalid signed pair")
    return a,b


def zero(x):
    return type(x) is int and x==0


def positive(x):
    return type(x) is int and x>0


def yes(conditions):
    return type(conditions) is dict and bool(conditions) and all(v is True for v in conditions.values())


def select(rows,name):
    row,=[r for r in rows if r["case_id"]==name]
    return row


def group_hash(order):
    return sha(json.dumps(order,separators=(",",":"),sort_keys=True).encode())


def weighted(abs_rows, columns):
    need(all(len(row)==len(columns) for row in abs_rows) and all(v>=0 for v in columns)
         and all(v>=0 for row in abs_rows for v in row),"typed nonnegative complete majorant")
    return [sum((v*c for v,c in zip(row,columns,strict=True)),Q(0)) for row in abs_rows]


def reserve(z,r):
    need(z>=0 and r>0,"positive tested radius and nonnegative derivative")
    if z>=1:
        return None
    h=1-z
    b=h/2
    f=h*r/4
    need(f+(z+b)*r==r-f<r,"exact reserved planning inequality")
    return {"headroom_exact":str(h),"b_cap_exact":str(b),"signed_force_cap_exact":str(f)}


def old_box_and_trial(e,tr,links):
    src=read(e["box"],e["box_sha256"],links)
    dest=read(tr["trial_box"],tr["trial_box_sha256"],links)
    permitted={"coordinate_radii_exact","coordinate_hulls_exact","corrector_radius_exact"}
    need(set(src)==set(dest) and all(src[k]==v for k,v in dest.items() if k not in permitted)
         and src["root_membership_assumed"] is dest["root_membership_assumed"] is False,
         "trial moved relative to old admitted exact center")
    r=q(src["corrector_radius_exact"])
    need(r==q(e["planned_geometry_radius_exact"]) and q(dest["corrector_radius_exact"])==R34,
         "same archived and frozen r34 trial radii")
    if r==R36:
        need(tr["trial_box_sha256"]!=e["box_sha256"] and tr["derivative_policy"]==SMALL,
             "old derivative bound silently enlarged")
    else:
        need(r==R34 and tr["trial_box_sha256"]==e["box_sha256"] and tr["derivative_policy"]==FRESH,
             "mixed prepared trial")
    for b, radius in ((src,r),(dest,R34)):
        free=b["free"]
        need(len(free)==len(set(free))==N and all(type(i) is int and 0<=i<650 for i in free),
             "invalid ordered projection")
        free=set(free)
        need(len(b["current_centers_exact"])==325,"whole 650-coordinate trial required")
        for i,(Y,de,ends) in enumerate(zip(b["current_centers_exact"],b["coordinate_radii_exact"],
                                          b["coordinate_hulls_exact"],strict=True)):
            need(len(Y)==len(de)==len(ends)==2,"wrong coordinate dimension")
            for j,(v,d,h) in enumerate(zip(map(q,Y),map(q,de),ends,strict=True)):
                need(d==(radius if 2*i+j in free else 0) and list(map(q,h))==[v-d,v+d],
                     "missing full matched trial hull")
    return dest


def cold(links):
    docs={p:read(p,h,links,p.suffix==".json") for p,h in PINS.items()}
    p=docs[D/"preparation/protocol.json"]
    m=docs[T/"preparation/manifest.json"]
    top_review=docs[T/"independent_prelaunch/review.json"]
    configured=docs[G/"independent_prelaunch/review.json"]
    expected={"source_sha256":PINS[LAUNCHER],
              "target_manifest_sha256":PINS[T/"preparation/manifest.json"],
              "target_context_review_sha256":PINS[T/"independent_prelaunch/review.json"],
              **{prefix+"_"+n+"_sha256":PINS[base/"preparation"/(n+".json")]
                 for prefix,base in (("geometry",G),("variation",D)) for n in ("protocol","preflight")}}
    gate=docs[D/"preparation/preflight.json"]
    need(p["stage"]==gate["stage"]==STAGE and gate["status"]=="PASS" and yes(gate["checks"])
         and gate["source_sha256"]==p["source_sha256"]==PINS[LAUNCHER]
         and gate["protocol_sha256"]==expected["variation_protocol_sha256"]
         and configured["status"]=="PASS_INDEPENDENT_REMAINING32_GEOMETRY_VARIATION_PRELAUNCH"
         and yes(configured["checks"]) and all(configured[k]==v for k,v in expected.items())
         and configured["exact_spatial_bridge_authorized"] is False
         and top_review["status"]=="PASS_INDEPENDENT_REMAINING32_TARGET_CONTEXTS" and yes(top_review["checks"])
         and top_review["manifest_sha256"]==expected["target_manifest_sha256"],
         "current cold stage and separately pinned prelaunch gates required")
    for name,h in p["kernel_source_sha256"].items():
        read(name,h,links,False)
    order=p["case_order"]
    need(len(order)==len(set(order))==10 and order==[r["case_id"] for r in p["trials"]]
         ==configured["permitted_case_ids"][STAGE]
         and Counter(t["derivative_policy"] for t in p["trials"])=={FRESH:7,SMALL:3}
         and all(t["native_cached_bound_usable_at_trial_radius"] is False for t in p["trials"])
         and {t["case_id"] for t in p["trials"] if q(t["archived_center_shift_exact"])>0}==RECENTERED,
         "must replay exactly the ten fresh cases only")
    for tr in p["trials"]:
        e=select(m["targets"],tr["case_id"])
        saved=select(top_review["per_case"],tr["case_id"])
        need(all(tr[k]==v for k,v in e.items()) and all(tr[k]==saved[k] for k in (
            "model_contract_sha256","box_sha256","rule_sha256","routes_sha256")),
            "cold same-case independent review not retained")
        old_box_and_trial(e,tr,links)
    return p,m,top_review,{**expected,"launcher_review_sha256":PINS[G/"independent_prelaunch/review.json"]}


def admission(producer,expected,digest,links):
    admit=read(ACTUAL/"root_terminal_admission.json",digest,links)
    order=producer["case_order"]
    g=group_hash(order)
    need(admit["status"]=="PASS_ROOT_OBSERVED_REMAINING32_VARIATION_TERMINALS"
         and positive(admit["pid"]) and positive(admit["unified_exec_session"])
         and zero(admit["observed_exit_code"]) and admit["launch_group"]==g,
         "genuine external OS exit0 for this full producer variation group required")
    group=ACTUAL/"launches"/g
    parent={n:group/(n+".json") for n in ("started","terminal","manifest")}
    fields=("launched","started","terminal","results","observed_os_exit")
    files=admit["sha256"]
    names={canonical(path)[1] for path in parent.values()}|{
        canonical(ACTUAL/"cases"/case/(f+".json"))[1] for case in order for f in fields}
    need(set(files)==names and len(files)==53,"incomplete or extra externally frozen producer terminal closure")
    parent_docs={k:read(v,files[canonical(v)[1]],links) for k,v in parent.items()}
    start,end,result=(parent_docs[k] for k in ("started","terminal","manifest"))
    identity={"pid":admit["pid"],"bindings":expected,"stage":STAGE,"launch_group":g,"case_order":order}
    need(all(all(x[k]==v for k,v in identity.items()) for x in (start,end,result))
         and start["status"]=="RUNNING" and end["status"]=="FINISHED" and zero(end["exit_code"])
         and end["manifest_sha256"]==files[canonical(parent["manifest"])[1]]
         and result["status"]=="COMPLETE_ACTUAL_WORK_AWAITS_INDEPENDENT_REPLAY" and result["root_asserted"] is False
         and end["workers"]==result["workers"] and [x["case_id"] for x in end["workers"]]==order
         and len({x["pid"] for x in end["workers"]})==10,
         "producer group manifests, observed status or frozen sources differ")
    proof_hashes={}
    for w in end["workers"]:
        case=w["case_id"]
        folder=ACTUAL/"cases"/case
        assert set(fields)=={"launched","started","terminal","results","observed_os_exit"}
        sub={k:folder/(k+".json") for k in fields}
        hashes={k:files[canonical(v)[1]] for k,v in sub.items()}
        closed={k:read(sub[k],hashes[k],links) for k in fields if k!="results"}
        start0,launch0,end0,obs=(closed[k] for k in ("started","launched","terminal","observed_os_exit"))
        need(positive(w["pid"]) and {k:v for k,v in w.items() if k not in (
            "os_exit_sha256","results_sha256","candidate_status")}==obs and zero(obs["observed_exit_code"])
            and w["os_exit_sha256"]==hashes["observed_os_exit"]
            and obs["terminal_sha256"]==hashes["terminal"] and end0["results_sha256"]==hashes["results"]==w["results_sha256"]
            and end0["status"]=="FINISHED" and zero(end0["exit_code"]) and start0["status"]=="RUNNING"
            and all(all(x[k]==v for k,v in {"case_id":case,"pid":w["pid"],"bindings":expected,
                    "supervisor_pid":admit["pid"],"stage":STAGE,"launch_group":g}.items())
                    for x in (obs,launch0,start0,end0)), "worker is not a fully matched observed terminal")
        proof_hashes[case]=(sub["results"],hashes["results"],w)
    # No numeric payload is opened until every child has been admitted.
    return admit,proof_hashes


def prior_point(contract,protocol_review,links):
    native=read(contract["native_case"],contract["native_case_sha256"],links)
    r=q(contract["source_cached_radius_exact"])
    if "own_radius_attempts" in native:
        source,=[leaf for leaf in native["own_radius_attempts"] if q(leaf["radius_exact"])==r]
        need(q(native["chosen_radius_exact"])==r and source["whole_cube_AD_recomputed_for_this_radius"] is True,
             "matching independently passing native cached radius required")
    else:
        source=native
        need(q(source["radius_exact"])==r,"native cached radius drift")
    key=native.get("case_id",native.get("pilot_id"))
    need(key==contract["pilot_id"] and native["record_id"]==contract["physical_source_record_id"],
         "prior point provenance is not the current exact case")
    eta=q(source.get("eta_upper_exact",source.get("independent_eta_upper_exact")))
    Z0=q(source.get("center_Z0_upper_exact",source.get("independent_Z0_upper_exact")))
    z=q(source.get("ZQ_upper_exact",source.get("independent_ZQ_upper_exact")))
    err=q(source.get("posterior_error_upper_exact",source.get("independent_posterior_error_upper_exact")))
    need(z==q(contract["native_cached_ZQ_upper_exact"]) and 0<=Z0<=z<1
         and 0<=eta and eta+z*r<r and err==eta/(1-z)<r
         and r==q(protocol_review["native_cached_radius_exact"])
         and err==q(protocol_review["native_cached_posterior_upper_exact"]),
         "retained point and exact native whole-cube premises are not coherent")
    if key in RECENTERED:
        need(native["recentered_cache_packet_sha256"]==contract["point_cache_sha256"]
             and native["original_binary_R_sha256"]==contract["binary_R_sha256"]
             and native["replayed_u1_binary_sha256"]==contract["center_binary_sha256"]
             and native["archived_field_sha256"]==contract["archived_field_sha256"]
             and q(native["actual_binary_center_shift_sigma_exact"])==q(contract["archived_shift_exact"])>0,
             "recentered derivative must use new u1, same R and old archived observer")
    else:
        need(native["cache_packet_sha256"]==contract["point_cache_sha256"]
             and native["field_sha256"]==contract["archived_field_sha256"] and q(contract["archived_shift_exact"])==0,
             "cached point bounds for different archive or operator")
    rows=source.get("signed_center_force_intervals")
    need(type(rows) is list and len(rows)==N and max(max(map(abs,interval(r))) for r in rows)<=eta
         and source.get("precision_bits",256)>=256,
         "full independent signed cached anchor required")
    return native,eta,Z0,[interval(r) for r in rows]


def frozen_arrays(tr,ctx,packet,box):
    payload=packet["frozen_binary_cache"]
    raw=b64decode(payload["data"],validate=True)
    need(sha(raw)==payload["sha256"]==ctx["frozen_binary_cache_sha256"],"embedded original/current cache hash drift")
    with np.load(io.BytesIO(raw),allow_pickle=False) as z:
        a={k:z[k] for k in z.files}
    need({k:list(v.shape) for k,v in a.items()}==payload["array_shapes"]
         and a["preconditioner"].dtype.str==a["u"].dtype.str=="<f8"
         and a["preconditioner"].shape==(N,N) and a["u"].shape==(N,)
         and a["reference"].shape==(325,2) and a["triangles"].shape==(144,6)
         and a["shape"].shape==(64,6) and a["shape_grad"].shape==(144,64,6,2)
         and a["free"].tolist()==box["free"] and a["triangles"].tolist()==box["triangles"]
         and sha(a["u"].tobytes(order="C"))==tr["center_binary_sha256"]==ctx["center_binary_sha256"]
         and sha(a["preconditioner"].tobytes(order="C"))==tr["R_binary_sha256"]==ctx["binary_R_sha256"]
         and Q(float(a["energy_scale"][0]))==q(ctx["energy_scale_binary_exact"])==q(tr["energy_scale_binary_exact"])>0,
         "explicit frozen owner cache changed")
    x=[[Q(float(v)) for v in row] for row in a["reference"]]
    inc=dict(zip(box["free"],(Q(float(v)) for v in a["u"]),strict=True))
    same=[[str(v+inc.get(2*i+j,Q(0))) for j,v in enumerate(row)] for i,row in enumerate(x)]
    need(same==box["current_centers_exact"] and [[str(v) for v in row] for row in x]==box["reference_nodes_exact"],
         "new derivative arrays do not reconstruct the exact enlarged trial box center")
    for v in a.values():
        v.flags.writeable=False
    return a


def recalc(ad, arrays,params,drive):
    ad.ctx.prec=BITS
    ad.A.cache_clear()
    cache=[]
    n=Counter()
    site_functions={"bulk":ad.bulk_AD,"surface":ad.surface_AD}
    for t in ad.terms_from_arrays(arrays,params,drive):
        need(t["weight"]>0 and t["indices"],"strictly positive projected cached term required")
        _,H,_,_=site_functions[t["kind"]](t["q"],*t["parameters"])
        t["Hcenter"]=H
        t["absB"]=[[ad.UP(v) for v in row] for row in t["B"]]
        t["rowsumB"]=[sum(row,ad.arb(0)).upper() for row in t["absB"]]
        n[t["kind"]]+=1
        cache.append(t)
    need(dict(n)==SITE and len(cache)==9504,"replay did not reconstruct all cached sites exactly")
    vr,domain=ad.variation(cache,N,ad.A(R34))
    need(domain["counts"]==SITE and set(domain["strict_lower_margins"])==MARGINS
         and all(interval(v)[0]>0 for v in domain["strict_lower_margins"].values())
         and interval(domain["strict_lower_margins"]["bulk_J_min"])[0]>Q(1,1000000),
         "failed independent whole cached analytic domain")
    return [Q(str(vr[i,0].upper().fmpq())) for i in range(N)],domain


def entry(ad,tr,e,top,source,links):
    name=tr["case_id"]
    path,ph,worker=source
    producer=read(path,ph,links)
    need(all(producer[k]==worker[k] for k in ("case_id","pid","supervisor_pid","stage","launch_group","bindings"))
         and producer["trial"]==tr and producer["status"]==worker["candidate_status"]
         and producer["root_asserted"] is False and producer["exact_spatial_derivatives_computed"] is False,
         "producer variation result and observed actual child do not match")
    if producer["status"]!="FRESH_CACHED_VARIATION_AWAITS_INDEPENDENT_REVIEW":
        return {"status":"UNRESOLVED_PRODUCER_ANALYTIC_DOMAIN","case_id":name,"producer_result_sha256":ph,
                "archive_reason":producer.get("reason"),"trial_box_sha256":tr["trial_box_sha256"]}
    box=old_box_and_trial(e,tr,links)
    ctx=read(e["model_contract"],e["model_contract_sha256"],links)
    need(ctx["pilot_id"]==name and ctx["record_id"]==e["record_id"] and ctx["point_cache_sha256"]==tr["point_cache_sha256"]
         and ctx["archived_shift_exact"]==tr["archived_center_shift_exact"]==top["observer_shift_sigma_exact"],
         "trial and its independent exact context authority differ")
    pk=read(ctx["point_cache"],ctx["point_cache_sha256"],links)
    a=frozen_arrays(tr,ctx,pk,box)
    native,eta,z0,prior=prior_point(ctx,top,links)
    activation=q(tr["exact_drive"])
    strength=ctx["exact_requested_drive"]
    need(activation==q(strength["activation_strength_exact"])==q(strength["optical_dose_exact"])*q(strength["common_gain_exact"])
         and activation!=q(strength["activation_strength_binary_exact"])
         and activation==q(pk["activation_strength_exact"])==q(native.get("exact_requested_strength",native.get("exact_strength"))),
         "requested exact physical drive replaced by rounded solver activation")
    rows=[[abs(Q(v)) for v in row] for row in a["preconditioner"].tolist()]
    norm=max(sum(row) for row in rows)
    need(norm==q(producer["R_inf_exact"])>0 and producer["site_counts"]==SITE,"producer bound from mixed R or cache")
    f=[interval(t) for t in producer["signed_RFQ_center_rows"]]
    pe=max(max(map(abs,z)) for z in f)
    nrows=list(map(q,producer["center_defect_row_bounds_exact"]))
    need(len(f)==len(nrows)==N and pe==q(producer["eta_Q_upper_exact"]) and min(nrows)>=0
         and max(nrows)==q(producer["center_Z0_upper_exact"])
         and all(max(a[0],b[0])<=min(a[1],b[1]) for a,b in zip(f,prior,strict=True)),
         "independent prior signed anchor differs from producer point")
    pv=list(map(q,producer["unmapped_variation_rows_exact"]))
    pm=list(map(q,producer["preconditioned_variation_rows_exact"]))
    exact_p=weighted(rows,pv)
    pz=max(nrows)+max(pm)
    need(len(pm)==N and all(a<=b for a,b in zip(exact_p,pm,strict=True))
         and max(pm)==q(producer["variation_upper_exact"]) and pz==q(producer["ZQ_upper_exact"]),
         "producer cached derivative arithmetic does not reproduce")
    if producer["whole_analytic_domain"]["term_counts"]!=SITE:
        raise ValueError("producer whole-cube term counts changed")
    try:
        v,dom=recalc(ad,a,pk["parameters"],activation)
    except ValueError as error:
        return {"status":"UNRESOLVED_INDEPENDENT_WHOLE_CUBE_DOMAIN","case_id":name,
                "producer_result_sha256":ph,"trial_box_sha256":tr["trial_box_sha256"],"reason":str(error)}
    mapped=weighted(rows,v)
    own_z=z0+max(mapped)
    combined=max(own_z,pz)
    past_sigma=q(ctx["archived_shift_exact"])
    lhs=eta+own_z*R34
    selfmap=own_z<1 and lhs<R34
    err=eta/(1-own_z) if selfmap else None
    return {"status":"PASS_REVIEWED_FRESH_R34_CACHED_DERIVATIVE","case_id":name,"record_id":ctx["record_id"],
            "producer_result_sha256":ph,"trial_box_sha256":tr["trial_box_sha256"],
            "model_contract_sha256":e["model_contract_sha256"],"cache_packet_sha256":ctx["point_cache_sha256"],
            "native_cached_anchor_sha256":ctx["native_case_sha256"],"R_binary_sha256":ctx["binary_R_sha256"],
            "center_binary_sha256":ctx["center_binary_sha256"],"proof_drive_exact":str(activation),"radius_exact":str(R34),
            "archived_center_shift_exact":str(past_sigma),"archived_xi_exact":ctx["archived_addition_gap_exact"],
            "R_inf_exact":str(norm),"site_counts":SITE,"analytic_domain":dom,"precision_bits":BITS,
            "signed_native_cached_center_rows_exact":[[str(a),str(b)] for a,b in prior],
            "native_eta_upper_exact":str(eta),"native_center_Z0_upper_exact":str(z0),
            "independent_unmapped_rows_exact":list(map(str,v)),"independent_preconditioned_rows_exact":list(map(str,mapped)),
            "independent_variation_upper_exact":str(max(mapped)),
            "independent_ZQ_upper_exact":str(own_z),"producer_ZQ_upper_exact":str(pz),
            "reviewed_ZQ_upper_exact":str(combined),"bridge_budget_if_positive_headroom":reserve(combined,R34),
            "diagnostic_cached_r34_selfmap_lhs_upper_exact":str(lhs),"diagnostic_cached_r34_selfmap_pass":selfmap,
            "diagnostic_cached_r34_posterior_upper_exact":str(err) if err is not None else None,
            "diagnostic_archived_shift_plus_error_below_observer_goal":selfmap and past_sigma+err<q(ctx["observer_goal_exact"]),
            "root_asserted":False,"exact_spatial_derivative_claimed":False}


def prepare():
    links={}
    p,_,_,expected=cold(links)
    exact=weighted([[Q(3),Q(1,2)],[Q(1,4),Q(2)]],[Q(1,8),Q(1,32)])
    tests={"all_ten_same_center_only_seven_shrink_three_enlarge":True,
           "exact_R_positive_row_weighted_sums":exact==[Q(25,64),Q(3,32)],
           "headroom_strict_and_no_positive_zero_slack_claim":reserve(Q(1,3),R34) is not None and reserve(Q(1),R34) is None,
           "zero_exit_type_rejects_boolean":not zero(False) and zero(0),
           "no_live_producer_packets_or_AD_in_cold_phase":not any("actual/cases" in f for f in links)}
    need(yes(tests),"cold independent fixtures failed")
    h=sha(SOURCE.read_bytes())
    protocol=save(HERE/"protocol.json",{"status":"FROZEN_INDEPENDENT_TEN_REMAINING32_VARIATION_REPLAY","source_sha256":h,
                   "stage_input_sha256":links,"case_order":p["case_order"],"producer_bindings":expected,
                   "fresh_cases":10,"old_matched_not_replayed":22,"precision_bits":BITS,"radius_exact":str(R34),
                   "site_counts":SITE,"scope":SCOPE,"cold_actual_payloads_opened":0,"cold_AD_jobs":0})
    pf=save(HERE/"preflight.json",{"status":"PASS","checks":tests,"source_sha256":h,"protocol_sha256":protocol,
                                 "cold_actual_payloads_opened":0,"cold_AD_jobs":0})
    print(json.dumps({"source_sha256":h,"protocol_sha256":protocol,"preflight_sha256":pf}))


def run(args):
    linked={}
    h=sha(SOURCE.read_bytes())
    p=read(HERE/"protocol.json",args.protocol_sha256,linked)
    pf=read(HERE/"preflight.json",args.preflight_sha256,linked)
    approved=read(HERE/"root_source_review.json",args.source_review_sha256,linked)
    phase={"source_sha256":h,"protocol_sha256":args.protocol_sha256,"preflight_sha256":args.preflight_sha256}
    need(pf["status"]=="PASS" and approved["status"]=="PASS" and yes(pf["checks"]) and yes(approved["checks"])
         and all(approved[k]==v for k,v in phase.items())
         and p["source_sha256"]==pf["source_sha256"]==h and pf["protocol_sha256"]==args.protocol_sha256,
         "independent frozen source/math prelaunch missing")
    current={}
    producer,m,top,expected=cold(current)
    need(p["stage_input_sha256"]==current and p["producer_bindings"]==expected
         and p["case_order"]==producer["case_order"],"cold phase changed")
    linked.update(current)
    a,rows=admission(producer,expected,args.admission_sha256,linked)
    folder=HERE/"actual_replay"
    folder.mkdir(exist_ok=False)
    phase.update(pid=os.getpid(),producer_admission_sha256=args.admission_sha256,
                 root_source_review_sha256=args.source_review_sha256)
    save(folder/"started.json",{"status":"RUNNING",**phase})
    result={"status":"FAILED","exit_code":1,**phase}
    start=time.perf_counter()
    try:
        sys.path.insert(0,str(ROOT/"scripts"))
        import review_algebraic_root as ad
        need(Path(ad.__file__).resolve()==AD,"changed independent Jet runtime module")
        summary=[]
        for t in producer["trials"]:
            name=t["case_id"]
            print(json.dumps({"phase":"independent_variation_started","case_id":name}),flush=True)
            record=entry(ad,t,select(m["targets"],name),select(top["per_case"],name),rows[name],linked)
            path=folder/(name+".json")
            hash_=save(path,{**record,"phase":phase})
            summary.append({"case_id":name,"native_file":canonical(path)[1],"native_sha256":hash_,**{
                k:record.get(k) for k in ("status","reviewed_ZQ_upper_exact","independent_ZQ_upper_exact","producer_ZQ_upper_exact",
                                          "producer_result_sha256","trial_box_sha256","cache_packet_sha256","center_binary_sha256",
                                          "R_binary_sha256","proof_drive_exact","radius_exact","bridge_budget_if_positive_headroom")}})
            print(json.dumps(summary[-1]),flush=True)
        complete=sum(r["status"]=="PASS_REVIEWED_FRESH_R34_CACHED_DERIVATIVE" for r in summary)
        receipt={"status":"PASS_INDEPENDENT_ACTUAL_TEN_FRESH_CACHED_DERIVATIVES" if complete==10
                   else "INDEPENDENT_VARIATION_REPLAY_WITH_UNRESOLVED_CASES",
                 **phase,"scope":SCOPE,"case_results":summary,"fresh_valid_derivative_count":complete,
                 "root_observed_producer_pid":a["pid"],"root_observed_producer_exec":a["unified_exec_session"],
                 "old_native22_included_in_fresh_authority":False,
                 "spatial_derivative_root_or_integration_count":0,"independent_terminal_observation_required":True,
                 "input_sha256":linked,"seconds":time.perf_counter()-start}
        rh=save(folder/"receipt.json",receipt)
        result.update(status="FINISHED",exit_code=0,receipt_sha256=rh)
        print(json.dumps({"status":receipt["status"],"receipt_sha256":rh}),flush=True)
    except BaseException as error:
        save(folder/"failure.json",{"error":str(error),"traceback":traceback.format_exc(),**phase})
        raise
    finally:
        save(folder/"terminal.json",result)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    sub=cli.add_subparsers(dest="phase",required=True)
    sub.add_parser("prepare")
    act=sub.add_parser("run")
    for flag in ("admission-sha256","protocol-sha256","preflight-sha256","source-review-sha256"):
        act.add_argument("--"+flag,required=True)
    args=cli.parse_args()
    if args.phase=="prepare":
        prepare()
    else:
        run(args)

