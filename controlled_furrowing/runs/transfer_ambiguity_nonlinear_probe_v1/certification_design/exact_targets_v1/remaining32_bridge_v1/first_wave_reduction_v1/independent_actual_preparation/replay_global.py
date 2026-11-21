"""Independent first-wave global endpoint replay, cold preparation only by default."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement
import os
from pathlib import Path
import sys
import traceback
from types import SimpleNamespace
from unittest.mock import patch

from flint import arb, arb_mat, ctx, fmpq

SOURCE=Path(__file__).resolve()
HERE=SOURCE.parent
ROOT=HERE.parents[6]
sys.path.insert(0,str(ROOT/"scripts"))
import transfer_ambiguity_remaining32_first_wave_reduce as D

S,W,V=D.S,D.W,D.V
N=611
BITS=768
OUT=HERE/"actual_replay"
PRODUCER_OS=D.ACTUAL/"root_terminal_admission.json"
SOURCE_REVIEW=HERE/"root_source_review.json"
RAW_SHA="bb9f059699e6d98fe626491f252edddb1664fbf83d79cb74258084ffd399c880"
PINS={
 D.SOURCE:"58de9f5e96b43005b0ba83dd0f4d8cd1b35f0690eef58290a75ef5abce692d92",
 D.PREP/"protocol.json":"4b9ecc9794f84e00d62f1fbc090b3a96739af044a2a426194d8fba07967fe487",
 D.PREP/"preflight.json":"2f1c69917e32d25d22a28577f80dc6a2c8c7fab025d93967e73fcf48cffeb79c",
 D.REVIEW_FILE:"7040ba3e30a276df43e41a61378b86f3223633b47fd3b7c5dce810e5abf58127",
 D.RAW:RAW_SHA,
 D.SPD/"actual_replay/os_exit_observation.json":"fd8ac1cb3b920b57096649544bba105442404fe19576aae99c53f849fe15a8d5",
}
OS_STATUS="PASS_ROOT_OBSERVED_FIRST_WAVE_SIGNED_REDUCTION_TERMINAL"
REVIEW="PASS_ROOT_FIRST_WAVE_INDEPENDENT_GLOBAL_ENDPOINT_SOURCE"
SCOPE="Independent exact endpoint scatter and outward positive/negative R split for all four spatial singleton hull targets."
ZERO={"actual_hull_leaves_opened":0,"new_spatial_integrals_or_AD_jobs":0,"actual_global_reductions":0}
read,save,sha,require=D.read,D.save,D.sha,D.require
rel=D.rel


def q(x):
    require(type(x) in (str,int,Fraction),"exact signed rational")
    return Fraction(x)


def pair(x):
    require(type(x) in (list,tuple) and len(x)==2,"full endpoint pair")
    ans=tuple(map(q,x))
    require(ans[0]<=ans[1],"ordered exact endpoints")
    return ans


def add(x,y):
    return x[0]+y[0],x[1]+y[1]


def pack(row):
    return list(map(str,row))


def static():
    require(all(sha(path)==h for path,h in PINS.items()),"static reviewed producer or observed metadata changed")
    a,w=D.static_phase()
    return {"status":"FROZEN_FIRST_WAVE_INDEPENDENT_GLOBAL_ENDPOINT_REPLAY","scope":SCOPE,
            "source_sha256":sha(SOURCE),"static_pins":{rel(p):h for p,h in PINS.items()},
            "worker_binding":w,"case_order":list(W.FIRST),"precision_bits":BITS,**ZERO,
            "recorded_spatial_producer_admission_sha256":RAW_SHA,
            "point_and_whole_roles":[S.POINT,S.BOX],"root_no_heldout_or_continuum_claim":True}


def validated_channel(leaf,name,indices,raw):
    channel=leaf["signed_channels"][name]
    pairs=[]
    misses=[]
    require(len(channel["rows"])==len(indices)==len(raw),"complete local signed channel")
    for x,idx,producer in zip(channel["rows"],indices,raw,strict=True):
        i=x["index"]
        require((type(i) is int and type(idx) is int
                or type(i) is list and type(idx) is list and all(type(k) is int for k in i))
                and i==idx,"typed signed channel row")
        a,b=pair(x["independent_exact"])
        c,d=pair(x["producer_exact"])
        e,f=pair(x["hull_exact"])
        separate=b<c or d<a
        contained=c<=a and b<=d
        require((c,d)==pair(producer) and not separate and (e,f)==(min(a,c),max(b,d))
                and x["producer_contains_independent"] is contained and x["disjoint"] is False,
                "independent valid enclosure full exact union and no disjoint required")
        if not contained:
            misses.append(idx)
        pairs.append((e,f))
    require(channel["disjoint_indices"]==[] and channel["producer_containment_misses"]==misses,
            "retained true nonnesting diagnostics")
    return pairs,len(misses)


def scatter_add(force,upper,ids,F,H):
    require(len(ids)==len(F) and len(H)==len(ids)*(len(ids)+1)//2 and
            len(set(ids))==len(ids) and all(type(i) is int and 0<=i<len(force) for i in ids),
            "complete typed local projection")
    for index,value in zip(ids,F,strict=True):
        force[index]=add(force[index],value)
    for (i,j),value in zip(combinations_with_replacement(ids,2),H,strict=True):
        label=min(i,j),max(i,j)
        upper[label]=add(upper.get(label,(Fraction(0),Fraction(0))),value)


def bound(r,force,H):
    """Independent tight endpoint algebra with a fixed signed R and box-valued entries."""
    n=len(r)
    require(n>0 and len(force)==n and all(len(row)==n for row in r),"exact fixed square R and force")
    require(all(type(i) is type(j) is int and 0<=i<=j<n for i,j in H),"typed symmetric global H map")
    rr=[[q(v) for v in row] for row in r]
    def a(x):
        return arb(fmpq(x.numerator,x.denominator))
    ctx.prec=BITS
    pos=arb_mat(n,n,[a(max(0,v)) for row in rr for v in row])
    neg=arb_mat(n,n,[a(max(0,-v)) for row in rr for v in row])
    L=arb_mat(n,n+1)
    U=arb_mat(n,n+1)
    for (i,j),x in H.items():
        lo,hi=pair(x)
        for s,t in ((i,j),(j,i)):
            L[s,t],U[s,t]=a(lo),a(hi)
    for i,x in enumerate(force):
        lo,hi=pair(x)
        L[i,n],U[i,n]=a(lo),a(hi)
    # R = R_+ - R_-, so endpoint monotonicity applies before roundout.
    low=pos*L-neg*U
    high=pos*U-neg*L
    def ends(i,j):
        a,b=q(str(low[i,j].lower().fmpq())),q(str(high[i,j].upper().fmpq()))
        require(a<=b,"outward ordered linear image bounds")
        return a,b
    hrows=[sum(max(map(abs,ends(i,j))) for j in range(n)) for i in range(n)]
    return [ends(i,n) for i in range(n)],hrows


def classify(force,H,anchor,R,zq,radius):
    pre,hrows=bound(R,force,H)
    zq,radius=q(zq),q(radius)
    require(zq>=0 and radius>0 and len(anchor)==len(pre),"same whole cached derivative and signed point anchor")
    full=[add(x,pair(y)) for x,y in zip(pre,anchor,strict=True)]
    eta=max(max(abs(a),abs(b)) for a,b in full)
    lower=max(Fraction(0) if a<=0<=b else min(abs(a),abs(b)) for a,b in full)
    b=max(hrows)
    z=zq+b
    lhs=eta+z*radius
    selfmap=z<1 and lhs<radius
    exclusion=lower>(1+z)*radius
    require(not (selfmap and exclusion),"contradictory diagnostics")
    posterior=eta/(1-z) if selfmap else None
    return {"eta_upper_exact":str(eta),"RF_lower_exact":str(lower),"ZQ_upper_exact":str(zq),
            "b_H_upper_exact":str(b),"z_upper_exact":str(z),"selfmap_lhs_upper_exact":str(lhs),
            "strict_selfmap":selfmap,"strict_cube_exclusion":exclusion,
            "selfmap_margin_lower_exact":str(radius-lhs),"exclusion_margin_lower_exact":str(lower-(1+z)*radius),
            "posterior_error_upper_exact":None if posterior is None else str(posterior),
            "signed_R_DF_rows_exact":list(map(pack,pre)),"signed_R_F_I_rows_exact":list(map(pack,full)),
            "RH_absolute_row_bounds_exact":list(map(str,hrows))}


def observer_minimum(arith,target,certified_spd,norm):
    shift=q(target["archived_center_shift_exact"])
    xi=q(target["original_addition_gap_exact"])
    goal=q(target["observer_goal_exact"])
    require(type(certified_spd) is bool and 0<=shift<goal and xi>=0 and norm>0,"same archive and finite norm")
    if not arith["strict_selfmap"]:
        return {"archived_shift_plus_posterior":None,"archived_state_budget_proved":False,
                "physical_nodal_error_um_upper_exact":None,"normalized_H_euclidean_lower_exact":None,
                "strict_minimum_proved":False}
    e=q(arith["posterior_error_upper_exact"])
    require(0<=e<S.RADIUS,"same-cube Banach posterior")
    e0=shift+e
    m=q(arith["z_upper_exact"])
    return {"archived_shift_plus_posterior":str(e0),"archived_state_budget_proved":e0<goal,
            "physical_nodal_error_um_upper_exact":str(Fraction(67,2)*(e0+xi)),
            "strict_minimum_proved":certified_spd,
            "normalized_H_euclidean_lower_exact":str((1-m)/norm) if certified_spd else None}

def observed_producer(observed_h,admitted_h,bindings):
    obs=read(PRODUCER_OS,observed_h)
    require(obs["status"]==OS_STATUS and V.positive_int(obs["pid"]) and V.positive_int(obs["unified_exec_session"])
            and V.zero(obs["observed_exit_code"]) and obs["producer_reduction_input_admission_sha256"]==admitted_h,
            "full externally observed global producer process, separately bound")
    expected={rel(D.ACTUAL/(k+".json")) for k in ("started","terminal","results")}
    expected|={rel(D.ACTUAL/c/(k+".json")) for c in W.FIRST for k in ("assembled_operands","results")}
    files=obs["sha256"]
    require(len(files)==11 and set(files)==expected,"exact four-case global producer file closure")
    def load(path):
        return read(path,files[rel(path)])
    start,end,receipt=(load(D.ACTUAL/(k+".json")) for k in ("started","terminal","results"))
    require(all(V.positive_int(x["pid"]) and x["pid"]==obs["pid"] and x["bindings"]==bindings
                and x["case_order"]==list(W.FIRST) for x in (start,end,receipt))
            and start["status"]=="RUNNING" and end["status"]=="FINISHED" and V.zero(end["exit_code"])
            and end["results_sha256"]==files[rel(D.ACTUAL/"results.json")]
            and receipt["status"]=="COMPLETE_FIRST_WAVE_SIGNED_REDUCTION_AWAITS_INDEPENDENT_GLOBAL_REVIEW"
            and receipt["root_asserted"] is False and receipt["full_five_gate_complete"] is False
            and tuple(x["case_id"] for x in receipt["cases"])==W.FIRST, "actual unpromoted global producer terminal")
    return obs,receipt,load


def assembled(a,case,phases,raw_parts):
    folder,summary,phase,native_receipt=phases[case]
    opts,d,tail,proofs=W.reviewed_context(a,case)
    ctx=d["ctx"]
    t=opts["target"]
    plan=read(S.PLAN/case/"plan.json",opts["plan_sha256"])
    cat={(r["kind"],r["slot"]):r for r in d["entry"]["component_catalogue"]}
    F=[(Fraction(0),Fraction(0)) for _ in range(N)]
    H={}
    totals=Counter()
    maps=[]
    nativepins={}
    for meta in summary["case_components"]:
        job=(meta["kind"],meta["slot"])
        kind,slot=job
        catrow=cat[job]
        proof=proofs[job]
        name=f"{kind}_{slot:03d}.json"
        raw_folder,manifest=raw_parts[case,kind]
        rawpath=raw_folder/name
        leafpath=folder/name
        ids=catrow["global_indices"]
        require(meta["global_indices"]==ids==proof["global_indices"]
                and type(slot) is int and meta["native_file"]==rel(leafpath) and meta["packet_file"]==rel(rawpath)
                and meta["packet_sha256"]==manifest["components"][slot]["sha256"],
                "same signed point/H maps from exact admitted local order")
        own=read(leafpath,meta["native_sha256"])
        p=read(rawpath,meta["packet_sha256"])
        require(own["phase"]==phase and own["status"]==D.LEAF_PASS and own["packet_sha256"]==meta["packet_sha256"]
                and own["selector_leaf_content_sha256"]==p["proof_sha256"]==W.B.digest(proof)
                and own["record_id"]==p["record_id"]==t["record_id"] and own["global_indices"]==p["global_indices"]==ids
                and own["kind"]==p["kind"]==kind and type(own["slot"]) is type(p["slot"]) is int
                and own["slot"]==p["slot"]==slot and own["root_asserted"] is p["root_asserted"] is False
                and own["only_independently_valid_enclosures_used_for_hull_truth"] is True
                and own["overlap_not_used_as_proof"] is True and own["checked_cache_sites"]==(64 if kind=="bulk" else 12)
                and p["selector_plan_sha256"]==opts["plan_sha256"] and p["trial_box_sha256"]==t["trial_box_sha256"]
                and p["worker_identity"]=={k:manifest[k] for k in ("pid","supervisor_pid","stage","case_id","part","bindings")}
                and q(own["radius_exact"])==q(p["whole_box_hessian"]["radius_exact"])==S.RADIUS
                and own["point_force_scope"]==p["point_force"]["scope"]==S.POINT
                and own["whole_H_scope"]==p["whole_box_hessian"]["scope"]==S.BOX,
                "freshly matched signed own mathematical proof and original raw component provenance")
        ij=list(map(list,combinations_with_replacement(range(len(ids)),2)))
        raw_F=p["point_force"]
        raw_H=p["whole_box_hessian"]["H_I_minus_H_cache_upper_exact"]
        own_t=own["independent_point_and_whole_tails_exact"]
        require(len(own_t["F"])==len(ids) and own["producer_force_tails_exact"]==raw_F["force_tail_upper_exact"]
                and all(0<=q(x)<=tail[0] for x in own_t["F"]+raw_F["force_tail_upper_exact"])
                and [v[:2] for v in own_t["H"]]==[v[:2] for v in raw_H]==ij
                and all(0<=q(v[2])<=tail[1] and all(type(x) is int for x in v[:2]) for v in own_t["H"])
                and all(len(v)==4 and all(type(x) is int for x in v[:2]) for v in raw_H),
                "separate point and whole-box independent tail and typed H")
        fr,nf=validated_channel(own,"point_force",list(range(len(ids))),raw_F["F_I_minus_F_cache_intervals_exact"])
        hr,nh=validated_channel(own,"whole_box_H",ij,[v[2:] for v in raw_H])
        _,nc=validated_channel(own,"point_same_Y_cache_correction",list(range(len(ids))),raw_F["Qstar_minus_cache_correction_exact"])
        for k,n in zip(D.CHANNELS,(nf,nh,nc),strict=True):
            totals[k]+=n
        scatter_add(F,H,ids,fr,hr)
        maps.append({"kind":kind,"slot":slot,"global_indices":ids})
        nativepins[rel(leafpath)]=meta["native_sha256"]
    require(len(maps)==168 and dict(totals)==summary["producer_containment_miss_counts"],
            "all signed hulls and original diagnostics independently reconstructed")
    R=[[V.A.T.rational(x) for x in row] for row in ctx.R]
    norm=max(sum(abs(x) for x in row) for row in R)
    anchor=list(map(pair,plan["signed_cached_RF_center_rows_exact"]))
    zq=q(plan["sizing"]["ZQ_exact"])
    require(V.A.a_sha(ctx.R)==t["R_binary_sha256"] and norm==q(plan["sizing"]["R_inf_exact"])
            and V.A.a_sha(ctx.u)==t["center_binary_sha256"] and ctx.drive==q(t["exact_drive"])
            and ctx.row["point_cache_sha256"]==t["point_cache_sha256"] and len(anchor)==N,
            "same frozen exact binary R, point cached anchor and derivative")
    return t,ctx,plan,R,norm,anchor,zq,F,H,maps,nativepins,native_receipt,dict(totals)


def own_spd(spd,t,ctx,plan,norm):
    if spd is None:
        return None
    records,phase,sha_receipt=spd
    row=records[t["case_id"]]
    path=D.SPD/"actual_replay"/(t["case_id"]+".json")
    require(row["receipt"]==rel(path),"canonical native SPD leaf")
    leaf=read(path,row["sha256"])
    if leaf["status"]!=D.SPD_PASS:
        require(row["status"]==leaf["status"],"SPD status mismatch")
        return None
    ident=leaf["target_identity"]
    targets={"pilot_id":t["case_id"],"record_id":t["record_id"],"cache_packet_sha256":t["point_cache_sha256"],
        "R_binary_sha256":t["R_binary_sha256"],"center_binary_sha256":t["center_binary_sha256"],
        "field_sha256":t["original_observer_field_sha256"],"drive_exact":t["exact_drive"],
        "energy_scale_binary_exact":t["energy_scale_binary_exact"],"model_contract_sha256":t["model_contract_sha256"],
        "embedded_cache_sha256":ctx.row["frozen_binary_cache_sha256"],"parameter_cell_sha256":ctx.row["parameter_cell_sha256"],
        "dose_exact":str(ctx.dose),"gain_exact":str(ctx.base["common_gain"])}
    bounds=list(map(q,leaf["independent_row_absolute_sum_upper_exact"]))
    require(row["status"]==D.SPD_PASS and leaf["replay_phase"]==phase and ident==row["target_identity"]==targets
            and row["trial_box_sha256"]==leaf["trial_box_sha256"]==t["trial_box_sha256"]
            and leaf["R_file_sha256"]==t["point_cache_sha256"] and leaf["R_file"]==ctx.row["point_cache_file"]
            and leaf["native_cached_derivative_case_sha256"]==plan["bindings"]["cached_derivative_case_sha256"]
            and leaf["native_geometry_case_sha256"]==plan["bindings"]["geometry_replay_case_sha256"]
            and q(leaf["archived_shift_exact"])==q(t["archived_center_shift_exact"])
            and q(leaf["archived_addition_gap_exact"])==q(t["original_addition_gap_exact"])
            and q(leaf["exact_R_infinity_norm"])==norm
            and len(bounds)==N and min(bounds)>=0 and max(bounds)==q(leaf["independent_q_infinity_upper_exact"])<1
            and leaf["strict_congruence_bound"] is True and leaf["root_or_minimum_claimed"] is False,
            "same-target strict exact Sym(R) premise required")
    return {"exact_Sym_R_native_receipt_sha256":sha_receipt,"exact_Sym_R_native_case_sha256":row["sha256"]}


def compare_producer(a,b,old):
    x=old["arithmetic"]
    require(not (a["strict_selfmap"] and x["strict_exclusion_of_this_cube"])
            and not (a["strict_cube_exclusion"] and x["strict_contraction_and_selfmap"]),
            "global classifications mathematically contradict")
    # Outputs are alternative valid enclosures. Overlap is only diagnostic.
    for n,k in (("signed_R_F_I_rows_exact","signed_RF_I_rows_exact"),
                ("signed_R_DF_rows_exact","signed_RF_defect_rows_exact")):
        xx=list(map(pair,a[n]))
        yy=list(map(pair,x[k]))
        require(len(xx)==len(yy)==N and all(u[0]<=v[1] and v[0]<=u[1] for u,v in zip(xx,yy,strict=True)),
                "global signed interval systems disagree")
    return {"producer_strict_selfmap":x["strict_contraction_and_selfmap"],"own_strict_selfmap":a["strict_selfmap"],
            "producer_strict_exclusion":x["strict_exclusion_of_this_cube"],"own_strict_exclusion":a["strict_cube_exclusion"],
            "producer_z_upper_exact":x["z_upper_exact"],"own_z_upper_exact":a["z_upper_exact"],
            "producer_H_row_norm_upper_exact":x["b_upper_exact"],"own_H_row_norm_upper_exact":a["b_H_upper_exact"],
            "actual_observer_transfer_still_required":b["archived_state_budget_proved"]}
def actual_case(a,case,phases,parts,spd,old_rows,load,mark):
    t,c,plan,R,norm,anchor,zq,F,H,maps,pins,nh,misses=assembled(a,case,phases,parts)
    flat_H=[[i,j,*pack(v)] for (i,j),v in sorted(H.items())]
    pdir=D.ACTUAL/case
    operands=load(pdir/"assembled_operands.json")
    old=load(pdir/"results.json")
    old_row=V.A.pick(old_rows,case)
    require(sha(pdir/"assembled_operands.json")==old_row["assembled_operands_sha256"]==old["assembled_operands_sha256"]
            and sha(pdir/"results.json")==old_row["results_sha256"]
            and operands["target"]==old["target"]==t and old["status"]=="FIRST_WAVE_SPATIAL_REDUCTION_AWAITS_INDEPENDENT_GLOBAL_REVIEW"
            and operands["root_asserted"] is old["root_asserted"] is False and old["exact_spatial_root_certified"] is False
            and all(V.zero(old[k]) for k in ("new_spatial_integrals","new_AD_H_moments_or_solves"))
            and operands["cached_anchor_already_preconditioned"] is True and operands["maps"]==maps
            and operands["signed_native_R_F_Q_exact"]==list(map(pack,anchor))
            and operands["F_I_minus_F_Q_hull_global_exact"]==list(map(pack,F))
            and operands["H_I_minus_H_Q_hull_global_upper_exact"]==flat_H
            and operands["native_component_sha256"]==pins and operands["whole_case_hull_receipt_sha256"]==nh
            and operands["producer_containment_miss_counts"]==misses
            and operands["selector_and_cached_geometry_bindings"]==plan["bindings"]
            and q(operands["ZQ_whole_cube_upper_exact"])==zq and q(operands["exact_R_infinity_norm"])==norm,
            "independent exact global scatter must equal the same produced signed array, not merely its norm")
    arithmetic=classify(F,H,anchor,R,zq,S.RADIUS)
    premise=own_spd(spd,t,c,plan,norm)
    transfer=observer_minimum(arithmetic,t,premise is not None,norm)
    diagnostics=compare_producer(arithmetic,transfer,old)
    root=arithmetic["strict_selfmap"]
    only_cube_exclusion=arithmetic["strict_cube_exclusion"]
    doc={"status":("PASS_INDEPENDENT_FIRST_WAVE_EXACT_SPATIAL_ROOT" if root else
            "SUFFICIENT_INDEPENDENT_EXCLUSION_OF_ONLY_THIS_CUBE" if only_cube_exclusion else
            "UNRESOLVED_INDEPENDENT_FIRST_WAVE_EXACT_SPATIAL_TEST"),
         **mark,"case_id":case,"target":t,"radius_exact":str(S.RADIUS),"normalized_energy_H":True,
         "scope":SCOPE,"precision_bits":BITS,"independent_trusted_hull_native_receipt":nh,
         "independent_actual_hull_leaves_sha256":pins,"exact_endpoint_scatter_equal_all":True,
         "own_signed_hulls_global_F_exact":list(map(pack,F)),"own_signed_hulls_global_H_upper_exact":flat_H,
         "same_signed_cached_RF_rows_exact":list(map(pack,anchor)),"same_R_norm_exact":str(norm),
         "arithmetic":arithmetic,"producer_comparison":diagnostics,"observer_and_optional_minimum":transfer,
         "same_R_positive_symmetric_part_premise":premise,"closed_cube_invertible_R_from_contraction":root,
         "new_spatial_integrations_H_moments_AD_or_solves":0,"root_branch_or_continuum_or_identified_set_claimed":False,
         "original_observer_field_replaced":False,"full_five_gate_complete":False}
    h=save(OUT/(case+".json"),doc)
    return {"case_id":case,"native_file":rel(OUT/(case+".json")),"native_sha256":h,
            "status":doc["status"],"strict_local_minimum":transfer["strict_minimum_proved"],
            "observer_state_budget_proved":transfer["archived_state_budget_proved"]}


def run(args):
    p=read(HERE/"protocol.json",args.protocol_sha256)
    f=read(HERE/"preflight.json",args.preflight_sha256)
    require(args.review.resolve()==SOURCE_REVIEW.resolve(),"canonical independent global source approval")
    r=read(args.review,args.review_sha256)
    expected={"source_sha256":sha(SOURCE),"protocol_sha256":args.protocol_sha256,
              "preflight_sha256":args.preflight_sha256,
              "producer_root_terminal_admission_sha256":args.producer_observed_sha256,
              "producer_reduction_input_admission_sha256":args.producer_admission_sha256}
    require(r["status"]==REVIEW and V.yes(r["checks"]) and all(r[k]==v for k,v in expected.items())
            and r["new_spatial_integrations_authorized"] is False and p==static() and f["status"]=="PASS"
            and V.yes(f["checks"]) and all(V.zero(f[k]) for k in ZERO) and f["source_sha256"]==sha(SOURCE)
            and f["protocol_sha256"]==args.protocol_sha256 and sha(HERE/"replay_frozen.py")==sha(SOURCE),
            "current cold source, actual observed global producer and no other phase approved")
    a,w,j,b=D.gate(SimpleNamespace(protocol_sha256=PINS[D.PREP/"protocol.json"],
         preflight_sha256=PINS[D.PREP/"preflight.json"],review=D.REVIEW_FILE,review_sha256=PINS[D.REVIEW_FILE],
         admission=D.ADMIT_FILE,admission_sha256=args.producer_admission_sha256))
    require(j["root_producer_admission_sha256"]==RAW_SHA,"same 715-artifact externally observed first wave")
    obs,old,load=observed_producer(args.producer_observed_sha256,args.producer_admission_sha256,b)
    raw,parts=D.observed_raw(RAW_SHA,a,w)
    phases=D.open_all_hull_metadata(a,j,raw)
    spd_spec=j["optional_spd"]
    require(spd_spec is None or spd_spec["observed_os_exit_sha256"]==PINS[D.SPD/"actual_replay/os_exit_observation.json"],
            "optional exact same-R proof must match root-observed 32 native proof")
    spd=D.optional_spd(spd_spec)
    OUT.mkdir(exist_ok=False)
    mark={**expected,"pid":os.getpid(),"root_source_review_sha256":args.review_sha256,"producer_observed_pid":obs["pid"]}
    save(OUT/"started.json",{"status":"RUNNING",**mark})
    terminal={"status":"FAILED","exit_code":1,**mark}
    rows=[]
    try:
        with patch.object(D,"case_hulls",side_effect=AssertionError("producer endpoint scatter forbidden")), \
             patch.object(D,"reduce_case",side_effect=AssertionError("producer mathematical execution forbidden")), \
             patch.object(D.K,"assemble",side_effect=AssertionError("producer scatter kernel forbidden")), \
             patch.object(D.K,"contract",side_effect=AssertionError("producer contraction kernel forbidden")):
            for case in W.FIRST:
                row=actual_case(a,case,phases,parts,spd,old["cases"],load,mark)
                rows.append(row)
                print(row,flush=True)
        require(tuple(x["case_id"] for x in rows)==W.FIRST,"complete independent four-case end")
        result={"status":"COMPLETE_INDEPENDENT_FIRST_WAVE_SIGNED_ENDPOINT_REPLAY",
                "scope":SCOPE,**mark,"case_results":rows,
                "exact_spatial_roots":sum(x["status"]=="PASS_INDEPENDENT_FIRST_WAVE_EXACT_SPATIAL_ROOT" for x in rows),
                "strict_local_minima":sum(x["strict_local_minimum"] for x in rows),
                "new_spatial_integrations_H_moments_AD_or_solves":0,"archived_observers_replaced":False,
                "full_five_gate_complete":False}
        h=save(OUT/"receipt.json",result)
        terminal.update(status="FINISHED",exit_code=0,receipt_sha256=h)
        print({"receipt_sha256":h},flush=True)
    except BaseException as error:
        save(OUT/"failure.json",{"error":str(error),"traceback":traceback.format_exc(),"completed_cases":rows,**mark})
        raise
    finally:
        save(OUT/"terminal.json",terminal)


def fixtures():
    own={"signed_channels":{"test":{"rows":[{"index":0,"independent_exact":["-2","1"],"producer_exact":["0","3"],
           "hull_exact":["-2","3"],"producer_contains_independent":False,"disjoint":False}],
           "producer_containment_misses":[0],"disjoint_indices":[]}}}
    good,n=validated_channel(own,"test",[0],[["0","3"]])
    require(good==[(Fraction(-2),Fraction(3))] and n==1,"nonnested union valid")
    wrong=dict(own["signed_channels"]["test"]["rows"][0],hull_exact=["0","1"])
    bad={"signed_channels":{"test":{**own["signed_channels"]["test"],"rows":[wrong]}}}
    require(V.A.reject(lambda:validated_channel(bad,"test",[0],[["0","3"]])),"reject mere intersection")
    F=[(Fraction(0),Fraction(0)) for _ in range(2)]
    H={}
    scatter_add(F,H,[1,0],[pair(("3","4")),pair(("-2","1"))],[pair(("0","1")),pair(("-1","0")),pair(("2","3"))])
    require(F==[(-2,1),(3,4)] and H=={(1,1):(0,1),(0,1):(-1,0),(0,0):(2,3)},"permuted signed scatter")
    r=[[Fraction(2),Fraction(-3)],[Fraction(-5),Fraction(7)]]
    ff,hh=bound(r,F,H)
    def exact_oracle(invals):
        return [(sum(min(q(k)*a,q(k)*b) for k,(a,b) in zip(row,invals,strict=True)),
                 sum(max(q(k)*a,q(k)*b) for k,(a,b) in zip(row,invals,strict=True))) for row in r]
    left=exact_oracle((H[0,0],H[0,1]))
    right=exact_oracle((H[0,1],H[1,1]))
    require(ff==exact_oracle(F) and hh==[max(map(abs,a))+max(map(abs,b)) for a,b in zip(left,right,strict=True)],
            "independent positive/negative R endpoint arithmetic matches exact signed oracle")
    a=classify([[Fraction(-1,8)]*2],{(0,0):[Fraction(1,8)]*2},[[Fraction(1,4)]*2],[[2]],Fraction(1,4),S.RADIUS)
    b=observer_minimum(a,{"archived_center_shift_exact":"1/16","original_addition_gap_exact":"1/32",
                          "observer_goal_exact":"1/8"},True,Fraction(2))
    require(a["strict_selfmap"] and q(a["eta_upper_exact"])==0 and q(a["z_upper_exact"])==Fraction(1,2)
            and q(b["physical_nodal_error_um_upper_exact"])==Fraction(201,64)
            and q(b["normalized_H_euclidean_lower_exact"])==Fraction(1,4),"native cached anchor once and original observer origin")
    boundary=classify([(S.RADIUS,S.RADIUS)],{},[(0,0)],[[1]],0,S.RADIUS)
    require(boundary["strict_selfmap"] is False and boundary["strict_cube_exclusion"] is False,
            "no promotion on strict boundary")
    return {"same_local_union_recomputed_and_intersection_rejected":True,
            "signed_positive_negative_split_matches_exact_fraction_oracle":True,
            "already_preconditioned_anchor_once_and_archive_shift_xi":True,
            "strict_boundary_retained_unresolved":True,
            "optional_same_R_positive_symmetric_part_supported_without_new_SPD_job":True}


def prepare():
    p=static()
    checks=fixtures()
    require(V.yes(checks) and not OUT.exists(),"cold independent no actual global closure opened")
    ph=save(HERE/"protocol.json",p)
    fh=save(HERE/"preflight.json",{"status":"PASS",**ZERO,"source_sha256":sha(SOURCE),"protocol_sha256":ph,"checks":checks})
    with (HERE/"replay_frozen.py").open("xb") as f:
        f.write(SOURCE.read_bytes())
    print({"source_sha256":sha(SOURCE),"protocol_sha256":ph,"preflight_sha256":fh},flush=True)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    sub=cli.add_subparsers(dest="phase",required=True)
    sub.add_parser("prepare")
    actual=sub.add_parser("run")
    for k in ("protocol","preflight","review","producer-observed","producer-admission"):
        actual.add_argument("--"+k+"-sha256",required=True)
        if k=="review":
            actual.add_argument("--review",required=True,type=lambda x:Path(x).resolve())
    args=cli.parse_args()
    prepare() if args.phase=="prepare" else run(args)
