"""Independent exact endpoint readback of the cached-quadrature observer transfer."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
APP = HERE.parent
DESIGN = APP.parents[1]
RUN = DESIGN.parent
ROOT = next(p for p in RUN.parents if (p/"PLAN.md").is_file())
SOURCE = Path(__file__).resolve()
LEN = Q(67,2)
SIDES = ("max_contrast_change","max_negative_contrast_change")
OBS = ("maximum_shortening_um","region_averaged_shortening_um","apical_max_shortening_um","basal_max_shortening_um")
VARIANT = "current_disk__spatial_max"
RECASES = {"max_contrast_change_small_13","max_negative_contrast_change_small_16"}
PINS = {
    APP/"results.json":"1df2bc372edd12687fad17edf63013509994e6bbe41c02150cefb25952eab9f6",
    APP/"root_source_review.json":"5914ec9e27df980fe403fe5cfa727ae9aa54e8c0aa79637745ac27422df93ac2",
    APP/"apply.py":"2a7adc3e1dc8e05c6df568efa048b04c5574baaafa451a32cdb97b3af4a5b990",
    APP/"protocol.json":"3f8dc8d5018cfcdae9fe0362be196730686d7ff6f208a8fcdebd2c6fde1f4613",
    APP/"preflight.json":"1a31927c679951f1c996fb09c5a86ce76bad434a10cf2a23f5d6e541a96c0261",
    DESIGN/"budgets_v1/actual/results.json":"689948766db950ee4ddf399e6731aea08890cc7a1f84c5615d19e7ded68b0872",
    DESIGN/"budgets_v1/independent_actual/receipt.json":"978a2a4095bb8735103a5ee2c11ad266a6b0e052e22af2b7a36ef29f1f38265c",
    RUN/"field_manifest.json":"308ca4c05e476ba6a14cbc7b8d34018983b8f11263e2c224488672b35eb38ac0",
    RUN/"preparation/protocol.json":"9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94",
    RUN/"independent_actual_preparation/actual_replay/receipt.json":"c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd",
}
RECDIR = APP.parent/"remaining_cached_30_v1/recenter_v1"
PINS[RECDIR/"independent_actual_preparation/actual_replay/receipt.json"] = "ac6013e52166c434a0133219ce45db36daf8d6f32e663f8c214d3d48756cd025"
PINS[RECDIR/"independent_actual_preparation/actual_replay/os_exit_observation.json"] = "251d443883f9225f160a43d3894c49d506bf99a0597224134acd4dc607869322"
LINKS = {}


def require(ok, why):
    if not ok:
        raise ValueError(why)


def read(path, expected, parse=True):
    path = (ROOT/path).resolve()
    name = path.relative_to(ROOT).as_posix()
    raw = path.read_bytes()
    hash_ = sha256(raw).hexdigest()
    require(hash_ == expected == LINKS.get(name,expected),f"digest mismatch {name}")
    LINKS[name]=hash_
    return json.loads(raw) if parse else None


def strict_pair(x):
    a,b=map(Q,x)
    require(a<=b,"inverted interval")
    return a,b


def unique(rows, **kwargs):
    leaf, = [r for r in rows if all(r[k]==v for k,v in kwargs.items())]
    return leaf


def edges_at(reference, elements, z):
    rim_edges = {}
    for a,b,c,m1,m2,m3 in elements:
        for i,j,m in ((a,b,m1),(b,c,m2),(c,a,m3)):
            if reference[i][1] == reference[j][1] == z:
                start,end=sorted((i,j),key=lambda k:reference[k][0])
                require(start!=end and reference[start][0]<reference[end][0] and
                        rim_edges.get((start,end),m)==m,"nonconforming reference trace edge")
                rim_edges[start,end]=m
    ordered=[[i,m,j] for (i,j),m in rim_edges.items()]
    ordered.sort(key=lambda e:reference[e[0]][0])
    require(len(ordered)==12 and all(x[-1]==y[0] for x,y in zip(ordered,ordered[1:])),
            "need twelve conforming, ordered P2 edges")
    return ordered


def original_curves(state,b):
    ref=[tuple(map(Q,n)) for n in state["reference_nodes"]]
    cur=[tuple(map(Q,n)) for n in state["current_nodes"]]
    u=list(map(Q,state["free_displacement"]))
    p=state["parameters"]
    rim,height=Q(p["outer_radius"]),Q(p["height"])
    free=[2*i+j for i,(x,_) in enumerate(ref) for j in (0,1)
          if x!=rim and (j==1 or x!=0)]
    require(len(ref)==len(cur)==325 and len(free)==len(u)==611 and len(state["triangles"])==144
            and free==b["center_and_model"]["free_map"],"verified topology/free dimension differs")
    correction=dict(zip(free,u,strict=True))
    exact=[[x+correction.get(2*i+j,Q(0)) for j,x in enumerate(r)] for i,r in enumerate(ref)]
    wrapped=[tuple(Q(float(x)) for x in node) for node in exact]
    require(wrapped==cur,"serialized coordinates fail to match the archived sum")
    xi=max(abs(x-y) for e,c in zip(exact,cur,strict=True) for x,y in zip(e,c,strict=True))
    sides={}
    for name,z in (("top",height),("bottom",Q(0))):
        edges=edges_at(ref,state["triangles"],z)
        i,j=edges[0][0],edges[-1][-1]
        anchor=[2*i,2*j]
        require(not (set(anchor)&set(free)) and [exact[i][0],exact[j][0]]==[Q(0),rim],
                "both current endpoint radii must remain fixed")
        rec=b["center_and_model"]["fixed_trace_endpoints"][name]
        require(rec["radial_dofs"]==anchor and rec["trace_nodes"]==edges,
                "archived fixed trace binding mismatch")
        sides[name]=[[tuple(x*LEN for x in cur[k]) for k in e] for e in edges]
    require(Q(b["center_and_model"]["reference_height_um_exact"])==LEN*height,"changed original reference")
    return sides, xi, u


def curve(sides,epsilon,roi):
    require(sides[0][0][0]<=roi[0]<roi[1]<=sides[-1][-1][0],"common fixed current-radius domain missing")
    # Derivatives of each quadratic Lagrange coordinate at t=0 and t=1.
    ds=[[(tuple(-3*p+4*m-q for p,m,q in zip(*edge,strict=True))),
          (tuple(p-4*m+3*q for p,m,q in zip(*edge,strict=True)))] for edge in sides]
    dr=min(d[0] for edge in ds for d in edge)
    require(dr>8*epsilon,"radial monotonicity exhausts verified error")
    kappa=max(max(abs(d[1]) for d in edge)/min(d[0] for d in edge) for edge in ds)
    return dr, kappa, Q(5,4)*(1+kappa)*epsilon


def sigma_check(leaf, root, archive, u, rec):
    k=archive["case_id"]
    if k not in RECASES:
        require(Q(root["sigma"])==0,"shift assigned to non-recentered state")
        return Q(0)
    meta=unique(rec["case_results"],case_id=k)
    require(leaf["status"]==meta["status"]=="PASS_INDEPENDENT_RECENTERED_CACHED_LOCAL_ROOT_AND_OLD_OBSERVER_BUDGET"
            and meta["native_sha256"]==root["cached_proof_leaf_sha256"] and leaf["original_field_or_scores_replaced"] is False,
            "missing native independently proved recenter")
    f=RECDIR/"actual"/k/"exact_binary_step.json"
    key=f.relative_to(ROOT).as_posix()
    require(rec["producer_metadata_and_packet_sha256"][key]==leaf["exact_step_sha256"],"recorded binary shift is mismatched")
    step=read(f,leaf["exact_step_sha256"])
    cols=step["ordered_rows"]["columns"]
    data=step["ordered_rows"]["values"]
    require(cols==["old_u","native_signed_RF_midpoint","exact_unrounded_u1","rounded_binary_u1","exact_binary_shift"]
            and len(data)==len(u)==611,"wrong shift input")
    actual=[]
    for start,row in zip(u,data,strict=True):
        old,mid,ideal,new,shift=map(Q,row)
        require(old==start and old-mid==ideal and Q(float(ideal))==new and new-old==shift,
                "binary recenter round-once row inconsistent")
        actual.append(abs(new-old))
    sigma=max(actual)
    require(sigma==Q(step["sigma_exact"])==Q(root["sigma"])==Q(leaf["actual_binary_center_shift_sigma_exact"])>0,
            "reported sigma does not equal observed maximum displacement")
    return sigma


def radial_error(source, budget, root, bfield, state, rec, leaf):
    if "own_radius_attempts" in leaf:
        active=unique(leaf["own_radius_attempts"],radius_exact=leaf["chosen_radius_exact"])
        eta,z,r,stored_post=map(Q,(active["eta_upper_exact"],active["ZQ_upper_exact"],
                                    active["radius_exact"],active["posterior_error_upper_exact"]))
        lhs=Q(active["strict_selfmap_lhs_upper_exact"])
        goal=Q(leaf["observer_radius_budget_exact"])
    else:
        eta=Q(leaf.get("eta_upper_exact",leaf.get("independent_eta_upper_exact")))
        z=Q(leaf.get("ZQ_upper_exact",leaf.get("independent_ZQ_upper_exact",leaf.get("independent_z_upper_exact"))))
        r=Q(leaf["radius_exact"])
        stored_post=Q(leaf.get("posterior_error_upper_exact",leaf.get("independent_posterior_error_upper_exact")))
        lhs=Q(leaf.get("strict_selfmap_lhs_upper_exact",leaf.get("independent_selfmap_lhs_upper_exact")))
        goal=Q(leaf.get("archived_observer_goal_exact",leaf.get("science_radius_budget_exact",leaf.get("observer_goal_exact"))))
    require(0<=z<1 and 0<=eta and lhs==eta+z*r<r==Q(root["radius"]) and r<goal==Q(budget["e_goal_dimensionless_exact"]),
            "strict contraction differs from independently verified leaf")
    e=eta/(1-z)
    require(e==stored_post==Q(root["posterior"]),"posterior arithmetic drift")
    key="archived_field_sha256" if source["case_id"] in RECASES else "field_sha256"
    mech=budget["center_and_model"]
    drive=Q(leaf.get("exact_requested_strength",leaf.get("exact_strength")))
    require(leaf[key]==source["field_sha256"]==bfield["field_sha256"] and leaf["record_id"]==source["record_id"]
            and drive==Q(mech["dose_exact"])*Q(mech["gain_exact"])==Q(bfield["step"]["requested_strength_exact"]),
            "root and frozen archived drive do not refer to same target")
    curves,xi,u=original_curves(state,budget)
    sigma=sigma_check(leaf,root,source,u,rec)
    epsilon=LEN*(sigma+e+xi)
    require(sigma+e<goal and xi==Q(root["archived_xi"])==Q(mech["node_addition_gap_dimensionless_exact"])
            and epsilon==Q(root["epsilon_um"]),"shift, posterior and archived storage must each count")
    return curves,epsilon,(sigma,e,xi)


def execute():
    data={p:read(p,h,p.suffix==".json") for p,h in PINS.items()}
    actual=data[APP/"results.json"]
    require(actual["status"]=="CACHED_QUADRATURE_PRIMARY_OBSERVER_APPLICATION_AWAITING_INDEPENDENT_READBACK"
            and actual["primary_observer"]==VARIANT and actual["cached_root_fields"]==36
            and actual["primary_calibration_comparisons"]==136 and actual["large_Q_values"]==2
            and actual["exact_spatial_root_observer_results"]==0 and actual["experimental_compatibility_or_PDE_claimed"] is False,
            "result scope or primary field counts changed")
    authority=actual["source_and_actual_sha256"]
    for name, h in authority.items():
        read(name,h,parse=False)
    cold=data[APP/"root_source_review.json"]
    require(cold["status"]=="PASS" and cold["checks"] and all(v is True for v in cold["checks"].values())
            and all(cold[k]==PINS[APP/file] for k,file in (("source_sha256","apply.py"),("protocol_sha256","protocol.json"),
                                                          ("preflight_sha256","preflight.json"))),"root source gate mismatch")
    recdir=RECDIR/"independent_actual_preparation/actual_replay"
    rec=data[recdir/"receipt.json"]
    ro=data[recdir/"os_exit_observation.json"]
    e=read(recdir/"terminal.json",ro["terminal_sha256"])
    require(ro["status"]=="PASS_ROOT_OBSERVED_INDEPENDENT_RECENTER_REPLAY_TERMINAL"
            and type(ro["observed_exit_code"]) is type(e["exit_code"]) is int
            and ro["observed_exit_code"]==e["exit_code"]==0 and ro["pid"]==e["pid"]==rec["pid"]
            and PINS[recdir/"receipt.json"]==ro["receipt_sha256"]==e["receipt_sha256"]
            and rec["cached_local_roots_certified"]==2,"native recenter scope or OS proof differs")
    plan=data[RUN/"preparation/protocol.json"]
    archived=data[RUN/"field_manifest.json"]["fields"]
    primary=data[RUN/"independent_actual_preparation/actual_replay/receipt.json"]
    bd=data[DESIGN/"budgets_v1/actual/results.json"]
    delta=Q(plan["selection"]["conditional_source_allowance_um_exact"])
    require(delta==Q(actual["frozen_allowance_um_exact"])==Q(bd["conditional_source_allowance_um_exact"])>0,
            "original allowance not preserved")
    by_sha={h:name for name,h in authority.items()}
    seen=set(); small=[]; large=[]; transfer=[]
    for a in actual["cases"]:
        endpoint,step=a["endpoint_id"],a["checkpoint"]
        demand_id=f"{endpoint}_{step}"
        require(a["case_id"]==demand_id and demand_id not in seen,"duplicate case or wrong selection")
        seen.add(demand_id)
        b=unique(bd["field_budgets"],endpoint_id=endpoint,checkpoint=step)
        field=unique(archived,endpoint_id=endpoint,checkpoint=step)
        rt=a["root"]
        state=read(field["field"],field["field_sha256"])
        verifier=read(field["verification"],field["verification_sha256"])
        require(field["field_sha256"]==a["field_sha256"]==b["field_sha256"] and
                unique(primary["endpoints"],endpoint_id=endpoint)["endpoint_id"]==endpoint,
                "input changed from original field authority")
        leaf=read(by_sha[rt["cached_proof_leaf_sha256"]],rt["cached_proof_leaf_sha256"])
        faces,eps,terms=radial_error(a,b,rt,field,state,rec,leaf)
        roi=strict_pair(b["fixed_roi_um_exact"])
        require(roi==strict_pair(a["ROI_um_exact"]),"changed ROI")
        E=[]
        for face,side in (("top","apical"),("bottom","basal")):
            dr,slope,error=curve(faces[face],eps,roi)
            got=a["apical_and_basal_graph_error_bounds"][side]
            require(dr-8*eps==Q(got["radial_derivative_margin_exact"])>0 and
                    slope==Q(got["base_graph_Lipschitz_upper_exact"]) and error==Q(got["Eulerian_vertical_error_um_exact"]),
                    "independent P2 derivative or graph error differs")
            E.append(error)
        transfer.append({"case_id":demand_id,"posterior_exact":str(terms[1]),"sigma_exact":str(terms[0]),
                         "archived_xi_exact":str(terms[2]),"epsilon_um_exact":str(eps)})
        if step=="large_03":
            q0=strict_pair(b["stored_Q_um_exact"])
            end=unique(primary["endpoints"],endpoint_id=endpoint)
            require(roi[0]==0 and roi[1]==Q(plan["large_observable"]["shoulder_radius_um_exact"])
                    and q0==strict_pair(end["large_Q_um_exact"])==strict_pair(verifier["surface"]["Q_um_exact"])
                    ==strict_pair(a["stored_Q_um_exact"]),"different original Q or shoulder")
            extra=eps+E[0]
            q1=(q0[0]-extra,q0[1]+extra)
            require(Q(a["Q_absolute_error_um_exact"])==extra and strict_pair(a["cached_root_Q_um_exact"])==q1,
                    "axis plus shoulder error mismatch")
            large.append({"endpoint_id":endpoint,"Q_cached_um_exact":list(map(str,q1)),
                          "error_um_exact":str(extra),"Q_stored_um_exact":list(map(str,q0))})
        else:
            require(b["role"]=="small_record" and step in {f"small_{i:02d}" for i in range(1,18)}
                    and roi==(Q(0),LEN*Q(state["parameters"]["patch_radius"]))
                    and verifier["variant"]==VARIANT and len(a["comparisons"])==4,"wrong small observer convention")
            for obs,extra in zip(OBS,(sum(E),sum(E),E[0],E[1]),strict=True):
                pred=strict_pair(verifier["observer_intervals_um_exact"][obs])
                y=Q(field["step"]["observed_exact"][obs])
                row=unique(a["comparisons"],observable=obs)
                lo,hi=pred[0]-extra,pred[1]+extra
                res=max(Q(0),lo-y,y-hi),max(abs(lo-y),abs(hi-y))
                status="IN" if res[1]<=delta else "OUT" if res[0]>delta else "UNRESOLVED"
                require(pred==strict_pair(row["stored_interval_um_exact"]) and y==Q(row["observed_um_exact"])
                        and extra==Q(row["extra_observer_error_um_exact"]) and (lo,hi)==strict_pair(row["cached_root_interval_um_exact"])
                        and res==strict_pair(row["absolute_residual_interval_um_exact"]) and status==row["comparison_at_original_delta"],
                        "independent signed residual comparison differs")
                small.append({"endpoint_id":endpoint,"case_id":demand_id,"observable":obs,"standing":status,
                              "absolute_residual_um_exact":list(map(str,res)),"strict_slack_um_exact":str(delta-res[1])})
    require(len(seen)==36 and len(small)==136 and len(large)==2
            and all(sum(r["endpoint_id"]==end for r in small)==68 for end in SIDES),"coverage incomplete")
    pos,neg=(strict_pair(unique(large,endpoint_id=end)["Q_cached_um_exact"]) for end in SIDES)
    gap=(pos[0]-neg[1],pos[1]-neg[0])
    counts={end:{s:sum(r["endpoint_id"]==end and r["standing"]==s for r in small)
                 for s in ("IN","OUT","UNRESOLVED")} for end in SIDES}
    require(gap==strict_pair(actual["cached_root_positive_minus_negative_large_Q_um_exact"])
            and (gap[0]>0) is actual["strict_positive_cached_gap"] and counts==actual["per_endpoint_counts"],
            "cached Q gap or counts differ from independent readback")
    minslack=min(small,key=lambda r:Q(r["strict_slack_um_exact"]))
    out={"status":"PASS_INDEPENDENT_ACTUAL_36_CACHED_ROOT_OBSERVERS","scope":VARIANT,
         "cached_fixed_parameter_roots":36,"primary_signed_residuals":136,"large_Q_values":2,"counts":counts,
         "frozen_model_discrepancy_um_exact":str(delta),"minimum_slack":minslack,
         "cache_root_signed_large_gap_um_exact":list(map(str,gap)),"large_Q_cached":large,
         "observer_transfer_terms":transfer,"all_rows":small,"readback_sha256":sha256(SOURCE.read_bytes()).hexdigest(),
         "premises_reused":"Independently reviewed cached roots and archived observer intervals; no observer solver replay. Separate P2 trace and exact signed interval implementation.",
         "current_input_sha256":LINKS,"source_compatibility_or_PDE":False,
         "new_novelty_claim":False,"alternate_observer_robustness_claimed":False,
         "exact_spatial_observer_applications_from_this_readback":0,"full_five_gate_goal_achieved":False}
    f=HERE/"receipt.json"
    with f.open("x",encoding="utf-8",newline="\n") as handle:
        json.dump(out,handle,indent=2)
        handle.write("\n")
    print(json.dumps({"receipt_sha256":sha256(f.read_bytes()).hexdigest(),"bindings":len(LINKS),
                      "counts":counts,"min_slack_um":float(Q(minslack["strict_slack_um_exact"])),
                      "gap_display_um":list(map(float,gap)),"large_Q_display_um":
                      {r["endpoint_id"]:list(map(lambda v:float(Q(v)),r["Q_cached_um_exact"])) for r in large}}))


if __name__=="__main__":
    execute()

