"""Explicit archived-field authority and canonical proof centers for the remaining 32."""
import argparse
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

import transfer_ambiguity_exact_targets as X

SOURCE = Path(__file__).resolve()
ROOT = X.ROOT
OUT = X.OUT / "remaining32_targets_v1"
PREP = OUT / "preparation"
DESIGN = X.OUT / "remaining32_bridge_design_v1"
INV = X.OUT.parent / "full_outcome_coverage_v1"
WAVES = INV / "remaining_cached_30_v1"
M = INV / "remaining_cached_wave1"
R = WAVES / "recenter_v1"
Q, require, sha, doc = X.exact, X.require, X.sha, X.doc
T, E = X.T, X.E
R36, R34, R32 = (Fraction(1, 2**k) for k in (36, 34, 32))
RECENTERED = {"max_negative_contrast_change_small_16", "max_contrast_change_small_13"}
IDS = tuple(sorted(f"{e}_small_{i:02}" for e in X.sweep.ENDPOINTS for i in range(1,17)))
MP = M/"preparation/protocol.json"
WP = WAVES/"preparation/protocol.json"
RP = R/"preparation/protocol.json"
PINS = {
    DESIGN/"manifest.json": "5f043b6a2fac20bdbb02dd3a5c9e58a860b2a4453b246373e28178f99ab49d53",
    DESIGN/"receipt.json": "a3a44dda65b1467088c63808bbb91906bc98009c968a4b7e13e5865dbfe9902b",
    MP: "c1543cbbd2520664313b56114b66a885faaf4d75869f56dbb1aa83a43a78fe3a",
    WP: "e31870b0c6ac61e95a7b2657eab73d8231657af86c65f7a1d1172ff0b45a1e4c",
    RP: "bb5f0aef19bf5a7d5a6f954d7b41f03219bd19056be4e4eb4672454a588bec63",
    ROOT/"scripts/transfer_ambiguity_exact_targets.py": "53ff92a87d164871e04323759b27e015fc8e19d039dc0dd15c2d6080acda7cc8",
}
STATUS = "COLD_REMAINING32_EXACT_SPATIAL_SINGLETON_TARGETS"
NO_MATH = {"new_geometry_certificates":0,"new_AD_derivatives":0,"new_integrations":0,"new_solves":0}


def rel(p):
    return p.resolve().relative_to(ROOT).as_posix()


def canonical(name):
    p = (ROOT / name).resolve()
    require(p.is_relative_to(ROOT) and rel(p)==name, "canonical study-relative path required")
    return p


def read(p, h, opened=None, purpose=None):
    value = X.pinned(p, h)
    if opened is not None:
        opened[rel(p)] = {"sha256":h,"purpose":purpose}
    return value


def a_sha(a):
    return X.array_sha(a)


def pick(items, case):
    hits = [item for item in items if item["case_id"]==case]
    require(len(hits)==1, "one authoritative case required")
    return hits[0]


def scope_for(r, z):
    require(r in (R36,R34,R32) and 0<=z<1, "admitted finite native cached radius required")
    if r == R34:
        return r,"EXISTING_MATCHED_WHOLE_CUBE_CACHED_ZQ"
    if r == R36:
        return r,"EXISTING_MATCHED_SMALL_CUBE_LARGER_TRIAL_REQUIRES_NEW_AD"
    require(z>Fraction(4,5), "new headroom policy requires review for this radius group")
    return R34,"FRESH_VARIATION_REQUIRED_BEFORE_SPATIAL_LAUNCH"


def load_authority():
    X.bootstrap()
    for p,h in PINS.items():
        require(sha(p)==h,"exact frozen authority bytes required")
    d=read(DESIGN/"manifest.json",PINS[DESIGN/"manifest.json"])
    receipt=read(DESIGN/"receipt.json",PINS[DESIGN/"receipt.json"])
    require(receipt["status"]=="PASS_READ_ONLY_REMAINING32_BRIDGE_DESIGN"
            and receipt["actual_consumer_launch_authorized"] is False
            and receipt["bindings"]["manifest.json"]==PINS[DESIGN/"manifest.json"]
            and d["root_or_worker_launch_authorized"] is False and d["full_five_gate_complete"] is False
            and d["actual_fields_caches_or_step_vectors_opened"]==0
            and all(x is True for x in d["checks"].values()),"complete cold census authority required")
    links={**d["input_metadata_and_sources_sha256"], **{rel(p):h for p,h in PINS.items()}}
    for name,h in links.items():
        require(sha(canonical(name))==h,"reviewed census dependency changed")
    rows={row["case_id"]:row for row in d["cases"]}
    require(tuple(rows)==IDS and len(d["cases"])==len(IDS)==32,"complete ordered case set required")
    m=read(MP,PINS[MP])
    w=read(WP,PINS[WP])
    r=read(RP,PINS[RP])
    rec={entry["case"]["case_id"]:entry for entry in r["cases"]}
    require(set(rec)==RECENTERED and set(w["cases"])|set(m["ordered_case_ids"])==set(IDS)
            and not set(w["cases"]).intersection(m["ordered_case_ids"]),"32 unique original case protocols")
    return rows,w["cases"],m["cases"],rec,links


def case_authority(case,authority):
    rows,w,m,rec,links=authority
    require(case in rows,"unknown remaining32 case")
    row=rows[case]
    c=w[case] if case in w else pick(m,case)
    require(all(row[k]==c[k] for k in ("case_id","endpoint_id","checkpoint","record_id","field","field_sha256",
                                     "verification","verification_sha256"))
            and Q(row["dose_exact"])==Q(c["drive"]["optical_dose_exact"])
            and Q(row["requested_drive_exact"])==Q(c["drive"]["activation_strength_exact"])
            and Q(row["solver_binary_drive_exact"])==Q(c["drive"]["activation_strength_binary_exact"])
            and Q(row["node_addition_gap_exact"])==Q(c["xi_dimensionless_exact"])
            and Q(row["observer_goal_exact"])==Q(c["e_goal_dimensionless_exact"])
            and c["observer_variant"]=="current_disk__spatial_max","same original case, exact drive and scoring authority")
    native=read(canonical(row["native_case"]),row["native_case_sha256"])
    phase=read(canonical(row["native_phase"]),row["native_phase_sha256"])
    require(all(row[key]==links[row[path]] for key,path in (("native_case_sha256","native_case"),
                                                           ("native_phase_sha256","native_phase")))
            and native["case_id"]==case and native["record_id"]==c["record_id"], "same reviewed native case")
    if case in rec:
        pinned=rec[case]
        produced=R/"actual"/case
        require(c==pinned["case"] and native["archived_field_sha256"]==c["field_sha256"]
                and native["old_native_signed_point_sha256"]==pinned["native_proof_sha256"]
                and native["original_cache_packet_sha256"]==pinned["cache_packet_sha256"],
                "same old case and signed anchor for verified recenter")
        files=phase["producer_metadata_and_packet_sha256"]
        cache=produced/"recentered_cache_packet.json"
        step=produced/"exact_binary_step.json"
        old=canonical(pinned["cache_packet"])
        digests=(files[rel(cache)],pinned["cache_packet_sha256"],files[rel(step)])
        require(digests==(row["operator_cache_packet_sha256"],native["original_cache_packet_sha256"],
                          native["exact_step_sha256"]),"current recentered packet and old packet hash")
        return c,row,native,(cache,old,step),digests,pinned
    if case in w:
        folder=WAVES/"actual"/phase["wave_id"]/case
    else:
        folder=M/"actual"/case
    cache=folder/"cache_packet.json"
    return c,row,native,(cache,cache,None),(row["operator_cache_packet_sha256"],)*2+(None,),None


def check_hashes(a, expected):
    require(set(a)==set(expected) and all(a_sha(a[k])==v for k,v in expected.items()),
            "native audited per-tensor bytes required")


def compare_roles(old,new,expected,c,row,witness,original):
    check_hashes(old,expected)
    require(set(old)==set(new) and all(a.dtype==new[k].dtype and a.shape==new[k].shape
            and (k=="u" or a_sha(a)==a_sha(new[k])) for k,a in old.items()),
            "verified recenter preserves all operator tensors")
    require(witness["case_id"]==c["case_id"] and witness["round_once_ties_to_even"] is True
            and witness["native_proof_sha256"]==row["original_unresolved_native_sha256"]
            and witness["u0_sha256"]==a_sha(old["u"])
            and witness["u1_sha256"]==a_sha(new["u"])==row["unsigned_binary_center_sha256"]
            and a_sha(old["u"])!=a_sha(new["u"]), "exact step and observed archived/new centers required")
    at,=[a for a in original["own_radius_attempts"] if Q(a["radius_exact"])==R36]
    forces=at["signed_center_force_intervals"]
    evidence=witness["ordered_rows"]["values"]
    deltas=[]
    for x,y,f,w in zip(old["u"],new["u"],forces,evidence,strict=True):
        before,after=map(T.rational,(x,y))
        mid=(Q(f["lower_exact"])+Q(f["upper_exact"]))/2
        target=before-mid
        delta=after-before
        require(np.float64(float(target)).tobytes()==y.tobytes() and
                list(map(str,(before,mid,target,after,delta)))==w,
                "single nearest-even binary step from the exact signed native midpoint")
        deltas.append(abs(delta))
    sigma=max(deltas)
    require(0<sigma==Q(witness["sigma_exact"])==Q(row["archived_field_center_shift_sigma_exact"]),
            "exact observer center shift")


def objects(case,authority,opened):
    X.bootstrap()
    c,row,native,paths,hashes,rec=case_authority(case,authority)
    current,old_p,step_p=paths
    current_h,old_h,step_h=hashes
    archived=read(canonical(c["field"]),c["field_sha256"],opened,"archived observer geometry and parameters")
    read(canonical(c["verification"]),c["verification_sha256"],opened,"archived field verification")
    old=read(old_p,old_h,opened,"independently verified original cache")
    packet=old if rec is None else read(current,current_h,opened,"independently verified recentered cache")
    require(old.get("case_id",old.get("pilot_id"))==case and old["record_id"]==c["record_id"]
            and old["input_state"]==c["field"] and old["input_state_sha256"]==c["field_sha256"]
            and old["field_verification_sha256"]==c["verification_sha256"]
            and old["free_dof"]==611 and old["integration_policy"]["volume_duffy_gauss_order"]==8
            and old["integration_policy"]["surface_gauss_order"]==12
            and old["status"]=="FROZEN_FRESH_EXACT_BINARY_QUADRATURE_CACHE",
            "authoritative original Q8/S12 field/cache")
    base,dose,drive=X.exact_parameters(old,c)
    cell={k:(base[k],base[k]) for k in T.ORDER}
    require(doc(cell)==c["closed_parameter_cell_exact"],"exact singleton cell match")
    u0=T.Target(old)
    X.validate_cache_shapes(u0.a)
    X.verify_state(u0,archived,old,c,base,drive)
    for header,tgt in ((old,u0),):
        require(header["frozen_binary_cache"]["array_shapes"]=={k:list(v.shape) for k,v in tgt.a.items()},
                "exact shape header match")
    if rec is None:
        target=u0
        expected=native.get("per_tensor_sha256",native.get("cache_tensor_content_sha256"))
        if expected is None:
            path=canonical(native["reused_signed_cached_anchor_file"])
            anchor=read(path,native["reused_signed_cached_anchor_sha256"])
            require(anchor["cache_packet_sha256"]==old_h,"same cached M signed anchor")
            expected=anchor["cache_tensor_content_sha256"]
        check_hashes(target.a,expected)
        require(Q(row["archived_field_center_shift_sigma_exact"])==0,"unmoved observer center")
    else:
        target=T.Target(packet)
        X.validate_cache_shapes(target.a)
        require(packet["frozen_binary_cache"]["array_shapes"]==old["frozen_binary_cache"]["array_shapes"]
                and set(packet)==set(old)|{"recenter"}
                and all(packet[k]==v for k,v in old.items() if k not in ("status","frozen_binary_cache"))
                and packet["status"]=="FROZEN_SAME_CACHED_OPERATOR_AT_NEW_PROOF_CENTER"
                and packet["recenter"]["exact_step_sha256"]==step_h
                and packet["recenter"]["original_point_cache_sha256"]==old_h,
                "same exact operator and score archive header")
        proof=read(canonical(rec["native_proof"]),rec["native_proof_sha256"])
        witness=read(step_p,step_h,opened,"independently verified exact binary step")
        compare_roles(u0.a,target.a,rec["per_tensor_sha256"],c,row,witness,proof)
    require(a_sha(target.a["preconditioner"])==row["binary_R_sha256"]
            and a_sha(target.a["u"])==row["unsigned_binary_center_sha256"],
            "exact current R and center hashes required")
    radius,stage=scope_for(Q(row["cached_proof_radius_exact"]),Q(row["cached_ZQ_upper_exact"]))
    names={k:row[k] for k in ("native_case","native_case_sha256","native_phase","native_phase_sha256",
                             "native_observed_os_exit_sha256","native_terminal_sha256")}
    flags={**NO_MATH,"root_claimed":False,"whole_box_spatial_geometry_verified":False}
    contract={"status":"FROZEN_EXPLICIT_REMAINING32_SINGLETON","record_id":"transfer_exact_"+case,
        "pilot_id":case,"physical_source_record_id":c["record_id"],"archived_field":c["field"],
        "archived_field_sha256":c["field_sha256"],"archived_observer_variant":c["observer_variant"],
        "archived_observer_allowance_um_exact":c["frozen_conditional_allowance_um_exact"],
        "archived_verification":c["verification"],"archived_verification_sha256":c["verification_sha256"],
        "observer_goal_exact":row["observer_goal_exact"],"archived_addition_gap_exact":row["node_addition_gap_exact"],
        "archived_shift_exact":row["archived_field_center_shift_sigma_exact"],"verification_center_kind":row["spatial_verification_center"],
        "center_binary_sha256":row["unsigned_binary_center_sha256"],"original_point_cache":rel(old_p),
        "original_point_cache_sha256":old_h,"point_cache":rel(current),"point_cache_sha256":current_h,
        "frozen_binary_cache_sha256":packet["frozen_binary_cache"]["sha256"],"binary_R_sha256":row["binary_R_sha256"],
        "recenter_exact_step_sha256":step_h,"exact_requested_drive":c["drive"],
        "material_parameters_binary_exact":old["material_parameter_storage_exact"],"parameter_cell":doc(cell),
        "parameter_cell_sha256":E._digest(doc(cell)),"energy_scale_binary_exact":str(target.energy_scale),
        "source_cached_radius_exact":row["cached_proof_radius_exact"],"native_cached_ZQ_upper_exact":row["cached_ZQ_upper_exact"],
        "planned_geometry_radius_exact":str(radius),"derivative_policy":stage,
        "matched_cached_ZQ_available":stage.startswith("EXISTING_MATCHED"),**names,**flags}
    predictor=np.zeros((611,5),dtype="<f8")
    predictor.flags.writeable=False
    for a in target.a.values():
        a.flags.writeable=False
    return target,predictor,base,dose,drive,radius,contract

def context(parts, paths):
    target,S,base,dose,drive,radius,contract=parts
    c={k:(v,v) for k,v in base.items()}
    row={"record_id":contract["record_id"],"source_record_id":contract["physical_source_record_id"],
         "dose_exact":str(dose),"parameter_order":list(T.ORDER),"parameter_cell":contract["parameter_cell"],
         "parameter_cell_sha256":contract["parameter_cell_sha256"],
         "base_parameters_exact":{k:str(v) for k,v in base.items()},
         "point_cache_file":contract["point_cache"],"point_cache_sha256":contract["point_cache_sha256"],
         "frozen_binary_cache_sha256":contract["frozen_binary_cache_sha256"],
         "model_observation_contract_sha256":paths["model_contract_sha256"],"predictor_sha256":paths["predictor_sha256"]}
    ctx=X.SingletonContext(contract["pilot_id"],contract["record_id"],row,target,S,c,base,dose,drive)
    check_context(ctx)
    return ctx


def check_context(ctx):
    require(type(ctx) is X.SingletonContext and ctx.pilot_id in IDS and ctx.record_id=="transfer_exact_"+ctx.pilot_id
            and ctx.S.shape==(611,5) and np.all(ctx.S==0)
            and set(ctx.cell)==set(ctx.base)==set(T.ORDER)
            and all(ctx.cell[k]==(ctx.base[k],ctx.base[k]) for k in T.ORDER)
            and ctx.dose*ctx.base["common_gain"]==ctx.drive==Q(ctx.target.packet["activation_strength_exact"])
            and ctx.target.packet.get("case_id",ctx.target.packet.get("pilot_id"))==ctx.pilot_id
            and ctx.target.packet["record_id"]==ctx.row["source_record_id"]
            and ctx.target.packet["frozen_binary_cache"]["sha256"]==ctx.row["frozen_binary_cache_sha256"],
            "explicit same-field exact singleton context")
    return E._context_identity(ctx)


def node_box(ctx, radius):
    identity=check_context(ctx)
    require(type(radius) is Fraction and radius in (R34,R36),"explicit planned positive exact radius required")
    Xs=[list(map(T.rational,a)) for a in ctx.target.a["reference"]]
    Y=ctx.target.nodes()
    spread=[[Fraction(0),Fraction(0)] for _ in Y]
    for k in ctx.free:
        n,comp=divmod(int(k),2)
        spread[n][comp]=radius
    require(len(Y)==325 and sum(x==radius for row in spread for x in row)==611,"complete 650 exact-coordinate hull")
    return doc({**identity,"reference_nodes_exact":Xs,"current_centers_exact":Y,
        "coordinate_radii_exact":spread,"coordinate_hulls_exact":[
        [[v-e,v+e] for v,e in zip(row,errs,strict=True)] for row,errs in zip(Y,spread,strict=True)],
        "predictor_hulls_exact":[[Fraction(0),Fraction(0)] for _ in range(611)],
        "triangles":ctx.target.a["triangles"].tolist(),"free":ctx.free.tolist(),
        "corrector_radius_exact":str(radius),"root_membership_assumed":False,
        "any_spatial_derivative_bound_assumed":False})


def routes(ctx,rule):
    identity=check_context(ctx)
    target=ctx.target
    a=target.a
    def projection(ids):
        dofs=[2*int(i)+k for i in ids for k in (0,1)]
        selected=[(j,target.index[dof]) for j,dof in enumerate(dofs) if dof in target.index]
        require(bool(selected) and len(set(g for _,g in selected))==len(selected),"valid projected local map")
        return {"free_local":[j for j,_ in selected],"global_indices":[g for _,g in selected]}
    bulk=[]
    for i,ids in enumerate(a["triangles"].tolist()):
        n=sum(T.rational(a["reference"][j][0])==0 for j in ids[:3])
        require(n in (0,1,2),"positive material radial span")
        bulk.append({"element":i,"nodes":ids,"chart":("off_axis","axis_vertex","axis_edge_physical")[n],
                     "active":str(T.rational(a["active"][i,0])),**projection(ids),
                     "auxiliary_element_rule_sha256":E._bulk_rule(target,rule,i)})
    surfaces=[]
    for face in (0,1):
        for edge,ids in enumerate(a[f"surface_{face}_edges"].tolist()):
            E._surface_rule(target,rule,face,edge)
            surfaces.append({"face":face,"edge":edge,"nodes":ids,**projection(ids),
                "Rlo_exact":str(T.rational(a["reference"][ids[0]][0])),
                "Rhi_exact":str(T.rational(a["reference"][ids[-1]][0])),
                "trace_rule_sha256":E._digest(doc([rule["surface_s"][face],rule["surface_weights"][face][edge]]))})
    require(Counter(r["chart"] for r in bulk)==X.CHARTS and len(surfaces)==24
            and sum(r["active"]=="1" for r in bulk)==24
            and {i for row in bulk for i in row["global_indices"]}==set(range(611)),
            "all 144 bulk, 24 trace routes and free coordinates required")
    return {"status":"COLD_EXACT_REFERENCE_ROUTE_METADATA_AWAITS_GEOMETRY",
            "identity":identity,"bulk":bulk,"surface":surfaces,"bulk_chart_counts":X.CHARTS,
            "derived_from_frozen_reference_and_exact_rule":True,"positive_deformed_geometry_proved":False}


def predictor_doc(parts):
    contract=parts[-1]
    return {"status":"ZERO_AFFINE_PREDICTOR_FOR_FIXED_SINGLETON","record_id":contract["record_id"],
            "parameter_cell_sha256":contract["parameter_cell_sha256"],"shape":[611,5],
            "little_endian_zero_binary64_sha256":a_sha(parts[1]),"identically_zero":True}


def calculations(ctx,radius):
    rule=ctx.target.auxiliary_rule()
    return {"row":ctx.row,"rule":doc(rule),"routes":routes(ctx,rule),"box":node_box(ctx,radius)}


def store_case(folder,parts):
    paths={}
    def add(k,value):
        p=folder/(k+".json")
        paths[k]=rel(p)
        paths[k+"_sha256"]=X.save(p,value)
    add("model_contract",parts[-1])
    add("predictor",predictor_doc(parts))
    ctx=context(parts,paths)
    for k,v in calculations(ctx,parts[-2]).items():
        add(k,v)
    return {"case_id":ctx.pilot_id,"record_id":ctx.record_id,
            "planned_geometry_radius_exact":str(parts[-2]),"derivative_policy":parts[-1]["derivative_policy"],**paths}


def open_context(case_id,manifest_sha256):
    """Read one current context after an external prep manifest pin. No launch approval."""
    manifest=read(PREP/"manifest.json",manifest_sha256)
    require(manifest["status"]==STATUS and manifest["source_sha256"]==sha(SOURCE)
            and tuple(e["case_id"] for e in manifest["targets"])==IDS,
            "exact current cold manifest and complete ordered scope")
    authority=load_authority()
    e=pick(manifest["targets"],case_id)
    opened={}
    with X.no_integrals():
        parts=objects(case_id,authority,opened)
        recorded=manifest["opened_actual_input_files"]
        require(all(recorded[k]==v for k,v in opened.items()),"same proved source and payload files")
        ctx=context(parts,e)
        expected={"model_contract":parts[-1],"predictor":predictor_doc(parts),**calculations(ctx,parts[-2])}
        for k,v in expected.items():
            require(v==read(canonical(e[k]),e[k+"_sha256"]),"regenerated exact case context mismatch")
    return ctx,ctx.target.auxiliary_rule(),e


def reject(f):
    try:
        f()
    except (ValueError,KeyError):
        return True
    return False


def fixtures(authority):
    first=IDS[0]
    a,b=IDS[:2]
    rows,w,m,rec,links=authority
    changed=dict(rows)
    changed[first]=dict(rows[first],requested_drive_exact="0")
    bad=(changed,w,m,rec,links)
    mixed=dict(rows)
    mixed[a]=dict(rows[b])
    mix=(mixed,w,m,rec,links)
    assert reject(lambda:case_authority(first,bad)) and reject(lambda:case_authority(a,mix))
    copies=np.array([1.],dtype="<f8")
    zero=np.array([0.],dtype="<f8")
    assert reject(lambda:check_hashes({"u":zero},{"u":a_sha(copies)}))
    assert scope_for(R34,Fraction(1,3))[0]==R34 and scope_for(R36,Fraction(1,16))[0]==R36
    assert scope_for(R32,Fraction(7,8))[1]=="FRESH_VARIATION_REQUIRED_BEFORE_SPATIAL_LAUNCH"
    assert reject(lambda:scope_for(R32,Fraction(1))) and reject(lambda:case_authority("foreign_case",authority))
    return {"foreign_case_wrong_exact_drive_and_changed_center_rejected":True,
            "mixed_proof_source_rejected":True,"no_enlargement_of_native_narrow_cube":True,
            "containing_cube_loose_bound_requires_fresh_plan":True}


def prepare():
    authority=load_authority()
    checks=fixtures(authority)
    PREP.mkdir(parents=True,exist_ok=False)
    closure=X.local_source_closure()
    protocol={"status":"FROZEN_REMAINING32_CONTEXT_PREPARATION","source_sha256":sha(SOURCE),
              "ordered_cases":list(IDS),"input_metadata_and_source_sha256":authority[-1],
              "local_source_closure_sha256":{**closure,rel(SOURCE):sha(SOURCE)},
              "geometry_or_variation_actual_launch_authorized":False,"full_five_gate_complete":False,**NO_MATH}
    ph=X.save(PREP/"protocol.json",protocol)
    opened={}
    targets=[]
    with X.no_integrals():
        for case in IDS:
            parts=objects(case,authority,opened)
            folder=PREP/case
            folder.mkdir(exist_ok=False)
            targets.append(store_case(folder,parts))
            print(case,flush=True)
    require(len(opened)==100 and sum(p["purpose"].endswith("cache") for p in opened.values())==34,
            "all 64 original field/verification, 34 caches and two binary step files opened")
    counts=Counter(r["derivative_policy"] for r in targets)
    require(counts=={"EXISTING_MATCHED_WHOLE_CUBE_CACHED_ZQ":22,
                     "FRESH_VARIATION_REQUIRED_BEFORE_SPATIAL_LAUNCH":7,
                     "EXISTING_MATCHED_SMALL_CUBE_LARGER_TRIAL_REQUIRES_NEW_AD":3},
            "22 matched large, 7 planned fresh variation, 3 true narrow proofs")
    manifest={"status":STATUS,"source_sha256":sha(SOURCE),"protocol_sha256":ph,
              "targets":targets,"derivative_policy_counts":dict(counts),"opened_actual_input_files":opened,
              "ranges":"650 exact hulls per case with 611 free radii and zero predictor",
              "geometry_or_variation_actual_launch_authorized":False,
              "spatial_bridge_actual_launch_authorized":False,"full_five_gate_complete":False,**NO_MATH}
    mh=X.save(PREP/"manifest.json",manifest)
    pf={"status":"PASS","checks":{**checks,"all_32_exact_cached_and_archive_contexts":True,
         "no_field_scores_replaced":True,"all_locations_and_explicit_source_hashes_bound":True,
         "complete_box_rules_projection_and_zero_predictor":True},
         "source_sha256":sha(SOURCE),"protocol_sha256":ph,"manifest_sha256":mh,**NO_MATH}
    pfh=X.save(PREP/"preflight.json",pf)
    (PREP/"transfer_ambiguity_remaining32_exact_targets_frozen.py").write_bytes(SOURCE.read_bytes())
    print(doc({"source_sha256":sha(SOURCE),"protocol_sha256":ph,"manifest_sha256":mh,"preflight_sha256":pfh}))


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("mode",choices=("prepare",))
    cli.parse_args()
    prepare()
