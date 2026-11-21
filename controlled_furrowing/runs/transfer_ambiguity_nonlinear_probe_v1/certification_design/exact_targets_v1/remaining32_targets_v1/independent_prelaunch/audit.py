"""Cold independent review with direct binary cache parsing and exact coordinate reconstruction."""
import base64
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import io
import json
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE.parent
EX=OUT.parent
ROOT=next(p for p in EX.parents if (p/"PLAN.md").is_file())
PREP=OUT/"preparation"
SOURCE=ROOT/"scripts/transfer_ambiguity_remaining32_exact_targets.py"
CENSUS=EX/"remaining32_bridge_design_v1/manifest.json"
REQUIRED={
    "source_sha256": (SOURCE,"d1a9c9e4839ff85a9926b6fd5ecc11e4cd8e0831c7e6263ecfbe068ef22ea3f4"),
    "protocol_sha256": (PREP/"protocol.json","a12669fd1e50cbd7435d016ba60223791c8a70119f622048c28e24ef5a9f456b"),
    "manifest_sha256": (PREP/"manifest.json","dfea134640f99b45e9ca75d1ce718043d67f2917375506997d2465abda8c1242"),
    "preflight_sha256": (PREP/"preflight.json","86ddb06c75c2deedea29d62bb7e03da6142d55b41a0291e1567c9e4dba6aab4f"),
}
META={
    CENSUS:"5f043b6a2fac20bdbb02dd3a5c9e58a860b2a4453b246373e28178f99ab49d53",
    OUT/"context_readback.json":"f4f4b3366a2bd412d582deaa9bc8958e81c401cb7621dade4a427dd654316649",
}
SIDES=("max_contrast_change","max_negative_contrast_change")
IDS=sorted(f"{sign}_small_{i:02d}" for sign in SIDES for i in range(1,17))
RECENT={"max_contrast_change_small_13","max_negative_contrast_change_small_16"}
ORDER=("apical_area_modulus","basal_area_modulus","apical_surface_tension","basal_surface_tension","common_gain")
R36,R34,R32=(Q(1,2**n) for n in (36,34,32))
MATCH="EXISTING_MATCHED_WHOLE_CUBE_CACHED_ZQ"
FRESH="FRESH_VARIATION_REQUIRED_BEFORE_SPATIAL_LAUNCH"
SMALL="EXISTING_MATCHED_SMALL_CUBE_LARGER_TRIAL_REQUIRES_NEW_AD"
NO_MATH=("new_geometry_certificates","new_AD_derivatives","new_integrations","new_solves")
SHAPES={"reference":(325,2),"free":(611,),"u":(611,),"triangles":(144,6),
        "shape":(64,6),"shape_grad":(144,64,6,2),"qr":(144,64),"volume_weights":(144,64),
        "active":(144,1),"energy_scale":(1,),"preconditioner":(611,611)}
for i in range(2):
    for k,v in {"edges":(12,3),"N":(12,3),"dN":(12,3),"r0":(12,12),"l0":(12,),"weights":(12,12)}.items():
        SHAPES[f"surface_{i}_{k}"]=v
LINKS={}


def verify(ok,why):
    if not ok:
        raise ValueError(why)


def digest(raw):
    return sha256(raw).hexdigest()


def file(path,want,parse=True):
    p=(ROOT/path).resolve()
    name=p.relative_to(ROOT).as_posix()
    raw=p.read_bytes()
    g=digest(raw)
    verify(g==want==LINKS.get(name,g),f"changed input {name}")
    LINKS[name]=g
    return json.loads(raw) if parse else raw


def index_doc(ob,name,k,h):
    return file(ob[k],ob[h]) if name else file(ob[k],ob[h],False)


def aq(a):
    return digest(a.tobytes(order="C"))


def jsha(x):
    return digest(json.dumps(x,sort_keys=True,separators=(",",":"),ensure_ascii=True).encode())


def qtable(array):
    return [[Q(float(v)) for v in row] for row in array]


def require_zero(x):
    return type(x) is int and x==0


def unpack(packet):
    desc=packet["frozen_binary_cache"]
    raw=base64.b64decode(desc["data"],validate=True)
    verify(digest(raw)==desc["sha256"],"embedded byte digest does not match file")
    with np.load(io.BytesIO(raw),allow_pickle=False) as handle:
        verify(set(handle.files)==set(SHAPES) and len(handle.files)==len(SHAPES),"binary cache tensor names fail")
        arr={name:handle[name] for name in handle.files}
    verify(desc["array_shapes"]=={name:list(a.shape) for name,a in arr.items()},"binary array shape metadata drift")
    for key,a in arr.items():
        verify(a.shape==SHAPES[key],"binary cache exact dimensions fail")
        if key in ("free","triangles","surface_0_edges","surface_1_edges"):
            verify(a.dtype.kind in ("i","u"),"integer topology tensor required")
        else:
            verify(a.dtype.str=="<f8" and np.isfinite(a).all(),"finite little-endian binary64 tensor required")
    return arr


def radius_native(native,c,r):
    st=native["status"]
    if "own_radius_attempts" in native:
        verify(st=="PASS_INDEPENDENT_CACHED_LOCAL_ROOT" and Q(native["chosen_radius_exact"])==r,
               "no matched passing original wave box")
        leaf,=[x for x in native["own_radius_attempts"] if Q(x["radius_exact"])==r]
        verify(leaf["whole_cube_AD_recomputed_for_this_radius"] is True,"whole-box bound missing")
    else:
        verify(st in ("PASS_INDEPENDENT_M_CACHED_LOCAL_ROOT","PASS_INDEPENDENT_NEGATIVE_M_SMALLER_CACHED_ROOT",
                      "PASS_INDEPENDENT_RECENTERED_CACHED_LOCAL_ROOT_AND_OLD_OBSERVER_BUDGET")
               and Q(native["radius_exact"])==r,"no matched passing native radius")
        leaf=native
    eta=Q(leaf.get("eta_upper_exact",leaf.get("independent_eta_upper_exact")))
    z=Q(leaf.get("ZQ_upper_exact",leaf.get("independent_ZQ_upper_exact")))
    lhs=Q(leaf["strict_selfmap_lhs_upper_exact"])
    error=Q(leaf.get("posterior_error_upper_exact",leaf.get("independent_posterior_error_upper_exact")))
    verify(lhs==eta+z*r<r and eta>=0 and 0<=z<1 and error==eta/(1-z)<r,
           "case authority does not prove claimed cached fixed-point bound")
    verify((eta,z,error)==tuple(Q(c[k]) for k in ("cached_eta_upper_exact","cached_ZQ_upper_exact","cached_posterior_upper_exact")),
           "native exact bound differs from census")
    return z,error


def reconstruct_step(contract, native, census, authority, old, new, pk0, pk1):
    name=contract["pilot_id"]
    verify(set(pk1)==set(pk0)|{"recenter"} and pk1["status"]=="FROZEN_SAME_CACHED_OPERATOR_AT_NEW_PROOF_CENTER"
           and all(pk0[k]==pk1[k] for k in pk0 if k not in ("status","frozen_binary_cache")),
           "recenter altered a physical header")
    verify(set(old)==set(new) and all(old[k].shape==new[k].shape and old[k].dtype==new[k].dtype for k in old)
           and {k for k in old if aq(old[k])!=aq(new[k])}=={"u"},
           "exactly one array is allowed to differ in a recentered cache")
    phase=file(contract["native_phase"],contract["native_phase_sha256"])
    rp=phase["producer_metadata_and_packet_sha256"]
    chosen=contract["recenter_exact_step_sha256"]
    wp,=[name for name,h in rp.items() if h==chosen and name.endswith("/exact_binary_step.json")]
    w=file(wp,chosen)
    old_native=native["old_native_signed_point_sha256"]
    loc,=[name for name,h in authority.items() if h==old_native]
    orig=file(loc,old_native)
    a,=[r for r in orig["own_radius_attempts"] if Q(r["radius_exact"])==R36]
    corners=a["signed_center_force_intervals"]
    verify(len(corners)==len(w["ordered_rows"]["values"])==len(old["u"])==611
           and w["case_id"]==name and w["round_once_ties_to_even"] is True
           and w["native_proof_sha256"]==old_native and orig["chosen_radius_exact"] is None
           and pk1["recenter"]["exact_step_sha256"]==chosen
           and pk1["recenter"]["original_point_cache_sha256"]==contract["original_point_cache_sha256"]
           and w["u0_sha256"]==aq(old["u"]) and w["u1_sha256"]==aq(new["u"]),
           "recenter step provenance is not the archived unresolved point")
    shifts=[]
    for x,y,c,row in zip(old["u"],new["u"],corners,w["ordered_rows"]["values"],strict=True):
        old_x, new_x=Q(float(x)),Q(float(y))
        mid=(Q(c["lower_exact"])+Q(c["upper_exact"]))/2
        val=old_x-mid
        shift=new_x-old_x
        verify(np.float64(float(val)).tobytes()==y.tobytes()
               and list(map(Q,row))==[old_x,mid,val,new_x,shift],"bad archived signed midpoint binary step")
        shifts.append(abs(shift))
    sigma=max(shifts)
    verify(sigma>0 and sigma==Q(w["sigma_exact"])==Q(contract["archived_shift_exact"])
           ==Q(census["archived_field_center_shift_sigma_exact"])==Q(native["actual_binary_center_shift_sigma_exact"]),
           "recenter observer transport distance drift")
    return sigma


def rules_and_routes(A,ref,tri,map_free,identity,rule,routes,scale,mask,top_z):
    bary=[]
    for row in qtable(A["shape"]):
        r=row[1]+(row[3]+row[4])/2
        z=row[2]+(row[4]+row[5])/2
        bary.append([str(1-r-z),str(r),str(z)])
    verify(all(min(map(Q,row))>0 for row in bary) and bary==rule["volume_barycentric"],
           "expected 64 exact interior binary-reference auxiliary sites")
    vw=[[str(v/scale) for v in row] for row in qtable(A["volume_weights"])]
    sw=[]
    ss=[]
    for f in range(2):
        sw.append([[str(v/scale) for v in row] for row in qtable(A[f"surface_{f}_weights"])])
        ss.append([str(row[2]+row[1]/2) for row in qtable(A[f"surface_{f}_N"])])
    verify(vw==rule["volume_weights"] and sw==rule["surface_weights"] and ss==rule["surface_s"]
           and all(Q(v)>0 for row in vw for v in row) and all(Q(v)>0 for w in sw for row in w for v in row)
           and all(0<Q(v)<1 for row in ss for v in row),"same exact normalized rule required at all 9504 sites")
    def projection(ids):
        candidates=[2*i+j for i in ids for j in (0,1)]
        loc=[i for i,n in enumerate(candidates) if n in map_free]
        glob=[map_free[candidates[i]] for i in loc]
        verify(loc and len(glob)==len(set(glob)),"invalid local projection")
        return {"free_local":loc,"global_indices":glob}
    bulk=[]
    chart=("off_axis","axis_vertex","axis_edge_physical")
    for i,ids in enumerate(tri):
        axis=sum(ref[j][0]==0 for j in ids[:3])
        verify(axis<3,"nonpositive reference radial measure")
        bulk.append({"element":i,"nodes":ids,"chart":chart[axis],"active":str(mask[i]),
                     **projection(ids),"auxiliary_element_rule_sha256":jsha([bary,vw[i]])})
    surfaces=[]
    for f,zref in enumerate((Q(0),top_z)):
        ref_edges={}
        for a,b,c,ab,bc,ca in tri:
            for u,v,m in ((a,b,ab),(b,c,bc),(c,a,ca)):
                if ref[u][1]==ref[v][1]==zref:
                    u,v=sorted((u,v),key=lambda i:ref[i][0])
                    verify(ref[u][0]<ref[v][0] and ref_edges.get((u,v),m)==m,"nonconforming trace reference edge")
                    ref_edges[(u,v)]=m
        ordered=sorted(([u,m,v] for (u,v),m in ref_edges.items()),key=lambda edge:ref[edge[0]][0])
        verify(len(ordered)==12 and ordered==A[f"surface_{f}_edges"].tolist()
               and all(x[-1]==y[0] for x,y in zip(ordered,ordered[1:])),"surface route not derived from its triangle trace")
        for i,ids in enumerate(ordered):
            surfaces.append({"face":f,"edge":i,"nodes":ids,**projection(ids),
                             "Rlo_exact":str(ref[ids[0]][0]),"Rhi_exact":str(ref[ids[-1]][0]),
                             "trace_rule_sha256":jsha([ss[f],sw[f][i]])})
    counts=dict(Counter(row["chart"] for row in bulk))
    verify(counts=={"axis_vertex":6,"axis_edge_physical":6,"off_axis":132}
           and routes["bulk_chart_counts"]==counts and routes["identity"]==identity and routes["bulk"]==bulk
           and routes["surface"]==surfaces and routes["positive_deformed_geometry_proved"] is False
           and set.union(*(set(b["global_indices"]) for b in bulk))==set(range(611)),
           "stored full projections and reference charts disagree with independent reconstruction")


def check_one(e,c,manifest,premises):
    name=e["case_id"]
    six=("model_contract","row","predictor","box","rule","routes")
    docs={}
    for k in six:
        location=PREP/name/(k+".json")
        verify((ROOT/e[k]).resolve()==location,"foreign context artifact path")
        docs[k]=file(location,e[k+"_sha256"])
    con,row,pre,box,rule,routes=(docs[k] for k in six)
    verify(con["pilot_id"]==name==c["case_id"] and e["record_id"]==con["record_id"]=="transfer_exact_"+name
           and row["record_id"]==con["record_id"] and con["physical_source_record_id"]==c["record_id"]
           and all(require_zero(con[k]) for k in NO_MATH) and con["root_claimed"] is False
           and con["whole_box_spatial_geometry_verified"] is False,"different physical state or premature proof claim")
    for a,b in (("archived_field","field"),("archived_verification","verification"),("point_cache_sha256","operator_cache_packet_sha256"),
                ("binary_R_sha256","binary_R_sha256"),("observer_goal_exact","observer_goal_exact"),
                ("archived_addition_gap_exact","node_addition_gap_exact"),
                ("native_case","native_case"),("native_case_sha256","native_case_sha256"),
                ("native_phase","native_phase"),("native_phase_sha256","native_phase_sha256"),
                ("native_terminal_sha256","native_terminal_sha256"),("native_observed_os_exit_sha256","native_observed_os_exit_sha256")):
        verify(con[a]==c[b],f"census identity drift {a}")
    for path_key,hash_key in (("archived_field","archived_field_sha256"),("archived_verification","archived_verification_sha256"),
                             ("point_cache","point_cache_sha256"),("original_point_cache","original_point_cache_sha256")):
        verify(manifest["opened_actual_input_files"][con[path_key]]["sha256"]==con[hash_key],"missing admitted physical artifact")
    verify(con["archived_field_sha256"]==c["field_sha256"] and con["archived_verification_sha256"]==c["verification_sha256"],
           "source field/verification digest drift")
    original=file(con["archived_field"],con["archived_field_sha256"])
    file(con["archived_verification"],con["archived_verification_sha256"],False)
    native=file(con["native_case"],con["native_case_sha256"])
    verify(native["case_id"]==name and native["record_id"]==c["record_id"],"native proof from foreign case")
    oldp=file(con["original_point_cache"],con["original_point_cache_sha256"])
    cached=file(con["point_cache"],con["point_cache_sha256"])
    old=unpack(oldp)
    A=unpack(cached) if name in RECENT else old
    mat=oldp["parameters"]
    verify(original["parameters"]==mat and oldp["input_state_sha256"]==con["archived_field_sha256"]
           and oldp["field_verification_sha256"]==con["archived_verification_sha256"]
           and oldp["input_state"]==con["archived_field"] and oldp["record_id"]==c["record_id"]
           and oldp["case_id"]==name and Q(mat["mu"])==1 and Q(mat["bulk"])==5000
           and oldp["integration_policy"]["volume_duffy_gauss_order"]==8
           and oldp["integration_policy"]["surface_gauss_order"]==12,
           "wrong physical parameters or original quadrature state")
    ph=con["exact_requested_drive"]
    verify(all(oldp[k]==v==cached[k] for k,v in ph.items()) and con["material_parameters_binary_exact"]==oldp["material_parameter_storage_exact"],
           "independent unchanged exact drive and binary material storage required")
    dose,gain,strength,rounded=map(lambda k:Q(ph[k]),("optical_dose_exact","common_gain_exact","activation_strength_exact","activation_strength_binary_exact"))
    verify(0<dose==Q(c["dose_exact"]) and gain>0 and dose*gain==strength==Q(c["requested_drive_exact"])
           and strength!=rounded==Q(c["solver_binary_drive_exact"])==Q(original["activation_strength"]),
           "rounded solver activation substituted for exact requested drive")
    expected={k:Q(mat[k]) for k in ORDER[:-1]}|{"common_gain":gain}
    cell={k:[str(v),str(v)] for k,v in expected.items()}
    verify(con["parameter_cell"]==row["parameter_cell"]==cell and con["parameter_cell_sha256"]==jsha(cell)==row["parameter_cell_sha256"]
           and row["parameter_order"]==list(ORDER) and row["base_parameters_exact"]=={k:str(v) for k,v in expected.items()},
           "exact material singleton identity mismatch")
    scale=Q(float(A["energy_scale"][0]))
    verify(scale>0 and scale==Q(original["reference_volume"])==Q(oldp["normalization_binary_exact"])==Q(con["energy_scale_binary_exact"])
           and scale==Q(float(A["volume_weights"].sum())),"binary reference volume normalization drift")
    r=Q(con["source_cached_radius_exact"])
    verify(r==Q(c["cached_proof_radius_exact"]),"native radius drift")
    z,error=radius_native(native,c,r)
    verify(z==Q(con["native_cached_ZQ_upper_exact"]),"native derivative differs")
    expected_policy={R36:(R36,SMALL),R34:(R34,MATCH),R32:(R34,FRESH)}
    verify(r in expected_policy and (r!=R32 or z>Q(4,5)),"not an admitted planned radius")
    rad,policy=expected_policy[r]
    verify(Q(e["planned_geometry_radius_exact"])==Q(con["planned_geometry_radius_exact"])==rad
           and e["derivative_policy"]==con["derivative_policy"]==policy
           and con["matched_cached_ZQ_available"] is (r==rad),"mislabeled whole-box cached eligibility")
    if name in RECENT:
        sigma=reconstruct_step(con,native,c,premises,old,A,oldp,cached)
    else:
        sigma=Q(0)
        verify(con["point_cache_sha256"]==con["original_point_cache_sha256"]
               and con["recenter_exact_step_sha256"] is None and Q(con["archived_shift_exact"])==0,
               "nonrecentered case has an unapproved source change")
        tensor_hashes=native.get("per_tensor_sha256",native.get("cache_tensor_content_sha256"))
        if tensor_hashes is None:
            prior=file(native["reused_signed_cached_anchor_file"],native["reused_signed_cached_anchor_sha256"])
            tensor_hashes=prior["cache_tensor_content_sha256"]
        verify(tensor_hashes=={k:aq(a) for k,a in A.items()},"root proof and current cache tensor bytes differ")
    verify(aq(A["u"])==con["center_binary_sha256"]==c["unsigned_binary_center_sha256"]
           and aq(A["preconditioner"])==con["binary_R_sha256"]==c["binary_R_sha256"]
           and cached["frozen_binary_cache"]["sha256"]==con["frozen_binary_cache_sha256"]==row["frozen_binary_cache_sha256"]
           and sigma+error<Q(con["observer_goal_exact"]),"recenter transport, current R, center or cache SHA drift")
    ref=qtable(A["reference"])
    tri=A["triangles"].tolist()
    verify(ref==[list(map(Q,v)) for v in original["reference_nodes"]] and tri==original["triangles"]
           and old["u"].tolist()==original["free_displacement"],"old/new field topology or saved u0 mismatch")
    rim,height,patch,z0,z1=map(lambda k:Q(mat[k]),("outer_radius","height","patch_radius","activation_z_min","activation_z_max"))
    verify(rim==3 and height==1 and Q(mat["reference_curvature"])==0 and 0<patch<rim and 0<=z0<z1<=height,
           "unexpected reference body or active patch")
    mask=[]
    for a,b,g,*_ in tri:
        pts=[ref[a],ref[b],ref[g]]
        cr,cz=(sum(p[j] for p in pts)/3 for j in (0,1))
        active=int(cr<patch and z0<cz<z1)
        if active:
            verify(all(x<=patch and z0<=y<=z1 for x,y in pts),"active reference triangle crosses patch boundary")
        mask.append(active)
    verify(mask==A["active"][:,0].tolist()==original["activation_elements"] and sum(mask)==24,"physical active mask differs")
    free=[2*i+j for i,(r0,_) in enumerate(ref) for j in (0,1) if r0!=rim and (r0!=0 or j==1)]
    verify(free==A["free"].tolist()==box["free"] and len(free)==611 and tri==box["triangles"],
           "wrong exact reduced physical coordinate map")
    m={i:j for j,i in enumerate(free)}
    fixed=dict.fromkeys(set(range(650))-set(m),Q(0))
    disp0=fixed|dict(zip(free,map(Q,original["free_displacement"]),strict=True))
    disp1=fixed|dict(zip(free,(Q(float(v)) for v in A["u"]),strict=True))
    oldY=[[ref[i][j]+disp0[2*i+j] for j in range(2)] for i in range(325)]
    stored=[[Q(v) for v in n] for n in original["current_nodes"]]
    xi=max(abs(x-y) for a,b in zip(oldY,stored,strict=True) for x,y in zip(a,b,strict=True))
    verify([[Q(float(v)) for v in n] for n in oldY]==stored and xi==Q(con["archived_addition_gap_exact"]),
           "observer must still refer to original saved binary addition gap")
    Y=[[ref[i][j]+disp1[2*i+j] for j in range(2)] for i in range(325)]
    spread=[[rad if 2*i+j in m else Q(0) for j in range(2)] for i in range(325)]
    hull=[[[str(v-g),str(v+g)] for v,g in zip(node,err,strict=True)] for node,err in zip(Y,spread,strict=True)]
    verify(box["reference_nodes_exact"]==[[str(v) for v in n] for n in ref]
           and box["current_centers_exact"]==[[str(v) for v in n] for n in Y]
           and box["coordinate_radii_exact"]==[[str(v) for v in n] for n in spread]
           and box["coordinate_hulls_exact"]==hull and box["predictor_hulls_exact"]==[["0","0"]]*611
           and Q(box["corrector_radius_exact"])==rad and box["root_membership_assumed"] is False
           and box["any_spatial_derivative_bound_assumed"] is False,
           "independent 650-coordinate exact nodal box does not match")
    ident={"record_id":con["record_id"],"dose_exact":str(dose),"parameter_cell_sha256":jsha(cell),
           "model_observation_contract_sha256":e["model_contract_sha256"],"point_cache_sha256":con["point_cache_sha256"],
           "predictor_sha256":e["predictor_sha256"],"precision_bits":640}
    verify(all(box[k]==v for k,v in ident.items()) and all(row[k]==v for k,v in ident.items() if k!="precision_bits")
           and row["point_cache_file"]==con["point_cache"] and row["source_record_id"]==c["record_id"]
           and pre["identically_zero"] is True and pre["shape"]==[611,5]
           and pre["record_id"]==con["record_id"] and pre["parameter_cell_sha256"]==jsha(cell)
           and pre["little_endian_zero_binary64_sha256"]==digest(bytes(611*5*8)),
           "hash-linked context identity or zero predictor differs")
    rules_and_routes(A,ref,tri,m,ident,rule,routes,scale,mask,height)
    return {**{k:e[k] for k in ("case_id","record_id","model_contract_sha256","box_sha256","routes_sha256","rule_sha256",
                               "planned_geometry_radius_exact","derivative_policy")},
            "native_cached_radius_exact":str(r),"native_ZQ_upper_exact":str(z),
            "planned_radius_uses_native_matched_ZQ":r==rad,"native_cached_posterior_upper_exact":str(error),
            "center_sha256":con["center_binary_sha256"],"binary_R_sha256":con["binary_R_sha256"],
            "exact_drive":str(strength),"observer_archived_xi_exact":str(xi),"observer_shift_sigma_exact":str(sigma)}


def main():
    d={k:file(path,h,k!="source_sha256") for k,(path,h) in REQUIRED.items()}
    p,m,f=d["protocol_sha256"],d["manifest_sha256"],d["preflight_sha256"]
    pins={k:h for k,(_,h) in REQUIRED.items()}
    census=file(CENSUS,META[CENSUS])
    sample=file(OUT/"context_readback.json",META[OUT/"context_readback.json"])
    verify(sample["status"]=="PASS_PUBLIC_CONTEXT_READBACK" and sample["manifest_sha256"]==pins["manifest_sha256"]
           and len(sample["cases"])==6 and sample["wrong_external_digest_and_foreign_case_rejected"] is True,
           "factory public API test receipt missing")
    verify(f["status"]=="PASS" and f["checks"] and all(v is True for v in f["checks"].values())
           and all(f[k]==pins[k] for k in ("source_sha256","protocol_sha256","manifest_sha256"))
           and m["source_sha256"]==p["source_sha256"]==pins["source_sha256"] and m["protocol_sha256"]==pins["protocol_sha256"]
           and all(require_zero(x[k]) for x in (p,m,f) for k in NO_MATH)
           and m["geometry_or_variation_actual_launch_authorized"] is False
           and m["spatial_bridge_actual_launch_authorized"] is False,"not a cold and scope-limited frozen context")
    verify(file(PREP/"transfer_ambiguity_remaining32_exact_targets_frozen.py",pins["source_sha256"],False)==d["source_sha256"],
           "producer bytes and frozen snapshot differ")
    for k,h in (p["input_metadata_and_source_sha256"]|p["local_source_closure_sha256"]).items():
        file(k,h,False)
    for k,item in m["opened_actual_input_files"].items():
        file(k,item["sha256"],False)
    verify(len(m["opened_actual_input_files"])==100 and p["ordered_cases"]==IDS
           and [x["case_id"] for x in m["targets"]]==IDS==[x["case_id"] for x in census["cases"]],
           "incomplete or mixed sorted identities")
    out=[check_one(e,c,m,p["input_metadata_and_source_sha256"]) for e,c in zip(m["targets"],census["cases"],strict=True)]
    totals=dict(Counter(x["derivative_policy"] for x in out))
    verify(totals==m["derivative_policy_counts"]=={MATCH:22,FRESH:7,SMALL:3},"mixed derivative planning scope")
    verify({x["case_id"] for x in out if Q(x["observer_shift_sigma_exact"])>0}==RECENT,"wrong recentered subset")
    checks={
        "source_protocol_preflight_manifest_current_and_factory_negative_fixtures":True,
        "actual_cache_bytes_32_binary_R_and_current_centers_independently_read":True,
        "all_32_physical_material_exact_drive_and_normalization":True,
        "two_exact_binary_recenters_change_only_u_with_original_observer":True,
        "all_20800_exact_coordinate_hulls_and_zero_predictors_independently_rebuilt":True,
        "all_reference_masks_local_free_maps_and_auxiliary_rules_independently_rebuilt":True,
        "actual_native_radius_recomputed_no_cross_radius_promotion":True,
        "strict_cold_geometry_context_only_scope":True,
    }
    payload={"status":"PASS_INDEPENDENT_REMAINING32_TARGET_CONTEXTS",**pins,"ordered_case_ids":IDS,
             "checks":checks,"derivative_policy_counts":totals,"per_case":out,
             "review_source_sha256":digest(Path(__file__).read_bytes()),"input_sha256":LINKS,
             "independent_reconstruction":"Read binary cache archives directly; no imports or calls to producer. The six-case public context API test is a pinned complementary premise.",
             "scope":"Exact center/box, reference route and cached premise admission only. No deformed geometry or new contraction.",
             "new_geometry_certificates":0,"new_AD_derivatives":0,"new_integrations":0,"new_solves":0,
             "actual_mechanics_launch_authorized":False,"full_five_gate_complete":False}
    path=HERE/"review.json"
    with path.open("x",encoding="utf-8",newline="\n") as stream:
        json.dump(payload,stream,indent=2)
        stream.write("\n")
    print(json.dumps({"review_sha256":digest(path.read_bytes()),"cases":len(out),"input_pins":len(LINKS),"counts":totals}))


if __name__=="__main__":
    main()

