"""Exact active coefficient producer for the frozen negative K–L continuation pilot."""

import argparse
import base64
from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
import io
import json
from math import factorial, lcm
import os
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PILOT = ROOT / "runs/exact_loading_tube_pilot_preparation_v1"
PILOT_ID = "negative_KL_exact_spatial_loading_secant_v1"
K, L = (f"max_negative_contrast_change_small_{i:02}" for i in (16,17))
EX = "runs/transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1/"
CF = "runs/transfer_ambiguity_nonlinear_probe_v1/certification_design/"
K_CACHE = CF + "full_outcome_coverage_v1/remaining_cached_30_v1/recenter_v1/actual/" + K + "/recentered_cache_packet.json"
L_CACHE = CF + "cached_pilots_v1/actual/" + L + "/cache_packet.json"
PILOT_PINS = {
    "source_sha256": "d39ae6c3bfa3d932eaba03ddd59fc275fd5b3b1c094e7f9a544810d8c10dc8bc",
    "protocol_sha256": "fe64074d1a169a7672039c73dc07ca35c728cc1c299cf2a28025e0a0e7a97e43",
    "preflight_sha256": "911f6839bd9cc49dafe626e688cc30d88a283ad87b5aaf92fecc5b629a684881",
    "plan_sha256": "c3b1870e5ad815bc5df5ec94d4731b0228e529ba2e0651f847102b8e0ecaa277",
    "envelope_sha256": "f50c04b1802176abf4f5494898c8128018d4db8f8d0bb410af7774e5a3fff925",
}
PINNED = {
    K_CACHE: "c7ef4cb944fb73ec0b498f97dbebf3727d4a7e5773957a1b573a57a8b77b6549",
    L_CACHE: "7c94b3f989cab5fb31d979851b3ee94a155b96e2251a34092f1cab2b6a0c96a3",
    "runs/exact_loading_tube_design_v1/receipt.json": "dbfb9871da2b384fcd57e44015a35bfeebbfccf25e23a6f45178f3aef488ff3c",
    EX + "loading_affine_identity_v1/independent_readback/review_receipt.json":
        "2e77f4dc2a6cd369f6368b3028a408fdd970cf0c9e3f892f08f9e60727a1be9a",
    EX + "negative_KL_active_energy_v1/independent_actual/root_acceptance.json":
        "e2c36606b3911ac51ed99b56e0fb6477235cc5666ca78d152b85473b7dc6c3e1",
}
EDGES = ((0,1),(1,2),(2,0))
N = 611
DIMENSIONS = {"reference":(325,2),"free":(N,),"triangles":(144,6),"active":(144,1),
              "energy_scale":(1,),"u":(N,),"preconditioner":(N,N)}
INT_ARRAYS = {"free","triangles"}
REQUIRED_CHECKS = ("same_frozen_K_R_cache_reference_material_gain_and_normalization",
                   "same_KL_exact_secant_and_signed_projection",
                   "positive_orientation_and_exact_P2_axial_cubic_moments",
                   "common_directed_factor_and_exact_signed_sparse_R_product",
                   "same_old_K_cube_short_load_scope_and_strict_tests",
                   "run_is_write_once_and_creates_no_subprocesses")


def require(condition,message):
    if not condition:
        raise ValueError(message)


def sha(payload):
    return sha256(payload).hexdigest()


def digest(path):
    return sha(path.read_bytes())


def canonical(name):
    path = (ROOT/name).resolve()
    require(path.is_relative_to(ROOT) and path.relative_to(ROOT).as_posix()==name,
            "canonical study-relative input required")
    return path


def read(path,expected):
    payload = path.read_bytes()
    require(sha(payload)==expected,f"pinned input changed: {path.name}")
    return json.loads(payload)


def write(path,value):
    with path.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(value,f,indent=2,default=str,allow_nan=False)
        f.write("\n")
    return digest(path)


def exact(x):
    require(type(x) in (str,int,Q),"exact rational required")
    return Q(x)


def dot(a,b):
    return sum((x*y for x,y in zip(a,b,strict=True)),Q(0))


def triangle_geometry(ref):
    (r0,z0),(r1,z1),(r2,z2) = ref
    ds_r,dt_r = r1-r0,r2-r0
    determinant = ds_r*(z2-z0)-dt_r*(z1-z0)
    radii = (r0,r1,r2)
    require(determinant>0 and min(radii)>=0 and max(radii)>0,
            "positive material orientation and nonnegative nontrivial reference radius required")
    return determinant,radii,((dt_r-ds_r)/determinant,-dt_r/determinant,ds_r/determinant)


def cubic_moment(i,j,k):
    """Integral of lambda_i lambda_j lambda_k on the unit simplex."""
    return Q(factorial((i,j,k).count(0))*factorial((i,j,k).count(1))
             *factorial((i,j,k).count(2)),120)


def active_gram(ref):
    """Exact six-node axial Gram matrix, before the common 2*pi*g/S factor."""
    determinant,radii,dz = triangle_geometry(ref)
    grad = [tuple(dz[i]*(3 if j==i else -1) for j in range(3)) for i in range(3)]
    grad.extend(tuple(4*((dz[i] if m==j else 0)+(dz[j] if m==i else 0)) for m in range(3))
                for i,j in EDGES)
    radial_moments = [[sum(r*cubic_moment(i,j,k) for k,r in enumerate(radii))
                       for j in range(3)] for i in range(3)]
    return tuple(tuple(determinant*sum(a[i]*b[j]*radial_moments[i][j]
                                        for i in range(3) for j in range(3)) for b in grad) for a in grad)


def checked_map(nodes,local,indices,projection):
    expected = [(2*i+c,projection[2*n+c]) for i,n in enumerate(nodes) for c in (0,1)
                if 2*n+c in projection]
    require(type(local) is list and type(indices) is list and
            all(type(j) is int for j in local+indices) and list(zip(local,indices,strict=True))==expected,
            "ordered exact free-coordinate projection mismatch")


def scatter_active(ref,triangles,free,slots):
    """C0 and d0 where C=alpha*C0 and d=alpha*d0, with full reference offsets."""
    projection = {i:j for j,i in enumerate(free)}
    n = len(free)
    rows = [defaultdict(Q) for unused in range(n)]
    offset = [Q(0)]*n
    leaves = []
    for slot in slots:
        nodes = triangles[slot]
        points = [ref[i] for i in nodes]
        local = active_gram(points[:3])
        for comp in (0,1):
            gather = [projection.get(2*j+comp) for j in nodes]
            ref_comp = [point[comp] for point in points]
            for index, r in zip(gather,local,strict=True):
                if index is None:
                    continue
                offset[index] += dot(r,ref_comp)
                for j,a in zip(gather,r,strict=True):
                    if j is not None and a:
                        rows[index][j] += a
        leaves.append({"element":slot,"nodes":nodes,"axial_gram_exact":local})
    out = [{j:x for j,x in sorted(row.items()) if x} for row in rows]
    require(all(row.get(i,0)>=0 and all(x==out[j].get(i,0) for j,x in row.items()) for i,row in enumerate(out)),
            "the exact projected active Gram matrix must be symmetric with nonnegative diagonal")
    return out,offset,leaves


def multiply_signed(rows,offset,R,u,v):
    """Exact signed R*A0 and R*d0 actions, summing one sparse product row at a time."""
    n = len(rows)
    require(n>0 and len(offset)==len(R)==len(u)==len(v)==n,"all projected inputs must have the same dimension")
    common = lcm(*(a.denominator for row in rows for a in row.values()),
                 *(a.denominator for a in offset))
    integers = [{j:a.numerator*(common//a.denominator) for j,a in row.items()} for row in rows]
    b_int = [a.numerator*(common//a.denominator) for a in offset]
    source_rows = [i for i,row in enumerate(integers) if row or b_int[i]]
    a_out,k_out,active_norms,R_norms = [],[],[],[]
    for row in R:
        require(len(row)==n,"square binary R required")
        ratios = [r.as_integer_ratio() for r in row]
        require(all(d>0 and (d & (d-1))==0 for _,d in ratios),"exact binary coefficients required in R")
        rden = max(d for _,d in ratios)
        r_int = [m*(rden//d) for m,d in ratios]
        den = rden*common
        signed = defaultdict(int)
        signed_offset = 0
        for i in source_rows:
            x = r_int[i]
            signed_offset += x*b_int[i]
            for j,m in integers[i].items():
                signed[j] += x*m
        a_out.append((sum((m*u[j] for j,m in signed.items()),Q(signed_offset)))/den)
        k_out.append(sum((m*v[j] for j,m in signed.items()),Q(0))/den)
        active_norms.append(Q(sum(abs(x) for x in signed.values()),den))
        R_norms.append(Q(sum(abs(x) for x in r_int),rden))
    return a_out,k_out,active_norms,R_norms


def directed_factor(gain,scale):
    from flint import arb, ctx, fmpq

    require(gain>0 and scale>0,"positive exact gain and binary normalization required")
    old = ctx.prec
    try:
        ctx.prec = 640
        multiplier = 2*gain/scale
        ball = arb.pi()*arb(fmpq(multiplier.numerator,multiplier.denominator))
        lower,upper = (Q(str(x.fmpq())) for x in (ball.lower(),ball.upper()))
        require(0<lower<=upper,"common factor must be positive")
        return lower,upper
    finally:
        ctx.prec = old


def scale_signed(value,factor):
    products = sorted(value*q for q in factor)
    return tuple(products)


def old_K_extension(eta,z,r,Rnorm,gap,B,M):
    """Sufficient certificate arithmetic for an unchanged-geometry old-cube subinterval."""
    require(0<=eta and 0<=z<1 and r>0 and Rnorm>0 and gap>0 and B>=0 and M>=0,
            "valid nonnegative same-old-cube estimates required")
    margin = (1-z)*r-eta
    require(margin>0,"strict original exact K cube bound required")
    drift = B+r*M
    h = min(gap/2,margin/(2*drift)) if drift else gap/2
    new_z,new_eta = z+h*M,eta+h*B
    new_margin = (1-new_z)*r-new_eta
    require(0<h<gap and new_z<1 and new_margin>0,"strict positive extension margin required")
    return {"status":"CANDIDATE_SAME_EXACT_K_CUBE_SHORT_LOAD_EXTENSION_AWAITS_INDEPENDENT_REVIEW",
            "positive_dose_increment_exact":h,"original_K_radius_exact":r,
            "eta_upper_exact":new_eta,"Z_upper_exact":new_z,"strict_selfmap_margin_lower_exact":new_margin,
            "posterior_upper_exact":new_eta/(1-new_z),"normalized_H_curvature_lower_exact":(1-new_z)/Rnorm,
            "predictor":"constant exact K proof center","links_only_to_original_K_root":True,
            "reaches_L":False,"whole_KL_path_certified":False}


def cache_arrays(packet):
    """Decode only arrays needed by this operator, preserving little-endian byte identity."""
    import numpy as np

    block = packet["frozen_binary_cache"]
    blob = base64.b64decode(block["data"],validate=True)
    require(type(block["bytes"]) is int and len(blob)==block["bytes"] and sha(blob)==block["sha256"],
            "exact embedded frozen cache required")
    with np.load(io.BytesIO(blob),allow_pickle=False) as archive:
        require(len(archive.files)==len(set(archive.files)) and set(DIMENSIONS).issubset(archive.files),
                "complete cache array catalogue required")
        arrays = {name:archive[name] for name in DIMENSIONS}
    for name, dimensions in DIMENSIONS.items():
        a = arrays[name]
        require(a.shape==dimensions and a.dtype.str==("<i8" if name in INT_ARRAYS else "<f8")
                and np.isfinite(a).all(),f"literal shape, dtype and finite entries required for {name}")
    return arrays


def binary(x):
    return Q.from_float(float(x))


def bytes_sha(a):
    return sha(a.tobytes(order="C"))


def decode_anchor(plan,box,letter,packet,a):
    case = (K,L)[letter]
    dose = exact(plan["dose_interval_exact"][letter])
    gain,scale = map(exact,(plan["gain_exact"],plan["energy_scale_binary_exact"]))
    check_status = ("FROZEN_SAME_CACHED_OPERATOR_AT_NEW_PROOF_CENTER",
                    "FROZEN_FRESH_EXACT_BINARY_QUADRATURE_CACHE")[letter]
    require(packet["status"]==check_status and packet.get("case_id",packet.get("pilot_id"))==case
            and type(packet["free_dof"]) is int and packet["free_dof"]==N
            and packet["material_parameter_storage_exact"]==plan["material_parameters_binary_exact"]
            and exact(packet["optical_dose_exact"])==exact(box["dose_exact"])==dose
            and exact(packet["common_gain_exact"])==gain
            and exact(packet["activation_strength_exact"])==dose*gain
            and binary(a["energy_scale"][0])==exact(packet["normalization_binary_exact"])==scale
            and box["record_id"]=="transfer_exact_"+case,
            "frozen material, normalization, optical dose and requested drive must agree")
    free = a["free"].tolist()
    ref = [tuple(map(binary,row)) for row in a["reference"]]
    u = list(map(binary,a["u"]))
    require(free==box["free"] and a["triangles"].tolist()==box["triangles"]
            and ref==[tuple(map(exact,row)) for row in box["reference_nodes_exact"]]
            and u==list(map(exact,plan[("K","L")[letter]+"_center_displacement_exact"])),
            "same ordered projection, binary reference and exact displacement centers required")
    projection = dict(zip(free,u,strict=True))
    reconstructed = [tuple(ref[n][c]+projection.get(2*n+c,Q(0)) for c in (0,1)) for n in range(325)]
    require(reconstructed==[tuple(map(exact,row)) for row in box["current_centers_exact"]],
            "proof center must be literal X+Pu arithmetic")
    return ref,u


def load_inputs():
    for filename,key in (("pilot.py","source"),("protocol.json","protocol"),("preflight.json","preflight"),
                         ("plan.json","plan"),("envelope.json","envelope")):
        require(digest(PILOT/filename)==PILOT_PINS[key+"_sha256"],"closed pilot source and phase required")
    plan = read(PILOT/"plan.json",PILOT_PINS["plan_sha256"])
    env = read(PILOT/"envelope.json",PILOT_PINS["envelope_sha256"])
    inputs = dict(PINNED)

    def linked(name,h):
        require(inputs.setdefault(name,h)==h,"admission path has inconsistent digests")
        return read(canonical(name),h)

    crosswalk = linked("runs/exact_loading_tube_design_v1/receipt.json",PINNED["runs/exact_loading_tube_design_v1/receipt.json"])
    law = EX+"loading_affine_identity_v1/independent_readback/review_receipt.json"
    energy_order = EX+"negative_KL_active_energy_v1/independent_actual/root_acceptance.json"
    require(linked(law,PINNED[law])["status"]=="PASS_INDEPENDENT_AFFINE_LOADING_SOURCE_AND_MATHEMATICS"
            and linked(energy_order,PINNED[energy_order])["certified_loading_path"] is False
            and plan["pilot_id"]==PILOT_ID and plan["case_ids"]==[K,L]
            and plan["path_certified"] is False and plan["exact_active_direction"]=="e_Z"
            and crosswalk["pilot"]["active_slots"]==plan["active_slots"]
            and crosswalk["pilot"]["same_gain_exact"]==plan["gain_exact"],"same affine active K–L target required")
    kb,lb=(linked(p,h) for p,h in zip(plan["anchor_box_files"],plan["anchor_box_sha256"],strict=True))
    kn,ln=(linked(p,h) for p,h in zip(plan["anchor_spatial_proof_files"],plan["anchor_spatial_proof_sha256"],strict=True))
    routes = linked(plan["K_routes_file"],plan["K_routes_sha256"])
    pk,pl=(linked(name,PINNED[name]) for name in (K_CACHE,L_CACHE))
    ka,la=map(cache_arrays,(pk,pl))
    ref,uk=decode_anchor(plan,kb,0,pk,ka)
    rl,ul=decode_anchor(plan,lb,1,pl,la)
    require(ref==rl and all(bytes_sha(ka[k])==bytes_sha(la[k]) for k in DIMENSIONS if k not in {"u","preconditioner"})
            and kn["target"]==plan["K_binding"]
            and kn["target"]["point_cache_sha256"]==PINNED[K_CACHE]
            and ln["target_identity"]["cache_packet_sha256"]==PINNED[L_CACHE]
            and bytes_sha(ka["preconditioner"])==plan["K_binding"]["R_binary_sha256"]
            and bytes_sha(ka["u"])==plan["K_binding"]["center_binary_sha256"]
            and bytes_sha(la["preconditioner"])==ln["target_identity"]["R_binary_sha256"],
            "same physical arrays, K proof center and both authentic binary preconditioners required")
    low,high=map(exact,plan["dose_interval_exact"])
    h=high-low
    v=list(map(exact,plan["predictor_slope_exact"]))
    require(exact(plan["nominal_dose_exact"])==low and exact(plan["dose_gap_exact"])==h>0
            and exact(plan["nominal_drive_exact"])==exact(plan["gain_exact"])*low
            and all(a+h*b==c for a,b,c in zip(uk,v,ul,strict=True)),"actual frozen K–L secant required")
    free=ka["free"].tolist()
    triangles=kb["triangles"]
    projection={i:j for j,i in enumerate(free)}
    outer=exact(plan["material_parameters_binary_exact"]["outer_radius"])
    expected_free=[2*n+c for n,(r,_) in enumerate(ref) for c in (0,1) if r!=outer and (r!=0 or c==1)]
    require(free==expected_free==env["free"] and triangles==env["triangles"]
            and kb["reference_nodes_exact"]==env["reference_nodes_exact"],"same exact fixed axis/rim and nodal injection required")
    slots=plan["active_slots"]
    mask=ka["active"][:,0].tolist()
    require(type(slots) is list and len(slots)==24 and all(type(i) is int and 0<=i<144 for i in slots)
            and all(w in (0.,1.) for w in mask) and slots==[i for i,w in enumerate(mask) if w==1.]
            and [x["element"] for x in routes["bulk"]]==list(range(144)),"complete exact reference mask and routes required")
    counts=Counter()
    for i,row in enumerate(routes["bulk"]):
        nodes=triangles[i]
        require(type(row["element"]) is int and len(nodes)==len(set(nodes))==6
                and all(type(n) is int and 0<=n<len(ref) for n in nodes)
                and nodes==row["nodes"] and exact(row["active"])==int(mask[i]),"canonical P2 node order and active route required")
        checked_map(nodes,row["free_local"],row["global_indices"],projection)
        _,radii,_=triangle_geometry([ref[n] for n in nodes[:3]])
        kind=("positive_radius","axis_vertex","axis_edge")[radii.count(0)]
        require(row["chart"]=={"positive_radius":"off_axis","axis_vertex":"axis_vertex",
                               "axis_edge":"axis_edge_physical"}[kind],
                f"changed chart at slot {i}: {row['chart']}")
        for mid,(a,b) in zip(nodes[3:],EDGES,strict=True):
            represented=(ka["reference"][nodes[a]]+ka["reference"][nodes[b]])/2
            require(represented.tolist()==ka["reference"][mid].tolist(),"changed frozen P2 edge ordering")
        if i in slots:
            counts[kind]+=1
    require(exact(kn["radius_exact"])==exact(kb["corrector_radius_exact"])==Q(1,2**34)
            and kn["status"]=="PASS_INDEPENDENT_FIRST_WAVE_EXACT_SPATIAL_ROOT"
            and kn["arithmetic"]["strict_selfmap"] is True
            and kn["observer_and_optional_minimum"]["strict_minimum_proved"] is True,
            "original K exact cube proof required for optional short extension")
    return plan,kn,inputs,ref,triangles,free,slots,counts,ka["preconditioner"],uk,v


def phase_bindings():
    return {"source_sha256":digest(Path(__file__)),"protocol_sha256":digest(HERE/"protocol.json"),
            "preflight_sha256":digest(HERE/"preflight.json")}


def stage_contract(protocol,preflight,phase,inputs):
    require(protocol["status"]=="COLD_KL_EXACT_ACTIVE_OPERATOR_PROTOCOL"
            and protocol["source_sha256"]==phase["source_sha256"]==preflight["source_sha256"]
            and preflight["protocol_sha256"]==phase["protocol_sha256"]
            and protocol["pilot_bindings"]==preflight["pilot_bindings"]==PILOT_PINS
            and protocol["input_sha256"]==inputs
            and preflight["status"]=="PASS_KL_EXACT_ACTIVE_OPERATOR_COLD_FIXTURES"
            and all(type(preflight[k]) is int and preflight[k]==0 for k in
                    ("actual_active_elements_integrated","actual_R_A_products","new_AD_derivatives","new_solves"))
            and protocol["producer_authorized"] is False,"same complete cold source, inputs and phase required")


def launch_gate(review,phase):
    require(review["status"]=="PASS_INDEPENDENT_KL_ACTIVE_OPERATOR_PRELAUNCH"
            and review["pilot_id"]==PILOT_ID and review["producer_bindings"]==phase
            and review["pilot_bindings"]==PILOT_PINS
            and isinstance(review["checks"],dict) and all(review["checks"].get(k) is True for k in REQUIRED_CHECKS)
            and review["active_operator_only_authorized"] is True
            and review["passive_or_geometry_jobs_authorized"] is False,
            "externally pinned active-operator-only cold-review release required")


def material_packet(plan,kn,ref,triangles,free,slots,counts,R,uk,v):
    rows,offset,leaves=scatter_active(ref,triangles,free,slots)
    a0,k0,row0,R_rows=multiply_signed(rows,offset,R,uk,v)
    gain,scale=map(exact,(plan["gain_exact"],plan["energy_scale_binary_exact"]))
    factor=directed_factor(gain,scale)
    a=[scale_signed(x,factor) for x in a0]
    k=[scale_signed(x,factor) for x in k0]
    M= factor[1]*max(row0)
    B=max(abs(y) for x in a for y in x)
    Rnorm=max(R_rows)
    require(Rnorm==exact(plan["same_R_infinity_norm_exact"]),"same exact infinity norm of authenticated K R required")
    ar=kn["arithmetic"]
    gap=exact(plan["dose_gap_exact"])
    r=exact(kn["radius_exact"])
    extension=old_K_extension(exact(ar["eta_upper_exact"]),exact(ar["z_upper_exact"]),r,Rnorm,gap,B,M)
    extension["exact_load_interval"]=[exact(plan["nominal_dose_exact"]),
                                    exact(plan["nominal_dose_exact"])+extension["positive_dose_increment_exact"]]
    extension["exact_original_K_spatial_proof_sha256"]=plan["anchor_spatial_proof_sha256"][0]
    sparse={"local_axial_grams":leaves,"global_exact_A0_upper":
            [[i,j,a] for i,row in enumerate(rows) for j,a in row.items() if i<=j],
            "global_exact_reference_b0":offset,"dimension":len(rows),"component_cross_coupling":False,
            "active_gram_and_offset_exclude_common_factor":True}
    result={"status":"CANDIDATE_KL_EXACT_ACTIVE_OPERATOR_BOUNDS",
            "R_binary_sha256":plan["K_binding"]["R_binary_sha256"],
            "point_cache_sha256":plan["K_binding"]["point_cache_sha256"],
            "center_binary_sha256":plan["K_binding"]["center_binary_sha256"],
            "nominal_drive_exact":plan["nominal_drive_exact"],
            "energy_scale_binary_exact":scale,"fiber_radius_exact":plan["fiber_radius_exact"],
            "common_gain_exact":gain,"optical_dose_interval_exact":plan["dose_interval_exact"],
            "exact_R_infinity_norm":Rnorm,"common_2_pi_gain_over_energy_scale_interval_exact":factor,
            "signed_R_active_center_rows_exact":a,"signed_R_C_v_rows_exact":k,
            "R_C_norm_upper_exact":M,"R_active_center_norm_upper_exact":B,
            "pre_factor_R_A0_absolute_row_sums_exact":row0,
            "pre_factor_signed_R_active_center_exact":a0,"pre_factor_signed_R_A0_v_exact":k0,
            "active_element_slots":slots,"active_element_chart_counts":counts,"active_terms_only":True,
            "common_factor_precision_bits":640,"new_AD_derivatives":0,"new_passive_integrations":0,
            "new_solves":0,"whole_KL_path_certified":False,
            "separate_same_cube_short_load_candidate":extension}
    return result,sparse


def run(review_sha):
    phase=phase_bindings()
    protocol=read(HERE/"protocol.json",phase["protocol_sha256"])
    preflight=read(HERE/"preflight.json",phase["preflight_sha256"])
    review=read(HERE/"independent_prelaunch.json",review_sha)
    launch_gate(review,phase)
    plan,kn,inputs,*data=load_inputs()
    stage_contract(protocol,preflight,phase,inputs)
    out=HERE/"actual"
    out.mkdir(exist_ok=False)
    ident={"pilot_id":PILOT_ID,"pid":os.getpid(),"producer_bindings":phase,"bindings":PILOT_PINS,
           "independent_active_prelaunch_sha256":review_sha}
    write(out/"start.json",ident)
    try:
        result,sparse=material_packet(plan,kn,*data)
        gram_sha=write(out/"exact_sparse_operator.json",{**ident,**sparse})
        res=write(out/"result.json",{**ident,**result,"exact_sparse_operator_sha256":gram_sha})
    except Exception as error:
        failed=write(out/"failure.json",{**ident,"type":type(error).__name__,"error":str(error)})
        write(out/"terminal.json",{**ident,"status":"KL_LOADING_ACTIVE_OPERATOR_FAILED","failure_sha256":failed})
        raise
    write(out/"terminal.json",{**ident,"status":"KL_LOADING_OPERAND_STAGE_FINISHED",
                                "mathematical_status":result["status"],"result_sha256":res})


def rejected(fn):
    try:
        fn()
    except (ValueError,TypeError):
        return True
    return False


def fixtures():
    ref=[(Q(0),Q(0)),(Q(2),Q(0)),(Q(0),Q(3))]
    points=ref+[
        tuple((a+b)/2 for a,b in zip(ref[i],ref[j],strict=True)) for i,j in EDGES]
    gram=active_gram(ref)
    profiles=([z*z for _,z in points],[r*z for r,z in points])
    oracle=sum((dot(y,[dot(row,y) for row in gram])/2 for y in profiles),Q(0))
    one=[Q(1)]*6
    axial=[z for _,z in points]
    shift_matrix=active_gram([(r,z+Q(23,7)) for r,z in ref])
    # At positive radius the integrand is 2*(1+s)*t^2 and its integral is 1/5.
    off_ref=[(Q(1),Q(0)),(Q(2),Q(0)),(Q(1),Q(1))]
    off_pts=off_ref+[tuple((a+b)/2 for a,b in zip(off_ref[i],off_ref[j],strict=True)) for i,j in EDGES]
    off_gram=active_gram(off_ref)
    off_z=[z*z for _,z in off_pts]
    ax1=[(Q(0),Q(0)),(Q(1),Q(0)),(Q(1),Q(1))]
    ax1_gram=active_gram(ax1)
    ax1_u=[Q(0),Q(0),Q(1),Q(0),Q(1,4),Q(1,4)]
    checks={"exact_axis_polynomial_energy_and_translation_invariance":oracle==Q(24,5) and shift_matrix==gram,
            "cubic_rule_preserves_axis_partition_and_axial_linear_energy":
            all(dot(row,one)==0 for row in gram) and dot(axial,[dot(r,axial) for r in gram])==2,
            "independent_positive_radius_polynomial_oracle":dot(off_z,[dot(row,off_z) for row in off_gram])/2==Q(1,5),
            "independent_axis_vertex_polynomial_oracle":dot(ax1_u,[dot(row,ax1_u) for row in ax1_gram])/2==Q(2,15),
            "negative_orientation_and_bad_radius_fail_before_moments":all((
                rejected(lambda:active_gram([ref[0],ref[2],ref[1]])),
                rejected(lambda:active_gram([ref[0],ref[1],(Q(-1),Q(3))]))))}
    free=[1,2,3,6,7,8,9,11]   # Includes fixed nonzero axial reference coordinate at full index 5.
    rows,offset,_=scatter_active(points,[list(range(6))],free,[0])
    reference_free=[points[i//2][i%2] for i in free]
    wrong=[sum((x*reference_free[j] for j,x in row.items()),Q(0)) for row in rows]
    checks["fixed_reference_offset_retained_and_component_cross_terms_zero"] = (
        offset!=wrong and all((free[i]-free[j])%2==0 for i,row in enumerate(rows) for j in row))
    proj={i:j for j,i in enumerate(free)}
    checked_map(list(range(6)),free,list(range(len(free))),proj)
    checks["global_projection_rejects_bool_float_and_wrong_permutation"] = all((
        rejected(lambda:checked_map(list(range(6)),[True,*free[1:]],list(range(len(free))),proj)),
        rejected(lambda:checked_map(list(range(6)),free,[0.0,*range(1,len(free))],proj)),
        rejected(lambda:checked_map(list(range(6)),free,list(reversed(range(len(free)))),proj))))
    a,k,ar,rr=multiply_signed([{0:Q(2),1:Q(-2)},{0:Q(-2),1:Q(2)}],[Q(1),Q(-1)],
                             [[1.,1.],[-.5,1.]],[Q(2),Q(-1)],[Q(1),Q(0)])
    checks["signed_sparse_products_cancel_before_taking_norm"]=(
        a==[0,Q(-21,2)] and k==[0,-3] and ar==[0,6] and rr==[2,Q(3,2)]
        and rejected(lambda:multiply_signed([{0:Q(1)}],[Q(1)],[[Q(1,3)]],[Q(0)],[Q(0)])))
    pi=directed_factor(Q(1),Q(2))
    checks["shared_positive_directed_pi_factor_respects_negative_rows"]=(
        Q(333,106)<pi[0]<=pi[1]<Q(355,113) and scale_signed(-2,pi)==(-2*pi[1],-2*pi[0]))
    short=old_K_extension(Q(1,64),Q(1,4),Q(1,8),Q(2),Q(1,2),Q(1,8),Q(1,4))
    stationary=old_K_extension(Q(0),Q(0),Q(1),Q(1),Q(1,2),Q(0),Q(0))
    checks["short_extension_preserves_strict_margin_and_cannot_reach_L"]=(
        short["positive_dose_increment_exact"]==Q(1,4) and short["strict_selfmap_margin_lower_exact"]==Q(5,128)
        and short["reaches_L"] is False and stationary["positive_dose_increment_exact"]==Q(1,4)
        and rejected(lambda:old_K_extension(Q(1,2),Q(1,2),Q(1),Q(1),Q(1),Q(0),Q(0))))
    phase={k+"_sha256":"0"*64 for k in ("source","protocol","preflight")}
    review={"status":"PASS_INDEPENDENT_KL_ACTIVE_OPERATOR_PRELAUNCH","pilot_id":PILOT_ID,
            "producer_bindings":phase,"pilot_bindings":PILOT_PINS,"checks":dict.fromkeys(REQUIRED_CHECKS,True),
            "active_operator_only_authorized":True,"passive_or_geometry_jobs_authorized":False}
    launch_gate(review,phase)
    checks["external_review_allows_active_only_and_rejects_wrong_phase"] = rejected(
        lambda:launch_gate(review,{**phase,"preflight_sha256":"f"*64})) and rejected(
        lambda:launch_gate({**review,"passive_or_geometry_jobs_authorized":True},phase))
    require(all(checks.values()),f"failed cold fixtures: {checks}")
    return checks


def prepare():
    checks=fixtures()
    plan,kn,inputs,ref,tri,free,slots,counts,_,u,v=load_inputs()
    checks["all_actual_identities_maps_and_secant_matched_without_active_assembly"]=True
    source=digest(Path(__file__))
    protocol=write(HERE/"protocol.json",{"status":"COLD_KL_EXACT_ACTIVE_OPERATOR_PROTOCOL",
        "pilot_id":PILOT_ID,"pilot_bindings":PILOT_PINS,"source_sha256":source,"input_sha256":inputs,
        "hardcoded_law":"axial P2 reference-weighted Gram and full reference offset, exact state independent",
        "review_path":"independent_prelaunch.json","required_review_checks":REQUIRED_CHECKS,
        "producer_authorized":False,"passive_or_geometry_jobs_authorized":False,
        "short_extension_separate_and_conditional":True,"full_KL_path_certified":False})
    write(HERE/"preflight.json",{"status":"PASS_KL_EXACT_ACTIVE_OPERATOR_COLD_FIXTURES",
        "pilot_bindings":PILOT_PINS,"source_sha256":source,"protocol_sha256":protocol,"checks":checks,
        "active_slots":slots,"active_chart_counts":counts,"coordinates":2*len(ref),"triangles":len(tri),
        "free":len(free),"actual_active_elements_integrated":0,"actual_R_A_products":0,
        "new_AD_derivatives":0,"new_passive_integrations":0,"new_solves":0})


def verify_preparation():
    phase=phase_bindings()
    protocol=read(HERE/"protocol.json",phase["protocol_sha256"])
    preflight=read(HERE/"preflight.json",phase["preflight_sha256"])
    _,_,inputs,*_=load_inputs()
    stage_contract(protocol,preflight,phase,inputs)
    print(json.dumps({"producer_bindings":phase,"checks":fixtures()},indent=2))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_subparsers(dest="command",required=True)
    modes.add_parser("prepare")
    modes.add_parser("verify-preparation")
    modes.add_parser("run").add_argument("--review-sha",required=True)
    args=parser.parse_args()
    {"prepare":prepare,"verify-preparation":verify_preparation,"run":lambda:run(args.review_sha)}[args.command]()


if __name__=="__main__":
    main()
