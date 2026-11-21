"""Cold/actual scalar active-energy test at the certified negative K and L roots."""
import argparse
import base64
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
import itertools
import json
from math import factorial
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EX = HERE.parent
BR = EX / "remaining32_bridge_v1"
ROOT = EX.parents[3]
K = "max_negative_contrast_change_small_16"
L = "max_negative_contrast_change_small_17"
RADIUS = Q(1, 2**34)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def load(path, sha, pins):
    require(type(sha) is str and len(sha) == 64 and digest(path) == sha, f"incorrect pin {path}")
    pins[path.relative_to(ROOT).as_posix()] = sha
    return json.loads(path.read_bytes())


def save(path, doc):
    with path.open("x", encoding="utf8", newline="\n") as stream:
        json.dump(doc, stream, indent=2)
        stream.write("\n")
    return digest(path)


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def scale(c, a):
    return (c*a[0], c*a[1]) if c >= 0 else (c*a[1], c*a[0])


def product(a, b):
    corners = [x*y for x in a for y in b]
    return min(corners), max(corners)


def square(a):
    lo, hi = a
    return (Q(0) if lo <= 0 <= hi else min(lo*lo, hi*hi)), max(lo*lo, hi*hi)


def simplex_moment(det, i, j, k):
    powers = [0, 0, 0]
    for pos in (i, j, k):
        powers[pos] += 1
    return abs(det) * Q(factorial(powers[0])*factorial(powers[1])*factorial(powers[2]), 120)


def element(vertices, Y):
    """Enclose (1/2) int R_ref |d_Z Y|² on an affine P1 material triangle."""
    R, Z = zip(*vertices, strict=True)
    det = (R[1]-R[0])*(Z[2]-Z[0]) - (R[2]-R[0])*(Z[1]-Z[0])
    require(det != 0 and min(R) >= 0 and max(R) > 0, "nondegenerate nonnegative-radius material triangle")
    ld = [(R[(i+2)%3]-R[(i+1)%3])/det for i in range(3)]
    deriv = []
    for v in range(3):
        grad = [(Q(3) if i == v else Q(-1))*ld[i] for i in range(3)]
        grad += [4*(ld[j] if v == i else ld[i] if v == j else Q(0)) for i,j in ((0,1),(1,2),(2,0))]
        deriv.append([sum((scale(c, y[axis])[k] for c, y in zip(grad,Y,strict=True)),Q(0))
                      for axis in (0,1) for k in (0,1)])
    energy = Q(0), Q(0)
    for i,j in itertools.combinations_with_replacement(range(3),2):
        moment = sum(R[k]*simplex_moment(det,i,j,k) for k in range(3))
        for axis in (0,1):
            q = tuple(deriv[i][2*axis:2*axis+2])
            h = tuple(deriv[j][2*axis:2*axis+2])
            energy = add(energy, scale(moment*(Q(1,2) if i == j else 1),
                                       square(q) if i == j else product(q,h)))
    return max(Q(0), energy[0]), energy[1]


def classification(k, l):
    require(0 <= k[0] <= k[1] and 0 <= l[0] <= l[1], "ordered nonnegative energy enclosures")
    delta = l[0]-k[1], l[1]-k[0]
    status = ("EXCLUDES_SMOOTH_STABLE_SAME_OPERATOR_CONNECTION" if delta[0] > 0
              else "NECESSARY_ENERGY_ORDER_VERIFIED_ONLY" if delta[1] <= 0
              else "UNRESOLVED_NECESSARY_ENERGY_ORDER")
    return status, delta


def fixtures():
    vertices = [(Q(0),Q(0)),(Q(1),Q(0)),(Q(0),Q(1))]
    points = vertices+[(Q(1,2),Q(0)),(Q(1,2),Q(1,2)),(Q(0),Q(1,2))]
    def interval(x):
        return x,x
    affine = [[interval(R+2*Z),interval(3*Z)] for R,Z in points]
    # |d_Z Y|²=13 and integral R=1/6.
    assert element(vertices,affine)==(Q(13,12),)*2
    quadratic = [[interval(Z*Z),interval(R*Z)] for R,Z in points]
    # (1/2) int R[4 Z²+R²] = 2/60+1/40.
    assert element(vertices,quadratic)==(Q(7,120),)*2
    scale_interval = (Q(-1,2), Q(3,2))
    family = [[interval(Q(0)),scale(Z,scale_interval)] for R,Z in points]
    # Includes the exact endpoint examples, -1/2, 0, +3/2.
    enclosure = element(vertices,family)
    assert enclosure[0]==0 and enclosure[1]>=Q(3,16)
    order=[0,2,1,5,4,3]
    assert element([vertices[i] for i in order[:3]],[quadratic[i] for i in order])==(Q(7,120),)*2
    assert simplex_moment(Q(1),0,0,0)==Q(1,20)
    assert simplex_moment(Q(1),0,1,2)==Q(1,120)
    assert classification((Q(2),Q(3)),(Q(4),Q(5)))[0].startswith("EXCLUDES")
    assert classification((Q(2),Q(3)),(Q(1),Q(2)))[0].startswith("NECESSARY")
    assert classification((Q(2),Q(4)),(Q(1),Q(3)))[0].startswith("UNRESOLVED")
    return dict(exact_affine_unit_triangle=True,exact_quadratic_weighted_unit_triangle=True,
                interval_with_sign_crossing_contains_all_three_witnesses=True,
                reversed_material_vertex_order_same_integral=True,
                strict_exclusion_closed_pass_overlap_unresolved=True)




AUTHORITY = {
    EX/"remaining32_geometry_v1/preparation/protocol.json": "754da64bb96cd8a19707181315e6c3e3f014e8d37a898b83c1693342c82c9850",
    EX/"preparation/target_manifest.json": "cb25fbeaef753cea700c187cd8b377b7dcf4aa1fe2b0b246ed17134a07a89bc3",
    EX/"loading_affine_identity_v1/READBACK.md": "61e1a3a9362c8ce066588719c002cd6b66cba75840c8412b8fd5f9b3188f2b28",
    ROOT/"scripts/exact_fe_target_safe.py": "889b61cfa794d432021fbd4670b781169033d4b8451554b18e7f7a569d9d3954",
    ROOT/"scripts/exact_cohort_integrator.py": "092dc51161e708cdba3791ec47a7f8e843dd0e1afc866b6f65f286b8657d1745",
    ROOT/"scripts/quadratic_quadrature_study.py": "0eac4bdcfd2bb0b725940d15874fc762c16065d2a29f2a10d0da70c5725ac5dc",
}
KPHASE = BR/"first_wave_reduction_v1/independent_actual_preparation/actual_replay"
LPHASE = EX/"bridge_reduction_v1/independent_actual_preparation"
REVIEW_STATUS = "PASS_INDEPENDENT_COLD_NEGATIVE_KL_ACTIVE_ENERGY_TEST"


def find(rows, key, value):
    named, = (x for x in rows if x[key] == value)
    return named


def byte_pin(path, sha, pins):
    require(digest(path) == sha, f"stale source {path}")
    pins[path.relative_to(ROOT).as_posix()] = sha


def phase(dir, expected, status, pins, receipt_name="receipt.json"):
    o = load(dir/"os_exit_observation.json",expected,pins)
    t = load(dir/"terminal.json",o["terminal_sha256"],pins)
    r = load(dir/receipt_name,o["receipt_sha256"],pins)
    require(o["status"] == status and t["status"] == "FINISHED" and t["receipt_sha256"] == o["receipt_sha256"]
        and type(o["pid"]) is int and o["pid"] > 0 and type(t["pid"]) is int and t["pid"] == o["pid"]
        and type(r["pid"]) is int and r["pid"] == o["pid"] and type(o["observed_exit_code"]) is int
        and type(t["exit_code"]) is int and 0 == o["observed_exit_code"] == t["exit_code"]
        and "unified_exec_session" in o, "bound genuinely observed root proof completion")
    return r


def open_proofs(pins):
    ks = phase(KPHASE, "aa3c90101c64ccefe1403af2c296bcb22cafad3b2e35d73539aa7e3c5227f281",
               "PASS_ROOT_OBSERVED_FIRST_WAVE_INDEPENDENT_GLOBAL_EXIT", pins)
    kr = find(ks["case_results"],"case_id",K)
    k = load(ROOT/kr["native_file"], kr["native_sha256"],pins)
    require(ks["status"] == "COMPLETE_INDEPENDENT_FIRST_WAVE_SIGNED_ENDPOINT_REPLAY"
            and kr["status"] == k["status"] == "PASS_INDEPENDENT_FIRST_WAVE_EXACT_SPATIAL_ROOT"
            and k["pid"] == ks["pid"] and k["case_id"] == K and kr["strict_local_minimum"] is True
            and k["observer_and_optional_minimum"]["strict_minimum_proved"] is True
            and Q(k["observer_and_optional_minimum"]["normalized_H_euclidean_lower_exact"]) > 0,
            "proper negative K exact-spatial root and minimum")
    ls = load(LPHASE/"all_four_readback.json","859f507ad3762cfef44d53672bbe15a43592ce54d6e08ba8828bf50a925ec340",pins)
    lr = find(ls["cases"],"pilot_id",L)
    l = phase(LPHASE/"actual_replay"/L,lr["independent_os_observation_sha256"],
              "PASS_OBSERVED_INDEPENDENT_GLOBAL_REDUCTION_TRUE_OS_EXIT",pins)
    require(ls["status"] == "PASS_FOUR_INDEPENDENT_FIXED_INPUT_EXACT_SPATIAL_ROOTS_AND_MINIMA"
            and lr["independent_proof_sha256"] == digest(ROOT/lr["independent_proof_file"])
            and lr["independent_proof_sha256"] == pins[(LPHASE/"actual_replay"/L/"receipt.json").relative_to(ROOT).as_posix()]
            and lr["target_identity"] == l["target_identity"] and l["pilot_id"] == L
            and l["status"] == "PASS_INDEPENDENT_EXACT_SPATIAL_ROOT_AND_STRICT_LOCAL_MINIMUM"
            and Q(l["whole_box_normalized_H_euclidean_lower_exact"]) > 0, "proper negative L exact-spatial root and minimum")
    return k,l


def error(proof, a):
    r = Q(proof["radius_exact"])
    eta, z, e = (Q(a[k]) for k in ("eta_upper_exact","z_upper_exact","posterior_error_upper_exact"))
    require(r == RADIUS and a["strict_selfmap"] is True and eta >= 0 and 0 <= z < 1
            and eta+z*r == Q(a["selfmap_lhs_upper_exact"]) < r
            and 0 <= e == eta/(1-z) < r, "exact same-box spatial posterior")
    return e


def decode(packet):
    embedded = packet["frozen_binary_cache"]
    raw=base64.b64decode(embedded["data"],validate=True)
    require(sha256(raw).hexdigest() == embedded["sha256"], "exact binary state and mask bytes")
    names=("reference","u","triangles","free","active","energy_scale","surface_0_edges","surface_1_edges")
    with np.load(BytesIO(raw),allow_pickle=False) as data:
        result={k:data[k] for k in names}
    shapes={"reference":(325,2),"u":(611,),"triangles":(144,6),"free":(611,),
            "active":(144,1),"energy_scale":(1,),"surface_0_edges":(12,3),"surface_1_edges":(12,3)}
    require(all(result[k].shape==shapes[k] and np.isfinite(result[k]).all() for k in names),"archive array shapes")
    require(all(result[k].dtype == np.float64 for k in ("reference","u","active","energy_scale"))
            and all(result[k].dtype.kind in "iu" for k in ("free","triangles","surface_0_edges","surface_1_edges"))
            and set(result["active"].ravel()) == {0,1} and int(result["active"].sum()) == 24, "reference-bound binary flags")
    require(int(result["triangles"].min())>=0 and int(result["triangles"].max())<325
            and all(len(set(t))==6 for t in result["triangles"]), "valid P2 connectivity")
    return result


def input_target(name, fields, contract, row, box, packet, posterior):
    a = decode(packet)
    ref = [[Q(float(v)) for v in pos] for pos in a["reference"]]
    x = list(itertools.chain.from_iterable(ref))
    free = [int(j) for j in a["free"]]
    rim=Q(contract["material_parameters_binary_exact"]["outer_radius"])
    expected_free=[j for j in range(650) if x[j-j%2] != rim and (j%2 or x[j] != 0)]
    require(len(free)==611 and expected_free==free==box["free"] and box["triangles"]==a["triangles"].tolist(),
            "whole P2 geometry and fixed axis/rim restriction")
    center=x.copy()
    radii=[Q(0)]*650
    for j,u in zip(free,a["u"],strict=True):
        center[j]+=Q(float(u)); radii[j]=RADIUS
    require([[Q(v) for v in pos] for pos in box["reference_nodes_exact"]] == ref
            and [Q(v) for pos in box["current_centers_exact"] for v in pos] == center
            and [Q(v) for pos in box["coordinate_radii_exact"] for v in pos] == radii
            and Q(box["corrector_radius_exact"]) == RADIUS and 0<=posterior<RADIUS, "exact X+Pu proof center and radius")
    require(contract["pilot_id"]==name and contract["record_id"]==row["record_id"]==box["record_id"]==fields["record_id"]
            and contract["point_cache_sha256"]==row["point_cache_sha256"]==box["point_cache_sha256"]==fields["cache_sha"]
            and row["model_observation_contract_sha256"]==box["model_observation_contract_sha256"]==fields["contract_sha"]
            and contract["frozen_binary_cache_sha256"]==packet["frozen_binary_cache"]["sha256"],"same root proof and coordinate authority")
    d=contract["exact_requested_drive"]
    i,g,v=map(Q,(d["optical_dose_exact"],d["common_gain_exact"],d["activation_strength_exact"]))
    require(i==Q(row["dose_exact"])==Q(box["dose_exact"])==Q(packet["optical_dose_exact"])>0
         and g==Q(row["base_parameters_exact"]["common_gain"])==Q(packet["common_gain_exact"])>0
         and v==i*g==Q(fields["drive"])==Q(packet["activation_strength_exact"]),"same exact requested gain and dose")
    normalization=Q(float(a["energy_scale"][0]))
    require(normalization==Q(fields["scale"])==Q(contract["energy_scale_binary_exact"])
         ==Q(packet["normalization_binary_exact"])>0,"same exact positive normalization")
    require(Q(contract["material_parameters_binary_exact"]["reference_curvature"])==0,"flat material axial law")
    require(contract["material_parameters_binary_exact"]==packet["material_parameter_storage_exact"]
         and row["source_record_id"]==contract["physical_source_record_id"],"same case physical material and source record")
    return {"case_id":name,"id":row["source_record_id"],"dose":i,"gain":g,"normalization":normalization,
            "material":contract["material_parameters_binary_exact"],"effective_material":row["base_parameters_exact"],
            "reference":ref,"center":center,"radii":radii,"posterior":posterior,"a":a,
            "contract_sha256":fields["contract_sha"],"box_sha256":fields["box_sha"],
            "cache_sha256":fields["cache_sha"],"exact_drive":v}


def same_law(k,l):
    names=("reference","triangles","free","active","energy_scale","surface_0_edges","surface_1_edges")
    require(k["case_id"]==K and l["case_id"]==L and (k["id"],l["id"])==("optoGEF2_small_col_K","optoGEF2_small_col_L")
            and 0 < k["dose"] < l["dose"] and k["gain"] == l["gain"] > 0
            and k["normalization"] == l["normalization"] > 0 and k["material"]==l["material"]
            and k["effective_material"]==l["effective_material"],"same fixed material operator and increasing dose")
    require(all(k["a"][n].dtype==l["a"][n].dtype and k["a"][n].shape==l["a"][n].shape
                and k["a"][n].tobytes()==l["a"][n].tobytes() for n in names),
            "same mesh, flags, surfaces and boundary injection")
    return {n:sha256(k["a"][n].tobytes()).hexdigest() for n in names}


def inputs():
    pins={}
    for f,h in AUTHORITY.items():
        byte_pin(f,h,pins)
    k,l=open_proofs(pins)
    kp=EX/"remaining32_geometry_v1/preparation/protocol.json"
    lp=EX/"preparation/target_manifest.json"
    kt=find(load(kp,AUTHORITY[kp],pins)["trials"],"case_id",K)
    lt=find(load(lp,AUTHORITY[lp],pins)["targets"],"pilot_id",L)
    answers=[]
    for case,proof,t,box_key in ((K,k,kt,"trial_box"),(L,l,lt,"box_1")):
        cf,bf,rf=(ROOT/t[key] for key in ("model_contract",box_key,"row"))
        c,b,row=(load(path,t[key+"_sha256"],pins) for path,key in ((cf,"model_contract"),(bf,box_key),(rf,"row")))
        if case==K:
            q=proof["target"]
            require(q["trial_box_sha256"]==t["trial_box_sha256"],"K exact-spatial box")
            cache_sha, contract_sha, scale, drive=(q[m] for m in ("point_cache_sha256","model_contract_sha256",
                   "energy_scale_binary_exact","exact_drive"))
            error_bound=error(proof,proof["arithmetic"])
            proof_path=KPHASE/(K+".json")
        else:
            q=proof["target_identity"]
            for path,key in ((cf,"model_contract"),(bf,box_key)):
                require(proof["source_and_actual_sha256"][path.relative_to(ROOT).as_posix()]==t[key+"_sha256"],
                        "L independently checked exact-spatial box and contract")
            cache_sha, contract_sha, scale, drive=(q[m] for m in ("cache_packet_sha256","model_contract_sha256",
                  "energy_scale_binary_exact","drive_exact"))
            error_bound=error(proof,proof["independent_bounds"])
            proof_path=LPHASE/"actual_replay"/L/"receipt.json"
        data=load(ROOT/c["point_cache"],cache_sha,pins)
        fields={"cache_sha":cache_sha,"box_sha":t[box_key+"_sha256"],"contract_sha":contract_sha,
                "record_id":q["record_id"],"scale":scale,"drive":drive}
        a=input_target(case,fields,c,row,b,data,error_bound)
        a["exact_spatial_proof_sha256"]=pins[proof_path.relative_to(ROOT).as_posix()]
        answers.append(a)
    shared=same_law(*answers)
    return answers,shared,pins




def summarized(a):
    fields=("dose","gain","normalization","posterior","exact_drive")
    identities=("case_id","id","contract_sha256","box_sha256","cache_sha256","exact_spatial_proof_sha256")
    return {**{k:a[k] for k in identities}, **{k+"_exact":str(a[k]) for k in fields}}


def protocol(cases, shared, pins):
    return {"status":"FROZEN_TWO_CERTIFIED_NEGATIVE_ROOTS_ACTIVE_ENERGY_ORDER",
            "source_sha256":digest(Path(__file__)), "inputs_sha256":pins,
            "same_model_array_bytes_sha256":shared,"cases":[summarized(a) for a in cases],
            "same_effective_binary_material_and_exact_gain":True,
            "scope":"One necessary energy ordering between the certified negative K and L exact-spatial roots.",
            "output_units":"Scalar V divided by the shared positive 2*pi/S; normalized P2 material coordinates.",
            "cold_actual_energy_evaluations":0}


def evaluate(a):
    """Integrate only the active polynomial, on the actual posterior coordinate box."""
    center, e, free = a["center"],a["posterior"],a["radii"]
    Y=[[(center[j]-e,center[j]+e) if free[j] else (center[j],center[j])
        for j in (2*i,2*i+1)] for i in range(325)]
    result=Q(0),Q(0)
    entries=[]
    for e, (mask, nodes) in enumerate(zip(a["a"]["active"],a["a"]["triangles"],strict=True)):
        if int(mask[0]) == 0:
            continue
        coords=[a["reference"][int(i)] for i in nodes[:3]]
        enclosure=element(coords,[Y[int(i)] for i in nodes])
        result=add(result,enclosure)
        entries.append({"element":e,"unit_drive_energy_over_normalization_exact":list(map(str,enclosure))})
    require(len(entries)==24,"complete frozen support")
    return result,entries


def prepare():
    checks=fixtures()
    cases,shared,pins=inputs()
    checks.update(authenticated_two_same_operator_exact_spatial_roots=True,
                  exact_Y_centers_fixed_coordinates_ownership=True,
                  successful_true_OS_both_independent_proofs=True,
                  froze_before_actual_active_energy_evaluation=True)
    h=save(HERE/"protocol.json",protocol(cases,shared,pins))
    f=save(HERE/"preflight.json",{"status":"PASS","source_sha256":digest(Path(__file__)),
        "protocol_sha256":h,"checks":checks,"actual_energy_evaluations":0,
        "new_AD_FE_Hessian_force_solves_or_spatial_moment_jobs":0})
    print(json.dumps({"source_sha256":digest(Path(__file__)),"protocol_sha256":h,
                      "preflight_sha256":f,"current_evidence_pins":len(pins)}))


def run(args):
    source=digest(Path(__file__))
    p=load(HERE/"protocol.json",args.protocol_sha256,{})
    f=load(HERE/"preflight.json",args.preflight_sha256,{})
    g=load(HERE/"independent_prelaunch_review.json",args.review_sha256,{})
    expected={"source_sha256":source,"protocol_sha256":args.protocol_sha256,"preflight_sha256":args.preflight_sha256}
    for x in (f,g):
        require(type(x["checks"]) is dict and x["checks"] and all(v is True for v in x["checks"].values()),
                "complete literal true checks")
    require(g["status"]==REVIEW_STATUS and all(g[k]==v for k,v in expected.items())
            and f["source_sha256"]==source and f["status"]=="PASS" and f["protocol_sha256"]==args.protocol_sha256
            and type(f["actual_energy_evaluations"]) is int and f["actual_energy_evaluations"]==0,
            "matched independently reviewed cold phase")
    cases,shared,pins=inputs()
    require(p==protocol(cases,shared,pins),"exact unchanged input protocol")
    dest=HERE/"actual"
    dest.mkdir(exist_ok=False)
    stamp={**expected,"review_sha256":args.review_sha256,"pid":os.getpid()}
    save(dest/"started.json",{"status":"RUNNING",**stamp})
    terminal={"status":"FAILED","exit_code":1,**stamp}
    try:
        pairs=[evaluate(a) for a in cases]
        state, delta=classification(pairs[0][0],pairs[1][0])
        rows=[{**summarized(a),"unit_drive_energy_over_normalization_exact":list(map(str,b)),"active_elements":data}
              for a,(b,data) in zip(cases,pairs,strict=True)]
        h=save(dest/"results.json",{"status":"CANDIDATE_ACTIVE_ENERGY_ORDER_AWAITS_INDEPENDENT_READBACK",
                 **stamp,"classification":state,"V_L_minus_V_K_over_shared_normalization_exact":list(map(str,delta)),
                 "two_certified_same_operator_cases":rows,"actual_energy_evaluations":2,
                 "certified_loading_path":False,"source_fit_or_biological_dynamics_claimed":False})
        terminal.update(status="FINISHED",exit_code=0,results_sha256=h)
        print(json.dumps({"classification":state,"results_sha256":h}))
    finally:
        save(dest/"terminal.json",terminal)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="mode",required=True)
    sub.add_parser("prepare")
    run_parser=sub.add_parser("run")
    for arg in ("protocol","preflight","review"):
        run_parser.add_argument("--"+arg+"-sha256",required=True)
    args=parser.parse_args()
    prepare() if args.mode=="prepare" else run(args)

