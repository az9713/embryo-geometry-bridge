"""Independent fixed-record mathematical audit of C's old cube exclusion."""
from base64 import b64decode
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
from itertools import product
import json
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
COHORT=STAGE.parents[1]
ROOT=COHORT.parents[1]
CENTER=STAGE.parent
C="optoGEF2_small_col_C"
GLOBAL=CENTER/"independent_global_actual"/C
GLOBAL_PRE=CENTER/"independent_global_preparation"
DIRECT=COHORT/"remaining15_reduction_preparation"
ORIG=COHORT/"records"/C/"original_cube_v1"
PINS={
 STAGE/"candidate.json":"820fab9bf66259253624743791c0ddf013b684518b1424992b261521c8ee08cd",
 STAGE/"PROOF.md":"7bf8801b74a48d56cb4f487842d3491f0f69a40f1571ba75afd97960df0018d3",
 STAGE/"arithmetic_terminal.json":"ec3d0163f9d7f0b66a142bfebde23c225ea26596b18d3521dae5ced3703fac2e",
 GLOBAL/"receipt.json":"7cb318cd6dc5f879848b12c3dcdb86f668c7ad7227c7c58761dba6f37950a02d",
}


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_bytes())


def digest(value):
    return sha256(json.dumps(value,sort_keys=True,separators=(",",":")).encode()).hexdigest()


def p(x):
    lo,hi=map(Q,x)
    assert lo<=hi
    return lo,hi


def add(*xs):
    return tuple(sum(x) for x in zip(*xs))


def mul(a,b):
    values=[q*r for q,r in product(a,b)]
    return min(values),max(values)


def intersect(a,b):
    lo,hi=max(a[0],b[0]),min(a[1],b[1])
    assert lo<=hi
    return lo,hi


def main():
    pins={}
    def bind(path,h=None):
        h=h or sha(path)
        rel=path.relative_to(ROOT).as_posix()
        assert sha(path)==h==pins.setdefault(rel,h),rel
        return h
    for path,h in PINS.items():
        bind(path,h)
    d=read(STAGE/"candidate.json")
    t=read(STAGE/"arithmetic_terminal.json")
    assert t["status"]=="EXACT_COMPARISON_FINISHED_REVIEW_PENDING"
    assert t["candidate_sha256"]==PINS[STAGE/"candidate.json"] and t["written_proof_sha256"]==PINS[STAGE/"PROOF.md"]
    for f,k in (("protocol.json","protocol_sha256"),("preflight.json","preflight_sha256"),
                ("exclude.py","source_sha256"),("exclude_frozen.py","source_sha256")):
        bind(STAGE/f,d[k])
    frozen=read(STAGE/"protocol.json")
    pre=read(STAGE/"preflight.json")
    assert frozen["source_sha256"]==pre["source_sha256"]==d["source_sha256"]
    assert pre["protocol_sha256"]==d["protocol_sha256"] and all(x is True for x in pre["checks"].values())
    assert frozen["input_pins"]==d["input_pins"]
    for k,h in (d["input_pins"]|d["source_specific_predictor_cache_and_proof_files"]).items():
        bind(ROOT/k,h)
    protocol=read(ORIG/"protocol.json")
    box=read(ORIG/"coordinate_box.json")
    orig=protocol["row"]
    cell=orig["parameter_cell"]
    base=orig["base_parameters_exact"]
    assert orig["record_id"]==d["record_id"]==box["record_id"]==C
    assert cell==d["full_closed_parameter_cell_exact"] and base==d["base_parameters_exact"]
    ch=digest(cell)
    assert ch==d["closed_cell_sha256"]==protocol["parameter_cell_sha256"]==box["parameter_cell_sha256"]
    assert d["branch_domain_id"]==orig["branch_domain_id"] and box["domain_id"]==d["coordinate_domain_id"]
    assert orig["point_cache_sha256"]==d["same_R_binary64_target_cache_binding"]
    assert orig["predictor_sha256"]==d["same_predictor_sha256"]
    assert Q(orig["dose_exact"])==Q(d["dose_exact"])==Q(4089863,12500000)
    assert sha(ORIG/"coordinate_box.json")==protocol["coordinate_box_sha256"]
    assert box["root_membership_assumed"] is False and all(Q(x)==0 for x in box["shift_exact"])
    ig=read(ORIG/"integration_prelaunch_independent_review.json")
    bind(ORIG/"integration_prelaunch_independent_review.json")
    bind(ORIG/"domains.json",ig["domains_result_sha256"])
    bind(ORIG/"body_geometry.json",ig["body_sha256"])
    dom=read(ORIG/"domains.json")
    body=read(ORIG/"body_geometry.json")
    assert ig["status"]=="PASS" and all(v is True for v in ig["checks"].values())
    assert ig["protocol_sha256"]==sha(ORIG/"protocol.json")==dom["protocol_sha256"]
    assert ig["coordinate_box_sha256"]==dom["coordinate_box_sha256"]==protocol["coordinate_box_sha256"]
    assert ig["body_sha256"]==dom["body_sha256"] and body["global_injectivity_for_entire_box"] is True
    assert body["all_element_J_positive"] is True and Q(body["J_bound"]["lower_exact"])>Q(1,10**6)
    assert len(dom["surface_cases"])==24 and all(x["domain"]["status"]=="CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN" for x in dom["surface_cases"])

    inp=read(GLOBAL/"signed_inputs.json")
    review=read(GLOBAL/"receipt.json")
    os=read(GLOBAL/"os_exit_observation.json")
    for f,expected in ((GLOBAL_PRE/"protocol.json",review["review_protocol_sha256"]),
                       (GLOBAL_PRE/"preflight.json",review["review_preflight_sha256"]),
                       (GLOBAL_PRE/"review_centered_cohort_global_frozen.py",review["review_source_sha256"]),
                       (CENTER/(C+"_input_admission.json"),review["actual_admission_sha256"]),
                       (CENTER/"independent_actual"/(C+"_packet_review.json"),review["actual_packet_review_sha256"]),
                       (CENTER/"independent_actual"/(C+"_packet_review_os_observation.json"),review["actual_packet_review_os_sha256"]),
                       (CENTER/C/"results.json",review["producer_results_sha256"]),
                       (CENTER/C/"terminal.json",review["producer_terminal_sha256"])):
        bind(f,expected)
    assert os["observed_exit_code"]==0 and os["review_sha256"]==PINS[GLOBAL/"receipt.json"]
    assert review["future_signed_inputs_sha256"]==os["signed_inputs_sha256"]==sha(GLOBAL/"signed_inputs.json")
    assert review["status"]=="PASS_INDEPENDENT_CENTERED_REDUCTION_UNRESOLVED" and review["root_certified"] is False
    assert inp["same_closed_C_sha256"]==ch and inp["raw_gain_operand_proven_to_cover_exact_box"] is True
    assert review["gain_slots"]==review["producer_gain_contains_exact_rational_nodal_gain_slots"]==1656
    assert review["producer_gain_overlapping_only_slots"]==0 and inp["noncontaining_raw_gain_diagnostics"]==[]
    pi_lo,pi_hi=p(inp["gain_rational_pi_enclosure"])
    assert 3<pi_lo<=pi_hi<4 and pi_hi-pi_lo<=Q(1,2**1000)
    assert inp["declared_outward_dyadic_bits"]==800 and review["first_144_bulk_and_retained_24_surfaces"] is True

    old=read(DIRECT/"independent_actual"/(C+".json"))
    assert old["status"]=="PASS_DIRECT_REDUCTION_UNRESOLVED" and old["same_closed_C_sha256"]==ch
    prior_pins=old["producer_pins"]
    bind(DIRECT/C/"results.json",prior_pins["results_sha256"])
    bind(DIRECT/C/"terminal.json",prior_pins["terminal_sha256"])
    bind(ROOT/prior_pins["observation_file"],prior_pins["observation_sha256"])
    prior=read(DIRECT/C/"results.json")
    direct=[p(x) for x in prior["preconditioned_RF_intervals_exact"]]
    assert len(direct)==611
    order=inp["parameter_order"]
    assert len(order)==5 and order[-1]=="common_gain" and set(order)==set(cell)
    offsets=[tuple(Q(x)-Q(base[k]) for x in cell[k]) for k in order]
    assert offsets==list(map(p,inp["closed_five_offsets_exact"]))
    operands=(inp["point_bulk_RF_intervals_exact"],inp["RH_bulk_S_signed_intervals_exact"],
              inp["independent_affine_gain_RF_intervals_exact"],inp["retained_surface_RF_intervals_exact"])
    assert all(len(x)==611 for x in operands)
    alt=[]
    for point,hcol,rg,surf in zip(*operands,strict=True):
        assert len(hcol)==5
        coefficients=list(map(p,hcol))
        coefficients[-1]=add(coefficients[-1],p(rg))
        alt.append(add(p(point),p(surf),*(mul(c,a) for c,a in zip(coefficients,offsets,strict=True))))
    assert alt==list(map(p,inp["independent_coupled_predictor_RF_intervals_exact"]))
    final=[intersect(x,y) for x,y in zip(alt,direct,strict=True)]
    assert final==list(map(p,inp["independent_final_predictor_RF_intervals_exact"]))
    assert digest(inp["independent_final_predictor_RF_intervals_exact"])==d["signed_predictor_array_sha256"]

    cache=read(ROOT/orig["point_cache_file"])
    data=b64decode(cache["frozen_binary_cache"]["data"],validate=True)
    assert sha256(data).hexdigest()==cache["frozen_binary_cache"]["sha256"]==orig["frozen_binary_cache_sha256"]
    with np.load(BytesIO(data),allow_pickle=False) as f:
        a={k:f[k] for k in ("preconditioner","reference","free","u","triangles")}
    with np.load(ROOT/orig["predictor_file"],allow_pickle=False) as f:
        S=f["S"]
    R=a["preconditioner"]
    assert S.shape==(611,5) and R.shape==(611,611)
    assert np.isfinite(S).all() and np.isfinite(R).all()
    norm=max(sum(abs(Q.from_float(float(x))) for x in row) for row in R)
    assert norm==Q(protocol["R_infinity_norm_exact"])==Q(inp["normalization_norm_exact"])
    free=a["free"].tolist()
    assert free==box["free"] and len(set(free))==611
    assert a["triangles"].tolist()==box["triangles"]
    r=Q(box["corrector_radius_exact"])
    assert r==Q(d["corrector_radius_in_normalized_coordinate_exact"])==Q(protocol["normal_radius_exact"])==Q(1,10**11)
    assert r==Q(inp["same_radius_exact"])
    centers=[[Q.from_float(float(x)) for x in row] for row in a["reference"]]
    assert centers==[list(map(Q,row)) for row in box["reference_nodes_exact"]]
    radii=[[Q(0),Q(0)] for _ in centers]
    for i,f in enumerate(free):
        node,axis=divmod(f,2)
        lo,hi=add(*(mul((Q(float(S[i,k])),)*2,offsets[k]) for k in range(5)))
        assert (lo,hi)==p(box["predictor_hulls_exact"][i])
        centers[node][axis]+=Q(float(a["u"][i]))
        radii[node][axis]=r+max(-lo,hi)
    assert centers==[list(map(Q,row)) for row in box["current_centers_exact"]]
    assert radii==[list(map(Q,row)) for row in box["coordinate_radii_exact"]]
    source=protocol["cached_bound_provider"]
    c_audit=read(ROOT/source)
    capsule=read(ROOT/orig["independent_root_capsule_file"])
    (wrapper,)=[e for e in capsule["record_audits"] if e["record_id"]==C]
    assert wrapper["audit"]==c_audit and wrapper["review_file"]==source and wrapper["review_sha256"]==sha(ROOT/source)
    root=read(ROOT/orig["root_certificate_file"])
    assert c_audit["status"]=="PASS" and c_audit["uniform_quantifier_verified"] is True
    assert all(x is True for x in c_audit["checks"].values())
    assert c_audit["certificate_sha256"]==orig["root_certificate_sha256"]
    assert root["point_cache_file"]==orig["point_cache_file"] and root["predictor_file"]==orig["predictor_file"]
    assert c_audit["branch_domain_id"]==root["branch_domain_id"]==d["branch_domain_id"]
    assert c_audit["parameter_cell_sha256"]==root["parameter_cell_sha256"]==ch
    oldb=c_audit["bounds"]
    assert Q(oldb["radius"]["lower_exact"])<=r<=Q(oldb["radius"]["upper_exact"])
    eta=max(Q(oldb["eta"]["upper_exact"]),Q(root["uniform_cell_result"]["eta_uniform"]["upper_exact"]))
    zQ=Q(oldb["Z"]["upper_exact"])
    s=review["scalar_exact_readback"]
    assert eta==Q(protocol["cached_eta_upper_exact"])==Q(s["old_cached_eta_exact"])==Q(old["exact_scalars"]["old_cached_eta_upper_exact"])
    assert zQ==Q(protocol["cached_Z_upper_exact"])==Q(s["old_cached_Z_exact"])==Q(old["exact_scalars"]["old_cached_Z_upper_exact"])
    b=Q(s["same_fullbox_b_exact"])
    z=zQ+b
    assert b==Q(old["exact_scalars"]["b_upper_exact"])==Q(inp["same_full_box_b_exact"])
    assert z==Q(old["exact_scalars"]["total_Z_upper_exact"])==Q(inp["same_total_Z_upper_exact"])
    index=max(range(611),key=lambda j:max(final[j][0],-final[j][1],Q(0)))
    witness=final[index]
    distance=max(witness[0],-witness[1],Q(0))
    lower=distance-eta
    cap=(1+z)*r
    margin=lower-cap
    expected={"preconditioned_defect_norm_lower":distance,"cached_predictor_residual_upper":eta,
              "integral_predictor_residual_lower":lower,"cached_derivative_defect_upper":zQ,
              "whole_box_H_defect_upper":b,"integral_derivative_defect_upper":z,
              "predictor_residual_upper_necessary_if_old_cube_root":cap,"strict_exclusion_margin":margin}
    assert all(Q(d["comparison_exact"][k])==v for k,v in expected.items())
    assert d["comparison_exact"]["witness_free_row_zero_based"]==index==6
    assert p(d["comparison_exact"]["witness_signed_defect_interval_exact"])==witness
    assert witness[0]>0 and margin>0 and d["strict_exact_comparison_passed"] is True
    assert d["existing_sufficient_test_unchanged"] is True
    assert d["existence_in_any_other_domain_asserted"] is d["biological_source_incompatibility_asserted"] is False
    assert all(sha(ROOT/k)==v for k,v in pins.items())
    receipt={
      "status":"PASS_INDEPENDENT_UNIFORM_C_OLD_CUBE_EXCLUSION",
      "record_id":C,"closed_cell_sha256":ch,"scope":d["scope"],
      "checks":{"whole_same_C_binary_S_R_dose_cached_residual_and_geometry_bound":True,
                "all_611_signed_predictor_defect_intersections_recomposed":True,
                "original_full_H_bound_and_geometry_apply_on_same_theta_segment":True,
                "strict_exact_contradiction_of_necessary_root_bound":True,
                "finite_FE_and_only_original_matched_cube_nonexistence":True},
      "old_matched_cube_integral_root_excluded":True,"new_root_certified":False,
      "conclusion":"For every theta in the full closed cell, no zero of normalized exact-integral F_I lies in the original matched closed radius-r predictor cube.",
      "math":"At any fixed theta, the positive signed predictor defect lower bound and cached predictor force upper imply ||RF_I(p)||∞ >= d_lower - eta_Q. If F_I(u*)=0 in the matched cube, the same-theta segment from predictor to u* stays in the matched cube and the reviewed whole nodal box. The derivative norm is at most 1+z because RH_I=I-(I-RH_I); segment integration gives ||RF_I(p)||∞ <= (1+z)r. Strict contradiction. Neither SPD nor invertibility nor failure of any sufficient fixed-point test is a premise.",
      "reuse":"Reviewed point bulk, bulk-only derivative path, full whole-box derivative including 24 surfaces, centered alternate gain via 1024-bit pi and exact nodal box. No geometry kernel or mechanical solve re-executed.",
      "witness_free_row_zero_based":index,"witness_interval_exact":list(map(str,witness)),
      "comparison_exact":{k:str(v) for k,v in expected.items()},
      "comparison_display":{k:float(v) for k,v in expected.items()},
      "new_force_or_integral_calculation":False,"global_absence_claim":False,"physical_or_source_effect_claim":False,
      "input_and_source_sha256":pins,"audit_source_sha256":sha(Path(__file__)),
    }
    o=HERE/"receipt.json"
    with o.open("x",encoding="utf-8") as out:
        json.dump(receipt,out,indent=2)
        out.write("\n")
    print(sha(o),len(pins),float(lower),float(cap),float(margin))


if __name__=="__main__":
    main()

