"""Exact original-cube exclusion candidate for record C, using reviewed endpoints only."""
import argparse
from collections import Counter
from datetime import datetime,timezone
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
BASE="runs/exact_integration_cohort_v1"
REC="optoGEF2_small_col_C"
ORIG=f"{BASE}/records/{REC}/original_cube_v1"
GLOBAL=f"{BASE}/remaining15_centered_reduction_preparation/independent_global_actual/{REC}"
DIRECT=f"{BASE}/remaining15_reduction_preparation/independent_actual/{REC}.json"
CELL="64f93fefa034f4c004b93451d98b2fd04e0072855be8c5e3bba9babb8f394dbf"
PINS={
    f"{GLOBAL}/signed_inputs.json":"1e7842b0b58c3298e65edbaf4db21247acaa9144ba330e7d9a483b1bed9149eb",
    f"{GLOBAL}/receipt.json":"7cb318cd6dc5f879848b12c3dcdb86f668c7ad7227c7c58761dba6f37950a02d",
    f"{GLOBAL}/os_exit_observation.json":"857b6b6f1de596501bbb92d236471cfbc892f251881b03d182f2ab3765a637c3",
    DIRECT:"b66b821666c74ee68d2c0105df3087c4bec9e3a408d92f14cfa88231c9fc7d08",
    f"{ORIG}/protocol.json":"775b87588ae463396a407102b4672a17e9946957b7a07779097fccdc487a634e",
    f"{ORIG}/coordinate_box.json":"8d1595f2ed493a2a3d90b559a2690005073943ce84b1b58a00c19c41c2617bf1",
}
CONCLUSION="All parameters in the same closed cell and the original matched predictor cube of record C, for the normalized exact-spatial-integration fixed P2 force."


def require(value,msg):
    if not value:
        raise ValueError(msg)


def hash_data(raw):
    return sha256(raw).hexdigest()


def file_hash(f):
    return hash_data(f.read_bytes())


def canonical(obj):
    return hash_data(json.dumps(obj,sort_keys=True,separators=(",",":"),default=str).encode())


def bound(rel,pin,parse=True):
    path=(ROOT/rel).resolve()
    require(path.is_relative_to(ROOT) and path.relative_to(ROOT).as_posix()==rel,"exact canonical study path")
    raw=path.read_bytes()
    require(hash_data(raw)==pin,f"changed proof input {rel}")
    return json.loads(raw) if parse else raw


def write(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False,default=str)+"\n").encode()
    with path.open("xb") as out:
        out.write(raw)
    return hash_data(raw)


def exact(value):
    require(type(value) in (Q,int,str),"exact rational text or integer required")
    return Q(value)


def interval(value):
    require(len(value)==2,"two exact interval endpoints required")
    lo,hi=map(exact,value)
    require(lo<=hi,"ordered interval required")
    return lo,hi


def zero_distance(value):
    lo,hi=interval(value)
    return lo if lo>0 else -hi if hi<0 else Q(0)


def comparison(signed_defect,eta_Q,z,r):
    require(signed_defect and eta_Q>=0 and 0<=z<1 and r>0,"reviewed nonnegative constants and positive normal radius")
    row=max(range(len(signed_defect)),key=lambda j:zero_distance(signed_defect[j]))
    d=zero_distance(signed_defect[row])
    lower=d-eta_Q
    maximum_at_anchor_if_root=(1+z)*r
    return row,d,lower,maximum_at_anchor_if_root,lower-maximum_at_anchor_if_root,lower>maximum_at_anchor_if_root


def tests():
    checks={}
    def ck(key,yes):
        require(yes,key)
        checks[key]=True
    ck("positive_negative_and_zero_crossing_distance",
       [zero_distance(v) for v in ((1,4),(-5,-2),(-2,4),(0,2))]==[1,2,0,0])
    result=comparison([(Q(-1),Q(2)),(Q(-7),Q(-4)),(Q(2),Q(3))],Q(1,4),Q(1,2),Q(2))
    ck("largest_component_lower_bound_valid_for_both_signs",result==(1,Q(4),Q(15,4),Q(3),Q(3,4),True))
    ck("threshold_equality_is_not_exclusion",
       comparison([(Q(2),Q(2))],Q(1,2),Q(1,2),Q(1))[-2:]==(Q(0),False))
    ck("failed_self_map_sufficiency_alone_is_no_exclusion",
       comparison([(Q(-2),Q(2))],Q(0),Q(1,4),Q(1))[-1] is False)
    ck("same_row_bound_beats_signed_cancellation",
       comparison([(Q(2),Q(3))],Q(1,4),Q(1,2),Q(1))[-2:]==(Q(1,4),True))
    try:
        comparison([(Q(2),Q(1))],Q(0),Q(1,4),Q(1))
    except ValueError:
        checks["reversed_intervals_rejected"]=True
    else:
        raise AssertionError("accepted reversed interval")
    for eta,z,r in ((Q(-1),Q(0),Q(1)),(Q(0),Q(1),Q(1)),(Q(0),Q(0),Q(0))):
        try:
            comparison([(0,1)],eta,z,r)
        except ValueError:
            continue
        raise AssertionError("accepted invalid scalar")
    checks["invalid_eta_z_r_rejected"]=True
    return checks


def prepare():
    for path,h in PINS.items():bound(path,h,False)
    checks=tests()
    own=Path(__file__).read_bytes()
    with (HERE/"exclude_frozen.py").open("xb") as f:f.write(own)
    p=write(HERE/"protocol.json",{
        "status":"FROZEN_EXACT_C_OLD_CUBE_EXCLUSION_PREPARATION",
        "source_sha256":hash_data(own),"input_pins":PINS,"record_id":REC,"closed_cell_sha256":CELL,
        "scope":CONCLUSION,"independent_predictor_defect_operand":"independent_final_predictor_RF_intervals_exact",
        "cached_bound":"reviewed upper bound on ||R F_Q(p(theta),theta)|| on same family",
        "derivative_bound":"reviewed upper bound on ||I-R H_I|| on entire old predictor cube family",
        "comparison":"max_j distance(0,[RD(p)]_j) - eta_Q > (1+z)r, strict rational endpoints",
        "assumptions":"Original P2 free coordinates, fixed original cached predictor u0+S(theta-theta0), reviewed same signed R, all 611 free coordinates, same branch and dose.",
        "classification_policy":"Exact comparison only. Even if strict, candidate pending separate independent proof review. Existing sufficient-test result is left intact.",
        "new_force_cache_polynomial_integral_or_parameter_search":False})
    pf=write(HERE/"preflight.json",{"status":"PASS_C_OLD_CUBE_EXCLUSION_PREPARATION_ONLY","source_sha256":hash_data(own),
          "protocol_sha256":p,"checks":checks,"actual_comparison_executed":False})
    print(json.dumps({"source_sha256":hash_data(own),"protocol_sha256":p,"preflight_sha256":pf}))


def apply(args):
    own=Path(__file__).read_bytes()
    pd=bound((HERE/"protocol.json").relative_to(ROOT).as_posix(),args.protocol_sha256)
    pf=bound((HERE/"preflight.json").relative_to(ROOT).as_posix(),args.preflight_sha256)
    require(pd["status"]=="FROZEN_EXACT_C_OLD_CUBE_EXCLUSION_PREPARATION"
            and pf["status"]=="PASS_C_OLD_CUBE_EXCLUSION_PREPARATION_ONLY"
            and pf["protocol_sha256"]==args.protocol_sha256 and pf["source_sha256"]==pd["source_sha256"]==hash_data(own)
            and pf["checks"]==tests() and pd["input_pins"]==PINS and pd["record_id"]==REC
            and (HERE/"exclude_frozen.py").read_bytes()==own,"same independently readable frozen comparison")
    data={n:bound(n,h) for n,h in PINS.items()}
    inp=data[f"{GLOBAL}/signed_inputs.json"]
    rev=data[f"{GLOBAL}/receipt.json"]
    os=data[f"{GLOBAL}/os_exit_observation.json"]
    old=data[DIRECT]
    protocol=data[f"{ORIG}/protocol.json"]
    box=data[f"{ORIG}/coordinate_box.json"]
    require(inp["status"]=="INDEPENDENT_EXACT_CENTERED_INPUT_ARRAYS"
            and rev["status"]=="PASS_INDEPENDENT_CENTERED_REDUCTION_UNRESOLVED" and rev["root_certified"] is False
            and rev["future_signed_inputs_sha256"]==PINS[f"{GLOBAL}/signed_inputs.json"]
            and os["status"]=="OBSERVED_ZERO_EXIT_INDEPENDENT_GLOBAL_CENTERED_REPLAY"
            and type(os["observed_exit_code"]) is int and os["observed_exit_code"]==0
            and os["review_sha256"]==PINS[f"{GLOBAL}/receipt.json"]
            and os["signed_inputs_sha256"]==PINS[f"{GLOBAL}/signed_inputs.json"]
            and old["status"]=="PASS_DIRECT_REDUCTION_UNRESOLVED",
            "same independent antecedent checks and true global readback outcome")
    row=protocol["row"]
    cell=row["parameter_cell"]
    base=row["base_parameters_exact"]
    require(all(x["record_id"]==REC for x in (inp,rev,os,old,protocol,row,box))
            and canonical(cell)==CELL==inp["same_closed_C_sha256"]==old["same_closed_C_sha256"]
            ==protocol["parameter_cell_sha256"]==box["parameter_cell_sha256"]==row["parameter_cell_sha256"]
            and protocol["coordinate_box_sha256"]==PINS[f"{ORIG}/coordinate_box.json"]
            and protocol["coordinate_box_file"]==f"{ORIG}/coordinate_box.json"
            and all(exact(x)==0 for x in box["shift_exact"])
            and box["root_membership_assumed"] is False and protocol["root_asserted"] is False
            and box["domain_id"]==f"exactFE:{REC}:original_cube_v1","one closed unchanged C and original cube only")
    require(len(inp["parameter_order"])==len(cell)==len(base)==5 and set(inp["parameter_order"])==set(cell)
            and all([str(exact(t)-exact(base[k])) for t in cell[k]]==ds
                    for k,ds in zip(inp["parameter_order"],inp["closed_five_offsets_exact"],strict=True)),
            "same fixed theta offsets for all predictor rows")
    scal=rev["scalar_exact_readback"]
    old_scal=old["exact_scalars"]
    eta,z0,b,r=map(exact,(scal["old_cached_eta_exact"],scal["old_cached_Z_exact"],
                        scal["same_fullbox_b_exact"],scal["same_radius_exact"]))
    z=exact(scal["total_Z_upper_exact"])
    require(r==Q(1,10**11)==exact(box["corrector_radius_exact"])==exact(protocol["normal_radius_exact"])
            ==exact(inp["same_radius_exact"])==exact(old_scal["old_radius_exact"])
            and eta==exact(old_scal["old_cached_eta_upper_exact"])==exact(protocol["cached_eta_upper_exact"])
            and z0==exact(old_scal["old_cached_Z_upper_exact"])==exact(protocol["cached_Z_upper_exact"])
            and b==exact(old_scal["b_upper_exact"])==exact(inp["same_full_box_b_exact"])
            and z==z0+b==exact(inp["same_total_Z_upper_exact"])==exact(old_scal["total_Z_upper_exact"])
            and exact(scal["exact_R_infinity_norm"])==exact(old_scal["R_infinity_norm_exact"])
            ==exact(inp["normalization_norm_exact"])==exact(protocol["R_infinity_norm_exact"])
            and (rev["global_F_entries"],rev["H_upper_entries"])==(611,6676)
            and len(box["free"])==len(set(box["free"]))==611,"same original derivative, projection, radius and cached predictor residual")
    require(old["root_certified"] is False and inp["raw_gain_operand_proven_to_cover_exact_box"] is True
            and rev["producer_gain_contains_exact_rational_nodal_gain_slots"]==rev["gain_slots"]==1656
            and rev["producer_gain_overlapping_only_slots"]==0,"independent exact-box gain and signed endpoint scope")
    cached_source=protocol["cached_bound_provider"]
    bound_files={}
    for prefix in ("predictor","point_cache","proof_center_state","root_certificate","independent_root_capsule"):
        path,expected=row[prefix+"_file"],row[prefix+"_sha256"]
        require(protocol["source_and_input_hashes"][path]==expected,"canonical original predictor/cache/premise provenance")
        bound(path,expected,False)
        bound_files[path]=expected
    cached_digest=protocol["source_and_input_hashes"][cached_source]
    ca=bound(cached_source,cached_digest)
    capsule=bound(row["independent_root_capsule_file"],row["independent_root_capsule_sha256"])
    wrapper,=[v for v in capsule["record_audits"] if v["record_id"]==REC]
    require(wrapper["review_file"]==cached_source and wrapper["review_sha256"]==cached_digest
            and wrapper["audit"]==ca,"same independent C cached-bound provider")
    require(ca["status"]=="PASS" and ca["uniform_quantifier_verified"] is True
            and all(x is True for x in ca["checks"].values())
            and ca["certificate_sha256"]==row["root_certificate_sha256"]
            and ca["branch_domain_id"]==row["branch_domain_id"] and ca["parameter_cell_sha256"]==CELL,
            "same independently reviewed cached root premises, never substitution of quadrature root")
    root=bound(row["root_certificate_file"],row["root_certificate_sha256"])
    require(root["record_id"]==REC and root["branch_domain_id"]==row["branch_domain_id"]
            and root["predictor_file"]==row["predictor_file"] and root["point_cache_file"]==row["point_cache_file"]
            and root["parameter_cell_sha256"]==CELL and exact(root["dose_exact"])==exact(row["dose_exact"]),
            "one original cached predictor family")
    ca_lo,ca_hi=interval([ca["bounds"]["radius"][k] for k in ("lower_exact","upper_exact")])
    require(ca_lo<=r<=ca_hi and exact(ca["bounds"]["Z"]["upper_exact"])==z0
            and max(exact(ca["bounds"]["eta"]["upper_exact"]),exact(root["uniform_cell_result"]["eta_uniform"]["upper_exact"]))==eta,
            "exact original audited radius inclusion and cached residual/norm bounds")
    values=inp["independent_final_predictor_RF_intervals_exact"]
    require(len(values)==611 and rev["independent_gain_predictor_passes"] is False,"full independent predictor-only signed rows")
    j,d,lower,cap,margin,ok=comparison(values,eta,z,r)
    require(ok and margin>0,"strict old-cube exclusion condition did not pass")
    witness=interval(values[j])
    status="CANDIDATE_UNIFORM_OLD_CUBE_EXCLUSION_INDEPENDENT_REVIEW_PENDING"
    candidate={"status":status,"record_id":REC,"closed_cell_sha256":CELL,
          "dose_exact":row["dose_exact"],"branch_domain_id":row["branch_domain_id"],"coordinate_domain_id":box["domain_id"],
          "full_closed_parameter_cell_exact":cell,"base_parameters_exact":base,"corrector_radius_in_normalized_coordinate_exact":r,
          "profile":"p(theta) = u0 + S(theta - theta0); ||u-p(theta)||_infinity <= r; exact original unshifted family",
          "same_R_binary64_target_cache_binding":row["point_cache_sha256"],"same_predictor_sha256":row["predictor_sha256"],
          "comparison_exact":{"witness_free_row_zero_based":j,"witness_signed_defect_interval_exact":witness,
              "preconditioned_defect_norm_lower":d,"cached_predictor_residual_upper":eta,"integral_predictor_residual_lower":lower,
              "cached_derivative_defect_upper":z0,"whole_box_H_defect_upper":b,"integral_derivative_defect_upper":z,
              "predictor_residual_upper_necessary_if_old_cube_root":cap,"strict_exclusion_margin":margin},
          "strict_exact_comparison_passed":ok,
          "source_sha256":hash_data(own),"protocol_sha256":args.protocol_sha256,"preflight_sha256":args.preflight_sha256,
          "input_pins":PINS,"source_specific_predictor_cache_and_proof_files":{**bound_files,cached_source:cached_digest},
          "signed_predictor_array_sha256":canonical(values),
          "scope":CONCLUSION,"independent_exclusion_review_complete":False,"existing_sufficient_test_unchanged":True,
          "existence_in_any_other_domain_asserted":False,"material_law_rejected":False,"biological_source_incompatibility_asserted":False,
          "new_mechanics_or_force_computation":False}
    h=write(HERE/"candidate.json",candidate)
    proof = rf"""# Candidate exclusion for record C’s original predictor cube

**Exact comparison passed. Independent review is pending.**

Let \(F_I(u,\theta)\) be the normalized exact-spatial-integration force of the fixed P2 mechanical model and \(F_Q(u,\theta)\) the matching cached quadrature force. Both use the selected source dose \(4089863/12500000\). The original predictor is \(p_\theta=u_0+S(\theta-\theta_0)\) in 611 normalized free coordinates. For every parameter \(\theta\) in the same closed cell, the original cube is \(\|u-p_\theta\|_\infty\le r\), with \(r=10^{{-11}}\).

The independently replayed signed predictor intervals bound \(R(F_I-F_Q)(p_\theta,\theta)\) over that entire parameter cell. Row {j} in zero-based free-coordinate order lies in \([{float(witness[0]):.15g},{float(witness[1]):.15g}]\) to display precision, strictly away from zero. Its distance from zero implies \(\|R(F_I-F_Q)(p_\theta,\theta)\|_\infty\ge d_*\) for every \(\theta\). The cached predictor residual has the uniform upper bound \(\|RF_Q(p_\theta,\theta)\|_\infty\le\eta_Q\). The reverse triangle inequality gives
\[
\|RF_I(p_\theta,\theta)\|_\infty\ge d_* - \eta_Q.
\]

Suppose \(u_\star\) is a zero of \(F_I\) in that same cube at the same \(\theta\). The line from \(p_\theta\) to \(u_\star\) stays inside the reviewed domain. The full integral Hessian \(H_I\) includes the retained 24 surfaces, and the reviewed whole-box bound is \(\|I-RH_I\|_\infty\le z\). Hence \(\|RH_I\|_\infty\le1+z\), and integrating the derivative on the line yields
\[
\|RF_I(p_\theta,\theta)\|_\infty
=\|R(F_I(p_\theta,\theta)-F_I(u_\star,\theta))\|_\infty
\le(1+z)\,r.
\]

The exact comparison in [candidate.json](candidate.json) gives \(d_* - \eta_Q\approx {float(lower):.15g}\), while \((1+z)r\approx {float(cap):.15g}\). The exact difference is strictly positive, approximately \({float(margin):.15g}\). Thus the two inequalities are incompatible for every \(\theta\) in the full originally declared closed cell.

This exclusion argument uses a positive residual lower bound. The separate centered sufficient test keeps its unresolved status. The radius is a tolerance in normalized free coordinates. This argument gives no existence or exclusion result elsewhere, no source-data or biological incompatibility, and no rejection of the material law. The candidate uses the prior reviewed force and derivative enclosures and performs no new field computation.
"""
    with (HERE/"PROOF.md").open("x",encoding="utf8",newline="\n") as p:p.write(proof)
    term=write(HERE/"arithmetic_terminal.json",{"status":"EXACT_COMPARISON_FINISHED_REVIEW_PENDING","candidate_sha256":h,
                "written_proof_sha256":file_hash(HERE/"PROOF.md"),"utc":datetime.now(timezone.utc).isoformat()})
    print(json.dumps({"status":status,"candidate_sha256":h,"written_proof_sha256":file_hash(HERE/"PROOF.md"),
                      "arithmetic_terminal_sha256":term,"witness_row":j,"lower":float(lower),"necessary_upper":float(cap),
                      "margin":float(margin)}))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="mode",required=True)
    sub.add_parser("prepare")
    actual=sub.add_parser("compare")
    for flag in ("protocol-sha256","preflight-sha256"):actual.add_argument("--"+flag,required=True)
    args=parser.parse_args()
    prepare() if args.mode=="prepare" else apply(args)

