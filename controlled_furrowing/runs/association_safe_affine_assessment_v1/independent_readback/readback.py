"""Independent exact affine and metadata readback. Never imports assessor or an optimizer."""
import argparse
import ast
import csv
import json
import math
import os
from pathlib import Path
from fractions import Fraction
from hashlib import sha256
import traceback

HERE=Path(__file__).resolve().parent
ASSESS=HERE.parent
ROOT=ASSESS.parents[1]
SOURCE=Path(__file__).resolve()
PINS={
 ASSESS/"assess.py":"065759b33106c8909ed9d7b1658dc422e1f1125628da7fd78684d9f5c106b598",
 ASSESS/"protocol.json":"e409bd1d039c6cefbf87c67d98c8de433bfbbdc0b9b1caa7bacbf0fcb0da69d3",
 ASSESS/"inputs.json":"116a19b185ff743a2ccc4f9d137bfe2d95ab98fef09af54185f438a77cbfbad9",
 ASSESS/"results.json":"9e2d67bd0eba02b300f9a8936513c046b7e170abdb51422eaff77edcf764d517",
 ASSESS/"ASSESSMENT.md":"fe18a94d058627890d3fca6bafd3fd0073007e64547f61ddf21d9d1df3f3892d",
 ROOT/"evidence/joint_observation_contract_audit.json":"3fa58c360b04bdadb718586511d75bcc25a8292113455e7f6d95615c929ce2c1",
 ROOT/"scripts/material_cohort_descent.py":"8343f4908a76a47f8c83159e7befbbaa23afdc7d4a5f4689bda9ec901e3957c1",
 ROOT/"scripts/transfer_ambiguity_nonlinear_probe.py":"87ca09b46f9569903276923d937bef63d00bba9ad1b99d78b7e1708168f7e631",
}
SMALL="runs/material_cohort_descent_v1/results.json"
SP="runs/material_cohort_descent_v1/protocol.json"
LARGE="runs/transfer_ambiguity_linear_v1/actual/results.json"
LP="runs/transfer_ambiguity_linear_v1/preparation/protocol.json"
OLD_REVIEW="runs/transfer_ambiguity_linear_v1/actual_independent/receipt.json"
CSV="data/derived/calibration_small_18um.csv"
CHANNELS=("maximum_shortening_um","region_averaged_shortening_um","apical_max_shortening_um","basal_max_shortening_um")
OBS=("current_disk__spatial_max","current_disk__regional_mean","current_radial_line__spatial_max","current_radial_line__regional_mean")
COLS=("log_apical_area","log_basal_area","apical_tension_additive","log_basal_tension","log_common_gain")
LETTERS=tuple("BCDEFGHIJKLMNOPQR")
NAMES=("max_contrast_change","max_negative_contrast_change")
EPS=Fraction(1,1000)
SCOPE="Exact rational LP verification on frozen binary64-rounded tangents, and source-contract summary linkage audit."


def need(ok,msg):
    if not ok:
        raise ValueError(msg)


def rel(path):
    return path.resolve().relative_to(ROOT).as_posix()


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def raw(path, expected=None):
    path=path.resolve()
    need(path.is_relative_to(ROOT),"path outside study")
    b=path.read_bytes()
    if expected is not None:
        need(sha256(b).hexdigest()==expected,f"authority digest mismatch {rel(path)}")
    return b


def read(path,expected=None):
    return json.loads(raw(path,expected))


def persist(path,value):
    with path.open("x",encoding="utf8",newline="\n") as f:
        json.dump(value,f,indent=2,allow_nan=False)
        f.write("\n")
    return digest(path)


def q(x):
    need(type(x) in (str,int,Fraction),"exact rationals require typed endpoints")
    return Fraction(x)


def pair(values):
    need(len(values)==2,"two endpoints")
    a,b=map(q,values)
    need(a<=b,"reversed endpoints")
    return a,b


def midpoint(values):
    return sum(pair(values),Fraction(0))/2


def binary(x):
    return Fraction.from_float(float(x))


def dot(a,b):
    return sum((x*y for x,y in zip(a,b,strict=True)),Fraction(0))


def exact_rank(matrix):
    rows=[list(r) for r in matrix]
    rank=0
    for k in range(len(rows[0])):
        best=next((i for i in range(rank,len(rows)) if rows[i][k]),None)
        if best is None:
            continue
        rows[rank],rows[best]=rows[best],rows[rank]
        pivot=rows[rank][k]
        for i in range(rank+1,len(rows)):
            if rows[i][k]:
                factor=rows[i][k]/pivot
                rows[i]=[x-factor*y for x,y in zip(rows[i],rows[rank],strict=True)]
        rank+=1
    return rank


def input_protocol():
    for path,h in PINS.items():
        raw(path,h)
    p=read(ASSESS/"protocol.json",PINS[ASSESS/"protocol.json"])
    need(p["status"]=="FROZEN_ARCHIVED_AFFINE_COMPARISON_PROTOCOL" and
         p["source_sha256"]==PINS[ASSESS/"assess.py"] and "completed read-only comparison" in p["design_timing"],
         "retrospective reproducibility protocol, no new preregistration")
    actual=p["study_relative_sha256"]
    for path,h in actual.items():
        raw(ROOT/path,h)
    return p,actual


def source_table(pins):
    audit=read(ROOT/"evidence/joint_observation_contract_audit.json",PINS[ROOT/"evidence/joint_observation_contract_audit.json"])
    need(audit["mean_duplicate_exact"] is True and audit["dose_max_duplicate_in_Fig4h_small_exact"] is True
         and type(audit["mean_values_unique_count"]) is int and audit["mean_values_unique_count"]==17,
         "verified exact duplicate scalar keys in archived small-block audit")
    pairs=audit["pairing"]
    need(pairs["Fig2c_to_2d_verified_by_unique_mean"] is True
         and pairs["Fig2e_withincolumn_apical_basal_pair_author_supported"] is True
         and "not independently verified specimen IDs" in pairs["Fig2e_to_otherpanels_canonical_order"],
         "only published scalar joins, component pairing is within-panel")
    h=audit["workbook_headers"]
    need(h["Fig 2c"]["A4"]=="Activation magnitude" and h["Fig 2d"]["A5"].startswith("Region-averaged")
         and h["Fig 4h"]["A3"].startswith("18-") and len(h["Fig 2e"])==3
         and all("dose" not in x.lower() and "specimen" not in x.lower() for x in h["Fig 2e"].values()),
         "18um small panel duplicated, component panel supplies no dose or specimen ID")
    values={(x["record"],x["observable"]):(x["locator"],q(x["canonical_decimal"]))
            for x in audit["source_crosschecks"]}
    need(len(values)==len(audit["source_crosschecks"])==85,"five saved scalar cells for each of 17 archive records")
    rows=list(csv.DictReader(raw(ROOT/CSV,pins[CSV]).decode("utf8").splitlines()))
    need(len(rows)==17,"exact archived 17-row section")
    mapping=[]
    means=[]
    for letter,row in zip(LETTERS,rows,strict=True):
        rid="optoGEF2_small_col_"+letter
        need(row["record_id"]==rid and row["identity_status"]=="Source-column record; no explicit biological embryo ID",
             "source column label must not be treated as a specimen ID")
        loc={"activation_magnitude":f"Fig 2c!{letter}4","maximum_shortening_um":f"Fig 2d!{letter}4",
             "region_averaged_shortening_um":f"Fig 2c!{letter}5","apical_max_shortening_um":f"Fig 2e!{letter}4",
             "basal_max_shortening_um":f"Fig 2e!{letter}5"}
        need(all(values[(rid,key)]==(loc[key],q(row[key])) for key in loc)
            and row["duplicate_maximum_source"]==f"Fig 4h!{letter}5","same exact archived small scalar locators")
        means.append(q(row["region_averaged_shortening_um"]))
        for channel in CHANNELS:
            mapping.append({"record_id":rid,"observable":channel,"source":loc[channel],
                "dose_association":("EXACT_PUBLISHED_SCALAR_JOIN" if channel in CHANNELS[:2] else
                                     "COMMON_B_R_ORDER_ASSUMPTION_PRESERVING_COMPONENT_PAIR")})
    need(len(set(means))==17,"sampled small source mean key is unique")
    return mapping

def matrices(protocol,pins,labels):
    src=read(ROOT/SMALL,pins[SMALL])
    smallp=read(ROOT/SP,pins[SP])
    oldp=read(ROOT/LP,pins[LP])
    old=read(ROOT/LARGE,pins[LARGE])
    review=read(ROOT/OLD_REVIEW,pins[OLD_REVIEW])
    need(src["protocol_sha256"]==pins[SP] and old["protocol_sha256"]==pins[LP] and
         old["input_hashes"][SMALL]==pins[SMALL] and review["status"]=="PASS_INDEPENDENT_ACTUAL_LOCAL_AFFINE_DIAGNOSTIC"
         and all(x is True for x in review["checks"].values()) and all(review["bindings"][name]==pins[name]
                                                            for name in (SMALL,SP,LARGE,LP)),
         "active original tangent source receipts match")
    archived_box=[tuple(map(binary,x)) for x in smallp["trust_box"]]
    need(archived_box==list(map(pair,oldp["trust_box_exact"]))==list(map(pair,protocol["trust_box_exact"]))
         and len(archived_box)==5 and archived_box[2]==(0,binary(.02))
         and all(ends==(-binary(.1),binary(.1)) for i,ends in enumerate(archived_box) if i!=2)
         and smallp["parameters"]["apical_surface_tension"]==0 and smallp["parameters"]["mu"]==1
         and smallp["parameters"]["bulk"]==5000 and Fraction(smallp["fixed_length_um"])==Fraction(67,2),
         "original exact rounding, one-sided Ta and fixed scale")
    need(tuple(protocol["parameter_columns"])==tuple(old["columns"])==tuple(oldp["columns"])
         ==tuple(smallp["active_LP_columns"])==COLS and tuple(protocol["observer_variants"])==OBS,
         "all five unchanged original parameter directions")
    need(q(oldp["shoulder_um_exact"])==Fraction("24.94698765")
         and old["heldout_response_tables_read"] is False,
         "original fixed radius, archive only model large objective")
    new_J={variant:[] for variant in OBS}
    recordpins={}
    order=[{"record_id":x["record_id"],"observable":x["observable"]} for x in labels]
    for n,(row,sp) in enumerate(zip(src["records"],smallp["records"],strict=True)):
        expect=f"optoGEF2_small_col_{LETTERS[n]}"
        need(row["record_id"]==sp["record_id"]==expect
             and review["bindings"][row["file"]]==row["sha256"],"same archive record in original independent audit")
        r=read(ROOT/row["file"],row["sha256"])
        recordpins[row["file"]]=row["sha256"]
        need(r["index"]==n and r["record_id"]==expect and r["input_field_sha256"]==sp["field_sha256"]
             and tuple(x["coordinate"] for x in r["coordinate_diagnostics"])==COLS
             and all(len(x["two_steps"])==2 for x in r["coordinate_diagnostics"]),
             "same saved state, tangent direction and selected original second step")
        for variant in OBS:
            for obs in CHANNELS:
                cols=[binary(midpoint(x["two_steps"][1]["central_difference_intervals"][variant][obs]))
                      for x in r["coordinate_diagnostics"]]
                new_J[variant].append(cols)
    need(len(recordpins)==17,"all17 original per-record archives")
    for v,row in zip(OBS,src["candidates"],strict=True):
        aff=row["affine_LP"]
        need(v==row["variant"] and row["all68_constraint_order"]==order
             and list(map(pair,aff["trust_box_exact"]))==archived_box
             and new_J[v]==[list(map(q,line)) for line in aff["affine_derivative_coefficients_exact"]]
             and any(q(x)!=0 for x in aff["affine_residual_intercepts_exact"]),
             "archived model rounded coefficients reproduced from individual second-scale observer FD intervals")
    objective=[binary(midpoint(p)) for p in old["contrast_gradient_intervals_exact"]]
    need(objective==list(map(q,old["rounded_affine_contrast_gradient_exact"])) and len(objective)==5
         and q(protocol["allowance_um_exact"])==EPS and EPS in list(map(q,oldp["small_prediction_change_budgets_um_exact"])),
         "original large archived signed objective and exact unchanged tolerance")
    old_range,=[pair(r["contrast_range_width_um_exact"]) for r in old["rows"]
                if r["mode"]=="small_prediction_change" and q(r["budget_exact"])==q(r["allowance_um_exact"])==EPS]
    keep=[i for i,x in enumerate(labels) if x["dose_association"]=="EXACT_PUBLISHED_SCALAR_JOIN"]
    need(len(keep)==34 and len(order)==68 and all(
         new_J[OBS[k]][i]==new_J[OBS[k+1]][i] for k in (0,2) for i in keep),
         "same within-weight two total observer channels")
    return {"order":order,"subset":{"source_joined34":keep,"ordered68":list(range(68))},
            "J":new_J,"box":archived_box,"objective":objective,"old_range":old_range,"record_sha256":recordpins}


def certificate(c,J,B,w,eps):
    step=list(map(q,w["step_exact"]))
    dual=list(map(q,w["multipliers_exact"]))
    need(len(step)==len(c)==len(B)==5 and len(dual)==2*len(J),"all exact LP directions and both multipliers per row")
    need(w["exact_primal_feasible"] is True and w["exact_weak_duality"] is True
         and 0<=q(w["repair_alpha_exact"])<=1,"witness flags and recorded feasibility repair")
    need(all(lo<=v<=hi for v,(lo,hi) in zip(step,B,strict=True)),"primal point outside frozen box")
    signed=J+[[-x for x in row] for row in J]
    slack=[eps-dot(row,step) for row in signed]
    need(min(slack)>=0,"primal violates a centered signed source constraint")
    need(min(dual)>=0,"negative dual certificate invalid")
    remainder=[c[j]-sum(row[j]*l for row,l in zip(signed,dual,strict=True)) for j in range(5)]
    support=sum(x*(hi if x>=0 else lo) for x,(lo,hi) in zip(remainder,B,strict=True))
    upper=eps*sum(dual)+support
    lower=dot(c,step)
    need(0<=lower<=upper and (lower,upper,upper-lower)==tuple(q(w[k]) for k in
            ("primal_lower_exact","dual_upper_exact","duality_gap_exact")),"stored exact objective/weak duality mismatch")
    return (lower,upper),step,min(slack)


def need_reject(f):
    try:
        f()
    except ValueError:
        return True
    return False


def source_law(pins):
    parent=raw(ROOT/"scripts/axisymmetric_mixed.py",pins["scripts/axisymmetric_mixed.py"]).decode()
    child=raw(ROOT/"scripts/axisymmetric_quadratic.py",pins["scripts/axisymmetric_quadratic.py"]).decode()
    nonlinear=raw(ROOT/"scripts/transfer_ambiguity_nonlinear_probe.py",PINS[ROOT/"scripts/transfer_ambiguity_nonlinear_probe.py"]).decode()
    def func(text,name):
        return next(n for n in ast.walk(ast.parse(text)) if isinstance(n,ast.FunctionDef) and n.name==name)
    def assignment(node,s):
        for x in ast.walk(node):
            if isinstance(x,ast.Assign) and any(ast.unparse(t)==s for t in x.targets):
                return ast.unparse(x.value)
        return None
    bulk=func(parent,"energy_components")
    full=func(child,"energy_components")
    drive=func(nonlinear,"strength_receipt")
    terms={"shear":assignment(bulk,"shear"),"pointwise_log_volume":assignment(full,"full"),
           "inherited_active":assignment(bulk,"active"),"projected_term_replacement":assignment(full,"difference"),
           "pointwise_total":assignment(full,"result['total']"),"trial_penalty":assignment(bulk,"invalid"),
           "exact_loading_product":assignment(drive,"requested")}
    need("p.mu / 2" in terms["shear"] and "self.parameters.bulk / 2" in terms["pointwise_log_volume"]
         and "activation_strength / 2" in terms["inherited_active"]
         and terms["exact_loading_product"]=="dose * gain"
         and terms["projected_term_replacement"]=="full - result['projected_volumetric']"
         and terms["pointwise_total"]=="result['total'] + difference"
         and "torch.relu(1e-06 - J)" in terms["trial_penalty"],
         "same pinned inherited bulk, active and pointwise replacement")
    increments={ast.unparse(x.targets[0]):ast.unparse(x.value) for x in ast.walk(bulk) if isinstance(x,ast.Assign)}
    surface=increments["surface"]
    tension=increments["tension_energy"]
    need("weight * modulus / 2 * (Js - 1) ** 2" in surface and "weight * tension * Js" in tension,
         "same material area and tension linearity")
    terms["surface"]=surface
    terms["tension"]=tension
    return {"source_expressions":terms,
            "mechanical_claim":"Homogeneous in all seven force-scale coefficients on the smooth positive physical domain, fixed geometry/support/constraints, fixed inputs, no independent unscaled force and scale-invariant branch selection.",
            "literal_numerical_analytic_branch":"J > 1e-6 keeps the inherited clamp inactive and the unscaled inversion penalty zero. This is a physical law and exact algebra check, not a floating-run equivalence claim.",
            "mechanical_ratios":["kappa/mu","K_a/(mu L0)","K_b/(mu L0)","T_a/(mu L0)","T_b/(mu L0)","g/mu"],
            "optional_common_conversion_ratio":"g*beta/mu; g and beta reciprocally exchange only when the conversion is common.",
            "L0_um_exact":"67/2","tangent_fixed_kappa_over_mu_exact":"5000"}

def compare_exact(data,inp,result):
    import numpy as np
    need(inp["protocol_sha256"]==PINS[ASSESS/"protocol.json"]
         and result["status"]=="COMPLETE_ROUNDED_AFFINE_CHANNEL_COMPARISON"
         and result["source_sha256"]==PINS[ASSESS/"assess.py"]
         and result["input_snapshot_sha256"]==PINS[ASSESS/"inputs.json"]
         and result["protocol_sha256"]==PINS[ASSESS/"protocol.json"]
         and result["larger_observation_data_read"] is False and
         type(result["new_force_derivative_or_equilibrium_jobs"]) is int and result["new_force_derivative_or_equilibrium_jobs"]==0,
         "same frozen input/producer/result with no scientific jobs")
    ignore=inp["residual_offsets_exact"]
    need(len(ignore)==68 and all(type(x) is str and q(x)==0 for x in ignore)
         and inp["observation_order"]==data["order"] and inp["zero_based_row_indices"]==data["subset"]
         and q(inp["allowance_um_exact"])==EPS and tuple(inp["parameter_columns"])==COLS
         and list(map(pair,inp["trust_box_exact"]))==data["box"]
         and list(map(q,inp["large_contrast_c_exact"]))==data["objective"]
         and pair(inp["original_disk_spatial_max_width_um_exact"])==data["old_range"]
         and {v:[list(map(q,r)) for r in J] for v,J in inp["small_J_by_variant_exact"].items()}==data["J"],
         "same full 68, certified 34 subset, exact rounding and zero baseline offsets")
    expected=[(v,k) for v in OBS for k in ("source_joined34","ordered68")]
    need([(r["variant"],r["group"]) for r in result["rows"]]==expected,"all eight ordered comparisons")
    verified=[]
    widths={}
    intervals={}
    negative_case=None
    for row,(v,group) in zip(result["rows"],expected,strict=True):
        J=[data["J"][v][i] for i in data["subset"][group]]
        B=data["box"]
        need(len(J)==row["row_count"] and tuple(row["endpoints"])==NAMES,"correct exact signed comparison")
        exacts=[]
        endpoint_details={}
        for sign,name in zip((1,-1),NAMES,strict=True):
            c=[sign*x for x in data["objective"]]
            w=row["endpoints"][name]
            bounds,d,slack=certificate(c,J,B,w,EPS)
            need(all(abs(dot(data["J"][v][i],d))<=EPS for i in data["subset"]["source_joined34"]),
                 "full-model primal always remains feasible on the subset, dropping constraints only")
            exacts.append(bounds)
            intervals[v,group,name]=bounds
            endpoint_details[name]={"primal_and_dual_exact":list(map(str,bounds)),
                                    "minimum_signed_primal_slack_exact":str(slack),
                                    "gap_exact":str(bounds[1]-bounds[0])}
            if negative_case is None:
                negative_case=c,J,B,w
        rr=tuple(sum(values) for values in zip(*exacts,strict=True))
        need(rr==pair(row["width_um_exact"]) and list(map(float,rr))==row["width_um_display"]
             and rr[0]>0 and rr[1]-rr[0]<Fraction(2,10**15),"same exact two-endpoint width")
        rank=exact_rank(J)
        small_svd=np.linalg.svd(np.array([[float(q) for q in x] for x in J]),compute_uv=False).tolist()
        need(rank==row["numerical_rank_diagnostic"]==5 and len(small_svd)==5 and all(
            math.isclose(a,b,rel_tol=2e-12,abs_tol=1e-15)
            for a,b in zip(small_svd,row["singular_values_floating_diagnostic"],strict=True)),
            "rounded exact row rank and numerical SVD check, neither a physical identifiability proof")
        widths[v,group]=rr
        verified.append({"variant":v,"group":group,"rows":len(J),"rank_of_exact_stored_matrix":rank,
                         "singular_values_binary64_diagnostic":small_svd,"width_um_exact":list(map(str,rr)),
                         "endpoints":endpoint_details})
    growth={}
    for v in OBS:
        a,b=widths[v,"source_joined34"]
        c,d=widths[v,"ordered68"]
        rr=(a/d-1,b/c-1)
        expected_growth=result["widening"][v]
        need(a>d and rr==pair(expected_growth["fractional_width_increase_exact"])
             and list(map(float,[100*x for x in rr]))==expected_growth["percent_display"],
             "strict widening and safe interval ratio")
        # The exact inclusion F68 ⊆ F34 is independent of the numerical
        # certificate. An endpoint making the opposite strict claim is impossible.
        for name in NAMES:
            need(intervals[v,"ordered68",name][0]<=intervals[v,"source_joined34",name][1],
                 "endpoint bounds contradict mathematical feasible set inclusion")
        growth[v]={"ratio_increase_exact":list(map(str,rr)),"percent_interval_display":list(map(float,[100*x for x in rr]))}
    need(widths[OBS[0],"source_joined34"]==widths[OBS[1],"source_joined34"] and
         widths[OBS[2],"source_joined34"]==widths[OBS[3],"source_joined34"],"same dropped-channel widths")
    a,b=widths[OBS[0],"ordered68"]
    c,d=data["old_range"]
    need(a<=d and c<=b,"same original baseline 68-row numerical certificate agrees")
    c,J,B,w=negative_case
    bad_step={**w,"step_exact":["1",*w["step_exact"][1:]]}
    bad_dual={**w,"multipliers_exact":["-1",*w["multipliers_exact"][1:]]}
    altered_upper={**w,"dual_upper_exact":str(q(w["dual_upper_exact"])+EPS)}
    expanded_constraint=J+[[Fraction(1000000)]*5]
    expanded_dual={**w,"multipliers_exact":[*w["multipliers_exact"][:len(J)],"0",
                                           *w["multipliers_exact"][len(J):],"0"]}
    neg={"box_infeasible_step":need_reject(lambda:certificate(c,J,B,bad_step,EPS)),
         "negative_dual_weight":need_reject(lambda:certificate(c,J,B,bad_dual,EPS)),
         "optimistic_or_corrupt_saved_objective":need_reject(lambda:certificate(c,J,B,altered_upper,EPS)),
         "wrong_objective_sign":need_reject(lambda:certificate([-x for x in c],J,B,w,EPS)),
         "centered_primal_observation_violation":need_reject(lambda:certificate(c,expanded_constraint,B,expanded_dual,EPS))}
    need(all(neg.values()),"altered witness not rejected")
    return verified,growth,neg


def frozen_protocol():
    p,pins=input_protocol()
    return {"status":"FROZEN_INDEPENDENT_ASSOCIATION_AFFINE_READBACK","scope":SCOPE,"source_sha256":digest(SOURCE),
            "assessment_sha256":{rel(p):h for p,h in PINS.items()},"original_source_sha256":pins,
            "method":"Independent reconstruction from archived per-record intervals and original observable tangent intervals, exact rational endpoint inequalities, no optimizer.",
            "posthoc_reproducibility_audit_not_preregistration":True,"max_source_summary_joins":34,"no_actual_new_mechanics":True}


def prepare():
    need(not (HERE/"actual").exists(),"retain past audits, no repeat")
    p=frozen_protocol()
    need(exact_rank([[Fraction(i==j) for j in range(5)] for i in range(5)])==5
         and exact_rank([[Fraction(1)]*5]*4)==1 and need_reject(lambda:pair(("1","0"))),
         "cold exact rank and endpoint order fixtures")
    h=persist(HERE/"protocol.json",p)
    print({"source_sha256":digest(SOURCE),"protocol_sha256":h},flush=True)


def run(ph):
    frozen=read(HERE/"protocol.json",ph)
    need(frozen==frozen_protocol(),"changed independent exact protocol or source")
    out=HERE/"actual"
    out.mkdir(exist_ok=False)
    mark={"pid":os.getpid(),"source_sha256":digest(SOURCE),"protocol_sha256":ph}
    persist(out/"started.json",{"status":"RUNNING",**mark})
    terminal={"status":"FAILED","exit_code":1,**mark}
    try:
        p,pins=input_protocol()
        a=read(ASSESS/"inputs.json",PINS[ASSESS/"inputs.json"])
        b=read(ASSESS/"results.json",PINS[ASSESS/"results.json"])
        need(a["study_relative_sha256"]==pins,"same input snapshot source closure")
        labels=source_table(pins)
        rebuilt=matrices(p,pins,labels)
        rows,growth,neg=compare_exact(rebuilt,a,b)
        functions=source_law(pins)
        h=persist(out/"receipt.json",{"status":"PASS_INDEPENDENT_ASSOCIATION_SAFE_AFFINE_READBACK",**mark,
              "scope":SCOPE,"all_protocol_input_hashes_revalidated":True,"protocol_input_hash_count":len(pins),
              "original_observer_tangent_records":rebuilt["record_sha256"],"row_labels_and_dose_evidence":labels,
              "exact_ordered68_constraint_rows":68,"published_summary_dose_linked_rows":34,
              "component_rows_allowed_by_order_only":34,"verified_specimen_ID_joins":0,
              "all_16_endpoint_certificates_verified":True,"all_8_ranges_verified":True,
              "retrospective_reproducibility_not_preregistration":True,"rows":rows,
              "widening":growth,"negative_witnesses_rejected":neg,"source_law":functions,
              "physical_or_experimental_identifiability_proved":False,"raw_workbook_or_heldout_response_payloads_opened":0,
              "new_force_AD_FE_spatial_or_equilibrium_jobs":0,"scipy_optimization_or_assessor_import":False,
              "full_five_gate_complete":False})
        terminal.update(status="FINISHED",exit_code=0,receipt_sha256=h)
        print({"receipt_sha256":h,"endpoint_count":16,"width_count":8,**mark},flush=True)
    except BaseException as e:
        persist(out/"failure.json",{"error":str(e),"traceback":traceback.format_exc(),**mark})
        raise
    finally:
        persist(out/"terminal.json",terminal)


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("mode",choices=("prepare","run"))
    a.add_argument("--protocol-sha256")
    args=a.parse_args()
    prepare() if args.mode=="prepare" else run(args.protocol_sha256)
