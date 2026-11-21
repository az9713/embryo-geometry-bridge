"""Frozen prelaunch review of observer-derived certification budgets."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
PRE=HERE.parent
ROOT=HERE.parents[5]
PRIMARY=ROOT/"runs/transfer_ambiguity_nonlinear_probe_v1"
NATIVE=PRIMARY/"independent_actual_preparation/actual_replay/receipt.json"
MANIFEST=PRIMARY/"field_manifest.json"
PRIMARY_PROTOCOL=PRIMARY/"preparation/protocol.json"
SOURCE=ROOT/"scripts/transfer_ambiguity_error_budgets.py"
PINS={
 SOURCE:"410a4944798b1f42f0360cbdf236d9da68f8de2615d93f860bccdd8be54d8859",
 PRE/"protocol.json":"332b66e28a9c990b0e6c59953d36c841d5c047e8bb9a5b9d47467e97598a3d27",
 PRE/"preflight.json":"4588b3dd1779d00e85f80e8dec62d3862a9605c883f73ff57c1ada78374b0e3f",
 PRIMARY_PROTOCOL:"9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94",
 NATIVE:"c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd",
 MANIFEST:"308ca4c05e476ba6a14cbc7b8d34018983b8f11263e2c224488672b35eb38ac0",
}
L=Q(67,2)
OBS=("maximum_shortening_um","region_averaged_shortening_um","apical_max_shortening_um","basal_max_shortening_um")


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def write(path,data):
    with path.open("x",encoding="utf-8") as f:
        json.dump(data,f,indent=2)
        f.write("\n")


def true(checks):
    return bool(checks) and all(value is True for value in checks.values())


def synthetic_formula(m):
    # Distinct face slope bounds, exact axis end nodes and an interior ROI.
    R=(Q(0),Q(1,2),Q(1))
    a=[[(t,2*t+10) for t in R]]
    b=[[(t,-t) for t in R]]
    faces,t=m.geometry_targets({"top":a,"bottom":b},(Q(0),Q(3,4)))
    A=Q(5,4)*L*3
    B=Q(5,4)*L*2
    assert t==1/(16*L)
    assert Q(faces["top"]["C_um_per_dimensionless_coordinate_exact"])==A
    assert Q(faces["bottom"]["C_um_per_dimensionless_coordinate_exact"])==B
    # First two share total height, last two get one face; disk mean is
    # invariant under changing the origin of the reference height.
    xi=Q(1,2**54)
    margins=[Q(1,1000),Q(1,900),Q(1,1100),Q(1,800)]
    rows={k:{"margin_um_exact":str(e),"stored_worst_residual_um_exact":str(2-e)}
          for k,e in zip(OBS,margins,strict=True)}
    got=m.small_budget(rows,faces,xi,t,Q(2))
    expected_coeff=(A+B,A+B,A,B)
    e=min([t]+[marg/(2*c) for marg,c in zip(margins,expected_coeff,strict=True)])-xi
    assert Q(got["e_goal_dimensionless_exact"])==e>0
    check=m.validate_design_transfer({"top":a,"bottom":b},(Q(0),Q(3,4)),faces,got)
    assert check["status"]=="CONDITIONAL_EXACT_ERROR_TRANSFER"
    caps=[Q(got["conditional_observer_error_caps_um_exact"][k]) for k in OBS]
    assert all(q<=marg/2 for q,marg in zip(caps,margins,strict=True))
    assert tuple(map(Q,got["conditional_observer_coefficient_um_exact"].values()))==expected_coeff
    for g in (Q(1,100),Q(2)):
        b0=m.large_budget((Q(0),Q(1)),faces,xi,t,g)
        assert Q(b0["balanced_error_allocation_um_exact"])==g/4
        assert Q(b0["e_goal_dimensionless_exact"])==min(t,g/(4*(L+A)))-xi>0
        assert Q(b0["conditional_observer_error_caps_um_exact"]["Q"])<=g/4
    assert m.strict_goal({"Q":Q(0)},{"Q":A},xi,t)["e_goal_dimensionless_exact"] is None
    assert m.strict_goal({"Q":Q(3)},{"Q":A},t,t)["e_goal_dimensionless_exact"] is None
    x,gaps=m.addition_gaps([1.0,0.0,3.0],[1.0,0.0,3.0],[0],[2.**-54])
    assert x==[1+xi,0,3] and gaps==[xi,0,0]
    return {"trace_radial_cap_exact":str(t),"synthetic_shifted_small_goal_exact":str(e),
            "expected_coefficients_exact":list(map(str,expected_coeff))}


def main():
    assert all(sha(f)==h for f,h in PINS.items())
    p,pre,primary,native,manifest=[read(f) for f in (PRE/"protocol.json",PRE/"preflight.json",PRIMARY_PROTOCOL,NATIVE,MANIFEST)]
    assert pre["status"]=="PASS" and true(pre["checks"])
    assert pre["source_sha256"]==p["source_sha256"]==PINS[SOURCE] and pre["protocol_sha256"]==PINS[PRE/"protocol.json"]
    assert pre["primary_protocol_sha256"]==PINS[PRIMARY_PROTOCOL]
    assert pre["primary_actual_fields_opened"]==pre["actual_budgets_computed"]==0
    assert pre["new_force_evaluations"]==pre["new_equilibrium_solves"]==0
    assert p["root_error_bounds_proved"]==p["new_force_evaluations"]==p["new_equilibrium_solves"]==0
    assert not (PRE.parent/"actual").exists()
    assert sha(PRE/"transfer_ambiguity_error_budgets_frozen.py")==PINS[SOURCE]
    for file,h in (p["fixed_design_hashes"]|p["dependencies"]|p["terminal_authority_hashes"]).items():
        assert sha(ROOT/file)==h
    assert p["observer"]==primary["observable_variant"]=="current_disk__spatial_max"
    assert Q(p["length_unit_um_exact"])==L
    assert Q(p["sigma_exact"])==Q(p["radial_margin_spending_fraction_exact"])==Q(1,2)
    assert Q(p["frozen_conditional_allowance_um_exact"])==Q(primary["selection"]["conditional_source_allowance_um_exact"])
    assert len(p["frozen_44_checkpoint_plan"])==native["checkpoints"]==len(manifest["fields"])==44
    assert native["exact_equilibrium_position_bounded"] is native["source_compatible_witness"] is False
    assert native["actual_results_sha256"]==manifest["primary_results_sha256"]
    assert manifest["independent_primary_review_sha256"]==PINS[NATIVE]
    assert manifest["independent_replay_observed_os_exit_code"]==0
    for x,y in zip(p["frozen_44_checkpoint_plan"],manifest["fields"],strict=True):
        assert all(x[k]==y[k] for k in ("endpoint_id","checkpoint","role","step","field","verification"))
        for k in ("field","verification"):
            assert native["all_hashes"][y[k]]==y[k+"_sha256"]
    selected=[x for x in p["frozen_44_checkpoint_plan"] if x["role"] in ("small_record","large_final")]
    assert p["selected_36_fields"]==selected
    assert len(selected)==36 and sum(x["role"]=="small_record" for x in selected)==34
    # Audit exact source paths and synthetic formulas. Do not call admit or
    # calculate and do not open any field/verification payload.
    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_error_budgets as m
    def forbidden(*args,**kw):
        raise AssertionError("actual field/budget path disabled")
    with (patch.object(m,"admit",side_effect=forbidden),patch.object(m,"calculate",side_effect=forbidden),
          patch.object(m,"model_receipt",side_effect=forbidden)):
        assert m.protocol()==p
        fixture=synthetic_formula(m)
    assert not (PRE.parent/"actual").exists()
    binds={"source_sha256":PINS[SOURCE],"protocol_sha256":PINS[PRE/"protocol.json"],
           "preflight_sha256":PINS[PRE/"preflight.json"],"primary_protocol_sha256":PINS[PRIMARY_PROTOCOL],
           "primary_actual_review_sha256":PINS[NATIVE],"whole_44_field_manifest_sha256":PINS[MANIFEST]}
    checks={
      "all_frozen_sources_protocol_preflight_and_native_terminal_premises_match":True,
      "all_44_metadata_and_36_selected_nonzero_field_manifest_rows_match":True,
      "exact_binary_injection_gap_and_fixed_611_map_axis_rim_guards_read":True,
      "reviewed_same_radius_P2_slope_and_derivative_margin_composition_correct":True,
      "small_half_slack_and_large_quarter_gap_allocations_reproduced":True,
      "observer_and_radial_cap_min_then_subtract_rounding_gap_and_reject_nonpositive":True,
      "no_actual_field_payload_open_or_budget_calculation_or_force_or_solve":True,
    }
    support=HERE/"support.json"
    write(support,{"status":"PASS_OBSERVER_ERROR_BUDGET_MATH_AND_SOURCE_REVIEW","checks":checks,
                   "bindings":{f.relative_to(ROOT).as_posix():h for f,h in PINS.items()},
                   "fixture":fixture,"scope":"Mathematical and source review of conditional dimensionless target budgets only.",
                   "actual_field_payloads_opened":0,
                   "root_error_bound_or_root_proved":False,
                   "source_uncertainty_supplied":False,
                   "source_notes":"model_receipt reconstructs free map independently from reference axis/rim conditions, demands exact fixed graph radii, checks binary addition/gain rounding and freshly implies geometry-local constants. validate_design_transfer crosschecks transfer at half radial reserve. Geometry, positive thickness and equilibrium of future roots remain separate.",
                   "audit_source_sha256":sha(Path(__file__))})
    gate=PRE/"prelaunch_independent_review.json"
    write(gate,{"status":"PASS","checks":checks,**binds,
                "support_file":support.relative_to(ROOT).as_posix(),"support_sha256":sha(support),
                "scope":"Conditional observer budget launch only, no new force, solve or root proof."})
    print(sha(gate),sha(support))


if __name__=="__main__":
    main()

