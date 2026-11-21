"""Read the revised bridge authority and exercise cold fixtures only."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import sys
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[7]
sys.path.insert(0,str(ROOT/"scripts"))
import transfer_ambiguity_exact_bridge as J

SOURCE = "58ae92f849414f185c4e28bbb2ec3b9d05361391c4d6e3fa71f519162e30e5e7"
PROTOCOL = "4d38e999565b452a276777b66690d9b8ed4df8f3576cd8ed26d8757b48900581"
PREFLIGHT = "8d7acca719a3f9466abc828c2732c5cb825a4012d2a7c1c30445cd306ba91484"
OLD_HOLD = "054a387ebd1b5e739c08cf92eb155626e35d09a20673cc0ca54669fc2b2aabf7"

def reject(fn):
    try:
        fn()
    except ValueError:
        return True
    return False


def main():
    p = J.PREP
    assert J.sha(J.SOURCE) == J.sha(p/"transfer_ambiguity_exact_bridge_frozen.py") == SOURCE
    assert J.sha(J.OUT/"preparation/independent_prelaunch/initial_revision/receipt.json") == OLD_HOLD
    protocol = J.pinned(p/"protocol.json",PROTOCOL)
    pf = J.pinned(p/"preflight.json",PREFLIGHT)
    assert pf["status"] == "PASS" and len(pf["checks"]) == 13 and J.sweep.all_true(pf["checks"])
    assert pf["selector_plans_run"] == pf["new_spatial_or_cached_AD_jobs"] == 0
    assert not J.PLAN.exists() and not J.ACTUAL.exists()
    stop = AssertionError("actual selector, spatial moment or cache AD entered the cold audit")
    with patch.object(J,"build_bulk",side_effect=stop), patch.object(J,"build_surface",side_effect=stop), \
            patch.object(J,"evaluate",side_effect=stop), patch.object(J.C,"local_terms",side_effect=stop), \
            patch.object(J.C,"assemble_point",side_effect=stop):
        state = J.static_inputs()
        assert protocol == J.authority_pins(state)
        assert pf["checks"] == J.fixtures(state)
        accepted = J.derivative_inputs(J.DERIV_ADMISSION,J.DERIV_ADMISSION_SHA,state[1])

    native = J.pinned(J.NATIVE/"receipt.json",J.NATIVE_SHA)
    wrapper = J.pinned(J.DERIV_ADMISSION,J.DERIV_ADMISSION_SHA)
    assert native["reviewed_cases"] == wrapper["reviewed_cases"]
    for pid,case,n in zip(J.EX.PILOTS,state[1],native["reviewed_cases"],strict=True):
        assert n["pilot_id"] == case["pilot_id"] == pid
        proof = J.pinned(J.ROOT/n["native_case_file"],n["native_case_sha256"])
        h = accepted[pid]
        zi = Q(proof["reviewed_ZQ_upper_exact"])
        zp = Q(proof["producer_ZQ_upper_exact"])
        z0 = Q(case["center_Z0_upper_exact"])
        eta = Q(case["point_eta_upper_exact"])
        assert zi == z0+Q(proof["own_variation_upper_exact"])
        assert eta+zi*J.RADIUS == Q(proof["own_cached_selfmap_lhs_upper_exact"]) < J.RADIUS
        z = max(zi,zp)
        assert z == Q(h["ZQ_upper_exact"]) < 1
        q = h["planning"]
        assert q["strict_positive_headroom"] is True and q["quarter_reference_is_gate"] is False
        assert Q(q["headroom_exact"]) == 1-z
        b = (1-z)/2
        f = ((1-z)-b)*J.RADIUS/2
        assert b == Q(q["b_star_exact"]) and f == Q(q["preconditioned_signed_total_force_allowance_exact"])
        assert f+b*J.RADIUS < (1-z)*J.RADIUS
        rn = Q(case["R_inf_exact"])
        old = case["prior_optional_tail_policy"]
        assert Q(q["per_component_force_entry_tail_cap_exact"]) == min(f/(10_000*rn*6),Q(old["tail_force_budget_exact"]))
        assert Q(q["per_component_H_entry_tail_cap_exact"]) == min(b/(10_000*rn*72),Q(old["tail_H_entry_budget_exact"]))

    # Mutate reviewed evidence after the real byte/path check, proving semantic rejection.
    original = J.pinned
    def replay():
        return J.derivative_inputs(J.DERIV_ADMISSION,J.DERIV_ADMISSION_SHA,state[1])
    def bad(path,edit):
        def read(file,h):
            value = original(file,h)
            return edit(deepcopy(value)) if file == path else value
        with patch.object(J,"pinned",side_effect=read):
            return reject(replay)
    assert bad(J.NATIVE/"receipt.json",lambda d:{**d,"status":"PASS_INDEPENDENT_UNRELATED","reviewed_cases":[]})
    assert bad(J.NATIVE/"terminal.json",lambda d:{**d,"pid":d["pid"]+1})
    assert bad(J.NATIVE/"os_exit_observation.json",lambda d:{**d,"observed_exit_code":False})
    first = native["reviewed_cases"][0]
    assert bad(J.ROOT/first["native_case_file"],lambda d:{**d,"preconditioner_infinity_norm_exact":"1"})
    assert bad(J.ROOT/first["native_case_file"],lambda d:{**d,"reviewed_ZQ_upper_exact":"0"})
    assert reject(lambda: J.derivative_inputs(J.DERIV_ADMISSION,"f"*64,state[1]))
    assert reject(lambda: J.derivative_inputs(J.NATIVE/"receipt.json",J.DERIV_ADMISSION_SHA,state[1]))

    def fixture(slot,ix,force,upper):
        return {"kind":"bulk","slot":slot,"global_indices":ix,
                "point_force":{"scope":J.POINT,"F_I_minus_F_cache_intervals_exact":force},
                "whole_box_hessian":{"scope":J.BOX,"radius_exact":str(J.RADIUS),
                                      "H_I_minus_H_cache_upper_exact":upper}}
    a = fixture(0,[2,0],[(1,2),(-3,1)],[(0,0,4,5),(0,1,-2,3),(1,1,5,6)])
    b = fixture(1,[0,2],[(2,4),(-2,-1)],[(0,0,-4,1),(0,1,1,2),(1,1,-5,-4)])
    order = [("bulk",0),("bulk",1)]
    f,h = J.scatter([a,b],order,3)
    assert f == [(-1,5),(0,0),(-1,1)] and h == {(0,0):(1,7),(0,2):(-1,5),(2,2):(-1,1)}
    assert reject(lambda: J.scatter([a,{**b,"slot":True}],order,3))
    assert reject(lambda: J.scatter([b,a],order,3))
    assert reject(lambda: J.scatter([a,a],order,3))
    src = J.SOURCE.read_text(encoding="utf8")
    assert "('integral','node_power','site_power','defect')" in src
    assert src.count("error.component_scope=") == 2
    assert not J.PLAN.exists() and not J.ACTUAL.exists()

    gate = {
        "status":"PASS","source_sha256":SOURCE,"protocol_sha256":PROTOCOL,"preflight_sha256":PREFLIGHT,
        "geometry_review_sha256":J.DX.PINS[J.DX.GEOMETRY_REVIEW],
        "cached_root_review_sha256":J.EX.PINS[J.EX.CACHE_REVIEW],
        "audit_source_sha256":J.sha(Path(__file__)),"initial_hold_receipt_sha256":OLD_HOLD,
        "checks":{
            "current_source_protocol_preflight_and_all_cold_input_pins":True,
            "native_review_and_producer_then_each_case_OS_authority_and_exact_z_links":True,
            "reject_mixed_native_producer_case_preconditioner_and_nonliteral_exit":True,
            "whole_radius_exact_cache_domain_and_optional_headroom_recomputed":True,
            "fresh_Arb640_factory_and_distinct_point_and_whole_box_route_scope":True,
            "explicit_Qstar_and_same_Y_cache_correction_available_across_routes":True,
            "independent_signed_scatter_permuted_columns_and_complete_upper_H_fixture":True,
            "point_moment_cache_cleared_and_structured_order_errors_covered":True,
            "independent_actual_plan_and_OS_receipts_are_required_before_workers":True,
            "cold_only_no_actual_selectors_integrals_or_AD":True},
        "scope":[
            "Cold source and real bound-metadata review only. No selector, point-force, Hessian or cache matrix was recomputed.",
            "Approval permits the plan phase. No worker launch before separate independent actual per-component plan review and the observed plan process exit.",
            "A future reducer must freshly admit whole 144+24 part manifests, worker and supervisor true OS exits, exact packet declarations, and independently reviewed packet endpoints before signed assembly."
        ],
        "pinned_inputs":len(protocol["input_sha256"]),
        "optional_planning_bounds":{pid:accepted[pid]["ZQ_upper_exact"] for pid in J.EX.PILOTS},
    }
    path = p/"prelaunch_independent_review.json"
    h = J.save(path,gate)
    state2,links = J.plan_gate(path,h)
    assert J.authority_pins(state2) == protocol and links["source_prelaunch_review_sha256"] == h
    print(json.dumps({"gate_sha256":h,"pin_count":gate["pinned_inputs"],"source":SOURCE}))


if __name__ == "__main__":
    main()
