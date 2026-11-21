"""Independent inert generic-wiring fixtures for later spatial observer transfer."""
from fractions import Fraction as Q
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent))
import apply as A


def rejected(fn):
    try:
        fn()
    except (KeyError,ValueError):
        return True
    return False


def scalar(M):
    r=A.RADIUS
    s,xi,e,z=r/16,r/64,r/8,Q(1,2)
    eta=e*(1-z)
    leaf={"radius_exact":str(r),"status":A.ROOT_PASS,"root_asserted":True,
          "target":{"archived_center_shift_exact":str(s),"original_addition_gap_exact":str(xi),
                    "observer_goal_exact":str(r)},
          "arithmetic":{"eta_upper_exact":str(eta),"z_upper_exact":str(z),"RF_lower_exact":"0",
                        "strict_selfmap":True,"strict_cube_exclusion":False,
                        "selfmap_lhs_upper_exact":str(eta+z*r),"posterior_error_upper_exact":str(e)},
          "observer_and_optional_minimum":{"archived_state_budget_proved":True,
                       "archived_shift_plus_posterior":str(s+e),
                       "physical_nodal_error_um_upper_exact":str(M.LENGTH*(s+e+xi))}}
    row={"status":A.ROOT_PASS,"exact_spatial_root_certified":True,"shift_inclusive_observer_budget_proved":True}
    assert A.locate(M,row,leaf)==(e,s,xi,M.LENGTH*(s+e+xi))
    excluded={**leaf,"status":A.EXCLUDE,"root_asserted":False,
        "arithmetic":{**leaf["arithmetic"],"eta_upper_exact":str(2*r),"RF_lower_exact":str(2*r),
                      "strict_selfmap":False,"strict_cube_exclusion":True,
                      "selfmap_lhs_upper_exact":str(Q(5,2)*r),"posterior_error_upper_exact":None}}
    assert A.locate(M,{"status":A.EXCLUDE,"exact_spatial_root_certified":False,
                      "shift_inclusive_observer_budget_proved":False},excluded) is None
    assert rejected(lambda:A.locate(M,{**row,"exact_spatial_root_certified":False},leaf))
    assert rejected(lambda:A.locate(M,row,{**leaf,"root_asserted":1}))
    assert rejected(lambda:A.locate(M,row,{**leaf,"target":{**leaf["target"],"original_addition_gap_exact":"0"}}))
    assert M.compare((Q(0),Q(1)),Q(3),Q(1,2),Q(2))[2]=="UNRESOLVED"
    assert M.compare((Q(0),Q(1)),Q(3),Q(1,2),Q(1))[2]=="OUT"


def wire(M,wave,profile):
    with TemporaryDirectory(dir=A.HERE) as scratch:
        work=Path(scratch)
        root=work/"global_reduction_v1/independent"
        root.mkdir(parents=True)
        pins={}
        def put(path,obj):
            h=A.save(path,obj)
            pins[A.rel(path)]=h
            return h
        identity={**A.GLOBAL_PHASE,"mode":"independent","wave_id":wave,"pid":91,
              **{k:str(i)*64 for i,k in enumerate(("hull_admission_sha256","raw_wave_admission_sha256",
                          "raw_root_observed_exit_sha256","producer_observed_os_sha256"),1)}}
        native=[]
        for c in profile["case_order"]:
            p=root/c
            p.mkdir()
            op=put(p/"operands.json",{"inert_fixture":True})
            rh=put(p/"results.json",{**identity,"case_id":c,"status":A.ROOT_PASS,"operands_sha256":op,
                   "target":profile["raw_worker_targets"][c],"radius_exact":str(A.RADIUS),
                   "full_five_gate_complete":False,"archived_observer_field_replaced":False})
            native.append({"case_id":c,"status":A.ROOT_PASS,"operands_sha256":op,"case_results_sha256":rh})
        started=put(root/"started.json",{"status":"RUNNING",**identity})
        rh=put(root/"results.json",{"status":A.G.INDEP_STATUS,**identity,"cases":native,"signed_hull_components":672,
                                  "new_force_AD_or_spatial_integrals":0,"full_five_gate_complete":False})
        terminal=put(root/"terminal.json",{"status":"FINISHED",**identity,"exit_code":0,"results_sha256":rh})
        outer={"status":A.OBSERVED,"mode":"independent","wave_id":wave,"pid":91,"observed_exit_code":0,
               "unified_exec_chunk_id":"a31cee","started_sha256":started,"terminal_sha256":terminal,"results_sha256":rh,
               "sha256":pins}
        path=root/"root_observed_exit.json"
        def attempt(row):
            if path.is_file():
                path.unlink()
            h=A.save(path,row)
            return A.native_gate(M,wave,profile,h,{})
        with patch.object(A.X,"wave_folder",return_value=work), patch.object(A.G,"producer_observed",return_value={}) as prior:
            assert tuple(attempt(outer))==tuple(profile["case_order"]) and prior.call_count==1
            assert rejected(lambda:attempt({**outer,"mode":"producer"}))
            assert rejected(lambda:attempt({**outer,"observed_exit_code":False}))
            assert rejected(lambda:attempt({**outer,"pid":91.0}))
            assert rejected(lambda:attempt({**outer,"sha256":{}}))
            (root/profile["case_order"][-1]/"operands.json").write_text('{"mutated":true}',encoding="utf-8")
            assert rejected(lambda:attempt(outer))


def main():
    m,profiles,_,_,_=A.source_inputs()
    scalar(m)
    wave=A.X.WAVES[0]
    wire(m,wave,profiles[wave])
    print(A.save(A.HERE/"synthetic_current_schema_readback.json",{"status":"PASS","source_sha256":A.sha(A.HERE/"apply.py"),
        "fixture_source_sha256":A.sha(Path(__file__)),"checks":{
        "generic_independent_true_exit_and_eleven_native_files":True,"changed_last_operand_float_pid_and_producer_mode_rejected":True,
        "fixed_sigma_e_xi_exact_bound_and_mock_cube_exclusion":True,
        "conditional_residual_boundary_signed_cases_checked":True},"opened_actual_later_global_or_observer_outputs":0,
        "actual_new_scientific_jobs":0}))


if __name__=="__main__":
    main()
