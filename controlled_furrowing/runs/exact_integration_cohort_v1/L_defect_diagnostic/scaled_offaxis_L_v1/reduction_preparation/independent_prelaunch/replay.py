"""Audit frozen scaled-L reducer admission and exact arithmetic before live inputs."""
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
import sys
from contextlib import redirect_stdout
from unittest.mock import patch

OUT = Path(__file__).resolve().parent
JOB = OUT.parent
ROOT = JOB.parents[4]
GATE = JOB / "prelaunch_independent_review.json"
EXPECTED = {
    "source_sha256": "0ebaf0c1dfd65a3962841be51b714f3540b6a7f07e226dab130591bd4b7f70e9",
    "protocol_sha256": "d5dae93ac5c3fc6d37dd24eac251f5873f2534fb121ec8e230fd31ebc50b3f1d",
    "preflight_sha256": "c82febf0c55b1c19e208e7f6d503dd80969e4c3920f8c469d6f24ae51be41007",
}


def reject(action, text):
    try:
        action()
    except ValueError as error:
        assert text in str(error), str(error)
    else:
        raise AssertionError("expected a binding or numerical rejection")


def packet_fixture(R,pins):
    pack = lambda a,b:{"lower_exact":str(a),"upper_exact":str(b)}
    make_h = lambda pairs:[{"i":i,"j":j,"enclosure":pack(a,b)} for (i,j),(a,b) in zip(((0,0),(0,1),(1,1)),pairs)]
    shared = {"record_id":R.RECORD,"dose_exact":"3/2","parameter_cell_sha256":"c",
              "model_observation_contract_sha256":"m","point_cache_sha256":"q","predictor_sha256":"s",
              "precision_bits":640,"element":6,"chart":"off_axis","global_indices":[3,0],
              "free_local":[1,5],"active":"0","nodal_hulls_sha256":"y","auxiliary_element_rule_sha256":"aux"}
    p = {**shared,"unchanged_degree":12,"unchanged_rho_exact":"1/4","force_tail_budget_exact":"1/8",
         "H_entry_tail_budget_exact":"1/16","force_tail_max_exact":"1/16","H_entry_tail_max_exact":"1/32"}
    o = {**shared,"degree":12,"rho_exact":"1/4","force_integral_minus_cache":[pack(-4,2), pack(1,4)],
         "H_integral_minus_cache_upper":make_h([(-3,-1),(-1,3),(0,4)])}
    force = [pack(-1,2),pack(2,3)]
    H = make_h([(-2,-1),(-1,1),(2,4)])
    n = {**shared,**{k:v for k,v in pins.items() if k != "prelaunch_sha256"},
         "status":"SCALED_OFFAXIS_COMPONENT_UNREVIEWED","degree":12,"rho_exact":"1/4",
         "selector_proof_sha256":R.D._digest(p),"original_packet_sha256":"a"*64,
         "force_integral_minus_cache":[pack(-1,3),pack(2,3)],"force_intersection_with_original":force,
         "H_intersection_with_original":H,"H_integral_minus_cache_upper":make_h([(-2,0),(-2,1),(2,4)]),
         "force_two_mass_tail_upper_exact":["1/16","1/32"],
         "same_Y_cache_consistency":{"all64_convex_strain_lines_positive":True,
         "minimum_cache_convex_line_J_lower":pack("1/2","3/4"),"max_coefficientwise_delta_q_upper":pack(0,1)}}
    for row,value in zip(n["H_integral_minus_cache_upper"],["1/64","1/32","0"]):
        row["tail_upper_exact"] = value
    verify = lambda x:R._scaled_packet(o,x,p,pins,"a"*64,[3,0])
    assert verify(n) == (force,H)
    for key in ("record_id","global_indices","free_local","nodal_hulls_sha256",
                "auxiliary_element_rule_sha256","rho_exact","predictor_sha256"):
        reject(lambda:verify({**n,key:"wrong"}),"provenance")
    reject(lambda:verify({**n,"force_two_mass_tail_upper_exact":["1/16","1/4"]}),"unchanged budget")
    reject(lambda:verify({**n,"force_two_mass_tail_upper_exact":["1/32","1/32"]}),"match unchanged")
    changed = deepcopy(n)
    changed["same_Y_cache_consistency"]["minimum_cache_convex_line_J_lower"] = pack(0,1)
    reject(lambda:verify(changed),"analytic cache")
    changed = deepcopy(n)
    changed["H_intersection_with_original"][1]["enclosure"] = pack(0,1)
    reject(lambda:verify(changed),"Hessian intersection")
    plus = lambda seq_or_map:[max(abs(Q(y["lower_exact"])),abs(Q(y["upper_exact"]))) for y in seq_or_map]
    assert max(plus(force)) <= min(max(plus(o["force_integral_minus_cache"])),max(plus(n["force_integral_minus_cache"])))
    old_g,old_H,inter_g,inter_H = [R.D.ZERO]*4,{},[R.D.ZERO]*4,{}
    R.D._scatter(old_g,old_H,[3,0],o["force_integral_minus_cache"],o["H_integral_minus_cache_upper"])
    R.D._scatter(inter_g,inter_H,[3,0],force,H)
    # An unchanged incident element is signed-scattered into both versions.
    incident_g = [pack(-1,0),pack(-2,1)]
    incident_H = make_h([(0,1),(2,4),(-1,0)])
    for out_g,out_H in ((old_g,old_H),(inter_g,inter_H)):
        R.D._scatter(out_g,out_H,[0,2],incident_g,incident_H)
    matrix = [[1.0,-0.25,0.5,-0.125],[-1.0,2.0,1.0,0.],[0.,0.,2.,0.],[-0.5,0.25,-0.5,1.]]
    rg0,rh0,n0 = R.D._dyadic_reduce(matrix,old_g,R.D._h_rows(4,old_H))
    rg1,rh1,n1 = R.D._dyadic_reduce(matrix,inter_g,R.D._h_rows(4,inter_H))
    assert n0 == n1 and all(a[0] <= b[0] <= b[1] <= a[1] for a,b in zip(rg0,rg1))
    assert all(b <= a for a,b in zip(rh0,rh1))


def admission_fixture(R,scale,pins):
    names = scale["actual_expected_files"]
    rows = [{"element":i,"file":name,"sha256":sha256(name.encode()).hexdigest()}
            for i,name in zip(R.OFFAXIS,names)]
    manifest = {**pins,"record_id":R.RECORD,"status":"ALL_132_SCALED_OFFAXIS_PACKETS_UNREVIEWED",
                "unchanged_axis_elements":R.AXIS,"original_surface_packets_retained":24,
                "global_reduction_run":False,"files":rows}
    terminal_sha,manifest_sha,review_sha = "a"*64,"b"*64,"c"*64
    terminal = {"status":"FINISHED","pid":424242,"exit_code":0,"results_sha256":manifest_sha}
    touched,outputs = [],[]

    def bind(path,h,dict_,parse=True):
        assert path in (f"{R.STAGE}/actual/terminal.json",f"{R.STAGE}/actual/results.json")
        terminal_file = path.endswith("terminal.json")
        assert h == (terminal_sha if terminal_file else manifest_sha)
        touched.append(path)
        dict_[path] = h
        return terminal if terminal_file else manifest

    def write(path,value):
        assert path.name == "input_admission.json"
        outputs.append(value)
        return "d"*64

    def reviewer(expected,dict_):
        assert expected == review_sha
        dict_["integrity_gate"] = expected

    with patch.object(R,"_bound_preparation",return_value=({},scale,[],{}, {},pins)), \
            patch.object(R,"_matches_prelaunch",side_effect=reviewer),patch.object(R.D,"_bound",side_effect=bind), \
            patch.object(R.D,"_write",side_effect=write),redirect_stdout(StringIO()):
        R.admit(review_sha,terminal_sha,manifest_sha)
        assert len(touched) == 2 and len(outputs) == 1
        assert outputs[0]["expected_scaled_packets"] == rows and outputs[0]["local_scaled_payloads_opened"] is False
        terminal["exit_code"] = 1
        reject(lambda:R.admit(review_sha,terminal_sha,manifest_sha),"terminal")
        terminal["exit_code"] = 0
        terminal["results_sha256"] = "0"*64
        reject(lambda:R.admit(review_sha,terminal_sha,manifest_sha),"binding")
        terminal["results_sha256"] = manifest_sha
        manifest["files"] = rows[::-1]
        reject(lambda:R.admit(review_sha,terminal_sha,manifest_sha),"canonical")
        manifest["files"] = rows
        manifest["prelaunch_sha256"] = "0"*64
        reject(lambda:R.admit(review_sha,terminal_sha,manifest_sha),"reviewed preparation")
        assert len(outputs) == 1


def run():
    if not __debug__:
        raise RuntimeError("Review requires assertions")
    bindings = {}

    def raw(name,expected=None):
        target = ROOT/name
        assert target.resolve().is_relative_to(ROOT)
        assert "scaled_offaxis_L_v1/actual/" not in name, name
        data = target.read_bytes()
        value = sha256(data).hexdigest()
        assert expected is None or value == expected, name
        assert bindings.setdefault(name,value) == value, name
        return data

    prefix = JOB.relative_to(ROOT).as_posix()
    src = raw(f"{prefix}/reduce_scaled_L.py",EXPECTED["source_sha256"])
    assert src == raw(f"{prefix}/reduce_scaled_L_frozen.py",EXPECTED["source_sha256"])
    proto,pre = [json.loads(raw(f"{prefix}/{key}.json",EXPECTED[f"{key}_sha256"])) for key in ("protocol","preflight")]
    assert proto["source_sha256"] == pre["source_sha256"] == EXPECTED["source_sha256"]
    assert pre["protocol_sha256"] == EXPECTED["protocol_sha256"]
    assert pre["status"] == "PASS_SCALED_REDUCTION_PREPARATION_ONLY"
    assert len(pre["checks"]) == pre["check_count"] == 19 and all(v is True for v in pre["checks"].values())
    assert proto["actual_scaled_reduction_run"] is proto["actual_scaled_packet_payloads_opened"] is False
    for path,h in proto["input_hashes"].items():
        raw(path,h)
    sys.path.insert(0,str(JOB))
    import reduce_scaled_L as R
    original_open = Path.open

    def guarded_open(path,*args,**kwargs):
        assert not path.resolve().is_relative_to(ROOT/R.STAGE/"actual"), path
        return original_open(path,*args,**kwargs)

    with patch.object(Path,"open",guarded_open):
        inputs,scale,selectors,old_admit,old_review,pins = R._bound_preparation()
        assert inputs == proto["input_hashes"]
        axis = list(range(6))+list(range(72,78))
        offaxis = set(range(144))-set(axis)
        assert proto["offaxis_intersections_required"] == scale["offaxis_element_indices"] == sorted(offaxis)
        assert proto["original_axis_packets_retained"] == scale["original_axis_element_indices"] == axis
        assert [v["element"] for v in selectors] == sorted(offaxis)
        assert proto["original_surface_packets_retained"] == len(old_admit["expected_packets"]["surface"]) == 24
        assert proto["expected_force_entries"] == old_review["force_entries"] == 611
        assert proto["expected_upper_Hessian_entries"] == old_review["upper_Hessian_entries"] == 6676
        assert len(R.tests()) == 18
        packet_fixture(R,pins)
        admission_fixture(R,scale,pins)
        path = f"{R.ORIGINAL}/protocol.json"
        p = json.loads(raw(path,old_admit["source_and_manifest_hashes"][path]))
        c = R.D.C.load_context(R.RECORD)
        assert c.row == p["row"] and c.dose == Q(scale["dose_exact"])
        assert c.old_radius == Q(scale["original_cube_radius_exact"])
        prior = {}
        eta0,z0,r,provider = R.D._old_premises(c,p,prior)
        assert (eta0,z0,r) == tuple(Q(old_review[k]) for k in ("cached_eta_upper_exact","cached_Z_upper_exact","normal_radius_exact"))
        norm = Q(old_review["R_infinity_norm_exact"])
        test_z = z0+(1-z0)/4
        stable = R.D._minimum(c,prior,z0,test_z,r,norm,provider)
        assert stable["strict_local_minimum"] is True
        assert Q(stable["normalized_curvature_lower_exact"]) == (1-test_z)/norm > 0
        radius_fraction_room = 1-z0-eta0/r
        assert 0 < radius_fraction_room < Q(old_review["b_upper_exact"])
        assert not R.D._strict_test(eta0,z0,Q(0),radius_fraction_room,r)[3]
        assert R.D._strict_test(eta0,z0,Q(0),radius_fraction_room/2,r)[3]
        for path,h in c.input_hashes.items():
            raw(path,h)
        for path,h in prior.items():
            raw(path,h)

    own = Path(__file__).relative_to(ROOT).as_posix()
    raw(own)
    for path,h in tuple(bindings.items()):
        raw(path,h)
    report = {
        **EXPECTED,"status":"PASS",
        "scope":"Prelaunch audit of complete L scaled/off-axis replacement admission, signed interval reduction, and conditional original-cube test. Shared reviewed moment, derivative, cache and endpoint arithmetic premises. No live scaled manifest, terminal or payload opened and no actual reduction.",
        "checks":{
            "frozen_reducer_protocol_preflight_and_parent_preparation_bindings_match":True,
            "source_admits_external_terminal_and_complete_digest_manifest_before_scaled_reads":True,
            "mock_terminal_and_admission_gate_reject_wrong_exit_digest_context_or_order":True,
            "original_132_scaled_candidates_plus_12_axis_and_24_surface_pins":True,
            "synthetic_provenance_same_geometry_rule_phi_and_cache_signs_enforced":True,
            "signed_F_and_symmetric_H_intersection_narrows_exact_R_contractions":True,
            "strict_original_cube_and_conditional_same_R_cached_SPD_homotopy_checked":True,
            "live_actual_read_guard_held_for_entire_prelaunch_replay":True,
        },
        "record_id":R.RECORD,"admitted_scope":{"offaxis":132,"axis_retained":12,"surface_retained":24},
        "old_cube_available_normalized_margin_exact":str(radius_fraction_room),
        "actual_scaled_packet_payloads_opened":False,"actual_scaled_reduction_run":False,
        "review_source_file":own,"review_source_sha256":bindings[own],
        "input_hashes":dict(sorted(bindings.items())),
        "actual_followthrough":"Use terminal/manifest digests observed only after true process exit. Independently review exact same-Y cache equality, stored tails, and all aggregated endpoint/scalar conclusions after actual reduction.",
    }
    if "--freeze" in sys.argv:
        receipt = OUT/"replay.json"
        assert not receipt.exists() and not GATE.exists()
        with receipt.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(report,indent=2)+"\n")
        report["support_file"] = receipt.relative_to(ROOT).as_posix()
        report["support_sha256"] = sha256(receipt.read_bytes()).hexdigest()
        with GATE.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":report["status"],"checks":len(report["checks"]),
                      "old_cube_available_normalized_margin":float(radius_fraction_room)}),flush=True)


if __name__ == "__main__":
    run()
