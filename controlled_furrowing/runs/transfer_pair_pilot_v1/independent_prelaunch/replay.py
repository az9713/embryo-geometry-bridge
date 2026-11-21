"""Independent transfer-pair prelaunch replay. Nonlinear entry points are blocked."""
from contextlib import redirect_stdout
from fractions import Fraction as Q
from hashlib import sha256
from io import StringIO
import csv
import json
from math import inf, nextafter
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

OUT = Path(__file__).resolve().parent
JOB = OUT.parent
ROOT = JOB.parents[1]
PREP = JOB / "preparation"
DRIVER = "scripts/transfer_pair_pilot.py"
GATE = JOB / "prelaunch_independent_review.json"
EXPECTED = {
    "driver_sha256": "a05d677ac7d919c2f012b20040c663e6f181563bc008fd4c216a81e6a30b7d18",
    "protocol_sha256": "a2745a7e50e9a64ee6e5e597f04b30981c7910ad5e3ca2d3a16e6d4d0f3de927",
    "preflight_sha256": "c6defad69457de5d28cffb1c1a657ddf36b71575b75a7078a2dff7df875a816d",
}
CASE_IDS = {"baseline_small", "baseline_large", "joint_small", "joint_large"}
LENGTH = Q(67, 2)


def rejects(action):
    try:
        action()
    except ValueError:
        return
    raise AssertionError("invalid bound fixture accepted")


def verify_rounding(exact, candidate):
    value = float(candidate)
    assert Q.from_float(value) == candidate
    lower = Q.from_float(nextafter(value, -inf))
    upper = Q.from_float(nextafter(value, inf))
    assert (lower+candidate)/2 <= exact <= (candidate+upper)/2
    assert candidate == Q.from_float(float(exact))


def independent_pair(p, docs):
    keys = {"record_id", "reported_diameter_um", "activation_magnitude"}
    small = [{k: r[k] for k in keys} for r in csv.DictReader(StringIO(docs["data/derived/calibration_small_18um.csv"]))]
    large = json.loads(docs["data/derived/spatial_inputs_only.json"], parse_float=str)
    assert len(small) == 17 and len(large) == 31
    assert all(set(r) == keys for r in large)
    assert len({r["record_id"] for r in large}) == 31 and len({r["record_id"] for r in small}) == 17
    median = sorted(Q(r["activation_magnitude"]) for r in small)[8]
    ranked = sorted(large, key=lambda r: (-Q(r["reported_diameter_um"]), abs(Q(r["activation_magnitude"]) - median), r["record_id"]))
    chosen = ranked[0]
    assert chosen["record_id"] == "fig4f_col_I"
    dose = Q(chosen["activation_magnitude"])
    small_d, = {Q(r["reported_diameter_um"]) for r in small}
    large_d = Q(chosen["reported_diameter_um"])
    assert small_d == Q(935512037, 50000000)
    assert large_d == Q(498939753, 10000000)
    assert dose == Q(44883687, 50000000)
    sel = p["input_selection"]
    assert sel["large_record_id"] == chosen["record_id"]
    assert tuple(Q(sel[k]) for k in ("large_reported_diameter_um_exact", "matched_small_reported_diameter_um_exact",
                                     "matched_dose_exact", "small_training_median_dose_exact")) == (large_d, small_d, dose, median)
    assert sel["heldout_response_tables_read"] is False
    obs = p["observable"]
    assert Q(obs["shoulder_radius_um_exact"]) == large_d/2 == Q("24.94698765")
    assert Q(obs["length_unit_um_exact"]) == LENGTH
    assert obs["units"] == "um" and obs["formula"] == "z_apical(current R=0) - z_apical(current R=R_shoulder)"

    basep = json.loads(docs["runs/joint_gain_extension_v2/protocol.json"])
    base = json.loads(docs["runs/joint_gain_extension_v2/direct_results.json"])
    jointp = json.loads(docs["runs/material_cohort_direct_v1/protocol.json"])
    joint = json.loads(docs["runs/material_cohort_direct_v1/results.json"])
    basefit, = [item for item in base["variants"] if item["variant"] == "current_disk__spatial_max"]
    jointfit, = [item for item in joint["variants"] if item["variant"] == "current_disk__spatial_max"]
    for fit in (basefit, jointfit):
        assert len(fit["all68_residuals"]["constraint_rows"]) == 68 and fit["all68_residuals"]["records"] == 17
    params = {"baseline": basep["physical_parameters"], "joint": jointp["new_parameters"]}
    gains = {"baseline": Q(basefit["selected_exact_SURROGATE_gain"]), "joint": Q(jointp["new_gain_exact"])}
    assert gains["baseline"] != gains["joint"]
    assert all(gain > 0 for gain in gains.values())
    assert p["numerical_policy"]["solve_quadrature"]["volume_duffy_gauss_order"] == 8
    assert p["numerical_policy"]["cross_force_quadrature"]["volume_duffy_gauss_order"] == 12
    assert p["numerical_policy"]["solve_quadrature"]["surface_gauss_order"] == 12
    assert p["numerical_policy"]["cross_force_quadrature"]["surface_gauss_order"] == 16
    assert {case["case_id"] for case in p["cases"]} == CASE_IDS and len(p["cases"]) == 4
    planned = [j*dose/3 for j in range(4)]
    final = {}
    for case in p["cases"]:
        material, size = case["case_id"].split("_")
        assert (case["candidate"], case["diameter_level"]) == (material,size)
        d = small_d if size == "small" else large_d
        assert Q(case["diameter_um_exact"]) == d
        assert Q(case["patch_radius_model_exact"]) == d/(2*LENGTH)
        actual_r = Q.from_float(case["parameters"]["patch_radius"])
        assert actual_r == Q(case["patch_radius_model_binary64_exact"])
        verify_rounding(d/(2*LENGTH), actual_r)
        expected = dict(params[material], patch_radius=float(d/(2*LENGTH)))
        assert case["parameters"] == expected
        assert Q(case["gain_exact"]) == gains[material]
        assert [Q(step["dose_exact"]) for step in case["schedule"]] == planned
        for d0, step in zip(planned, case["schedule"]):
            request = d0*gains[material]
            binary = Q(step["solver_strength_binary64_exact"])
            assert Q(step["requested_strength_exact"]) == request
            verify_rounding(request, binary)
            assert Q(step["binary64_rounding_difference_exact"]) == binary-request
        final[case["case_id"]] = case["schedule"][-1]["solver_strength_binary64_exact"]
    assert final["baseline_small"] == final["baseline_large"]
    assert final["joint_small"] == final["joint_large"]
    return {"chosen_large_record": chosen["record_id"], "shared_dose_schedule_exact": list(map(str,planned)),
            "shoulder_um_exact": str(large_d/2), "final_binary64_strength_exact": final}


def independent_observer(T):
    # Curved radial trace: at t=1/2, current r=5/8 and z=9/10.
    # The material radius at the shoulder is therefore 1/2, not 5/8.
    ref = [(0,0),(1,0),(1,1),(0,1),(Q(1,2),0),(1,Q(1,2)),
           (Q(1,2),Q(1,2)),(Q(1,2),1),(0,Q(1,2))]
    cur = [(r+r*(1-r)/2, Q(4,5)+Q(r,5)) if z == 1 else (r,z) for r,z in ref]
    tri = [[0,1,2,4,5,6],[0,2,3,6,7,8]]
    state = {"parameters":{"outer_radius":1,"height":1,"reference_curvature":0},
             "reference_nodes":ref,"current_nodes":cur,"triangles":tri}
    receipt = T.surface_receipt(state, Q(5,8)*LENGTH, full_profile=True)
    assert receipt["usable_stored_graph"] is True and receipt["height_units"] == receipt["radius_units"] == "um"
    q = list(map(Q, receipt["Q_um_exact"]))
    assert q[0] <= -LENGTH/10 <= q[1] < 0 and q[1]-q[0] <= Q(1, 2**60)
    assert Q(receipt["minimum_current_thickness_um_lower_exact"]) > 0
    profile = receipt["sampled_full_current_profiles"]
    assert Q(profile[0]["radius_um_exact"]) == 0 and Q(profile[-1]["radius_um_exact"]) == LENGTH
    assert not T.surface_receipt(state, LENGTH)["usable_stored_graph"]


def zero_replay(T, p, pre):
    out = {}
    fail_nonlinear = AssertionError("prelaunch entered nonlinear solve or continuation")
    with patch.object(T.AxisymmetricModel, "advance", side_effect=fail_nonlinear), \
            patch.object(T.AxisymmetricModel, "solve", side_effect=fail_nonlinear), \
            patch.object(T.CrossModel, "solve", side_effect=fail_nonlinear):
        for case in p["cases"]:
            model, cross = T.build_models(T.Parameters(**case["parameters"]))
            assert np.array_equal(model.reference, cross.reference)
            assert np.array_equal(model.triangles, cross.triangles) and np.array_equal(model.free, cross.free)
            assert np.array_equal(model.activation_elements, cross.activation_elements)
            r = model.reference[:,0]
            assert np.array_equal(model.fixed[:,0], (r == 0) | (r == 3))
            assert np.array_equal(model.fixed[:,1], r == 3)
            assert model.reference.shape == (325,2) and model.triangles.shape == (144,6)
            assert len(model.free) == 611 and int(sum(model.activation_elements)) == 24
            assert case["parameters"]["patch_radius"] in model.radial
            assert all(case["parameters"][k] in model.layers for k in ("activation_z_min", "activation_z_max"))
            assert model.integration_policy.volume_duffy_gauss_order == 8 and cross.integration_policy.volume_duffy_gauss_order == 12
            assert model.integration_policy.surface_gauss_order == 12 and cross.integration_policy.surface_gauss_order == 16
            zero = [0.] * len(model.free)
            assert np.array_equal(model.coordinates(T.torch.tensor(zero, dtype=T.torch.float64)).detach().numpy(),
                                  model.reference)
            r8 = float(np.max(np.abs(model.objective(zero, 0.)[1])))
            r12 = float(np.max(np.abs(cross.objective(zero, 0.)[1])))
            init = pre["actual_initial_meshes"][case["case_id"]]
            assert r8 == init["Q8_residual"] and r12 == init["Q12_residual"]
            assert max(r8, r12) < 1e-7 == p["numerical_policy"]["solve"]["force_tolerance"]
            state = T.identity_state(model)
            geom = T.exact_geometry(state)
            assert geom["global_embedding_certified"] is True and geom["all_exact_J_positive"] is True
            assert 0 < Q(geom["J_lower_exact"]) <= 1 <= Q(geom["J_upper_exact"])
            q = T.surface_receipt(state, Q(p["observable"]["shoulder_radius_um_exact"]))
            assert q["usable_stored_graph"] is True and q["Q_um_exact"] == ["0","0"]
            out[case["case_id"]] = {"Q8_zero_residual":r8,"Q12_zero_residual":r12}
    return out


def supervisor_fixture(T):
    with TemporaryDirectory(dir=OUT) as temp, patch.object(T, "OUT", Path(temp)), \
            patch.object(T, "relative", side_effect=lambda p: p.as_posix()), redirect_stdout(StringIO()):
        rel = Path(temp)
        case = {"case_id": "baseline_small"}
        path = rel/"jobs"/case["case_id"]
        calls = []

        class Child:
            pid = 999

            def wait(self):
                running = json.loads((path/"job.json").read_bytes())
                assert running["started"] is True and running["terminal"] is False and running["exit_code"] is None
                assert running["pid"] == self.pid
                # A completion-looking worker file cannot override an OS failure.
                T.save(path/"worker.json", {"status": "COMPLETE_CONDITIONAL_STORED_FIELD"})
                calls.append("wait")
                return 2

        with patch.object(T.subprocess, "Popen", return_value=Child()) as launch:
            decision = T.launch_case(rel/"synthetic_gate.json", case, dict(EXPECTED))
            assert launch.call_count == 1 and calls == ["wait"]
            assert launch.call_args.args[0][-2:] == ["--case", "baseline_small"]
        assert decision["status"] == "PROCESS_EXITED" and decision["terminal"] is True and decision["exit_code"] == 2
        assert decision["worker_status"] == "COMPLETE_CONDITIONAL_STORED_FIELD"
        assert T.summarize([decision])["status"] == "INCOMPLETE_PILOT"
        with patch.object(T.subprocess,"Popen",side_effect=OSError("synthetic launch failure")):
            failed = T.launch_case(rel/"synthetic_gate.json", {"case_id":"baseline_large"}, dict(EXPECTED))
        assert failed["status"] == "LAUNCH_FAILED" and failed["terminal"] is True
        assert failed["started"] is False and failed["exit_code"] is None
        assert T.summarize([failed])["status"] == "INCOMPLETE_PILOT"


def run():
    if not __debug__:
        raise RuntimeError("Review requires enabled assertions")
    inputs = {}

    def raw(path, h=None):
        file = ROOT/path
        data = file.read_bytes()
        observed = sha256(data).hexdigest()
        assert h is None or observed == h, path
        assert inputs.setdefault(path, observed) == observed, path
        return data

    prep = PREP.relative_to(ROOT).as_posix()
    raw(DRIVER, EXPECTED["driver_sha256"])
    p = json.loads(raw(f"{prep}/protocol.json",EXPECTED["protocol_sha256"]))
    pre = json.loads(raw(f"{prep}/preflight.json",EXPECTED["preflight_sha256"]))
    assert pre["status"] == "PASS" and len(pre["checks"]) == 40
    assert all(v is True for v in pre["checks"].values())
    assert pre["nonlinear_jobs_launched"] == 0
    assert pre["heldout_response_tables_read"] is False
    for key in ("driver_sha256", "protocol_sha256"):
        assert pre[key] == EXPECTED[key]
    docs = {file:raw(file,h).decode("utf8") for file,h in p["input_hashes"].items()}
    for file,h in p["source_hashes"].items():
        raw(file,h)
    pair = independent_pair(p, docs)
    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_pair_pilot as T
    T.torch.set_num_threads(1)

    allowed = {(ROOT/path).resolve() for path in docs | p["source_hashes"]}
    underlying = Path.open
    seen = set()

    def checked_open(path, *args, **kwargs):
        assert path.resolve() in allowed, f"unfrozen input read: {path}"
        seen.add(path.resolve())
        return underlying(path, *args, **kwargs)

    with patch.object(Path,"open",checked_open):
        assert T.create_protocol() == p
    assert seen == allowed
    for key in EXPECTED:
        bad = {**EXPECTED, key: "0"*64, "status": "PASS", "checks":{"fixture":True}}
        rejects(lambda: T.verify_review(bad,EXPECTED))
    rejects(lambda: T.verify_review({**EXPECTED,"status":"PASS","checks":{"fixture":1}},EXPECTED))
    independent_observer(T)
    zeros = zero_replay(T,p,pre)
    supervisor_fixture(T)
    script = Path(__file__).relative_to(ROOT).as_posix()
    raw(script)
    for path, h in tuple(inputs.items()):
        raw(path,h)
    report = {**EXPECTED,"status":"PASS",
              "scope":"Independent prelaunch of four conditional stored-field calculations only. No nonlinear jobs or larger response score. Exact observer intervals apply to serialized P2 curves.",
              "checks":{
                  "all_frozen_source_protocol_preflight_and_input_hashes_current":True,
                  "input_only_largest_then_nearest_median_selection_replayed":True,
                  "archived_material_gain_provenance_and_four_unmodified_cases":True,
                  "exact_diameter_radius_dose_gain_and_binary64_rounding":True,
                  "only_frozen_paths_opened_reproducing_protocol_no_larger_responses":True,
                  "negative_Q_uses_current_radius_and_physical_micrometres":True,
                  "same_P2_space_node_BC_and_Q8_Q12_zero_residuals_with_positive_stored_geometry":True,
                  "true_OS_wait_required_and_nonzero_exit_blocks_comparison":True,
                  "mixed_or_nonliteral_prelaunch_hash_binding_rejected":True,
              },"selection_and_schedule":pair,"identity_replay":zeros,
              "input_hashes":dict(sorted(inputs.items())),
              "nonlinear_jobs_launched":0,"heldout_response_tables_read":False,
              "postselection_used":False}
    if "--freeze" in sys.argv:
        receipt = OUT/"replay.json"
        assert not receipt.exists() and not GATE.exists()
        with receipt.open("x", encoding="utf8",newline="\n") as file:
            file.write(json.dumps(report,indent=2)+"\n")
        report["review_receipt_file"] = receipt.relative_to(ROOT).as_posix()
        report["review_receipt_sha256"] = sha256(receipt.read_bytes()).hexdigest()
        with GATE.open("x",encoding="utf8",newline="\n") as file:
            file.write(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status": report["status"],"checks": len(report["checks"]),
                      "selection_and_schedule":pair}),flush=True)


if __name__ == "__main__":
    run()
