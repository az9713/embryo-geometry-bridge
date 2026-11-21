"""Metadata-only census and exact radius advice. This is no launch admission."""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
EX = HERE.parent
ROOT = EX.parents[3]
DESIGN = EX.parent
INV = DESIGN / "full_outcome_coverage_v1"
W30 = INV / "remaining_cached_30_v1"
M = INV / "remaining_cached_wave1"
R36, R34 = Q(1, 2**36), Q(1, 2**34)
EP = ("max_contrast_change", "max_negative_contrast_change")
CASE_IDS = {f"{ep}_small_{k:02}" for ep in EP for k in range(1, 17)}
REC = {"max_negative_contrast_change_small_16", "max_contrast_change_small_13"}
PILOTS = {f"{ep}_{step}" for ep in EP for step in ("small_17", "large_03")}
pins = {}
BOUNDS = {"status": "PASS_INDEPENDENT_CACHED_LOCAL_ROOT"}


def rel(path):
    return path.resolve().relative_to(ROOT).as_posix()


def checked_bytes(path, expect=None):
    # Deliberately forbid all actual mechanics/cache/field payloads here.
    assert path.name not in {"field.json", "verification.json", "cache_packet.json", "recentered_cache_packet.json",
                             "exact_binary_step.json", "assembled_operands.json"}
    assert path.suffix in (".json", ".py", ".md")
    blob = path.read_bytes()
    h = sha256(blob).hexdigest()
    key = rel(path)
    assert pins.get(key, h) == h and (expect is None or expect == h), key
    pins[key] = h
    return blob


def read(path, expect=None):
    return json.loads(checked_bytes(path, expect))


def posint(x):
    return type(x) is int and x > 0


def zero(x):
    return type(x) is int and x == 0


def observed(folder, expected_os):
    os = read(folder / "os_exit_observation.json")
    assert os["status"] == expected_os and posint(os["pid"]) and posint(os["unified_exec_session"])
    assert zero(os["observed_exit_code"])
    end = read(folder / "terminal.json", os["terminal_sha256"])
    assert end["status"] == "FINISHED" and zero(end["exit_code"]) and end["pid"] == os["pid"]
    assert os["receipt_sha256"] == end["receipt_sha256"]
    return read(folder / "receipt.json", os["receipt_sha256"])


def phase_link(folder):
    return {"native_phase": rel(folder / "receipt.json"),
            "native_phase_sha256": pins[rel(folder / "receipt.json")],
            "native_terminal_sha256": pins[rel(folder / "terminal.json")],
            "native_observed_os_exit_sha256": pins[rel(folder / "os_exit_observation.json")]}


def bind_original(row, native):
    assert native["case_id"] == row["case_id"] and native["record_id"] == row["record_id"]
    assert native["field_sha256"] == row["field_sha256"]
    assert Q(native.get("exact_requested_strength", native.get("exact_strength"))) == Q(row["requested_drive_exact"])


def proposal(r, z, sigma, goal):
    assert 0 < r and 0 <= z and 0 <= sigma < goal
    if z >= 1 or sigma + r >= goal:
        return "REVIEW_EXISTING_CUBE_HEADROOM_BEFORE_SPATIAL"
    if r == R34:
        return "TRY_SAME_CENTER_RADIUS_2^-34_AFTER_ALL_GATES"
    if r == R36:
        return "2^-36_PROVED_CACHE_FALLBACK_FRESH_2^-34_DERIVATIVE_FIRST"
    return "ADMIT_CURRENT_RADIUS_OR_FRESH_MATCHED_2^-34_TRIAL"


def entry(row, proof, r, z, eta, posterior, sig, cache, tensor, uh, name, links):
    goal = Q(row["observer_goal_exact"])
    assert 0 < r and 0 <= eta and 0 <= z < 1 and eta + z*r < r
    assert posterior == eta/(1-z) and 0 <= posterior < r and 0 <= sig and sig+posterior < goal
    return {**{k:row[k] for k in ("case_id", "record_id", "endpoint_id", "checkpoint", "material_group", "field",
            "field_sha256", "verification", "verification_sha256", "dose_exact", "requested_drive_exact",
            "solver_binary_drive_exact", "observer_goal_exact", "node_addition_gap_exact")},
            "operator_cache_packet_sha256":cache, "binary_R_sha256":tensor,
            "unsigned_binary_center_sha256":uh,
            "archived_field_center_shift_sigma_exact":str(sig),
            "observer_center": "ARCHIVED_FIELD_SEPARATE_FROM_PROOF_CENTER",
            "spatial_verification_center": "VERIFIED_RECENTERED_U1" if sig else "VERIFIED_ARCHIVED_BINARY_U0",
            "cached_proof_radius_exact":str(r), "cached_ZQ_upper_exact":str(z),
            "cached_eta_upper_exact":str(eta), "cached_posterior_upper_exact":str(posterior),
            "cached_cube_trial_advice":proposal(r,z,sig,goal),
            "native_case":rel(name), "native_case_sha256":pins[rel(name)], **links}


def make():
    inv = read(INV / "receipt.json")
    assert inv["status"] == "PASS_READ_ONLY_36_FIELD_COVERAGE_INVENTORY" and inv["full_five_gate_complete"] is False
    all_rows = {r["case_id"]:r for r in inv["outcome_states"]}
    assert len(all_rows) == len(inv["outcome_states"]) == 36 and set(all_rows) == CASE_IDS | PILOTS
    native = {}
    plans = read(W30 / "preparation/protocol.json")
    assert set(plans["cases"]) == CASE_IDS - {f"{ep}_small_09" for ep in EP}
    for i in range(1,9):
        wave = f"wave_{i:02}"
        folder = W30 / "independent_actual_preparation/actual_replay" / wave
        phase = observed(folder, "PASS_OBSERVED_INDEPENDENT_CACHED_WAVE_TRUE_OS_EXIT")
        assert phase["wave_id"] == wave
        assert phase["pid"] == read(folder/"terminal.json")["pid"]
        for summary in phase["cases"]:
            case = summary["case_id"]
            assert case not in native
            row = all_rows[case]
            name = folder / (case + ".json")
            assert summary["native_file"] == rel(name)
            raw = read(name, summary["native_sha256"])
            bind_original(row, raw)
            chosen = summary["chosen_radius_exact"]
            if case in REC:
                assert raw["status"] == summary["status"] == "UNRESOLVED_INDEPENDENT_CACHED_ROOT" and chosen is None
                native[case] = {"original_unresolved_native_file":rel(name),
                                "original_unresolved_native_sha256":pins[rel(name)]}
                continue
            assert raw["status"] == summary["status"] == "PASS_INDEPENDENT_CACHED_LOCAL_ROOT"
            assert chosen == raw["chosen_radius_exact"]
            at, = [a for a in raw["own_radius_attempts"] if a["radius_exact"] == chosen]
            assert at["whole_cube_AD_recomputed_for_this_radius"] is True
            native[case] = entry(row,raw,Q(chosen),Q(at["ZQ_upper_exact"]),Q(at["eta_upper_exact"]),
              Q(raw["posterior_error_upper_exact"]),Q(0),raw["cache_packet_sha256"],
              raw["per_tensor_sha256"]["preconditioner"],raw["per_tensor_sha256"]["u"],name,phase_link(folder))
    assert set(native) == set(plans["cases"])
    v = read(M / "VERIFIED-CACHED-COVERAGE.json")
    assert v["status"] == "PASS_VERIFIED_TWO_ADDITIONAL_M_CACHED_LOCAL_ROOTS"
    for summary in v["case_results"]:
        case = summary["case_id"]
        assert case not in native and case in {f"{ep}_small_09" for ep in EP}
        name = ROOT / summary["native_proof_file"]
        folder = name.parent
        raw = read(name,summary["native_proof_sha256"])
        if case.startswith(EP[0] + "_"):
            assert raw["status"] == "PASS_INDEPENDENT_M_CACHED_LOCAL_ROOT"
            phase = observed(folder,"PASS_OBSERVED_INDEPENDENT_M_REPLAY_TRUE_OS_EXIT")
            proof, = [r for r in phase["results"] if r["case_id"] == case]
            assert proof["native_sha256"] == summary["native_proof_sha256"]
            eta = raw["independent_eta_upper_exact"]
            u = raw["cache_tensor_content_sha256"]["u"]
        else:
            assert raw["status"] == "PASS_INDEPENDENT_NEGATIVE_M_SMALLER_CACHED_ROOT"
            assert observed(folder,"PASS_OBSERVED_NEGATIVE_M_NARROW_TRUE_OS_EXIT") == raw
            anchor = read(ROOT/raw["reused_signed_cached_anchor_file"],raw["reused_signed_cached_anchor_sha256"])
            bind_original(all_rows[case], anchor)
            eta = raw["eta_upper_exact"]
            u = anchor["cache_tensor_content_sha256"]["u"]
        bind_original(all_rows[case], raw)
        native[case] = entry(all_rows[case],raw,Q(summary["radius_exact"]),Q(summary["cached_ZQ_upper_exact"]),
            Q(eta),Q(summary["posterior_error_upper_exact"]),Q(0),summary["cache_packet_sha256"],
            summary["preconditioner_binary_sha256"],u,name,phase_link(folder))
    folder = W30/"recenter_v1/independent_actual_preparation/actual_replay"
    phase = observed(folder,"PASS_ROOT_OBSERVED_INDEPENDENT_RECENTER_REPLAY_TERMINAL")
    assert phase["status"] == "PASS_INDEPENDENT_TWO_RECENTERED_CACHED_ROOTS_AND_OLD_OBSERVER_BUDGETS"
    assert {c["case_id"] for c in phase["case_results"]} == REC
    for summary in phase["case_results"]:
        case = summary["case_id"]
        row = all_rows[case]
        original = native[case]
        name = folder/(case+".json")
        assert summary["native_file"] == rel(name)
        raw = read(name,summary["native_sha256"])
        assert raw["status"] == "PASS_INDEPENDENT_RECENTERED_CACHED_LOCAL_ROOT_AND_OLD_OBSERVER_BUDGET"
        assert raw["archived_field_sha256"] == row["field_sha256"] and row["record_id"] == raw["record_id"]
        assert raw["old_native_signed_point_sha256"] == original["original_unresolved_native_sha256"]
        assert Q(raw["exact_requested_strength"]) == Q(row["requested_drive_exact"])
        assert Q(raw["archived_addition_gap_xi0_exact"]) == Q(row["node_addition_gap_exact"])
        assert Q(raw["archived_observer_goal_exact"]) == Q(row["observer_goal_exact"])
        native[case] = {**entry(row,raw,Q(raw["radius_exact"]),Q(raw["ZQ_upper_exact"]),Q(raw["eta_upper_exact"]),
                    Q(raw["posterior_error_upper_exact"]),Q(raw["actual_binary_center_shift_sigma_exact"]),
                    raw["recentered_cache_packet_sha256"],raw["original_binary_R_sha256"],
                    raw["replayed_u1_binary_sha256"],name,phase_link(folder)),**original}
    assert set(native) == CASE_IDS
    for pilot in sorted(PILOTS):
        folder = EX/"bridge_reduction_v1/independent_actual_preparation/actual_replay"/pilot
        phase = observed(folder,"PASS_OBSERVED_INDEPENDENT_GLOBAL_REDUCTION_TRUE_OS_EXIT")
        assert phase["status"] == "PASS_INDEPENDENT_EXACT_SPATIAL_ROOT_AND_STRICT_LOCAL_MINIMUM"
        assert phase["pilot_id"] == pilot and phase["full_five_gate_goal_achieved"] is False
    for name in ("transfer_ambiguity_exact_targets","transfer_ambiguity_exact_bridge",
                 "transfer_ambiguity_exact_bridge_reduction","transfer_ambiguity_exact_bridge_reduce_driver_v3",
                 "transfer_ambiguity_cached_derivative_extension","transfer_ambiguity_remaining30_cached",
                 "transfer_ambiguity_recenter_cached","transfer_ambiguity_preconditioner_spd"):
        checked_bytes(ROOT / ("scripts/"+name+".py"))
    checks = {"exact_32_distinct_state_keys_and_separate_four_spatial_pilots":True,
              "all_native_cached_roots_proofs_and_observed_terminals_read":True,
              "two_recentered_states_retain_archived_observer_and_failed_native_attempt":True,
              "all_exact_cached_selfmap_and_posterior_bounds_verified":True,
              "optional_radius_does_not_promote_new_cube":True,
              "same_center_contained_cube_logic":
                  proposal(R34,Q(3,8),Q(0),Q(1,100)) == "TRY_SAME_CENTER_RADIUS_2^-34_AFTER_ALL_GATES"
                  and proposal(R36,Q(1,16),Q(1,2**35),Q(1,100)).startswith("2^-36")
                  and proposal(R34,Q(1),Q(0),Q(1,100)).startswith("REVIEW"),
              "metadata_only_no_actual_mechanics":True}
    assert all(v is True for v in checks.values())
    return {"status":"COLD_METADATA_ONLY_REMAINING32_BRIDGE_DESIGN",
            "root_or_worker_launch_authorized":False,"full_five_gate_complete":False,
            "existing_spatial_roots":4,"existing_cached_roots_among_remaining":32,
            "target_spatial_states_remaining":32,"actual_fields_caches_or_step_vectors_opened":0,
            "actual_AD_geometry_force_H_moments_or_integration_jobs":0,
            "cached_radius_counts":dict(Counter(r["cached_proof_radius_exact"] for r in native.values())),
            "source_sha256":sha256(Path(__file__).read_bytes()).hexdigest(),"checks":checks,
            "cases":list(native[k] for k in sorted(native)), "input_metadata_and_sources_sha256":dict(sorted(pins.items()))}


if __name__ == "__main__":
    result = make()
    with (HERE/"manifest.json").open("x",encoding="utf-8",newline="\n") as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write("\n")
    print(result["cached_radius_counts"])
