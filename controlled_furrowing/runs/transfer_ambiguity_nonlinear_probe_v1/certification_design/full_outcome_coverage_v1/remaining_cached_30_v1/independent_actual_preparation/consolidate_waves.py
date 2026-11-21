"""Consolidate recorded slots; use only exact comparisons, with no numerical replay."""
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
STUDY = next(parent for parent in HERE.parents if (parent / "PLAN.md").is_file())
PRODUCER_INDEX_SHA = "610910bb12e706897a6f42bd0fe457d2f33fdad8a77a7efd090176a21551e7c3"
WAVE_NATIVE_AND_OS = (
    ("872e395a29afb5389bb1168b79674016b1224e78af462db86376b947f7210f60", "c5e6b2e36dd005243ebb812b3b825f8aebf1bfc634b5d54737ee2705587bf036"),
    ("485b34ad2c838b56e091b14823c663de931797739569080ab8a2444850ea38ac", "42fe685976a2eaf41c5e5f3d6b4f418c06fc74093c751cf1c390729305da1fb5"),
    ("9a9d40746c777fe6e48776ba4e9a14aa2853d1ebacfb732bb2d9319861aeeef2", "5b7636f8ec7bedd6fd7a3502c17f9e32cc723a84dbf8e273650d7049d1df35bf"),
    ("9283218568a78812f83cda513a47f1cc34d55314caad8b57b6db7eba7e6f8fe6", "4f7e0968f8462570e505b98135192d4dbb2bb52c7208538fec81c32c237439e8"),
    ("77ea54a7f2ddf8ea18d8056d48cdbdec80e988a3399301319dd1ecc96fdfd5a7", "4103bfa454b42e7c18450d0d2613d98ae115072e5c4612920112ad129b165ce1"),
    ("11b383b758da997661445bc8aa08ab06606d666479c1630aad5b808ed916e055", "93942ac8e9a12367116872411ffc376ed338fed9d4687638178ecc4473e80122"),
    ("87d8d4308964be3b7bd4e3b51c768aba770cdb7e70f5580171c453b538d18d01", "e39ceecf09c3116b623cd8fb182d8496cb08ce43edb3778e2b07c3891e0e4b41"),
    ("b0adfa31be66227ef6de06a32740cfbe9b4727d32972e2550c7b77e29af493a7", "722704f03a0b3498584057cf39488174023464f67cc893da065cb4d197aabfe8"),
)
bindings = {}


def read(path, expected=None, parse=True):
    path = (STUDY / path).resolve()
    payload = path.read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    assert expected is None or digest == expected, path
    bindings[path.relative_to(STUDY).as_posix()] = digest
    return json.loads(payload) if parse else digest


waves = []
seen = {}
read(HERE.parent / "actual" / "producer_wave_completion_index.json", PRODUCER_INDEX_SHA)
for i, (native_sha, os_sha) in enumerate(WAVE_NATIVE_AND_OS, 1):
    wave_id = f"wave_{i:02d}"
    base = HERE / "actual_replay" / wave_id
    observed = read(base / "os_exit_observation.json", os_sha)
    assert observed["wave_id"] == wave_id
    assert type(observed["observed_exit_code"]) is int and observed["observed_exit_code"] == 0
    assert observed["receipt_sha256"] == native_sha
    terminal = read(base / "terminal.json", observed["terminal_sha256"])
    assert terminal["wave_id"] == wave_id and terminal["pid"] == observed["pid"]
    assert type(terminal["exit_code"]) is int and terminal["exit_code"] == 0
    assert terminal["status"] == "FINISHED" and terminal["receipt_sha256"] == native_sha
    native = read(base / "receipt.json", native_sha)
    assert native["wave_id"] == wave_id and native["pid"] == observed["pid"]
    assert native["exact_spatial_root_count"] == 0 and native["new_equilibrium_solves_or_preconditioners"] == 0
    for key, file in (("source_sha256", "replay.py"), ("protocol_sha256", "protocol.json"), ("preflight_sha256", "preflight.json")):
        assert terminal[key] == native[key]
        read(HERE / file, native[key], parse=False)
    admission_sha = native["root_terminal_admission_sha256"]
    assert terminal["root_terminal_admission_sha256"] == admission_sha
    read(HERE.parent / "actual" / wave_id / "root_terminal_admission.json", admission_sha)
    resolved = 0
    for row in native["cases"]:
        case = row["case_id"]
        assert case not in seen and Path(row["native_file"]).resolve().parent == base
        point = read(row["native_file"], row["native_sha256"])
        for key in ("case_id", "record_id", "status", "field_sha256", "cache_packet_sha256", "chosen_radius_exact", "posterior_error_upper_exact"):
            assert row[key] == point[key], (case, key)
        read(point["field"], point["field_sha256"], parse=False)
        assert point["rounded_strength_used_in_proof"] is False and point["new_R_or_solver_evaluations"] == 0
        goal = Q(point["observer_radius_budget_exact"])
        failed = []
        accepted = []
        for step in point["own_radius_attempts"]:
            r = Q(step["radius_exact"])
            z = Q(step["ZQ_upper_exact"])
            eta = Q(step["eta_upper_exact"])
            lhs = Q(step["strict_selfmap_lhs_upper_exact"])
            assert r > 0 and z >= 0 and eta >= 0 and lhs == eta + z*r
            assert step["whole_cube_AD_recomputed_for_this_radius"] is True
            if z < 1 and lhs < r:
                posterior = Q(step["posterior_error_upper_exact"])
                assert 0 <= posterior == eta/(1-z) < min(goal, r)
                accepted.append((r, posterior))
            else:
                assert step["posterior_error_upper_exact"] is None
                failed.append({key: step[key] for key in ("radius_exact", "eta_upper_exact", "ZQ_upper_exact", "strict_selfmap_lhs_upper_exact", "display")})
        is_resolved = point["status"] == "PASS_INDEPENDENT_CACHED_LOCAL_ROOT"
        assert is_resolved == bool(accepted) and len(accepted) <= 1
        if is_resolved:
            resolved += 1
            assert accepted == [(Q(row["chosen_radius_exact"]), Q(row["posterior_error_upper_exact"]))]
        else:
            assert failed and row["chosen_radius_exact"] is None and row["posterior_error_upper_exact"] is None
        seen[case] = {**row, "wave_id": wave_id, "exact_requested_strength": point["exact_requested_strength"],
                      "observer_radius_budget_exact": point["observer_radius_budget_exact"],
                      "original_unresolved_attempts": failed}
    assert resolved == native["cached_local_root_count"]
    waves.append({"wave_id": wave_id, "cached_local_roots": resolved, "cases": len(native["cases"]),
                  "independent_native_sha256": native_sha, "observed_exit_record_sha256": os_sha})
unresolved = sorted(case for case, row in seen.items() if row["chosen_radius_exact"] is None)
assert len(seen) == 30 and unresolved == ["max_contrast_change_small_13", "max_negative_contrast_change_small_16"]
read(Path(__file__), parse=False)
packet = {
    "status": "PASS_CACHED_WAVE_ACCOUNTING_WITH_TWO_UNRESOLVED",
    "scope": "Independent cached quadrature roots on the fixed singleton parameters, using the frozen binary cache and preconditioner.",
    "covered_candidate_count": 30,
    "strict_cached_local_root_count": 28,
    "unresolved_case_ids": unresolved,
    "exact_spatial_integral_root_count_from_this_receipt": 0,
    "full_five_gate_goal_achieved": False,
    "waves": waves,
    "cases": list(seen.values()),
    "source_and_receipt_sha256": bindings,
}
dest = HERE / "all_waves_readback.json"
with dest.open("x", encoding="utf-8", newline="\n") as stream:
    stream.write(json.dumps(packet, indent=2, ensure_ascii=False) + "\n")
print(f"{dest.relative_to(STUDY)} {hashlib.sha256(dest.read_bytes()).hexdigest()} {len(bindings)} digests")

