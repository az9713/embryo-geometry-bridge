"""Independent prelaunch checks for the cohort exact endpoint reducer."""
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement, product
import json
from pathlib import Path
from random import Random
import sys
from types import SimpleNamespace
from unittest.mock import patch

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
PREP = OUT.parent / "reduction_preparation"
EVIDENCE = ROOT / "evidence/exact_cohort_reducer_prelaunch_independent_review.json"
SOURCE = "scripts/review_cohort_exact_record.py"
FROZEN = {
    SOURCE: "67db748a5d40f31471ddd3d6c84fae280ff7ead650e6eb4e950a4766e54e0285",
    "protocol.json": "cfa83fd31db3a052337abad734d529b94d1e24bf223ef36fc981afd6838b0680",
    "preflight.json": "4cfbd0699c25d827e40635748a5d5af04c124e6aa0b3198f1ebea1c9b3287c8a",
}


def pack(lo, hi=None):
    return {"lower_exact": str(lo), "upper_exact": str(lo if hi is None else hi)}


def h_entry(i, j, lo, hi=None):
    return {"i": i, "j": j, "enclosure": pack(lo, hi)}


def rejects(call, message=None):
    try:
        call()
    except ValueError as error:
        assert message is None or message in str(error), str(error)
    else:
        raise AssertionError("invalid fixture accepted")


def dyadic_oracle(E):
    rng = Random(6172026)
    count = 0
    for n in range(1, 6):
        for _ in range(12):
            g, H = [E.ZERO] * n, {}
            lo_g, hi_g = [Q(0)] * n, [Q(0)] * n
            lo_h = [[Q(0)] * n for _ in range(n)]
            hi_h = [row.copy() for row in lo_h]
            affected = set()
            for size in range(n, 0, -1):
                ids = rng.sample(range(n), size)
                f = []
                for k in ids:
                    a = Q(rng.randint(-19, 4), 16)
                    b = a + Q(rng.randrange(20), 32)
                    lo_g[k] += a
                    hi_g[k] += b
                    f.append(pack(a, b))
                entries = []
                for i, j in combinations_with_replacement(range(size), 2):
                    a = Q(rng.randint(-35, 8), 32)
                    b = a + Q(rng.randrange(35), 16)
                    k, l = ids[i], ids[j]
                    affected.add(tuple(sorted((k, l))))
                    lo_h[k][l] += a
                    hi_h[k][l] += b
                    if k != l:
                        lo_h[l][k] += a
                        hi_h[l][k] += b
                    entries.append(h_entry(i, j, a, b))
                rng.shuffle(entries)
                E._scatter(g, H, ids, f, entries)
            assert g == list(zip(lo_g, hi_g))
            assert H == {(i, j): (lo_h[i][j], hi_h[i][j]) for i, j in affected}
            h_abs = [sum(max(abs(a), abs(b)) for a, b in zip(l, r)) for l, r in zip(lo_h, hi_h)]
            assert E._h_rows(n, H) == h_abs
            R = [[rng.randint(-30, 30) / 2**rng.randrange(12) for i in range(n)] for j in range(n)]
            as_q = [[Q.from_float(v) for v in row] for row in R]
            vertices = list(product(*g))
            evalu = [[sum(v * x for v, x in zip(row, point)) for point in vertices] for row in as_q]
            expected_g = [(min(vals), max(vals)) for vals in evalu]
            expected_h = [sum(abs(v)*h for v, h in zip(row, h_abs)) for row in as_q]
            expected_norm = max(sum(map(abs, row)) for row in as_q)
            assert E._dyadic_reduce(R, g, h_abs) == (expected_g, expected_h, expected_norm)
            count += 1
    rejects(lambda: E._dyadic_reduce([[1.0]], [(Q(1,3), Q(1,3))], [Q(0)]), "dyadic")
    rejects(lambda: E._scatter([E.ZERO], {}, [0], [pack(0)], [h_entry(0, 0, 0), h_entry(0, 0, 0)]), "unique")
    rejects(lambda: E._scatter([E.ZERO], {}, [0], [pack(0)], []), "coverage")
    rejects(lambda: E._packed(pack(True)), None)
    return count


def producer_fixture(E):
    row = [(Q(-1,8), Q(1,16))]
    H = {(0,0):(Q(-1,16), Q(1,16))}
    data = {
        "force_integral_minus_cache": [pack(-Q(1,4), Q(1,8))],
        "preconditioned_force_defect": [pack(-Q(1,4), Q(1,8))],
        "H_absolute_row_bounds": [pack(Q(1,8), Q(1,8))],
        "H_integral_minus_cache_upper": [h_entry(0,0,-Q(1,8),Q(1,8))],
        "RF_defect_uniform_upper_exact": "1/4", "b_upper_exact": "1/8",
        "eta_upper_exact": "3/8", "z_upper_exact": "3/8", "selfmap_lhs_exact": "3/4",
        "r_exact": "1", "sufficient_conditions_pass": True,
        "status": "CANDIDATE_UNIFORM_EXACT_FE_ROOT_UNREVIEWED",
        "root_assertion_independently_reviewed": False, "root_error_upper_exact": "3/5",
    }
    test = lambda p: E._verify_producer(p, row, H, row, [Q(1,16)], Q(1,8), Q(1,16), Q(1,8), Q(1,4), Q(1))
    assert test(data)
    rejects(lambda: test({**data, "H_absolute_row_bounds": [pack(Q(1,32))]}), "domination")
    rejects(lambda: test({**data, "H_integral_minus_cache_upper": [h_entry(0,0,-Q(1,32), Q(1,8))]}), "containment")
    rejects(lambda: test({**data, "preconditioned_force_defect": [pack(-Q(1,16), Q(1,8))]}), "containment")
    rejects(lambda: test({**data, "RF_defect_uniform_upper_exact": "1/16"}), "domination")
    rejects(lambda: test({**data, "root_error_upper_exact": "1/2"}), "arithmetic")
    assert E._strict_test(Q(1,4), Q(1,4), Q(1,8), Q(0), Q(1,2))[-2:] == (False, None)
    assert E._strict_test(Q(0), Q(0), Q(0), Q(1), Q(1))[-2:] == (False, None)


def admission_fixture(E):
    record = "optoGEF2_small_col_L"
    folder = E._record_dir(record)
    source = "scripts/cohort_exact_pilot.py"
    support = "synthetic_support.json"
    h = lambda path: sha256(path.encode()).hexdigest()
    pins = {k: h(f"{folder}/{v}") for k, v in E.PRIMARY_PINS.items()}
    steps = []
    documents = {}

    def success(results_digest):
        return {"status": "FINISHED", "exit_code": 0, "pid": 42, "results_sha256": results_digest}

    def fake_read(path, digest, bindings, parse=True):
        assert all(not path.startswith(f"{folder}/{part}/") for part in ("bulk", "surface"))
        assert digest == h(path), path
        bindings[path] = digest
        steps.append(path)
        return documents[path] if parse else None

    def write(path, value):
        assert path.name == record + "_input_admission.json"
        steps.append(value)
        return h(path.name)

    def put(name, value):
        documents[f"{folder}/{name}"] = value

    put("protocol.json", {"record_id": record, "source_and_input_hashes": {source: h(source)}})
    put("domains.json", {"record_id": record, "status": "PASS_FRESH_DOMAINS_AND_SELECTORS_UNREVIEWED",
                         "actual_force_H_integrals_run": False, "protocol_sha256": pins["protocol"],
                         "body_sha256": h(f"{folder}/body_geometry.json")})
    put("domains_terminal.json", success(pins["domains"]))
    put("integration_prelaunch_independent_review.json", {
        "record_id": record, "status": "PASS", "checks": {"fixture": True}, "domains_result_sha256": pins["domains"],
        "protocol_sha256": pins["protocol"], "body_sha256": h(f"{folder}/body_geometry.json"),
        "preflight_sha256": h(f"{folder}/preflight.json"), "support_file": support,
        "support_sha256": h(support), "source_sha256": h(source),
    })
    packets = {}
    for part, names in (
        ("bulk", [f"element_{i:03d}.json" for i in range(144)]),
        ("surface", [f"face{i}_edge{j:02d}.json" for i in (0,1) for j in range(12)]),
    ):
        files = [f"{folder}/{part}/{name}" for name in names]
        packets[part] = [{"file": name, "sha256": h(name)} for name in files]
        put(f"{part}_terminal.json", success(h(f"{folder}/{part}_results.json")))
        put(f"{part}_results.json", {"record_id": record, "part": part,
                                    "domains_sha256": pins["domains"], "status": "COMPONENT_PART_COMPLETE_UNREVIEWED",
                                    "files": packets[part]})
    put("assembled.json", {"record_id": record, "protocol_sha256": pins["protocol"],
                           "domains_sha256": pins["domains"], "component_packets": packets["bulk"] + packets["surface"]})
    with patch.object(E, "_prepared", side_effect=lambda: ({}, {})), \
            patch.object(E, "_bound", side_effect=fake_read), patch.object(E, "_write", side_effect=write):
        E.admit(record, pins)
        assert isinstance(steps[-1], dict) and steps[-1]["expected_packets"] == packets
        assert steps[-1]["local_component_payloads_opened"] is False
        n = len([x for x in steps if isinstance(x, dict)])
        rejects(lambda: E.admit(record, {k:v for k,v in pins.items() if k != "domains"}), "external")
        review = documents[f"{folder}/integration_prelaunch_independent_review.json"]
        review["record_id"] = "optoGEF2_small_col_F"
        rejects(lambda: E.admit(record, pins), "independently reviewed")
        review["record_id"] = record
        review["checks"] = {"non_boolean": 1}
        rejects(lambda: E.admit(record, pins), "independently reviewed")
        review["checks"] = {"fixture": True}
        manifest = documents[f"{folder}/bulk_results.json"]
        saved = manifest["files"]
        manifest["files"] = saved[:-1] + saved[:1]
        rejects(lambda: E.admit(record, pins), "canonical")
        manifest["files"] = saved
        assembled = documents[f"{folder}/assembled.json"]
        assembled["component_packets"] = assembled["component_packets"][::-1]
        rejects(lambda: E.admit(record, pins), "exactly match")
        assert len([x for x in steps if isinstance(x, dict)]) == n


def all_cached_energy_premises(E, raw):
    bindings = {}
    authority = json.loads(raw(E.C.AUTHORITY, E.C.AUTHORITY_SHA256))
    energy = json.loads(raw(E.ENERGY_RESULT, E.PREPARATION_PINS[E.ENERGY_RESULT]))
    energy_review = json.loads(raw(E.ENERGY_REVIEW, E.PREPARATION_PINS[E.ENERGY_REVIEW]))
    records = []
    for row in authority["records"]:
        idx = row["record_id"]
        c = E.C.load_context(idx)
        for name, expected in c.input_hashes.items():
            raw(name, expected)
        root = json.loads(raw(row["root_certificate_file"], row["root_certificate_sha256"]))
        capsules = json.loads(raw(row["independent_root_capsule_file"], row["independent_root_capsule_sha256"]))
        original, = [a for a in capsules["record_audits"] if a["record_id"] == idx]
        audit = json.loads(raw(original["review_file"], original["review_sha256"]))
        assert audit == original["audit"]
        if "bounds" in audit:
            source, bounds = original["review_file"], audit["bounds"]
            eta_key, z_key = "eta", "Z"
        else:
            source = audit["reuse_review_file"]
            nested = json.loads(raw(source, audit["reuse_review_sha256"]))
            bounds = nested["independent_bounds"]
            eta_key, z_key = "eta_uniform", "Z_uniform"
        eta = max(Q(bounds[eta_key]["upper_exact"]), Q(root["uniform_cell_result"]["eta_uniform"]["upper_exact"]))
        z = Q(bounds[z_key]["upper_exact"])
        r = Q(root["uniform_cell_result"]["selected_radius_exact"])
        assert 0 <= z < 1 and eta + z*r < r
        p = {"cached_bound_provider": source, "cached_eta_upper_exact": str(eta),
             "cached_Z_upper_exact": str(z), "normal_radius_exact": str(r)}
        assert E._old_premises(c, p, bindings) == (eta, z, r, source)
        match, = [a for a in energy["records"] if a["record_id"] == idx]
        result = json.loads(raw(match["file"], match["sha256"]))
        independent, = [a for a in energy_review["record_audits"] if a["record_id"] == idx]
        compared = json.loads(raw(independent["review_file"], independent["review_sha256"]))
        assert compared == independent["audit"] and all(v is True for v in compared["checks"].values())
        norm = Q(compared["exact_R_infinity_norm"])
        binary_R = sha256(c.R.tobytes(order="C")).hexdigest()
        assert binary_R == result["root_binding"]["preconditioner_binary_sha256"]
        assert result["root_binding"]["independent_Z_source"]["file"] == source
        assert Q(compared["uniform_Hessian_lower_bound_exact"]) > 0
        choice = z + (1-z)/3
        reused = E._minimum(c, bindings, z, choice, r, norm, source)
        assert reused["strict_local_minimum"] is True
        assert Q(reused["normalized_curvature_lower_exact"]) == (1-choice)/norm
        records.append({"record_id": idx, "cached_bound_provider": source, "full_tube_Z_exact": str(z),
                        "same_binary_R_sha256": binary_R, "cached_SPD_review_sha256": independent["review_sha256"]})
    assert len(records) == 17 and len({item["record_id"] for item in records}) == 17
    for name, expected in bindings.items():
        raw(name, expected)
    return records


def run():
    if not __debug__:
        raise RuntimeError("Replay requires assertions")
    bindings = {}

    def raw(path, expected=None):
        text = (ROOT / path).read_bytes()
        digest = sha256(text).hexdigest()
        assert expected is None or digest == expected, path
        assert bindings.setdefault(path, digest) == digest, path
        return text

    prep = PREP.relative_to(ROOT).as_posix()
    p = json.loads(raw(f"{prep}/protocol.json", FROZEN["protocol.json"]))
    preflight = json.loads(raw(f"{prep}/preflight.json", FROZEN["preflight.json"]))
    assert raw(SOURCE, FROZEN[SOURCE]) == raw(f"{prep}/review_cohort_exact_record_frozen.py", FROZEN[SOURCE])
    for name, expected in p["input_hashes"].items():
        raw(name, expected)
    for path, expected in p["superseded_preparations"].items():
        raw(f"{prep}/{path}", expected)
    assert preflight["status"] == "PASS_REDUCTION_PREPARATION_ONLY"
    assert preflight["protocol_sha256"] == FROZEN["protocol.json"] and preflight["source_sha256"] == FROZEN[SOURCE]
    assert len(preflight["checks"]) == preflight["check_count"] and all(x is True for x in preflight["checks"].values())
    assert p["real_component_payloads_opened"] is p["real_record_review_run"] is False
    assert preflight["real_record_review_run"] is preflight["component_integration_run"] is False

    sys.path.insert(0, str(ROOT / "scripts"))
    import review_cohort_exact_record as E
    assert E._prepared()[0] == p and E.PRIMARY_PINS == p["actual_admission_requires_external_primary_pins"]
    replay_file = Path(__file__).relative_to(ROOT).as_posix()
    raw(replay_file)
    count = dyadic_oracle(E)
    producer_fixture(E)
    admission_fixture(E)
    rows = all_cached_energy_premises(E, raw)
    for name, expected in tuple(bindings.items()):
        raw(name, expected)
    report = {
        "status": "PASS",
        "scope": "Prelaunch endpoint reducer and cached SPD premise review. No actual local component payloads opened, endpoint reduction, assembly, or new root result.",
        "source_sha256": FROZEN[SOURCE],
        "protocol_sha256": FROZEN["protocol.json"],
        "preflight_sha256": FROZEN["preflight.json"],
        "checks": {
            "frozen_reducer_and_author_protocol_preflight_current": True,
            "exact_scatter_and_signed_dyadic_products_agree_with_independent_corner_oracle": True,
            "producer_signed_containment_and_widened_upper_bound_domination": True,
            "strict_selfmap_and_contraction_equalities_rejected": True,
            "mock_admission_binds_complete_manifests_without_component_payload_reads": True,
            "mock_admission_rejects_cross_record_nonboolean_missing_duplicate_or_reordered_inputs": True,
            "all17_cached_SPD_premises_bind_same_context_R_Z_radius_and_curvature": True,
            "homotopy_and_inverse_norm_argument_reviewed_for_original_cube": True,
        },
        "dyadic_corner_fixtures": count,
        "cached_SPD_records": len(rows),
        "caller_requirements": [
            "Use externally observed current primary hashes and retain the exact admission file digest produced by admit before invoking review. The separately pinned review here assesses the frozen reducer.",
            "The reducer trusts previously reviewed local integration and domain results. Containment is rechecked at the assembly layer, with no claim to independently replay the shared integral kernels.",
        ],
        "actual_local_payloads_opened": False, "actual_endpoint_reduction_run": False,
        "review_source_file": replay_file, "review_source_sha256": bindings[replay_file],
        "input_hashes": dict(sorted(bindings.items())),
    }
    if "--freeze" in sys.argv:
        receipt = OUT / "replay.json"
        assert not receipt.exists() and not EVIDENCE.exists()
        with receipt.open("x", encoding="utf8", newline="\n") as stream:
            stream.write(json.dumps({**report, "cached_SPD_detail": rows}, indent=2) + "\n")
        report["review_receipt_file"] = receipt.relative_to(ROOT).as_posix()
        report["review_receipt_sha256"] = sha256(receipt.read_bytes()).hexdigest()
        with EVIDENCE.open("x", encoding="utf8", newline="\n") as stream:
            stream.write(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "checks": len(report["checks"]),
                      "dyadic_corner_fixtures": count, "cached_SPD_records": len(rows)}), flush=True)


if __name__ == "__main__":
    run()
