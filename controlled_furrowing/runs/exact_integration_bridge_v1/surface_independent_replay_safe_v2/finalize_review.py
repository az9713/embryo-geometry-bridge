"""Hash/interval readback of the one completed SAFE six-edge replay."""
from pathlib import Path
from fractions import Fraction as Q
import datetime
import hashlib
import json

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]


def read(p):
    return json.loads(p.read_text(encoding="utf8"))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write(p, obj):
    if p.exists():
        raise RuntimeError("write-once output exists: " + str(p))
    p.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="utf8")


def interval(d):
    lo, hi = Q(d["lower_exact"]), Q(d["upper_exact"])
    if lo > hi:
        raise AssertionError("unordered interval")
    return lo, hi


def main():
    p = read(OUT / "protocol.json")
    results = read(OUT / "results.json")
    job = read(OUT / "job.json")
    tests = read(OUT / "selfchecks.json")
    checks = {}
    bindings = {}

    def check(k, v):
        if v is not True:
            raise AssertionError(k)
        checks[k] = True

    def bind(file, expected=None):
        digest = sha(file)
        if expected is not None:
            check("hash_" + file.relative_to(ROOT).as_posix(), digest == expected)
        bindings[file.relative_to(ROOT).as_posix()] = digest

    for name, digest in p["input_hashes"].items():
        bind(ROOT / name, digest)
    for name in ("protocol.json", "results.json", "job.json", "selfchecks.json", "replay.py"):
        bind(OUT / name)
    check("true_terminal_retained_job", job["status"] == "FINISHED" and job["exit_code"] == 0)
    check("completed_six_safe_replays", results["status"] == "PASS_SIX_SAFE_V2_INDEPENDENT_POINT_SURFACE_REPLAYS" and len(results["outputs"]) == 6)
    check("source_and_protocol_current", results["source_sha256"] == sha(OUT / "replay.py") and results["protocol_sha256"] == sha(OUT / "protocol.json"))
    check("precision_and_order_frozen", p["precision_bits"] == 640 and p["truncation_order"] == 28)
    check("own58_zero_crossing_maps_domain_checks", len(tests["checks"]) == 58 and all(v is True for v in tests["checks"].values()))
    check("reused113_safe_synthetic_checks", len(tests["explicitly_reused_safe_surface_synthetic_checks"]) == 113 and all(v is True for v in tests["explicitly_reused_safe_surface_synthetic_checks"].values()))
    check("selfcheck_source_current", tests["source_sha256"] == sha(OUT / "replay.py"))
    seen = []
    force_count = 0
    hessian_count = 0
    for item in results["outputs"]:
        file = OUT / item["file"]
        bind(file, item["sha256"])
        record = read(file)
        face, edge = record["face"], record["edge"]
        seen.append((face, edge))
        producer = ROOT / f"runs/exact_integration_bridge_v1/surface_assembly_B_safe_v2/edges/face{face}_edge{edge:02d}.json"
        bind(producer, record["producer_file_sha256"])
        q = read(producer)
        label = f"{face}_{edge}"
        check("record_PASS_" + label, record["status"] == "PASS_SAFE_V2_INDEPENDENT_POINT_SURFACE_EDGE_REPLAY" and bool(record["checks"]) and all(v is True for v in record["checks"].values()))
        check("same_global_map_" + label, record["global_indices"] == q["global_indices"])
        check("same_force_count_" + label, len(record["force_enclosures"]) == len(q["force_integral_minus_cache"]))
        for i, (a, b) in enumerate(zip(record["force_enclosures"], q["force_integral_minus_cache"])):
            x, y = interval(a), interval(b)
            check(f"rational_force_overlap_{label}_{i}", max(x[0], y[0]) <= min(x[1], y[1]))
            force_count += 1
        H = {(x["i"], x["j"]): x["enclosure"] for x in q["H_integral_minus_cache_upper"]}
        independent_H = {(x["i"], x["j"]): x["enclosure"] for x in record["Hessian_enclosures"]}
        check("same_Hessian_keys_" + label, set(H) == set(independent_H))
        for (i, j), enclosure in independent_H.items():
            x, y = interval(enclosure), interval(H[i, j])
            check(f"rational_Hessian_overlap_{label}_{i}_{j}", max(x[0], y[0]) <= min(x[1], y[1]))
            hessian_count += 1
    check("same_six_predeclared_cases", seen == [(f, e) for f in (0, 1) for e in (0, 1, 11)])
    check("all30force92Hessian_overlaps", force_count == 30 and hessian_count == 92)
    check("not_a_root_transfer", results["root_transfer_proved"] is False)
    execution = {
        "status": "FINISHED",
        "unified_exec_session": 49333,
        "observed_pid": 6296,
        "observed_creation_time": "2026-11-19T13:48:32.214266-08:00",
        "actual_unified_exec_exit_code": 0,
        "result_sha256": sha(OUT / "results.json"),
        "elapsed_seconds": results["elapsed_seconds"],
    }
    write(OUT / "execution_receipt.json", execution)
    bind(OUT / "execution_receipt.json")
    bind(Path(__file__))
    receipt = {
        "status": "PASS_SAFE_V2_SIX_INDEPENDENT_POINT_SURFACE_REPLAYS",
        "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "checks": checks,
        "bound_files": bindings,
        "own_preflight_checks": 58,
        "explicitly_reused_safe_surface_synthetic_checks": 113,
        "runtime_case_checks": sum(item["checks"] for item in results["outputs"]),
        "independent_rational_overlap_readback": {"force": force_count, "Hessian_upper": hessian_count},
        "actual_reintegration_coverage": {"edges": 6, "all_B_surface_edges": 24, "faces": [0, 1], "edge_indices_each_face": [0, 1, 11]},
        "precision_bits": 640,
        "truncation_order": 28,
        "reused_arithmetic": "Reviewed SAFE exact-zero polynomial arithmetic, rational moments, scalar constitutive primitives, SAFE binomial tails. Own raw P2 construction, generic jet AD, Horner/full-product contractions and normalized-cache correction. This is not a second independent arithmetic implementation.",
        "old_unsafe_replay_remains_withdrawn": True,
        "complete_point_certificate_restored_by_this_receipt": False,
        "uniform_root_transfer_proved": False,
        "new_equilibrium_or_optimization_run": False,
        "scope": "Six predeclared B point surface edges only; no full24 independent reintegration, parameter/corrector-tube/root, continuum, biological fit or novelty claim.",
    }
    dest = ROOT / "evidence/exact_surface_safe_independent_review.json"
    write(dest, receipt)
    print(json.dumps({"status": receipt["status"], "checks": len(checks), "bindings": len(bindings), "results_sha256": sha(OUT / "results.json"), "receipt_sha256": sha(dest), "elapsed_seconds": results["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
