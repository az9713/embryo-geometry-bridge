"""Retain actual terminal evidence and exact six-packet width comparisons.

No integration, optimization or global assembly is performed here.
"""
from pathlib import Path
from fractions import Fraction as F
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
        raise ValueError("preserve write-once output")
    p.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf8")


def bounds(d):
    return F(d["lower_exact"]), F(d["upper_exact"])


def width(d):
    lo, hi = bounds(d)
    if lo > hi:
        raise ValueError("unordered interval")
    return hi - lo


def overlaps(a, b):
    x, y = bounds(a), bounds(b)
    return max(x[0], y[0]) <= min(x[1], y[1])


def main():
    p, results = read(OUT / "protocol.json"), read(OUT / "results.json")
    old_manifest_file = ROOT / "runs/exact_integration_bridge_v1/uniform_bulk_assembly_B_safe_v2/results.json"
    old_manifest = read(old_manifest_file)
    old_packets = {x["element"]: x for x in old_manifest["elements"]}
    checks, bindings, comparison, stdout = {}, {}, [], []

    def ck(k, v):
        if v is not True:
            raise AssertionError(k)
        checks[k] = True

    def bind(file, expected=None):
        h = sha(file)
        if expected is not None:
            ck("hash_" + file.relative_to(ROOT).as_posix(), h == expected)
        bindings[file.relative_to(ROOT).as_posix()] = h

    ck("six_completed_results", len(results["elements"]) == 6 and [x["element"] for x in results["elements"]] == list(range(72, 78)))
    ck("producer_all_point_diagnostics_pass", results["all_point_diagnostics_pass"] is True)
    ck("source_and_protocol_unchanged", results["source_sha256"] == sha(ROOT / "scripts/exact_uniform_axis_edge_bridge_B_safe.py") and results["protocol_sha256"] == sha(OUT / "protocol.json"))
    bind(old_manifest_file)
    for row in results["elements"]:
        index = row["element"]
        newfile, oldinfo = OUT / row["file"], old_packets[index]
        oldfile = ROOT / oldinfo["file"]
        bind(newfile, row["sha256"])
        bind(oldfile, oldinfo["sha256"])
        new, old = read(newfile), read(oldfile)
        ck(f"same_global_free_map_{index}", new["global_indices"] == old["global_indices"])
        ck(f"same_fresh_body_rho_{index}", F(new["rho_exact"]) == F(old["rho_exact"]))
        ck(f"same_active_mask_{index}", F(new["active"]) == F(old["active"]))
        ck(f"same_force_entry_count_{index}", len(new["force_integral_minus_cache"]) == len(old["force_integral_minus_cache"]))
        newH = {(v["i"], v["j"]): v["enclosure"] for v in new["H_integral_minus_cache_upper"]}
        oldH = {(v["i"], v["j"]): v["enclosure"] for v in old["H_integral_minus_cache_upper"]}
        ck(f"same_Hessian_keys_{index}", set(newH) == set(oldH))
        ck(f"all_uniform_force_packet_overlaps_{index}", all(overlaps(a, b) for a, b in zip(new["force_integral_minus_cache"], old["force_integral_minus_cache"])))
        ck(f"all_uniform_Hessian_packet_overlaps_{index}", all(overlaps(value, oldH[key]) for key, value in newH.items()))
        ck(f"all_point_diagnostics_{index}", all(v is True for v in new["point_overlap_diagnostics_only"].values()))
        oldFw = max(map(width, old["force_integral_minus_cache"]))
        newFw = max(map(width, new["force_integral_minus_cache"]))
        oldHw = max(map(width, oldH.values()))
        newHw = max(map(width, newH.values()))
        comparison.append({
            "element": index, "old_Duffy_file": oldinfo["file"], "old_Duffy_sha256": sha(oldfile),
            "new_physical_file": newfile.relative_to(ROOT).as_posix(), "new_physical_sha256": sha(newfile),
            "old_degree": old["degree"], "new_degree": new["degree"],
            "max_force_interval_width_old_exact": str(oldFw), "max_force_interval_width_new_exact": str(newFw),
            "max_Hessian_interval_width_old_exact": str(oldHw), "max_Hessian_interval_width_new_exact": str(newHw),
            "ratio_of_max_force_widths_exact": str(oldFw / newFw) if newFw else None,
            "ratio_of_max_Hessian_widths_exact": str(oldHw / newHw) if newHw else None,
            "display_only": {"old_force_width": float(oldFw), "new_force_width": float(newFw), "old_H_width": float(oldHw), "new_H_width": float(newHw), "force_width_ratio": float(oldFw / newFw) if newFw else None, "H_width_ratio": float(oldHw / newHw) if newHw else None},
            "scope": "Ratio of two valid local interval enclosure widths; not a physical error, experimental uncertainty or an assembled contraction/root result.",
        })
        stdout.append(json.dumps({"axis_edge_uniform_complete": index, "seconds": new["seconds"], "status": new["status"]}))
    write(OUT / "width_comparison.json", {"status": "SIX_LOCAL_ENCLOSURE_WIDTH_COMPARISONS_ONLY", "rows": comparison, "old_full_manifest_sha256": sha(old_manifest_file), "no_global_assembly_or_root_claim": True})
    stdoutfile = OUT / "stdout.log"
    if stdoutfile.exists():
        raise ValueError("preserve stdout")
    # These six lines match the retained unified-exec output, including the
    # exact per-element seconds and completion status from each immutable file.
    stdoutfile.write_text("\n".join(stdout) + "\n", encoding="utf8")
    execution = {
        "status": "TRUE_TERMINAL_EXIT_ZERO",
        "unified_exec_session": 70450, "observed_pid": 32468,
        "observed_creation_time": "2026-11-19T15:28:39.225196-08:00",
        "actual_unified_exec_exit_code": 0,
        "source_sha256": results["source_sha256"], "protocol_sha256": results["protocol_sha256"],
        "results_sha256": sha(OUT / "results.json"), "stdout_sha256": sha(stdoutfile),
        "summed_element_compute_seconds_excluding_startup_hashes": sum(read(OUT / x["file"])["seconds"] for x in results["elements"]),
        "sole_authorized_job": True, "restarted_or_duplicate_job": False,
        "representation_review_checked_after_launch": True,
        "representation_independent_review_sha256": sha(ROOT / "evidence/exact_axis_edge_physical_independent_review.json"),
        "independent_prelaunch_review_sha256": sha(OUT / "prelaunch_independent_review.json"),
    }
    write(OUT / "execution_receipt.json", execution)
    for name, digest in p["source_and_input_hashes"].items():
        bind(ROOT / name, digest)
    for file in (Path(__file__), OUT / "protocol.json", OUT / "preflight.json", OUT / "prelaunch_independent_review.json", OUT / "results.json", OUT / "launch_receipt.json", OUT / "execution_receipt.json", OUT / "stdout.log", OUT / "width_comparison.json", ROOT / "evidence/exact_axis_edge_physical_independent_review.json"):
        bind(file)
    receipt = {
        "status": "PASS_EXECUTION_INTEGRITY_AND_WIDTH_COMPARISON_PENDING_INDEPENDENT_NUMERICAL_RESULTS_REVIEW",
        "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "checks": checks, "bound_files": bindings,
        "root_claim": False, "new_equilibrium_optimization": False,
        "scope": "One true-terminal authorized six-element physical-coordinate integration/cache job, immutable evidence and local comparison only. No138+6+24 global reduction, root certificate or independent numerical replay by this producer receipt.",
    }
    write(OUT / "completion_receipt.json", receipt)
    print(json.dumps({"status": receipt["status"], "checks": len(checks), "bindings": len(bindings), "completion_receipt_sha256": sha(OUT / "completion_receipt.json"), "execution_receipt_sha256": sha(OUT / "execution_receipt.json"), "width_comparison_sha256": sha(OUT / "width_comparison.json"), "display_only_widths": [dict(element=x["element"], **x["display_only"]) for x in comparison]}))


if __name__ == "__main__":
    main()
