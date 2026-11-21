"""Freeze draft/figures and their receipt bindings; no publication or science."""
from pathlib import Path
import datetime
import hashlib
import json
import re

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[3]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text(encoding="utf8"))


def write(p, obj):
    if p.exists():
        raise ValueError("preserve frozen authoring receipt")
    p.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf8")


def main():
    draft = OUT / "EXACT-INTEGRAL-ROOT.md"
    primary = ROOT / "evidence/exact_integral_root_B_safe_independent_review.json"
    direct = ROOT / "evidence/exact_uniform_direct_root_B_crosscheck.json"
    result = ROOT / "runs/exact_integration_bridge_v1/integral_root_B_safe_v2/results.json"
    doc = draft.read_text(encoding="utf8")
    a, b, plots = read(primary), read(direct), read(OUT / "figure_generation.json")
    checks = {
        "final_primary1500_review_exact_hash": sha(primary) == "4500e5b34039eafe2e8af742cad739b1133ab414521fd6ed604cf4787d80d561",
        "primary1500_nonempty_literal_PASS": a["status"] == "PASS" and len(a["checks"]) == 1500 and all(v is True for v in a["checks"].values()),
        "primary_review_binds_actual_result": a["certificate_sha256"] == sha(result),
        "direct_crosscheck_exact_hash": sha(direct) == "f28b458b66d718c0fa7b4ef35a29ee463fc59612d887dee52f7e661d0d41ecfe",
        "direct_crosscheck_nonempty_PASS": b["status"] == "PASS_DIRECT_UNIFORM_FORCE_ROOT_CROSSCHECK" and bool(b["checks"]) and all(v is True for v in b["checks"].values()),
        "sixteen_exact_figure_arithmetic_checks": len(plots["checks"]) == 16 and all(v is True for v in plots["checks"].values()),
        "word_budget_under2000": len(doc.split()) <= 2000,
        "three_targets_explicit": all(x in doc for x in ("E_Q", "E_I", "continuum energy")),
        "independent_and_coupled_domains_distinguished": "independent domain" in doc and "coupled cubes" in doc,
        "parameter_variation_not_dropped": "V:=bs+A" in doc,
        "SPDH0_homotopy_explicit": "H_t=(1-t)H_0+tH" in doc,
        "exclusion_not_from_point_escape_alone": "existence plus common-cube uniqueness" in doc,
        "second_route_shared_kernels_disclosed": "sharing the validated uniform force/Hessian kernels" in doc,
        "predictor_not_root_legend_and_figure_scope": "diamonds mark predictors, not equilibria" in doc,
        "renderable_determinant_definition": "\\(J=\\det F\\)" in doc,
        "fixed_error_sentence_clarified": "Its error relative to the preconditioned integration defect is bounded by m, defined below." in doc,
        "no_private_windows_path_in_draft": "C:\\Users" not in doc and "C:/Users" not in doc,
        "final_review_pending_prose_removed": "final independent review is pending" not in doc,
    }
    local_links = []
    # Resolve against the FINAL publication destination, not this run folder.
    for target in re.findall(r"!?\[[^\]]*\]\(([^)]+)\)", doc):
        if target.startswith(("https://", "http://", "#")):
            continue
        resolved = (ROOT / "docs" / target).resolve()
        checks["final_docs_link_" + target] = resolved.is_file()
        local_links.append({"link": target, "target_relative_to_project": resolved.relative_to(ROOT.resolve()).as_posix(), "sha256": sha(resolved) if resolved.is_file() else None})
    for name, h in plots["figures"].items():
        checks["figure_exact_hash_" + name] = sha(OUT / name) == h
    if not all(v is True for v in checks.values()):
        raise AssertionError([k for k, v in checks.items() if v is not True])
    native = {
        "status": "PASS_NATIVE_TWO_FIGURE_INSPECTION",
        "checks": {"both_figures_viewed": True, "all_text_labels_fit": True, "strict_threshold_and_both_metrics_readable": True, "coordinate588_Q_I_and_both_cubes_distinguished": True, "predictor_diamonds_explicitly_not_equilibria": True, "figures_original_and_not_observations": True},
        "PNG_hashes": {name: h for name, h in plots["figures"].items() if name.endswith(".png")},
        "scope": "Native images inspected; HTML rendering and publication are parent-owned later gates.",
    }
    write(OUT / "native_visual_inspection.json", native)
    files = [Path(__file__), draft, primary, direct, result, OUT / "plot_root_result.py", OUT / "figure_data_exact.json", OUT / "figure_generation.json", OUT / "native_visual_inspection.json"] + [OUT / name for name in plots["figures"]]
    receipt = {
        "status": "FROZEN_DRAFT_AND_ORIGINAL_VISUALS_READY_FOR_PARENT_PUBLICATION_REVIEW",
        "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "checks": checks,
        "bound_files": {file.relative_to(ROOT).as_posix(): sha(file) for file in files},
        "local_links_resolved_at_final_docs_destination": local_links,
        "draft_words_whitespace_count": len(doc.split()),
        "figure_generation_status_is_historical": "Figures were generated while final scalar review was pending; this receipt binds the now-final1500PASS and direct-route receipts without altering earlier evidence.",
        "public_docs_builders_or_Space_modified": False,
        "science_or_numeric_certificate_replayed": False,
        "scope": "Authoring/figure arithmetic/receipt binding and native image inspection only; full publication/HTML/browser gates remain parent-owned. No new all17/PDE/fullphysicaldomain/dynamics/biology/novelty claim.",
    }
    write(OUT / "authoring_binding.json", receipt)
    print(json.dumps({"status": receipt["status"], "checks": len(checks), "draft_sha256": sha(draft), "binding_sha256": sha(OUT / "authoring_binding.json"), "native_inspection_sha256": sha(OUT / "native_visual_inspection.json"), "words": len(doc.split()), "local_links": len(local_links)}))


if __name__ == "__main__":
    main()
