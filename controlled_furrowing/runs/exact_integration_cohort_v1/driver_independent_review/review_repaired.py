"""Combine repaired L driver checks and issue the narrowly scoped gate."""

import inspect
import json
from fractions import Fraction as Q
from pathlib import Path
from unittest.mock import patch

import check_driver_gate as gate_check
import check_selector_smoke as surface_check
import review_prelaunch as exact_check

P = surface_check.P
HERE = Path(__file__).resolve().parent
OUT = exact_check.PILOT / "original_cube_v1"
SOURCE = "b73fb8bbaecdb7ed17e8be986c9ad74acc3abe1ed4dabe1843bc9cbed11f8213"
PROTOCOL = "79118e68dafcc869d12a6b5ee399e4eb342ef77433d24a24fd8a57a3f0021478"
PREFLIGHT = "9f6e8dbb7699b8eb140b57d465933d878356dda4fe9dfa03a33f4dad7d204f4e"


def failure_receipt_check():
    captured = {}

    def retain(path, value):
        captured[Path(path).name] = value

    with (patch.object(P, "gate"), patch.object(P, "write", side_effect=retain),
          patch.object(P, "domains_work", side_effect=RuntimeError("synthetic failure"))):
        try:
            P.domains(exact_check.RECORD)
        except RuntimeError as error:
            assert str(error) == "synthetic failure"
        else:
            raise AssertionError("Driver did not propagate the synthetic failure")
    failure = captured["domains_failure.json"]
    terminal = captured["domains_terminal.json"]
    return (set(captured) == {"domains_started.json", "domains_failure.json", "domains_terminal.json"} and
            failure["error"] == "synthetic failure" and failure["source_sha256"] == SOURCE and
            failure["protocol_sha256"] == PROTOCOL and terminal["status"] == "FAILED" and
            terminal["exit_code"] == 1 and terminal["pid"] == captured["domains_started.json"]["pid"])


def exact_rounding_check(protocol):
    context = P.C.load_context(exact_check.RECORD)
    box = P.read(P.ROOT / protocol["coordinate_box_file"])
    nodes = P.S.common_nodes(context.target, box)
    rule = context.target.auxiliary_rule()
    results = []
    for face, edge in surface_check.EDGES:
        fields = context.target.surface(face, edge, nodes=nodes)
        domain = P.S.uniform_G_bounds(fields, max_depth=8)
        low, high = map(Q, domain["G_bounds_exact"])
        midpoint = Q(domain["q0_exact"])
        raw = Q(domain["rho_exact"])
        before = dict(domain)
        order, attempts = P.surface_selector(context, fields, domain, rule, face, edge, protocol)
        rho = Q(domain["rho_exact"])
        assert Q(domain["raw_rho_exact"]) == raw == (high - low) / (high + low)
        assert rho * 2**48 == -(-raw * 2**48 // 1) and raw <= rho < 1
        assert midpoint == (low + high) / 2 and Q(domain["q0_exact"]) == midpoint
        assert all(value == domain[key] for key, value in before.items() if key != "rho_exact")
        assert order == attempts[-1]["degree"]
        results.append({"face": face, "edge": edge, "degree": order,
                        "raw_rho_display": float(raw), "outward_rho_exact": str(rho)})
    return results


def write_new(path, value):
    with path.open("x", encoding="utf8") as output:
        json.dump(value, output, indent=2)
        output.write("\n")


def main():
    for name in ("domains_started.json", "domains.json", "domains_terminal.json",
                 "domains_prelaunch_independent_review.json"):
        assert not (OUT / name).exists(), f"Preserve existing {name}"
    exact = exact_check.review(OUT, SOURCE, PROTOCOL, PREFLIGHT)
    protocol = P.read(OUT / "protocol.json")
    guards = gate_check.verify(OUT)
    smoke = surface_check.smoke(OUT)
    rounded = exact_rounding_check(protocol)
    domain_source = inspect.getsource(P.domains_work)
    checks = {
        **exact["checks"], **guards["checks"],
        "reviewed_body_kernel_and_both_transitive_geometry_sources_bound":
            protocol["source_and_input_hashes"]["scripts/root_body_geometry.py"] ==
            "18013c6c50c5f679eb6c4abc04e3cf9ce447ee18d0bf1efbf6dcfbd6806f49c2" and
            all(o["expected_sha256"] == o["checked_current_sha256"] == o["driver_protocol_sha256"]
                for o in exact["geometry_dependency_bindings"].values()),
        "old_L_cached_bound_provider_direct_record_review":
            protocol["cached_bound_provider"].endswith("/optoGEF2_small_col_L.json"),
        "actual_L_repaired_surface_selectors_and_36_unweighted_site_domains_pass":
            smoke["status"] == "PASS" and len(smoke["cases"]) == 3 and
            sum(row["sites_checked"] for row in smoke["cases"]) == 36,
        "each_smoke_rho_rounds_outward_and_keeps_G_and_center_unchanged": len(rounded) == 3,
        "domain_work_loads_requested_L_exact_box_and_fresh_all_element_body":
            all(t in domain_source for t in ("C.load_context(record)", "G.certify_box(",
                                             "S.common_nodes(context.target,box)", "range(144)")),
        "domain_work_selects_all_closed_edges_and_288_coefficientwise_cache_lines":
            all(t in domain_source for t in ("for face in(0,1)", "for edge in range(12)",
                                             "for iq in range(12)", "S.local_site_domain(")),
        "domain_worker_has_no_actual_force_moment_or_weighted_derivative_calls":
            not any(t in domain_source for t in ("bulk_component(", "surface_component(",
                                                 "cache_edge_correction(", "arb_mat(", ".defect(")),
        "synthetic_domain_failure_writes_both_failure_and_terminal_receipts": failure_receipt_check(),
        "no_full_domain_or_force_job_launched_by_this_review":
            not any((OUT / name).exists() for name in ("domains_started.json", "body_geometry.json",
                                                       "bulk_started.json", "surface_started.json")),
    }
    assert all(v is True for v in checks.values()), [k for k, v in checks.items() if not v]
    support_path = HERE / "repaired_L_prelaunch_checks.json"
    write_new(support_path, {"exact_box_and_budget_review": exact, "approval_guard_checks": guards,
                             "repaired_bounded_surface_smoke": smoke, "exact_outward_rounding_checks": rounded})
    receipt = {
        "status": "PASS", "checks": checks, "source_sha256": SOURCE,
        "protocol_sha256": PROTOCOL, "preflight_sha256": PREFLIGHT, "record_id": exact_check.RECORD,
        "support_file": support_path.relative_to(exact_check.ROOT).as_posix(),
        "support_sha256": P.sha(support_path),
        "review_source_hashes": {Path(m.__file__).name: P.sha(m.__file__)
                                 for m in (exact_check, gate_check, surface_check)}
                                | {Path(__file__).name: P.sha(__file__)},
        "scope": "Prelaunch for the fixed L original cube fresh body and local domain/selector phase only. "
                 "No full body job, force/Hessian integral, moment-defect or weighted cache sum, root proof, "
                 "or other cohort record is approved by this review.",
        "other_record_followup": "Equality to a radius upper enclosure does not generalize to non-dyadic "
                                 "exact selected radii. Future cohort preparation must check enclosure "
                                 "containment while retaining the exact authoritative selected radius.",
    }
    gate_path = OUT / "domains_prelaunch_independent_review.json"
    write_new(gate_path, receipt)
    P.gate(OUT, "domains")
    print(json.dumps({"status": "PASS", "checks": len(checks), "receipt_sha256": P.sha(gate_path),
                      "support_sha256": P.sha(support_path), "rounded": rounded}, indent=2))


if __name__ == "__main__":
    main()
