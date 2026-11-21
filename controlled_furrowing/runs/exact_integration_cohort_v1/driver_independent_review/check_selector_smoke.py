"""Smoke-check actual L surface domains and scalar selectors without integrations."""

from fractions import Fraction as Q
from pathlib import Path
import argparse
import json
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import cohort_exact_pilot as P

RECORD = "optoGEF2_small_col_L"
EDGES = ((1, 11), (1, 1), (0, 0))


def smoke(out):
    protocol_file = out / "protocol.json"
    protocol = P.read(protocol_file)
    P.assert_inputs(protocol)
    P.precision.prec = 640
    P.core.A.cache_clear()
    context = P.C.load_context(RECORD)
    box = P.read(ROOT / protocol["coordinate_box_file"])
    nodes = P.S.common_nodes(context.target, box)
    rule = context.target.auxiliary_rule()
    report = {"record_id": RECORD, "protocol_sha256": P.sha(protocol_file), "cases": []}
    for face, edge in EDGES:
        started = time.perf_counter()
        fields = context.target.surface(face, edge, nodes=nodes)
        domain = P.S.uniform_G_bounds(fields, max_depth=8)
        rho = Q(domain["rho_exact"])
        row = {
            "face": face, "edge": edge, "raw_rho_upper_exact": str(rho),
            "raw_rho_numerator_bits": rho.numerator.bit_length(),
            "raw_rho_denominator_bits": rho.denominator.bit_length(),
        }
        try:
            order, attempts = P.surface_selector(context, fields, domain, rule, face, edge, protocol)
            lower = attempts[-1]
            assert order == lower["degree"] and 2 <= order <= protocol["max_surface_order"]
            assert Q(lower["force_tail_upper_exact"]) <= Q(protocol["force_tail_budget_exact"])
            assert Q(lower["H_entry_tail_upper_exact"]) <= Q(protocol["H_entry_tail_budget_exact"])
            sites = [P.S.local_site_domain(context.target, box, face, edge, iq) for iq in range(12)]
            assert all(site["global_indices"] == fields["global_indices"] and
                       site["selected_local"] == fields["free_local"] and
                       Q(site["convex_line_exact"][0][0]) > 0 and
                       Q(site["tangent_norm_squared_lower_exact"]) > 0 and
                       Q(site["auxiliary_normalized_weight_exact"]) == rule["surface_weights"][face][edge][iq]
                       for iq, site in enumerate(sites))
            row.update(status="PASS", selected_order=order, attempt_count=len(attempts),
                       final_attempt=lower, sites_checked=len(sites))
        except Exception as failure:
            row.update(status="FAIL", error=repr(failure), traceback=traceback.format_exc())
        row["seconds"] = time.perf_counter() - started
        report["cases"].append(row)
    report["status"] = "PASS" if all(row["status"] == "PASS" for row in report["cases"]) else "FAIL"
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run_directory")
    args = parser.parse_args()
    out = P.BASE / RECORD / args.run_directory
    result = smoke(out)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "PASS" else 1)
