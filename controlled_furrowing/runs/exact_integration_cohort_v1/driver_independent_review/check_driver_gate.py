"""Exercise gate rejection rules without writing a gate or launching the driver."""

import argparse
import json
import sys
from pathlib import Path
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import cohort_exact_pilot as P


def verify(out):
    source = P.sha(P.SOURCE)
    protocol = P.sha(out / "protocol.json")
    preflight = P.sha(out / "preflight.json")
    accepted_example = {
        "status": "PASS", "checks": {"substantive_check": True},
        "source_sha256": source, "protocol_sha256": protocol, "preflight_sha256": preflight,
    }
    rejected_examples = {
        "missing_pass": {"status": "FAIL"},
        "empty_checks": {"checks": {}},
        "truthy_number_check": {"checks": {"substantive_check": 1}},
        "truthy_string_check": {"checks": {"substantive_check": "true"}},
        "explicit_false_check": {"checks": {"substantive_check": False}},
        "changed_source_hash": {"source_sha256": "0" * 64},
        "changed_protocol_hash": {"protocol_sha256": "0" * 64},
        "changed_preflight_hash": {"preflight_sha256": "0" * 64},
    }
    reader = P.read
    file_path = out / "domains_prelaunch_independent_review.json"

    def call_gate(input_receipt):
        def read(path):
            return input_receipt if Path(path) == file_path else reader(path)
        with patch.object(P, "read", side_effect=read):
            P.gate(out, "domains")

    call_gate(accepted_example)
    checks = {"complete_match_accepted": True}
    for name, override in rejected_examples.items():
        try:
            call_gate({**accepted_example, **override})
        except ValueError:
            checks[name + "_rejected"] = True
        else:
            raise AssertionError(name + " unexpectedly passed")
    return {"status": "PASS", "checks": checks, "source_sha256": source,
            "protocol_sha256": protocol, "preflight_sha256": preflight}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run_directory")
    args = parser.parse_args()
    path = P.BASE / "optoGEF2_small_col_L" / args.run_directory
    print(json.dumps(verify(path), indent=2))
