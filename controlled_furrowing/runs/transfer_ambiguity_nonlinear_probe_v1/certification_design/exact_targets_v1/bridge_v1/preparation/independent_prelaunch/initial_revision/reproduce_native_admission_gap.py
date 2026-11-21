import sys
from pathlib import Path
from unittest.mock import patch

ROOT = Path.cwd()
sys.path.insert(0,str(ROOT/"scripts"))
import transfer_ambiguity_exact_bridge as B

assert B.sha(B.SOURCE) == "fd79fea21ed8f33863e0db9e5ff9ba47f20ead4def0295db13e92bba8ceabeca"
A = B.DX.ACTUAL
N = B.DX.OUT / "independent_actual/FOREIGN.json"
WRAP = ROOT/"mock_wrapper.json"
native = {"status":"PASS_INDEPENDENT_UNRELATED","checks":{"positive":True},"reviewed_cases":[]}
candidates = []
u = {"status":"PASS_REVIEWED_LARGER_CACHED_DERIVATIVE_BRIDGE_INPUTS","checks":{"true":True},
     "radius_exact":str(B.RADIUS),"native_independent_review_file":B.EX.relative(N),
     "native_independent_review_sha256":"native","derivative_results_sha256":"top","derivative_terminal_sha256":"end",
     "reviewed_cases":[]}
bs = {"source_sha256":B.PINS[B.DERIV_SOURCE]}
top = {"status":"COMPLETE_FOUR_CACHED_DERIVATIVE_ATTEMPTS_UNREVIEWED","bindings":bs,"case_results":[]}
end = {"status":"FINISHED","bindings":bs,"results_sha256":"top","pid":7,"workers":[]}
files = {WRAP:u,N:native,A/"results.json":top,A/"terminal.json":end}
for i,p in enumerate(B.EX.PILOTS):
    h = f"producer-{i}"
    c = {"pilot_id":p,"box_sha256":f"box-{i}","cache_packet_sha256":f"cache-{i}",
         "R_inf_exact":"10","prior_optional_tail_policy":{"max_force_incidence":6,
             "max_H_entry_row_incidence":72,"tail_force_budget_exact":"1/1000","tail_H_entry_budget_exact":"1/1000"}}
    candidates.append(c)
    top["case_results"].append({"pilot_id":p,"results_sha256":h,"ZQ_upper_exact":"1/5"})
    u["reviewed_cases"].append({**c,"producer_result_sha256":h,"reviewed_ZQ_upper_exact":"1/7"})
    identity = {"pilot_id":p,"pid":100+i,"supervisor_pid":7,"bindings":bs}
    end["workers"].append({**identity,"observed_exit_code":0,"terminal_sha256":"child","os_exit_receipt_sha256":"os"})
    files[A/p/"results.json"] = {"case":c,"whole_box_cache_radius_exact":str(B.RADIUS),
                                "whole_box_controls_only_cached_quadrature":True,"ZQ_upper_exact":"1/5"}
    files[A/p/"terminal.json"] = {**identity,"status":"FINISHED","results_sha256":h}
    files[A/p/"observed_os_exit.json"] = {"observed_exit_code":0}
with patch.object(B,"pinned",side_effect=lambda path,digest: files[path]):
    accepted = B.derivative_inputs(WRAP,"wrapper",candidates)
assert len(accepted)==4
print("UNRELATED_NATIVE_REVIEW_ACCEPTED",[(k,v["ZQ_upper_exact"]) for k,v in accepted.items()])
