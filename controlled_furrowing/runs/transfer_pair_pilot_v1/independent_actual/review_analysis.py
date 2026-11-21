"""Check downstream comparison data/figure bindings against the independent four-field replay."""

import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path

HERE=Path(__file__).resolve().parent
PILOT=HERE.parent
ROOT=HERE.parents[2]
ANALYSIS=PILOT/"analysis"
FILES={
    "analyze.py":"fdec9194a9c666dba665f180485015cfc3e3c4f2b69c80962b91c4af7b8282aa",
    "protocol.json":"3e260b0d1acc069d22eade9c187b0ddf71dda45687aa18daebe38c77e4195349",
    "receipt.json":"42be29376602c5ae9c08d9d580324928600a35b7201cb452ba501a936257e234",
    "current_radius_profiles.png":"e17288fa5b1cd0824778bac5de522fc0a372ad70f02f8eb446200e003dd4160e",
    "shape_contrast_and_diameter_effect.png":"59746706b826af8fc8d6469849846664ba1b0b43249ca4bf83f57866b790900b",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def bounds(record):
    return tuple(map(Q,record["exact"]))


def main():
    checks={}
    def ck(name, condition):
        assert condition,name
        checks[name]=True
    ck("fixed_analyzer_protocol_receipt_and_both_figure_digests",all(sha(ANALYSIS/k)==v for k,v in FILES.items()))
    protocol=read(ANALYSIS/"protocol.json")
    analysis=read(ANALYSIS/"receipt.json")
    replay_path=HERE/"actual_review.json"
    replay=read(replay_path)
    binding=replay["all_replayed_hashes"]
    ck("producer_current_AR_receipts_and_analysis_source_protocol_bound",
       replay["status"]=="PASS_ACTUAL_FOUR_SAVED_FIELD_REPLAY" and
       analysis["results_sha256"]==protocol["result_sha256"]==replay["producer_results_sha256"]==sha(PILOT/"results.json") and
       analysis["source_sha256"]==protocol["analyzer_sha256"]==FILES["analyze.py"] and
       analysis["protocol_sha256"]==FILES["protocol.json"] and
       analysis["supervisor_receipt_sha256"]==sha(PILOT/"job.json"))
    ck("analysis_carries_same_frozen_inputs_and_sources",
       all(binding[p]==h==sha(ROOT/p) for p,h in
           {**analysis["checked_current_source_hashes"],**analysis["checked_current_input_hashes"]}.items()))
    ck("all_four_Q_exact_intervals_match_independent_Lagrange_inversion",
       all(bounds(analysis["signed_Q_um"][k])==bounds(v) for k,v in replay["Q_um"].items()))
    ck("both_diameter_changes_and_large_gap_exact",
       all(bounds(analysis["diameter_effect_large_minus_small_Q_um"][k])==bounds(v)
           for k,v in replay["large_minus_small_Q_um"].items()) and
       bounds(analysis["joint_minus_baseline_Q_um"]["large"])==bounds(replay["joint_minus_baseline_large_Q_um"]) and
       analysis["large_pair_disjoint_stored_curve_intervals"] is True)
    ck("independent_difference_of_differences_exact",
       bounds(analysis["difference_of_differences_joint_minus_baseline_effect_um"])==
       bounds(replay["joint_minus_baseline_diameter_change_um"]))
    numerical_review=[]
    for name,case in analysis["cases"].items():
        own=replay["cases"][name]
        for p,h in ((case["profiles"],case["profiles_sha256"]),
                    (case["worker_receipt"],case["worker_receipt_sha256"]),
                    (case["terminal_sidecar"],case["terminal_sidecar_sha256"])):
            assert binding[p]==h==sha(ROOT/p)
        assert case["samples"]==own[-1]["profile_sample_count"]
        assert case["plot_connects_retained_samples"] is True and case["root_error_enclosure"] is None
        assert case["minimum_graph_thickness_um_lower_exact"]==own[-1]["current_thickness_lower_um_exact"]
        assert len(case["checkpoints"])==len(own)==4
        for row,local in zip(case["checkpoints"],own):
            f8,f12=row["Q8_normalized_force_linf"],row["Q12_normalized_force_linf"]
            assert abs(f8-local["Q8_force_linf"]) <= 3e-14 and abs(f12-local["Q12_force_linf"]) <= 3e-14
            assert max(f8,f12)<analysis["normalized_force_limit"]
            assert binding[row["field"]]==row["field_sha256"] and binding[row["verification"]]==row["verification_sha256"]
            numerical_review.append((f8,f12))
    ck("all_four_OS_receipts_profile_samples_and_16_archived_residuals_match",len(numerical_review)==16)
    ck("plot_provenance_uses_predeclared_both_PNGs",
       set(analysis["visualizations"])=={x for x in FILES if x.endswith(".png")} and
       all(str((ROOT/row["path"]).resolve())==str((ANALYSIS/name).resolve()) and
           row["sha256"]==FILES[name] for name,row in analysis["visualizations"].items()))
    ck("scope_leaves_equilibrium_full_profile_supremum_and_source_claims_open",
       analysis["root_position_error_enclosure_um"] is None and analysis["full_profile_supremum_enclosure_um"] is None and
       not any(analysis[k] for k in ("root_validated","source_compatible_witness","larger_source_outcomes_read_or_scored","novelty_established")) and
       protocol["heldout_outcomes_read"] is False and "sample" in protocol["profile_plot"].lower())
    d={"status":"PASS_TRANSFER_ANALYSIS_AND_FIGURE_PROVENANCE_READBACK","checks":checks,
       "reviewer_source_sha256":sha(Path(__file__)),"actual_review_sha256":sha(replay_path),
       "analysis_hashes":FILES,"max_sampled_Q8_normalized_force_linf":max(x[0] for x in numerical_review),
       "max_sampled_Q12_normalized_force_linf":max(x[1] for x in numerical_review),
       "scope":"Static data binding and exact interval comparison. The figures connect stored samples. No "
               "full-profile supremum, equilibrium-position, integration, continuum or measurement uncertainty bound."}
    path=HERE/"analysis_readback.json"
    with path.open("x",encoding="utf8") as file:
        json.dump(d,file,indent=2)
        file.write("\n")
    print(json.dumps({"status":d["status"],"checks":len(checks),"sha256":sha(path),
                      "max_Q8":d["max_sampled_Q8_normalized_force_linf"],
                      "max_Q12":d["max_sampled_Q12_normalized_force_linf"]},indent=2))


if __name__=="__main__":
    main()
