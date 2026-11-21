"""Independent pre-matrix admission of the four observed first-wave hull replays."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "PLAN.md").is_file())
sys.path.insert(0, str(ROOT / "scripts"))
import transfer_ambiguity_remaining32_first_wave_reduce as reducer

R = reducer
SOURCE = "58de9f5e96b43005b0ba83dd0f4d8cd1b35f0690eef58290a75ef5abce692d92"
PROTOCOL = "4b9ecc9794f84e00d62f1fbc090b3a96739af044a2a426194d8fba07967fe487"
PREFLIGHT = "2f1c69917e32d25d22a28577f80dc6a2c8c7fab025d93967e73fcf48cffeb79c"
SOURCE_REVIEW = "7040ba3e30a276df43e41a61378b86f3223633b47fd3b7c5dce810e5abf58127"
HULL_REVIEW = "987577c3b4eaa0c7dce0655a04a88ece28bf0e891fbc84fb573acd26db95edec"
RAW = "bb9f059699e6d98fe626491f252edddb1664fbf83d79cb74258084ffd399c880"
SPD = {
    "observed_os_exit_sha256":"fd8ac1cb3b920b57096649544bba105442404fe19576aae99c53f849fe15a8d5",
    "replay_phase":{
        "source_sha256":"28c511fc6087d130f353f67c3cb655aa729031b557ae0bd84c9714b632435a13",
        "protocol_sha256":"a740f6a3305421de51a8fa7ce1e5da94bf308ae35ab7226c0ef3336e1beed8c7",
        "preflight_sha256":"c8150950f934719f6212b579bcac066573c94b09ac5943cd7f4c9b5d489339d7",
    },
}


def require(state, message):
    if not state:
        raise ValueError(message)


def ep(pair):
    require(type(pair) is list and len(pair) == 2 and all(type(x) is str for x in pair), "explicit exact endpoints")
    x, y = map(Fraction, pair)
    require(x <= y, "reversed exact endpoints")
    return x, y


def equal_index(x, y):
    if type(y) is int:
        return type(x) is int and x == y
    return type(x) is list and len(x) == 2 and all(type(v) is int for v in x) and tuple(x) == y


def channel(native, original, ids):
    rows = native["rows"]
    require(type(rows) is list and len(rows) == len(original) == len(ids), "complete signed endpoint channel")
    misses = []
    for row, raw, index in zip(rows, original, ids, strict=True):
        require(equal_index(row["index"], index), "typed local index mismatch")
        ind = ep(row["independent_exact"])
        cast = ep(raw)
        require(ep(row["producer_exact"]) == cast, "original producer interval altered")
        require(ind[0] <= cast[1] and cast[0] <= ind[1] and row["disjoint"] is False,
                "true disjointness or fabricated disposition")
        require(ep(row["hull_exact"]) == (min(ind[0], cast[0]), max(ind[1], cast[1])),
                "exact complete endpoint hull mismatch")
        contains = cast[0] <= ind[0] <= ind[1] <= cast[1]
        require(row["producer_contains_independent"] is contains, "nesting diagnostic altered")
        if not contains:
            misses.append(row["index"])
    reported = native["producer_containment_misses"]
    require(type(reported) is list and len(reported) == len(misses) and
            all(equal_index(x, y if type(y) is int else tuple(y)) for x,y in zip(reported, misses, strict=True))
            and native["disjoint_indices"] == [], "missing literal containment diagnostics or disjoint HOLD")
    return len(misses), len(rows)


def inspect_case(case, auth, native, raw_parts, joined_spd):
    folder, receipt, phase, receipt_pin = native
    opts, d, cap, proofs = R.W.reviewed_context(auth, case)
    plan = R.read(R.S.PLAN / case / "plan.json", opts["plan_sha256"])
    require(plan["bindings"]["cached_derivative_case_sha256"] and plan["bindings"]["geometry_replay_case_sha256"],
            "same target requires independent derivative and geometry")
    require(len(plan["signed_cached_RF_center_rows_exact"]) == 611
            and all(ep(x) for x in plan["signed_cached_RF_center_rows_exact"]), "all exact signed cached anchor rows")
    matches = Counter()
    scalar_rows = Counter()
    links = {}
    for selected in receipt["case_components"]:
        k, slot = selected["kind"], selected["slot"]
        proof = proofs[k,slot]
        path, part = raw_parts[case,k]
        raw = part["components"][slot]
        filename = f"{k}_{slot:03d}.json"
        leafpath, rawpath = folder/filename, path/filename
        require(selected["native_file"] == R.rel(leafpath) and selected["packet_file"] == R.rel(rawpath)
                and selected["packet_sha256"] == raw["sha256"], "mixed, reordered or substituted proof leaf")
        leaf = R.read(leafpath, selected["native_sha256"])
        source = R.read(rawpath, raw["sha256"])
        links[R.rel(leafpath)] = selected["native_sha256"]
        R.W.validate_packet(source, proof, opts["target"], slot, k, cap)
        require(leaf["phase"] == phase and leaf["packet_sha256"] == raw["sha256"]
                and leaf["selector_leaf_content_sha256"] == R.W.B.digest(proof)
                and leaf["status"] == R.LEAF_PASS and leaf["kind"] == k
                and type(leaf["slot"]) is int and leaf["slot"] == slot
                and leaf["record_id"] == opts["target"]["record_id"]
                and source["trial_box_sha256"] == opts["target"]["trial_box_sha256"]
                and source["selector_plan_sha256"] == opts["plan_sha256"]
                and source["worker_identity"] == {x:part[x] for x in ("pid","supervisor_pid","stage","case_id","part","bindings")}
                and leaf["only_independently_valid_enclosures_used_for_hull_truth"] is True
                and leaf["overlap_not_used_as_proof"] is True and leaf["root_asserted"] is False,
                "same independently proved signed local scopes and matching producer leaf")
        g = proof["global_indices"]
        require(type(g) is list and all(type(x) is int and 0 <= x < 611 for x in g)
                and leaf["global_indices"] == selected["global_indices"] == g
                and all(type(i) is int for i in leaf["global_indices"] + selected["global_indices"])
                and (leaf["point_force_scope"],leaf["whole_H_scope"]) == (R.S.POINT,R.S.BOX)
                and Fraction(leaf["radius_exact"]) == R.S.RADIUS
                and leaf["checked_cache_sites"] == (64 if k == "bulk" else 12),
                "exact free projection, separate scopes and complete sites")
        F = source["point_force"]
        H = source["whole_box_hessian"]
        upper = list(combinations_with_replacement(range(len(g)), 2))
        originals = (F["F_I_minus_F_cache_intervals_exact"], [row[2:] for row in H["H_I_minus_H_cache_upper_exact"]],
                     F["Qstar_minus_cache_correction_exact"])
        channels = leaf["signed_channels"]
        require(set(channels) == set(R.CHANNELS), "complete separate three-channel proof")
        for key, signed, indices in zip(R.CHANNELS, originals, (list(range(len(g))), upper, list(range(len(g)))), strict=True):
            miss, count = channel(channels[key], signed, indices)
            matches[key] += miss
            scalar_rows[key] += count
        tails = leaf["independent_point_and_whole_tails_exact"]
        force_tail = tails["F"]
        hessian_tail = tails["H"]
        c = leaf["independent_same_Y_whole_H_correction_exact"]
        require(len(force_tail) == len(g) and all(type(x) is str and 0 <= Fraction(x) <= cap[0] for x in force_tail)
                and leaf["producer_force_tails_exact"] == F["force_tail_upper_exact"]
                and len(hessian_tail) == len(c) == len(upper), "exact independently sufficient tail budget and local derivative coverage")
        for h, value, pair in zip(hessian_tail, c, upper, strict=True):
            require(len(h) == 3 and len(value) == 4 and equal_index(h[:2],pair) and equal_index(value[:2],pair)
                    and type(h[2]) is str and 0 <= Fraction(h[2]) <= cap[1] and ep(value[2:]),
                    "same typed upper Hessian tail and correction")
    require(len(links) == 168 and dict(matches) == receipt["producer_containment_miss_counts"],
            "all 168 native endpoint leaves with truthful counts")
    rnorm = Fraction(plan["sizing"]["R_inf_exact"])
    established = R.min_premise(joined_spd, {"opts":opts,"d":d,"plan":plan}, rnorm)
    require(established is not None, "same R SPD supported for this first wave")
    return {"case_id":case,"native_receipt_sha256":receipt_pin,"leaf_sha256":links,"signed_scalar_counts":dict(scalar_rows),
            "producer_containment_miss_counts":dict(matches),"disjoint_signed_scalar_counts":dict.fromkeys(R.CHANNELS,0),
            "same_R_SPD_bindings":established}


def cold():
    a, w = R.static_phase()
    p = R.read(R.PREP/"protocol.json",PROTOCOL)
    f = R.read(R.PREP/"preflight.json",PREFLIGHT)
    review = R.read(R.REVIEW_FILE,SOURCE_REVIEW)
    R.hull_static(HULL_REVIEW)
    bind = {"source_sha256":SOURCE,"protocol_sha256":PROTOCOL,"preflight_sha256":PREFLIGHT,
            "first_wave_worker_admission_sha256":R.W_ADMIT,
            "joined_bridge_input_admission_sha256":a["joined_bridge_input_admission_sha256"]}
    require(p["source_sha256"] == f["source_sha256"] == R.sha(R.SOURCE) == SOURCE
            and f["protocol_sha256"] == PROTOCOL and p["worker_runtime_binding"] == w
            and review["status"] == R.REVIEW and R.V.yes(review["checks"])
            and all(review[k] == v for k,v in bind.items()), "same approved first-wave reducer and native proof method")
    return a,w,review,bind


def run(os_pins):
    a,w,review,bind = cold()
    raw,parts = R.observed_raw(RAW,a,w)
    refs = [{"case_id":case,"observed_os_exit_sha256":h} for case,h in zip(R.W.FIRST,os_pins,strict=True)]
    base = {"status":R.ADMISSION,**bind,"reducer_review_sha256":SOURCE_REVIEW,
            "root_producer_admission_sha256":RAW,"hull_source_review_sha256":HULL_REVIEW,"hull_case_observations":refs,
            "all_four_cases_required_before_arithmetic":True,"new_spatial_integrals_authorized":False,"optional_spd":SPD}
    # All four root-observed true exits and phase/target metadata before any leaf.
    native = R.open_all_hull_metadata(a,base,raw)
    spd = R.optional_spd(SPD)
    cases = [inspect_case(case,a,native[case],parts,spd) for case in R.W.FIRST]
    require(len({leaf for item in cases for leaf in item["leaf_sha256"]}) == 672, "complete independent leaf evidence")
    doc = {**base,"checks":{
               "all_eight_raw_workers_and_four_independent_hull_replays_genuine_successful_exits":True,
               "all_672_raw_and_672_independently_proved_local_packets_checked":True,
               "all_three_signed_channels_exact_rational_full_hull_endpoints_replayed_no_disjoint":True,
               "all_nonnested_diagnostics_retained_with_matching_analytic_tail_and_selector":True,
               "all_four_same_R_SPD_premises_independently_observed":True,
               "independent_global_endpoint_replay_still_required":True},
           "admission_source_sha256":R.sha(Path(__file__)), "case_readback":cases, "matrix_arithmetic_jobs":0,
           "new_spatial_integrals":0,"root_or_minimum_claimed":False}
    R.check_gates(review,doc,bind,SOURCE_REVIEW)
    print({"admission_sha256":R.save(R.ADMIT_FILE,doc),"independent_local_leaves":672,"roots_claimed":0},flush=True)


def fixtures():
    cold()
    r = {"index":[0,1],"independent_exact":["-3/7","2/3"],"producer_exact":["-1/5","9/7"],
         "hull_exact":["-3/7","9/7"],"producer_contains_independent":False,"disjoint":False}
    wrap = {"rows":[r],"producer_containment_misses":[[0,1]],"disjoint_indices":[]}
    require(channel(wrap,[r["producer_exact"]],[(0,1)]) == (1,1), "safe nonnested exact hull")
    variants = ({"index":[False,1]}, {"hull_exact":["-1/5","2/3"]}, {"producer_contains_independent":True},
                {"independent_exact":["-1","-1/2"], "hull_exact":["-1","9/7"],"disjoint":True})
    for bad in variants:
        require(R.V.A.reject(lambda:channel({**wrap,"rows":[{**r,**bad}]},[r["producer_exact"]],[(0,1)])),
                "independent exact guard accepted corrupted data")
    print(R.save(HERE/"preflight.json",{"status":"PASS","source_sha256":R.sha(Path(__file__)),"checks":{
        "cold_reducer_and_hull_sources_match":True,"exact_nonnested_hull_accepted_after_valid_independent_enclosure":True,
        "typed_index_shrunk_hull_false_nesting_and_real_disjoint_rejected":True,"no_actual_payloads_opened":True},
        "same_root_producer_admission_sha256":RAW,"root_claimed":False}),flush=True)


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("phase", choices=("prepare","run"))
    args.add_argument("--observed-os-sha256",nargs=4)
    params = args.parse_args()
    if params.phase == "prepare":
        fixtures()
    else:
        require(params.observed_os_sha256 is not None, "four external true-observation pins required")
        run(params.observed_os_sha256)
