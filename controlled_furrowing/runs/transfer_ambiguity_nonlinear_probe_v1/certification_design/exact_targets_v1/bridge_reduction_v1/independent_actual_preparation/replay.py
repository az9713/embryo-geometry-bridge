"""Independent rational scatter and directed endpoint contraction for reviewed leaves."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace

from flint import arb, arb_mat, ctx as arb_context, fmpq

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
ROOT = next(p for p in OUT.parents if (p / "PLAN.md").is_file())
EXACT = OUT.parent
CALLER = ROOT / "scripts/transfer_ambiguity_exact_bridge_reduce_driver_v3.py"
GATE = OUT / "caller_preparation_v3/independent_prelaunch/review.json"
SPD = EXACT / "preconditioner_spd_v1/independent_actual_preparation/actual_replay"
GEOMETRY = EXACT / "independent_actual/replay"
PINS = {
    CALLER: "445fef269a49f5af4e4a32cf491995f87172bcda214520ed4574d355fb8405d8",
    GATE: "a0174119fa4eb6a25bff0257482708f22507726d33bad967f7bf848049a077e5",
    SPD / "os_exit_observation.json": "b5a4a1ab78dc56c7a1147bb75b657dc141d3c2ed4a08ea31074bb5c5a2a0a3a0",
    SPD / "receipt.json": "769a65a128aaa79018c050c7bd61d38999166621742946cd0f6a13876854ba83",
    EXACT / "preconditioner_spd_v1/METHOD.md": "0d6ca85eb910c616e6ccac15d67d63b37082d86818c1797b0304a744bf3e5750",
    GEOMETRY / "receipt.json": "884e0bc76a27465c86b62a202cab491d764aedf030623592238e09adee86d77d",
    GEOMETRY / "os_exit_observation.json": "c2be4be8b403e170e20ec54c34f17a69b4aec4299e6f9f1ab7ecb3d4afa337d9",
}
PILOTS = tuple(f"{p}_{q}" for p in ("max_contrast_change", "max_negative_contrast_change")
               for q in ("small_17", "large_03"))
N = 611
RADIUS = Q(1, 2**34)
PREC = 1024
POINT = "EXACT_CENTER_POINT_ONLY"
CUBE = "ENTIRE_CLOSED_FIXED_PARAMETER_CUBE"
FKEY = "F_I_minus_F_cache_intervals_exact"
HKEY = "H_I_minus_H_cache_upper_exact"
ORDER = [(k, i) for k, n in (("bulk", 144), ("surface", 24)) for i in range(n)]
ZERO = (Q(0), Q(0))


def demand(ok, message):
    if not ok:
        raise ValueError(message)


def relative(path):
    return path.resolve().relative_to(ROOT).as_posix()


def literal_checks(d):
    return type(d) is dict and bool(d) and all(v is True for v in d.values())


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path, expected, retained, *, json_payload=True):
    if isinstance(path, str):
        path = ROOT / path
    name = relative(path)
    payload = path.read_bytes()
    actual = sha256(payload).hexdigest()
    demand(actual == expected and (name not in retained or retained[name] == actual), f"changed bytes: {name}")
    retained[name] = actual
    return json.loads(payload) if json_payload else None


def save(path, payload):
    with path.open("x", encoding="utf-8", newline="\n") as destination:
        destination.write(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    return digest(path)


def rational(x):
    demand(type(x) in (str, Q, int), "strict exact rational encoding required")
    return Q(x)


def bounds(row):
    demand(type(row) in (list, tuple) and len(row) == 2, "exact signed interval required")
    a, b = map(rational, row)
    demand(a <= b, "unordered interval")
    return a, b


def plus(a, b):
    return a[0] + b[0], a[1] + b[1]


def contains(parent, child):
    return parent[0] <= child[0] <= child[1] <= parent[1]


def absolute_cap(row):
    return max(abs(row[0]), abs(row[1]))


def absolute_floor(row):
    lo, hi = row
    return max(lo, -hi, Q(0))


def exact_maps(ids, n):
    demand(type(ids) is list and ids and all(type(i) is int and 0 <= i < n for i in ids)
           and len(ids) == len(set(ids)), "invalid literal projected free map")
    return ids


def local_hessian(rows, m):
    expected = list(combinations_with_replacement(range(m), 2))
    demand(type(rows) is list and len(rows) == len(expected), "incomplete local upper Hessian")
    answer = []
    for row, loc in zip(rows, expected, strict=True):
        demand(type(row) is list and len(row) == 4 and type(row[0]) is type(row[1]) is int
               and tuple(row[:2]) == loc, "permuted or malformed literal local Hessian index")
        answer.append((row[0], row[1], bounds(row[2:])))
    return answer


def scatter(n, components):
    """Full sparse matrix, with each off-diagonal mirrored exactly once."""
    force = [ZERO for _ in range(n)]
    H = {}
    for ids, local_F, local_H in components:
        exact_maps(ids, n)
        demand(len(ids) == len(local_F), "unmatched local force")
        for global_row, endpoints in zip(ids, local_F, strict=True):
            force[global_row] = plus(force[global_row], endpoints)
        for a, b, endpoints in local_H:
            i, j = ids[a], ids[b]
            H[i, j] = plus(H.get((i, j), ZERO), endpoints)
            if i != j:
                H[j, i] = plus(H.get((j, i), ZERO), endpoints)
    return force, H


def ball(q):
    return arb(fmpq(q.numerator, q.denominator))


def product_endpoints(R, interval_entries, columns):
    """Apply R+ and R- to directed endpoint matrices, never to Arb radius balls."""
    n = len(R)
    positive = arb_mat(n, n, [ball(max(x, Q(0))) for row in R for x in row])
    negative = arb_mat(n, n, [ball(min(x, Q(0))) for row in R for x in row])
    lo, hi = arb_mat(n, columns), arb_mat(n, columns)
    for (i, j), (a, b) in interval_entries.items():
        lo[i, j], hi[i, j] = ball(a), ball(b)
    return positive * lo + negative * hi, positive * hi + negative * lo


def directed(lo, hi, i, j):
    a, b = lo[i, j], hi[i, j]
    demand(a.is_finite() and b.is_finite(), "nonfinite directed matrix bound")
    return bounds((Q(str(a.lower().fmpq())), Q(str(b.upper().fmpq()))))


def row_bounds(R, F, H, anchor):
    n = len(R)
    lo_f, hi_f = product_endpoints(R, {(i, 0): row for i, row in enumerate(F)}, 1)
    RFdef = [directed(lo_f, hi_f, i, 0) for i in range(n)]
    full = [plus(row, cd) for row, cd in zip(RFdef, anchor, strict=True)]
    lo_h, hi_h = product_endpoints(R, H, n)
    norm_H = [sum((absolute_cap(directed(lo_h, hi_h, i, j)) for j in range(n)), Q(0))
              for i in range(n)]
    return RFdef, full, norm_H


def classifications(eta, lower, z, r):
    lhs = eta + z * r
    exclusion = lower - (1+z)*r
    exists = z < 1 and lhs < r
    demand(not (exists and exclusion > 0), "inconsistent sufficient tests")
    return {
        "eta_upper_exact": str(eta), "RF_lower_exact": str(lower), "z_upper_exact": str(z),
        "selfmap_lhs_upper_exact": str(lhs), "selfmap_margin_lower_exact": str(r-lhs),
        "strict_selfmap": exists, "posterior_error_upper_exact": str(eta/(1-z)) if exists else None,
        "strict_exclusion_of_this_cube": exclusion > 0, "exclusion_margin_lower_exact": str(exclusion),
    }


def identity(ctx, case, entry):
    return {
        "pilot_id": ctx.pilot_id, "record_id": ctx.record_id, "field_sha256": case["field_sha256"],
        "cache_packet_sha256": case["cache_packet_sha256"], "embedded_cache_sha256": case["embedded_cache_sha256"],
        "R_binary_sha256": case["preconditioner_binary_sha256"], "model_contract_sha256": entry["model_contract_sha256"],
        "parameter_cell_sha256": ctx.row["parameter_cell_sha256"], "dose_exact": str(ctx.dose),
        "gain_exact": str(ctx.base["common_gain"]), "drive_exact": str(ctx.drive),
        "energy_scale_binary_exact": str(ctx.target.energy_scale),
    }


def prior_verified_geometry(pilot, case, plan, retained):
    raw = read(GEOMETRY / "receipt.json", PINS[GEOMETRY / "receipt.json"], retained)
    observed = read(GEOMETRY / "os_exit_observation.json", PINS[GEOMETRY / "os_exit_observation.json"], retained)
    end = read(GEOMETRY / "terminal.json", observed["terminal_sha256"], retained)
    demand(raw["status"] == "PASS_INDEPENDENT_ACTUAL_EIGHT_SINGLETON_GEOMETRY_BOXES"
           and literal_checks(raw["checks"])
           and observed["status"] == "PASS_OBSERVED_INDEPENDENT_GEOMETRY_REPLAY_TRUE_OS_EXIT"
           and type(observed["observed_exit_code"]) is int and observed["observed_exit_code"] == 0
           and observed["receipt_sha256"] == end["receipt_sha256"] == PINS[GEOMETRY / "receipt.json"]
           and end["pid"] == observed["pid"] and type(end["exit_code"]) is int and end["exit_code"] == 0,
           "unfinished independent geometry")
    row, = [g for g in raw["geometry_records"] if g["pilot_id"] == pilot and rational(g["radius_exact"]) == RADIUS]
    leaf = read(row["receipt"], row["sha256"], retained)
    demand(case["whole_box_geometry_review"] == row["receipt"]
           and case["whole_box_geometry_review_sha256"] == row["sha256"]
           and leaf["pilot_id"] == pilot and leaf["record_id"] == case["record_id"]
           and rational(leaf["radius_exact"]) == RADIUS
           and leaf["box_sha256"] == row["box_sha256"] == case["box_sha256"] == plan["context_box_sha256"]
           and leaf["proof_center_field_sha256"] == case["field_sha256"]
           and leaf["all_144_body_and_24_trace_domains"] is True
           and min(rational(leaf[key]) for key in ("strong_monotonicity_lower_exact", "min_surface_G_exact",
                                                    "min_trace_hoop_exact", "min_site_tangent_squared_exact")) > 0,
           "missing whole-cube, strictly admissible same-body/trace domains")
    return row["sha256"]


def prior_verified_spd(pilot, ctx, case, entry, plan, Rnorm, retained):
    raw = read(SPD / "receipt.json", PINS[SPD / "receipt.json"], retained)
    o = read(SPD / "os_exit_observation.json", PINS[SPD / "os_exit_observation.json"], retained)
    end = read(SPD / "terminal.json", o["terminal_sha256"], retained)
    start = read(SPD / "started.json", o["started_sha256"], retained)
    demand(raw["status"] == "PASS_INDEPENDENT_ACTUAL_FOUR_EXACT_SYMMETRIC_PRECONDITIONER_PARTS_SPD"
           and literal_checks(raw["checks"]) and raw["root_or_minimum_certified"] is False
           and o["status"] == "PASS_OBSERVED_INDEPENDENT_ACTUAL_SPD_TRUE_OS_EXIT"
           and type(o["observed_exit_code"]) is int and o["observed_exit_code"] == 0
           and o["receipt_sha256"] == end["receipt_sha256"] == PINS[SPD / "receipt.json"]
           and raw["pid"] == o["pid"] == end["pid"] == start["pid"]
           and type(end["exit_code"]) is int and end["exit_code"] == 0,
           "unfinished independent SPD proof")
    row, = [c for c in raw["ordered_cases"] if c["pilot_id"] == pilot]
    leaf = read(row["case_receipt_file"], row["case_receipt_sha256"], retained)
    demand(leaf["target_identity"] == row["target_identity"] == identity(ctx, case, entry)
           and leaf["pilot_id"] == pilot and leaf["record_id"] == case["record_id"]
           and leaf["exact_target_symmetric_part_spd"] is True and leaf["R_symmetry_assumed"] is False
           and row["status"] == "PASS_INDEPENDENT_EXACT_SYMMETRIC_PART_SPD"
           and 0 <= rational(leaf["fresh_congruence_q_upper_exact"]) <= rational(leaf["max_reviewed_q_upper_exact"]) < 1
           and leaf["max_reviewed_q_upper_exact"] == row["max_reviewed_q_upper_exact"]
           and rational(leaf["exact_R_infinity_norm"]) == rational(row["exact_R_infinity_norm"]) == Rnorm > 0
           and leaf["cache_packet_sha256"] == case["cache_packet_sha256"]
           and leaf["cached_point_proof_sha256"] == case["verified_point_sha256"]
           and leaf["larger_radius_box_sha256"] == case["box_sha256"]
           and leaf["derivative_case_sha256"] == plan["derivative"]["native_case_sha256"],
           "SPD receipt is independent but does not match this exact target and R")
    return row["case_receipt_sha256"]


def observe(pilot, admission_hash, retained):
    folder = OUT / "actual_v3" / pilot
    ext = read(folder / "root_terminal_admission.json", admission_hash, retained)
    demand(ext["status"] == "PASS_ROOT_OBSERVED_SINGLETON_GLOBAL_REDUCTION_TERMINAL"
           and ext["pilot_id"] == pilot and ext["source_sha256"] == PINS[CALLER]
           and type(ext["pid"]) is int and ext["pid"] > 0
           and type(ext["unified_exec_session"]) is int and ext["unified_exec_session"] > 0
           and type(ext["observed_exit_code"]) is int and ext["observed_exit_code"] == 0,
           "external true OS exit-0 gate required before reading any actual data")
    files = {name: folder / (name+".json") for name in ("started", "terminal", "results", "assembled_operands")}
    paths = {name: relative(p) for name, p in files.items()}
    records = {name: read(files[name], ext["sha256"][path], retained) for name, path in paths.items()}
    start, end, result, assembled = (records[k] for k in files)
    demand(start["status"] == "RUNNING" and end["status"] == "FINISHED"
           and result["status"] == "SPATIAL_REDUCTION_CANDIDATE_AWAITING_INDEPENDENT_REVIEW"
           and assembled["status"] == "ADMITTED_TWO_PROOF_SIGNED_SPATIAL_OPERANDS"
           and type(end["exit_code"]) is int and end["exit_code"] == 0
           and start["pid"] == end["pid"] == result["pid"] == ext["pid"]
           and all(blob["pilot_id"] == pilot for blob in records.values())
           and end["results_sha256"] == ext["sha256"][paths["results"]]
           and end["assembled_operands_sha256"] == result["assembled_operands_sha256"] == ext["sha256"][paths["assembled_operands"]]
           and all(blob["bindings"] == start["bindings"] for blob in (end, result, assembled))
           and result["exact_spatial_root_certified"] is False
           and result["independent_actual_reduction_review_required"] is True
           and type(result["force_H_AD_or_equilibrium_jobs_run"]) is int and result["force_H_AD_or_equilibrium_jobs_run"] == 0,
           "mismatched or mixed producer actual reductions")
    return result, assembled


def admit(pilot, result, retained):
    # Use the cold peer-reviewed structural and mathematical proof-admission logic
    # as a premise. Actual numerical assembly and contraction below are separate.
    sys.path.insert(0, str(ROOT / "scripts"))
    import transfer_ambiguity_exact_bridge_reduce_driver_v3 as reviewed
    links = result["bindings"]
    args = SimpleNamespace(pilot_id=pilot)
    for stem, prefix in (("admission", "producer_admission"), ("raw_observation", "raw_os"),
                         ("companion_observation", "companion_os")):
        path = (ROOT / links[prefix+"_file"]).resolve()
        sha = links[prefix+"_sha256"]
        read(path, sha, retained)
        setattr(args, stem, path)
        setattr(args, stem+"_sha256", sha)
    protocol, bindings = reviewed.gate(GATE, PINS[GATE])
    r0, r1, plan, plinks, entries, checked = reviewed.admitted(args, bindings)
    demand(checked == links, "external proof sources differ from producer admission")
    for source, sha in protocol["source_pins"].items():
        read(source, sha, retained, json_payload=False)
    for source, sha in protocol["base_authority"]["source_and_prerequisite_sha256"].items():
        read(source, sha, retained, json_payload=False)
    ctx, case, entry, routes, anchor, z, proof_digest, cached_eta = reviewed.base.context(
        pilot, r0, plan, plinks, protocol["base_authority"])
    for path, key in ((case["box"], "box_sha256"), (case["verified_point"], "verified_point_sha256"),
                      (case["cache_packet"], "cache_packet_sha256"), (entry["routes"], "routes_sha256"),
                      (entry["model_contract"], "model_contract_sha256"), (case["field"], "field_sha256")):
        h = case[key] if key in case else entry[key]
        read(path, h, retained, json_payload=False)
    for path, h in r0["checked_packet_sha256"].items():
        read(path, h, retained, json_payload=False)
    rawdir, newdir = reviewed.RAW / pilot, reviewed.V3 / pilot
    for directory, key in ((rawdir, "raw_native"), (newdir, "companion")):
        for phase in ("started", "terminal", "receipt"):
            read(directory / (phase+".json"), links[key+"_"+phase+"_sha256"], retained)
    demand(r1["point_force_scope"] == POINT if "point_force_scope" in r1 else True,
           "wrong companion point domain")
    return reviewed, r0, r1, plan, plinks, entries, ctx, case, entry, routes, anchor, z, proof_digest, cached_eta


def local_packet(component, raw, fresh, route, plan_leaf, expected_kind, expected_slot):
    records = (component, raw, fresh, plan_leaf)
    ids = exact_maps(route["global_indices"], N)
    for blob in records:
        demand(blob["kind"] == expected_kind and type(blob["slot"]) is int and blob["slot"] == expected_slot
               and exact_maps(blob["global_indices"], N) == ids, "component map disagrees with factory route")
    f, h = component["point_force"], component["whole_box_hessian"]
    demand(f["scope"] == POINT and h["scope"] == CUBE and rational(h["radius_exact"]) == RADIUS
           and raw["point_force_scope"] == fresh["point_force_scope"] == POINT
           and raw["whole_box_H_scope"] == fresh["whole_box_H_scope"] == CUBE
           and raw["containment_failures"]["direct_point_cache"] == []
           and set(fresh["containment_failures"]) == {"direct_point_cache", "point_force", "whole_box_H"}
           and all(not misses for misses in fresh["containment_failures"].values()),
           "missing same-point cache/direct proof or wrong derivative domain")
    F = list(map(bounds, f[FKEY]))
    f_fresh = list(map(bounds, fresh["independent_point_F_exact"]))
    H = local_hessian(h[HKEY], len(ids))
    h_fresh = local_hessian(fresh["independent_whole_box_H_exact"], len(ids))
    demand(len(F) == len(f_fresh) == len(ids)
           and all(contains(a, b) for a, b in zip(F, f_fresh, strict=True))
           and all(contains(a[2], b[2]) for a, b in zip(H, h_fresh, strict=True)),
           "independent leaf is not strictly enclosed")
    if expected_kind == "bulk":
        demand(raw["containment_failures"]["whole_box_H"] == []
               and fresh["independent_whole_box_H_exact"] == raw["independent_whole_box_H_exact"]
               and fresh["fresh_surface_H_tails"] == [], "bulk H did not retain its successful native leaf")
    else:
        demand(len(fresh["fresh_surface_H_tails"]) == len(h_fresh), "missing freshly verified surface H entries")
    return ids, F, H


def validate_producer(result, zq, own_force, own_derivative):
    r = result["arithmetic"]
    old_full = list(map(bounds, r["signed_RF_I_rows_exact"]))
    old_rows = list(map(rational, r["RH_absolute_row_sum_bounds_exact"]))
    demand(len(old_full) == len(old_rows) == len(own_force) == N
           and all(max(a[0], b[0]) <= min(a[1], b[1]) for a, b in zip(old_full, own_force, strict=True)),
           "independent and producer force enclosures are disjoint")
    old_eta = max(map(absolute_cap, old_full))
    old_lower = max(map(absolute_floor, old_full))
    old_b = max(old_rows)
    old_z = zq + old_b
    checks = classifications(old_eta, old_lower, old_z, RADIUS)
    demand(rational(r["eta_I_upper_exact"]) == old_eta and rational(r["RF_I_norm_lower_exact"]) == old_lower
           and rational(r["b_upper_exact"]) == old_b and rational(r["ZQ_upper_exact"]) == zq
           and rational(r["z_upper_exact"]) == old_z and rational(r["radius_exact"]) == RADIUS
           and all(r[k] == checks[k] for k in ("selfmap_lhs_upper_exact", "selfmap_margin_lower_exact",
                                                 "posterior_error_upper_exact", "strict_exclusion_of_this_cube",
                                                 "exclusion_margin_lower_exact"))
           and r["strict_contraction_and_selfmap"] is checks["strict_selfmap"],
           "producer-derived scalars or strict decisions do not recompute")
    own_eta = max(map(absolute_cap, own_force))
    own_lower = max(map(absolute_floor, own_force))
    own_b = max(own_derivative)
    proof = classifications(max(own_eta, old_eta), min(own_lower, old_lower), zq+max(own_b, old_b), RADIUS)
    # Leave an insufficient producer test unresolved even if independent arithmetic could otherwise improve it.
    demand(not proof["strict_selfmap"] or checks["strict_selfmap"], "cannot upgrade producer HOLD")
    demand(not proof["strict_exclusion_of_this_cube"] or checks["strict_exclusion_of_this_cube"],
           "cannot upgrade producer exclusion HOLD")
    return proof, classifications(own_eta, own_lower, zq+own_b, RADIUS), str(own_b)


def run(args):
    retained = {}
    protocol = read(HERE / "protocol.json", args.protocol_sha256, retained)
    pf = read(HERE / "preflight.json", args.preflight_sha256, retained)
    self_sha = digest(Path(__file__))
    demand(self_sha == protocol["source_sha256"] == pf["source_sha256"] and pf["status"] == "PASS"
           and pf["protocol_sha256"] == args.protocol_sha256 and literal_checks(pf["checks"]),
           "current cold independent mathematical preflight required")
    for path, pin in PINS.items():
        read(path, pin, retained, json_payload=False)
    result, assembled = observe(args.pilot_id, args.admission_sha256, retained)
    folder = HERE / "actual_replay" / args.pilot_id
    folder.mkdir(parents=True, exist_ok=False)
    phase = {"pid": os.getpid(), "pilot_id": args.pilot_id, "source_sha256": self_sha,
             "protocol_sha256": args.protocol_sha256, "preflight_sha256": args.preflight_sha256,
             "root_admission_sha256": args.admission_sha256}
    save(folder / "started.json", {"status": "RUNNING", **phase})
    end = {"status": "FAILED", "exit_code": 1, **phase}
    try:
        data = admit(args.pilot_id, result, retained)
        reviewed, r0, r1, plan, links, entries, ctx, case, entry, routes, anchor, zq, proof_digest, qeta = data
        adj = []
        for (kind, slot), pack in zip(ORDER, entries, strict=True):
            row, new, prod, proof, old_path, new_path = pack
            production = read(prod, row["packet_sha256"], retained)
            selection = read(proof, row["plan_file_sha256"], retained)
            old = read(old_path, row["replay_sha256"], retained)
            leaf = read(new_path, new["replay_sha256"], retained)
            reviewed.leaf_refinement(row, new, old, leaf, production, tuple(map(rational, plan["tail_caps_exact"])))
            reviewed.base.component_operands(row, leaf, production, selection, routes[kind][slot],
                                             case["record_id"], links, proof_digest(selection))
            adj.append(local_packet(production, old, leaf, routes[kind][slot], selection, kind, slot))
        demand(len(adj) == 168 and len(anchor) == N == ctx.R.shape[0] == ctx.R.shape[1]
               and ctx.R.dtype.str == "<f8" and ctx.R.flags.writeable is False
               and case["preconditioner_binary_sha256"] == sha256(ctx.R.tobytes(order="C")).hexdigest()
               and qeta >= max(map(absolute_cap, anchor)) and 0 <= zq < 1,
               "incomplete signed assembly or wrong binary preconditioner, anchor or whole-box ZQ")
        R = [list(map(Q.from_float, row)) for row in ctx.R.tolist()]
        RN = max(sum(map(abs, row), Q(0)) for row in R)
        geomsha = prior_verified_geometry(args.pilot_id, case, plan, retained)
        spdsha = prior_verified_spd(args.pilot_id, ctx, case, entry, plan, RN, retained)
        F, H = scatter(N, adj)
        demand(list(map(bounds, assembled["signed_F_defect_exact"])) == F
               and assembled["whole_box_signed_H_defect_upper_exact"] ==
               [[i,j,str(a),str(b)] for (i,j),(a,b) in sorted(H.items()) if i<=j]
               and assembled["maps"] == [{"kind": kind, "slot": slot, "global_indices": pkt[0]}
                                          for (kind, slot), pkt in zip(ORDER, adj, strict=True)]
               and list(map(bounds, assembled["signed_R_F_cache_exact"])) == anchor
               and rational(assembled["ZQ_upper_exact"]) == zq
               and rational(assembled["radius_exact"]) == RADIUS
               and assembled["closed_cell_exact"] == ctx.row["parameter_cell"]
               and assembled["cache_packet_sha256"] == case["cache_packet_sha256"]
               and assembled["binary_R_sha256"] == case["preconditioner_binary_sha256"]
               and rational(assembled["same_exact_drive"]) == ctx.drive,
               "independent exact scatter/anchor disagrees with assembled producer operands")
        arb_context.prec = PREC
        fd, full, hb = row_bounds(R, F, H, anchor)
        final, own, own_b = validate_producer(result, zq, full, hb)
        success = final["strict_selfmap"]
        posterior = rational(final["posterior_error_upper_exact"]) if success else None
        z = rational(final["z_upper_exact"])
        goal = rational(case["observer_goal_exact"])
        demand(goal == rational(result["observer_design_goal_exact"]), "observer design goal changed")
        # D(E) gives F, D(F) is symmetric on the C^2 admissible neighborhood,
        # and the prior same-matrix SPD-part certificate yields H_I positive definite.
        lower_H = (1-z)/RN if success else None
        status = ("PASS_INDEPENDENT_EXACT_SPATIAL_ROOT_AND_STRICT_LOCAL_MINIMUM" if success else
                  "PASS_INDEPENDENT_EXACT_SPATIAL_EXCLUSION_OF_MATCHED_CUBE" if final["strict_exclusion_of_this_cube"] else
                  "UNRESOLVED_INDEPENDENT_EXACT_SPATIAL_ROOT_TEST")
        packet = {"status": status, **phase, "target_identity": identity(ctx, case, entry),
                  "scope": "Single fixed-input locally embedded finite-element model with exact spatial integration and chosen normalization.",
                  "radius_exact": str(RADIUS), "precision_bits": PREC, "counts": {"bulk": 144, "surface": 24, "free": N},
                  "source_and_actual_sha256": dict(sorted(retained.items())),
                  "same_geometry_admissibility_sha256": geomsha, "same_R_spd_part_sha256": spdsha,
                  "ZQ_upper_exact": str(zq), "independent_RH_defect_norm_upper_exact": own_b,
                  "conservative_compared_bounds": final, "independent_bounds": own,
                  "whole_box_normalized_H_euclidean_lower_exact": str(lower_H) if success else None,
                  "exact_R_infinity_norm": str(RN), "observer_state_design_goal_exact": str(goal),
                  "posterior_below_observer_state_design_goal": success and posterior < goal,
                  "signed_R_force_defect_rows_exact": [[str(x) for x in a] for a in fd],
                  "signed_R_full_force_rows_exact": [[str(x) for x in a] for a in full],
                  "RH_defect_absolute_row_sum_upper_exact": list(map(str, hb)),
                  "premises": "Shared cold-reviewed structural and leaf calculus validation (both proof phases); exact scatter, directed interval products, rational norms and theorem application independently evaluated.",
                  "new_AD_Hessian_force_moment_or_solve_jobs": 0,
                  "full_five_gate_goal_achieved": False}
        h = save(folder / "receipt.json", packet)
        end.update(status="FINISHED", exit_code=0, receipt_sha256=h, classification=status)
        print(json.dumps({"status": status, "receipt_sha256": h, "pid":os.getpid()}), flush=True)
    finally:
        save(folder / "terminal.json", end)


def fixtures():
    arb_context.prec = PREC
    # Nontrivial permutations, cancellation and distinct local upper indexing.
    comps = [([1, 0], [bounds((2,3)),bounds((-1,1))],
              local_hessian([[0,0,"1/7","3/7"],[0,1,"-2/5","1/5"],[1,1,"1","2"]], 2)),
             ([0], [bounds((1,2))], local_hessian([[0,0,"-1","0"]], 1))]
    F, H = scatter(2, comps)
    check_scatter = F == [(Q(0),Q(3)),(Q(2),Q(3))] and H == {
        (0,0):(Q(0),Q(2)), (1,1):(Q(1,7),Q(3,7)),
        (0,1):(Q(-2,5),Q(1,5)), (1,0):(Q(-2,5),Q(1,5))}
    R = [[Q(-1,3),Q(1,5)],[Q(4,7),Q(-2,3)]]
    def brute(entries, j):
        for weights in R:
            terms = [sorted(w*x for x in entries[i,j]) for i,w in enumerate(weights)]
            yield tuple(sum(t[k] for t in terms) for k in range(2))
    fd, full, rows = row_bounds(R,F,H,[(Q(-1),Q(3)),(Q(1),Q(2))])
    f_ref = list(brute(dict(zip(((0,0),(1,0)),F)),0))
    h_ref = list(zip(*(brute(H,j) for j in range(2)), strict=True))
    direct = all(contains(a,b) for a,b in zip(fd,f_ref,strict=True))
    single = all(contains(out, plus(fi, ai)) for out,fi,ai in zip(full,f_ref,[(Q(-1),Q(3)),(Q(1),Q(2))],strict=True))
    norms = all(a >= sum(map(absolute_cap,b)) for a,b in zip(rows,h_ref,strict=True))
    strict = (classifications(Q(1,8), Q(0), Q(1,2), Q(1))["strict_selfmap"]
              and not classifications(Q(1,2), Q(0), Q(1,2), Q(1))["strict_selfmap"]
              and classifications(Q(2), Q(2), Q(1,2), Q(1))["strict_exclusion_of_this_cube"])
    def reject(fn):
        try:
            fn()
        except (ValueError, TypeError):
            return True
        return False
    return {
        "signed_permuted_exact_scatter_mirror_unshared_diagonal": check_scatter,
        "directed_sign_splits_contain_exact_rational_products": direct and norms,
        "preconditioned_cache_anchor_added_exactly_once": single,
        "strict_root_and_exclusion_boundary": strict,
        "typed_int_only_and_signed_interval_order": reject(lambda: exact_maps([False],2))
          and reject(lambda: local_hessian([[0.0,0,"1","2"]],1)) and reject(lambda: bounds((2,1)))
          and reject(lambda: bounds((0.5,2))),
    }


def prepare():
    retained = {}
    for path, value in PINS.items():
        read(path, value, retained, json_payload=False)
    gates = fixtures()
    demand(literal_checks(gates), "independent cold exact arithmetic fixture failed")
    h = digest(Path(__file__))
    proto_hash = save(HERE / "protocol.json", {
        "status": "FROZEN_INDEPENDENT_ACTUAL_SIGNED_REDUCTION",
        "source_sha256": h, "source_and_prerequisite_sha256": retained,
        "radius_exact": str(RADIUS), "precision_bits": PREC, "pilot_ids": PILOTS,
        "full_actual_authentication": "Per-pilot fresh external root observation of OS exit 0, matching started/terminal/results/assembled; coherent producer/plan/raw/companion leaves.",
        "cold_actual_components_read": 0, "cold_611_matrix_jobs": 0,
    })
    preflight = save(HERE / "preflight.json", {
        "status": "PASS", "checks": gates, "source_sha256":h,"protocol_sha256":proto_hash,
        "cold_actual_components_read": 0, "cold_611_matrix_jobs": 0})
    print(json.dumps({"source_sha256":h,"protocol_sha256":proto_hash,"preflight_sha256":preflight}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    stages = parser.add_subparsers(dest="phase", required=True)
    stages.add_parser("prepare")
    act = stages.add_parser("run")
    act.add_argument("--pilot-id", required=True, choices=PILOTS)
    for flag in ("admission-sha256", "protocol-sha256", "preflight-sha256"):
        act.add_argument("--"+flag, required=True)
    options = parser.parse_args()
    if options.phase == "prepare":
        prepare()
    else:
        run(options)

