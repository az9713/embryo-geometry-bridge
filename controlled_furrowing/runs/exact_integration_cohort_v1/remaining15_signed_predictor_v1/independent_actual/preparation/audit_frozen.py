"""Independently replay D's signed, original-cube predictor bound."""
import argparse
import ast
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
import os
from pathlib import Path
import sys
from unittest.mock import patch

from flint import arb, arb_mat, ctx

SOURCE = Path(__file__).resolve()
OUT = SOURCE.parent
STAGE = OUT.parent
ROOT = STAGE.parents[2]
PREP = OUT / "preparation"
sys.path.insert(0, str(ROOT / "scripts"))
import exact_signed_cohort_predictor as J

RECORD = "optoGEF2_small_col_D"
PID = 11408
PREC = 256
N = 611
ORDER = tuple(J.P.NAMES)
USE = (J.SHARED[0], J.SHARED[1], J.SHARED[2], J.GAIN)
PRODUCER_OS = {"observer": "root", "unified_exec_session": 95017, "pid": PID, "observed_exit_code": 0}
OLD = ROOT / "runs/exact_integration_cohort_v1/L_defect_diagnostic/signed_predictor_L_v1/independent_actual/replay.py"
OLD_SHA = "42d3512869b18fd686390274dbf524b4b26667ab167c0712b4e4ab28e2fbaae4"
READBACK_HELPERS = {"pair", "ball_pair", "aloft", "add", "mul", "inside", "overlap",
                    "site_gradient_columns", "build_independently", "predictor_enclosures", "synthetic"}
PINS = {
    STAGE / "preparation/protocol.json": "314b71743c1f90cdeb7bd1e8b54f312a1b2c6b264eef08ba8a9421c00d3f79ac",
    STAGE / "preparation/preflight.json": "ff5d044978df37ced65ab16cb599790ef87b5d8bb8b038c97987515f9c0e91fc",
    STAGE / "preparation/prelaunch_independent_review.json": "ee73604f3d202b81c2d63a4f006ea51eb90a980f8d66595981dc0a1c7de6e044",
    J.SOURCE: "5a7e6bdffa829335d1c88129d21f73f284a886053ab85d7bb4dbad3c8a6d2b6f",
    STAGE / "actual/started.json": "72fc7b81585a3a6529a0c21051e474752df1567a4178cd5c5d8c1a8dc3462c34",
    STAGE / "actual/terminal.json": "cf1588cf9b546e762a2796932e6736038faa180c1988bd98c66ff383de05e12c",
    STAGE / "actual/results.json": "cbd1cac250722c99c7cdf5b8100b3ba53352cadd49c220399420919af2704fde",
    OLD: OLD_SHA,
}


def now():
    return datetime.now(timezone.utc).isoformat()


def load(path, expected, parse=True):
    data = path.read_bytes()
    assert sha256(data).hexdigest() == expected, path
    return json.loads(data) if parse else data


def write(path, value):
    with path.open("x", encoding="utf-8", newline="\n") as f:
        json.dump(value, f, indent=2, default=str, allow_nan=False)
        f.write("\n")
    return J.sha(path)


def helper_namespace():
    tree = ast.parse(load(OLD, OLD_SHA, False))
    selected = [item for item in tree.body if isinstance(item, ast.FunctionDef) and item.name in READBACK_HELPERS]
    assert {item.name for item in selected} == READBACK_HELPERS
    # Evaluate only pure function definitions and explicit record-neutral constants.
    # No top-level statements, path tables, records, or data from the earlier L job.
    namespace = {key: globals()[key] for key in ("arb", "arb_mat", "Q", "product", "Counter", "ORDER", "USE", "N")}
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(OLD), "exec"), namespace)
    return namespace


def classify(rows, z, r):
    assert 0 <= z < 1 and r > 0 and rows
    lows = [max(Q(0), a, -b) for a,b in rows]
    highs = [max(abs(a), abs(b)) for a,b in rows]
    ilow = max(range(len(rows)), key=lows.__getitem__)
    ihigh = max(range(len(rows)), key=highs.__getitem__)
    lower, upper = lows[ilow], highs[ihigh]
    rhs_excl = (1+z)*r
    rhs_exists = (1-z)*r
    exists = upper < rhs_exists
    excluded = lower > rhs_excl
    assert not (exists and excluded)
    return {
        "eta_I_upper_exact": upper, "uniform_RF_I_lower_exact": lower,
        "lower_witness_row_zero_based": ilow, "upper_witness_row_zero_based": ihigh,
        "strict_eta_budget_exact": rhs_exists, "selfmap_lhs_exact": upper+z*r,
        "exclusion_threshold_exact": rhs_excl, "exclusion_margin_exact": lower-rhs_excl,
        "selfmap_slack_exact": rhs_exists-upper,
        "sufficient_existence_test_passes": exists,
        "sufficient_old_cube_exclusion_passes": excluded,
        "conditional_root_error_exact": upper/(1-z) if exists else None,
        "zero_excluding_rows": sum(v > 0 for v in lows),
    }


def incoming():
    files = STAGE / "preparation"
    afiles = STAGE / "actual"
    prep = load(files / "protocol.json", PINS[files / "protocol.json"])
    preflight = load(files / "preflight.json", PINS[files / "preflight.json"])
    gate = load(files / "prelaunch_independent_review.json", PINS[files / "prelaunch_independent_review.json"])
    started = load(afiles / "started.json", PINS[afiles / "started.json"])
    terminal = load(afiles / "terminal.json", PINS[afiles / "terminal.json"])
    assert all(item["status"] == "PASS" and item["checks"] and all(v is True for v in item["checks"].values())
               for item in (preflight, gate))
    assert J.sha(J.SOURCE) == PINS[J.SOURCE] == gate["source_sha256"]
    assert prep["source_sha256"] == preflight["source_sha256"] == gate["source_sha256"]
    assert preflight["protocol_sha256"] == gate["protocol_sha256"] == PINS[files / "protocol.json"]
    assert gate["preflight_sha256"] == PINS[files / "preflight.json"]
    assert terminal["status"] == "FINISHED" and type(terminal["exit_code"]) is int and terminal["exit_code"] == 0
    assert terminal["record_id"] == started["record_id"] == prep["record_id"] == gate["record_id"] == RECORD
    assert terminal["pid"] == started["pid"] == PRODUCER_OS["pid"] and PRODUCER_OS["observed_exit_code"] == 0
    assert terminal["results_sha256"] == PINS[afiles / "results.json"]
    for key in ("source_sha256", "protocol_sha256", "preflight_sha256",
                "global_reduction_review_sha256", "independent_signed_inputs_sha256", "global_reduction_results_sha256"):
        assert gate[key] == started[key] == terminal[key]
    assert started["prelaunch_review_sha256"] == terminal["prelaunch_review_sha256"] == PINS[files / "prelaunch_independent_review.json"]
    assert prep["input_hashes"][J.REVIEW] == gate["global_reduction_review_sha256"]
    assert prep["input_hashes"][J.INPUTS] == gate["independent_signed_inputs_sha256"]
    assert prep["input_hashes"][J.PRODUCER] == gate["global_reduction_results_sha256"]
    links = dict(prep["input_hashes"])
    for path, h in PINS.items():
        links[path.relative_to(ROOT).as_posix()] = h
        if path != afiles / "results.json":
            load(path, h, False)
    return prep, links, terminal


def source_plan(source_hash, links):
    return {
        "status": "FROZEN_D_INDEPENDENT_SIGNED_REPLAY_PREPARATION",
        "record_id": RECORD, "precision_bits": PREC, "audit_source_sha256": source_hash,
        "input_sha256": links, "producer_os": PRODUCER_OS, "cache_site_count": 9504,
        "helpers_extracted_as_function_definitions_from_prior_review": sorted(READBACK_HELPERS),
        "scope": [
            "D original closed parameter cell and parameter-matched original numerical cube only.",
            "Use the selected D binary cache coordinates, free order, R, and S. Rebuild cached gR, signed affine columns via local H*S, and rowwise zero-radius path remainder from 9504 sites.",
            "The same five offsets multiply combined signed slopes once. The alternate independent exact bulk gain array and reviewed centered force arrays are numerical premises.",
            "Shared with previously reviewed methods: binary cache restoration, element maps, local energy derivatives, full whole-box derivative and previously reviewed centered integral arrays. No integration or solve.",
            "An independent bound is constructed and classified on its own. Any interval overlap between implementations is a diagnostic only."
        ],
    }


def prepare():
    _, links, _ = incoming()
    helpers = helper_namespace()
    ctx.prec = PREC
    J.P.core.A.cache_clear()
    stop = AssertionError("preparation entered actual cached site machinery")
    with patch.object(J.P, "build_point", side_effect=stop), \
            patch.object(J.P, "variation", side_effect=stop), \
            patch.object(J.P.core, "restore_cache", side_effect=stop), \
            patch.object(J.P.core, "local_terms", side_effect=stop), \
            patch.object(J.P.core, "assemble_point", side_effect=stop):
        helpers["synthetic"](J.P.core)
    c = Q(1,2)
    assert classify([(-c,c)],c,1)["sufficient_existence_test_passes"] is False
    assert classify([(-2,-Q(3,2))],c,1)["sufficient_old_cube_exclusion_passes"] is False
    assert classify([(-2,-Q(13,8))],c,1)["sufficient_old_cube_exclusion_passes"] is True
    assert classify([(-Q(1,4),Q(1,4))],c,1)["conditional_root_error_exact"] == c
    assert not (OUT / "started.json").exists()
    PREP.mkdir(exist_ok=False)
    source_hash = J.sha(SOURCE)
    with (PREP / "audit_frozen.py").open("xb") as h:
        h.write(SOURCE.read_bytes())
    protocol_sha = write(PREP / "protocol.json", source_plan(source_hash, links))
    preflight_sha = write(PREP / "preflight.json", {
        "status": "PASS_PREPARATION_ONLY", "audit_source_sha256": source_hash, "protocol_sha256": protocol_sha,
        "checks": {"root_observed_producer_OS_and_hash_links_match": True, "bulk_gain_surface_gradient_fixtures": True,
                   "theta_first_state_path_and_cancellation_fixtures": True, "strict_two_thresholds_and_equality": True,
                   "producer_result_payload_not_opened_and_all_cached_sites_guarded": True},
    })
    print(json.dumps({"source_sha256": source_hash, "protocol_sha256": protocol_sha, "preflight_sha256": preflight_sha}))


def replay(protocol_sha, preflight_sha):
    protocol = load(PREP / "protocol.json", protocol_sha)
    pf = load(PREP / "preflight.json", preflight_sha)
    assert pf["status"] == "PASS_PREPARATION_ONLY" and pf["checks"] and all(v is True for v in pf["checks"].values())
    assert J.sha(SOURCE) == J.sha(PREP / "audit_frozen.py") == pf["audit_source_sha256"] == protocol["audit_source_sha256"]
    assert pf["protocol_sha256"] == protocol_sha
    production, links, term = incoming()
    assert protocol == source_plan(J.sha(SOURCE), links)
    ppath = STAGE / "actual/results.json"
    producer = load(ppath, PINS[ppath])
    assert producer["record_id"] == RECORD and producer["status"] == "UNRESOLVED_SIGNED_PREDICTOR"
    for key in ("source_sha256", "protocol_sha256", "preflight_sha256", "prelaunch_review_sha256",
                "global_reduction_review_sha256", "independent_signed_inputs_sha256", "global_reduction_results_sha256"):
        assert producer[key] == term[key]
    assert producer["input_hashes"] == production["input_hashes"]
    assert producer["actual_independent_scope_review_required"] is True
    assert producer["independently_certified_root_or_exclusion"] is False and producer["global_model_absence_claimed"] is False

    forbid = AssertionError("producer point, combination, or variation called by independent replay")
    with patch.object(J.P, "build_point", side_effect=forbid), \
            patch.object(J.P, "variation", side_effect=forbid), \
            patch.object(J.P.core, "assemble_point", side_effect=forbid), \
            patch.object(J.signed, "affine_rows", side_effect=forbid), \
            patch.object(J.signed, "combine_slopes", side_effect=forbid):
        live, record, (bulk, surface, gain, matrix, delta), _, oldz, radius, b, z, norm, provider = J.input_bundle()
        assert record.record_id == RECORD and z == oldz+b and radius == record.old_radius == Q(1,10**11)
        assert live == production["input_hashes"]
        assert list(map(lambda a: tuple(map(Q,a)), producer["closed_parameter_offsets_exact"])) == delta
        assert producer["same_closed_C_sha256"] == production["same_closed_C_sha256"] == record.row["parameter_cell_sha256"]
        for key, value in (("r_exact",radius), ("z_upper_exact",z), ("old_cached_Z_upper_exact",oldz),
                           ("b_upper_exact",b), ("exact_R_infinity_norm",norm)):
            assert Q(producer[key]) == value
        ctx.prec = PREC
        J.P.core.A.cache_clear()
        import torch
        torch.set_num_threads(1)
        rt = {"utc":now(), "pid":os.getpid(), "protocol_sha256":protocol_sha, "preflight_sha256":preflight_sha}
        write(OUT / "started.json", rt)
        try:
            h = helper_namespace()
            center = dict(zip(USE, (bulk, surface, matrix, gain), strict=True))
            cached, affine, rem, domains = h["build_independently"](J.P.core, record, delta)
            anchors, coeffs, rows = h["predictor_enclosures"](cached, affine, rem, center, delta)
            assert len(rows) == len(cached) == len(affine) == len(rem) == N
            assert all(len(row) == 5 for row in affine)

            pair = h["pair"]
            pc = [pair(x) for x in producer["cached_signed_anchor_gR"]]
            pm = [[pair(x) for x in row] for row in producer["cached_signed_affine_M"]]
            pe = list(map(Q, producer["nonnegative_cached_predictor_remainder_upper_exact"]))
            assert len(pc) == len(pm) == len(pe) == N and min(pe) >= 0
            assert all(len(row) == 5 for row in pm)
            pa, ps, pr = h["predictor_enclosures"](pc, pm, pe, center, delta)
            assert pa == [pair(x) for x in producer["signed_combined_origins"]]
            assert ps == [[pair(x) for x in row] for row in producer["signed_combined_five_slopes"]]
            assert pr == [pair(x) for x in producer["signed_full_integral_RF_predictor_intervals"]]

            own_decision = classify(rows, z, radius)
            prod_decision = classify(pr, z, radius)
            for key, value in prod_decision.items():
                if key not in producer:
                    continue
                got = Q(producer[key]) if isinstance(value, Q) else producer[key]
                assert got == value, key
            assert producer["conditional_curvature"] is None
            flags = ("sufficient_existence_test_passes", "sufficient_old_cube_exclusion_passes")
            separately_unresolved = all(not d[k] for d in (own_decision, prod_decision) for k in flags)
            own_overlap = h["overlap"]
            overlap_counts = {
                "independent_and_producer_cached_anchor_overlap": sum(own_overlap(x,y) for x,y in zip(cached,pc,strict=True)),
                "independent_and_producer_cached_affine_overlap": sum(own_overlap(x,y) for a,b0 in zip(affine,pm,strict=True) for x,y in zip(a,b0,strict=True)),
                "independent_and_producer_full_RF_overlap": sum(own_overlap(x,y) for x,y in zip(rows,pr,strict=True)),
            }

            for name, digest in links.items():
                load(ROOT / name, digest, False)
            result = {
                "status": "PASS_INDEPENDENT_D_SIGNED_REPLAY_UNRESOLVED" if separately_unresolved else "INDEPENDENT_D_SIGNED_REPLAY_DECISION_REQUIRES_SCOPE_REVIEW",
                "record_id": RECORD, "producer_os": PRODUCER_OS, "same_closed_C_sha256": producer["same_closed_C_sha256"],
                "root_or_exclusion_certified": False, "precision_bits": PREC,
                "whole_box_radius_exact": radius, "whole_box_z_exact": z, "whole_box_b_exact": b, "old_cached_z_exact": oldz,
                "replayed_domains": domains, "shared_exact_five_offsets": delta,
                "independent_cached_gR": cached, "independent_cached_affine_M": affine,
                "independent_cached_remainder_upper": rem, "independent_full_origins": anchors,
                "independent_full_slopes": coeffs, "independent_full_RF": rows,
                "independent_classification": own_decision, "producer_classification_readback": prod_decision,
                "producer_endpoint_composition_exactly_reconstructed": True, "diagnostic_overlap_counts": overlap_counts,
                "scope": protocol["scope"], "input_sha256": links, "audit_source_sha256": J.sha(SOURCE), "runtime": rt,
            }
            digest = write(OUT / "results.json", result)
            pd = write(OUT / "terminal.json", {**rt, "utc":now(), "results_sha256":digest, "status":"FINISHED","exit_code":0})
            print(json.dumps({"status":result["status"], "results_sha256":digest, "terminal_sha256":pd,
                              **{k: float(v) if isinstance(v,Q) else v for k,v in own_decision.items()},
                              "z":float(z), "r":float(radius), **overlap_counts}))
        except BaseException as error:
            write(OUT / "terminal.json", {**rt, "utc":now(),"status":"FAILED", "exit_code":1, "message":str(error)})
            raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    stages = parser.add_subparsers(dest="stage", required=True)
    stages.add_parser("prepare")
    actual = stages.add_parser("replay")
    actual.add_argument("--protocol-sha", required=True)
    actual.add_argument("--preflight-sha", required=True)
    cli = parser.parse_args()
    if cli.stage == "prepare":
        prepare()
    else:
        replay(cli.protocol_sha, cli.preflight_sha)
