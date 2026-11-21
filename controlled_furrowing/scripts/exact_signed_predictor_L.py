"""Signed full exact-FE predictor force on the frozen record-L original cube."""
import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
import os
from pathlib import Path
from time import monotonic
from unittest.mock import patch

import numpy as np
from flint import arb, ctx as precision

import parametric_root_cells as P
import review_cohort_exact_record as D
import exact_cohort_context as C


ROOT = Path(__file__).resolve().parents[1]
BASE = "runs/exact_integration_cohort_v1/L_defect_diagnostic"
STAGE = ROOT / BASE / "signed_predictor_L_v1"
PREP = STAGE / "preparation"
CENTERED = BASE + "/centered_bulk_L_v1"
RECORD = "optoGEF2_small_col_L"
ORIGINAL = D._record_dir(RECORD)
PIN = {
    "scripts/parametric_root_cells.py": "362a19548aa60c4234017316e22deef4a80e93b73c725f47132573a1a9103b60",
    "scripts/algebraic_root_feasibility.py": "7321d38bf17cf354cc03be2871f3fbaa45bcebbcbfca7be90d0f867b80078acf",
    "scripts/exact_cohort_context.py": "c22612a1de1ccfd34e8b0f0859f6413c4da09c765ba4945d060340e93e928f9f",
    "scripts/review_cohort_exact_record.py": "67db748a5d40f31471ddd3d6c84fae280ff7ead650e6eb4e950a4766e54e0285",
    ORIGINAL + "/protocol.json": "66bfdb7dbe49d3ece787e7ecc9f16e004f7611a8b946be0e078f8a6806a4c937",
    CENTERED + "/actual_independent_review.json": "89df6cbaaaae0237b07e39c7cef05fbbeba20117d73a67153a4adcb459b6d6d0",
    CENTERED + "/reduction_preparation/results.json": "5f07f0aea44bc176d615de9e593e16695cc478082b7cf9fc9847aada2a0604e8",
    CENTERED + "/reduction_preparation/actual_independent_review.json": "cce6f01147b6ed436e3f3a8620e94bbd9df84a8fc1aa45e9bbed9cba8c1f8c4a",
    CENTERED + "/reduction_preparation/input_admission.json": "1b0b09ed4e11746eba26485f3296dc9201a68e6fcb13c98d64ff78335374deb3",
    CENTERED + "/reduction_preparation/terminal.json": "e6d0e6e35ee1457adbd9b6519ba6cd65ac0bc2aaaef4a472a722c48cef4dab33",
}
USED_ARRAYS = (
    "separated_point_bulk_RF_intervals_exact",
    "separated_retained_surface_RF_intervals_exact",
    "coupled_RH_bulk_S_intervals_exact",
    "separated_gain_derivative_RF_intervals_exact",
    "exact_closed_parameter_offsets",
)


def _relative(path):
    return Path(path).relative_to(ROOT).as_posix()


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _bound(path, hashes, parse=True):
    return D._bound(path, PIN[path], hashes, parse=parse)


def _interval(pair, dyadic=False):
    lo, hi = map(D._fraction, pair)
    D.require(lo <= hi, "finite ordered closed exact interval required")
    if dyadic:
        D._exp(lo)
        D._exp(hi)
    return lo, hi


def _add(*values):
    return sum((a for a, _ in values), Q(0)), sum((b for _, b in values), Q(0))


def _mul(a, b):
    products = [x * y for x in a for y in b]
    return min(products), max(products)


def affine_rows(base, slopes, offsets, remainder):
    """Exact signed affine enclosure, then rowwise symmetric nonnegative remainder."""
    D.require(len(base) == len(slopes) == len(remainder) and bool(base)
              and len(offsets) == 5 and all(len(row) == 5 for row in slopes),
              "complete five-parameter indexed residual rows required")
    D.require(min(remainder) >= 0 and all(lo <= 0 <= hi for lo, hi in offsets),
              "nonnegative remainder and base-in-cell exact shifts required")
    return [_add(a, (-e, e), *(_mul(j, d) for j, d in zip(row, offsets)))
            for a, row, e in zip(base, slopes, remainder)]


def combine_slopes(cached, bulk, gain, gain_index):
    D.require(len(cached) == len(bulk) == len(gain) and 0 <= gain_index < 5
              and all(len(a) == len(b) == 5 for a, b in zip(cached, bulk)),
              "matching cached and defect slope columns required")
    return [
        [_add(a, b, g if j == gain_index else D.ZERO) for j, (a, b) in enumerate(zip(c, h))]
        for c, h, g in zip(cached, bulk, gain)
    ]


def strict_bound(eta, z, radius):
    D.require(eta >= 0 and z >= 0 and radius > 0, "finite nonnegative sufficient-test scalars and positive radius required")
    lhs = eta + z * radius
    passed = z < 1 and lhs < radius
    return lhs, passed, eta / (1 - z) if passed else None


def _inputs():
    hashes = {}
    for module, p in ((P, "scripts/parametric_root_cells.py"), (P.core, "scripts/algebraic_root_feasibility.py"),
                      (D, "scripts/review_cohort_exact_record.py"), (C, "scripts/exact_cohort_context.py")):
        _bound(p, hashes, parse=False)
        D.require(_sha(module.__file__) == hashes[p], "wrong imported frozen source")
    protocol = _bound(ORIGINAL + "/protocol.json", hashes)
    ctx = C.load_context(RECORD)
    D.require(ctx.row == protocol["row"] and tuple(P.NAMES) == tuple(C.ORDER)
              and ctx.R.shape == (611, 611) and ctx.u.shape == (611,) and ctx.S.shape == (611, 5)
              and np.isfinite(ctx.R).all() and np.isfinite(ctx.S).all(),
              "same authoritative old cached point, R, S and ordered material coordinates required")
    for name, digest in ctx.input_hashes.items():
        D._bound(name, digest, hashes, parse=False)
    box = D._bound(protocol["coordinate_box_file"], protocol["coordinate_box_sha256"], hashes)
    canonical = C.exact_node_box(ctx, ctx.old_radius)
    D.require(all(D._digest(box[k]) == D._digest(v) for k, v in canonical.items()),
              "unchanged original declared cube and full C required")
    old_eta, old_Z, r, old_provider = D._old_premises(ctx, protocol, hashes)
    resultpath = CENTERED + "/reduction_preparation/results.json"
    final_review = _bound(CENTERED + "/reduction_preparation/actual_independent_review.json", hashes)
    component_review = _bound(CENTERED + "/actual_independent_review.json", hashes)
    admission = _bound(CENTERED + "/reduction_preparation/input_admission.json", hashes)
    terminal = _bound(CENTERED + "/reduction_preparation/terminal.json", hashes)
    D._terminal(terminal, PIN[resultpath])
    result = _bound(resultpath, hashes)
    D.require(final_review["status"] == "PASS_ACTUAL_CENTERED_ORIGINAL_CUBE_REDUCTION_REVIEW"
              and component_review["status"] == "PASS_ACTUAL_CENTERED_PACKET_AND_COMPOSITION_REPLAY"
              and D._yes(final_review["checks"]) and D._yes(component_review["checks"])
              and final_review["component_replay_sha256"] == component_review["support_sha256"]
              and final_review["reducer_results_sha256"] == PIN[resultpath]
              and final_review["admission_sha256"] == PIN[CENTERED + "/reduction_preparation/input_admission.json"]
              and final_review["terminal_sha256"] == PIN[CENTERED + "/reduction_preparation/terminal.json"]
              and len(admission["expected_centered_packets"]) == 144 and result["component_count"] == 168
              and result["retained_surface_components"] == 24 and result["same_point_cell_and_cube"] is True
              and result["force_entries"] == 611 and result["upper_Hessian_entries"] == 6676
              and result["predictor_parameter_order"] == list(P.NAMES), "same reviewed centered signed-defect arrays required")
    D._bound(component_review["support_file"], component_review["support_sha256"], hashes, parse=False)
    for key in USED_ARRAYS:
        D.require(D._digest(result[key]) == component_review["expected_reducer_array_digests"][key],
                  f"centered array differs from independent endpoint replay: {key}")
    for f in (0, 1):
        for i in range(12):
            path = f"{ORIGINAL}/surface/face{f}_edge{i:02d}.json"
            D._bound(path, result["source_inputs_and_packets"][path], hashes, parse=False)
    b = Q(result["b_upper_exact"])
    total_z = old_Z + b
    D.require(Q(result["normal_radius_exact"]) == r and Q(result["cached_Z_upper_exact"]) == old_Z
              and Q(result["cached_eta_upper_exact"]) == old_eta
              and Q(result["z_upper_exact"]) == total_z and 0 <= total_z < 1,
              "same original radius and reviewed full-box derivative bound required")
    bulk, surface, gain = [[_interval(v, dyadic=True) for v in result[key]] for key in (
        USED_ARRAYS[0], USED_ARRAYS[1], USED_ARRAYS[3]
    )]
    bulk_S = [[_interval(v, dyadic=True) for v in row] for row in result[USED_ARRAYS[2]]]
    delta = [_interval(v) for v in result[USED_ARRAYS[4]]]
    cell = P.check_cell(ctx.row["parameter_cell"])
    D.require(len(bulk) == len(surface) == len(gain) == len(bulk_S) == 611 and all(len(row) == 5 for row in bulk_S)
              and len(delta) == 5 and cell == ctx.cell
              and delta == [(cell[key][0] - ctx.base[key], cell[key][1] - ctx.base[key]) for key in P.NAMES]
              and all(a <= 0 <= b for a, b in delta) and ctx.drive == ctx.base["common_gain"] * ctx.dose,
              "same selected exact offsets, dose and Signed RF arrays required")
    pinned_arguments = {
        "point_certificate_file": ctx.row["point_cache_file"],
        # The old builder joins P.OUT with this value. An absolute reference
        # names the one selected affine predictor file without changing globals.
        "predictor_file": str((ROOT / ctx.row["predictor_file"]).resolve()),
        "dose_exact": str(ctx.dose), "base_parameters_exact": dict(ctx.row["base_parameters_exact"]),
    }
    return hashes, ctx, protocol, result, bulk, surface, gain, bulk_S, delta, b, total_z, r, old_Z, old_provider, pinned_arguments


def _arb_bounds(value):
    return D._packed(P.core.pack(value))


def _point_and_remainder(context, arguments, offsets):
    point = P.build_point(arguments)
    D.require(point["n"] == 611 and np.array_equal(point["S"], context.S)
              and len(point["cache"]) == 9504
              and Counter(x["kind"] for x in point["cache"]) == {"bulk": 9216, "surface": 288}
              and point["base"] == [context.base[k] for k in P.NAMES]
              and D._contains(_arb_bounds(point["dose"]), (context.dose, context.dose))
              and _arb_bounds(point["Z0"])[1] < 1, "same frozen decoded point object, S and nonsingular anchor required")
    cells = P.cell_context(point, context.row["parameter_cell"], "affine")
    D.require(all(D._contains(_arb_bounds(ball), exact) for ball, exact in zip(cells["delta"], offsets)),
              "cached variation's Arb cell must enclose the exact offsets")
    remainder, min_J = P.variation(point, cells, arb(0), True)
    left = [_arb_bounds(point["gR"][i, 0]) for i in range(611)]
    slopes = [[_arb_bounds(point["affine_M"][i, a]) for a in range(5)] for i in range(611)]
    bounds = [_arb_bounds(x)[1] for x in remainder]
    D.require(len(bounds) == 611 and min(bounds) >= 0 and _arb_bounds(min_J)[0] > Q(1, 10**6),
              "full valid cached predictor path and remainder required")
    return left, slopes, bounds, _arb_bounds(min_J), _arb_bounds(point["Z0"])


def tests():
    checks = {}
    deltas = [(Q(-1), Q(2)), (Q(-2), Q(0)), (Q(0), Q(1)), (Q(-1), Q(1)), (Q(-1), Q(2))]
    g = [(Q(2), Q(2))]
    point = [(Q(-2), Q(-2))]
    cache = [[(Q(x), Q(x)) for x in (1, -1, 2, -2, 3)]]
    defect = [[(Q(x), Q(x)) for x in (0, 2, -2, 1, -4)]]
    gain = [(Q(1), Q(1))]
    slopes = combine_slopes(cache, defect, gain, 4)
    checks["shared_gain_slope_cancels_exactly"] = slopes[0][4] == D.ZERO
    bounds = affine_rows([_add(g[0], point[0])], slopes, deltas, [Q(1, 8)])[0]
    corners = [sum(v[0] * q for v, q in zip(slopes[0], corner)) for corner in product(*deltas)]
    checks["mixed_asymmetric_five_coordinate_corners_and_sym_remainder"] = (
        bounds == (min(corners) - Q(1, 8), max(corners) + Q(1, 8)))
    checks["midpoint_is_not_used_for_interval_slopes"] = _mul((Q(-2), Q(3)), (Q(-4), Q(1))) == (Q(-12), Q(8))
    eps = Q(1, 1024)
    direct = affine_rows([_add((Q(1), Q(1) + eps), (Q(-1), -1 + eps))],
                         [[D.ZERO] * 5], deltas, [eps])[0]
    checks["signed_cancellation_remains_valid_after_outward_endpoint_rounding"] = (
        direct == (Q(-1, 1024), Q(3, 1024)) and strict_bound(D._abs_upper(direct), Q(1, 2), Q(1))[1])
    checks["triangle_norm_bound_can_fail_while_signed_passes"] = not strict_bound(Q(2), Q(1, 2), Q(1))[1]
    checks["strict_selfmap_equality_is_unresolved"] = strict_bound(Q(1, 8), Q(3, 4), Q(1, 2)) == (Q(1, 2), False, None)
    checks["strict_pass_yields_exact_radius_error"] = strict_bound(Q(1, 8), Q(1, 2), Q(1))[-2:] == (True, Q(1, 4))
    for name, action in (
        ("mismatched_five_column_predictor", lambda: affine_rows([D.ZERO], [[D.ZERO] * 4], deltas, [Q(0)])),
        ("negative_scalar_remainder", lambda: affine_rows([D.ZERO], [[D.ZERO] * 5], deltas, [Q(-1)])),
        ("disjoint_base_gain_offsets", lambda: affine_rows([D.ZERO], [[D.ZERO] * 5],
                                                        deltas[:4] + [(Q(1), Q(2))], [Q(0)])),
        ("nondyadic_bound_from_claimed_Arb", lambda: _interval(("1/3", "1/2"), dyadic=True)),
        ("wrong_gain_column", lambda: combine_slopes(cache, defect, gain, 5)),
    ):
        try:
            action()
        except ValueError:
            checks[name] = True
        else:
            checks[name] = False
    D.require(D._yes(checks), "all exact signed predictor arithmetic checks required")
    return checks


def prepare():
    forbidden = AssertionError("signed predictor preparation entered actual cached quadrature")
    with (
        patch.object(P, "build_point", side_effect=forbidden),
        patch.object(P, "variation", side_effect=forbidden),
        patch.object(P.core, "assemble_point", side_effect=forbidden),
        patch.object(P.core, "local_terms", side_effect=forbidden),
    ):
        hashes, c, _, result, _, _, _, _, delta, b, z, r, old_Z, _, point = _inputs()
        checks = tests()
    checks.update({"all611_same_record_5_column_arrays": True, "signed_arrays_match_independent_endpoint_replay": True,
                   "old_contraction_and_original_radius_retained": z == old_Z + b and r == c.old_radius,
                   "same_base_exact_cell_and_unchanged_24_surface_endpoint_packets": True,
                   "no_actual_9504_cached_sites_assembled_or_varied": True})
    D.require(D._yes(checks), "frozen signed-predictor prerequisites required")
    names = ("exact_signed_predictor_L_frozen.py", "protocol.json", "preflight.json")
    D.require(not any((PREP / n).exists() for n in names), "preserve write-once signed predictor preparation")
    PREP.mkdir(parents=True, exist_ok=True)
    with (PREP / names[0]).open("xb") as output:
        output.write(Path(__file__).read_bytes())
    p = {
        "status": "FROZEN_SIGNED_ORIGINAL_L_PREDICTOR_PREPARATION_REVIEW_PENDING", "record_id": RECORD,
        "source_file": _relative(__file__), "source_sha256": _sha(__file__), "input_hashes": hashes,
        "precision_bits": 640, "base_parameters_exact": c.base, "point_kernel_arguments": point,
        "parameter_cell_sha256": c.row["parameter_cell_sha256"], "parameter_order": list(P.NAMES),
        "same_exact_closed_offsets": delta, "radius_exact": r, "unchanged_z_upper_exact": z,
        "previous_separate_norm_eta_exact": result["eta_upper_exact"], "strict_signed_force_budget_exact": r * (1 - z),
        "scope": "Full signed normalized integral force only at p(theta) on the original closed C. The full-box "
                 "derivative bound z is unchanged. No new point integration such as F_I, material fit, or nonlinear solve.",
        "cached_force": "Rebuild a rigorous 9504-term cached-Q point using the exact frozen point cache and selected "
                        "affine predictor S. Take gR and affine_M from audited PR.build_point. "
                        "Use its affine variation with radius zero for a rigorous rowwise predictor-only remainder.",
        "cached_path_derivation": "At fixed Y0, cached force is affine in the five material parameters. Move theta "
                                  "there first, then vary Y along the predictor at that theta. Subtract H_Q0*S*delta "
                                  "from that nodal segment; audited variation bounds the resulting H_Q-H_Q0 remainder.",
        "signed_composition": "Add cached gR, signed point-bulk RF and retained whole-box surface RF rowwise. "
                              "Add cached affine M, RH_bulk*S and RG in the same gain column. Multiply once by the "
                              "same five closed parameter shifts before taking any row magnitude. Add the cached remainder symmetrically.",
        "actual_cached_point_or_variation_run": False,
        "failed_attempts": {_relative(p): _sha(p) for p in sorted(PREP.glob("development_revision_*/*")) if p.is_file()},
    }
    psha = D._write(PREP / names[1], p)
    report = {"status": "PASS_SIGNED_PREDICTOR_PREPARATION_ONLY", "source_sha256": p["source_sha256"],
              "protocol_sha256": psha, "checks": checks, "actual_cached_9504_site_sums": False}
    presha = D._write(PREP / names[2], report)
    print(json.dumps({"status": report["status"], "source_sha256": p["source_sha256"], "protocol_sha256": psha,
                      "preflight_sha256": presha, "checks": len(checks)}), flush=True)


def run(review_sha):
    gate = D._bound(_relative(PREP / "prelaunch_independent_review.json"), review_sha, {})
    current = {"source_sha256": _sha(__file__), "protocol_sha256": _sha(PREP / "protocol.json"),
               "preflight_sha256": _sha(PREP / "preflight.json")}
    D.require(gate["status"] == "PASS" and D._yes(gate["checks"]) and all(gate[k] == v for k, v in current.items()),
              "bound independent signed predictor prelaunch required")
    protocol = json.loads((PREP / "protocol.json").read_bytes())
    hashes, c, _, reviewed, bulk, surface, gain, bulk_S, delta, b, z, r, old_Z, provider, args = _inputs()
    D.require(hashes == protocol["input_hashes"] and protocol["point_kernel_arguments"] == args
              and D._digest(delta) == D._digest(protocol["same_exact_closed_offsets"])
              and Q(protocol["unchanged_z_upper_exact"]) == z, "frozen signed predictor inputs changed")
    precision.prec = 640
    P.core.A.cache_clear()
    import torch
    torch.set_num_threads(1)
    out = STAGE / "actual"
    out.mkdir(parents=True, exist_ok=True)
    utc = lambda: datetime.now(timezone.utc).isoformat()
    start = monotonic()
    D._write(out / "started.json", {**current, "prelaunch_sha256": review_sha, "pid": os.getpid(), "utc": utc()})
    try:
        g, cached_M, rem, min_J, Z0 = _point_and_remainder(c, args, delta)
        origins = [_add(x, y, t) for x, y, t in zip(g, bulk, surface)]
        slopes = combine_slopes(cached_M, bulk_S, gain, P.NAMES.index("common_gain"))
        intervals = affine_rows(origins, slopes, delta, rem)
        eta = max(map(D._abs_upper, intervals))
        lhs, passed, error = strict_bound(eta, z, r)
        norm = Q(reviewed["R_infinity_norm_exact"])
        candidate = D._minimum(c, hashes, old_Z, z, r, norm, provider) if passed else None
        worst = max(range(611), key=lambda i: D._abs_upper(intervals[i]))
        result = {
            **current, "prelaunch_sha256": review_sha,
            "status": "CONDITIONAL_UNREVIEWED_SIGNED_PREDICTOR_SUFFICIENT_TEST_PASSES" if passed
                      else "UNREVIEWED_SIGNED_PREDICTOR_SUFFICIENT_TEST_UNRESOLVED",
            "record_id": RECORD, "scope": "Signed RF_I at parameter-matched original predictors only; full-box z reviewed previously.",
            "same_closed_C_sha256": c.row["parameter_cell_sha256"], "base_parameters_exact": c.base,
            "signed_cached_anchor_gR": g, "signed_cached_affine_M": cached_M,
            "nonnegative_cached_predictor_remainder_upper_exact": rem,
            "cached_predictor_path_J_lower_exact": min_J[0], "cached_anchor_Z0_upper_exact": Z0[1],
            "signed_exact_integral_anchor_including_wholebox_surface": origins, "signed_combined_affine_columns": slopes,
            "same_closed_parameter_offsets": delta, "signed_full_integral_RF_predictor_intervals": intervals,
            "worst_row_index_zero_based": worst, "eta_I_upper_exact": eta, "b_upper_exact": b,
            "cached_Z_upper_exact": old_Z, "z_upper_exact": z, "r_exact": r,
            "strict_eta_budget_exact": r * (1 - z), "selfmap_lhs_exact": lhs, "previous_separate_norm_eta_exact": reviewed["eta_upper_exact"],
            "sufficient_test_pass": passed, "conditional_root_error_exact": error, "conditional_curvature": candidate,
            "independently_reviewed_root_assertion": False,
            "shared_primitives": ["frozen parametric cached-Q exact point assembly", "audited cached predictor variation",
                                  "independently reviewed centered bulk and retained whole-box surface endpoint arrays"],
            "input_hashes": hashes, "seconds": monotonic() - start,
        }
        digest = D._write(out / "results.json", result)
        D._write(out / "terminal.json", {"status": "FINISHED", "exit_code": 0,
                                        "results_sha256": digest, "pid": os.getpid(), "utc": utc()})
        print(json.dumps({"status": result["status"], "results_sha256": digest}), flush=True)
    except BaseException as error:
        D._write(out / "failure.json", {"error": str(error), **current})
        D._write(out / "terminal.json", {"status": "FAILED", "exit_code": 1, "pid": os.getpid(), "utc": utc()})
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest="command", required=True)
    command.add_parser("prepare")
    work = command.add_parser("run")
    work.add_argument("--review-sha256", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare()
    else:
        run(args.review_sha256)
