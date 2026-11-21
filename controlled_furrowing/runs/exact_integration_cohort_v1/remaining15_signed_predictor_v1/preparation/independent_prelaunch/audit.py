from contextlib import ExitStack
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(ROOT / "scripts"))
import exact_signed_cohort_predictor as J

SOURCE_PIN = "5a7e6bdffa829335d1c88129d21f73f284a886053ab85d7bb4dbad3c8a6d2b6f"
PROTOCOL_PIN = "314b71743c1f90cdeb7bd1e8b54f312a1b2c6b264eef08ba8a9421c00d3f79ac"
PREFLIGHT_PIN = "ff5d044978df37ced65ab16cb599790ef87b5d8bb8b038c97987515f9c0e91fc"
GLOBAL_PIN = "ac2e091abaf8175d3d50b952b30bda8ebfcfed1cdb73ae3297bb8813df62ae4c"
INPUTS_PIN = "f6414069a462400cf32b55ddad59fce4625a06a1231b45027d6af95c643fb176"
RESULTS_PIN = "a30067ae17c8d85418bf6002822a665a1737bfd6103aa92fa7ceef2c785e450e"


def read(path, digest):
    data = path.read_bytes()
    assert sha256(data).hexdigest() == digest, path
    return json.loads(data)


def rejected(action):
    try:
        action()
    except ValueError:
        return True
    return False


def pair(value):
    a, b = map(F, value)
    assert a <= b
    return a, b


def serial_digest(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                             allow_nan=False).encode()).hexdigest()


def run():
    assert J.sha(J.SOURCE) == J.sha(J.PREP / "exact_signed_cohort_predictor_frozen.py") == SOURCE_PIN
    protocol = read(J.PREP / "protocol.json", PROTOCOL_PIN)
    preflight = read(J.PREP / "preflight.json", PREFLIGHT_PIN)
    assert preflight["status"] == "PASS" and len(preflight["checks"]) == 20
    assert all(value is True for value in preflight["checks"].values())
    assert preflight["actual_cached_sites_assembled"] == preflight["added_exact_integral_component_jobs"] == 0
    assert protocol["expected_finite_cache_site_count"] == 9504
    assert protocol["actual_cached_site_assembly_executed"] is False
    assert not J.ACTUAL.exists()
    pins = protocol["input_hashes"]
    assert len(pins) == 26
    for path, digest in pins.items():
        assert J.sha(ROOT / path) == digest, path
    assert pins[J.REVIEW] == GLOBAL_PIN and pins[J.INPUTS] == INPUTS_PIN and pins[J.PRODUCER] == RESULTS_PIN
    review = read(ROOT / J.REVIEW, GLOBAL_PIN)
    payload = read(ROOT / J.INPUTS, INPUTS_PIN)
    producer = read(ROOT / J.PRODUCER, RESULTS_PIN)
    terminal = read(ROOT / J.PRODUCER_TERMINAL, pins[J.PRODUCER_TERMINAL])
    observed = read(ROOT / J.OS_EXIT, pins[J.OS_EXIT])
    signed_fields = ["point_bulk_RF_intervals_exact", "retained_surface_RF_intervals_exact",
                     "RH_bulk_S_signed_intervals_exact", "closed_five_offsets_exact", "parameter_order"]
    for key in signed_fields:
        assert payload[key] == producer[key]
    assert payload["reducer_exact_array_digest"] == review["exact_array_sha256"]
    assert all(serial_digest(producer[k]) == h for k, h in review["exact_array_sha256"].items())
    assert review["producer_gain_contains_exact_rational_nodal_gain_slots"] == 1656
    assert review["producer_gain_overlapping_only_slots"] == 0
    assert payload["raw_gain_operand_proven_to_cover_exact_box"] is True
    assert payload["noncontaining_raw_gain_diagnostics"] == []

    targets = [(J.P, "build_point"), (J.P, "variation"), (J.P.core, "assemble_point"), (J.P.core, "local_terms")]
    with ExitStack() as guards:
        for module, name in targets:
            guards.enter_context(patch.object(module, name, side_effect=AssertionError(name)))
        bundle = J.input_bundle()
        assert json.loads(json.dumps(J.prep_doc(bundle), default=str)) == protocol
        links, ctx, (bulk, surf, gain, bulk_s, offsets), args, old_z, r, b, z, norm, provider = bundle
        curvature = J.D._minimum(ctx, dict(links), old_z, z, r, norm, provider)
        assert curvature["strict_local_minimum"] is True
        assert F(curvature["normalized_curvature_lower_exact"]) == (1-z)/norm > 0
        assert J.precision.prec == 640
        restored, u, R = J.P.core.restore_cache(read(ROOT / ctx.row["point_cache_file"], pins[ctx.row["point_cache_file"]]))
        assert J.np.array_equal(R, ctx.R) and J.np.array_equal(u, ctx.u)
        assert len(restored.triangles) == 144 and len(restored.shape) == 64
        assert len(restored.surfaces) == 2
        assert [len(v[0]) * len(v[1]) for v in restored.surfaces] == [144, 144]

    assert ctx.record_id == protocol["record_id"] == "optoGEF2_small_col_D"
    assert args["point_certificate_file"] == ctx.row["point_cache_file"]
    assert Path(args["predictor_file"]) == (ROOT / ctx.row["predictor_file"]).resolve()
    assert F(args["dose_exact"]) * ctx.base["common_gain"] == ctx.drive
    assert args["base_parameters_exact"] == ctx.row["base_parameters_exact"]
    assert z == old_z + b and 0 <= z < 1 and r == F(1, 10**11)
    assert tuple(payload["parameter_order"]) == J.P.NAMES == J.C.ORDER
    assert offsets == [(ctx.cell[k][0] - ctx.base[k], ctx.cell[k][1] - ctx.base[k]) for k in J.P.NAMES]
    assert len(bulk) == len(surf) == len(gain) == len(bulk_s) == 611
    assert all(len(v) == 5 for v in bulk_s)
    assert [*map(pair, payload[J.GAIN])] == gain
    assert [*map(pair, payload["verified_producer_gain_RF_intervals_exact"])] != gain
    assert all(x.denominator.bit_count() == 1 for row in (*bulk, *surf, *gain, *(v for row in bulk_s for v in row)) for x in row)

    for key, change in (
        ("same_closed_C_sha256", "f"*64), ("record_id", "optoGEF2_small_col_C"),
        ("same_radius_exact", str(r/2)), ("parameter_order", list(reversed(J.P.NAMES))),
    ):
        assert rejected(lambda: J.parse_operand({**payload, key: change}, ctx))
    wrong = deepcopy(payload)
    wrong[J.GAIN][600][1] = "1/3"
    assert rejected(lambda: J.parse_operand(wrong, ctx))
    wrong = deepcopy(payload)
    wrong[signed_fields[2]][610][4] = ["0", "0"]
    assert rejected(lambda: J.check_authority(review, observed, terminal, producer, wrong, ctx.record_id))
    assert rejected(lambda: J.check_authority(review, {**observed, "observed_exit_code": False},
                                              terminal, producer, payload, ctx.record_id))
    assert rejected(lambda: J.check_authority({**review, "producer_gain_overlapping_only_slots": 1},
                                              observed, terminal, producer, payload, ctx.record_id))
    assert rejected(lambda: J.check_authority(review, observed, {**terminal, "pid": terminal["pid"]+1},
                                              producer, payload, ctx.record_id))

    # Exact corner enumeration gives an independent interval oracle for the signed combination.
    h = F(1, 16)
    offset = [(-h, 2*h), (-2*h, 3*h), (0, h), (-h, 0), (-3*h, 2*h)]
    cached = [[(1, 2), (-4, -2), (2, 3), (-3, 2), (-9, -8)]]
    bulk_s = [[(-3, -2), (1, 2), (-1, 3), (-2, 1), (4, 5)]]
    gain = [(4, 5)]
    slope = J.signed.combine_slopes(cached, bulk_s, gain, 4)
    pooled = [(-2,0), (-3,0), (1,6), (-5,3), (-1,2)]
    assert slope == [pooled]
    values = [a + e + sum(j*d for j,d in zip(s, t)) for a,e in product((h, 3*h), (-h,h))
              for s in product(*pooled) for t in product(*offset)]
    exact = (min(values), max(values))
    assert J.signed.affine_rows([(h, 3*h)], slope, offset, [h]) == [exact] and exact == (-F(11,8), F(27,16))
    assert J.signed.combine_slopes([[(0, 0)]*4 + [(-5, -5)]],
                                  [[(0, 0)]*4 + [(2, 2)]], [(3, 3)], 4)[0][4] == (0, 0)
    cases = (
        ((-F(1, 4), F(1, 3)), True, False),
        ((-F(1, 2), F(1, 2)), False, False),
        ((-F(5, 2), -F(3, 2)), False, False),
        ((-F(5, 2), -F(25, 16)), False, True),
        ((F(25, 16), F(5, 2)), False, True),
    )
    for row, exists, excluded in cases:
        result = J.classify([row, (-F(1,4), F(1,4))], F(1,2), F(1))
        assert result["sufficient_existence_test_passes"] is exists
        assert result["sufficient_old_cube_exclusion_passes"] is excluded
    assert not J.ACTUAL.exists()

    return {
        "status": "PASS",
        "record_id": ctx.record_id,
        "source_sha256": SOURCE_PIN,
        "protocol_sha256": PROTOCOL_PIN,
        "preflight_sha256": PREFLIGHT_PIN,
        "global_reduction_review_sha256": GLOBAL_PIN,
        "independent_signed_inputs_sha256": INPUTS_PIN,
        "global_reduction_results_sha256": RESULTS_PIN,
        "audit_source_sha256": J.sha(__file__),
        "checks": {
            "all_26_byte_pins_same_record_and_producer_observed_OS_link": True,
            "exact_611_rows_and_5_columns_original_closed_offsets": True,
            "independent_exact_gain_used_after_lossless_authority_checks": True,
            "selected_cache_embeds_exact_same_U_and_R_and_explicit_selected_S": True,
            "Arb640_before_cache_and_context_with_explicit_D_point_arguments": True,
            "no_L_inputs_in_selected_pure_helper_call_graph": True,
            "signed_affine_arithmetic_corner_oracle_and_exact_gain_cancellation": True,
            "whole_box_Z_and_original_radius_fixed_with_cached_SPD_premise": True,
            "foreign_cell_wrong_terminal_overlap_only_malformed_row_rejected": True,
            "strict_existence_and_exclusion_boundaries_and_negative_witness": True,
            "all_actual_point_variation_local_terms_and_assembly_guarded": True,
        },
        "binding_scalars": {"radius_exact": str(r), "cached_Z_upper_exact": str(old_z),
                            "whole_box_defect_b_exact": str(b), "whole_box_Z_upper_exact": str(z),
                            "cached_SPD_support": curvature},
        "scope": [
            "Preparation only. No actual cached-site sum, integration, or solve was performed.",
            "Explicit D point cache, selected S, dose, and base are passed to the reviewed helper. Its actual variation uses radius zero on the predictor segment with full parameter cell, then supplies componentwise nonnegative remainder bounds.",
            "At fixed Y0 the cached force is affine in the five parameters. The state-segment remainder uses the Hessian difference from the anchor and offsets S delta. Independently enclosed bulk RH*S and affine gain force slopes add before each shared offset; retained surface RF is a separate whole-box interval.",
            "Whole-cube derivative z = cached z plus full derivative defect b; it is distinct from predictor-only signed RF. At a root inside the convex matched cube the predictor force is at most (1+z)r by integrating RH along the state segment.",
            "A strict passing inequality gives a candidate only. Reviewed actual replay remains required. Scope is the original D closed cell and numerical predictor-matched cubes."
        ],
    }


if __name__ == "__main__":
    result = run()
    with (J.PREP / "prelaunch_independent_review.json").open("x", encoding="utf-8") as output:
        json.dump(result, output, indent=2)
        output.write("\n")
