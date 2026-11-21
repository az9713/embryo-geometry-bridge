"""Bounded preparation tests using frozen inputs. Run with --freeze to write receipts."""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys
from time import monotonic
from unittest.mock import patch

import numpy as np


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import exact_cohort_context as C


def _frozen_candidate_checks(ctx, authority, checks):
    """Fault injection is confined to in-memory reads and document copies."""
    def reject(label, action, message=None):
        try:
            action()
        except ValueError as error:
            if message is not None and message not in str(error):
                raise AssertionError(f"{label}: wrong rejection") from error
            checks[label] = True
        else:
            raise AssertionError(f"{label}: accepted malformed input")

    row = ctx.row
    root = json.loads((ROOT / row["root_certificate_file"]).read_bytes())
    other_row = next(r for r in authority["records"] if r["record_id"] != row["record_id"])
    other_root = json.loads((ROOT / other_row["root_certificate_file"]).read_bytes())
    reject("cross_record_root", lambda: C._validate_root_links(row, other_root, authority["input_hashes"]),
           "root identity")
    for kind in ("point_cache", "predictor"):
        key = f"{kind}_file"
        bad = deepcopy(root)
        bad[key] = other_row[key]
        reject(f"cross_record_{kind}_path", lambda: C._validate_root_links(row, bad, authority["input_hashes"]),
               f"path/digest mismatch: {key}")
        bad = deepcopy(root)
        bad["artifact_hashes"][row[key]] = other_row[f"{kind}_sha256"]
        reject(f"root_artifact_{kind}_digest_mismatch",
               lambda: C._validate_root_links(row, bad, authority["input_hashes"]),
               f"path/digest mismatch: {key}")
    for file in [C.AUTHORITY] + [row[k] for k in (
            "root_certificate_file", "point_cache_file", "predictor_file", "proof_center_state_file")]:
        real_read = Path.read_bytes
        def changed(path):
            raw = real_read(path)
            return raw + b"\n" if path == ROOT / file else raw

        with patch.object(Path, "read_bytes", changed):
            reject(f"changed_authoritative_bytes:{file}", lambda: C.load_context(ctx.record_id),
                   f"changed frozen artifact: {file}")

    recentered = next(r for r in authority["records"] if r["record_id"] == "optoGEF2_small_col_F")
    recentered_root = json.loads((ROOT / recentered["root_certificate_file"]).read_bytes())
    bad = deepcopy(recentered_root)
    bad["proof_center_state_file"] = bad["original_seed_state_file"]
    reject("recentered_seed_substitution",
           lambda: C._validate_root_links(recentered, bad, authority["input_hashes"]),
           "proof-center path/digest mismatch")

    packet = deepcopy(ctx.target.packet)
    encoded = packet["frozen_binary_cache"]["data"]
    packet["frozen_binary_cache"]["data"] = ("B" if encoded[0] != "B" else "A") + encoded[1:]
    state = json.loads((ROOT / row["proof_center_state_file"]).read_bytes())
    reject("embedded_cache_digest_mismatch", lambda: C._target(packet, row, state),
           "embedded cache bytes or digest changed")
    bad_state = deepcopy(state)
    bad_state["free_displacement"][0] += 1.0
    reject("embedded_displacement_not_from_proof_center",
           lambda: C._target(ctx.target.packet, row, bad_state),
           "cache does not match proof center")
    reject("unknown_exact_record", lambda: C.load_context("L"), "one authoritative record")
    reject("wrong_case_exact_record", lambda: C.load_context(ctx.record_id.lower()),
           "one authoritative record")
    reject("duplicate_exact_record", lambda: C._one_record([row, row], ctx.record_id),
           "one authoritative record")
    reject("boolean_record_id", lambda: C.load_context(True), "record ID string")

    for label, value in (
        ("negative", -ctx.old_radius), ("below_old", ctx.old_radius / 2), ("zero", 0),
        ("boolean", True), ("numpy_boolean", np.bool_(False)), ("infinity", float("inf")),
        ("nan", float("nan")), ("zero_denominator", "1/0"), ("nonnumeric", "old"),
    ):
        reject("radius_" + label, lambda: C.exact_node_box(ctx, value))
    for label, value in (
        ("short", [0] * 610), ("long", [0] * 612), ("column", np.zeros((611, 1))),
        ("scalar", Fraction(0)), ("text", "0" * 611), ("wrong_rank", [[0, 0]] * 611),
        ("boolean_entry", [True] + [0] * 610), ("nan_entry", [float("nan")] + [0] * 610),
        ("inf_entry", [float("-inf")] + [0] * 610), ("nonnumeric_entry", ["shift"] + [0] * 610),
    ):
        reject("shift_" + label, lambda: C.exact_node_box(ctx, ctx.old_radius, shift=value))


def _verify_boxes(contexts, checks):
    def check(label, condition):
        if not condition:
            raise AssertionError(label)
        checks[label] = True

    records = []
    for ctx in contexts:
        tag = ctx.record_id
        box = C.exact_node_box(ctx, ctx.old_radius)
        root = json.loads((ROOT / ctx.row["root_certificate_file"]).read_bytes())
        uniform = root["uniform_cell_result"]
        old_bounds = list(map(Fraction, uniform["full_nodal_predictor_bounds_exact"]))
        old_outer = Fraction(uniform["uniform_outer_nodal_radius_about_base"]["upper_exact"])
        old_root_error = Fraction(uniform["uniform_error_from_predictor"]["upper_exact"])
        ref = box["reference_nodes_exact"]
        centers = box["current_centers_exact"]
        radii = box["coordinate_radii_exact"]
        hulls = box["predictor_coordinate_hulls_exact"]
        offsets = box["predictor_hulls_exact"]
        free_map = {int(dof): i for i, dof in enumerate(ctx.free)}
        displacement = [Fraction.from_float(float(value)) for value in ctx.u]
        endpoint_offsets = {
            key: tuple(Fraction(bound) - Fraction(ctx.row["base_parameters_exact"][key])
                       for bound in ctx.row["parameter_cell"][key]) for key in C.ORDER
        }

        check(tag + ": cell/base/gain/dose verified exactly",
              all(tuple(map(Fraction, ctx.row["parameter_cell"][key])) == ctx.cell[key]
                  and Fraction(ctx.row["base_parameters_exact"][key]) == ctx.base[key]
                  and ctx.cell[key][0] <= ctx.base[key] <= ctx.cell[key][1] for key in C.ORDER)
              and ctx.dose == Fraction(root["dose_exact"])
              and ctx.drive == Fraction(root["base_strength_exact"]) == ctx.base["common_gain"] * ctx.dose
              and ctx.cached_binary_drive == Fraction.from_float(float(ctx.drive)))
        check(tag + ": original record dimensions and box scope",
              len(centers) == len(ref) == len(hulls) == len(radii) == 325
              and all(len(row) == 2 for rows in (ref, centers, radii, hulls) for row in rows)
              and len(offsets) == len(old_bounds) == len(free_map) == 611
              and len(box["triangles"]) == 144 and all(len(tri) == 6 for tri in box["triangles"])
              and box["corrector_radius_exact"] == ctx.old_radius
              and box["root_membership_assumed"] is False and box["record_id"] == tag
              and box["parameter_cell_sha256"] == C.CELL_SHA256)
        check(tag + ": old root error bound fits old selected corrector",
              0 <= old_root_error <= ctx.old_radius and max(old_bounds) + old_root_error <= old_outer)
        count = 0
        axis_count = 0
        for node, axis in product(range(325), range(2)):
            dof = 2 * node + axis
            center, r = centers[node][axis], radii[node][axis]
            minimum, maximum = hulls[node][axis]
            if dof not in free_map:
                check(f"{tag}: constrained coordinate {dof}", minimum == maximum == center == ref[node][axis]
                      and r == 0)
                axis_count += (axis == 0 and center == 0)
                continue
            i = free_map[dof]
            coefficients = [Fraction.from_float(float(v)) for v in ctx.S[i]]
            # Independent extremum witness: evaluate at the two signed corners of C.
            values = [
                sum(v * endpoint_offsets[key][side if v >= 0 else 1 - side]
                    for key, v in zip(C.ORDER, coefficients))
                for side in (0, 1)
            ]
            check(f"{tag}: free coordinate {dof}", tuple(values) == offsets[i]
                  and center == ref[node][axis] + displacement[i]
                  and minimum == center + values[0] and maximum == center + values[1]
                  and center - r <= minimum - ctx.old_radius
                  and maximum + ctx.old_radius <= center + r
                  and r == max(map(abs, values)) + ctx.old_radius
                  and max(map(abs, values)) <= old_bounds[i])
            count += 1

        geometry_bytes = json.dumps({
            key: box[key] for key in ("reference_nodes_exact", "triangles", "free")
        }, default=str, sort_keys=True, separators=(",", ":")).encode()
        records.append({
            "record_id": tag,
            "dose_exact": str(ctx.dose),
            "base_strength_exact": str(ctx.drive),
            "cached_binary_drive_exact": str(ctx.cached_binary_drive),
            "selected_old_radius_exact": str(ctx.old_radius),
            "largest_nodal_radius_at_old_r_exact": str(max(map(max, radii))),
            "exact_free_coordinate_enclosures_checked": count,
            "fixed_coordinate_enclosures_checked": 650 - count,
            "fixed_zero_axis_coordinates": axis_count,
            "geometry_sha256": sha256(geometry_bytes).hexdigest(),
            "proof_center_state_file": ctx.row["proof_center_state_file"],
        })
    check("all17 exact boxes use identical reference, connectivity and free map",
          len({row["geometry_sha256"] for row in records}) == 1)
    check("F and K both use distinct recentered proof-center states",
          all(ctx.row["proof_center_state_file"].endswith("/proof_center_state.json") for ctx in contexts
              if ctx.record_id in ("optoGEF2_small_col_F", "optoGEF2_small_col_K")))
    ctx = contexts[0]
    additive_shift = [Fraction(1, 10**8)] + [Fraction(0)] * 609 + [Fraction(-1, 10**10)]
    radius = ctx.old_radius + max(map(abs, additive_shift))
    shifted = C.exact_node_box(ctx, radius, shift=additive_shift)
    old = C.exact_node_box(ctx, ctx.old_radius)
    for i, dof in enumerate(ctx.free):
        node, axis = divmod(int(dof), 2)
        old_center = old["current_centers_exact"][node][axis]
        old_min, old_max = old["predictor_coordinate_hulls_exact"][node][axis]
        new_center = shifted["current_centers_exact"][node][axis]
        new_min, new_max = shifted["predictor_coordinate_hulls_exact"][node][axis]
        larger_r = shifted["coordinate_radii_exact"][node][axis]
        check(f"positive/negative shift and larger cube:{int(dof)}",
              new_center == old_center + additive_shift[i]
              and new_min == old_min + additive_shift[i] and new_max == old_max + additive_shift[i]
              and new_center - larger_r <= min(old_min, new_min) - ctx.old_radius
              and new_center + larger_r >= max(old_max, new_max) + ctx.old_radius)
    return records


def run(freeze):
    start = monotonic()
    authority = json.loads((ROOT / C.AUTHORITY).read_bytes())
    ids = [row["record_id"] for row in authority["records"]]
    if len(ids) != 17 or len(set(ids)) != 17:
        raise AssertionError("all 17 distinct frozen exact record IDs required")
    contexts = [C.load_context(record_id) for record_id in ids]
    checks = {}
    _frozen_candidate_checks(contexts[0], authority, checks)
    records = _verify_boxes(contexts, checks)

    # The F/K matrix archives carry seed lineage. Even reading them is unnecessary.
    for ctx in contexts[-2:]:
        archive_path = ROOT / ctx.target.packet["matrix_archive"]
        real_read = Path.read_bytes
        def guard_historical_archive(path):
            if path == archive_path:
                raise AssertionError("historical seed matrix archive accessed")
            return real_read(path)

        with patch.object(Path, "read_bytes", guard_historical_archive):
            C.load_context(ctx.record_id)
        checks[ctx.record_id + ": historical seed matrix archive never read"] = True

    source_path = Path(C.__file__).resolve()
    source = source_path.read_bytes()
    source_sha = sha256(source).hexdigest()
    suite_sha = sha256(Path(__file__).read_bytes()).hexdigest()
    report = {
        "status": "PASS_CONTEXT_PREPARATION_ONLY",
        "source_sha256": source_sha,
        "preflight_source_sha256": suite_sha,
        "elapsed_seconds": monotonic() - start,
        "check_count": len(checks),
        "checks": checks,
        "record_count": len(records),
        "records": records,
        "exact_integration_started": False,
        "matrix_assembly_started": False,
        "new_equilibrium_solves": False,
        "new_exact_FE_root_claims": False,
    }
    if freeze:
        names = ("exact_cohort_context_frozen.py", "protocol.json", "preflight.json")
        if any((OUT / name).exists() for name in names):
            raise ValueError("write-once preparation output exists")
        protocol = {
            "status": "FROZEN_CONTEXT_PREPARATION_REVIEW_PENDING",
            "version": 1,
            "utc": datetime.now(timezone.utc).isoformat(),
            "source_sha256": source_sha,
            "source_file": source_path.relative_to(ROOT).as_posix(),
            "source_snapshot": names[0],
            "preflight_source_file": Path(__file__).relative_to(ROOT).as_posix(),
            "preflight_source_sha256": suite_sha,
            "frozen_authority": {"file": C.AUTHORITY, "sha256": C.AUTHORITY_SHA256},
            "record_input_hashes": {ctx.record_id: ctx.input_hashes for ctx in contexts},
            "parameter_cell": authority["records"][0]["parameter_cell"],
            "parameter_cell_sha256": C.CELL_SHA256,
            "target_definition_from_frozen_authority": authority["target_definition"],
            "physical_scope": "The full original exact-FE physical domain and unchanged closed C, "
                              "fixed model, active mask and record doses apply.",
            "nodal_set": "For each theta in C, select any free error e with max norm <= supplied radius. "
                         "The box encloses reference plus injection of u + S(theta-base) + shift + e. "
                         "The radius is at least the frozen old corrector radius.",
            "center_policy": "Each exact free center is reference + proof displacement + explicit shift. "
                             "The coordinate radius uses both signed exact affine endpoints. "
                             "The default shift is the zero vector.",
            "box_requires_root_membership": False,
            "engine_boundary": "Binding and exact affine enclosure only. Independent root review required "
                               "before any cohort integration job.",
            "failed_attempts": {
                path.relative_to(OUT).as_posix(): sha256(path.read_bytes()).hexdigest()
                for path in sorted(OUT.glob("development_revision_*/*")) if path.is_file()
            },
        }
        with (OUT / names[0]).open("xb") as file:
            file.write(source)
        protocol_raw = (json.dumps(protocol, indent=2) + "\n").encode()
        with (OUT / names[1]).open("xb") as file:
            file.write(protocol_raw)
        report["protocol_sha256"] = sha256(protocol_raw).hexdigest()
        with (OUT / names[2]).open("x", encoding="utf8", newline="\n") as file:
            json.dump(report, file, indent=2)
            file.write("\n")
    summary = {key: report[key] for key in (
        "status", "source_sha256", "preflight_source_sha256", "elapsed_seconds", "check_count", "record_count")}
    if freeze:
        summary["protocol_sha256"] = report["protocol_sha256"]
        summary["preflight_sha256"] = sha256((OUT / names[2]).read_bytes()).hexdigest()
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", action="store_true")
    run(parser.parse_args().freeze)
