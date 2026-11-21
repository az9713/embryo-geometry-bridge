"""Independent source, context, nodal hull, incidence, and budget review for L."""

import argparse
import base64
import hashlib
import io
import itertools
import json
import sys
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RECORD = "optoGEF2_small_col_L"
SOURCE = ROOT / "scripts/cohort_exact_pilot.py"
INVENTORY = ROOT / "runs/exact_integration_cohort_v1/inventory/inventory.json"
AUTHORITY = ROOT / "runs/exact_integration_bridge_v1/target_preparation_safe_v2/protocol.json"
PILOT = ROOT / "runs/exact_integration_cohort_v1/records" / RECORD
ORDER = ("apical_area_modulus", "basal_area_modulus", "apical_surface_tension",
         "basal_surface_tension", "common_gain")
GEOMETRY = {
    "scripts/mixed_geometry_certificate_v2.py": "165d3b914b1ab9bad2e6591a7df63c09e004f33a742cfa5f4430584746edfbd2",
    "scripts/mixed_geometry_certificate.py": "fb0d49d680286ff49d0dc6713baca41fbebbf724442bff9731292e1441b03fd8",
}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def binary(value):
    return Q.from_float(float(value))


def ordered_lookup(rows):
    selected, = [row for row in rows if row["record_id"] == RECORD]
    return selected


def exact_row_norm(preconditioner):
    """Sum binary64 magnitudes as integers over one common power of two."""
    row_norms = []
    for row in preconditioner:
        parts = [abs(float(x)).as_integer_ratio() for x in row]
        power_of_two = max(denominator for _, denominator in parts)
        numerator = sum(value * (power_of_two // den) for value, den in parts)
        row_norms.append(Q(numerator, power_of_two))
    largest = max(row_norms)
    return largest, row_norms.index(largest)


def review(out, expected_source, expected_protocol, expected_preflight):
    protocol_file = out / "protocol.json"
    preflight_file = out / "preflight.json"
    protocol = read(protocol_file)
    preflight = read(preflight_file)
    hashes = protocol["source_and_input_hashes"]
    inventory = read(INVENTORY)
    authority = read(AUTHORITY)
    record = ordered_lookup(authority["records"])
    summary = ordered_lookup(inventory["records"])
    root = read(ROOT / record["root_certificate_file"])
    independent = ordered_lookup(read(ROOT / record["independent_root_capsule_file"])["record_audits"])
    audit = independent["audit"]
    row_review = ROOT / independent["review_file"]
    box_file = ROOT / protocol["coordinate_box_file"]
    box = read(box_file)
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)
        if not condition:
            raise ValueError(name)

    check("reviewed_source_protocol_and_preflight_bytes_exact",
          (sha(SOURCE.read_bytes()), sha(protocol_file.read_bytes()), sha(preflight_file.read_bytes()))
          == (expected_source, expected_protocol, expected_preflight))
    check("preflight_binds_source_protocol_and_preparation_only",
          preflight["status"] == "PASS_PREPARATION_ONLY" and bool(preflight["checks"]) and
          all(value is True for value in preflight["checks"].values()) and
          preflight["source_sha256"] == expected_source and preflight["protocol_sha256"] == expected_protocol and
          preflight["actual_body_or_force_job_run"] is False)
    check("all_declared_input_hashes_current", all(sha((ROOT / name).read_bytes()) == digest
                                                   for name, digest in hashes.items()))
    check("frozen_preparation_inventory_and_all_L_pointer_bindings",
          protocol["row"] == record and
          all(summary[key] == record[key] for key in summary.keys() & record.keys()) and
          inventory["all_input_sha256"][AUTHORITY.relative_to(ROOT).as_posix()] == sha(AUTHORITY.read_bytes()) and
          all(hashes[record[name + "_file"]] == record[name + "_sha256"] for name in
              ("root_certificate", "independent_root_capsule", "predictor", "proof_center_state", "point_cache")))
    check("selected_exact_parameter_cell_canonical_hash", sha(json.dumps(record["parameter_cell"], sort_keys=True,
                                                                          separators=(",", ":")).encode())
          == protocol["parameter_cell_sha256"] == root["parameter_cell_sha256"] == audit["parameter_cell_sha256"])
    check("selected_exact_original_branch_context",
          protocol["record_id"] == record["record_id"] == root["record_id"] == box["record_id"] == RECORD and
          root["branch_domain_id"] == record["branch_domain_id"] == audit["branch_domain_id"] and
          record["model_observation_contract_sha256"] == root["model_observation_contract_sha256"] ==
          audit["model_observation_contract_sha256"] == box["model_observation_contract_sha256"] and
          Q(root["dose_exact"]) == Q(record["dose_exact"]) == Q("1.51057344") and
          Q(root["base_strength_exact"]) == Q(record["dose_exact"]) * Q(record["base_parameters_exact"]["common_gain"]))
    check("cached_root_current_independently_audited", sha(row_review.read_bytes()) == independent["review_sha256"] and
          hashes[independent["review_file"]] == independent["review_sha256"] and
          audit["status"] == "PASS" and audit["uniform_quantifier_verified"] is True and
          bool(audit["checks"]) and all(v is True for v in audit["checks"].values()) and
          audit["certificate_sha256"] == record["root_certificate_sha256"])
    eta = max(Q(audit["bounds"]["eta"]["upper_exact"]),
              Q(root["uniform_cell_result"]["eta_uniform"]["upper_exact"]))
    z = Q(audit["bounds"]["Z"]["upper_exact"])
    r = Q(record["selected_corrector_radius_exact"])
    check("authoritative_eta_max_and_independent_Z_and_selected_normal_radius",
          eta == Q(protocol["cached_eta_upper_exact"]) and 0 <= z < 1 and
          z == Q(protocol["cached_Z_upper_exact"]) == Q(summary["cached_root_tube_Z_upper_exact"]) and
          r == Q(protocol["normal_radius_exact"]) == Q(box["corrector_radius_exact"]) ==
          Q(root["uniform_cell_result"]["selected_radius_exact"]) > 0)
    check("old_cached_certificate_strict_selfmap_continues_to_hold", eta + z * r < r)

    packet = read(ROOT / record["point_cache_file"])
    binary_packet = base64.b64decode(packet["frozen_binary_cache"]["data"], validate=True)
    check("binary_cache_digest_matches_frozen_authority", sha(binary_packet)
          == packet["frozen_binary_cache"]["sha256"] == record["frozen_binary_cache_sha256"])
    with np.load(io.BytesIO(binary_packet), allow_pickle=False) as archive:
        free, reference, u, triangle_nodes, preconditioner = (archive[name] for name in
                                                             ("free", "reference", "u", "triangles", "preconditioner"))
        edges = [archive[f"surface_{face}_edges"].tolist() for face in (0, 1)]
    with np.load(ROOT / record["predictor_file"], allow_pickle=False) as archive:
        predictor = archive["S"]
    check("binary64_predictor_preconditioner_and_free_map_shapes", predictor.shape == (611, 5) and
          preconditioner.shape == (611, 611) and len(set(free)) == len(u) == len(free) == 611 and
          np.isfinite(preconditioner).all() and np.isfinite(predictor).all() and
          predictor.dtype == preconditioner.dtype == np.dtype("float64"))
    check("fresh_box_bound_and_no_root_membership_assumed", sha(box_file.read_bytes())
          == protocol["coordinate_box_sha256"] == hashes[protocol["coordinate_box_file"]] and
          box["root_membership_assumed"] is False and protocol["root_asserted"] is False and
          len(box["shift_exact"]) == 611 and all(Q(v) == 0 for v in box["shift_exact"]))
    state = read(ROOT / record["proof_center_state_file"])
    check("proof_center_state_matches_binary_mesh_and_u", np.array_equal(state["free_displacement"], u) and
          np.array_equal(state["reference_nodes"], reference) and np.array_equal(state["triangles"], triangle_nodes))
    check("fresh_box_mesh_and_map_match_record", box["free"] == free.tolist() and
          box["triangles"] == triangle_nodes.tolist() and len(reference) == len(box["reference_nodes_exact"]) ==
          len(box["current_centers_exact"]) == len(box["coordinate_radii_exact"]) == 325)

    delta = [tuple(Q(endpoint) - Q(record["base_parameters_exact"][name])
                   for endpoint in record["parameter_cell"][name]) for name in ORDER]
    vertices = tuple(itertools.product(*delta))
    index = {int(g): i for i, g in enumerate(free)}
    check("all_32_vertices_of_unchanged_independent_parameter_cell", len(set(vertices)) == 32)
    mismatch = []
    for node, coordinate in itertools.product(range(325), range(2)):
        i = index.get(2 * node + coordinate)
        fixed_value = binary(reference[node, coordinate])
        center = Q(box["current_centers_exact"][node][coordinate])
        radius = Q(box["coordinate_radii_exact"][node][coordinate])
        hull = tuple(map(Q, box["predictor_coordinate_hulls_exact"][node][coordinate]))
        if i is None:
            good = (center, radius, hull) == (fixed_value, 0, (fixed_value, fixed_value))
        else:
            values = [sum((binary(x) * d for x, d in zip(predictor[i], v)), Q(0)) for v in vertices]
            lower, upper = min(values), max(values)
            good = (center == fixed_value + binary(u[i]) and
                    radius == max(abs(lower), abs(upper)) + r and
                    tuple(map(Q, box["predictor_hulls_exact"][i])) == (lower, upper) and
                    hull == (center + lower, center + upper) and
                    center - radius <= hull[0] - r <= hull[1] + r <= center + radius)
        if Q(box["reference_nodes_exact"][node][coordinate]) != fixed_value or not good:
            mismatch.append((node, coordinate))
    check("all_650_exact_centers_radii_and_32_corner_hulls_enclose_original_cube", not mismatch)
    axis = [i for i in range(325) if reference[i, 0] == 0]
    check("all_axis_radial_nodes_remain_exactly_fixed",
          all(2 * i not in index and Q(box["current_centers_exact"][i][0]) ==
              Q(box["coordinate_radii_exact"][i][0]) == 0 for i in axis))

    r_norm, row_index = exact_row_norm(preconditioner)
    check("exact_integer_row_reduction_matches_driver_R_norm_and_inventory",
          r_norm == Q(protocol["R_infinity_norm_exact"]) == Q(summary["cached_root_tube_R_infinity_norm_exact"]))
    selectors = [tuple(index[k] for k in sorted({2 * node + axis for node in element for axis in (0, 1)})
                       if k in index) for element in triangle_nodes.tolist() + edges[0] + edges[1]]
    check("144_bulk_and_24_surface_unique_local_free_maps", len(selectors) == 168 and
          {len(face) for face in edges} == {12} and all(len(s) == len(set(s)) > 0 for s in selectors))
    f_counts = Counter(itertools.chain.from_iterable(selectors))
    h_counts = Counter(i for selector in selectors for i, _ in itertools.product(selector, repeat=2))
    expected_f = [f_counts[i] for i in range(611)]
    expected_h = [h_counts[i] for i in range(611)]
    fc, hc = max(expected_f), max(expected_h)
    check("all_611_force_and_Hessian_row_incidence_counts_reconstructed",
          expected_f == protocol["force_incidence"] and expected_h == protocol["H_row_incidence"] and
          (fc, hc) == (protocol["max_force_incidence"], protocol["max_H_row_incidence"]))
    f_budget = Q(protocol["force_tail_budget_exact"])
    h_budget = Q(protocol["H_entry_tail_budget_exact"])
    check("exact_local_tail_accounting_and_infinity_norm_budget",
          protocol["margin_divisor"] == 10000 and f_budget == (1-z)*r / (10000*r_norm*fc) and
          h_budget == (1-z) / (10000*r_norm*hc) and (r_norm*fc*f_budget, r_norm*hc*h_budget) ==
          (Q(protocol["global_force_tail_budget_exact"]), Q(protocol["global_derivative_tail_budget_exact"])) ==
          ((1-z)*r/10000, (1-z)/10000))
    guard = protocol.get("surface_domain_depth") == 8 and protocol.get("max_surface_order") == 64
    check("bounded_selector_policy_640_precision_depth8_degree32_and64",
          guard and protocol["precision_bits"] == 640 and protocol["max_bulk_order"] == 32)
    return {
        "checks": checks,
        "record_id": RECORD,
        "source_sha256": expected_source,
        "protocol_sha256": expected_protocol,
        "preflight_sha256": expected_preflight,
        "geometry_dependency_bindings": {
            name: {"expected_sha256": digest, "checked_current_sha256": sha((ROOT / name).read_bytes()),
                   "driver_protocol_sha256": hashes.get(name)} for name, digest in GEOMETRY.items()
        },
        "exact_R_infinity_norm": str(r_norm),
        "exact_R_maximizing_row": row_index,
        "max_force_incidence": fc,
        "max_H_row_incidence": hc,
        "local_force_tail_exact": str(f_budget),
        "local_H_entry_tail_exact": str(h_budget),
        "global_preconditioned_force_tail_exact": protocol["global_force_tail_budget_exact"],
        "global_preconditioned_derivative_tail_exact": protocol["global_derivative_tail_budget_exact"],
        "remaining_old_Q_selfmap_slack_exact": str(r - eta - z*r),
        "source_and_input_binding_count": len(hashes),
        "scope": "Original L proof cube independent byte/metadata, nodal-box, row-norm and scalar-budget comparisons only.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run_directory")
    parser.add_argument("source_sha256")
    parser.add_argument("protocol_sha256")
    parser.add_argument("preflight_sha256")
    args = parser.parse_args()
    print(json.dumps(review(PILOT / args.run_directory, args.source_sha256,
                            args.protocol_sha256, args.preflight_sha256), indent=2))
