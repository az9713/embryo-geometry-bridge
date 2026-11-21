"""Verify existing cohort receipts and summarize their scopes without rerunning proofs."""

import base64
import hashlib
import io
import json
from fractions import Fraction
from pathlib import Path
import sys

# This script name is also the name of a standard-library module NumPy uses.
HERE = Path(__file__).resolve().parent
if sys.path and Path(sys.path[0]).resolve() == HERE:
    sys.path.pop(0)

import numpy as np

ROOT = HERE.parents[2]
PROTOCOL = "runs/exact_integration_bridge_v1/target_preparation_safe_v2/protocol.json"
JOINED = "evidence/shared_cell_cohort_inference.json"
ENERGY_REVIEW = "evidence/cohort_energy_stability_independent_review.json"
EXACT_B_REVIEW = "evidence/exact_integral_root_B_safe_independent_review.json"
EXACT_B_DELIVERY = "evidence/exact_integral_root_B_final_delivery.json"
EXACT_B = "runs/exact_integration_bridge_v1/integral_root_B_safe_v2/results.json"
EXACT_B_GEOMETRY = "runs/exact_integration_bridge_v1/retube_B_geometry_independent_review/results.json"


def require(condition, description):
    if not condition:
        raise ValueError(description)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def map_rows(document, key):
    rows = document[key]
    mapping = {row["record_id"]: row for row in rows}
    require(len(mapping) == len(rows), f"Duplicate record IDs in {key}")
    return mapping


def checked_json(filename, hashes, expected=None):
    raw = (ROOT / filename).read_bytes()
    actual = sha(raw)
    require(expected is None or actual == expected, f"Changed input {filename}")
    require(hashes.setdefault(filename, actual) == actual, f"Inconsistent input {filename}")
    return json.loads(raw)


def check_bound_bytes(filename, expected, hashes):
    actual = sha((ROOT / filename).read_bytes())
    require(actual == expected, f"Changed input {filename}")
    require(hashes.setdefault(filename, actual) == actual, f"Inconsistent input {filename}")


def matching_checks(record, status):
    return record["status"] == status and all(v is True for v in record["checks"].values())


def read_geometric_cache(packet, embedded_digest, predictor_file):
    embedded = packet["frozen_binary_cache"]
    raw = base64.b64decode(embedded["data"], validate=True)
    require(sha(raw) == embedded["sha256"] == embedded_digest, "Embedded cache hash")
    with np.load(io.BytesIO(raw), allow_pickle=False) as archive:
        # R is already bound by the compressed archive digest. Do not load or process it.
        arrays = {name: archive[name] for name in archive.files if name != "preconditioner"}
    with np.load(ROOT / predictor_file, allow_pickle=False) as archive:
        predictor = archive["S"]
    require(predictor.shape == (611, 5) and np.isfinite(predictor).all(), "Invalid predictor")
    return arrays


def common_geometry(arrays):
    reference = arrays["reference"]
    triangles = arrays["triangles"]
    free = arrays["free"]
    require(reference.shape == (325, 2) and triangles.shape == (144, 6), "Mesh shape")
    require(free.shape == (611,) and len(set(free.tolist())) == 611, "Free mapping")
    axis_vertex_count = (reference[triangles[:, :3], 0] == 0).sum(axis=1)
    mask = arrays["active"].reshape(-1)
    require(set(mask.tolist()) == {0, 1} and sum(mask) == 24, "Activation mask")
    fixed = sorted(set(range(650)) - set(free.tolist()))
    axis_radial = np.flatnonzero(reference[:, 0] == 0) * 2
    require(set(axis_radial.tolist()).issubset(fixed), "Unconstrained axis radius")
    classes = {}
    for count, name in enumerate(("off_axis", "axis_vertex", "axis_edge")):
        indices = np.flatnonzero(axis_vertex_count == count)
        classes[name] = {
            "elements": indices.tolist(),
            "active_elements": indices[mask[indices] == 1].tolist(),
        }
    require([len(row["elements"]) for row in classes.values()] == [132, 6, 6], "Axis classes")
    return {
        "nodes": 325,
        "triangles": 144,
        "geometry_classes_by_zero_radius_vertex_count": classes,
        "active_element_indices": np.flatnonzero(mask).tolist(),
        "free_to_full_interleaved_coordinate_indices": free.tolist(),
        "fixed_full_interleaved_coordinate_indices": fixed,
        "fixed_axis_radial_full_coordinate_indices": axis_radial.tolist(),
        "cached_surface_edge_node_indices": [arrays[f"surface_{k}_edges"].tolist() for k in (0, 1)],
        "cached_volume_sites_per_element": arrays["shape"].shape[0],
        "cached_surface_sites_per_edge": arrays["surface_0_N"].shape[0],
        "energy_scale_binary_exact": str(Fraction(float(arrays["energy_scale"][0]))),
    }


def build_inventory():
    hashes = {}
    protocol = checked_json(PROTOCOL, hashes)
    for filename, expected in protocol["input_hashes"].items():
        check_bound_bytes(filename, expected, hashes)
    joined = map_rows(checked_json(JOINED, hashes), "proof_provenance")
    energy_review = checked_json(ENERGY_REVIEW, hashes)
    require(matching_checks(energy_review, "PASS"), "Cohort energy review")
    energy_reviews = map_rows(energy_review, "record_audits")
    records = []
    first_geometry_arrays = None
    common_cell = None

    for metadata in sorted(protocol["records"], key=lambda row: row["record_id"]):
        record_id = metadata["record_id"]
        root = checked_json(metadata["root_certificate_file"], hashes, metadata["root_certificate_sha256"])
        capsule = checked_json(metadata["independent_root_capsule_file"], hashes,
                               metadata["independent_root_capsule_sha256"])
        root_row = map_rows(capsule, "record_audits")[record_id]
        audit = root_row["audit"]
        check_bound_bytes(root_row["review_file"], root_row["review_sha256"], hashes)
        require(matching_checks(audit, "PASS") and audit["uniform_quantifier_verified"] is True,
                f"Root review for {record_id}")
        require(root["record_id"] == record_id and root["status"] == "CERTIFIED_UNIFORM_LOCAL_ROOT",
                f"Root status for {record_id}")
        require(audit["certificate_sha256"] == metadata["root_certificate_sha256"],
                f"Root certificate digest for {record_id}")
        for name in ("branch_domain_id", "parameter_cell_sha256", "model_observation_contract_sha256"):
            require(root[name] == audit[name] == metadata[name], f"{name} for {record_id}")
        require(Fraction(root["dose_exact"]) == Fraction(metadata["dose_exact"]), f"Dose for {record_id}")
        require(root["parameter_cell"] == metadata["parameter_cell"], f"Parameter cell for {record_id}")
        cell = {name: metadata[name] for name in ("parameter_cell", "parameter_cell_sha256",
                                                 "model_observation_contract_sha256")}
        if common_cell is None:
            common_cell = cell
        require(cell == common_cell, f"Shared cell for {record_id}")

        point_file = metadata["point_cache_file"]
        predictor_file = metadata["predictor_file"]
        state_file = metadata["proof_center_state_file"]
        for name, filename in (("point_cache", point_file), ("predictor", predictor_file)):
            require(root[name + "_file"] == filename and
                    root["artifact_hashes"][filename] == metadata[name + "_sha256"],
                    f"Authoritative {name} for {record_id}")
            check_bound_bytes(filename, metadata[name + "_sha256"], hashes)
        packet = checked_json(point_file, hashes, metadata["point_cache_sha256"])
        require(packet["input_state"] == state_file and
                packet["input_state_sha256"] == metadata["proof_center_state_sha256"],
                f"Proof center binding for {record_id}")
        state = checked_json(state_file, hashes, metadata["proof_center_state_sha256"])
        arrays = read_geometric_cache(packet, metadata["frozen_binary_cache_sha256"], predictor_file)
        require(np.array_equal(arrays["u"], state["free_displacement"]), f"Proof center u for {record_id}")
        require(np.array_equal(arrays["reference"], state["reference_nodes"]) and
                np.array_equal(arrays["triangles"], state["triangles"]), f"Proof center mesh for {record_id}")
        geometry_arrays = {name: value for name, value in arrays.items() if name != "u"}
        if first_geometry_arrays is None:
            first_geometry_arrays = geometry_arrays
        require(geometry_arrays.keys() == first_geometry_arrays.keys() and
                all(np.array_equal(value, first_geometry_arrays[name])
                    for name, value in geometry_arrays.items()), f"Shared cache geometry for {record_id}")
        archive_file = packet["matrix_archive"]
        check_bound_bytes(archive_file, packet["matrix_archive_sha256"], hashes)

        uniform = root["uniform_cell_result"]
        radius = Fraction(metadata["selected_corrector_radius_exact"])
        require(Fraction(uniform["selected_radius_exact"]) == radius > 0, f"Radius for {record_id}")
        root_bounds = audit.get("bounds")
        reused = None
        if root_bounds is None:
            reused = audit["reuse_review_file"]
            reused_audit = checked_json(reused, hashes, audit["reuse_review_sha256"])
            root_bounds = reused_audit["independent_bounds"]

        observer_ref = joined[record_id]
        observer_file = observer_ref["observer_certificate_file"]
        observer = checked_json(observer_file, hashes, observer_ref["observer_certificate_sha256"])
        application = checked_json(observer_ref["application_review_file"], hashes,
                                   observer_ref["application_review_sha256"])
        application_row = map_rows(application, "record_audits")[record_id]
        require(matching_checks(application_row, "PASS_CONDITIONAL_UNIFORM_OBSERVER_APPLICATION") and
                application_row["root_certificate_sha256"] == metadata["root_certificate_sha256"] and
                application_row["observer_certificate_sha256"] == observer_ref["observer_certificate_sha256"],
                f"Observer/source binding for {record_id}")
        body = observer["body_geometry"]
        require(body["status"] == "CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING" and
                len(body["triangles"]) == 144 and body["global_injectivity_for_entire_box"] is True,
                f"Old whole-body geometry for {record_id}")

        energy_row = energy_reviews[record_id]
        energy_audit = energy_row["audit"]
        energy_file = f"runs/cohort_energy_stability_v1/records/{record_id}/result.json"
        energy = checked_json(energy_file, hashes, energy_audit["producer_result_sha256"])
        check_bound_bytes(energy_row["review_file"], energy_row["review_sha256"], hashes)
        energy_binding = energy["root_binding"]
        require(matching_checks(energy_audit, "PASS") and
                energy_audit["root_certificate_sha256"] == energy_binding["root_certificate_sha256"] ==
                metadata["root_certificate_sha256"] and energy["whole_tube_Hessian_positive"] is True and
                energy["strict_local_minimum_at_selected_root"] is True, f"Energy review for {record_id}")
        require(Fraction(energy_audit["normal_root_radius_exact"]) == radius, f"Energy tube for {record_id}")

        center_lineage = None
        if record_id.endswith(("_F", "_K")):
            refinement = checked_json(root["center_refinement_file"], hashes, root["center_refinement_sha256"])
            require(refinement["new_proof_center_state_file"] == state_file and
                    refinement["new_proof_center_state_sha256"] == metadata["proof_center_state_sha256"] and
                    refinement["S_sha256"] == metadata["predictor_sha256"] and
                    refinement["R_binary_unchanged"] is True, f"Recentered lineage for {record_id}")
            seed_file = refinement["original_seed_state_file"]
            seed = checked_json(seed_file, hashes, refinement["original_seed_state_sha256"])
            require(not np.array_equal(arrays["u"], seed["free_displacement"]), f"Original seed for {record_id}")
            center_lineage = {
                "center_refinement_file": root["center_refinement_file"],
                "original_seed_state_file": seed_file,
                "translation_linf_exact": refinement["translation_linf_exact"],
                "matrix_archive_is_old_seed_tangents": True,
                "R_and_S_unchanged_by_cached_root_recentering": True,
            }

        jacobian = body["J_bound"]
        rho = max(1 - Fraction(jacobian["lower_exact"]), Fraction(jacobian["upper_exact"]) - 1)
        entry = {name: value for name, value in metadata.items()
                 if name not in (*common_cell, "new_exact_target_root_certified")}
        entry.update({
            "reused_M_root_review_file": reused,
            "matrix_archive_file": archive_file,
            "proof_center_context_validated": True,
            "cached_root_recentering_lineage": center_lineage,
            "independent_root_record_review_file": root_row["review_file"],
            "cached_root_independent_bounds": root_bounds,
            "cached_root_predictor_error": uniform["uniform_error_from_predictor"],
            "cached_root_tube_Z_upper_exact": energy_audit["independent_fulltube_Z_upper_exact"],
            "cached_root_tube_quadrature_J_lower_exact": energy_binding["quadrature_J_lower_exact"],
            "cached_root_tube_hessian_lower_exact": energy_audit["uniform_Hessian_lower_bound_exact"],
            "cached_root_tube_R_infinity_norm_exact": energy_audit["exact_R_infinity_norm"],
            "cached_energy_file": energy_file,
            "cached_energy_independent_record_review_file": energy_row["review_file"],
            "old_nodal_box_observer_file": observer_file,
            "old_nodal_box_application_review_file": observer_ref["application_review_file"],
            "old_nodal_box_body_J_bound_exact": jacobian,
            "old_nodal_box_body_rho_upper_exact": str(rho),
            "old_nodal_box_body_strong_monotonicity_lower_exact": body["strong_monotonicity_lower_exact"],
            "old_nodal_box_surface_error_transfer": observer["surface_error_transfer"],
            "new_exact_integral_root_certified": record_id == "optoGEF2_small_col_B",
        })
        records.append(entry)

    require(len(records) == len(joined) == len(energy_reviews) == 17, "Cohort count")
    delivered = checked_json(EXACT_B_DELIVERY, hashes)
    require(matching_checks(delivered, "PASS_SAFE_UNIFORM_INTEGRAL_ROOT_COMPONENT_AND_PUBLICATION"),
            "Current B delivery")
    exact = checked_json(EXACT_B, hashes, delivered["artifact_hashes"][f"controlled_furrowing/{EXACT_B}"])
    exact_review = checked_json(EXACT_B_REVIEW, hashes,
                              delivered["artifact_hashes"][f"controlled_furrowing/{EXACT_B_REVIEW}"])
    require(exact_review["status"] == "PASS" and exact_review["certificate_sha256"] == hashes[EXACT_B] and
            exact["root_certified"] is True and exact["strict_local_minimum_certified"] is True,
            "Reviewed exact-integral B result")
    geometry = checked_json(EXACT_B_GEOMETRY, hashes)
    require(geometry["status"] == "PASS_FRESH_COMMON_NODAL_GEOMETRY_COMPONENT", "Fresh B geometry")
    worst_j = max(records, key=lambda row: Fraction(row["old_nodal_box_body_rho_upper_exact"]))
    worst_z = max(records, key=lambda row: Fraction(row["cached_root_tube_Z_upper_exact"]))
    return {
        "status": "PASS_READ_ONLY_RECEIPT_INVENTORY",
        "scope": "Recorded certificates and current file bindings. No proof, solve, integration or matrix calculation rerun.",
        "all_input_sha256": hashes,
        "common_cell_and_contract": common_cell,
        "common_cache_geometry": common_geometry(first_geometry_arrays),
        "records": records,
        "current_exact_B": {"certificate": EXACT_B, "review": EXACT_B_REVIEW, "delivery": EXACT_B_DELIVERY,
                            "fresh_geometry_review": EXACT_B_GEOMETRY},
        "ranking_from_existing_proofs": {"largest_old_body_rho": worst_j["record_id"],
                                         "largest_old_cached_tube_Z": worst_z["record_id"]},
        "recommended_generalization_pilot": "optoGEF2_small_col_L",
    }


if __name__ == "__main__":
    result = build_inventory()
    print(json.dumps(result, indent=2, sort_keys=True))
