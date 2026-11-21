"""Load one frozen selected record and construct a root-independent exact nodal box."""
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from io import BytesIO
from pathlib import Path, PurePosixPath
import base64
import json

import numpy as np

import exact_fe_target_safe as T


ROOT = Path(__file__).resolve().parents[1]
AUTHORITY = "runs/exact_integration_bridge_v1/target_preparation_safe_v2/protocol.json"
AUTHORITY_SHA256 = "68164b3a922854b05d55c57d6872e92451e176484d7f5bf2ba3127101819850b"
ORDER = T.ORDER
CELL_SHA256 = "64f93fefa034f4c004b93451d98b2fd04e0072855be8c5e3bba9babb8f394dbf"


@dataclass(frozen=True)
class CohortContext:
    record_id: str
    row: dict
    target: T.Target
    S: np.ndarray
    input_hashes: dict
    cell: dict
    base: dict
    dose: Fraction
    drive: Fraction
    old_radius: Fraction
    cached_binary_drive: Fraction

    @property
    def u(self):
        return self.target.a["u"]

    @property
    def R(self):
        return self.target.a["preconditioner"]

    @property
    def free(self):
        return self.target.a["free"]


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _verified_bytes(relative, expected_sha256, bindings):
    relative_path = PurePosixPath(relative)
    _require(relative_path.as_posix() == relative and not relative_path.is_absolute()
             and ".." not in relative_path.parts and "\\" not in relative and ":" not in relative,
             f"noncanonical artifact path: {relative}")
    path = ROOT / relative_path
    _require(path.resolve().is_relative_to(ROOT), f"artifact escapes study root: {relative}")
    _require(bindings.get(relative, expected_sha256) == expected_sha256,
             f"conflicting artifact digest: {relative}")
    raw = path.read_bytes()
    _require(sha256(raw).hexdigest() == expected_sha256, f"changed frozen artifact: {relative}")
    bindings[relative] = expected_sha256
    return raw


def _one_record(rows, record_id):
    matches = [row for row in rows if row["record_id"] == record_id]
    _require(len(matches) == 1, f"exactly one authoritative record required: {record_id!r}")
    return matches[0]


def _cell_digest(cell):
    canonical = {k: [str(value) for value in cell[k]] for k in ORDER}
    raw = json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode()
    return sha256(raw).hexdigest()


def _rational_cell(value):
    _require(set(value) == set(ORDER), "exact five-parameter cell required")
    cell = {key: tuple(map(T.rational, value[key])) for key in ORDER}
    _require(all(len(bounds) == 2 and bounds[0] <= bounds[1] for bounds in cell.values()),
             "closed ordered parameter intervals required")
    return cell


def _validate_root_links(row, root, authority_hashes):
    _require(root["status"] == "CERTIFIED_UNIFORM_LOCAL_ROOT"
             and root["record_id"] == row["record_id"] and root["free_dofs"] == 611,
             "selected root identity or status differs from authority")
    for field in ("parameter_cell_sha256", "model_observation_contract_sha256", "branch_domain_id"):
        _require(root[field] == row[field], f"selected root context mismatch: {field}")
    for prefix in ("point_cache", "predictor"):
        key = f"{prefix}_file"
        path = row[key]
        digest = row[f"{prefix}_sha256"]
        _require(root[key] == path and root["artifact_hashes"].get(path) == digest
                 and authority_hashes.get(path) == digest,
                 f"selected root artifact path/digest mismatch: {key}")
    _require(all(authority_hashes.get(path) == digest
                 for path, digest in root["artifact_hashes"].items()),
             "selected root artifact bindings differ from authority")

    recentered = "proof_center_state_file" in root
    center_key = "proof_center_state_file" if recentered else "saved_state"
    digest_key = "proof_center_state_sha256" if recentered else "saved_state_sha256"
    center_file = row["proof_center_state_file"]
    center_hash = row["proof_center_state_sha256"]
    _require(root[center_key] == center_file and root[digest_key] == center_hash
             and authority_hashes.get(center_file) == center_hash,
             "selected root proof-center path/digest mismatch")
    if recentered:
        _require(root["artifact_hashes"].get(center_file) == center_hash
                 and root["original_seed_state_file"] != center_file,
                 "recentered root requires its distinct proof center")
    return recentered


def _target(packet, row, state):
    embedded = packet["frozen_binary_cache"]
    raw = base64.b64decode(embedded["data"], validate=True)
    _require(len(raw) == embedded["bytes"] and embedded["sha256"] == row["frozen_binary_cache_sha256"]
             and sha256(raw).hexdigest() == embedded["sha256"], "embedded cache bytes or digest changed")
    target = T.Target(packet)
    arrays = target.a
    _require(embedded["array_shapes"] == {key: list(v.shape) for key, v in arrays.items()},
             "cache array manifest mismatch")
    _require(arrays["u"].shape == (611,) and arrays["preconditioner"].shape == (611, 611)
             and arrays["u"].dtype == arrays["preconditioner"].dtype == np.dtype("float64"),
             "binary64 proof displacement and preconditioner dimensions required")
    for state_key, cache_key in (("free_displacement", "u"), ("reference_nodes", "reference"),
                                 ("triangles", "triangles")):
        _require(np.array_equal(state[state_key], arrays[cache_key]),
                 f"embedded cache does not match proof center: {state_key}")
    for array in arrays.values():
        array.setflags(write=False)
    return target


def load_context(record_id):
    """Load only the exact-ID record, binding its current proof center, binary u/R and predictor S."""
    _require(type(record_id) is str and bool(record_id), "exact record ID string required")
    bindings = {}
    authority = json.loads(_verified_bytes(AUTHORITY, AUTHORITY_SHA256, bindings))
    row = _one_record(authority["records"], record_id)
    hashes = authority["input_hashes"]

    def selected_bytes(prefix):
        path = row[f"{prefix}_file"]
        digest = row[f"{prefix}_sha256"]
        _require(hashes.get(path) == digest, f"selected artifact is not bound by authority: {prefix}")
        return _verified_bytes(path, digest, bindings)

    for source in (Path(T.__file__), Path(T.core.__file__)):
        relative = source.resolve().relative_to(ROOT).as_posix()
        _verified_bytes(relative, authority["source_hashes"][relative], bindings)

    root = json.loads(selected_bytes("root_certificate"))
    recentered = _validate_root_links(row, root, hashes)
    review = json.loads(selected_bytes("independent_root_capsule"))
    audit = _one_record(review["record_audits"], record_id)["audit"]
    _require(audit["status"] == "PASS" and bool(audit["checks"])
             and all(value is True for value in audit["checks"].values())
             and audit["certificate_sha256"] == row["root_certificate_sha256"]
             and audit["uniform_quantifier_verified"] is True
             and all(audit[key] == row[key] for key in (
                 "parameter_cell_sha256", "model_observation_contract_sha256", "branch_domain_id")),
             "matching independent root capsule required")

    # Root hash binds the original run protocol for the base and source-dose checks.
    original_path = (PurePosixPath(row["root_certificate_file"]).parents[2] / "protocol.json").as_posix()
    original = json.loads(_verified_bytes(original_path, root["protocol_sha256"], bindings))
    source_row = _one_record(original["records"], record_id)
    _require(row["root_certificate_file"] == source_row["output_folder"] + "/result.json"
             and row["branch_domain_id"] == source_row["branch_domain_id"],
             "root protocol record and output path mismatch")
    contract = json.loads(_verified_bytes(original["common_model_observation_contract_file"],
                                          row["model_observation_contract_sha256"], bindings))
    _require(original["common_model_observation_contract_sha256"]
             == row["model_observation_contract_sha256"], "physical contract hash mismatch")

    uniform = root["uniform_cell_result"]
    cell = _rational_cell(row["parameter_cell"])
    _require(_cell_digest(cell) == CELL_SHA256 == row["parameter_cell_sha256"]
             == original["shared_parameter_cell_sha256"] == contract["parameter_cell_sha256"]
             and all(_rational_cell(candidate) == cell for candidate in (
                 root["parameter_cell"], uniform["cell_bounds_exact"], original["shared_parameter_cell"],
                 contract["parameter_cell"])), "selected record must retain the full shared closed C")
    base = {key: T.rational(row["base_parameters_exact"][key]) for key in ORDER}
    _require(base == {key: T.rational(original["base_parameters_exact"][key]) for key in ORDER}
             and all(cell[key][0] <= base[key] <= cell[key][1] for key in ORDER)
             and base["common_gain"] == T.rational(source_row["base_gain_exact"]),
             "selected exact base differs from the frozen base or leaves C")
    dose = T.rational(row["dose_exact"])
    drive = base["common_gain"] * dose
    old_radius = T.rational(row["selected_corrector_radius_exact"])
    _require(dose > 0 and all(T.rational(value) == dose for value in (
                 root["dose_exact"], uniform["dose_exact"], source_row["dose_exact"]))
             and drive > 0 and drive == T.rational(root["base_strength_exact"])
             == T.rational(source_row["exact_requested_strength"])
             and root["strength_relation"] == "common_gain_times_source_dose",
             "exact dose or gain-times-dose mismatch")
    _require(old_radius > 0 and old_radius == T.rational(uniform["selected_radius_exact"])
             and uniform["status"] == root["status"] and uniform["branch_domain_id"] == row["branch_domain_id"]
             and uniform["calibration_record_ids_covered"] == [record_id],
             "old branch radius or record mismatch")

    packet = json.loads(selected_bytes("point_cache"))
    cached_drive = T.rational(packet["activation_strength_binary_exact"])
    _require(packet["input_state"] == row["proof_center_state_file"]
             and packet["input_state_sha256"] == row["proof_center_state_sha256"]
             and packet.get("record_id", record_id) == record_id
             and T.rational(packet.get("dose_exact", dose)) == dose
             and T.rational(packet.get("exact_requested_strength", drive)) == drive
             and cached_drive == T.rational(source_row["solver_strength_binary_exact"]),
             "point cache proof-center, record or drive mismatch")
    state = json.loads(selected_bytes("proof_center_state"))
    variable_material_keys = set(ORDER[:-1])
    fixed = {key: value for key, value in packet["parameters"].items()
             if key not in variable_material_keys}
    _require(state["parameters"] == packet["parameters"] == original["fixed_full_parameters"]
             and fixed == contract["fixed_full_parameters"]
             and all(T.rational(packet["parameters"][key]) == base[key] for key in variable_material_keys),
             "cache changes the fixed physical model or base parameters")
    if recentered:
        _require(T.rational(state["proof_activation_strength_exact"]) == drive
                 == T.rational(packet["exact_base_gain_times_dose"])
                 and T.rational(state["original_solver_strength_binary_exact"]) == cached_drive
                 and state["original_seed_state_file"] == root["original_seed_state_file"] == source_row["saved_state"]
                 and state["original_seed_state_sha256"] == root["original_seed_state_sha256"]
                 == source_row["saved_state_sha256"], "recentered proof center mismatch")
    else:
        _require(T.rational(state["activation_strength"]) == cached_drive
                 and row["proof_center_state_file"] == source_row["saved_state"]
                 and row["proof_center_state_sha256"] == source_row["saved_state_sha256"],
                 "saved proof-center state mismatch")
    target = _target(packet, row, state)
    with np.load(BytesIO(selected_bytes("predictor")), allow_pickle=False) as archive:
        _require(archive.files == ["S"], "selected predictor archive must contain only S")
        predictor = archive["S"]
    _require(predictor.shape == (611, 5) and predictor.dtype == np.dtype("float64")
             and np.isfinite(predictor).all(), "finite binary64 611 by 5 predictor required")
    predictor.setflags(write=False)
    return CohortContext(record_id, row, target, predictor, bindings, cell, base,
                         dose, drive, old_radius, cached_drive)


def _exact_rational(value):
    # Floats denote their exact binary rational. Decimal strings denote the stated rational.
    try:
        return T.rational(value)
    except (TypeError, ValueError, OverflowError, ZeroDivisionError) as error:
        raise ValueError("finite non-Boolean rational required") from error


def exact_node_box(context, radius, shift=None):
    """Enclose X + injection(u + S·(theta-base) + shift + e), theta in C, ||e||∞ <= radius.

    The shift is a length-611 vector ordered like u, added to the old predictor. Each
    free coordinate is centered at X + u + shift. The symmetric radius encloses its
    exact affine interval plus the corrector radius. Constrained coordinates stay fixed.
    """
    radius = _exact_rational(radius)
    _require(radius >= context.old_radius, "radius must be at least the selected old radius")
    if shift is None:
        shift = (Fraction(0),) * 611
    else:
        _require(isinstance(shift, (tuple, list, np.ndarray)) and np.shape(shift) == (611,),
                 "shift must be an exact vector with shape (611,)")
        shift = tuple(map(_exact_rational, shift))

    reference = [list(map(T.rational, row)) for row in context.target.a["reference"]]
    centers = [row.copy() for row in reference]
    radii = [[Fraction(0), Fraction(0)] for _ in reference]
    coordinate_hulls = [[(v, v) for v in row] for row in reference]
    deltas = [tuple(endpoint - context.base[key] for endpoint in context.cell[key]) for key in ORDER]
    predictor_hulls = []
    for i, (coordinate, row) in enumerate(zip(context.free, context.S)):
        terms = [(value * lo, value * hi)
                 for value, (lo, hi) in zip(map(T.rational, row), deltas)]
        lower = sum(map(min, terms), Fraction(0))
        upper = sum(map(max, terms), Fraction(0))
        predictor_hulls.append((lower, upper))
        node, axis = divmod(int(coordinate), 2)
        center = reference[node][axis] + T.rational(context.u[i]) + shift[i]
        centers[node][axis] = center
        radii[node][axis] = max(-lower, upper) + radius
        coordinate_hulls[node][axis] = (center + lower, center + upper)
    return {
        "record_id": context.record_id,
        "parameter_cell_sha256": context.row["parameter_cell_sha256"],
        "model_observation_contract_sha256": context.row["model_observation_contract_sha256"],
        "root_membership_assumed": False,
        "reference_nodes_exact": reference,
        "current_centers_exact": centers,
        "coordinate_radii_exact": radii,
        "predictor_hulls_exact": predictor_hulls,
        "predictor_coordinate_hulls_exact": coordinate_hulls,
        "triangles": context.target.a["triangles"].tolist(),
        "free": context.free.tolist(),
        "corrector_radius_exact": radius,
        "shift_exact": shift,
    }
