"""Independently replay the frozen cohort context bindings and exact box arithmetic."""
from base64 import b64decode, b64encode
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
from itertools import product
import json
from pathlib import Path, PurePosixPath
import sys
from unittest.mock import patch

import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
PREP = OUT.parent / "context_preparation"
EVIDENCE = ROOT / "evidence/exact_cohort_context_independent_review.json"
AUTHORITY = "runs/exact_integration_bridge_v1/target_preparation_safe_v2/protocol.json"
SOURCE = "scripts/exact_cohort_context.py"
FROZEN = {
    SOURCE: "c22612a1de1ccfd34e8b0f0859f6413c4da09c765ba4945d060340e93e928f9f",
    "protocol.json": "c1834232a5f92ef1e3d997aeba3bc6a543ef2fe4a2e1b0d20995b4a4fd9f2409",
    "preflight.json": "0795ee2ce87e8bf37e3bfcb33c3fb8611b361f1909be9506f722564b184b886a",
    AUTHORITY: "68164b3a922854b05d55c57d6872e92451e176484d7f5bf2ba3127101819850b",
}
ORDER = ("apical_area_modulus", "basal_area_modulus", "apical_surface_tension",
         "basal_surface_tension", "common_gain")
RECENTERED = {"optoGEF2_small_col_F", "optoGEF2_small_col_K"}
B_GEOMETRY = "runs/exact_integration_bridge_v1/retube_B_geometry"
B_CANDIDATE = "runs/exact_integration_bridge_v1/retube_B_preparation/candidate_vectors.npz"


def digest(data):
    return sha256(data).hexdigest()


def q(value):
    return Q.from_float(float(value)) if isinstance(value, (float, np.floating)) else Q(value)


def binary_equal(a, b):
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


def select(rows, record_id):
    row, = [row for row in rows if row["record_id"] == record_id]
    return row


def read_inputs():
    bindings = {}

    def raw(name, expected=None):
        path = PurePosixPath(name)
        assert path.as_posix() == name and not path.is_absolute() and ".." not in path.parts
        assert (ROOT / path).resolve().is_relative_to(ROOT)
        data = (ROOT / path).read_bytes()
        value = digest(data)
        if expected is not None:
            assert value == expected, name
        assert bindings.setdefault(name, value) == value, name
        return data

    def doc(name, expected=None):
        return json.loads(raw(name, expected))

    return bindings, raw, doc


def arrays(raw):
    with np.load(BytesIO(raw), allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in archive.files}


def corner_extrema(predictor, corners):
    extrema = []
    for row in predictor:
        coefficients = tuple(map(q, row))
        values = [sum(v * t for v, t in zip(coefficients, corner)) for corner in corners]
        extrema.append((min(values), max(values)))
    return extrema


def compare_box(c, box, a, extrema, radius, shift):
    free = {int(k): i for i, k in enumerate(a["free"])}
    assert box["triangles"] == a["triangles"].tolist()
    assert box["free"] == list(free) and box["root_membership_assumed"] is False
    assert box["record_id"] == c.record_id and box["corrector_radius_exact"] == radius
    assert box["shift_exact"] == tuple(shift)
    assert box["parameter_cell_sha256"] == c.row["parameter_cell_sha256"]
    assert box["model_observation_contract_sha256"] == c.row["model_observation_contract_sha256"]
    axis_zeros = 0
    for node, coordinate in product(range(325), range(2)):
        x = q(a["reference"][node, coordinate])
        center = box["current_centers_exact"][node][coordinate]
        extent = box["coordinate_radii_exact"][node][coordinate]
        actual = box["predictor_coordinate_hulls_exact"][node][coordinate]
        assert box["reference_nodes_exact"][node][coordinate] == x
        i = free.get(2 * node + coordinate)
        if i is None:
            assert center == x and actual == (x, x) and extent == 0
            axis_zeros += coordinate == 0 and x == 0
            continue
        assert coordinate != 0 or x != 0
        lo, hi = extrema[i]
        u = x + q(a["u"][i]) + shift[i]
        assert box["predictor_hulls_exact"][i] == (lo, hi)
        assert center == u and actual == (u + lo, u + hi)
        assert extent == max(abs(lo), abs(hi)) + radius
        assert center - extent <= u + lo - radius <= u + hi + radius <= center + extent
    return axis_zeros


def record(C, row, protocol, authoritative, raw, doc):
    record_id = row["record_id"]
    paths = ("root_certificate", "independent_root_capsule", "point_cache",
             "predictor", "proof_center_state")
    expected = {row[f"{p}_file"]: row[f"{p}_sha256"] for p in paths}
    assert all(authoritative["input_hashes"][name] == h for name, h in expected.items())
    root = doc(row["root_certificate_file"], row["root_certificate_sha256"])
    assert root["record_id"] == record_id and root["status"] == "CERTIFIED_UNIFORM_LOCAL_ROOT"
    capsule = doc(row["independent_root_capsule_file"], row["independent_root_capsule_sha256"])
    entry = select(capsule["record_audits"], record_id)
    audit = doc(entry["review_file"], entry["review_sha256"])
    assert audit == entry["audit"] and audit["certificate_sha256"] == row["root_certificate_sha256"]
    assert audit["status"] == "PASS" and audit["uniform_quantifier_verified"] is True
    assert audit["checks"] and all(v is True for v in audit["checks"].values())
    assert len(root["uniform_cell_result"]["full_nodal_predictor_bounds_exact"]) == 611
    for key in ("parameter_cell_sha256", "branch_domain_id", "model_observation_contract_sha256"):
        assert root[key] == row[key] == audit[key]
    for prefix in ("point_cache", "predictor"):
        path = row[f"{prefix}_file"]
        assert root[f"{prefix}_file"] == path and root["artifact_hashes"][path] == expected[path]
    original_name = (PurePosixPath(row["root_certificate_file"]).parents[2] / "protocol.json").as_posix()
    original = doc(original_name, root["protocol_sha256"])
    selected = select(original["records"], record_id)
    assert selected["output_folder"] + "/result.json" == row["root_certificate_file"]
    contract_name = original["common_model_observation_contract_file"]
    contract = doc(contract_name, row["model_observation_contract_sha256"])
    cell = {key: tuple(map(Q, row["parameter_cell"][key])) for key in ORDER}
    base = {key: Q(row["base_parameters_exact"][key]) for key in ORDER}
    cell_hash = digest(json.dumps({k: list(map(str, v)) for k, v in cell.items()},
                                  sort_keys=True, separators=(",", ":")).encode())
    assert cell_hash == protocol["parameter_cell_sha256"] == row["parameter_cell_sha256"]
    for other in (root["parameter_cell"], root["uniform_cell_result"]["cell_bounds_exact"],
                  protocol["parameter_cell"], contract["parameter_cell"], original["shared_parameter_cell"]):
        assert {k: tuple(map(Q, v)) for k, v in other.items()} == cell
    assert {k: Q(v) for k, v in original["base_parameters_exact"].items()} == base
    assert all(cell[k][0] <= b <= cell[k][1] for k, b in base.items())
    dose = Q(row["dose_exact"])
    drive = dose * base["common_gain"]
    assert dose == Q(selected["dose_exact"]) == Q(root["dose_exact"]) > 0
    assert drive == Q(root["base_strength_exact"]) == Q(selected["exact_requested_strength"]) > 0
    packet = doc(row["point_cache_file"], row["point_cache_sha256"])
    state_path = row["proof_center_state_file"]
    state_hash = row["proof_center_state_sha256"]
    state = doc(state_path, state_hash)
    assert packet.get("record_id", record_id) == record_id and packet["input_state"] == state_path
    assert packet["input_state_sha256"] == state_hash
    assert packet["parameters"] == state["parameters"] == original["fixed_full_parameters"]
    assert {k: v for k, v in packet["parameters"].items() if k not in ORDER[:4]} == contract["fixed_full_parameters"]
    assert all(q(packet["parameters"][k]) == base[k] for k in ORDER[:4])
    assert Q(packet.get("dose_exact", dose)) == dose
    assert Q(packet.get("exact_requested_strength", drive)) == drive
    cached_drive = Q(packet["activation_strength_binary_exact"])
    assert cached_drive == Q(selected["solver_strength_binary_exact"]) == q(float(drive))
    data = b64decode(packet["frozen_binary_cache"]["data"], validate=True)
    assert digest(data) == packet["frozen_binary_cache"]["sha256"] == row["frozen_binary_cache_sha256"]
    a = arrays(data)
    assert len(data) == packet["frozen_binary_cache"]["bytes"]
    assert {k: list(v.shape) for k, v in a.items()} == packet["frozen_binary_cache"]["array_shapes"]
    assert a["u"].shape == (611,) and a["preconditioner"].shape == (611, 611)
    assert a["u"].dtype == a["preconditioner"].dtype == np.dtype("float64")
    for name, key in (("u", "free_displacement"), ("reference", "reference_nodes"), ("triangles", "triangles")):
        assert binary_equal(a[name], np.asarray(state[key], dtype=a[name].dtype))
    s, = arrays(raw(row["predictor_file"], row["predictor_sha256"])).values()
    assert s.dtype == np.dtype("float64") and s.shape == (611, 5) and np.isfinite(s).all()
    historical = packet["matrix_archive"]
    real_read = Path.read_bytes

    def reject_historical(path):
        assert path != ROOT / historical, "historical matrix archive access"
        return real_read(path)

    with patch.object(Path, "read_bytes", reject_historical):
        c = C.load_context(record_id)
    assert c.record_id == record_id and c.row == row
    assert c.input_hashes == protocol["record_input_hashes"][record_id]
    for name, h in c.input_hashes.items():
        raw(name, h)
    assert binary_equal(c.S, s)
    assert set(c.target.a) == set(a) and all(binary_equal(a[k], c.target.a[k]) for k in a)
    assert all(not value.flags.writeable for value in [c.S, *c.target.a.values()])
    assert c.cell == cell and c.base == base and c.dose == dose and c.drive == drive
    assert c.cached_binary_drive == cached_drive
    old = Q(root["uniform_cell_result"]["selected_radius_exact"])
    assert old == c.old_radius == Q(row["selected_corrector_radius_exact"]) > 0

    recentering = None
    if record_id in RECENTERED:
        assert root["proof_center_state_file"] == state_path and root["proof_center_state_sha256"] == state_hash
        assert root["artifact_hashes"][state_path] == state_hash
        seed_path, seed_hash = root["original_seed_state_file"], root["original_seed_state_sha256"]
        seed = doc(seed_path, seed_hash)
        assert seed_path != state_path and seed_path == state["original_seed_state_file"] == selected["saved_state"]
        assert seed_hash == state["original_seed_state_sha256"] == selected["saved_state_sha256"]
        assert Q(state["proof_activation_strength_exact"]) == drive == Q(packet["exact_base_gain_times_dose"])
        refinement = doc(root["center_refinement_file"], root["center_refinement_sha256"])
        vectors = arrays(raw(refinement["center_arrays_file"], refinement["center_arrays_sha256"]))
        assert binary_equal(vectors["new_u"], a["u"])
        assert binary_equal(vectors["old_u"], np.asarray(seed["free_displacement"], dtype=np.float64))
        difference = [q(y) - q(x) for x, y in zip(vectors["old_u"], vectors["new_u"])]
        assert any(difference)
        assert difference == list(map(Q, refinement["actual_center_translation_exact"]))
        assert max(map(abs, difference)) == Q(refinement["translation_linf_exact"])
        original_packet = doc(packet["original_point_cache_file"], packet["original_point_cache_sha256"])
        original_embedded = original_packet["frozen_binary_cache"]
        original_bytes = b64decode(original_embedded["data"], validate=True)
        assert digest(original_bytes) == original_embedded["sha256"]
        original_a = arrays(original_bytes)
        assert set(original_a) == set(a)
        assert all(binary_equal(a[key], original_a[key]) for key in a if key != "u")
        assert refinement["S_sha256"] == row["predictor_sha256"]
        assert raw(refinement["S_file"], refinement["S_sha256"]) == raw(row["predictor_file"], row["predictor_sha256"])
        recentering = {"seed_file": seed_path, "proof_center_file": state_path,
                       "changed_coordinates": sum(v != 0 for v in difference),
                       "translation_linf_exact": str(max(map(abs, difference))),
                       "all_non_u_cache_arrays_bitwise_unchanged": True}
    else:
        assert root["saved_state"] == state_path == selected["saved_state"]
        assert root["saved_state_sha256"] == state_hash == selected["saved_state_sha256"]
        assert q(state["activation_strength"]) == cached_drive

    corners = tuple(product(*(tuple(x - base[k] for x in cell[k]) for k in ORDER)))
    extrema = corner_extrema(s, corners)
    zero = [Q(0)] * 611
    unshifted = C.exact_node_box(c, old)
    axis = compare_box(c, unshifted, a, extrema, old, zero)
    shift = [Q((-1)**j * (j % 7), 10**12) for j in range(611)]
    bigger = old + max(map(abs, shift))
    shifted = C.exact_node_box(c, bigger, shift)
    assert compare_box(c, shifted, a, extrema, bigger, shift) == axis == 13
    indexes = {int(k): i for i, k in enumerate(a["free"])}
    for node, j in product(range(325), range(2)):
        i = indexes.get(2 * node + j)
        mid = shifted["current_centers_exact"][node][j]
        ext = shifted["coordinate_radii_exact"][node][j]
        lo, hi = unshifted["predictor_coordinate_hulls_exact"][node][j]
        e = old if i is not None else Q(0)
        d = shift[i] if i is not None else Q(0)
        assert mid - ext <= lo + min(Q(0), d) - e
        assert hi + max(Q(0), d) + e <= mid + ext
    info = {"record_id": record_id, "binary_cache_sha256": digest(data),
            "u_bytes_sha256": digest(a["u"].tobytes()), "R_bytes_sha256": digest(a["preconditioner"].tobytes()),
            "S_bytes_sha256": digest(s.tobytes()), "record_binding_count": len(c.input_hashes),
            "exact_corner_coordinates": 650, "exact_shift_and_sweep_coordinates": 650,
            "axis_zero_coordinates": axis, "recentered": recentering}
    return c, info


def regression(C, c, raw, doc):
    review = doc("evidence/exact_retube_B_geometry_independent_review.json")
    assert review["status"].startswith("PASS") and review["checks"]
    assert all(value is True for value in review["checks"].values())
    geometry = doc(f"{B_GEOMETRY}/protocol.json", review["protocol_sha256"])
    raw(f"{B_GEOMETRY}/certification_run/results.json", review["producer_results_sha256"])
    box = doc(f"{B_GEOMETRY}/coordinate_boxes.json", geometry["prepared_box_sha256"])
    assert review["body_output_sha256"] == digest(raw(f"{B_GEOMETRY}/certification_run/whole_body_geometry.json"))
    assert box["parameter_cell_sha256"] == c.row["parameter_cell_sha256"]
    r = Q(geometry["new_common_enclosing_radius_exact"])
    new = C.exact_node_box(c, r)
    vectors = arrays(raw(B_CANDIDATE, geometry["source_and_input_hashes"][B_CANDIDATE]))
    assert binary_equal(vectors["u0"], c.u)
    shifts = [q(y) - q(x) for x, y in zip(vectors["u0"], vectors["u1"])]
    assert r == c.old_radius + max(map(abs, shifts))
    index = {int(k): i for i, k in enumerate(c.free)}
    equal = 0
    for node, j in product(range(325), range(2)):
        k = 2 * node + j
        i = index.get(k)
        old_mid = Q(box["current_centers_exact"][node][j])
        old_extent = Q(box["coordinate_radii_exact"][node][j])
        new_mid = new["current_centers_exact"][node][j]
        new_extent = new["coordinate_radii_exact"][node][j]
        outer = (old_mid - old_extent, old_mid + old_extent)
        inner = (new_mid - new_extent, new_mid + new_extent)
        assert outer == tuple(map(Q, box["rows"][k]["common_box_exact"]))
        assert outer[0] <= inner[0] <= inner[1] <= outer[1]
        assert old_extent == new_extent + abs(new_mid - old_mid)
        equal += outer == inner
        lo, hi = new["predictor_coordinate_hulls_exact"][node][j]
        d = shifts[i] if i is not None else Q(0)
        e = c.old_radius if i is not None else Q(0)
        swept = (lo + min(Q(0), d) - e, hi + max(Q(0), d) + e)
        assert swept == tuple(map(Q, box["rows"][k]["swept_translation_cube_hull_exact"]))
        assert inner[0] <= swept[0] <= swept[1] <= inner[1]
    return {"coordinates_exactly_equal": equal, "coordinates_strictly_contained": 650 - equal,
            "coordinates_checked_against_original_swept_endpoints": 650,
            "explanation": "Finished B centers on serialized current nodes and adds the exact binary reconstruction gap. The new box centers on exact X + u and therefore nests inside that completed geometry box at the same radius."}


def negative(C, c, authority, raw, doc):
    checks = {}

    def reject(label, call, expected=None):
        try:
            call()
        except ValueError as error:
            assert expected is None or expected in str(error), (label, str(error))
        else:
            raise AssertionError(f"accepted malformed fixture: {label}")
        checks[label] = True

    other = next(row for row in authority["records"] if row["record_id"] != c.record_id)
    foreign_root = doc(other["root_certificate_file"])
    root = doc(c.row["root_certificate_file"])
    reject("foreign_root", lambda: C._validate_root_links(c.row, foreign_root, authority["input_hashes"]), "root identity")
    for kind in ("point_cache", "predictor", "proof_center_state", "root_certificate"):
        real_read = Path.read_bytes
        watched = ROOT / c.row[f"{kind}_file"]
        replacement = raw(other[f"{kind}_file"])
        with patch.object(Path, "read_bytes", lambda path: replacement if path == watched else real_read(path)):
            reject(f"foreign_{kind}_bytes", lambda: C.load_context(c.record_id), "changed frozen artifact")
    bad = deepcopy(root)
    bad["artifact_hashes"][c.row["predictor_file"]] = other["predictor_sha256"]
    reject("foreign_predictor_root_hash", lambda: C._validate_root_links(c.row, bad, authority["input_hashes"]), "path/digest mismatch")
    for record_id in sorted(RECENTERED):
        row = select(authority["records"], record_id)
        fk = doc(row["root_certificate_file"])
        fk["proof_center_state_file"] = fk["original_seed_state_file"]
        reject(f"seed_substituted_for_{record_id}", lambda: C._validate_root_links(row, fk, authority["input_hashes"]), "proof-center path/digest")
    for label, val in (("zero", 0), ("below_old", c.old_radius / 2), ("bool", True),
                       ("nan", float("nan")), ("bad_fraction", "1/0")):
        reject(f"radius_{label}", lambda: C.exact_node_box(c, val))
    for label, val in (("rank_two", np.zeros((611, 1))), ("ragged", [0] * 610 + [[1]]),
                       ("bool", [0] * 610 + [False]), ("infinity", [0] * 610 + [float("inf")])):
        reject(f"shift_{label}", lambda: C.exact_node_box(c, c.old_radius, val))

    state = doc(c.row["proof_center_state_file"])
    packet = deepcopy(c.target.packet)
    old_bytes = b64decode(packet["frozen_binary_cache"]["data"], validate=True)
    packet["frozen_binary_cache"]["data"] = b64encode(old_bytes + b"!").decode()
    reject("embedded_cache_byte_tamper", lambda: C._target(packet, c.row, state), "cache bytes or digest")
    packet = deepcopy(c.target.packet)
    a = arrays(old_bytes)
    axis_dof = next(2 * j for j, ref in enumerate(a["reference"]) if ref[0] == 0)
    a["free"][0] = axis_dof
    buf = BytesIO()
    np.savez(buf, **a)
    altered = buf.getvalue()
    h = digest(altered)
    packet["frozen_binary_cache"].update(data=b64encode(altered).decode(), sha256=h, bytes=len(altered))
    row = dict(c.row, frozen_binary_cache_sha256=h)
    reject("otherwise_well_hashed_nonexact_axis_injection", lambda: C._target(packet, row, state), "axis coordinates must be fixed")
    return checks


def write_once(path, report):
    with path.open("x", encoding="utf8", newline="\n") as handle:
        handle.write(json.dumps(report, indent=2) + "\n")


def main():
    if not __debug__:
        raise RuntimeError("Replay requires enabled assertions")
    bindings, raw, doc = read_inputs()
    prep = PREP.relative_to(ROOT).as_posix()
    protocol = doc(f"{prep}/protocol.json", FROZEN["protocol.json"])
    preflight = doc(f"{prep}/preflight.json", FROZEN["preflight.json"])
    assert raw(SOURCE, FROZEN[SOURCE]) == raw(f"{prep}/{protocol['source_snapshot']}", FROZEN[SOURCE])
    assert protocol["source_sha256"] == preflight["source_sha256"] == FROZEN[SOURCE]
    assert preflight["protocol_sha256"] == FROZEN["protocol.json"]
    assert preflight["status"] == "PASS_CONTEXT_PREPARATION_ONLY"
    assert preflight["checks"] and all(v is True for v in preflight["checks"].values())
    assert all(preflight[k] is False for k in ("exact_integration_started", "matrix_assembly_started",
                                             "new_equilibrium_solves", "new_exact_FE_root_claims"))
    raw(protocol["preflight_source_file"], protocol["preflight_source_sha256"])
    assert protocol["preflight_source_sha256"] == preflight["preflight_source_sha256"]
    for path, hash_value in protocol["failed_attempts"].items():
        raw(f"{prep}/{path}", hash_value)
    authority = doc(AUTHORITY, FROZEN[AUTHORITY])
    assert protocol["frozen_authority"] == {"file": AUTHORITY, "sha256": FROZEN[AUTHORITY]}
    assert protocol["target_definition_from_frozen_authority"] == authority["target_definition"]
    sys.path.insert(0, str(ROOT / "scripts"))
    import exact_cohort_context as C
    checks = {"frozen_source_protocol_preflight_and_authority_exact": True}
    ids = [row["record_id"] for row in authority["records"]]
    assert len(ids) == len(set(ids)) == 17
    assert set(ids) == set(protocol["record_input_hashes"]) == {row["record_id"] for row in preflight["records"]}
    rows = []
    first = None
    for row in authority["records"]:
        c, info = record(C, row, protocol, authority, raw, doc)
        if c.record_id == "optoGEF2_small_col_B":
            first = c
        rows.append(info)
    assert {row["record_id"] for row in rows if row["recentered"] is not None} == RECENTERED
    checks.update({"all17_independent_authority_and_binary_cache_replay": True,
                   "all17_exact_32_corner_unshifted_and_shifted_box_replay": True,
                   "all17_shifted_boxes_contain_original_and_translated_old_cubes_and_swept_segments": True,
                   "F_and_K_distinct_corrected_centers_bound_to_original_model_R_and_S": True,
                   "all17_loads_avoid_historical_matrix_archives": True})
    b = regression(C, first, raw, doc)
    checks["finished_B_exact_common_box_containment_regression"] = True
    checks.update(negative(C, first, authority, raw, doc))
    replay_file = Path(__file__).relative_to(ROOT).as_posix()
    replay_sha = digest(raw(replay_file))
    for name, expected in tuple(bindings.items()):
        raw(name, expected)
    evidence = {"status": "PASS", "scope": "Independent selected-record context and exact coordinate hull preparation only. No integral, Hessian assembly, equilibrium solve or new root certification.",
                "source_sha256": FROZEN[SOURCE], "protocol_sha256": FROZEN["protocol.json"],
                "preflight_sha256": FROZEN["preflight.json"], "checks": checks,
                "records": 17, "review_source_file": replay_file, "review_source_sha256": replay_sha,
                "B_regression": b, "input_hashes": dict(sorted(bindings.items()))}
    receipt = OUT / "replay.json"
    if "--freeze" in sys.argv:
        assert not receipt.exists() and not EVIDENCE.exists()
        write_once(receipt, {**evidence, "records": rows})
        evidence["review_receipt_file"] = receipt.relative_to(ROOT).as_posix()
        evidence["review_receipt_sha256"] = digest(receipt.read_bytes())
        write_once(EVIDENCE, evidence)
    print(json.dumps({"status": evidence["status"], "checks": len(checks), "records": len(rows), "B_regression": b}))


if __name__ == "__main__":
    main()
