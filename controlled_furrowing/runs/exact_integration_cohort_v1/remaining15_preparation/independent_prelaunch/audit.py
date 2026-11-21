"""Read-only math and provenance audit; commit domain-only receipts after all records pass."""

from base64 import b64decode
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np

OUT = Path(__file__).resolve().parent
BATCH = OUT.parent
ROOT = OUT.parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
import cohort_exact_pilot as pilot
import exact_cohort_context as context
import mixed_geometry_certificate as geom1
import mixed_geometry_certificate_v2 as geom2

PINS = {
    "prepare_records.py": "c36d8c72271a355def1281b1b2bb3a81b9ecc5aa435a81806f4b2adfc28eae1e",
    "protocol.json": "36f58fe4a7b4c100add489b1608b33a264fd35b3d2f30d23d9e03fe2f1a640c0",
    "preflight.json": "2c7a39c9d2e225afe4cd0f25e8172d5f2c8d12c8820e2a10d4f527bb19143f3d",
    "manifest.json": "93a2f5858a450ca9a8db9f9bc944b918e2a7446abba614f97de5e47fe434495d",
}
DRIVER = "a21b680117f236cc38efc46e086be1d9aff80cfb09ae613691b1642db8bed25e"
AUTH = "68164b3a922854b05d55c57d6872e92451e176484d7f5bf2ba3127101819850b"
INTEGRATOR_PROTOCOL = "runs/exact_integration_cohort_v1/integrator_preparation/protocol.json"
GEOM = {
    "scripts/root_body_geometry.py": "18013c6c50c5f679eb6c4abc04e3cf9ce447ee18d0bf1efbf6dcfbd6806f49c2",
    "scripts/mixed_geometry_certificate_v2.py": "165d3b914b1ab9bad2e6591a7df63c09e004f33a742cfa5f4430584746edfbd2",
    "scripts/mixed_geometry_certificate.py": "fb0d49d680286ff49d0dc6713baca41fbebbf724442bff9731292e1441b03fd8",
}
ORDER = ("apical_area_modulus", "basal_area_modulus", "apical_surface_tension", "basal_surface_tension", "common_gain")
PREFIX = "optoGEF2_small_col_"
LETTERS = set("BCDEFGHIJKLMNOPQR")
REMAINING = sorted(PREFIX + c for c in LETTERS - {"B", "L"})


def require(test, msg):
    if not test:
        raise ValueError(msg)


def digest(raw):
    return sha256(raw).hexdigest()


def raw(path):
    absolute = (ROOT / path).resolve()
    require(absolute.is_relative_to(ROOT), f"foreign path {path}")
    return absolute.read_bytes()


@lru_cache(None)
def h(path):
    return digest(raw(path))


def check_hash(path, expected):
    require(h(path) == expected, f"changed file {path}")


def load(path, expected=None):
    data = raw(path)
    if expected is not None:
        require(digest(data) == expected, f"changed JSON {path}")
    return json.loads(data)


def yes(d):
    return isinstance(d, dict) and bool(d) and all(v is True for v in d.values())


def interval(b):
    lo, hi = Q(b["lower_exact"]), Q(b["upper_exact"])
    require(lo <= hi, "invalid interval")
    return lo, hi


def fraction(x):
    return Q.from_float(float(x))


def relative(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def one(rows, record):
    (row,) = [v for v in rows if v["record_id"] == record]
    return row


def write(path, d):
    b = (json.dumps(d, indent=2, allow_nan=False) + "\n").encode()
    with path.open("xb") as stream:
        stream.write(b)
    return digest(b)


def rebuild_box(a, u, S, bases, cell, r, box):
    free = list(map(int, a["free"]))
    require(len(free) == len(set(free)) == len(box["shift_exact"]) == len(box["predictor_hulls_exact"]) == 611
            and all(0 <= k < 650 for k in free) and box["free"] == free
            and all(Q(t) == 0 for t in box["shift_exact"]), "selected free map and zero shift")
    # Enumerate parameter vertices independently of the preparation's per-column extremum formula.
    directions = [tuple(Q(t) - q for t in cell[k]) for k, q in zip(ORDER, bases)]
    dof = {i: j for j, i in enumerate(free)}
    n_axis = 0
    require(all(np.shape(box[k]) == (325, 2) for k in
                ("reference_nodes_exact", "current_centers_exact", "coordinate_radii_exact"))
            and np.shape(box["predictor_coordinate_hulls_exact"]) == (325, 2, 2), "box layout")
    for n in range(325):
        for k in range(2):
            x = fraction(a["reference"][n, k])
            y = Q(box["current_centers_exact"][n][k])
            width = Q(box["coordinate_radii_exact"][n][k])
            lo, hi = map(Q, box["predictor_coordinate_hulls_exact"][n][k])
            require(Q(box["reference_nodes_exact"][n][k]) == x, "reference coordinate")
            j = dof.get(2*n + k)
            if j is None:
                require(x == y == lo == hi and width == 0, "constrained coordinate")
                n_axis += (x == 0 and k == 0)
                continue
            require(x != 0 or k != 0, "axis radial coordinate must stay constrained")
            vertex_values = [Q(0)]
            for coefficient, (dl, dh) in zip(map(fraction, S[j]), directions):
                vertex_values = [q + coefficient * v for q in vertex_values for v in (dl, dh)]
            mn, mx = min(vertex_values), max(vertex_values)
            require((mn, mx) == tuple(map(Q, box["predictor_hulls_exact"][j]))
                    and y == x + fraction(u[j]) and (lo, hi) == (y + mn, y + mx)
                    and y - width <= lo - r <= hi + r <= y + width
                    and width == max(abs(mn), abs(mx)) + r, "exact affine vertex hull and original radius")
    return n_axis


def audit_record(expected, row, item, authority, deps, input_hashes):
    files = {}
    folder = pilot.outdir(expected)
    names = {"protocol.json", "preflight.json", "coordinate_box.json"}
    require({p.name for p in folder.iterdir()} == names and set(item["files"]) == names,
            "record directory must be preparation-only")
    for key, f in item["files"].items():
        require(f["file"] == relative(folder / key) and f["sha256"] == input_hashes[f["file"]],
                "manifest/file/preflight hash binding")
        files[key] = load(f["file"], f["sha256"])
    p, pre, box = [files[k] for k in ("protocol.json", "preflight.json", "coordinate_box.json")]
    ph = item["files"]["protocol.json"]["sha256"]
    box_pin = item["files"]["coordinate_box.json"]
    pi = p["source_and_input_hashes"]
    require(all(pi.get(path) == val for path, val in deps.items()), "row missing pinned reachable source/dependency")
    for path, val in pi.items():
        require(input_hashes[path] == val, "batch preflight versus record input hash")
        check_hash(path, val)
    require(pre["status"] == "PASS_PREPARATION_ONLY" and yes(pre["checks"])
            and pre["source_sha256"] == DRIVER and pre["protocol_sha256"] == ph
            and pre["actual_body_or_force_job_run"] is False
            and p["row"] == item["selected_context"] == row and p["record_id"] == box["record_id"] == expected
            and p["coordinate_box_file"] == box_pin["file"] and p["coordinate_box_sha256"] == box_pin["sha256"]
            and pi[box_pin["file"]] == box_pin["sha256"] and p["root_asserted"] is False
            and box["root_membership_assumed"] is False and row["new_exact_target_root_certified"] is False
            and box["domain_id"] == f"exactFE:{expected}:original_cube_v1", "preparation and phase binding")
    b = {}
    for key in ("root_certificate", "point_cache", "proof_center_state", "predictor", "independent_root_capsule"):
        path, sha = row[key + "_file"], row[key + "_sha256"]
        require(authority["input_hashes"][path] == pi[path] == sha, "selected path provenance")
        b[key] = raw(path) if key == "predictor" else load(path, sha)
    root, state, packet = [b[k] for k in ("root_certificate", "proof_center_state", "point_cache")]
    uresult = root["uniform_cell_result"]
    original_path = (Path(row["root_certificate_file"]).parents[2] / "protocol.json").as_posix()
    require(pi[original_path] == root["protocol_sha256"], "source protocol")
    original = load(original_path, root["protocol_sha256"])
    original_row = one(original["records"], expected)
    contract_file = original["common_model_observation_contract_file"]
    contract = load(contract_file, pi[contract_file])
    cell, base = row["parameter_cell"], row["base_parameters_exact"]
    bases = list(map(Q, (base[k] for k in ORDER)))
    canonical = {k: [str(Q(t)) for t in cell[k]] for k in ORDER}
    cell_hash = digest(json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode())
    require(cell_hash == row["parameter_cell_sha256"] == p["parameter_cell_sha256"] == box["parameter_cell_sha256"]
            == uresult.get("parameter_cell_sha256", cell_hash),
            "selected parameter cell digest")
    require(cell == contract["parameter_cell"] == root["parameter_cell"] == original["shared_parameter_cell"]
            == uresult["cell_bounds_exact"] and base == original["base_parameters_exact"] == item["exact_base"]
            and all(Q(cell[k][0]) <= Q(base[k]) <= Q(cell[k][1]) for k in ORDER)
            and box["model_observation_contract_sha256"] == pi[contract_file] == row["model_observation_contract_sha256"]
            and root["record_id"] == expected and root["status"] == "CERTIFIED_UNIFORM_LOCAL_ROOT"
            and original_row["branch_domain_id"] == root["branch_domain_id"] == row["branch_domain_id"]
            and root["point_cache_file"] == row["point_cache_file"] and root["predictor_file"] == row["predictor_file"]
            and all(authority["input_hashes"][k] == v for k, v in root["artifact_hashes"].items()),
            "same target and closed cell")
    for path, value in root["artifact_hashes"].items():
        check_hash(path, value)
    dose = Q(row["dose_exact"])
    drive = dose * bases[-1]
    require(dose > 0 and drive > 0 and dose == Q(root["dose_exact"]) == Q(uresult["dose_exact"])
            == Q(original_row["dose_exact"]) == Q(item["exact_dose"])
            and drive == Q(root["base_strength_exact"]) == Q(original_row["exact_requested_strength"])
            == Q(item["exact_drive"]) and bases[-1] == Q(original_row["base_gain_exact"])
            and Q(packet["activation_strength_binary_exact"]) == Q(original_row["solver_strength_binary_exact"]),
            "selected exact gain, dose and drive")
    recentered = expected[-1] in "FK"
    center_path = row["proof_center_state_file"]
    require(packet["input_state"] == center_path and packet["input_state_sha256"] == row["proof_center_state_sha256"]
            and expected == packet.get("record_id", expected), "same-record cache binds selected state")
    if recentered:
        require(root["proof_center_state_file"] == center_path != original_row["saved_state"]
                and root["original_seed_state_file"] == state["original_seed_state_file"] == original_row["saved_state"]
                and root["original_seed_state_sha256"] == state["original_seed_state_sha256"] == original_row["saved_state_sha256"]
                and Q(state["proof_activation_strength_exact"]) == drive == Q(packet["exact_base_gain_times_dose"]),
                "selected recentered F/K proof center and drive")
    else:
        require(root["saved_state"] == center_path == original_row["saved_state"], "selected original state")
    embedded = packet["frozen_binary_cache"]
    data = b64decode(embedded["data"], validate=True)
    require(len(data) == embedded["bytes"] and digest(data) == embedded["sha256"] == row["frozen_binary_cache_sha256"],
            "embedded cache integrity")
    with np.load(BytesIO(data), allow_pickle=False) as archive:
        a = {k: archive[k] for k in archive.files}
    with np.load(BytesIO(b["predictor"]), allow_pickle=False) as archive:
        require(archive.files == ["S"], "predictor payload")
        S = archive["S"]
    R, u = a["preconditioner"], a["u"]
    require(embedded["array_shapes"] == {k: list(v.shape) for k, v in a.items()}
            and a["reference"].shape == (325, 2) and a["free"].shape == (611,)
            and R.shape == (611, 611) and S.shape == (611, 5) and u.shape == (611,)
            and all(v.dtype == np.float64 and np.isfinite(v).all() for v in (R, S, u))
            and np.array_equal(u, state["free_displacement"])
            and a["triangles"].tolist() == state["triangles"] == box["triangles"]
            and np.array_equal(a["reference"], state["reference_nodes"]), "binary arrays match selected proof center")
    require(all(digest(v.tobytes()) == item[k + "_binary_sha256"]
                for k, v in (("R", R), ("S", S), ("u", u))), "manifest binary hashes")
    ctx = context.load_context(expected)
    require(ctx.row == row and ctx.drive == drive and ctx.dose == dose
            and all(np.array_equal(c, v) for c, v in ((ctx.R,R),(ctx.S,S),(ctx.u,u))),
            "loaded runtime context agrees with independent cache read")
    audit = one(b["independent_root_capsule"]["record_audits"], expected)
    source, audit_sha = audit["review_file"], audit["review_sha256"]
    require(pi[source] == audit_sha, "selected independent cached proof input")
    proof = load(source, audit_sha)
    require(proof == audit["audit"] and yes(proof["checks"]) and proof["uniform_quantifier_verified"] is True
            and proof["certificate_sha256"] == row["root_certificate_sha256"], "selected independent cached proof")
    if expected.endswith("_M"):
        source, nested_pin = proof["reuse_review_file"], proof["reuse_review_sha256"]
        require(source == "evidence/parametric_actual_M_independent_review.json" and pi[source] == nested_pin,
                "M actual nested provider")
        nested = load(source, nested_pin)
        nested_audit = nested["root_audit"]
        require(nested["status"].startswith("PASS") and yes(nested["checks"]) and yes(nested_audit["checks"])
                and nested_audit["branch_domain_id"] == row["branch_domain_id"]
                and nested_audit["parameter_cell_sha256"] == cell_hash and Q(nested["base_strength_exact"]) == drive,
                "same M nested proof domain and drive")
        n = nested["independent_bounds"]
        bounds = (n["radius"], n["eta_uniform"], n["Z_uniform"])
    else:
        bounds = tuple(proof["bounds"][k] for k in ("radius", "eta", "Z"))
    ball, eta, z = map(interval, bounds)
    eta = max(eta[1], interval(uresult["eta_uniform"])[1])
    z = z[1]
    r = Q(row["selected_corrector_radius_exact"])
    require(ball[0] <= r <= ball[1] and 0 < r and 0 <= eta and 0 <= z < 1 and eta + z*r < r
            and r == Q(uresult["selected_radius_exact"]) == Q(box["corrector_radius_exact"])
            == Q(p["normal_radius_exact"]) == Q(item["selected_corrector_radius_exact"])
            and eta == Q(p["cached_eta_upper_exact"]) == Q(item["exact_cached_eta_upper"])
            and z == Q(p["cached_Z_upper_exact"]) == Q(item["exact_cached_Z_upper"])
            and source == p["cached_bound_provider"] == item["audited_radius_and_cached_eta_Z_provider"],
            "audited radius interval containment and cached proof bounds")
    n_axis = rebuild_box(a, u, S, bases, cell, r, box)
    require(n_axis == item["fixed_axis_coordinates"] and item["exact_box_coordinate_checks"] == 650,
            "all exact box coordinates checked")
    incidences, h_inc = np.zeros(611, dtype=int), np.zeros(611, dtype=int)
    mesh = (a["triangles"], a["surface_0_edges"], a["surface_1_edges"])
    require(tuple(v.shape for v in mesh) == ((144,6),(12,3),(12,3)), "element and edge coverage")
    d = {i: j for j, i in enumerate(map(int, a["free"]))}
    for group in (v for m in mesh for v in m):
        projection = [d[i] for n in group for i in (2*int(n), 2*int(n)+1) if i in d]
        require(projection and len(set(projection)) == len(projection), "local free map")
        incidences[projection] += 1
        h_inc[projection] += len(projection)
    rn = max(sum(map(abs, map(fraction, row)), Q(0)) for row in R)
    require(rn == Q(p["R_infinity_norm_exact"]) == Q(item["exact_R_infinity_norm"]) > 0
            and incidences.tolist() == p["force_incidence"] and h_inc.tolist() == p["H_row_incidence"]
            and 0 < max(incidences) == p["max_force_incidence"] and 0 < max(h_inc) == p["max_H_row_incidence"],
            "all force and Hessian incidence entries and exact R norm")
    f, he = Q(p["force_tail_budget_exact"]), Q(p["H_entry_tail_budget_exact"])
    require(f == Q(item["force_tail_budget_exact"]) > 0 and he == Q(item["H_entry_tail_budget_exact"]) > 0
            and p["margin_divisor"] == 10000 and p["precision_bits"] == 640
            and (p["max_bulk_order"], p["max_surface_order"], p["surface_domain_depth"]) == (32,64,8)
            and rn * int(max(incidences)) * f == Q(p["global_force_tail_budget_exact"]) == (1-z)*r/10000
            and rn * int(max(h_inc)) * he == Q(p["global_derivative_tail_budget_exact"]) == (1-z)/10000,
            "positive exact tail budgets and domain selector policy")
    geom = {}
    for key, shape in {"shape":(64,6), "shape_grad":(144,64,6,2),"volume_weights":(144,64),
                       "surface_0_N":(12,3), "surface_1_N":(12,3),
                       "surface_0_weights":(12,12), "surface_1_weights":(12,12)}.items():
        require(a[key].shape == shape and np.isfinite(a[key]).all(), f"cached rule array {key}")
        if key.endswith("_weights"):
            require((a[key] > 0).all(), f"cached positive weights {key}")
    require(fraction(a["energy_scale"].item()) == Q(item["normalized_profile"]["energy_scale_exact"]) > 0
            and ctx.target.energy_scale == fraction(a["energy_scale"].item())
            and set(np.unique(a["active"])) == {0.0,1.0} and np.count_nonzero(a["active"]) == 24
            and packet["parameters"]["mu"] == 1 and packet["parameters"]["bulk"] == 5000,
            "positive normalization and fixed bulk material")
    for key, v in a.items():
        if key.startswith("surface_") or key in ("reference","free","triangles","shape","shape_grad","active","qr","energy_scale","volume_weights"):
            geom[key] = digest(v.tobytes())
    record = {
        "status":"PASS", "phase":"domains_only", "record_id":expected,
        "source_sha256":DRIVER, "protocol_sha256":ph,
        "preflight_sha256":item["files"]["preflight.json"]["sha256"],
        "coordinate_box_sha256":box_pin["sha256"],
        "checks": dict.fromkeys(("current_dependencies_and_phase_guard","authoritative_record_cell_dose_and_base",
                                 "selected_state_cache_u_R_and_predictor_S","cached_bound_provider_and_original_radius_contains",
                                 "independent_exact_vertex_box_650","tail_budgets_and_fixed_rule"), True),
        "checked_coordinate_count":650, "constrained_zero_axis_coordinates":n_axis,
        "r_exact":str(r), "cached_radius_ball_exact":list(map(str, ball)),
        "provider":source, "recentered_proof_center":recentered,
        "F_tail_exact":str(f), "H_entry_tail_exact":str(he),
        "runtime_record_preparation_files": item["files"], "batch_manifest_sha256":PINS["manifest.json"],
        "audit_source_sha256":digest(Path(__file__).read_bytes()), "integration_authorized":False,
        "numerical_jobs_run_by_review":False,
    }
    return folder, record, geom


def run():
    batch_rel = relative(BATCH)
    for name, pin in PINS.items():
        check_hash(f"{batch_rel}/{name}", pin)
    check_hash(f"{batch_rel}/prepare_records_frozen.py", PINS["prepare_records.py"])
    p, pre, m = [load(f"{batch_rel}/{n}") for n in ("protocol.json", "preflight.json", "manifest.json")]
    require(pre["status"] == "PASS_REMAINING15_PREPARATION_ONLY" and yes(pre["checks"])
            and p["source_sha256"] == pre["source_sha256"] == m["frozen_source_sha256"] == PINS["prepare_records.py"]
            and pre["protocol_sha256"] == m["protocol_sha256"] == PINS["protocol.json"]
            and pre["manifest_sha256"] == PINS["manifest.json"]
            and m["status"] == "READY_FOR_INDEPENDENT_RECORD_DOMAIN_PRELAUNCH_REVIEWS" and not m["rejected_attempts"]
            and pre["actual_new_body_or_integral_or_cache_jobs_run"] is m["actual_new_body_or_integral_or_cache_jobs_run"] is False,
            "write-once batch preparation state")
    require(context.AUTHORITY_SHA256 == AUTH, "authority selection")
    authority = load(context.AUTHORITY, AUTH)
    ids = [r["record_id"] for r in authority["records"]]
    require(len(ids) == len(set(ids)) == 17 and set(ids) == {PREFIX+c for c in LETTERS}
            and p["selection"] == m["expected_record_ids"] == [r["record_id"] for r in m["records"]] == REMAINING
            and p["excluded_completed"] == m["excluded_completed_records"] == sorted((PREFIX+"B",PREFIX+"L")),
            "all 17 minus the two certified IDs")
    dependencies = load(INTEGRATOR_PROTOCOL, "cced8eb0a3b7b15d96fe0c4f43da042c807ecd109bcaf1dd37e0a0de1e979096")["source_hashes_and_inputs"]
    deps = {**dependencies, **GEOM, relative(pilot.__file__): DRIVER, context.AUTHORITY: AUTH}
    for module in (pilot.G, geom1, geom2, pilot.E, pilot.S, context):
        check_hash(relative(module.__file__), deps[relative(module.__file__)])
    for path, pin in pre["source_input_and_record_file_hashes"].items():
        check_hash(path, pin)
    require(all(pre["source_input_and_record_file_hashes"][path] == val
                for path, val in p["source_and_authority_hashes"].items()), "bound prior reviewed sources")
    forbidden = AssertionError("domain, integral or cache operation during independent preparation review")
    with (patch.object(pilot, "domains_work", side_effect=forbidden),
          patch.object(pilot.G, "certify_box", side_effect=forbidden),
          patch.object(pilot.E, "prepare_bulk_selector", side_effect=forbidden),
          patch.object(pilot.E, "bulk_component", side_effect=forbidden),
          patch.object(pilot.E, "surface_component", side_effect=forbidden),
          patch.object(pilot.core, "local_terms", side_effect=forbidden)):
        audited = [audit_record(r["record_id"], one(authority["records"], r["record_id"]), r, authority, deps,
                                pre["source_input_and_record_file_hashes"]) for r in m["records"]]
    require(all(geom == audited[0][2] for _, _, geom in audited), "retained shared reference geometry/rule only")
    require(sum(rec["checked_coordinate_count"] for _, rec, _ in audited) == 9750, "exact coordinate count")
    approved = []
    for folder, receipt, _ in audited:
        gate = folder / "domains_prelaunch_independent_review.json"
        gate_sha = write(gate, receipt)
        pilot.gate(folder, "domains")
        require(not (folder/"integration_prelaunch_independent_review.json").exists(), "domains-only gate")
        approved.append({"record_id":receipt["record_id"],"gate_file":relative(gate),"sha256":gate_sha})
    batch = {
        "status":"PASS_REMAINING15_DOMAIN_PRELAUNCH_ONLY",
        "checks":dict.fromkeys(("bound_batch_inputs","all_17_minus_B_and_L","all_9750_exact_coordinates",
                                "all_original_radius_balls_contain_selected_r","same_C_and_selected_dose_gain_R_S_cache",
                                "recentered_F_K_and_nested_actual_M","exact_tail_allocations",
                                "current_bound_geometry_and_kernel_dependencies","domain_gate_readback_all_15"), True),
        "reviewed_batch_pins":PINS, "audit_source_sha256":digest(Path(__file__).read_bytes()),
        "approved_gates":approved, "integration_authorized":False, "numerical_jobs_run_by_review":False,
    }
    sha = write(OUT/"batch_independent_review.json", batch)
    print(json.dumps({"status":batch["status"],"receipt_sha256":sha,"gates":len(approved)}))


if __name__ == "__main__":
    run()

