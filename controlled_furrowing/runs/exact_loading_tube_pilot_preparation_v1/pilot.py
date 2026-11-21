"""Closed K–L secant trial with separate geometry and certificate stages."""

import argparse
import ast
from fractions import Fraction as Q
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
from unittest.mock import patch


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCRIPTS = ROOT / "scripts"
EX = "runs/transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1/"
SPD = EX + "remaining32_preconditioner_spd_v1/independent_actual_preparation/actual_replay/"
DESIGN = "runs/exact_loading_tube_design_v1/"
K, L = (f"max_negative_contrast_change_small_{i:02}" for i in (16, 17))
PILOT_ID = "negative_KL_exact_spatial_loading_secant_v1"
RADIUS = Q(1, 2**12)
N = 611
TARGET_SHA = "dfea134640f99b45e9ca75d1ce718043d67f2917375506997d2465abda8c1242"
INPUTS = {
    DESIGN + "receipt.json": "dbfb9871da2b384fcd57e44015a35bfeebbfccf25e23a6f45178f3aef488ff3c",
    DESIGN + "root_readback.json": "5d97278e0b8157abf33959ca6ad692c1c1907ec06dbee09c11027716efa7bbe6",
    EX + "loading_affine_identity_v1/independent_readback/review_receipt.json":
        "2e77f4dc2a6cd369f6368b3028a408fdd970cf0c9e3f892f08f9e60727a1be9a",
    EX + "remaining32_targets_v1/independent_prelaunch/review.json":
        "67d33bbc3171839b1c578b02b5cfbf9ea6ba5d6cb351d0655dd14aaf99df75f9",
    EX + "remaining32_targets_v1/preparation/manifest.json": TARGET_SHA,
    SPD + "os_exit_observation.json": "fd8ac1cb3b920b57096649544bba105442404fe19576aae99c53f849fe15a8d5",
    SPD + "receipt.json": "de31a30f3270829c7fd89ebc75db58e45d2d9a911cd02dfc6e564493748d8615",
    SPD + K + ".json": "799fb6d40b2fe8f882e4787a0c08b7f6aae22f255a352b405957925b95861a37",
    "scripts/transfer_ambiguity_remaining32_exact_targets.py":
        "d1a9c9e4839ff85a9926b6fd5ecc11e4cd8e0831c7e6263ecfbe068ef22ea3f4",
    "scripts/root_body_geometry.py": "18013c6c50c5f679eb6c4abc04e3cf9ce447ee18d0bf1efbf6dcfbd6806f49c2",
    "scripts/exact_cohort_integrator.py": "092dc51161e708cdba3791ec47a7f8e843dd0e1afc866b6f65f286b8657d1745",
    "scripts/exact_uniform_surface_bridge_B_safe.py":
        "61cf47d2b1a9e7ee46f9fa680c24b4fdf7b4e31e4889379c56951ab6a7455ee9",
}
METHOD_SOURCES = (
    "scripts/algebraic_root_feasibility.py",
    "scripts/exact_scaled_cohort_bulk.py",
    "scripts/exact_centered_bulk_L.py",
    "scripts/exact_signed_predictor_L.py",
    "scripts/transfer_ambiguity_exact_bridge_reduction.py",
    "docs/COHORT-ENERGY-STABILITY.md",
)
NO_NUMERICS = {"new_geometry_certificates": 0, "new_AD_derivatives": 0,
               "new_integrations": 0, "new_solves": 0}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def relative(path):
    return path.resolve().relative_to(ROOT).as_posix()


def canonical(name):
    path = (ROOT / name).resolve()
    require(path.is_relative_to(ROOT) and relative(path) == name,
            "canonical study-relative input path required")
    return path


def digest(path):
    with path.open("rb") as stream:
        return hashlib_file(stream)


def hashlib_file(stream):
    digestor = sha256()
    for blob in iter(lambda: stream.read(1024 * 1024), b""):
        digestor.update(blob)
    return digestor.hexdigest()


def hashed_json(path, expected):
    payload = path.read_bytes()
    require(sha256(payload).hexdigest() == expected, f"digest mismatch for {relative(path)}")
    return json.loads(payload)


def save(path, value):
    with path.open("x", encoding="utf-8", newline="\n") as output:
        json.dump(value, output, indent=2, default=str, allow_nan=False)
        output.write("\n")
    return digest(path)


def exact(value):
    require(type(value) in (str, int, Q), "exact rational required")
    return Q(value)


def pair(value):
    require(isinstance(value, (list, tuple)) and len(value) == 2,
            "two rational interval endpoints required")
    a, b = map(exact, value)
    require(a <= b, "reversed interval")
    return a, b


def read_member(hashes, expected, suffix):
    names = [name for name, value in hashes.items() if value == expected and name.endswith(suffix)]
    require(len(names) == 1, f"one pinned {suffix} required")
    return names[0], hashed_json(canonical(names[0]), expected)


def int_list(values, low, high, count):
    require(isinstance(values, list) and len(values) == count
            and all(type(i) is int and low <= i < high for i in values)
            and len(set(values)) == count, "strict integral projection map required")
    return values


def node_rows(raw, count):
    require(isinstance(raw, list) and len(raw) == count
            and all(isinstance(row, list) and len(row) == 2 for row in raw),
            "complete two-coordinate node array required")
    return [tuple(map(exact, row)) for row in raw]


def maxabs(values):
    return max(map(abs, values), default=Q(0))


def secant_box(ref, yk, yl, free, radius, gap):
    """Return exact nodal envelope and displacement-coordinate predictor."""
    require(type(radius) is Q and radius > 0 and type(gap) is Q and gap > 0,
            "positive exact corrector radius and dose gap required")
    require(len(ref) == len(yk) == len(yl), "matching node arrays required")
    active = set(free)
    centers, radii, first, last, slopes = [], [], [], [], []
    for node, (x, a, b) in enumerate(zip(ref, yk, yl, strict=True)):
        yc, dy = [], []
        for comp, (xj, aj, bj) in enumerate(zip(x, a, b, strict=True)):
            coord = 2*node + comp
            if coord in active:
                first.append(aj-xj)
                last.append(bj-xj)
                slopes.append((bj-aj)/gap)
                margin = radius
            else:
                require(aj == bj == xj, "constrained coordinate motion")
                margin = Q(0)
            yc.append((aj+bj)/2)
            dy.append(abs(bj-aj)/2+margin)
        centers.append(yc)
        radii.append(dy)
    require(len(first) == len(free) and list(free) == sorted(active),
            "ordered complete free coordinate map required")
    return centers, radii, first, last, slopes


def contains_root(center_a, error_a, center_b, radius_b):
    require(exact(error_a) >= 0 and exact(radius_b) > 0, "invalid containment radii")
    distance = maxabs(exact(a)-exact(b) for a,b in zip(center_a, center_b, strict=True))
    return exact(radius_b) - distance - exact(error_a)


def quad_bounds(offset, slope, quadratic, a, b):
    """Exactly enclose quadratic family on a closed nonnegative t interval."""
    require(0 <= a <= b, "ordered nonnegative offsets required")
    u, s, k = pair(offset), pair(slope), pair(quadratic)

    def candidates(i):
        times = [a, b]
        if k[i] != 0:
            turning = -s[i]/(2*k[i])
            if a < turning < b:
                times.append(turning)
        return (u[i] + t*(s[i]+t*k[i]) for t in times)

    return min(candidates(0)), max(candidates(1))


def conclude(q, hc_v, active_center, active_v, zc, active_norm, rnorm, gap, radius):
    require(zc >= 0 and active_norm >= 0 and rnorm > 0 and gap > 0 and radius > 0,
            "invalid derivative, matrix norm, or domain bounds")
    residual = []
    slope = []
    for row, ha, ac, k in zip(q, hc_v, active_center, active_v, strict=True):
        s = tuple(a+b for a,b in zip(pair(ha), pair(ac), strict=True))
        slope.append(s)
        residual.append(quad_bounds(row, s, k, Q(0), gap))
    require(bool(residual), "nonempty signed residual rows required")
    eta = maxabs(x for row in residual for x in row)
    z = zc + gap*active_norm
    lhs = eta+z*radius
    passed = z < 1 and lhs < radius
    return {"strict_fiber_selfmaps": passed, "eta_upper_exact": str(eta),
            "z_upper_exact": str(z), "strict_margin_exact": str(radius-lhs),
            "selfmap_lhs_upper_exact": str(lhs),
            "posterior_error_upper_exact": str(eta/(1-z)) if passed else None,
            "whole_product_normalized_H_euclidean_lower_exact": str((1-z)/rnorm) if passed else None,
            "signed_R_Hv_plus_active_center_rows_exact": slope,
            "signed_RF_predictor_all_loads_rows_exact": residual}


def local_source_tree():
    seeds = [SCRIPTS / "transfer_ambiguity_remaining32_exact_targets.py"]
    hashes = {}
    while seeds:
        path = seeds.pop().resolve()
        name = relative(path)
        if name in hashes:
            continue
        hashes[name] = digest(path)
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            modules = ([alias.name for alias in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) else [])
            for module in modules:
                candidate = SCRIPTS / f"{(module or '').split('.')[0]}.py"
                if candidate.is_file():
                    seeds.append(candidate)
    return dict(sorted(hashes.items()))


def cold_plan():
    documents = {name: hashed_json(canonical(name), value)
                 for name, value in INPUTS.items() if name.endswith(".json")}
    for name, value in INPUTS.items():
        require(digest(canonical(name)) == value, f"fixed input changed: {name}")
    cross = documents[DESIGN+"receipt.json"]
    root = documents[DESIGN+"root_readback.json"]
    aff = documents[EX+"loading_affine_identity_v1/independent_readback/review_receipt.json"]
    rev = documents[EX+"remaining32_targets_v1/independent_prelaunch/review.json"]
    require(cross["joint_path_certified"] is False and root["loading_path_certified"] is False
            and root["K_L_physical_crosswalk_independently_reviewed_separately"] is True
            and root["receipt_sha256"] == INPUTS[DESIGN+"receipt.json"]
            and aff["status"] == "PASS_INDEPENDENT_AFFINE_LOADING_SOURCE_AND_MATHEMATICS"
            and rev["status"] == "PASS_INDEPENDENT_REMAINING32_TARGET_CONTEXTS"
            and rev["manifest_sha256"] == TARGET_SHA, "accepted source and K/L crosswalk required")
    hashes = cross["input_sha256_by_study_relative_path"]
    for path, h in hashes.items():
        require(digest(canonical(path)) == h, f"original crosswalk dependency changed: {path}")
    p = cross["pilot"]
    require(p["case_order"] == [K,L] and p["exact_active_direction"] == "e_Z"
            and len(int_list(p["active_slots"],0,144,24)) == 24, "fixed material/load trial required")
    kpath, kb = read_member(hashes, p["K_box_sha256"], f"/{K}.json")
    lpath, lb = read_member(hashes, p["L_box_sha256"], "/coordinate_box_1.json")
    knative_path, kn = read_member(hashes, p["K_spatial_proof_sha256"], f"/{K}.json")
    lnative_path, ln = read_member(hashes, p["L_spatial_proof_sha256"], "/receipt.json")
    manifest = documents[EX+"remaining32_targets_v1/preparation/manifest.json"]
    held = [row for row in manifest["targets"] if row["case_id"] == K]
    require(len(held) == 1, "unique reviewed K context required")
    row = held[0]
    kd, ld = pair(p["exact_dose_interval"])
    gap, gain, energy = map(exact, (p["exact_dose_gap"], p["same_gain_exact"],
                                    p["common_binary_energy_normalization_exact"]))
    require(gap == ld-kd > 0 and gain > 0 and energy > 0
            and list(map(exact,p["exact_requested_drives"])) == [gain*kd,gain*ld]
            and [exact(b["dose_exact"]) for b in (kb,lb)] == [kd,ld]
            and [b["record_id"] for b in (kb,lb)] == ["transfer_exact_"+i for i in (K,L)],
            "same exact material, dose and strength identity required")
    ref = node_rows(kb["reference_nodes_exact"],325)
    free = int_list(kb["free"],0,650,N)
    tri = kb["triangles"]
    require(len(tri)==144 and all(int_list(t,0,325,6) for t in tri)
            and lb["triangles"]==tri and lb["free"]==free
            and node_rows(lb["reference_nodes_exact"],325)==ref,
            "common 325-node, 144-element, 611-free-coordinate P2 reference required")
    require(all(exact(b["corrector_radius_exact"]) == Q(1,2**34)
                and b["root_membership_assumed"] is False for b in (kb,lb)), "exact archived endpoint cubes")
    target = kn["target"]
    sym = kn["same_R_positive_symmetric_part_premise"]
    spd = documents[SPD+K+".json"]
    spdv = documents[SPD+"receipt.json"]
    spdo = documents[SPD+"os_exit_observation.json"]
    spdt = spd["target_identity"]
    srows = [x for x in spdv["case_results"] if x["case_id"]==K]
    errors = list(map(exact,p["endpoint_posterior_error_upper_exact"]))
    require(kn["status"] == "PASS_INDEPENDENT_FIRST_WAVE_EXACT_SPATIAL_ROOT"
            and kn["arithmetic"]["strict_selfmap"] is True
            and kn["observer_and_optional_minimum"]["strict_minimum_proved"] is True
            and target["trial_box_sha256"] == p["K_box_sha256"]
            and row["model_contract_sha256"] == target["model_contract_sha256"]
            and target["point_cache_sha256"] == kb["point_cache_sha256"]
            and exact(target["exact_drive"]) == gain*kd
            and exact(target["energy_scale_binary_exact"]) == energy
            and exact(kn["arithmetic"]["posterior_error_upper_exact"]) == errors[0]
            and ln["status"] == "PASS_INDEPENDENT_EXACT_SPATIAL_ROOT_AND_STRICT_LOCAL_MINIMUM"
            and ln["pilot_id"] == L and ln["independent_bounds"]["strict_selfmap"] is True
            and 0 <= exact(ln["independent_bounds"]["posterior_error_upper_exact"]) <= errors[1]
            and exact(ln["conservative_compared_bounds"]["posterior_error_upper_exact"]) == errors[1]
            and all(0 <= e < Q(1,2**34) for e in errors)
            and spdo["status"] == "PASS_ROOT_OBSERVED_REMAINING32_SPD_REPLAY_EXIT"
            and type(spdo["observed_exit_code"]) is int and spdo["observed_exit_code"]==0
            and positive_int(spdo["pid"]) and positive_int(spdo["unified_exec_session"])
            and spdo["receipt_sha256"] == sym["exact_Sym_R_native_receipt_sha256"] == INPUTS[SPD+"receipt.json"]
            and sym["exact_Sym_R_native_case_sha256"] == INPUTS[SPD+K+".json"]
            and spdv["status"] == "PASS_INDEPENDENT_ACTUAL_REMAINING32_EXACT_SYMMETRIC_PARTS"
            and len(srows)==1 and srows[0]["sha256"] == INPUTS[SPD+K+".json"]
            and spd["status"] == "PASS_INDEPENDENT_REMAINING32_EXACT_SYMMETRIC_PART_SPD"
            and spd["strict_congruence_bound"] is True and exact(spd["independent_q_infinity_upper_exact"])<1
            and spd["trial_box_sha256"] == target["trial_box_sha256"]
            and spd["exact_R_infinity_norm"] == kn["same_R_norm_exact"]
            and spdt["R_binary_sha256"] == target["R_binary_sha256"]
            and spdt["center_binary_sha256"] == target["center_binary_sha256"]
            and spdt["cache_packet_sha256"] == target["point_cache_sha256"]
            and exact(spdt["drive_exact"])==gain*kd and exact(spdt["energy_scale_binary_exact"])==energy,
            "same signed K point residual and independently SPD K binary R required")
    yk, yl = (node_rows(b["current_centers_exact"],325) for b in (kb,lb))
    centers, radii, uk, ul, slope = secant_box(ref,yk,yl,free,RADIUS,gap)
    require(all(a+gap*v==b for a,b,v in zip(uk,ul,slope,strict=True)), "secant must hit L proof center")
    q = [pair(v) for v in kn["arithmetic"]["signed_R_F_I_rows_exact"]]
    require(len(q)==N and maxabs(x for v in q for x in v) == exact(kn["arithmetic"]["eta_upper_exact"]),
            "all preconditioned exact point force rows retained")
    linear = {"reference_nodes_exact":kb["reference_nodes_exact"], "triangles":tri,"free":free,
              "current_centers_exact":centers,"coordinate_radii_exact":radii,
              "root_membership_assumed":False,"any_spatial_derivative_bound_assumed":False}
    plan = {"status":"PREPARED_CLOSED_KL_SECANT_PILOT", "pilot_id":PILOT_ID,
            "case_ids":[K,L], "material_parameters_binary_exact":p["exact_material_parameters"],
            "exact_active_direction":p["exact_active_direction"], "active_slots":p["active_slots"],
            "dose_interval_exact":[kd,ld],"nominal_dose_exact":kd, "dose_gap_exact":gap,
            "gain_exact":gain, "nominal_drive_exact":gain*kd,"energy_scale_binary_exact":energy,
            "fiber_radius_exact":RADIUS,"anchor_box_files":[kpath,lpath],
            "anchor_box_sha256":[p["K_box_sha256"],p["L_box_sha256"]],
            "anchor_spatial_proof_files":[knative_path,lnative_path],
            "anchor_spatial_proof_sha256":[p["K_spatial_proof_sha256"],p["L_spatial_proof_sha256"]],
            "anchor_posterior_error_upper_exact":errors,"K_binding":target,"K_same_R_SymR":sym,
            "same_R_infinity_norm_exact":kn["same_R_norm_exact"],"K_routes_file":row["routes"],
            "K_routes_sha256":row["routes_sha256"], "K_archive_shift_unchanged_exact":p["K_original_observer_shift_preserved_exact"],
            "K_center_displacement_exact":uk,"L_center_displacement_exact":ul,"predictor_slope_exact":slope,
            "signed_K_R_F_I_point_rows_exact":q,
            "max_center_motion_exact":gap*maxabs(slope), "max_slope_exact":maxabs(slope),
            "largest_envelope_coordinate_radius_exact":max(r for a in radii for r in a),
            "anchor_posterior_containment_margins_exact":[RADIUS-e for e in errors],
            "envelope_is_cartesian_hull_of_closed_fibers":True,"only_K_center_point_force_reused":True,
            "cached_or_exact_endpoint_H_used_on_envelope":False,"path_certified":False}
    return plan, linear, len(hashes)


GEOMETRY_CHECKS = (
    "same_exact_K_center_R_drive_and_KL_physical_reference",
    "secant_envelope_contains_all_closed_fibers_and_radial_bounds_checked",
    "fresh_body_charts_and_trace_geometry_only",
    "write_once_run_and_terminal_with_no_subprocesses",
)
MECHANICS_CHECKS = (
    "same_closed_nodal_envelope_exact_nominal_drive_material_and_normalization",
    "same_binary_preconditioner_and_exact_symmetric_part_premise",
    "fresh_independent_full_envelope_body_traces_and_cached_analytic_domain",
    "fresh_exact_integral_uniform_derivative_and_signed_H_directional_action",
    "exact_rational_axial_Gram_coefficients_and_directed_pi",
    "independently_replayed_actual_operands_with_observed_zero_OS_exits",
)


def bindings():
    return {label+"_sha256":digest(path) for label,path in {
        "source":Path(__file__),"protocol":HERE/"protocol.json","preflight":HERE/"preflight.json",
        "plan":HERE/"plan.json","envelope":HERE/"envelope.json"}.items()}


def yes(required, document):
    checks = document["checks"]
    return isinstance(checks,dict) and all(checks.get(flag) is True for flag in required)


def gate_geometry(document, current):
    require(document["status"] == "PASS_INDEPENDENT_KL_LOADING_GEOMETRY_PRELAUNCH"
            and document["pilot_id"] == PILOT_ID and document["bindings"] == current
            and yes(GEOMETRY_CHECKS,document)
            and document["geometry_only_authorized"] is True
            and document["mechanical_kernels_authorized"] is False,
            "exact current phase plus independent geometry-only release required")


def prepared():
    current = bindings()
    protocol = hashed_json(HERE/"protocol.json",current["protocol_sha256"])
    preflight = hashed_json(HERE/"preflight.json",current["preflight_sha256"])
    stable = ("source_sha256", "plan_sha256", "envelope_sha256")
    require(protocol["status"] == "COLD_KL_LOADING_TUBE_PILOT_PROTOCOL"
            and preflight["status"] == "PASS_KL_LOADING_COLD_FIXTURES"
            and all(protocol["bindings"][k]==preflight["bindings"][k]==current[k] for k in stable)
            and preflight["bindings"]["protocol_sha256"]==current["protocol_sha256"]
            and protocol["geometry_or_mechanics_launch_authorized"] is False
            and all(preflight[k]==v and type(preflight[k]) is int for k,v in NO_NUMERICS.items()),
            "frozen source and passed cold preflight required")
    for sources in ("geometry_source_tree_sha256", "read_only_method_source_sha256", "individual_authority_sha256"):
        for name,h in protocol[sources].items():
            require(digest(canonical(name))==h,f"changed pinned input {name}")
    plan = hashed_json(HERE/"plan.json",current["plan_sha256"])
    env = hashed_json(HERE/"envelope.json",current["envelope_sha256"])
    if protocol["geometry_source_tree_sha256"] != local_source_tree():
        raise ValueError("geometry import paths changed since freezing")
    return current, plan, env, protocol


def geometry_bounds(plan, env):
    """Use the same pure body, chart and trace functions as the reviewed r34 stage."""
    sys.path.insert(0, str(SCRIPTS))
    import transfer_ambiguity_remaining32_exact_targets as A

    require(Path(A.__file__).resolve() == SCRIPTS/"transfer_ambiguity_remaining32_exact_targets.py"
            and Path(A.X.G.__file__).resolve() == SCRIPTS/"root_body_geometry.py"
            and Path(A.E.__file__).resolve() == SCRIPTS/"exact_cohort_integrator.py"
            and Path(A.E.SUR.__file__).resolve() == SCRIPTS/"exact_uniform_surface_bridge_B_safe.py",
            "actual imported body, chart and trace module identity")
    A.X.bootstrap()
    with A.X.no_integrals():
        ctx, _, record = A.open_context(K,TARGET_SHA)
        target = plan["K_binding"]
        proof_center = A.node_box(ctx,A.R34)
        require(ctx.drive == exact(plan["nominal_drive_exact"])
                and ctx.target.energy_scale == exact(plan["energy_scale_binary_exact"])
                and A.a_sha(ctx.R) == target["R_binary_sha256"]
                and A.a_sha(ctx.u) == target["center_binary_sha256"]
                and ctx.row["point_cache_sha256"] == target["point_cache_sha256"]
                and proof_center == hashed_json(canonical(plan["anchor_box_files"][0]),plan["anchor_box_sha256"][0])
                and record["routes"] == plan["K_routes_file"]
                and record["routes_sha256"] == plan["K_routes_sha256"], "same exact K context required")
        return domain_components(A,ctx,plan,env)


def domain_components(A, ctx, plan, box):
    G, E = A.X.G, A.E
    routes = hashed_json(canonical(plan["K_routes_file"]),plan["K_routes_sha256"])
    ref = node_rows(box["reference_nodes_exact"],325)
    y = node_rows(box["current_centers_exact"],325)
    offsets = node_rows(box["coordinate_radii_exact"],325)
    # Match the literal binary64 trial radial bounds in axisymmetric_mixed.py.
    radial_margins = [y[i//2][0] - offsets[i//2][0] - ref[i//2][0]
                      - Q.from_float(-(1-1e-7)*float(ref[i//2][0]))
                      for i in box["free"] if i % 2 == 0]
    floor = min(radial_margins)
    body = G.certify_box(box["reference_nodes_exact"],box["current_centers_exact"],
                         box["triangles"],box["coordinate_radii_exact"])
    if body["status"] != "CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING" or floor <= 0:
        return "UNRESOLVED_KL_WHOLE_ENVELOPE_BODY_OR_RADIAL_BOUND", body, [], [], floor
    bulk, trace = [], []
    require([row["element"] for row in routes["bulk"]] == list(range(144))
            and [(row["face"],row["edge"]) for row in routes["surface"]]
            == [(f,e) for f in range(2) for e in range(12)], "complete canonical bulk and trace route set")
    nodes = E.SUR.common_nodes(ctx.target,box)
    for row in routes["bulk"]:
        slot = row["element"]
        fields, ids = E._bulk_fields(ctx.target,nodes,slot)
        require(ids==row["nodes"] and fields["chart"]["kind"]==row["chart"]
                and fields["free_local"]==row["free_local"] and fields["global_indices"]==row["global_indices"]
                and fields["active"]==exact(row["active"]), "whole-envelope bulk route changed")
        try:
            lo, hi, raw, rho = E._body_rho(body["triangles"][slot])
            bulk.append({**row,"J_bounds_exact":[lo,hi],"raw_rho_exact":raw,"rho_exact":rho})
        except ValueError as error:
            bulk.append({**row,"unresolved":str(error)})
    for row in routes["surface"]:
        f,e = row["face"],row["edge"]
        fields = ctx.target.surface(f,e,nodes=nodes)
        require(fields["free_local"]==row["free_local"] and fields["global_indices"]==row["global_indices"],
                "whole-envelope trace projection changed")
        try:
            bounds = E.SUR.uniform_G_bounds(fields,max_depth=8)
            require(bounds["status"]=="CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN",
                    "whole-envelope trace geometry unresolved")
            sites = [E.SUR.local_site_domain(ctx.target,box,f,e,k) for k in range(12)]
            trace.append({**row,"domain":bounds,"surface_cache_line_geometry":sites})
        except ValueError as error:
            trace.append({**row,"unresolved":str(error)})
    passed = all("rho_exact" in x for x in bulk) and all("domain" in x for x in trace)
    return ("PASS_KL_FRESH_WHOLE_ENVELOPE_GEOMETRY_CANDIDATE" if passed else
            "UNRESOLVED_KL_WHOLE_ENVELOPE_BULK_OR_TRACE_DOMAIN"), body, bulk, trace, floor


def run_geometry(review_sha):
    current, plan, box, _ = prepared()
    review = hashed_json(HERE/"independent_geometry_prelaunch.json",review_sha)
    gate_geometry(review,current)
    out = HERE/"actual_geometry"
    out.mkdir(exist_ok=False)
    identity = {"bindings":current,"pilot_id":PILOT_ID,"pid":os.getpid(),"phase":"geometry_only",
                "independent_geometry_prelaunch_sha256":review_sha}
    save(out/"start.json",identity)
    try:
        status, body, bulk, trace, floor = geometry_bounds(plan,box)
        bh = save(out/"body.json",body)
        rh = save(out/"result.json",{**identity,"status":status,"body_sha256":bh,"bulk":bulk,"traces":trace,
                   "radial_trial_bound_margin_normalized_exact":floor,"new_AD_derivatives":0,
                   "new_integrations":0,"new_solves":0,"root_branch_claimed":False})
    except Exception as error:
        failure = save(out/"failure.json",{**identity,"error_type":type(error).__name__,"error":str(error)})
        save(out/"terminal.json",{**identity,"status":"KL_GEOMETRY_PROCESS_FAILURE","failure_sha256":failure})
        raise
    save(out/"terminal.json",{**identity,"status":"KL_GEOMETRY_STAGE_FINISHED", "result_sha256":rh,
                              "mathematical_status":status})


def positive_int(value):
    return type(value) is int and value>0


def observed_zero(status, observation, result_sha, terminal_sha, pid):
    require(observation["status"] == status and type(observation["observed_exit_code"]) is int
            and observation["observed_exit_code"]==0 and positive_int(pid)
            and positive_int(observation["pid"]) and observation["pid"]==pid
            and (positive_int(observation.get("unified_exec_session"))
                 or (type(observation.get("unified_exec_chunk_id")) is str
                     and bool(observation["unified_exec_chunk_id"].strip())))
            and observation["result_sha256"]==result_sha and observation["terminal_sha256"]==terminal_sha,
            "externally observed exact OS terminal and digest chain required")


def mechanics_gate(review,current,geometry_sha,body_sha,active_sha,hessian_sha):
    require(review["status"]=="PASS_INDEPENDENT_KL_LOADING_ALL_OPERANDS"
            and review["pilot_id"]==PILOT_ID and review["bindings"]==current
            and review["geometry_result_sha256"]==geometry_sha and review["geometry_body_sha256"]==body_sha
            and review["active_result_sha256"]==active_sha and review["hessian_result_sha256"]==hessian_sha
            and yes(MECHANICS_CHECKS,review), "independent actual review of all same-envelope operands required")


def operand_block(doc,status,current,plan):
    require(doc["status"]==status and doc["pilot_id"]==PILOT_ID and doc["bindings"]==current
            and doc["R_binary_sha256"]==plan["K_binding"]["R_binary_sha256"]
            and doc["point_cache_sha256"]==plan["K_binding"]["point_cache_sha256"]
            and exact(doc["nominal_drive_exact"])==exact(plan["nominal_drive_exact"])
            and exact(doc["energy_scale_binary_exact"])==exact(plan["energy_scale_binary_exact"]) and
            exact(doc["fiber_radius_exact"])==exact(plan["fiber_radius_exact"]),
            "same exact slab, central operator, R and normalized energy operands required")


def read_provenance(block,result_sha):
    terminal_sha = block["terminal_sha256"]
    terminal = hashed_json(canonical(block["terminal"]),terminal_sha)
    require(terminal["status"]=="KL_LOADING_OPERAND_STAGE_FINISHED"
            and terminal["result_sha256"]==result_sha, "successful terminal for exact operand result")
    observed = hashed_json(canonical(block["os_exit"]),block["os_exit_sha256"])
    observed_zero("PASS_ROOT_OBSERVED_KL_LOADING_OPERAND_EXIT",observed,result_sha,terminal_sha,terminal["pid"])


def run_certificate(review_sha):
    current, plan, _, _ = prepared()
    review = hashed_json(HERE/"independent_mechanics_review.json",review_sha)
    gp = HERE/"actual_geometry"
    terminal_sha = review["geometry_terminal_sha256"]
    terminal = hashed_json(gp/"terminal.json",terminal_sha)
    gsha, bsha = terminal["result_sha256"], review["geometry_body_sha256"]
    observed = hashed_json(gp/"os_exit_observation.json",review["geometry_os_exit_sha256"])
    observed_zero("PASS_ROOT_OBSERVED_KL_LOADING_GEOMETRY_EXIT",observed,gsha,terminal_sha,terminal["pid"])
    active_link, hessian_link = review["active"],review["hessian"]
    ah,hh = active_link["result_sha256"],hessian_link["result_sha256"]
    mechanics_gate(review,current,gsha,bsha,ah,hh)
    geom = hashed_json(gp/"result.json",gsha)
    body = hashed_json(gp/"body.json",bsha)
    require(terminal["status"]=="KL_GEOMETRY_STAGE_FINISHED"
            and terminal["bindings"]==geom["bindings"]==current
            and positive_int(geom["pid"]) and geom["pid"]==terminal["pid"]
            and geom["body_sha256"]==bsha and geom["status"]=="PASS_KL_FRESH_WHOLE_ENVELOPE_GEOMETRY_CANDIDATE"
            and body["status"]=="CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING"
            and len(geom["bulk"])==144 and len(geom["traces"])==24,
            "complete independently admitted whole-envelope geometry required")
    for block,sha in ((active_link,ah),(hessian_link,hh)):
        read_provenance(block,sha)
    active = hashed_json(canonical(active_link["result"]),ah)
    hessian = hashed_json(canonical(hessian_link["result"]),hh)
    operand_block(active,"PASS_INDEPENDENT_KL_EXACT_ACTIVE_OPERATOR_BOUNDS",current,plan)
    operand_block(hessian,"PASS_INDEPENDENT_KL_FRESH_WHOLE_ENVELOPE_H_BOUNDS",current,plan)
    r,h = exact(plan["fiber_radius_exact"]),exact(plan["dose_gap_exact"])
    norm = exact(plan["same_R_infinity_norm_exact"])
    cols = (plan["signed_K_R_F_I_point_rows_exact"], hessian["signed_R_H_Y_v_rows_exact"],
            active["signed_R_active_center_rows_exact"], active["signed_R_C_v_rows_exact"])
    require(all(isinstance(col,list) and len(col)==N for col in cols),"of every operand all 611 signed rows required")
    result = conclude(*cols,exact(hessian["Zc_upper_exact"]),exact(active["R_C_norm_upper_exact"]),norm,h,r)
    uk,ul,v = (plan[field] for field in ("K_center_displacement_exact","L_center_displacement_exact","predictor_slope_exact"))
    end = [exact(a)+h*exact(b) for a,b in zip(uk,v,strict=True)]
    ea,eb = plan["anchor_posterior_error_upper_exact"]
    margins = (contains_root(uk,ea,uk,r), contains_root(ul,eb,end,r))
    linked = result["strict_fiber_selfmaps"] and min(margins)>0
    status = "CANDIDATE_CERTIFIED_ADJACENT_KL_LOADING_SEGMENT" if linked else "UNRESOLVED_KL_LOADING_SEGMENT"
    save(HERE/"certificate_candidate.json",{**result,"status":status,"pilot_id":PILOT_ID,"bindings":current,
         "endpoint_root_inclusion_margin_exact":margins,"independent_mechanics_review_sha256":review_sha,
         "linked_to_both_anchors":linked,"finite_element_equilibria_only":True,"full_loading_path_certified":False,
         "requires_external_readback_for_promotion":True})


def fixtures():
    check = {}

    def failure(fn):
        try:
            fn()
        except (ValueError,KeyError):
            return True
        return False

    ref = [(0,0),(1,1)]
    center, radius, k, l, v = secant_box(ref,ref,[(0,Q(1,5)),(1,1)],[1],Q(1,16),Q(1,3))
    check["swept_envelope_covers_both_fibers_and_preserves_fixed_axis"] = (
        center==[[0,Q(1,10)],[1,1]] and radius==[[0,Q(13,80)],[0,0]]
        and k==[0] and l==[Q(1,5)] and v==[Q(3,5)])
    check["axis_motion_and_boolean_radius_and_reversed_intervals_rejected"] = all((
        failure(lambda:secant_box(ref,ref,[(1,0),(1,1)],[1],Q(1,16),Q(1,3))),
        failure(lambda:secant_box(ref,ref,ref,[1],True,Q(1,3))),failure(lambda:pair([1,0]))))
    check["quadratic_range_uses_signed_interior_turning_point"] = (
        quad_bounds([0,0],[-2,-2],[2,2],Q(0),Q(1))==(Q(-1,2),0))
    # F(u,t)=(1+t)u-t, p(t)=t; R=1, so q=0, H(Y,0)v=1, f_I(0)=-1, and factor k=1.
    def linear_trial(r,g):
        return conclude([[0,0]],[[1,1]],[[-1,-1]],[[1,1]],Q(0),Q(1),Q(1),g,r)
    good=linear_trial(Q(1,8),Q(1,4))
    check["signed_affine_cancellation_yields_real_tube_and_posterior"] = (
        good["strict_fiber_selfmaps"] is True and good["posterior_error_upper_exact"]=="1/12"
        and good["whole_product_normalized_H_euclidean_lower_exact"]=="3/4")
    check["strict_selfmap_and_z_equal_one_boundaries_remain_unresolved"] = (
        linear_trial(Q(1,12),Q(1,4))["strict_fiber_selfmaps"] is False
        and linear_trial(Q(1,8),Q(1))["strict_fiber_selfmaps"] is False)
    check["posterior_inclusion_test_distinguishes_overlap_and_strict_touch"] = (
        contains_root([0],Q(1,2),[Q(3,4)],Q(1,2))<0
        and contains_root([0],Q(1,4),[Q(1,4)],Q(1,2))==0
        and contains_root([0],Q(1,4),[Q(1,8)],Q(1,2))>0)
    pre = {k+"_sha256":"0"*64 for k in ("source","protocol","preflight","plan","envelope")}
    p = {"status":"PASS_INDEPENDENT_KL_LOADING_GEOMETRY_PRELAUNCH","pilot_id":PILOT_ID,"bindings":pre,
         "checks":dict.fromkeys(GEOMETRY_CHECKS,True),"geometry_only_authorized":True,"mechanical_kernels_authorized":False}
    gate_geometry(p,pre)
    check["positive_external_geometry_phase_roundtrip_and_foreign_envelope_rejection"] = failure(
        lambda:gate_geometry(p,{**pre,"envelope_sha256":"f"*64}))
    os_doc = {"status":"PASS_ROOT_OBSERVED_KL_LOADING_GEOMETRY_EXIT","pid":42,"observed_exit_code":0,
              "unified_exec_chunk_id":"fake-fixture-only","result_sha256":"b","terminal_sha256":"t"}
    observed_zero(os_doc["status"],os_doc,"b","t",42)
    check["independent_actual_gate_rejects_float_OS_code_and_boolean_PID"] = all((
        failure(lambda:observed_zero(os_doc["status"],{**os_doc,"observed_exit_code":0.0},"b","t",42)),
        failure(lambda:observed_zero(os_doc["status"],{**os_doc,"pid":True},"b","t",True))))
    m = {"status":"PASS_INDEPENDENT_KL_LOADING_ALL_OPERANDS","pilot_id":PILOT_ID,"bindings":pre,
         "geometry_result_sha256":"g","geometry_body_sha256":"b","active_result_sha256":"a",
         "hessian_result_sha256":"h","checks":dict.fromkeys(MECHANICS_CHECKS,True)}
    mechanics_gate(m,pre,"g","b","a","h")
    check["absent_fresh_full_envelope_H_and_foreign_active_result_rejected"] = all((
        failure(lambda:mechanics_gate({**m,"checks":{}},pre,"g","b","a","h")),
        failure(lambda:mechanics_gate(m,pre,"g","b","foreign","h"))))
    require(all(check.values()),f"failed cold tests: {check}")
    return check


def prepare():
    tests = fixtures()
    plan, envelope, old_hashes = cold_plan()
    decoded_envelope = json.loads(json.dumps(envelope,default=str))
    with patch.dict(geometry_bounds.__globals__, {"domain_components": lambda _,ctx,p,b:
            ctx.drive==exact(p["nominal_drive_exact"]) and b==decoded_envelope}):
        require(geometry_bounds(plan,decoded_envelope) is True,"real K context cold-route fixture failed")
    tests["real_K_candidate_context_and_frozen_box_route_without_geometry"] = True
    inputs = {"source_sha256":digest(Path(__file__)),"plan_sha256":save(HERE/"plan.json",plan),
              "envelope_sha256":save(HERE/"envelope.json",envelope)}
    protocol_sha = save(HERE/"protocol.json",{"status":"COLD_KL_LOADING_TUBE_PILOT_PROTOCOL",
        "pilot_id":PILOT_ID,"bindings":inputs,"individual_authority_sha256":INPUTS,
        "read_only_method_source_sha256":{name:digest(canonical(name)) for name in METHOD_SOURCES},
        "geometry_source_tree_sha256":local_source_tree(),"geometry_source_review_path":"independent_geometry_prelaunch.json",
        "geometry_source_required_checks":GEOMETRY_CHECKS,"mechanical_release_required_checks":MECHANICS_CHECKS,
        "geometry_or_mechanics_launch_authorized":False,"pilot_trial_radius_normalized_exact":RADIUS,
        "whole_KL_path_certified":False,**NO_NUMERICS})
    save(HERE/"preflight.json",{"status":"PASS_KL_LOADING_COLD_FIXTURES","bindings":{**inputs,"protocol_sha256":protocol_sha},
         "checks":tests,"frozen_design_hashes_rechecked":old_hashes,"native_K_and_L_centers_and_R_premise_matched":True,
         "all_signed_K_center_rows_matched":N,"whole_KL_path_certified":False,**NO_NUMERICS})


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    commands=cli.add_subparsers(dest="command",required=True)
    for name in ("prepare","verify-preparation"):
        commands.add_parser(name)
    for name in ("geometry","certificate"):
        commands.add_parser(name).add_argument("--review-sha",required=True)
    args=cli.parse_args()
    if args.command=="prepare":
        prepare()
    elif args.command=="verify-preparation":
        current,_,_,_=prepared()
        print(json.dumps({"bindings":current,"checks":fixtures()},indent=2))
    elif args.command=="geometry":
        run_geometry(args.review_sha)
    else:
        run_certificate(args.review_sha)


if __name__ == "__main__":
    main()
