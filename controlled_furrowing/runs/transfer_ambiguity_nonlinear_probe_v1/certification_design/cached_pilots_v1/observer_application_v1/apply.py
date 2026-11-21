"""Apply reviewed cached-quadrature nodal errors to the fixed observers."""
import argparse
import base64
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import io
from itertools import product
import json
from pathlib import Path
import sys

import numpy as np

SOURCE = Path(__file__).resolve()
OUT = SOURCE.parent
DESIGN = OUT.parents[1]
RUN = DESIGN.parent
ROOT = RUN.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from nodal_observer_error_transfer import endpoint_bounds
from quadratic_endpoint_observer import _state_controls
from quadratic_mean_derivatives import _boundary_node_ids

PILOTS = tuple(f"{e}_{c}" for e,c in product(
    ("max_contrast_change","max_negative_contrast_change"), ("small_17","large_03")))
OBS = ("maximum_shortening_um","region_averaged_shortening_um",
       "apical_max_shortening_um","basal_max_shortening_um")
VARIANTS = tuple(f"{w}__{a}" for w,a in product(
    ("current_disk","current_radial_line"), ("spatial_max","regional_mean")))
CACHE = DESIGN / "cached_pilots_v1/independent_actual_preparation/actual_replay"
GEOM = DESIGN / "exact_targets_v1/independent_actual/replay"
TARGETS = DESIGN / "exact_targets_v1/preparation/target_manifest.json"
BUD = DESIGN / "budgets_v1/actual/results.json"
BUD_REVIEW = DESIGN / "budgets_v1/independent_actual/receipt.json"
PRIMARY = RUN / "independent_actual_preparation/actual_replay/receipt.json"
SWEEP = RUN / "observer_robustness/actual/results.json"
SWEEP_REVIEW = RUN / "observer_robustness/independent_actual/receipt.json"
FIELDS = RUN / "field_manifest.json"
PAIR = RUN / "preparation/protocol.json"
TRANSFER = ROOT / "evidence/nodal_observer_error_independent_review.json"
LENGTH = F(67, 2)
RADIUS = F(1, 2**36)
PINS = {
    CACHE/"receipt.json": "4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17",
    CACHE/"os_exit_observation.json": "65858f55c1281d3953b775f147c7475a8604891f37cf56ed0fdd47697e6df6e3",
    GEOM/"receipt.json": "884e0bc76a27465c86b62a202cab491d764aedf030623592238e09adee86d77d",
    GEOM/"os_exit_observation.json": "c2be4be8b403e170e20ec54c34f17a69b4aec4299e6f9f1ab7ecb3d4afa337d9",
    TARGETS: "cb25fbeaef753cea700c187cd8b377b7dcf4aa1fe2b0b246ed17134a07a89bc3",
    BUD: "689948766db950ee4ddf399e6731aea08890cc7a1f84c5615d19e7ded68b0872",
    BUD_REVIEW: "978a2a4095bb8735103a5ee2c11ad266a6b0e052e22af2b7a36ef29f1f38265c",
    PRIMARY: "c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd",
    SWEEP: "53a17296ddae08b4d1fe6b7e9b20ad0a377bf8dc962142bbfb9eb9877a381647",
    SWEEP_REVIEW: "703b574bc51a083b10c7328d9622862e939ddabe791b18ce9e1817f07718f69f",
    FIELDS: "308ca4c05e476ba6a14cbc7b8d34018983b8f11263e2c224488672b35eb38ac0",
    PAIR: "9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94",
    TRANSFER: "335626b91c72ac5ec37b1b1dd95b4825d3dc847df02d9cf993e5ed6daff5f35a",
    ROOT/"scripts/quadratic_endpoint_observer.py": "075dac2148ed26a3bfe83fadf1b97602e4b08536f2be38cb95c2054776c0b926",
    ROOT/"scripts/quadratic_mean_derivatives.py": "681b69335a80ae0468dd98856c5562b24524827809ece465747b95454d152bff",
}


def sha(value):
    return sha256(value).hexdigest()


def file_sha(path):
    return sha(path.read_bytes())


def put(path, value):
    with path.open("x",encoding="utf8",newline="\n") as f:
        json.dump(value,f,indent=2,default=str,allow_nan=False)
        f.write("\n")
    return file_sha(path)


def bound(path, digest, links, parse=True):
    path = (ROOT/path).resolve()
    rel = path.relative_to(ROOT).as_posix()
    data = path.read_bytes()
    assert sha(data) == digest and links.setdefault(rel,digest) == digest, rel
    return json.loads(data) if parse else data


def nonempty_true(values):
    return bool(values) and all(v is True for v in values.values())


def exact_pair(value):
    low,high = map(F,value)
    assert low <= high
    return low,high


def rat_json(v):
    return F(v["numerator"],v["denominator"])


def inflate(interval, e):
    assert e >= 0
    return interval[0]-e, interval[1]+e


def delta_class(interval, measured, delta):
    assert delta >= 0
    lo,hi = interval[0]-measured, interval[1]-measured
    mn = max(F(0),lo,-hi)
    mx = max(abs(lo),abs(hi))
    return ("IN" if mx <= delta else "OUT" if mn > delta else "UNRESOLVED"), (mn,mx)


def posterior(e, z0, v, r, claimed=None):
    z = z0+v
    assert e >= 0 and z0 >= 0 and v >= 0 and 0 <= z < 1 and r > 0
    assert e+z*r < r
    error = e/(1-z)
    assert error < r and (claimed is None or error == claimed)
    return error


def verify_exit(base, docs, links, status):
    exit_doc = docs[base/"os_exit_observation.json"]
    term = bound(base/"terminal.json",exit_doc["terminal_sha256"],links)
    h = PINS[base/"receipt.json"]
    assert type(term["exit_code"]) is int and term["exit_code"] == 0
    assert term["pid"] == exit_doc["pid"] and term["receipt_sha256"] == exit_doc["receipt_sha256"] == h
    assert term["status"] == "FINISHED" and exit_doc["status"] == status
    assert type(exit_doc["observed_exit_code"]) is int and exit_doc["observed_exit_code"] == 0


def authority():
    links = {}
    docs = {x:bound(x,h,links,x.suffix==".json") for x,h in PINS.items()}
    cache = docs[CACHE/"receipt.json"]
    geom = docs[GEOM/"receipt.json"]
    assert cache["status"] == "PASS_INDEPENDENT_FOUR_CACHED_LOCAL_ROOTS"
    assert cache["root_certificates_cached_quadrature_only"] == 4
    assert cache["root_certificates_exact_spatial_integral"] == geom["exact_spatial_root_count"] == 0
    assert geom["status"] == "PASS_INDEPENDENT_ACTUAL_EIGHT_SINGLETON_GEOMETRY_BOXES"
    assert geom["force_hessian_or_AD_jobs"] == geom["spatial_derivative_or_tail_bound_count"] == 0
    assert nonempty_true(cache["checks"]) and nonempty_true(geom["checks"])
    verify_exit(CACHE,docs,links,"PASS_OBSERVED_SERIAL_INDEPENDENT_REPLAY_TRUE_OS_EXIT")
    verify_exit(GEOM,docs,links,"PASS_OBSERVED_INDEPENDENT_GEOMETRY_REPLAY_TRUE_OS_EXIT")
    assert tuple(docs[TARGETS]["pilot_order"]) == PILOTS
    assert {f["pilot_id"] for f in cache["field_results"]} == set(PILOTS)
    assert {(x["pilot_id"],F(x["radius_exact"])) for x in geom["geometry_records"]} == {
        (p,r) for p,r in product(PILOTS,(RADIUS,4*RADIUS))}
    assert len(geom["geometry_records"]) == 8
    t = docs[TRANSFER]
    assert t["status"] == "PASS_CONDITIONAL_EXACT_NODAL_TO_OBSERVER_TRANSFER" and nonempty_true(t["checks"])
    for name,digest in t["reviewed_hashes"].items():
        bound(Path(name),digest,links,False)
    b,bi = docs[BUD],docs[BUD_REVIEW]
    assert len(b["field_budgets"]) == len(bi["fields"]) == 36 and bi["root_error_bounds_proved"] == 0
    assert nonempty_true(bi["all_literal_checks"])
    assert len(b["four_pilot_selection"]) == len(bi["four_pilot_fields"]) == 4
    sweep,si = docs[SWEEP],docs[SWEEP_REVIEW]
    assert si["status"] == "PASS_INDEPENDENT_SECONDARY_STORED_OBSERVER_SWEEP" and nonempty_true(si["checks"])
    assert si["input_and_output_sha256"][SWEEP.relative_to(ROOT).as_posix()] == PINS[SWEEP]
    assert sweep["small_constraint_evaluations"] == si["independent_interval_exact_matches"] == 544
    assert sweep["primary_replay_exactly_matches"] is si["primary_exact_result_kept"] is True
    assert docs[PRIMARY]["status"] == "PASS_INDEPENDENT_ACTUAL_NONLINEAR_STORED_FIELD_REPLAY"
    assert docs[FIELDS]["independent_primary_review_sha256"] == PINS[PRIMARY]
    delta = F(sweep["frozen_primary_source_allowance_um_exact"])
    assert delta == F(si["source_configured_delta_um_exact"]) == F(b["conditional_source_allowance_um_exact"])
    assert delta == F(docs[PAIR]["selection"]["conditional_source_allowance_um_exact"])
    return docs,links,delta


def plan(links, delta):
    return {
        "status": "FROZEN_FOUR_CACHED_ROOT_OBSERVER_APPLICATION",
        "source_sha256": file_sha(SOURCE),
        "inputs_sha256": links, "pilot_order": PILOTS, "box_radius_exact": RADIUS,
        "frozen_allowance_um_exact": delta, "length_scale_um_exact": LENGTH,
        "scope": "Apply prior independently certified finite cached-quadrature root errors on four exact fixed inputs. Bind coherent exact narrow root boxes, then apply the reviewed P2 graph theorem to the retained observer intervals. No AD, equilibrium solve, continuum or spatial-integral root claim."
    }


def prep():
    _,links,delta = authority()
    assert posterior(F(1,8),F(1,4),F(1,4),F(1)) == F(1,4)
    try:
        posterior(F(1,2),F(1,4),F(1,4),F(1))
    except AssertionError:
        pass
    else:
        raise AssertionError("equality must fail strict self-map")
    f = F(1,16)
    assert delta_class((-2,-1),0,2)[0] == "IN"
    assert delta_class((-3,-2-f),0,2)[0] == "OUT"
    assert delta_class((-3,-2),0,2)[0] == "UNRESOLVED"
    assert inflate((1,2),f) == (1-f,2+f)
    # Exact, synthetic steep-line example: Ea=5/4*(1+3)*epsilon.
    top = [[(F(t),3*F(t)+10) for t in (0,F(1,2),1)]]
    bot = [[(F(t),0) for t in (0,F(1,2),1)]]
    tr = endpoint_bounds(top,bot,f,fixed_radial_endpoints=True,roi_um=(0,1))
    assert tr["status"] == "CONDITIONAL_EXACT_ERROR_TRANSFER"
    assert F(tr["apical"]["Eulerian_vertical_error_um_exact"]) == 5*f
    assert F(tr["endpoint_error_um_exact"]["maximum_shortening_um"]) == F(25,4)*f
    assert not (OUT/"results.json").exists()
    h = put(OUT/"protocol.json",plan(links,delta))
    q = put(OUT/"preflight.json",{"status":"PASS_PREPARATION_ONLY","source_sha256":file_sha(SOURCE),
        "protocol_sha256":h,"checks":{"pinned_parents_and_both_review_OS_exits":True,
        "posterior_exact_and_strict_boundary":True,"closed_interval_classification_and_inflation":True,
        "reviewed_transfer_synthetic_slope":True,"no_actual_endpoint_arithmetic_during_preflight":True}})
    print(json.dumps({"source_sha256":file_sha(SOURCE),"protocol_sha256":h,"preflight_sha256":q}))


def select(rows, **pairs):
    row, = [r for r in rows if all(r[k] == v for k,v in pairs.items())]
    return row


def controls_and_gap(state, packet, proof, contract, box, budget):
    binf = lambda x: F.from_float(float(x))
    blob = base64.b64decode(packet["frozen_binary_cache"]["data"],validate=True)
    assert sha(blob) == packet["frozen_binary_cache"]["sha256"] == proof["cache_array_contents_sha256"]
    assert sha(blob) == contract["frozen_binary_cache_sha256"]
    with np.load(io.BytesIO(blob),allow_pickle=False) as a:
        meshes = a["reference"].tolist(), a["triangles"].tolist()
        free = list(map(int,a["free"]))
        u = a["u"].tolist()
        assert a["preconditioner"].shape == (611,611)
        assert sha(a["preconditioner"].tobytes(order="C")) == contract["binary_preconditioner_sha256"]
        trace = {"bottom":a["surface_0_edges"].tolist(),"top":a["surface_1_edges"].tolist()}
    ref,tris = meshes
    assert ref == state["reference_nodes"] and tris == state["triangles"] == box["triangles"]
    assert u == state["free_displacement"] and free == box["free"] == budget["center_and_model"]["free_map"]
    assert len(ref) == len(state["current_nodes"]) == 325 and len(tris) == 144
    u = list(map(binf,u))
    X = [binf(c) for row in ref for c in row]
    stored = [binf(c) for row in state["current_nodes"] for c in row]
    rim = binf(state["parameters"]["outer_radius"])
    fixed = {i for i in range(650) if X[2*(i//2)] == rim or (i%2 == 0 and X[i] == 0)}
    assert free == [i for i in range(650) if i not in fixed] and len(free) == 611
    Y = X.copy()
    for i,v in zip(free,u,strict=True):
        Y[i] += v
    assert list(map(binf,Y)) == stored
    assert all(Y[i] == X[i] == stored[i] for i in fixed)
    assert X == [F(x) for v in box["reference_nodes_exact"] for x in v]
    assert Y == [F(x) for v in box["current_centers_exact"] for x in v]
    radii = [F(0) if i in fixed else RADIUS for i in range(650)]
    assert radii == [F(x) for v in box["coordinate_radii_exact"] for x in v]
    assert [(a-r,a+r) for a,r in zip(Y,radii,strict=True)] == [exact_pair(x) for row in box["coordinate_hulls_exact"] for x in row]
    assert len(box["predictor_hulls_exact"]) == 611 and all(exact_pair(x) == (0,0) for x in box["predictor_hulls_exact"])
    assert box["root_membership_assumed"] is False
    gaps = [abs(a-b) for a,b in zip(Y,stored,strict=True)]
    xi = max(gaps)
    assert xi == F(budget["center_and_model"]["node_addition_gap_dimensionless_exact"])
    ids = _boundary_node_ids(state)
    controls,height = _state_controls(state,LENGTH)
    endpoints = {}
    assert ids == trace and height == binf(state["parameters"]["height"])*LENGTH
    for face in ("top","bottom"):
        ordered = ids[face]
        assert len(ordered) == 12 and all(a[-1] == b[0] for a,b in zip(ordered,ordered[1:]))
        pair = [2*ordered[0][0],2*ordered[-1][-1]]
        assert set(pair) <= fixed and [Y[j] for j in pair] == [0,rim]
        assert [[tuple(LENGTH * stored[2*n+j] for j in (0,1)) for n in e] for e in ordered] == controls[face]
        endpoints[face] = pair
    return controls,height,xi,endpoints


def apply_pilot(pid, docs, links, delta):
    cache,geom = docs[CACHE/"receipt.json"],docs[GEOM/"receipt.json"]
    cmeta = select(cache["field_results"],pilot_id=pid)
    field = select(docs[FIELDS]["fields"],field=cmeta["field"])
    budget = select(docs[BUD]["field_budgets"],field=field["field"])
    target = select(docs[TARGETS]["targets"],pilot_id=pid)
    contract = bound(Path(target["model_contract"]),target["model_contract_sha256"],links)
    assert contract["source_field"] == cmeta["field"] == field["field"] == budget["field"]
    assert contract["source_field_sha256"] == cmeta["field_sha256"] == field["field_sha256"] == budget["field_sha256"]
    assert contract["source_geometry_verification"] == field["verification"] and contract["source_geometry_verification_sha256"] == field["verification_sha256"]
    state = bound(Path(field["field"]),field["field_sha256"],links)
    ver = bound(Path(field["verification"]),field["verification_sha256"],links)
    assert docs[PRIMARY]["all_hashes"][field["field"]] == field["field_sha256"]
    assert docs[PRIMARY]["all_hashes"][field["verification"]] == field["verification_sha256"]
    proof = bound(Path(contract["cached_root_independent_receipt"]),cmeta["receipt_sha256"],links)
    assert contract["cached_root_independent_receipt_sha256"] == cmeta["receipt_sha256"]
    assert proof["status"] == cmeta["status"] == "PASS_INDEPENDENT_CACHED_LOCAL_ROOT"
    assert proof["field"] == field["field"] and proof["field_sha256"] == field["field_sha256"]
    assert proof["record_id"] == cmeta["record_id"] == budget["record_id"] == contract["physical_source_record_id"]
    assert proof["pilot_id"] == contract["pilot_id"] == pid and proof["site_counts"] == {"bulk":9216,"surface":288}
    assert proof["rounded_strength_used_in_proof"] is False and proof["npz_current_model_equal_all_except_R"] is True
    assert cache["producer_and_payload_sha256"][contract["point_cache"]] == contract["point_cache_sha256"] == proof["cache_packet_sha256"]
    packet = bound(Path(contract["point_cache"]),proof["cache_packet_sha256"],links)
    assert packet["input_state"] == field["field"] and packet["input_state_sha256"] == field["field_sha256"]
    assert packet["parameters"] == state["parameters"] == budget["center_and_model"]["material_parameters"]
    assert packet["material_parameter_storage_exact"] == contract["material_parameters_binary_exact"]
    assert packet["pilot_id"] == pid and packet["record_id"] == proof["record_id"]
    exact_drive = contract["exact_requested_drive"]
    assert all(packet[k] == v for k,v in exact_drive.items())
    dose,gain,drive,rounded = (F(exact_drive[k]) for k in ("optical_dose_exact","common_gain_exact","activation_strength_exact","activation_strength_binary_exact"))
    assert 0 < drive == dose*gain == F(proof["exact_strength"])
    assert drive != rounded == F.from_float(float(drive)) == F.from_float(state["activation_strength"])
    assert F(field["step"]["dose_exact"]) == dose and F(field["step"]["requested_strength_exact"]) == drive
    assert F(packet["normalization_binary_exact"]) == F(contract["energy_scale_binary_exact"]) > 0
    assert F(contract["cached_root_radius_exact"]) == F(proof["radius_exact"]) == RADIUS
    eta,z0,v,z = (F(proof[k]) for k in ("independent_eta_upper_exact","independent_Z0_upper_exact",
                                        "independent_variation_upper_exact","independent_z_upper_exact"))
    assert z == z0+v and eta+z*RADIUS == F(proof["independent_selfmap_lhs_upper_exact"])
    assert RADIUS-eta-z*RADIUS == F(proof["strict_selfmap_gap_lower_exact"])
    err = posterior(eta,z0,v,RADIUS,F(proof["independent_posterior_error_upper_exact"]))
    assert err == F(contract["cached_root_error_exact"])
    assert RADIUS < F(proof["science_radius_budget_exact"]) == F(budget["e_goal_dimensionless_exact"]) == F(contract["observer_design_goal_exact"])

    gm = select(geom["geometry_records"],pilot_id=pid,radius_exact=str(RADIUS))
    gg = bound(Path(gm["receipt"]),gm["sha256"],links)
    assert gg["status"] == "PASS_FRESH_SINGLETON_WHOLE_GEOMETRY_ONLY"
    assert gg["pilot_id"] == pid and gg["record_id"] == contract["record_id"]
    assert gg["proof_center_field_sha256"] == field["field_sha256"]
    assert F(gg["radius_exact"]) == F(gg["cached_proof_radius_exact"]) == RADIUS
    assert gg["all_144_body_and_24_trace_domains"] is True and gg["larger_cached_derivative_proof_reused"] is False
    assert F(gg["strong_monotonicity_lower_exact"]) > 0 and exact_pair(gg["body_J_bounds_exact"])[0] > 0
    assert gm["box_sha256"] == gg["box_sha256"] == target["box_0_sha256"]
    assert gm["body_geometry_sha256"] == gg["body_geometry_sha256"]
    geometry_dir = DESIGN/"exact_targets_v1/geometry"/pid/"box_0"
    body_name = geometry_dir/"body_geometry.json"
    assert geom["checked_file_hashes"][body_name.relative_to(ROOT).as_posix()] == gg["body_geometry_sha256"]
    body = bound(body_name,gg["body_geometry_sha256"],links)
    assert body["status"] == "CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING"
    box = bound(Path(target["box_0"]),target["box_0_sha256"],links)
    assert box["record_id"] == contract["record_id"]
    assert box["point_cache_sha256"] == proof["cache_packet_sha256"]
    assert box["model_observation_contract_sha256"] == target["model_contract_sha256"]
    assert box["parameter_cell_sha256"] == contract["parameter_cell_sha256"]
    assert box["predictor_sha256"] == target["predictor_sha256"]
    assert F(box["dose_exact"]) == dose and F(box["corrector_radius_exact"]) == RADIUS
    assert all(exact_pair(v)[0] == exact_pair(v)[1] for v in contract["parameter_cell"].values())
    assert exact_pair(contract["parameter_cell"]["common_gain"]) == (gain,gain)

    controls,height,xi,radial_dofs = controls_and_gap(state,packet,proof,contract,box,budget)
    roi = exact_pair(budget["fixed_roi_um_exact"])
    assert roi[0] == 0 and roi[1] > 0
    e_um = LENGTH*(err+xi)
    tr = endpoint_bounds(controls["top"],controls["bottom"],e_um,fixed_radial_endpoints=True,roi_um=roi)
    assert tr["status"] == "CONDITIONAL_EXACT_ERROR_TRANSFER"
    for side,face in (("apical","top"),("basal","bottom")):
        v = tr[side]
        f = budget["stored_face_derivatives"][face]
        assert F(v["radial_derivative_margin_exact"]) + 8*e_um == F(f["radial_min_um_exact"])
        assert F(v["base_graph_Lipschitz_upper_exact"]) == F(f["global_stored_slope_cap_exact"])
        assert F(v["Eulerian_vertical_error_um_exact"]) == F(f["C_um_per_dimensionless_coordinate_exact"])*(err+xi)
    Ea,Eb = (F(tr[k]["Eulerian_vertical_error_um_exact"]) for k in ("apical","basal"))
    assert [F(tr["endpoint_error_um_exact"][k]) for k in OBS] == [Ea+Eb,Ea+Eb,Ea,Eb]
    result = {"pilot_id":pid,"endpoint_id":field["endpoint_id"],"role":field["role"],"record_id":proof["record_id"],
              "root_operator":"exact-real frozen binary-cache quadrature",
              "cached_root_receipt_sha256":cmeta["receipt_sha256"],"geometry_receipt_sha256":gm["sha256"],
              "box_sha256":gm["box_sha256"],
              "eta_exact":eta,"z_exact":z,"certified_radius_exact":RADIUS,"posterior_radius_exact":err,
              "serialized_X_plus_u_gap_exact":xi,"physical_nodal_error_um_exact":e_um,
              "graph_transfer":tr,"fixed_radial_endpoint_dofs":radial_dofs,
              "same_cached_root_cube_is_inside_independently_embedded_geometry":True,
              "ROI_um_exact":roi,"exact_dose":dose,"exact_gain":gain,"exact_drive":drive}

    if field["role"] == "large_final":
        endpoint = select(docs[PRIMARY]["endpoints"],endpoint_id=field["endpoint_id"])
        stored = exact_pair(budget["stored_Q_um_exact"])
        assert stored == exact_pair(endpoint["large_Q_um_exact"]) == exact_pair(ver["surface"]["Q_um_exact"])
        assert F(ver["surface"]["shoulder_radius_um_exact"]) == roi[1] == F(docs[PAIR]["large_observable"]["shoulder_radius_um_exact"])
        # The fixed axis centre is a control node with error e_um, not 5/4 e_um.
        bound_q = e_um+Ea
        assert bound_q == F(budget["conditional_Q_error_coefficient_um_exact"])*(err+xi)
        result.update({"stored_Q_um_exact":stored,"Q_error_upper_um_exact":bound_q,
                       "cached_root_Q_um_exact":inflate(stored,bound_q)})
        return result

    assert field["role"] == "small_record" and proof["record_id"] == "optoGEF2_small_col_L"
    sweep,si = docs[SWEEP],docs[SWEEP_REVIEW]
    sf = select(sweep["field_observer_receipts"],
                path=(RUN/"observer_robustness/actual"/field["endpoint_id"]/"small_17.json").relative_to(ROOT).as_posix())
    assert si["input_and_output_sha256"][sf["path"]] == sf["sha256"]
    obs = bound(Path(sf["path"]),sf["sha256"],links)
    for key in ("endpoint_id","checkpoint","role","step","field","field_sha256","verification","verification_sha256"):
        assert obs[key] == field[key]
    assert obs["record_id"] == proof["record_id"]
    oc = obs["observer_contract"]
    assert tuple(map(rat_json,oc["current_radius_ROI_um"])) == roi
    assert rat_json(oc["reference_apical_um"]) == height and rat_json(oc["reference_basal_um"]) == 0
    assert roi[1] == F.from_float(state["parameters"]["patch_radius"])*LENGTH
    rows = []
    for variant in VARIANTS:
        source = sweep["by_endpoint_and_variant"][field["endpoint_id"]][variant]["all68_residuals"]["constraint_rows"]
        for name,error in zip(OBS,(Ea+Eb,Ea+Eb,Ea,Eb),strict=True):
            row = select(source,record_id=proof["record_id"],observable=name)
            saved = exact_pair(row["prediction_interval_exact"])
            assert saved == exact_pair(obs["variant_observables_um_exact"][variant][name])
            y = F(row["observed_exact"])
            assert y == F(field["step"]["observed_exact"][name]) and F(row["discrepancy_allowance_exact"]) == delta
            enclosed = inflate(saved,error)
            status,residual = delta_class(enclosed,y,delta)
            rows.append({"variant":variant,"observable":name,"recorded_value_um_exact":y,
                         "saved_um_exact":saved,"observer_error_upper_um_exact":error,"cached_root_um_exact":enclosed,
                         "absolute_residual_um_exact":residual,"status_at_frozen_delta":status})
    assert len(rows) == 16
    result["root_observer_comparisons"] = rows
    return result


def run(protocol_sha, preflight_sha):
    docs,links,delta = authority()
    bindings = {}
    p = bound(OUT/"protocol.json",protocol_sha,bindings)
    pf = bound(OUT/"preflight.json",preflight_sha,bindings)
    assert p == json.loads(json.dumps(plan(links,delta),default=str))
    assert pf["status"] == "PASS_PREPARATION_ONLY" and nonempty_true(pf["checks"])
    assert pf["source_sha256"] == p["source_sha256"] == file_sha(SOURCE)
    assert pf["protocol_sha256"] == protocol_sha
    cases = [apply_pilot(pid,docs,links,delta) for pid in PILOTS]
    qplus = cases[1]["cached_root_Q_um_exact"]
    qminus = cases[3]["cached_root_Q_um_exact"]
    gap = qplus[0]-qminus[1],qplus[1]-qminus[0]
    splus = cases[1]["stored_Q_um_exact"]
    sminus = cases[3]["stored_Q_um_exact"]
    stored_gap = exact_pair(docs[PRIMARY]["stored_large_Q_signed_difference_um_exact"])
    assert stored_gap == (splus[0]-sminus[1],splus[1]-sminus[0])
    assert gap == inflate(stored_gap,cases[1]["Q_error_upper_um_exact"]+cases[3]["Q_error_upper_um_exact"])
    assert cases[1]["exact_dose"] == cases[3]["exact_dose"] == F(docs[PAIR]["matched_input"]["matched_dose_exact"])
    sr = [r for c in (cases[0],cases[2]) for r in c["root_observer_comparisons"]]
    assert len(sr) == 32
    out = {
        "status":"CANDIDATE_CACHED_QUADRATURE_ROOT_OBSERVER_APPLICATION_AWAITING_REVIEW",
        "bindings":bindings,"source_sha256":file_sha(SOURCE),"inputs_sha256":links,
        "four_pilots":cases,"frozen_model_discrepancy_allowance_um_exact":delta,
        "stored_large_signed_gap_um_exact":stored_gap,"cached_root_large_signed_gap_um_exact":gap,
        "strict_positive_cached_root_gap":gap[0]>0,"L_static_observer_comparisons":len(sr),
        "comparison_counts":dict(Counter(r["status_at_frozen_delta"] for r in sr)),
        "total_cached_root_small_states":2,"full_68_per_endpoint_claimed":False,
        "exact_spatial_integral_root_observer_claimed":False,"PDE_or_source_fit_claimed":False,
        "no_new_force_or_AD":True,"no_new_equilibrium_solves":True,
        "scope":p["scope"],
    }
    h = put(OUT/"results.json",out)
    print(json.dumps({"results_sha256":h,"gap_limits_um":[float(x) for x in gap],
                      "comparison_counts":out["comparison_counts"]},default=str))


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    stages = cli.add_subparsers(dest="stage",required=True)
    stages.add_parser("prepare")
    actual = stages.add_parser("run")
    actual.add_argument("--protocol-sha",required=True)
    actual.add_argument("--preflight-sha",required=True)
    args = cli.parse_args()
    if args.stage == "prepare":
        prep()
    else:
        run(args.protocol_sha,args.preflight_sha)
