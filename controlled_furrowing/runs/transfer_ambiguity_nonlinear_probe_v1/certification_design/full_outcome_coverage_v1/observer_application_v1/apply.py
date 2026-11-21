"""Apply independently certified cached-root error to frozen primary observers."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
COVER = HERE.parent
DESIGN = COVER.parent
RUN = DESIGN.parent
ROOT = next(p for p in RUN.parents if (p / "PLAN.md").is_file())
FOUR = DESIGN / "cached_pilots_v1/independent_actual_preparation/actual_replay"
M = COVER / "remaining_cached_wave1/VERIFIED-CACHED-COVERAGE.json"
WAVES = COVER / "remaining_cached_30_v1/independent_actual_preparation/all_waves_readback.json"
RECENTER = COVER / "remaining_cached_30_v1/recenter_v1/independent_actual_preparation/actual_replay"
BUDGETS = DESIGN / "budgets_v1/actual/results.json"
BUDGET_REVIEW = DESIGN / "budgets_v1/independent_actual/receipt.json"
PRIMARY = RUN / "independent_actual_preparation/actual_replay/receipt.json"
FIELDS = RUN / "field_manifest.json"
PAIR = RUN / "preparation/protocol.json"
TRANSFER = ROOT / "evidence/nodal_observer_error_independent_review.json"
TRANSFER_SOURCE = ROOT / "scripts/nodal_observer_error_transfer.py"
PINS = {
    FOUR / "receipt.json": "4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17",
    FOUR / "os_exit_observation.json": "65858f55c1281d3953b775f147c7475a8604891f37cf56ed0fdd47697e6df6e3",
    M: "d5591fe4cfa66ef6e1d25eeaafde6df1000352c94225dd680584dab3018901b1",
    WAVES: "e9c8587930e6cf86c3317403657290d3a24d92cb5a46be53fc0a8472b86c9fc2",
    COVER / "receipt.json": "235e8336c735b1bf1a75889018ae6a18586fd0905a65d38cb683c5ed236df19a",
    BUDGETS: "689948766db950ee4ddf399e6731aea08890cc7a1f84c5615d19e7ded68b0872",
    BUDGET_REVIEW: "978a2a4095bb8735103a5ee2c11ad266a6b0e052e22af2b7a36ef29f1f38265c",
    PRIMARY: "c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd",
    FIELDS: "308ca4c05e476ba6a14cbc7b8d34018983b8f11263e2c224488672b35eb38ac0",
    PAIR: "9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94",
    TRANSFER: "335626b91c72ac5ec37b1b1dd95b4825d3dc847df02d9cf993e5ed6daff5f35a",
    TRANSFER_SOURCE: "001769f23ba386fe6bde4530c12e011f72f8733bb1a4a2923182796e54da35de",
    ROOT / "scripts/annulus_observer.py": "77bcd1bc10e557de4a4be9c304041f9f0dfe0fcbac47871cabfbb11ffb7346b9",
    ROOT / "scripts/quadratic_endpoint_observer.py": "075dac2148ed26a3bfe83fadf1b97602e4b08536f2be38cb95c2054776c0b926",
    ROOT / "scripts/quadratic_mean_derivatives.py": "681b69335a80ae0468dd98856c5562b24524827809ece465747b95454d152bff",
}
ENDPOINTS = ("max_contrast_change", "max_negative_contrast_change")
CALIBRATION_CHECKPOINTS = tuple(f"small_{i:02d}" for i in range(1,18))
RECENTER_CASES = {"max_contrast_change_small_13", "max_negative_contrast_change_small_16"}
OBS = ("maximum_shortening_um", "region_averaged_shortening_um",
       "apical_max_shortening_um", "basal_max_shortening_um")
PRIMARY_VARIANT = "current_disk__spatial_max"
LENGTH = Q(67, 2)


def demand(ok, message):
    if not ok:
        raise ValueError(message)


def yes(checks):
    return type(checks) is dict and bool(checks) and all(v is True for v in checks.values())


def integer_zero(value):
    return type(value) is int and value == 0


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def rel(path):
    return path.resolve().relative_to(ROOT).as_posix()


def read(path, expected_hash, links, *, parse=True):
    path = (ROOT / path).resolve()
    location = rel(path)
    raw = path.read_bytes()
    actual = sha256(raw).hexdigest()
    demand(actual == expected_hash and (location not in links or links[location] == actual),
           f"mismatched current input: {location}")
    links[location] = actual
    return json.loads(raw) if parse else None


def store(path, document):
    with path.open("x", encoding="utf-8", newline="\n") as output:
        output.write(json.dumps(document, indent=2, default=str, allow_nan=False) + "\n")
    return digest(path)


def rational(x):
    demand(type(x) in (int, str, Q), "exact rational string/integer expected")
    return Q(x)


def binary(x):
    demand(type(x) in (float, int), "archived field coordinate or scalar must use binary64 encoding")
    return Q.from_float(float(x))


def pair(x):
    demand(type(x) in (list, tuple) and len(x) == 2, "closed interval pair expected")
    a, b = map(rational, x)
    demand(a <= b, "unordered observer interval")
    return a, b


def select(rows, **match):
    answer, = [row for row in rows if all(row[k] == v for k,v in match.items())]
    return answer


def case_id(row):
    return f"{row['endpoint_id']}_{row['checkpoint']}"


def certified_error(eta, z, r, expected_lhs, expected_error):
    eta, z, r, lhs, out = map(rational, (eta,z,r,expected_lhs,expected_error))
    demand(0 <= z < 1 and eta >= 0 and r > 0 and lhs == eta+z*r < r
           and out == eta/(1-z) < r, "insufficient cached contraction or contradictory posterior")
    return out, r


def nodal_total(sigma, posterior, xi, goal):
    sigma, posterior, xi, goal = map(rational, (sigma,posterior,xi,goal))
    demand(min(sigma, posterior, xi) >= 0 and sigma + posterior < goal,
           "root transfer must fit original observer displacement budget, including center shift")
    return LENGTH * (sigma + posterior + xi)


def residual_interval(stored, err, measurement, delta):
    a, b = pair(stored)
    demand(err >= 0 and delta >= 0, "observer radius and allowance must be nonnegative")
    enclosed = (a-err, b+err)
    lo, hi = (x - measurement for x in enclosed)
    res = (max(Q(0), lo, -hi), max(abs(lo), abs(hi)))
    status = "IN" if res[1] <= delta else "OUT" if res[0] > delta else "UNRESOLVED"
    return enclosed, res, status


def signed_gap(positive, negative):
    a,b = pair(positive)
    c,d = pair(negative)
    return a-d, b-c


def observed(folder, status, receipt_hash, os_hash, links):
    os = read(folder / "os_exit_observation.json", os_hash, links)
    end = read(folder / "terminal.json", os["terminal_sha256"], links)
    doc = read(folder / "receipt.json", receipt_hash, links)
    demand(os["status"] == status and integer_zero(os["observed_exit_code"])
           and type(os["pid"]) is int and os["pid"] > 0
           and type(os["unified_exec_session"]) is int and os["unified_exec_session"] > 0
           and end["status"] == "FINISHED" and integer_zero(end["exit_code"])
           and doc.get("pid", os["pid"]) == os["pid"] == end["pid"]
           and receipt_hash == os["receipt_sha256"] == end["receipt_sha256"],
           "external OS0 and independently computed receipt must match")
    return doc


def fixed_authority(links):
    return {p:read(p,h,links,parse=p.suffix==".json") for p,h in PINS.items()}


def archive_fields(docs, links):
    budgets, review, primary, fields, plan = (docs[p] for p in (BUDGETS,BUDGET_REVIEW,PRIMARY,FIELDS,PAIR))
    demand(review["status"] == "PASS_INDEPENDENT_ACTUAL_CONDITIONAL_OBSERVER_BUDGETS"
           and yes(review["all_literal_checks"]) and review["root_error_bounds_proved"] == 0
           and len(review["fields"]) == len(budgets["field_budgets"]) == 36,
           "review of the 36 positive conditional budgets required")
    demand(primary["status"] == "PASS_INDEPENDENT_ACTUAL_NONLINEAR_STORED_FIELD_REPLAY"
           and fields["independent_primary_review_sha256"] == PINS[PRIMARY]
           and integer_zero(fields["independent_replay_observed_os_exit_code"]) and len(fields["fields"]) == 44,
           "same original independently reviewed nonlinear archive required")
    delta = rational(budgets["conditional_source_allowance_um_exact"])
    demand(delta == rational(plan["selection"]["conditional_source_allowance_um_exact"]) and delta > 0,
           "original discrepancy allowance changed")
    b = {case_id(row):row for row in budgets["field_budgets"]}
    wanted = {f"{end}_{step}" for end in ENDPOINTS for step in (*CALIBRATION_CHECKPOINTS, "large_03")}
    demand(set(b) == wanted and docs[COVER/"receipt.json"]["conditional_observer"] == PRIMARY_VARIANT,
           "wrong roles, count, or frozen observer convention")
    for name, row in b.items():
        frozen = select(fields["fields"], endpoint_id=row["endpoint_id"], checkpoint=row["checkpoint"])
        audited = select(review["fields"], endpoint_id=row["endpoint_id"], checkpoint=row["checkpoint"])
        for key in ("field","field_sha256","verification","verification_sha256","role","step"):
            demand(row[key] == frozen[key], f"{name} differs from its archived frozen field")
        demand(all(row[k] == audited[k] for k in ("field","field_sha256","record_id","e_goal_dimensionless_exact"))
               and primary["all_hashes"][row["field"]] == row["field_sha256"]
               and primary["all_hashes"][row["verification"]] == row["verification_sha256"]
               and row["status"] == "POSITIVE_CONDITIONAL_ERROR_TARGET",
               "observer numbers do not match existing independent review")
    theorem = docs[TRANSFER]
    demand(theorem["status"] == "PASS_CONDITIONAL_EXACT_NODAL_TO_OBSERVER_TRANSFER"
           and yes(theorem["checks"]), "reviewed P2 error transfer required")
    for path, h in theorem["reviewed_hashes"].items():
        read(path, h, links, parse=False)
    return b, delta


def root_certificates(docs, budgets, recenter_hash, recenter_os_hash, links):
    roots = {}
    def add(name, leaf, source, lead, keys, *, recentered=False):
        demand(name not in roots and name in budgets, "duplicate or unexpected root state")
        b = budgets[name]
        c = b["center_and_model"]
        error, r = certified_error(*(leaf[k] for k in keys[:5]))
        goal = rational(leaf[keys[5]])
        field_key = "archived_field_sha256" if recentered else "field_sha256"
        demand(leaf[field_key] == b["field_sha256"] and lead[field_key] == b["field_sha256"]
               and leaf["record_id"] == b["record_id"] == lead["record_id"]
               and 0 < r < goal == rational(b["e_goal_dimensionless_exact"])
               and rational(leaf[keys[6]]) == rational(c["exact_target_strength"])
               == rational(c["dose_exact"])*rational(c["gain_exact"]),
               f"{name}: root proof is not the same frozen material, field, exact drive or budget")
        xi = rational(c["node_addition_gap_dimensionless_exact"])
        sigma = Q(0)
        if recentered:
            demand(leaf["strict_old_observer_goal_pass"] is True
                   and leaf["fresh_center_and_entire_cube_recomputed"] is True
                   and leaf["original_field_or_scores_replaced"] is False
                   and leaf["exact_spatial_or_minimum_claimed"] is False
                   and leaf["no_new_preconditioner_or_equilibrium_solves"] is True
                   and rational(leaf["archived_addition_gap_xi0_exact"]) == xi,
                   f"{name}: recentered proof did not retain the original field and operator")
            sigma = rational(leaf["actual_binary_center_shift_sigma_exact"])
            demand(sigma > 0 and rational(leaf["shift_plus_posterior_exact"]) == sigma + error
                   and rational(leaf["physical_nodal_coordinate_error_um_exact"]) == nodal_total(sigma,error,xi,goal),
                   f"{name}: missing independently reviewed actual binary center displacement")
        eps = nodal_total(sigma, error, xi, goal)
        roots[name] = {"authority":source, "cached_proof_leaf_sha256":lead.get("native_sha256",
                         lead.get("native_proof_sha256",lead.get("receipt_sha256"))),
                       "posterior":error, "radius":r, "sigma":sigma, "archived_xi":xi, "epsilon_um":eps}

    four = observed(FOUR, "PASS_OBSERVED_SERIAL_INDEPENDENT_REPLAY_TRUE_OS_EXIT",
                    PINS[FOUR/"receipt.json"],PINS[FOUR/"os_exit_observation.json"],links)
    demand(four == docs[FOUR/"receipt.json"] and four["status"] == "PASS_INDEPENDENT_FOUR_CACHED_LOCAL_ROOTS"
           and yes(four["checks"]) and four["root_certificates_cached_quadrature_only"] == 4,
           "four-pilot root authority not current")
    for row in four["field_results"]:
        name = row["pilot_id"]
        leaf = read(FOUR / (name+".json"), row["receipt_sha256"], links)
        demand(leaf["status"] == row["status"] == "PASS_INDEPENDENT_CACHED_LOCAL_ROOT"
               and leaf["pilot_id"] == name and leaf["rounded_strength_used_in_proof"] is False,
               "four-pilot leaf not strict")
        add(name,leaf,"independent four cached pilot roots",row,(
            "independent_eta_upper_exact","independent_z_upper_exact","radius_exact","independent_selfmap_lhs_upper_exact",
            "independent_posterior_error_upper_exact","science_radius_budget_exact","exact_strength"))

    two = docs[M]
    demand(two["status"] == "PASS_VERIFIED_TWO_ADDITIONAL_M_CACHED_LOCAL_ROOTS"
           and yes(two["checks"]) and two["new_cached_local_root_count"] == len(two["case_results"]) == 2,
           "the two strict M roots must preserve the old unresolved negative radius")
    for row in two["case_results"]:
        name = row["case_id"]
        leaf = read(row["native_proof_file"], row["native_proof_sha256"], links)
        demand(leaf["case_id"] == name and leaf["status"] in ("PASS_INDEPENDENT_M_CACHED_LOCAL_ROOT",
                   "PASS_INDEPENDENT_NEGATIVE_M_SMALLER_CACHED_ROOT")
               and leaf["posterior_error_upper_exact" if "posterior_error_upper_exact" in leaf
                        else "independent_posterior_error_upper_exact"] == row["posterior_error_upper_exact"],
               "M bound is not the independently accepted root")
        narrow = leaf["status"] == "PASS_INDEPENDENT_NEGATIVE_M_SMALLER_CACHED_ROOT"
        if narrow:
            demand(leaf["original_unresolved_test_preserved"] is True, "missing original -M attempt")
        add(name,leaf,"independent M roots with separate narrowed -M",row,(
            "eta_upper_exact" if narrow else "independent_eta_upper_exact","independent_ZQ_upper_exact",
            "radius_exact","strict_selfmap_lhs_upper_exact",
            "posterior_error_upper_exact" if narrow else "independent_posterior_error_upper_exact",
            "observer_goal_exact","exact_strength"))

    waves = docs[WAVES]
    demand(waves["status"] == "PASS_CACHED_WAVE_ACCOUNTING_WITH_TWO_UNRESOLVED"
           and waves["covered_candidate_count"] == len(waves["cases"]) == 30
           and waves["strict_cached_local_root_count"] == 28 and set(waves["unresolved_case_ids"]) == RECENTER_CASES,
           "cached wave coverage and historical misses changed")
    rejected = {}
    for row in waves["cases"]:
        name = row["case_id"]
        leaf = read(row["native_file"], row["native_sha256"], links)
        demand(leaf["case_id"] == name and leaf["status"] == row["status"],
               "wave leaf and retained outcome differ")
        if name in RECENTER_CASES:
            demand(leaf["chosen_radius_exact"] is None and row["original_unresolved_attempts"],
                   "cannot replace an unretained original-radius outcome")
            rejected[name] = row
            continue
        demand(leaf["status"] == "PASS_INDEPENDENT_CACHED_LOCAL_ROOT"
               and leaf["rounded_strength_used_in_proof"] is False,
               "wave root is unresolved")
        chosen = select(leaf["own_radius_attempts"],radius_exact=leaf["chosen_radius_exact"])
        data = {**leaf,**chosen}
        demand(chosen["posterior_error_upper_exact"] == row["posterior_error_upper_exact"]
               and chosen["whole_cube_AD_recomputed_for_this_radius"] is True,
               "wave radius does not match its whole-cube proof")
        add(name,data,"independent cached wave",row,("eta_upper_exact","ZQ_upper_exact","radius_exact",
            "strict_selfmap_lhs_upper_exact","posterior_error_upper_exact","observer_radius_budget_exact","exact_requested_strength"))

    rec = observed(RECENTER, "PASS_ROOT_OBSERVED_INDEPENDENT_RECENTER_REPLAY_TERMINAL", recenter_hash, recenter_os_hash, links)
    demand(rec["status"] == "PASS_INDEPENDENT_TWO_RECENTERED_CACHED_ROOTS_AND_OLD_OBSERVER_BUDGETS"
           and rec["cached_local_roots_certified"] == len(rec["case_results"]) == 2
           and rec["exact_spatial_root_count"] == rec["new_refits_equilibrium_solves_or_R"] == 0
           and rec["archived_observers_predictions_or_scores_replaced"] is False,
           "both recentered independent cached roots are required")
    for row in rec["case_results"]:
        name = row["case_id"]
        leaf = read(row["native_file"],row["native_sha256"],links)
        prior = rejected[name]
        demand(leaf["case_id"] == name and leaf["status"] == row["status"]
               == "PASS_INDEPENDENT_RECENTERED_CACHED_LOCAL_ROOT_AND_OLD_OBSERVER_BUDGET"
               and leaf["old_native_signed_point_sha256"] == prior["native_sha256"]
               and leaf["original_cache_packet_sha256"] == prior["cache_packet_sha256"]
               and leaf["archived_field_sha256"] == prior["field_sha256"]
               and leaf["recentered_cache_packet_sha256"] == row["recentered_cache_packet_sha256"]
               and leaf["actual_binary_center_shift_sigma_exact"] == row["actual_binary_center_shift_sigma_exact"],
               "recentered root cannot silently switch target or erase the original unresolved receipt")
        add(name,leaf,"independent recentered cached root; original observer center unchanged",row,
            ("eta_upper_exact","ZQ_upper_exact","radius_exact","strict_selfmap_lhs_upper_exact",
             "posterior_error_upper_exact","archived_observer_goal_exact","exact_requested_strength"),
            recentered=True)
    demand(set(roots) == set(budgets), "need exactly one proved root error for each of 36 observer states")
    return roots


def frozen_trace(state, budget, root):
    from nodal_observer_error_transfer import endpoint_bounds
    from quadratic_endpoint_observer import _state_controls
    from quadratic_mean_derivatives import _boundary_node_ids
    c = budget["center_and_model"]
    X = [binary(v) for node in state["reference_nodes"] for v in node]
    stored = [binary(v) for node in state["current_nodes"] for v in node]
    vec = list(map(binary,state["free_displacement"]))
    rim = binary(state["parameters"]["outer_radius"])
    fixed = {i for i in range(650) if X[2*(i//2)] == rim or (i%2 == 0 and X[i] == 0)}
    free = [i for i in range(650) if i not in fixed]
    demand(len(X) == len(stored) == 650 and len(free) == len(vec) == 611 and c["free_map"] == free
           and state["parameters"] == c["material_parameters"] and len(state["triangles"]) == 144,
           "archived fixed trace and reduced free map changed")
    Y = X.copy()
    for i, u in zip(free,vec,strict=True):
        Y[i] += u
    demand(all(Q.from_float(float(y)) == s for y,s in zip(Y,stored,strict=True))
           and max(abs(y-s) for y,s in zip(Y,stored,strict=True)) == root["archived_xi"],
           "original storage gap must bound original exact X+u against archived curve")
    controls, height = _state_controls(state,LENGTH)
    ids = _boundary_node_ids(state)
    for face in ("top","bottom"):
        edges = ids[face]
        dofs = [2*edges[0][0],2*edges[-1][-1]]
        desc = c["fixed_trace_endpoints"][face]
        exact_controls = [[(LENGTH*stored[2*n],LENGTH*stored[2*n+1]) for n in edge] for edge in edges]
        demand(desc["trace_nodes"] == edges and desc["radial_dofs"] == dofs and len(edges)==12
               and all(a[-1] == b[0] for a,b in zip(edges,edges[1:]))
               and set(dofs) <= fixed and [Y[i] for i in dofs] == [Q(0),rim]
               and list(map(rational,desc["proved_fixed_radii_um_exact"])) == [Q(0),LENGTH*rim]
               and controls[face] == exact_controls,
               "both archived curves must be conforming and cover a fixed current-radius ROI")
    roi = pair(budget["fixed_roi_um_exact"])
    demand(roi[0] == 0 < roi[1] <= LENGTH*rim and height == rational(c["reference_height_um_exact"]),
           "observer reference and fixed ROI changed")
    tr = endpoint_bounds(controls["top"],controls["bottom"],root["epsilon_um"],
                         fixed_radial_endpoints=True,roi_um=roi)
    demand(tr["status"] == "CONDITIONAL_EXACT_ERROR_TRANSFER", "strict graph error transfer failed")
    for face, side in (("top","apical"),("bottom","basal")):
        baseline = budget["stored_face_derivatives"][face]
        got = tr[side]
        demand(rational(got["radial_derivative_margin_exact"]) ==
               rational(baseline["radial_min_um_exact"])-8*root["epsilon_um"] > 0
               and rational(got["base_graph_Lipschitz_upper_exact"]) == rational(baseline["global_stored_slope_cap_exact"])
               and rational(got["Eulerian_vertical_error_um_exact"]) ==
               rational(baseline["C_um_per_dimensionless_coordinate_exact"])*root["epsilon_um"]/LENGTH,
               "new observer enclosure disagrees with archived exact coefficients")
    return tr,roi


def apply(budget, root, docs, delta, links):
    state = read(budget["field"],budget["field_sha256"],links)
    verification = read(budget["verification"],budget["verification_sha256"],links)
    tr,roi = frozen_trace(state,budget,root)
    name = case_id(budget)
    if budget["role"] == "large_final":
        e = root["epsilon_um"]+rational(tr["apical"]["Eulerian_vertical_error_um_exact"])
        shoulder = rational(docs[PAIR]["large_observable"]["shoulder_radius_um_exact"])
        demand(roi[1] == shoulder == rational(verification["surface"]["shoulder_radius_um_exact"]),
               "large-ROI shoulder is not frozen common radius")
        source = select(docs[PRIMARY]["endpoints"],endpoint_id=budget["endpoint_id"])
        stored = pair(budget["stored_Q_um_exact"])
        demand(stored == pair(source["large_Q_um_exact"]) == pair(verification["surface"]["Q_um_exact"]),
               "large baseline Q is not the independently verified stored one")
        result = {"case_id":name,"root":root,"stored_Q_um_exact":stored,
                  "Q_absolute_error_um_exact":e,"cached_root_Q_um_exact":(stored[0]-e,stored[1]+e)}
    else:
        demand(budget["role"] == "small_record" and roi[1] == binary(state["parameters"]["patch_radius"])*LENGTH,
               "unexpected small-role ROI")
        demand(set(budget["stored_residuals"]) == set(OBS), "all four calibrated quantities are required")
        rows = []
        for obs in OBS:
            frozen = budget["stored_residuals"][obs]
            prediction = pair(frozen["stored_interval_um_exact"])
            measurement = rational(frozen["observed_um_exact"])
            demand(measurement == rational(budget["step"]["observed_exact"][obs]),
                   "measured value differs from frozen plan")
            e = rational(tr["endpoint_error_um_exact"][obs])
            before = residual_interval(prediction,Q(0),measurement,delta)[1][1]
            wide,err,status = residual_interval(prediction,e,measurement,delta)
            demand(before == rational(frozen["stored_worst_residual_um_exact"])
                   and delta-before == rational(frozen["margin_um_exact"]) > 0,
                   "original independent cached observer budget margin altered")
            rows.append({"observable":obs,"observed_um_exact":measurement,"stored_interval_um_exact":prediction,
                         "extra_observer_error_um_exact":e,"cached_root_interval_um_exact":wide,
                         "absolute_residual_interval_um_exact":err,"comparison_at_original_delta":status})
        result = {"case_id":name,"root":root,"comparisons":rows}
    return {**result,"record_id":budget["record_id"],"field_sha256":budget["field_sha256"],
            "endpoint_id":budget["endpoint_id"],"checkpoint":budget["checkpoint"],
            "ROI_um_exact":roi,"apical_and_basal_graph_error_bounds":tr}


def fixtures():
    sys.path.insert(0,str(ROOT/"scripts"))
    from nodal_observer_error_transfer import endpoint_bounds
    def trace(slope, offset):
        return [[(x,slope*x+offset) for x in (Q(0),Q(1,2),Q(1))]]
    e = Q(1,16)
    trans = endpoint_bounds(trace(3,10),trace(0,0),e,
                            fixed_radial_endpoints=True,roi_um=(0,1))
    unresolved = endpoint_bounds(trace(3,10),trace(0,0),Q(1,8),
                                 fixed_radial_endpoints=True,roi_um=(0,1))
    def rejects(fn):
        try:
            fn()
        except ValueError:
            return True
        return False
    tests = {
        "strict_selfmap_is_checked_not_just_claimed":
            certified_error("1/8","1/2","1","5/8","1/4") == (Q(1,4),Q(1))
            and rejects(lambda:certified_error("1/2","1/2","1","1","1")),
        "recenter_shift_and_old_addition_gap_both_count":
            nodal_total("1/10","1/5","1/100","2/5") == LENGTH*Q(31,100)
            and rejects(lambda:nodal_total("1/5","1/5","1/100","2/5")),
        "signed_apical_difference_is_not_unsigned_contrast":
            signed_gap(("-5","-4"),("-7","-6")) == (Q(1),Q(3)),
        "strict_radial_guard_and_both_surface_errors_plus_axis_term":
            trans["status"] == "CONDITIONAL_EXACT_ERROR_TRANSFER"
            and rational(trans["apical"]["Eulerian_vertical_error_um_exact"])+e == 6*e
            and rational(trans["endpoint_error_um_exact"]["maximum_shortening_um"]) == Q(25,4)*e
            and unresolved["status"] == "UNRESOLVED_ERROR_TRANSFER",
        "calibration_uses_worst_signed_endpoint_absolute_residual":
            residual_interval(("-2","-1"),Q(1,16),Q(0),Q(2))[2] == "UNRESOLVED"
            and residual_interval(("-2","-1"),Q(0),Q(0),Q(2))[2] == "IN"
            and residual_interval(("-4","-3"),Q(0),Q(0),Q(2))[2] == "OUT",
        "reject_float_or_bool_exact_inputs":rejects(lambda:rational(False)) and rejects(lambda:rational(0.5)),
    }
    return tests


def prepare():
    links = {}
    docs = fixed_authority(links)
    checks = fixtures()
    demand(yes(checks), "cold exact synthetic observer fixture fails")
    h = digest(Path(__file__))
    p = store(HERE/"protocol.json", {
        "status":"FROZEN_CACHED_ROOT_36_FIELD_OBSERVER_APPLICATION", "source_sha256":h,
        "fixed_authority_sha256":links,"external_recenter_authority":["independent_native_receipt_sha256","independent_os_exit_sha256"],
        "fields":36,"primary_calibration_comparisons":136,"large_observables":2,"comparison_convention":PRIMARY_VARIANT,
        "frozen_allowance_exact":docs[BUDGETS]["conditional_source_allowance_um_exact"],
        "cold_actual_fields_read":0,"cold_numerical_mechanics_jobs":0,
        "scope":"Fixed independent cached-quadrature root error-transfer only. Experimental and exact-spatial claims remain separate.",
    })
    pf = store(HERE/"preflight.json",{"status":"PASS","source_sha256":h,"protocol_sha256":p,
         "checks":checks,"cold_actual_fields_read":0,"cold_numerical_mechanics_jobs":0})
    print(json.dumps({"source_sha256":h,"protocol_sha256":p,"preflight_sha256":pf}))


def run(args):
    links = {}
    pf = read(HERE/"preflight.json",args.preflight_sha256,links)
    p = read(HERE/"protocol.json",args.protocol_sha256,links)
    gate = read(HERE/"root_source_review.json",args.source_review_sha256,links)
    h = digest(Path(__file__))
    demand(pf["status"] == "PASS" and yes(pf["checks"]) and gate["status"] == "PASS" and yes(gate["checks"])
           and h == gate["source_sha256"] == pf["source_sha256"] == p["source_sha256"]
           and gate["protocol_sha256"] == pf["protocol_sha256"] == args.protocol_sha256
           and gate["preflight_sha256"] == args.preflight_sha256,
           "cold source, protocol, preflight and root review must all match")
    fixed = {}
    docs = fixed_authority(fixed)
    demand(p["fixed_authority_sha256"] == fixed,"fixed cold pins changed")
    links.update(fixed)
    # Late-bound but independently observed recenter proofs are only opened in actual phase.
    budgets, delta = archive_fields(docs,links)
    roots = root_certificates(docs,budgets,args.recenter_sha256,args.recenter_os_sha256,links)
    sys.path.insert(0,str(ROOT/"scripts"))
    rows = [apply(b,roots[key],docs,delta,links) for key,b in sorted(budgets.items())]
    counts = {end:{k:0 for k in ("IN","OUT","UNRESOLVED")} for end in ENDPOINTS}
    for row in rows:
        for item in row.get("comparisons",[]):
            counts[row["endpoint_id"]][item["comparison_at_original_delta"]] += 1
    demand(all(sum(c.values()) == 68 for c in counts.values()),"exactly 68 comparisons per endpoint required")
    a,b = [select(rows,endpoint_id=end,checkpoint="large_03") for end in ENDPOINTS]
    gap = signed_gap(a["cached_root_Q_um_exact"],b["cached_root_Q_um_exact"])
    demand(signed_gap(a["stored_Q_um_exact"],b["stored_Q_um_exact"])
           == pair(docs[BUDGETS]["stored_large_signed_gap_um_exact"])
           == pair(docs[PRIMARY]["stored_large_Q_signed_difference_um_exact"]),
           "large cached root intervals changed archived signed gap definition")
    result = {
        "status":"CACHED_QUADRATURE_PRIMARY_OBSERVER_APPLICATION_AWAITING_INDEPENDENT_READBACK",
        "scope":"Static primary comparison at the two frozen fixed parameter points for the binary-cache quadrature operator.",
        "cached_root_fields":len(rows),"primary_calibration_comparisons":sum(sum(c.values()) for c in counts.values()),
        "large_Q_values":2,"primary_observer":PRIMARY_VARIANT,"frozen_allowance_um_exact":delta,
        "per_endpoint_counts":counts,"cases":rows,"cached_root_positive_minus_negative_large_Q_um_exact":gap,
        "strict_positive_cached_gap":gap[0]>0,"source_and_actual_sha256":links,
        "exact_spatial_root_observer_results":0,"alternate_observer_robustness_claimed":False,
        "experimental_compatibility_or_PDE_claimed":False,"full_five_gate_goal_achieved":False,
        "new_equilibrium_cache_AD_force_or_moment_jobs":0,
    }
    print(json.dumps({"results_sha256":store(HERE/"results.json",result)}))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    choices = p.add_subparsers(dest="mode",required=True)
    choices.add_parser("prepare")
    actual = choices.add_parser("run")
    for item in ("protocol-sha256","preflight-sha256","source-review-sha256","recenter-sha256","recenter-os-sha256"):
        actual.add_argument("--"+item,required=True)
    args=p.parse_args()
    if args.mode=="prepare":
        prepare()
    else:
        run(args)

