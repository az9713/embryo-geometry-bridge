"""Independent exact replay of conditional observer budgets. No force kernels."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
PRIMARY = BASE.parents[1]
ROOT = PRIMARY.parents[1]
PRE = BASE / "preparation"
RUN = BASE / "actual"
NATIVE = PRIMARY / "independent_actual_preparation/actual_replay/receipt.json"
MANIFEST = PRIMARY / "field_manifest.json"
PR = PRIMARY / "preparation/protocol.json"
SRC = ROOT / "scripts/transfer_ambiguity_error_budgets.py"
L = Q(67, 2)
OBS = ("maximum_shortening_um", "region_averaged_shortening_um",
       "apical_max_shortening_um", "basal_max_shortening_um")
PINS = {
    RUN / "results.json": "689948766db950ee4ddf399e6731aea08890cc7a1f84c5615d19e7ded68b0872",
    RUN / "started.json": "9e56d7be57cfc28ba8db78c3dfc4b0c5dfb0c0592ac84a56a738523ce38ade36",
    RUN / "terminal.json": "83abadac599cba377034a522aa5c160de176b4b8b6b91f6223883415495d9928",
    SRC: "410a4944798b1f42f0360cbdf236d9da68f8de2615d93f860bccdd8be54d8859",
    PRE / "protocol.json": "332b66e28a9c990b0e6c59953d36c841d5c047e8bb9a5b9d47467e97598a3d27",
    PRE / "preflight.json": "4588b3dd1779d00e85f80e8dec62d3862a9605c883f73ff57c1ada78374b0e3f",
    PRE / "prelaunch_independent_review.json": "5a1aed542a55c2a42f99252917d53bb7d8eeda5affba6086cb5e4b4fec6834ff",
    PR: "9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94",
    NATIVE: "c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd",
    MANIFEST: "308ca4c05e476ba6a14cbc7b8d34018983b8f11263e2c224488672b35eb38ac0",
}


def digest(file):
    return sha256(file.read_bytes()).hexdigest()


def read(file):
    return json.loads(file.read_bytes())


def qmap(row):
    return {key: Q(value) for key, value in row.items()}


def pair(values):
    left, right = map(Q, values)
    assert left <= right
    return left, right


def stop_radius(x, t, k):
    p0, pm, p1 = x
    # Differentiate the power basis p0 + b*t + a*t*t.
    b = 4 * (pm - p0) - (p1 - p0)
    a = 2 * (p1 - 2 * pm + p0)
    return b + 2 * a * t


def layout(state, material, budget):
    parameters = state["parameters"]
    assert parameters == material["parameters"]
    assert parameters["reference_curvature"] == 0
    p = [tuple(map(Q, x)) for x in state["reference_nodes"]]
    x = [tuple(map(Q, y)) for y in state["current_nodes"]]
    triangles = state["triangles"]
    assert len(p) == len(x) == 325 and len(triangles) == 144
    rim = Q(parameters["outer_radius"])
    height = Q(parameters["height"])
    flat_p = [v for coords in p for v in coords]
    flat_x = [v for coords in x for v in coords]
    free = [2 * i + k for i, (r, _) in enumerate(p) for k in (0, 1)
            if r != rim and (k != 0 or r != 0)]
    assert free == budget["center_and_model"]["free_map"]
    assert len(free) == len(state["free_displacement"]) == 611
    free_offsets = dict(zip(free, map(Q, state["free_displacement"]), strict=True))
    gaps = [abs(y - r - free_offsets.get(i, 0))
            for i, (r, y) in enumerate(zip(flat_p, flat_x, strict=True))]
    xi = max(gaps)
    provenance = budget["center_and_model"]
    assert xi == Q(provenance["node_addition_gap_dimensionless_exact"])
    assert provenance["largest_gap_flat_indices"] == [i for i, g in enumerate(gaps) if g == xi]
    assert all(g == 0 for i, g in enumerate(gaps) if i not in free_offsets)
    assert Q(provenance["reference_height_um_exact"]) == height * L
    assert state["reached_target"] is state["converged"] is True
    strength = Q(state["activation_strength"])
    gain = Q(material["gain_exact"])
    dose = Q(budget["step"]["dose_exact"])
    exact_strength = gain * dose
    assert strength == Q(float(exact_strength)) == Q(budget["step"]["solver_strength_binary64_exact"])
    assert Q(budget["step"]["requested_strength_exact"]) == exact_strength > 0
    assert Q(budget["step"]["binary64_rounding_difference_exact"]) == strength - exact_strength
    assert Q(provenance["dose_exact"]) == dose
    assert Q(provenance["gain_exact"]) == gain
    assert Q(provenance["exact_target_strength"]) == exact_strength
    assert Q(provenance["solver_binary64_strength_exact"]) == strength
    assert Q(provenance["binary64_strength_rounding_error_exact"]) == strength - exact_strength
    assert provenance["material_parameters"] == parameters
    assert provenance["material_parameter_storage_exact"] == material["parameter_storage_exact"]
    assert provenance["active_elements"] == [i for i, flag in enumerate(state["activation_elements"]) if flag]
    assert len(provenance["active_elements"]) == 24

    # Extract sole-incident top and bottom edges directly from the P2 mesh.
    uses = {}
    for a, b, c, d, e, f in triangles:
        for i, j, m in ((a, b, d), (b, c, e), (c, a, f)):
            if p[i][0] > p[j][0]:
                i, j = j, i
            uses.setdefault(frozenset((i, j)), []).append((i, m, j))
    values = {}
    for name, z in (("top", height), ("bottom", Q(0))):
        trace = sorted((edges[0] for edges in uses.values() if len(edges) == 1
                        and all(p[i][1] == z for i in edges[0])), key=lambda e: p[e[0]][0])
        assert len(trace) == 12 and all(a[2] == b[0] for a, b in zip(trace, trace[1:]))
        fixed = provenance["fixed_trace_endpoints"][name]
        assert fixed["trace_nodes"] == list(map(list, trace))
        ends = (trace[0][0], trace[-1][-1])
        assert fixed["radial_dofs"] == [2*i for i in ends]
        assert all(2*i not in free_offsets for i in ends)
        assert [p[i][0] for i in ends] == [x[i][0] for i in ends] == [0, rim]
        assert list(map(Q, fixed["proved_fixed_radii_um_exact"])) == [0, L * rim]
        values[name] = (tuple(tuple(tuple(c * L for c in x[i]) for i in edge) for edge in trace), rim * L)
    return xi, values


def face_checks(face, controls, expected_right, roi):
    trace, rim = controls
    assert 0 < expected_right <= rim
    assert trace[0][0][0] == 0 and trace[-1][-1][0] == rim
    mins, slopes = [], []
    for k, (edge, saved) in enumerate(zip(trace, face["edges"], strict=True)):
        assert saved["edge"] == k
        radial, vertical = map(tuple, zip(*edge, strict=True))
        dr = [stop_radius(radial, t, k) for t in (0, 1)]
        dz = [stop_radius(vertical, t, k) for t in (0, 1)]
        assert [radial[0], radial[-1]] == list(map(Q, saved["radial_start_end_um_exact"]))
        assert dr == list(map(Q, saved["Rprime_local_endpoints_um_exact"]))
        assert dz == list(map(Q, saved["Zprime_local_endpoints_um_exact"]))
        mn = min(dr)
        assert mn > 0
        ls = max(map(abs, dz)) / mn
        assert Q(saved["radial_derivative_min_um_exact"]) == mn
        assert Q(saved["stored_slope_cap_exact"]) == ls
        mins.append(mn)
        slopes.append(ls)
    minimum, slope = min(mins), max(slopes)
    coefficient = 5 * L * (1 + slope) / 4
    assert Q(face["radial_min_um_exact"]) == minimum
    assert Q(face["global_stored_slope_cap_exact"]) == slope
    assert Q(face["C_um_per_dimensionless_coordinate_exact"]) == coefficient
    check_transfer(face["zero_radius_transfer_check"], Q(0), Q(0), minimum, slope, roi)
    return minimum, slope, coefficient


def check_transfer(transfer, epsilon, error, radial_min, slope, roi):
    assert transfer["status"] == "CONDITIONAL_EXACT_ERROR_TRANSFER"
    assert transfer["fixed_radial_endpoints_assumed"] is True
    assert list(map(Q, transfer["ROI_um_exact"])) == list(roi)
    assert Q(transfer["nodal_coordinate_error_um_exact"]) == epsilon
    assert Q(transfer["parametric_coordinate_error_um_exact"]) == Q(5, 4) * epsilon
    assert Q(transfer["radial_derivative_margin_exact"]) == radial_min - 8 * epsilon >= radial_min / 2 > 0
    assert Q(transfer["base_graph_Lipschitz_upper_exact"]) == slope
    assert Q(transfer["Eulerian_vertical_error_um_exact"]) == error


def primary_terms(native, result, endpoints):
    outs = {row["endpoint_id"]: row for row in result["endpoints"]}
    refs = {row["endpoint_id"]: row for row in native["endpoints"]}
    assert list(outs) == list(refs) == list(endpoints)
    plus, minus = list(outs)
    endpoints_q = {eid: pair(row["large_Q_um_exact"]) for eid, row in outs.items()}
    assert endpoints_q == {eid: pair(row["large_Q_um_exact"]) for eid, row in refs.items()}
    for eid, row in outs.items():
        assert list(map(Q, row["small_residual"]["all68_residuals"]["required_allowance_um_exact"])) == list(
            map(Q, refs[eid]["required_allowance_um_exact"]))
    a, b = endpoints_q[plus]
    c, d = endpoints_q[minus]
    gap = (a-d, b-c)
    assert gap == pair(native["stored_large_Q_signed_difference_um_exact"])
    assert gap == pair(result["large_Q_positive_direction_minus_negative_direction_um_exact"])
    g = max(gap[0], -gap[1], Q(0))
    assert g > 0
    return outs, endpoints_q, gap, g


def main():
    for path, expected in PINS.items():
        assert digest(path) == expected, path
    budget, starting, terminal = (read(RUN / f"{name}.json") for name in ("results", "started", "terminal"))
    prep, preflight, gate = (read(PRE / f"{name}.json")
                             for name in ("protocol", "preflight", "prelaunch_independent_review"))
    manifest, primary, native = map(read, (MANIFEST, PR, NATIVE))
    assert terminal["status"] == "FINISHED" and starting["status"] == "RUNNING"
    assert starting["pid"] == terminal["pid"] == 35592
    assert terminal["results_sha256"] == PINS[RUN / "results.json"]
    assert starting["utc"] < terminal["utc"]
    bindings = dict(zip(
        ("source_sha256","protocol_sha256","preflight_sha256","primary_protocol_sha256",
         "primary_actual_review_sha256","whole_44_field_manifest_sha256","prelaunch_review_sha256"),
        (PINS[x] for x in (SRC,PRE/"protocol.json",PRE/"preflight.json",PR,NATIVE,MANIFEST,
                           PRE/"prelaunch_independent_review.json")), strict=True))
    assert budget["bindings"] == starting["bindings"] == terminal["bindings"] == bindings
    assert gate["status"] == preflight["status"] == "PASS"
    assert all(gate[k] == v for k, v in bindings.items() if k != "prelaunch_review_sha256")
    assert gate["checks"] and all(v is True for v in gate["checks"].values())
    assert preflight["checks"] and all(v is True for v in preflight["checks"].values())
    assert digest(PRE / "transfer_ambiguity_error_budgets_frozen.py") == PINS[SRC]
    assert native["exact_equilibrium_position_bounded"] is native["source_compatible_witness"] is False
    assert native["new_equilibrium_solves"] == 0
    assert manifest["independent_replay_observed_os_exit_code"] == 0
    assert manifest["independent_primary_review_sha256"] == PINS[NATIVE]
    assert native["actual_results_sha256"] == manifest["primary_results_sha256"]
    nr_terminal = read(NATIVE.with_name("terminal.json"))
    assert nr_terminal["exit_code"] == 0 and nr_terminal["receipt_sha256"] == PINS[NATIVE]
    hashes = (prep["fixed_design_hashes"] | prep["dependencies"] | prep["terminal_authority_hashes"]
              | primary["source_hashes"] | native["all_hashes"])
    hashes.update({f.relative_to(ROOT).as_posix(): h for f, h in PINS.items()})
    for file, h in hashes.items():
        assert digest(ROOT / file) == h, file

    delta = Q(budget["conditional_source_allowance_um_exact"])
    assert delta == Q(primary["selection"]["conditional_source_allowance_um_exact"])
    assert delta == Q(prep["frozen_conditional_allowance_um_exact"])
    assert prep["observer"] == primary["observable_variant"] == "current_disk__spatial_max"
    assert L == Q(prep["length_unit_um_exact"]) == Q(primary["large_observable"]["length_unit_um_exact"])
    assert Q(prep["sigma_exact"]) == Q(prep["radial_margin_spending_fraction_exact"]) == Q(1,2)
    ep = {row["endpoint_id"]: row for row in primary["endpoints"]}
    original = read(PRIMARY / "actual/results.json")
    outputs, endpoints_q, gap, g = primary_terms(native, original, ep)
    assert pair(budget["stored_large_signed_gap_um_exact"]) == gap

    expected = []
    for eid, design in ep.items():
        small = [(f"small_{k:02}", "small_record", v) for k, v in enumerate(design["small_records"], 1)]
        roles = ("large_zero","large_checkpoint","large_checkpoint","large_final")
        large = [(f"large_{k:02}", role, v) for k, (role, v) in
                 enumerate(zip(roles, design["large_checkpoints"], strict=True))]
        expected.extend([(eid,"small_zero","small_zero",design["small_zero"])] +
                        [(eid,*entry) for entry in small+large])
    headers = ("endpoint_id","checkpoint","role","step")
    assert len(expected) == len(manifest["fields"]) == native["checkpoints"] == 44
    assert expected == [tuple(row[k] for k in headers) for row in manifest["fields"]]
    for row, planned in zip(manifest["fields"], prep["frozen_44_checkpoint_plan"], strict=True):
        for k in (*headers,"field","verification"):
            assert row[k] == planned[k]
        for k in ("field","verification"):
            assert hashes[row[k]] == row[k+"_sha256"]

    selected = [r for r in manifest["fields"] if r["role"] in ("small_record","large_final")]
    assert len(selected) == len(budget["field_budgets"]) == 36
    margins, qcaps, pilot, audit_rows = {}, [], [], []
    for entry, row in zip(selected, budget["field_budgets"], strict=True):
        for k in (*headers,"field","verification","field_sha256","verification_sha256"):
            assert row[k] == entry[k]
        eid, checkpoint, role, step = (row[key] for key in headers)
        small = role == "small_record"
        rid = step["record_id"] if small else primary["matched_input"]["large_record_id"]
        assert row["record_id"] == rid
        mat = ep[eid]["materials"]["small" if small else "large"]
        state, verify = (read(ROOT / row[k]) for k in ("field","verification"))
        assert state["integration_policy"] == primary["numerical_policy"]["solve_quadrature"]
        assert Q(entry["actual_strength_binary64_exact"]) == Q(state["activation_strength"])
        xi, curves = layout(state, mat, row)
        right = Q(mat["parameters"]["patch_radius"]) * L if small else Q(
            primary["large_observable"]["shoulder_radius_um_exact"])
        roi = (Q(0), right)
        assert tuple(map(Q, row["fixed_roi_um_exact"])) == roi
        faces = {name: face_checks(row["stored_face_derivatives"][name], curve, right, roi)
                 for name, curve in curves.items()}
        ca, cb = faces["top"][2], faces["bottom"][2]
        rcap = min(values[0] for values in faces.values()) / 16 / L
        assert Q(row["monotonicity_total_dimensionless_cap_exact"]) == rcap
        if small:
            assert verify["variant"] == prep["observer"]
            coefficients = dict(zip(OBS,(ca+cb,ca+cb,ca,cb),strict=True))
            assert coefficients == qmap(row["conditional_observer_coefficient_um_exact"])
            assert Q(row["small_margin_fraction_exact"]) == Q(1, 2)
            fobs = verify["observer_intervals_um_exact"]
            raw_rows = outputs[eid]["small_residual"]["all68_residuals"]["constraint_rows"]
            raw = {k: next(r for r in raw_rows if r["record_id"] == rid and r["observable"] == k)
                   for k in OBS}
            allocations = {}
            assert set(fobs) == set(OBS) == set(row["stored_residuals"])
            for name, v in row["stored_residuals"].items():
                pred, o = pair(fobs[name]), Q(step["observed_exact"][name])
                lo, hi = pred
                res_min = 0 if lo <= o <= hi else min(abs(o-lo), abs(o-hi))
                res_max = max(abs(o-x) for x in pred)
                rr = raw[name]
                assert Q(rr["observed_exact"]) == Q(v["observed_um_exact"]) == o
                assert pair(rr["prediction_interval_exact"]) == pair(v["stored_interval_um_exact"]) == pred
                assert Q(rr["minimum_interval_residual_exact"]) == res_min
                assert Q(rr["maximum_interval_residual_exact"]) == Q(v["stored_worst_residual_um_exact"]) == res_max
                assert Q(rr["discrepancy_allowance_exact"]) == delta
                margin = delta - res_max
                assert Q(v["margin_um_exact"]) == margin > 0
                margins[eid,rid,name] = margin
                allocations[name] = margin / 2
            assert allocations == qmap(row["per_observer_allocations_um_exact"])
        else:
            assert endpoints_q[eid] == pair(verify["surface"]["Q_um_exact"]) == pair(row["stored_Q_um_exact"])
            assert Q(row["stored_gap_distance_from_zero_lower_um_exact"]) == g
            coefficients, allocations = {"Q": L+ca}, {"Q": g/4}
            assert Q(row["balanced_error_allocation_um_exact"]) == g/4
            assert Q(row["conditional_Q_error_coefficient_um_exact"]) == L+ca

        caps = {k: allocations[k] / c for k, c in coefficients.items()}
        winner = min(caps, key=caps.get)
        t = min(caps[winner], rcap)
        error = t - xi
        assert error == Q(row["e_goal_dimensionless_exact"]) > 0
        assert row["status"] == "POSITIVE_CONDITIONAL_ERROR_TARGET"
        assert row["e_goal_dimensionless_display"] == float(error)
        assert caps == qmap(row["observation_total_dimensionless_cap_by_observable_exact"])
        assert row["limiting_observable"] == winner
        assert row["limiting_guard"] == ("radial_monotonicity" if rcap < caps[winner] else "observer_error")
        assert t == Q(row["total_dimensionless_coordinate_cap_exact"])
        assert t * L == Q(row["physical_nodal_error_goal_um_exact"])
        spent = {k: c*t for k,c in coefficients.items()}
        assert spent == qmap(row["conditional_observer_error_caps_um_exact"])
        assert all(cap <= allocations[k] for k, cap in spent.items())
        transfer = row["conditional_transfer_at_design_cap"]
        assert transfer["status"] == "CONDITIONAL_EXACT_ERROR_TRANSFER"
        for side, name in (("apical","top"),("basal","bottom")):
            minimum, slope, c = faces[name]
            check_transfer(transfer[side], L*t, c*t, minimum, slope, roi)
        assert qmap(transfer["endpoint_error_um_exact"]) == {
            k:c*t for k,c in zip(OBS,(ca+cb,ca+cb,ca,cb),strict=True)}
        if small:
            remaining = {k: margins[eid,rid,k] - cap for k,cap in spent.items()}
            assert remaining == qmap(row["conditional_remaining_allowance_um_exact"])
            assert all(kept >= margins[eid,rid,k]/2 > 0 for k,kept in remaining.items())
        else:
            qcaps.append(spent["Q"])
        audit_rows.append({k:row[k] for k in ("record_id","endpoint_id","checkpoint","field","field_sha256",
                                              "e_goal_dimensionless_exact","limiting_guard","limiting_observable")}
                          | {"binary_gap_exact": str(xi), "e_goal_dimensionless_display": float(error),
                             "root_displacement_radius_um_display": float(L*error)})

    assert len(margins) == budget["stored_primary_constraint_count"] == 136
    worst = min(margins, key=margins.get)
    least_margin = margins[worst]
    assert least_margin == Q(budget["minimum_stored_margin_um_exact"])
    assert delta - least_margin == pair(original["both_endpoints_required_allowance_um_exact"])[1]
    for eid in ep:
        rows = [r for r in audit_rows if r["endpoint_id"] == eid]
        candidate = min((r for r in rows if r["checkpoint"].startswith("small_")),
                        key=lambda r: (Q(r["e_goal_dimensionless_exact"]),r["record_id"]))
        raw_record = min((k for k in margins if k[0] == eid), key=lambda k:(margins[k],k[1]))[1]
        final = next(r for r in rows if r["checkpoint"] == "large_03")
        pilot.extend([candidate | {"record_with_smallest_raw_stored_observer_margin": raw_record},
                      final])
    for expected_pilot, reported in zip(pilot,budget["four_pilot_selection"],strict=True):
        assert all(expected_pilot[k]==v for k,v in reported.items() if k != "role")
        assert reported["role"] == ("large_final" if expected_pilot["checkpoint"]=="large_03" else "small_record")
    lhs, rhs = gap[0]-sum(qcaps), gap[1]+sum(qcaps)
    reserve = max(lhs,-rhs,Q(0))
    claimed = budget["conditional_large_gap_reserve"]
    assert (lhs,rhs) == pair(claimed["gap_if_each_future_error_meets_its_cap_um_exact"])
    assert reserve == Q(claimed["remaining_distance_from_zero_lower_um_exact"]) >= g/2 > 0
    assert budget["status"] == "COMPLETE_POSITIVE_CONDITIONAL_36_FIELD_BUDGETS"
    assert budget["exact_root_error_bounds_proved"] == budget["new_force_evaluations"] == budget["new_equilibrium_solves"] == 0
    out = {
        "status": "PASS_INDEPENDENT_ACTUAL_CONDITIONAL_OBSERVER_BUDGETS",
        "scope": "Sufficient radii for future exact-integral finite-element root certificates for the fixed inputs, prescribed binary model and primary observer. No equilibrium, physical stability, full-body geometry, continuum, or experimental discrepancy certificate.",
        "observed_budget_terminal": {"pid":35592,"session":60829,"exit_code":0,"observer":"root"},
        "source_and_input_hashes": hashes,
        "audit_source_sha256": digest(Path(__file__)),
        "all_literal_checks": {
            "actual_and_primary_terminals_and_frozen_hashes": True,
            "native_44_field_and_88_field_verification_payload_bindings": True,
            "independent_36_core_injection_and_gap_calculations": True,
            "all_864_p2_trace_derivative_and_lipschitz_bounds": True,
            "independent_136_residual_margins_under_frozen_allowance": True,
            "strict_radial_half_reserve_and_conservative_observer_caps": True,
            "both_large_Q_caps_preserve_at_least_half_stored_gap": True,
            "four_minimax_pilot_paths_and_hashes_match": True,
        },
        "numerical_observer_engines_reused": [],
        "premise_reuse": "Existing independently verified exact observer-value intervals are inputs. Their geometry values were not re-integrated during this budget review.",
        "fields":audit_rows,
        "four_pilot_fields":pilot,
        "smallest_margin_constraint":list(worst),
        "smallest_stored_margin_um_exact": str(least_margin),
        "smallest_stored_margin_um_display": float(least_margin),
        "stored_large_gap_distance_um_display": float(g),
        "conditional_remaining_large_gap_um_display": float(reserve),
        "root_error_bounds_proved": 0,
        "new_equilibrium_solves": 0,
        "new_force_evaluations": 0,
    }
    file = HERE/"receipt.json"
    with file.open("x", encoding="utf-8") as stream:
        json.dump(out,stream,indent=2)
        stream.write("\n")
    print(digest(file), len(hashes), float(least_margin), float(g), float(reserve))
    for candidate in pilot:
        print(*(candidate[k] for k in ("record_id","e_goal_dimensionless_display",
                                        "root_displacement_radius_um_display","field","field_sha256")))


if __name__ == "__main__":
    main()
