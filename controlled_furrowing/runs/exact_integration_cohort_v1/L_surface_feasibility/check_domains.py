"""Check L's original normal-cube surface domains and two-mass scalar tails."""

import json
import os
import sys
import time
import traceback
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np
import exact_uniform_surface_bridge_B_safe as U

T, S, JET, core = U.T, U.S, U.JET, U.core
F_LIMIT = Q(1, 10**22)
H_LIMIT = Q(1, 10**16)
MAX_ORDER = 128
COMPARISON_ORDER = 28
MAX_DEPTH = 8
GOAL = Q(1, 8)
RECORD = "optoGEF2_small_col_L"
INVENTORY = "runs/exact_integration_cohort_v1/inventory/inventory.json"
PREPARATION = "runs/exact_integration_bridge_v1/target_preparation_safe_v2/protocol.json"
B_POLICY = "runs/exact_integration_bridge_v1/uniform_surface_B_safe_v2/protocol.json"
INVENTORY_SHA = "e658b641576d7966466a57d2b66b5d26c2e5913700af9d2b5d6aded9877c725b"
PREPARATION_SHA = "68164b3a922854b05d55c57d6872e92451e176484d7f5bf2ba3127101819850b"


def path_of(module_file):
    return Path(module_file).resolve().relative_to(ROOT).as_posix()


def require(condition, why):
    if not condition:
        raise ValueError(why)


def source_and_inputs():
    require(U.sha(ROOT / INVENTORY) == INVENTORY_SHA, "Verified inventory hash")
    require(U.sha(ROOT / PREPARATION) == PREPARATION_SHA, "Frozen target preparation hash")
    inventory = U.read(ROOT / INVENTORY)
    prepared = U.read(ROOT / PREPARATION)
    row, = [row for row in prepared["records"] if row["record_id"] == RECORD]
    summary, = [row for row in inventory["records"] if row["record_id"] == RECORD]
    require(summary["proof_center_context_validated"] is True and
            summary["new_exact_integral_root_certified"] is False, "L inventory scope")
    require(all(row[k] == summary[k] for k in summary.keys() & row.keys()
                if k != "new_exact_target_root_certified"), "Inventory/preparation identity")
    policy = U.read(ROOT / B_POLICY)
    expected = {
        **inventory["all_input_sha256"],
        **policy["source_and_input_hashes"],
        INVENTORY: INVENTORY_SHA,
        PREPARATION: PREPARATION_SHA,
    }
    names = [INVENTORY, PREPARATION, B_POLICY, path_of(__file__)]
    names += [path_of(module.__file__) for module in (U, S, T, JET, core, S.M, U.L)]
    names += [row[key + "_file"] for key in ("root_certificate", "independent_root_capsule",
                                            "point_cache", "predictor", "proof_center_state")]
    names += [summary[k] for k in ("old_nodal_box_observer_file",
                                   "old_nodal_box_application_review_file")]
    bindings = {name: U.sha(ROOT / name) for name in names}
    require(all(digest == expected[name] for name, digest in bindings.items() if name in expected),
            "Selected frozen input or source changed")
    packet = U.read(ROOT / row["point_cache_file"])
    target = T.Target(packet)
    require(packet["input_state_sha256"] == row["proof_center_state_sha256"] and
            packet["frozen_binary_cache"]["sha256"] == row["frozen_binary_cache_sha256"],
            "Frozen embedded cache or center binding")
    with np.load(ROOT / row["predictor_file"], allow_pickle=False) as archive:
        predictor = archive["S"]
    radius = Q(row["selected_corrector_radius_exact"])
    require(radius == Q(summary["selected_corrector_radius_exact"]) > 0, "Original selected normal radius")
    return bindings, row, target, predictor


def point_spread_witness(target, face, edge, depth):
    fields = target.surface(face, edge)
    samples = []
    for k in range(2**depth + 1):
        s = Q(k, 2**depth)
        a, b, c = T.evaluate_surface(fields, s)[0]
        samples.append((a * a * (b * b + c * c), s))
    lower, s_lower = min(samples)
    upper, s_upper = max(samples)
    rho_lower = (upper - lower) / (upper + lower)
    return {
        "G_low_value_exact": str(lower),
        "G_high_value_exact": str(upper),
        "s_at_low_exact": str(s_lower),
        "s_at_high_exact": str(s_upper),
        "lower_bound_on_global_centered_rho_exact": str(rho_lower),
        "no_global_center_can_achieve_rho_le_1_over_8": rho_lower > GOAL,
    }


def check_edge(target, nodes, row, face, edge, thresholds):
    started = time.perf_counter()
    fields = target.surface(face, edge, nodes=nodes)
    domain = U.uniform_G_bounds(fields, max_depth=MAX_DEPTH, rho_goal=GOAL)
    rho = Q(domain["rho_exact"])
    require(domain["status"] == "CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN" and
            Q(domain["hoop_lower_exact"]) > 0 and Q(domain["G_bounds_exact"][0]) > 0 and
            0 <= rho < 1, "Certified positive hoop and G with convergent global expansion")
    domain_seconds = time.perf_counter() - started
    jet = JET.surface_squared_jacobian_jet(fields["q_physical_R"], fields["B_physical_R"])
    bounds = fields["Rlo"], fields["Rhi"]
    norms = [U.poly_norm_upper(poly, *bounds) for poly in jet["gradient"]]
    h_norms = {ij: U.poly_norm_upper(poly, *bounds) for ij, poly in jet["hessian_upper"].items()}
    require(bool(norms) and bool(h_norms), "Edge has nonempty local free jet")

    # Analytic zeroth material mass and finite zeroth auxiliary mass only.
    # There are no weighted force/Hessian cache sums or moment defects here.
    mass_integral = target.normalization() * U.A((bounds[1]**2 - bounds[0]**2) / 2)
    cache_weights = target.a[f"surface_{face}_weights"][edge]
    require(all(T.rational(weight) > 0 for weight in cache_weights), "Positive auxiliary weights")
    rule_mass_exact = sum(map(T.rational, cache_weights)) / target.energy_scale
    mass = core.upper(mass_integral + U.A(rule_mass_exact))
    K, tension, K_exact, tension_exact = U.material_boxes(row, face)
    q0 = Q(domain["q0_exact"])
    root = U.A(q0).sqrt()
    dmax = max(map(abs, U.rational_hull(tension - K)))
    prefactor_force = mass * U.A(dmax) / (2 * root)
    prefactor_hessian = mass * U.A(dmax) / (4 * U.A(q0) * root)
    max_gradient = U.A(max(norms))
    upper_entries = [(i, j, U.A(v), U.A(norms[i] * norms[j]))
                     for (i, j), v in h_norms.items()]

    def tails(order):
        A = prefactor_force * U.A(S.binomial_tail(Q(-1, 2), order, rho))
        B = prefactor_hessian * U.A(S.binomial_tail(Q(-3, 2), order, rho))
        force = core.upper(A * max_gradient)
        hessian, entry = max(((core.upper(A * v + B * gu_gv), (i, j))
                              for i, j, v, gu_gv in upper_entries), key=lambda pair: pair[0])
        return force, hessian, entry

    def status_at(order, result):
        force, hessian, entry = result
        return {
            "order": order,
            "force_tail_upper_exact": core.pack(force)["upper_exact"],
            "hessian_entry_tail_upper_exact": core.pack(hessian)["upper_exact"],
            "hessian_maximizing_local_indices": entry,
            "force_passes": force <= thresholds[0],
            "hessian_passes": hessian <= thresholds[1],
        }

    orders = {}
    order = 0
    while order <= MAX_ORDER and len(orders) < 2:
        result = tails(order)
        for offset, name in enumerate(("force", "hessian")):
            if name not in orders and result[offset] <= thresholds[offset]:
                orders[name] = order
        order += 1
    selected = max(orders.values()) if len(orders) == 2 else None
    at_28 = status_at(COMPARISON_ORDER, tails(COMPARISON_ORDER))
    at_selected = status_at(selected, tails(selected)) if selected is not None else None
    at_predecessor = status_at(selected - 1, tails(selected - 1)) if selected else None
    require(at_predecessor is None or not (at_predecessor["force_passes"] and
                                           at_predecessor["hessian_passes"]), "Selected order minimal")
    return {
        "face": face, "face_name": ["basal", "apical"][face], "edge": edge,
        "global_free_indices": fields["global_indices"],
        "Rlo_exact": str(bounds[0]), "Rhi_exact": str(bounds[1]),
        "domain": {k: v for k, v in domain.items() if k != "leaves"},
        "depth_saved_certificate_leaf_count": len(domain["leaves"]),
        "rho_policy_goal_achieved": rho <= GOAL,
        "base_point_witness_if_policy_exceeded": (point_spread_witness(target, face, edge, MAX_DEPTH)
                                                  if rho > GOAL else None),
        "K_exact": U.pack_bounds(K_exact), "tension_exact": U.pack_bounds(tension_exact),
        "tension_minus_K_abs_upper_exact": str(dmax),
        "material_measure_mass": core.pack(mass_integral),
        "auxiliary_measure_mass_exact": str(rule_mass_exact), "two_mass_upper": core.pack(mass),
        "physical_G_gradient_norms_exact": list(map(str, norms)),
        "physical_G_hessian_upper_norms_exact": {f"{i},{j}": str(v) for (i, j), v in h_norms.items()},
        "minimum_force_passing_order": orders.get("force"),
        "minimum_H_entry_passing_order": orders.get("hessian"),
        "minimum_joint_passing_order": selected, "at_order_28": at_28,
        "at_minimum_joint_order": at_selected, "at_predecessor_order": at_predecessor,
        "domain_seconds": domain_seconds, "elapsed_seconds": time.perf_counter() - started,
    }


def run():
    U.set_precision(640)
    bindings, row, target, predictor = source_and_inputs()
    nodes = target.node_tube(predictor, row["parameter_cell"], row["base_parameters_exact"],
                             row["selected_corrector_radius_exact"])
    thresholds = U.A(F_LIMIT), U.A(H_LIMIT)
    scope = {
        "domain": "L's original exact binary center plus affine predictor image of unchanged C "
                  "and original selected axis-aligned normal correction cube, enclosed coefficientwise.",
        "no_retube": True, "no_spatial_series_expansion": True, "no_moment_defects": True,
        "no_force_or_H_cache_sums": True, "no_global_matrix_arithmetic": True, "no_root_claim": True,
    }
    protocol = {
        "source_and_input_hashes": bindings,
        "record_id": RECORD, "dose_exact": row["dose_exact"], "original_radius_exact": row["selected_corrector_radius_exact"],
        "parameter_cell_sha256": row["parameter_cell_sha256"], "branch_domain_id": row["branch_domain_id"],
        "precision_bits": 640, "depth_maximum": MAX_DEPTH, "rho_policy_goal_exact": str(GOAL),
        "order_cap": MAX_ORDER, "force_tail_upper_limit_exact": str(F_LIMIT),
        "hessian_entry_tail_upper_limit_exact": str(H_LIMIT), "scope": scope,
    }
    U.write(HERE / "protocol.json", protocol)
    U.write(HERE / "started.json", {"pid": os.getpid(), "utc": U.stamp(),
                                     "protocol_sha256": U.sha(HERE / "protocol.json")})
    rows = []
    failures = []
    for face in (0, 1):
        for edge in range(12):
            try:
                result = check_edge(target, nodes, row, face, edge, thresholds)
            except Exception as error:
                failures.append({"face": face, "edge": edge, "error": repr(error),
                                 "traceback": traceback.format_exc()})
                print(f"Failed face={face} edge={edge}: {error}", flush=True)
                continue
            rows.append(result)
            print(f"face={face} edge={edge} rho={float(Q(result['domain']['rho_exact'])):.8f} "
                  f"min_order={result['minimum_joint_passing_order']} "
                  f"seconds={result['elapsed_seconds']:.3f}", flush=True)
    complete = len(rows) == 24 and all(row["minimum_joint_passing_order"] is not None for row in rows)
    counts = Counter(row["minimum_joint_passing_order"] for row in rows
                     if row["minimum_joint_passing_order"] is not None)
    U.write(HERE / "results.json", {
        "status": "PASS_L_ORIGINAL_24_SURFACE_DOMAINS_AND_TAIL_SELECTORS" if complete else "INCOMPLETE",
        "protocol_sha256": U.sha(HERE / "protocol.json"), "records": rows, "failures": failures,
        "summary": {"orders": dict(sorted(counts.items())),
                    "edges_above_rho_policy": [[row["face"], row["edge"]] for row in rows
                                               if not row["rho_policy_goal_achieved"]],
                    "edges_whose_tail_fails_at_28": [[row["face"], row["edge"]] for row in rows
                                                     if not (row["at_order_28"]["force_passes"] and
                                                             row["at_order_28"]["hessian_passes"])]},
        "scope": scope,
    })
    return 0 if complete else 1


if __name__ == "__main__":
    if any((HERE / filename).exists() for filename in ("started.json", "terminal.json", "results.json", "protocol.json")):
        raise SystemExit("Existing run artifacts are write-once")
    start = time.perf_counter()
    try:
        code = run()
    except Exception:
        code = 1
        U.write(HERE / "failure.json", {"traceback": traceback.format_exc()})
        traceback.print_exc()
    U.write(HERE / "terminal.json", {"exit_code": code, "pid": os.getpid(), "utc": U.stamp(),
                                      "elapsed_seconds": time.perf_counter() - start,
                                      "result_sha256": U.sha(HERE / "results.json")
                                      if (HERE / "results.json").exists() else None})
    raise SystemExit(code)
