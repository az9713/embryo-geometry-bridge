"""Verify actual L element-21 signed endpoints and centered mean-value arithmetic."""

import hashlib
import json
import sys
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "scripts"))

import exact_cohort_context as C
import exact_cohort_integrator as E

RECORD = "optoGEF2_small_col_L"
SOURCE = "9a5f93bb3d4a39020e28bee165ba606d1d17f63e44a17ee758959aa05109d329"
PROTOCOL = "4b2f6063653ed89c1e2c0c5d4ba6e74077fae92c62ef09c38f35101f4559ef8b"
PREFLIGHT = "2b9ce8dfa112198c549e513e19ae1c34084e1a9d8bce9e784d689a50db2f4909"
PRELAUNCH = "38d476bfe981a74f8bb60362e96fb3a7508c000eeb1ba53c9b98e22a543497b2"
RESULT = "e1f0a45dc3e3a7554e5fc84b60a2e670a21169d282da8e6de5f21d52e9eecdf0"
FOLDER = ROOT / "runs/exact_integration_cohort_v1/records" / RECORD / "original_cube_v1"
TERMS = ("volumetric_truncation", "quadratic_strain_defect", "finite_same_Y_cache_correction")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_bytes())


def endpoints(obj):
    lo, hi = (Q(obj[key]) for key in ("lower_exact", "upper_exact"))
    assert lo <= hi
    return lo, hi


def add(*intervals):
    return sum(lo for lo, _ in intervals), sum(hi for _, hi in intervals)


def multiply(a, b):
    products = [x*y for x in a for y in b]
    return min(products), max(products)


def contains(outer, inner):
    return outer[0] <= inner[0] <= inner[1] <= outer[1]


def intersection(*intervals):
    lo, hi = max(a[0] for a in intervals), min(a[1] for a in intervals)
    assert lo <= hi
    return lo, hi


def width(a):
    return a[1]-a[0]


def magnitude(a):
    return max(abs(a[0]), abs(a[1]))


def pack(interval):
    return list(map(str, interval))


def write(path, data):
    with path.open("x", encoding="utf8") as f:
        json.dump(data, f, indent=2)
        f.write("\n")


def main():
    checks = {}

    def check(name, condition):
        assert condition, name
        checks[name] = True

    bound_files = {"element21_force.py": SOURCE, "element21_force_frozen.py": SOURCE,
                   "protocol.json": PROTOCOL, "preflight.json": PREFLIGHT,
                   "prelaunch_independent_review.json": PRELAUNCH, "results.json": RESULT}
    check("all_actual_producer_review_and_frozen_source_hashes_current",
          all(sha(OUT/name) == value for name, value in bound_files.items()))
    protocol, preflight, prelaunch, report, started, terminal = (
        load(OUT / name) for name in ("protocol.json", "preflight.json", "prelaunch_independent_review.json",
                                     "results.json", "started.json", "terminal.json"))
    binding = {"source_sha256":SOURCE, "protocol_sha256":PROTOCOL, "preflight_sha256":PREFLIGHT,
               "review_sha256":PRELAUNCH}
    check("actual_pid_run_and_terminal_receipt_bind_exact_result_and_inputs",
          all(report[k] == started[k] == h for k, h in binding.items()) and
          started["pid"] == terminal["pid"] == 7284 and terminal["exit_code"] == 0 and
          terminal["status"] == "FINISHED" and terminal["results_sha256"] == RESULT and
          started["utc"] < terminal["utc"] and report["seconds"] > 0)
    check("prelaunch_PASS_and_nonempty_literal_checks_match_frozen_protocol",
          prelaunch["status"] == "PASS" and prelaunch["checks"] and
          all(value is True for value in prelaunch["checks"].values()) and
          all(prelaunch[k] == v for k, v in binding.items() if k != "review_sha256") and
          preflight["source_sha256"] == SOURCE and preflight["protocol_sha256"] == PROTOCOL)
    check("all_current_prelaunch_and_actual_domain_inputs_hash_bound",
          protocol["input_hashes"] == report["domain_hashes"] and
          all(prelaunch["input_hashes"][name] == actual for name, actual in report["domain_hashes"].items()) and
          all(sha(ROOT / name) == expected for name, expected in prelaunch["input_hashes"].items()))
    check("one_actual_L_element_and_no_Hessian_or_global_root_job",
          all(packet["record_id"] == RECORD and packet["element"] == 21
              for packet in (protocol, prelaunch, report, started)) and
          report["status"] == "LOCAL_FORCE_DIAGNOSTIC_COMPLETE_UNREVIEWED" and
          report["precision_bits"] == protocol["precision_bits"] == 640 and
          report["fresh_Hessian_moments_or_assembly"] is False and report["new_root_test"] is False)

    box = load(FOLDER / "coordinate_box.json")
    domain = load(FOLDER / "domains.json")
    old = load(FOLDER / "bulk/element_021.json")
    case = domain["bulk_selectors"][21]
    E.arb_context.prec = 640
    E.core.A.cache_clear()
    context = C.load_context(RECORD)
    target = context.target
    ids = target.a["triangles"][21].tolist()
    index = {int(g): i for i,g in enumerate(context.free)}
    dofs = [2*node+component for node in ids for component in (0,1)]
    global_indices = [index[i] for i in dofs]
    check("same_active_off_axis_element_and_12_free_columns",
          ids == protocol["element_nodes"] and global_indices == report["global_indices"] ==
          case["global_indices"] == old["global_indices"] and case["free_local"] == list(range(12)) and
          old["chart"] == case["chart"] == protocol["chart"] == "off_axis" and case["active"] == "1")
    check("unchanged_same_cell_dose_order_and_relative_radius",
          context.row["parameter_cell_sha256"] == protocol["same_C_sha256"] == case["parameter_cell_sha256"] and
          context.dose == Q(protocol["same_dose_exact"]) and old["degree"] == case["degree"] == protocol["order"] == 12 and
          Q(old["rho_exact"]) == Q(case["rho_exact"]) == Q(protocol["rho_exact"]))
    centers = target.nodes()
    radii = [Q(box["coordinate_radii_exact"][node][axis]) for node in ids for axis in (0,1)]
    check("exact_binary_proof_center_matches_unshifted_declared_box",
          all(centers[i][k] == Q(box["current_centers_exact"][i][k]) for i in ids for k in (0,1)) and
          all(Q(x)==0 for x in box["shift_exact"]) and all(x>0 for x in radii))
    wide_nodes = E.SUR.common_nodes(target, box)
    wide = {node: [E.SUR.rational_hull(v) for v in wide_nodes[node]] for node in ids}
    check("outward_whole_box_nodes_contain_declared_interval_and_point",
          all(contains(wide[node][k], (centers[node][k]-r,centers[node][k]+r))
              for (node,k),r in zip(((n,k) for n in ids for k in (0,1)), radii)))

    material_mass_plus_rule = Q(case["geometry_two_mass_upper_exact"])
    coefficient_tail = Q(case["phi_prime_tail_upper_exact"])
    tail_budget = Q(case["force_tail_budget_exact"])
    weights = [E.T.rational(v) / target.energy_scale for v in target.a["volume_weights"][21]]
    assert len(weights) == 64 and min(weights)>0
    totals = {}
    force_part_summaries = {}
    for label, nodes in (("point", centers), ("whole_box_direct", wide_nodes)):
        fields = target.element(21, nodes=nodes)
        jet = E.JET.bulk_jacobian_jet(fields["q"], fields["B"])
        norms = [E._norm(col, fields) for col in jet["gradient"]]
        profile = report[label]
        tails = list(map(Q, profile["two_mass_tail_upper_exact"]))
        check(label+"_12_two_mass_tails_exact_by_local_gradient_norm",
              len(tails) == len(norms) == 12 and
              tails == [material_mass_plus_rule*coefficient_tail*x for x in norms] and
              0 < min(tails) <= max(tails) <= tail_budget and
              max(tails) <= Q(case["force_tail_upper_exact"]))
        check(label+"_64_weight_metadata_and_positive_retained_cache_J",
              Q(profile["largest_cached_weight_exact"]) == max(weights) and
              endpoints(profile["cache_min_J_bound"])[0] > Q(1,10**6))
        parts = [[endpoints(obj) for obj in profile[key]] for key in TERMS]
        total = [endpoints(obj) for obj in profile["normalized_force_total"]]
        recomposed = [add(a,b,c,(-t,t)) for a,b,c,t in zip(*parts,tails)]
        check(label+"_12_signed_totals_enclose_exact_endpoints_plus_symmetric_tails",
              len(total) == len(recomposed) == 12 and
              all(contains(a,b) for a,b in zip(total,recomposed)))
        totals[label] = total
        force_part_summaries[label] = {
            **{key: {"max_width_exact":str(max(map(width, rows))),
                     "max_abs_upper_exact":str(max(map(magnitude, rows)))} for key,rows in zip(TERMS,parts)},
            "max_tail_upper_exact":str(max(tails)),
            "cache_J_lower_exact":str(endpoints(profile["cache_min_J_bound"])[0]),
        }

    old_force = [endpoints(v) for v in old["force_integral_minus_cache"]]
    hessians = {(v["i"],v["j"]):endpoints(v["enclosure"]) for v in old["H_integral_minus_cache_upper"]}
    check("old_element_width_tail_and_complete_symmetric_Hess_packet_bound",
          set(hessians) == {(i,j) for i in range(12) for j in range(i,12)} and
          max(map(width,old_force)) == Q(protocol["existing_max_local_width_exact"]) and
          max(endpoints(v)[1] for v in old["force_two_mass_tail_upper"]) ==
          Q(protocol["existing_max_local_force_tail_exact"]))
    state = [endpoints(v) for v in report["whole_box_centered_state_increment"]]
    gain = [endpoints(v) for v in report["whole_box_gain_derivative"]]
    gain_effect = [endpoints(v) for v in report["whole_box_gain_increment"]]
    point, direct = totals["point"], totals["whole_box_direct"]
    centered = [endpoints(v) for v in report["whole_box_force_centered"]]
    max_rows = [sum(r*magnitude(hessians[min(i,j),max(i,j)]) for j,r in enumerate(radii)) for i in range(12)]
    check("12_state_mean_value_increments_cover_frozen_H_times_exact_symmetric_offsets",
          len(state)==12 and all(contains(a,(-b,b)) for a,b in zip(state,max_rows)))
    g0 = context.base["common_gain"]
    gain_offset = tuple(g-g0 for g in context.cell["common_gain"])
    check("12_affine_gain_increments_cover_exact_cell_offset_times_retained_derivative",
          gain_offset[0] <= 0 <= gain_offset[1] and len(gain)==len(gain_effect)==12 and
          all(contains(a,multiply(gain_offset,b)) for a,b in zip(gain_effect,gain)))
    check("12_centered_totals_enclose_exact_point_plus_state_and_gain_increment",
          len(centered)==12 and all(contains(a,add(b,c,d)) for a,b,c,d in zip(centered,point,state,gain_effect)))
    declared_intersection = [tuple(map(Q,v)) for v in report["whole_box_intersection_exact"]]
    check("12_declared_intersections_equal_exact_direct_intersect_centered",
          len(declared_intersection)==12 and declared_intersection ==
          [intersection(a,b) for a,b in zip(direct,centered)])
    check("all_12_centered_intervals_also_nested_in_retained_old_and_new_direct",
          all(contains(a,c) and contains(b,c) and contains(c,p)
              for a,b,c,p in zip(old_force,direct,centered,point)) and
          centered == declared_intersection)
    check("old_and_new_direct_intersect_for_all_columns_without_nesting_assumption",
          all(contains(intersection(a,b),c) for a,b,c in zip(old_force,direct,centered)))
    check("all_fixed_point_signs_certified_while_whole_box_centered_covers_zero",
          all(lo>0 or hi<0 for lo,hi in point) and all(lo<=0<=hi for lo,hi in centered))

    rows = [{"local_column":i,"global_index":g,"point_exact":pack(point[i]),
             "old_box_exact":pack(old_force[i]),"direct_box_exact":pack(direct[i]),
             "centered_box_exact":pack(centered[i]),"supporting_state_increment_exact":pack(state[i]),
             "supporting_gain_increment_exact":pack(gain_effect[i]),
             "width_reduction_factor_vs_old_exact":str(width(old_force[i])/width(centered[i]))}
            for i,g in enumerate(global_indices)]
    maximum_widths = {label:str(max(map(width, values))) for label,values in
                      (("point",point),("old",old_force),("direct",direct),("centered",centered))}
    support = {"checks":checks, "result_sha256":RESULT, "source_sha256":SOURCE,
               "review_source_sha256":sha(Path(__file__)),
               "producer_started_sha256":sha(OUT/"started.json"),"producer_terminal_sha256":sha(OUT/"terminal.json"),
               "frozen_and_prelaunch_hashes":prelaunch["input_hashes"],
               "row_conclusions":rows,"maximum_widths_exact":maximum_widths,
               "old_direct_containment_count":sum(contains(a,b) for a,b in zip(old_force,direct)),
               "maximum_signed_point_component_abs_upper_exact":str(max(map(magnitude,point))),
               "point_strict_positive_count":sum(lo>0 for lo,_ in point),
               "point_strict_negative_count":sum(hi<0 for _,hi in point),
               "whole_centered_contains_zero_count":sum(lo<=0<=hi for lo,hi in centered),
               "largest_state_increment_abs_exact":str(max(map(magnitude,state))),
               "largest_gain_increment_abs_exact":str(max(map(magnitude,gain_effect))),
               "force_part_summaries":force_part_summaries}
    support_path = HERE / "endpoint_replay.json"
    write(support_path,support)
    receipt = {"status":"PASS_ACTUAL_LOCAL_FORCE_DIAGNOSTIC_ENDPOINT_REVIEW", "checks":checks,
               "source_sha256":SOURCE,"protocol_sha256":PROTOCOL,"preflight_sha256":PREFLIGHT,
               "prelaunch_sha256":PRELAUNCH,"results_sha256":RESULT,
               "support_file":support_path.relative_to(ROOT).as_posix(),"support_sha256":sha(support_path),
               "scope":"One L element and one frozen point/box at unchanged C. Independent exact-rational endpoint "
                        "composition, signed sign/width comparisons and MVT enclosure checks. This review relies on "
                        "the prelaunch-reviewed actual moment/cache and gain derivative interval primitives and the "
                        "old local Hessian defect packet. The context, low-degree field jets and polynomial norm "
                        "primitives are shared. No force integral, Hessian moment, cache weighted derivative or "
                        "global assembly/root computation was repeated."}
    out = OUT / "actual_independent_review.json"
    write(out,receipt)
    print(json.dumps({"status":receipt["status"],"checks":len(checks),"result_sha256":RESULT,
                      "receipt_sha256":sha(out),"support_sha256":sha(support_path)},indent=2))


if __name__ == "__main__":
    main()
