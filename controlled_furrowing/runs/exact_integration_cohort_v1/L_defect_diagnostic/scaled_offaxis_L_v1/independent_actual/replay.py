"""Audit completed scaled-L packet metadata and independently reduce exact endpoints."""
import argparse
from datetime import datetime
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from pathlib import Path
import sys
from unittest.mock import patch

from flint import arb

OUT = Path(__file__).resolve().parent
STAGE = OUT.parent
ROOT = STAGE.parents[3]
PREP = STAGE / "preparation"
REDUCTION = STAGE / "reduction_preparation"
COMPONENT_RECEIPT = OUT / "component_and_endpoint_replay.json"
REVIEW = STAGE / "actual_independent_review.json"
FINAL_REVIEW = REDUCTION / "actual_independent_review.json"
ACTUAL_PINS = {
    "actual/terminal.json": "c768fe8cce6978ec2b6865c784dc29484121ed4e284bb84d7358df58ac4a2232",
    "actual/results.json": "338d7ba6806d5d9c3013d9559e0cee4333031545409fd679b799b125b6e72c95",
}
REDUCER_PRELAUNCH = "8de286be02dd50f43eef95f676c8491e3b2abb5e902a18b47b7fbc0817a1c0f4"
PID = 30368
N = 611
ZERO = (Q(0),Q(0))


def record_inputs():
    hashes = {}

    def read(filename, expected=None, parse=True):
        file = (ROOT/filename).resolve()
        assert file.is_relative_to(ROOT)
        path = file.relative_to(ROOT).as_posix()
        data = file.read_bytes()
        h = sha256(data).hexdigest()
        assert expected is None or h == expected, path
        assert hashes.setdefault(path,h) == h, path
        return json.loads(data) if parse else data

    return hashes,read


def canonical(obj):
    return sha256(json.dumps(obj,default=str,separators=(",",":"),sort_keys=True).encode()).hexdigest()


def pack(row):
    lo,hi = (Q(row[key]) for key in ("lower_exact","upper_exact"))
    assert lo <= hi and all(x.denominator & (x.denominator-1) == 0 for x in (lo,hi))
    return lo,hi


def inter(a,b):
    lo,hi = max(a[0],b[0]),min(a[1],b[1])
    assert lo <= hi
    return lo,hi


def maxabs(x):
    return max(abs(x[0]),abs(x[1]))


def matrix_entries(packet,n,key="H_integral_minus_cache_upper"):
    result = {}
    for row in packet[key]:
        i,j = row["i"],row["j"]
        assert type(i) is type(j) is int and (i,j) not in result
        result[i,j] = row
    assert set(result) == set(combinations_with_replacement(range(n),2))
    return result


def scatter(g,H,ids,force,hess):
    assert len(force) == len(ids) == len(set(ids))
    for i,(lo,hi) in zip(ids,force):
        a,b = g[i]
        g[i] = (a+lo,b+hi)
    for (i,j),(lo,hi) in hess.items():
        key = tuple(sorted((ids[i],ids[j])))
        a,b = H.get(key,ZERO)
        H[key] = (a+lo,b+hi)


def reduce_direct(R,g,H):
    Hrow = [Q(0)]*N
    for (i,j),pair in H.items():
        up = maxabs(pair)
        Hrow[i] += up
        if i != j:
            Hrow[j] += up
    rf,rows = [],[]
    norm = Q(0)
    for row in R:
        weights = [Q.from_float(float(v)) for v in row]
        projected = [((w*a,w*b) if w >= 0 else (w*b,w*a)) for w,(a,b) in zip(weights,g) if w]
        rf.append((sum((a for a,b in projected),Q(0)),sum((b for a,b in projected),Q(0))))
        rows.append(sum((abs(w)*v for w,v in zip(weights,Hrow) if w and v),Q(0)))
        norm = max(norm,sum(map(abs,weights)))
    return (rf,rows,norm)


def verify_scaled(A,ctx,nodes,rule,body,case,proof,old,new,validated_pin,original_sha):
    same = ("record_id","dose_exact","parameter_cell_sha256","model_observation_contract_sha256",
            "point_cache_sha256","predictor_sha256","precision_bits","element","chart","global_indices",
            "free_local","active","nodal_hulls_sha256","auxiliary_element_rule_sha256")
    assert all(new[key] == old[key] == proof[key] for key in same)
    assert all(new[k] == v for k,v in validated_pin.items() if k != "prelaunch_sha256")
    assert new["status"] == "SCALED_OFFAXIS_COMPONENT_UNREVIEWED"
    assert new["original_packet_sha256"] == original_sha
    assert new["selector_proof_sha256"] == canonical(proof)
    assert new["degree"] == old["degree"] == proof["unchanged_degree"] == case["degree"]
    assert new["rho_exact"] == old["rho_exact"] == proof["unchanged_rho_exact"] == case["rho_exact"]
    assert proof["body_element_sha256"] == canonical(body)
    assert proof["physical_selector_sha256"] == canonical(case)
    j_lower,j_upper = (Q(body["J_bound"][k]) for k in ("lower_exact","upper_exact"))
    assert 0 < Q(1,10**6) < j_lower <= j_upper and body["positive_J_for_entire_box"] is True
    cache = new["same_Y_cache_consistency"]
    assert cache == old["cache_consistency"] and cache["all64_convex_strain_lines_positive"] is True
    Jcached = pack(cache["minimum_cache_convex_line_J_lower"])[0]
    assert Jcached > Q(1,10**6) and pack(cache["max_coefficientwise_delta_q_upper"])[0] >= 0
    identity,field,_,mass,ids,g = A._geometry(ctx,nodes,rule,body,case,proof["element"])
    assert all(identity[k] == proof[k] for k in identity)
    q,B = A.S.numerators(g,[nodes[j] for j in ids],field["free_local"])
    jet = A.JET.bulk_jacobian_jet(q,B)

    def bound(poly):
        assert all(i >= 0 and j >= 0 for i,j in poly)
        def magnitude(v):
            if isinstance(v,arb):
                assert v.is_finite()
                return max(abs(Q(str(v.lower().fmpq()))),abs(Q(str(v.upper().fmpq()))))
            return abs(v)
        return sum(magnitude(v) for v in poly.values()) / g.Rhat_min

    g_max,h_max = (Q(case[k]) for k in ("max_DJ_norm_upper_exact","max_D2J_norm_upper_exact"))
    a,b = A.M.volumetric_derivative_tail_bounds(case["degree"],Q(case["rho_exact"]),1,5000)
    assert a == Q(case["phi_prime_tail_upper_exact"]) and b == Q(case["phi_second_tail_upper_exact"])
    norms = [min(g_max,bound(v)) for v in jet["gradient"]]
    ft = [mass*a*v for v in norms]
    Ht = {(i,j):mass*(b*norms[i]*norms[j]+a*min(h_max,bound(v)))
          for (i,j),v in jet["hessian_upper"].items()}
    assert ft == list(map(Q,new["force_two_mass_tail_upper_exact"]))
    assert 0 <= min(ft) <= max(ft) == Q(proof["force_tail_max_exact"]) <= Q(proof["force_tail_budget_exact"])
    old_f = list(map(pack,old["force_integral_minus_cache"]))
    new_f = list(map(pack,new["force_integral_minus_cache"]))
    f = [inter(x,y) for x,y in zip(old_f,new_f)]
    n = len(old["global_indices"])
    assert len(old_f) == len(new_f) == len(f) == n == len(ft)
    assert f == list(map(pack,new["force_intersection_with_original"]))
    old_h = {key:pack(row["enclosure"]) for key,row in matrix_entries(old,n).items()}
    new_h = matrix_entries(new,n)
    assert Ht == {key:Q(row["tail_upper_exact"]) for key,row in new_h.items()}
    assert 0 <= min(Ht.values()) <= max(Ht.values()) == Q(proof["H_entry_tail_max_exact"]) <= Q(proof["H_entry_tail_budget_exact"])
    H = {key:inter(old_h[key],pack(row["enclosure"])) for key,row in new_h.items()}
    assert H == {key:pack(row["enclosure"]) for key,row in matrix_entries(new,n,"H_intersection_with_original").items()}
    return f,H,old_f,old_h,{"element":proof["element"],"force_columns":n,"H_entries":len(H),
                           "minimum_same_Y_cache_J_lower_exact":str(Jcached),
                           "F_scaled_subset_of_original":sum(u == v for u,v in zip(f,new_f)),
                           "H_scaled_subset_of_original":sum(H[k] == pack(row["enclosure"]) for k,row in new_h.items())}


def write_once(path,doc):
    with path.open("x",encoding="utf8",newline="\n") as f:
        f.write(json.dumps(doc,indent=2,default=str)+"\n")
    return sha256(path.read_bytes()).hexdigest()


def compose():
    hashes,read = record_inputs()
    terminal,manifest = [read(STAGE/name,h) for name,h in ACTUAL_PINS.items()]
    started = read(STAGE/"actual/started.json")
    assert type(terminal["exit_code"]) is int and terminal["exit_code"] == 0 and terminal["status"] == "FINISHED"
    assert terminal["pid"] == started["pid"] == PID and terminal["results_sha256"] == ACTUAL_PINS["actual/results.json"]
    seconds = (datetime.fromisoformat(terminal["utc"])-datetime.fromisoformat(started["utc"])).total_seconds()
    assert 0 < manifest["seconds"] < seconds+1
    prelaunch = read(REDUCTION/"prelaunch_independent_review.json",REDUCER_PRELAUNCH)
    assert prelaunch["status"] == "PASS" and all(v is True for v in prelaunch["checks"].values())

    sys.path.insert(0,str(REDUCTION))
    sys.path.insert(0,str(ROOT/"scripts"))
    import reduce_scaled_L as V
    import exact_scaled_offaxis_cohort as A
    blocked = AssertionError("actual moment integration or cache sums attempted during replay")
    with patch.object(A.S.Functional,"__init__",side_effect=blocked), \
            patch.object(A.S,"_eval_parts",side_effect=blocked), \
            patch.object(A.U,"uniform_cache",side_effect=blocked), \
            patch.object(A.U.C,"exact_bulk_maps",side_effect=blocked):
        p_hash,scale,selectors,oa,original,pins = V._bound_preparation()
        assert all(started[k] == manifest[k] == value for k,value in pins.items())
        for path,h in p_hash.items():
            read(path,h,False)
        expected_names = scale["actual_expected_files"]
        ids = sorted(set(range(144))-set(V.AXIS))
        assert len(ids) == len(expected_names) == 132
        assert manifest["status"] == "ALL_132_SCALED_OFFAXIS_PACKETS_UNREVIEWED"
        assert manifest["record_id"] == V.RECORD and manifest["unchanged_axis_elements"] == V.AXIS
        assert manifest["original_surface_packets_retained"] == 24 and manifest["global_reduction_run"] is False
        assert [entry["file"] for entry in manifest["files"]] == expected_names
        assert [entry["element"] for entry in manifest["files"]] == ids
        assert len(set(expected_names)) == len(manifest["files"])
        actual = {entry["element"]:entry for entry in manifest["files"]}
        proofs = {entry["element"]:entry for entry in selectors}
        inp,ctx,nodes,rule,body,domains,original_files,offaxis,axis = A._pinned_L_inputs()
        assert offaxis == ids and axis == V.AXIS and scale["input_hashes"] == inp
        assert scale["parameter_cell_sha256"] == ctx.row["parameter_cell_sha256"]
        assert ctx.old_radius == Q(scale["original_cube_radius_exact"])
        original_p = read(f"{V.ORIGINAL}/protocol.json",oa["source_and_manifest_hashes"][f"{V.ORIGINAL}/protocol.json"])
        box = read(original_p["coordinate_box_file"],original_p["coordinate_box_sha256"])
        assert all(canonical(box[k]) == canonical(v) for k,v in V.D.C.exact_node_box(ctx,ctx.old_radius).items())
        assert original_p["row"] == ctx.row
        f0,f1 = [ZERO]*N,[ZERO]*N
        H0,H1 = {},{}
        local = []
        count,hr = [0]*N,[0]*N
        axis_retained,surf_retained = [],[]
        for part,cases in (("bulk",domains["bulk_selectors"]),("surface",domains["surface_cases"])):
            by_path = {item["file"]:item["sha256"] for item in oa["expected_packets"][part]}
            for i,case in enumerate(cases):
                face,edge = divmod(i,12)
                name = f"element_{i:03d}.json" if part == "bulk" else f"face{face}_edge{edge:02d}.json"
                path = f"{V.ORIGINAL}/{part}/{name}"
                old = read(path,by_path[path])
                assert original["source_input_and_component_hashes"][path] == by_path[path]
                idx = V.D._packet_meta(ctx,part,i,case,old,body,original_p)
                n = len(idx)
                if part == "bulk" and i in actual:
                    rec = actual[i]
                    new = read(rec["file"],rec["sha256"])
                    f,H,pf,pH,info = verify_scaled(A,ctx,nodes,rule,body["triangles"][i],case,proofs[i],old,new,pins,by_path[path])
                    local.append(info)
                else:
                    pf = f = list(map(pack,old["force_integral_minus_cache"]))
                    pH = H = {key:pack(entry["enclosure"]) for key,entry in matrix_entries(old,n).items()}
                    if part == "bulk":
                        axis_retained.append(i)
                    else:
                        surf_retained.append((face,edge))
                scatter(f0,H0,idx,pf,pH)
                scatter(f1,H1,idx,f,H)
                for j in idx:
                    count[j] += 1
                    hr[j] += n
        assert [info["element"] for info in local] == ids
        assert axis_retained == V.AXIS and len(surf_retained) == 24
        assert count == original_p["force_incidence"] and hr == original_p["H_row_incidence"]
        assert len(H0) == len(H1) == 6676

    hash_F = canonical(f0)
    hash_H = canonical([(i,j,H0[i,j]) for i,j in sorted(H0)])
    assert hash_F == original["F_endpoint_sha256"] and hash_H == original["H_endpoint_sha256"]
    def nested(a,b):
        return a[0] <= b[0] <= b[1] <= a[1]
    assert all(nested(x,y) for x,y in zip(f0,f1)) and all(nested(H0[k],H1[k]) for k in H0)
    rf,rows,norm = reduce_direct(ctx.R,f1,H1)
    assert norm == Q(original_p["R_infinity_norm_exact"]) == Q(original["R_infinity_norm_exact"])
    RF,b = max(map(maxabs,rf)),max(rows)
    assert RF <= Q(original["RF_uniform_upper_exact"]) and b <= Q(original["b_upper_exact"])
    eta0,z0,r,provider = V.D._old_premises(ctx,original_p,hashes)
    eta,z = eta0+RF,z0+b
    lhs = eta+z*r
    ok = z < 1 and lhs < r
    error = eta/(1-z) if ok else None
    assert (eta,z,lhs,ok,error) == V.D._strict_test(eta0,z0,RF,b,r)
    curvature = V.D._minimum(ctx,hashes,z0,z,r,norm,provider) if ok else None
    result_fields = {
        "status":"CONDITIONAL_UNREVIEWED_SCALED_ORIGINAL_CUBE_SUFFICIENT_TEST_PASSES" if ok else "UNREVIEWED_SCALED_ORIGINAL_CUBE_SUFFICIENT_TEST_UNRESOLVED",
        "record_id":ctx.record_id,"scale_correction_elements":ids,"original_axis_elements":axis_retained,
        "original_surface_packet_keys":surf_retained,"component_count":168,"force_entries":N,"upper_Hessian_entries":len(H1),
        "F_endpoint_sha256":canonical(f1),"H_endpoint_sha256":canonical([(i,j,H1[i,j]) for i,j in sorted(H1)]),
        "RF_endpoint_sha256":canonical(rf),"abs_R_H_row_sha256":canonical(rows),"R_infinity_norm_exact":norm,
        "original_R_force_upper_exact":original["RF_uniform_upper_exact"],"original_b_upper_exact":original["b_upper_exact"],
        "RF_uniform_upper_exact":RF,"b_upper_exact":b,"cached_eta_upper_exact":eta0,"cached_Z_upper_exact":z0,
        "eta_upper_exact":eta,"z_upper_exact":z,"selfmap_lhs_exact":lhs,"normal_radius_exact":r,
        "sufficient_conditions_pass":ok,"conditional_root_error_upper_exact":error,"conditional_curvature":curvature,
        "same_point_cell_and_cube":True,"root_assertion_independently_reviewed":False,
    }
    read(Path(__file__).relative_to(ROOT),parse=False)
    for path,h in tuple(hashes.items()):
        read(path,h,False)
    report = {"status":"PASS_ACTUAL_132_PACKET_METADATA_AND_INDEPENDENT_ENDPOINT_REPLAY",
              "record_id":ctx.record_id,"checks":{
                    "external_132_terminal_exit0_matching_start_PID_manifest_and_freeze":True,
                    "all132_scaled_packet_exact_context_and_bound_selector_metadata":True,
                    "all132_cache_witnesses_identical_to_original_and_positive_J":True,
                    "all132_element_force_and_Hessian_entry_tails_recomputed_with_same_budgets":True,
                    "all_reported_intersections_equal_direct_exact_interval_intersection":True,
                    "all132_intersections_and_original12_axis24_surface_scattered_independently":True,
                    "original_retained_168_endpoints_match_old_independent_reduction_digests":True,
                    "direct_Fraction_preconditioner_aggregation_and_strict_scalar_gate":True,
              },"complete_scaled_elements":len(local),"force_columns_across_scaled":sum(x["force_columns"] for x in local),
              "H_entries_across_scaled":sum(x["H_entries"] for x in local),
              "local_metadata":local,"expected_actual_reducer_fields":result_fields,
              "minimum_abs_RF_lower_exact":str(max((Q(0) if a <= 0 <= b else min(abs(a),abs(b)) for a,b in rf))),
              "zero_containing_RF_coordinates":sum(a <= 0 <= b for a,b in rf),
              "scope":"Full same-box endpoint/provenance and exact-assembly replay. Actual scaled moment and cache evaluation remain shared reviewed premises. This does not compare the new reducer final record until it is externally pinned.",
              "input_hashes":dict(sorted(hashes.items()))}
    receipt_sha = write_once(COMPONENT_RECEIPT,report)
    review = {k:v for k,v in report.items() if k not in ("local_metadata",)}
    review["support_file"] = COMPONENT_RECEIPT.relative_to(ROOT).as_posix()
    review["support_sha256"] = receipt_sha
    review_sha = write_once(REVIEW,review)
    print(json.dumps({"status":report["status"],"review_sha":review_sha,"receipt_sha":receipt_sha,
                      "RF":float(RF),"b":float(b),"eta":float(eta),"z":float(z),
                      "lhs_over_radius":float(lhs/r),"sufficient_conditions_pass":ok,
                      "zero_containing_RF_coordinates":report["zero_containing_RF_coordinates"]}),flush=True)


def finalize(args):
    hashes,read = record_inputs()
    replay = read(COMPONENT_RECEIPT,args.component_sha)
    for path,h in replay["input_hashes"].items():
        read(path,h,False)
    result = read(REDUCTION/"results.json",args.result_sha)
    terminal = read(REDUCTION/"terminal.json",args.terminal_sha)
    admission = read(REDUCTION/"input_admission.json",args.admission_sha)
    started = read(REDUCTION/"started.json")
    assert type(terminal["exit_code"]) is int and terminal["exit_code"] == 0 and terminal["status"] == "FINISHED"
    assert terminal["pid"] == started["pid"] and terminal["results_sha256"] == args.result_sha
    assert started["admission_sha256"] == result["admission_sha256"] == args.admission_sha
    assert admission["local_scaled_payloads_opened"] is False
    assert admission["scaled_manifest_sha256"] == ACTUAL_PINS["actual/results.json"]
    assert admission["scaled_terminal_sha256"] == ACTUAL_PINS["actual/terminal.json"]
    assert admission["review_sha256"] == REDUCER_PRELAUNCH
    actual_manifest = read(STAGE/"actual/results.json",ACTUAL_PINS["actual/results.json"])
    assert admission["expected_scaled_packets"] == actual_manifest["files"]
    for path,h in result["source_inputs_and_packets"].items():
        read(path,h,False)
    expected = replay["expected_actual_reducer_fields"]
    assert all(result[k] == v for k,v in expected.items())
    report = {"status":"PASS_ACTUAL_SCALED_ORIGINAL_CUBE_REDUCTION_REVIEW",
              "record_id":result["record_id"],"result_sha256":args.result_sha,"terminal_sha256":args.terminal_sha,
              "admission_sha256":args.admission_sha,"component_replay_sha256":args.component_sha,
              "checks":{"true_terminal_and_admission_pin_match":True,"all_exact_reducer_fields_match_independent_replay":True,
                        "all_bound_original_scaled_and_shared_premise_digests_current":True},
              "result":expected,"scope":"Exact endpoint and assembly cross-check on the frozen record L original cube. The result status retains the outcome of its sufficient test.",
              "input_hashes":dict(sorted(hashes.items()))}
    review_sha = write_once(FINAL_REVIEW,report)
    print(json.dumps({"status":report["status"],"review_sha256":review_sha,
                      "original_cube_sufficient_conditions_pass":expected["sufficient_conditions_pass"]}),flush=True)


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Replay requires enabled assertions")
    parser = argparse.ArgumentParser()
    actions = parser.add_subparsers(dest="action",required=True)
    actions.add_parser("compose")
    final = actions.add_parser("finalize")
    for pin in ("component","result","terminal","admission"):
        final.add_argument(f"--{pin}-sha",required=True)
    args = parser.parse_args()
    compose() if args.action == "compose" else finalize(args)
