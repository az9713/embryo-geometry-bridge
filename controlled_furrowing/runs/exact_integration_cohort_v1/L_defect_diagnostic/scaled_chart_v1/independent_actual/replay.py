"""Exact endpoint/provenance replay of the completed scaled L e21 component."""
from datetime import datetime
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from pathlib import Path
import sys
from unittest.mock import patch

from flint import arb, ctx, fmpq

OUT = Path(__file__).resolve().parent
JOB = OUT.parent
ROOT = JOB.parents[3]
GATE = JOB / "actual_independent_review.json"
PINS = {
    "scaled_element21.py": "a173b3759fb66fe66e5f3e442c1bcce9a5a185e67515649a94499f82d1d67f60",
    "protocol.json": "76e5086c374fe86f3659ef1d1b099c4a733211245fbb1fa70251b01462c3fc5c",
    "preflight.json": "a35814e39106e14514f2afa6134996781d0e08f5036bcd6e22b1aafd5587a505",
    "prelaunch_independent_review.json": "fe2b0e6f6b4fe32a380f045b8ca6467f0da1dd7960af8c77b28413b9e88403d8",
    "started.json": "91ddacd59f46e08047fb223495ea12c3a06c764ae50f0ad238be9f1a45af589c",
    "terminal.json": "3d99e5a292b47d919e00e815b41d596fcb5ae5dcb7c245704ea98b7813c47ab9",
    "results.json": "f537a5f635d6e1f0d124bb9e52e135ecf3d88cc854559e9a27725e80cfe9683a",
}


def ball_bounds(ball):
    assert ball.is_finite()
    return (Q(str(ball.lower().fmpq())), Q(str(ball.upper().fmpq())))


def interval(row):
    lo, hi = (Q(row[key]) for key in ("lower_exact", "upper_exact"))
    assert lo <= hi
    assert all(v.denominator & (v.denominator - 1) == 0 for v in (lo, hi))
    return (lo, hi)


def intersection(a,b):
    lo,hi = max(a[0],b[0]),min(a[1],b[1])
    assert lo <= hi
    return (lo,hi)


def rows_by_column(rows, n):
    out = {}
    for row in rows:
        key = row["i"],row["j"]
        assert all(type(x) is int for x in key)
        assert key not in out and 0 <= key[0] <= key[1] < n
        out[key] = row
    assert set(out) == set(combinations_with_replacement(range(n), 2))
    return out


def metrics(f,h):
    absrows = [Q(0)] * len(f)
    for (i,j),bounds in h.items():
        mag = max(map(abs,bounds))
        absrows[i] += mag
        if i != j:
            absrows[j] += mag
    return {
        "max_F_width":max(v-u for u,v in f),
        "max_H_entry_width":max(v-u for u,v in h.values()),
        "max_abs_H_row_bound":max(absrows),
        "max_F_magnitude_bound":max(max(map(abs,bounds)) for bounds in f),
    }, absrows


def endpoint_replay(old,new):
    n = 12
    f_old = list(map(interval,old["force_integral_minus_cache"]))
    f_new = list(map(interval,new["force_integral_minus_cache"]))
    f_recorded = new["force_intersections_with_original_exact"]
    assert len(f_old) == len(f_new) == len(f_recorded) == n
    f_inter = [intersection(a,b) for a,b in zip(f_old,f_new)]
    assert f_inter == [tuple(map(Q,x)) for x in f_recorded]
    h_old = {key:interval(row["enclosure"]) for key,row in rows_by_column(old["H_integral_minus_cache_upper"],n).items()}
    new_indexed = rows_by_column(new["H_integral_minus_cache_upper"],n)
    h_new = {key:interval(row["enclosure"]) for key,row in new_indexed.items()}
    h_recorded = new["H_intersections_in_packet_order_exact"]
    assert len(h_recorded) == len(h_old) == len(h_new) == 78
    h_inter = {key:intersection(h_old[key],h_new[key]) for key in h_new}
    assert list(h_inter.values()) == [tuple(map(Q,x)) for x in h_recorded]

    old_summary,old_bounds = metrics(f_old,h_old)
    new_summary,new_bounds = metrics(f_new,h_new)
    joint_summary,joint_bounds = metrics(f_inter,h_inter)
    assert all(a <= b and a <= c for a,b,c in zip(joint_bounds,old_bounds,new_bounds))
    assert all(joint_summary[key] <= min(old_summary[key],new_summary[key]) for key in joint_summary)
    assert all(0 < new_summary[key] < old_summary[key] for key in joint_summary)
    return new_indexed,{
        "intervals_with_scaled_subset_of_original":{"force":sum(a == b for a,b in zip(f_inter,f_new)),
                                                   "Hessian":sum(h_inter[key] == h_new[key] for key in h_new)},
        "old":old_summary,"scaled":new_summary,"intersection":joint_summary,
        "old_over_scaled":{key:old_summary[key]/new_summary[key] for key in old_summary},
        "abs_H_row_bounds":{"old":old_bounds,"scaled":new_bounds,"intersection":joint_bounds},
    }


def analytical_mass(g,rule,c):
    old_precision = ctx.prec
    try:
        ctx.prec = 1024
        to_ball = lambda r:arb(fmpq(r.numerator,r.denominator))
        ((R0,Z0),(R1,Z1),(R2,Z2)) = g.corners
        area = abs((R1-R0)*(Z2-Z0)-(R2-R0)*(Z1-Z0))/2
        weighted_area = area*(R0+R1+R2)/3
        material = 2*arb.pi()*to_ball(weighted_area/c.target.energy_scale)
        aux = sum(rule["volume_weights"][21],Q(0))
        assert weighted_area > 0 and aux > 0
        return ball_bounds(material+to_ball(aux))
    finally:
        ctx.prec = old_precision


def support_keys(S,g,jet,q,B,drive,order):
    # Reconstruct the requested monomial keys without evaluating any moment.
    e = S.T.sub(jet["value"],g.Rhat)
    a,b = S.series_numerators(e,g.Rhat,order)
    grad = jet["gradient"]
    hessian = jet["hessian_upper"]
    requests = {(u+i,v+j,order+1) for factor in [*grad,*hessian.values()]
                for u,v in factor for i,j in a}
    for i,j in hessian:
        product = S.T.mul(grad[i],grad[j])
        requests.update((u+k,v+l,order+2) for u,v in product for k,l in b)
    d = (S.A(1),S.A(1)+drive,S.A(1),S.A(1)+drive)
    columns = list(zip(*B))
    for left,right in [(q,column) for column in columns] + [(columns[i],columns[j]) for i,j in hessian]:
        poly = {}
        for scale,p,s in zip(d,left,right):
            poly = S.T.add(poly,S.T.scale(S.T.mul(p,s),scale))
        requests.update((i,j,0) for i,j in poly)
        requests.update((i,j,2) for i,j in S.T.mul(left[4],right[4]))
    return requests


def run():
    if not __debug__:
        raise RuntimeError("Replay requires enabled assertions")
    hashes = {}

    def raw(path,pin=None):
        file = ROOT/path
        assert file.resolve().is_relative_to(ROOT)
        data = file.read_bytes()
        sha = sha256(data).hexdigest()
        assert pin is None or pin == sha, path
        assert hashes.setdefault(path,sha) == sha, path
        return data

    rel = JOB.relative_to(ROOT).as_posix()
    current = raw(f"{rel}/scaled_element21.py",PINS["scaled_element21.py"])
    assert current == raw(f"{rel}/scaled_element21_frozen.py",PINS["scaled_element21.py"])
    proto,pre,gate,started,term,result = [json.loads(raw(f"{rel}/{name}",PINS[name])) for name in (
        "protocol.json","preflight.json","prelaunch_independent_review.json","started.json","terminal.json","results.json")]
    assert term["status"] == "FINISHED" and type(term["exit_code"]) is int and term["exit_code"] == 0
    assert type(term["pid"]) is int and term["pid"] == started["pid"] > 0
    assert term["results_sha256"] == PINS["results.json"]
    assert 0 < result["seconds"] <= (datetime.fromisoformat(term["utc"]) - datetime.fromisoformat(started["utc"])).total_seconds()
    assert pre["status"] == "PASS_SCALED_CHART_PREPARATION_ONLY" and gate["status"] == "PASS"
    assert all(d["checks"] and all(x is True for x in d["checks"].values()) for d in (pre,gate))
    links = {"source_sha256":PINS["scaled_element21.py"],"protocol_sha256":PINS["protocol.json"],
             "preflight_sha256":PINS["preflight.json"]}
    assert proto["source_sha256"] == pre["source_sha256"] == links["source_sha256"]
    assert pre["protocol_sha256"] == links["protocol_sha256"]
    for key,h in links.items():
        assert result[key] == gate[key] == started[key] == h
    assert result["prelaunch_sha256"] == started["prelaunch_sha256"] == PINS["prelaunch_independent_review.json"]
    raw(gate["review_receipt_file"],gate["review_receipt_sha256"])
    for obj in (proto,gate,result):
        for path,h in obj["input_hashes"].items():
            raw(path,h)
    assert result["input_hashes"] == proto["input_hashes"]

    sys.path.insert(0,str(JOB))
    import scaled_element21 as S
    forbidden = AssertionError("actual integration or cached quadrature entered during endpoint replay")
    with patch.object(S,"material_moment",side_effect=forbidden), \
            patch.object(S,"_eval_parts",side_effect=forbidden), \
            patch.object(S.U,"uniform_cache",side_effect=forbidden), \
            patch.object(S.U.C,"exact_bulk_maps",side_effect=forbidden):
        inputs,c,case,old,rule,ids,nodes,_,field,drive,g,ft,ht = S._certificate()
    assert inputs == proto["input_hashes"]
    assert result["record_id"] == proto["record_id"] == old["record_id"] == c.record_id == "optoGEF2_small_col_L"
    assert result["element"] == old["element"] == proto["element"] == 21
    assert result["degree"] == old["degree"] == proto["unchanged_order"] == 12
    assert result["rho_exact"] == old["rho_exact"] == proto["unchanged_rho_exact"] == case["rho_exact"]
    assert result["precision_bits"] == old["precision_bits"] == proto["precision_bits"] == ctx.prec == 640
    assert result["global_indices"] == old["global_indices"] == proto["global_indices"] == field["global_indices"]
    assert len(set(result["global_indices"])) == 12
    assert 0 < g.Rhat_min <= 1 and S.D._digest(result["coordinate_geometry"]) == S.D._digest(S._metadata(g))
    assert result["status"] == "SCALED_CHART_LOCAL_FORCE_H_COMPONENT_UNREVIEWED"
    assert result["new_global_assembly_or_root_test"] is False
    new_h, improvement = endpoint_replay(old,result)

    q,B = S.numerators(g,[nodes[j] for j in ids],field["free_local"])
    jet = S.JET.bulk_jacobian_jet(q,B)
    assert all(i >= 0 and j >= 0 for poly in [*jet["gradient"],*jet["hessian_upper"].values()] for i,j in poly)
    coeff_norm = lambda poly:sum((max(map(abs,ball_bounds(v))) if isinstance(v,arb) else abs(v) for v in poly.values()),Q(0))/g.Rhat_min
    gn = list(map(coeff_norm,jet["gradient"]))
    hnorm = {key:coeff_norm(poly) for key,poly in jet["hessian_upper"].items()}
    t1,t2 = S.M.volumetric_derivative_tail_bounds(12,Q(case["rho_exact"]),1,5000)
    geom_both = Q(case["geometry_two_mass_upper_exact"])
    theoretical = analytical_mass(g,rule,c)
    actual_mass = interval(result["material_plus_rule_mass_upper"])
    assert 0 < theoretical[0] <= theoretical[1] <= geom_both and theoretical[1] <= actual_mass[1]
    assert theoretical[0] <= actual_mass[0]
    # Both values enclose the same positive mass. Allow distinct associations.
    assert abs(actual_mass[1]-geom_both) < Q(1,2**600)
    assert all(t >= 0 for t in (t1,t2))
    exact_ft = [geom_both*t1*bound for bound in gn]
    exact_ht = {(i,j):geom_both*(t2*gn[i]*gn[j]+t1*bound) for (i,j),bound in hnorm.items()}
    assert exact_ft == list(map(Q,result["force_tail_upper_exact"])) == ft
    assert exact_ht == {key:Q(entry["tail_upper_exact"]) for key,entry in new_h.items()} == ht
    assert max(exact_ft) <= Q(proto["unchanged_force_tail_budget_exact"])
    assert max(exact_ht.values()) <= Q(proto["unchanged_H_entry_tail_budget_exact"])

    cache = result["positive_same_Y_cache"]
    assert cache == old["cache_consistency"] and cache["all64_convex_strain_lines_positive"] is True
    cache_J = interval(cache["minimum_cache_convex_line_J_lower"])
    assert cache_J[0] > Q(1,10**6)
    assert interval(cache["max_coefficientwise_delta_q_upper"])[0] >= 0
    keys = support_keys(S,g,jet,q,B,drive,12)
    assert len(keys) == result["signed_moment_defects"]
    assert {d for i,j,d in keys} == {0,2,13,14}
    crossing = arb(0,S.A(Q(1,16)))
    assert S.T.constant(crossing) and not S.T.constant(arb(0))
    first,second = S.series_numerators({(1,0):crossing},S.T.ONE,4)
    assert (1,0) in first and (1,0) in second and not first[1,0].is_zero() and not second[1,0].is_zero()

    own = Path(__file__).relative_to(ROOT).as_posix()
    raw(own)
    for path,h in tuple(hashes.items()):
        raw(path,h)
    report = {
        **links,"status":"PASS_ACTUAL_SCALED_E21_ENDPOINT_REVIEW",
        "results_sha256":PINS["results.json"],"terminal_sha256":PINS["terminal.json"],
        "started_sha256":PINS["started.json"],"prelaunch_sha256":PINS["prelaunch_independent_review.json"],
        "checks":{
            "pinned_true_terminal_exit_zero_and_same_PID_time_and_all_inputs":True,
            "same_original_L_element_chart_projection_box_C_degree_rho_and_precision":True,
            "all12_force_and78_Hessian_exact_intersections_match_both_original_enclosures":True,
            "all12_and78_local_tails_recomputed_and_fit_original_budgets":True,
            "positive_exact_two_measure_mass_independently_bounded":True,
            "same_64_site_cache_witness_and_positive_J_analytic_cutoff":True,
            "exact_zero_only_synthetic_and_all_actual_moment_keys_reconstructed_without_evaluation":True,
            "local_widths_and_symmetric_absolute_H_rows_recomputed":True,
        },
        "record_id":result["record_id"],"element":21,"force_entries":12,"H_upper_entries":78,
        "signed_moment_keys_independently_counted":len(keys),"cache_minimum_J_lower_exact":str(cache_J[0]),
        "local_only_metrics_exact":improvement,
        "scope":"Independent exact intersection, tail, provenance and local endpoint-width replay. Actual moment contractions, AD jet and same-Y cache routines remain shared prelaunch premises; no second whole integration or cache sum. No global L reassembly, new exact root, or claim about any other element.",
        "actual_integrals_replayed":False,"cache_sums_replayed":False,
        "input_hashes":dict(sorted(hashes.items())),"review_source_file":own,"review_source_sha256":hashes[own],
    }
    if "--freeze" in sys.argv:
        receipt = OUT/"endpoint_replay.json"
        assert not receipt.exists() and not GATE.exists()
        with receipt.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(report,indent=2,default=str)+"\n")
        report["support_file"] = receipt.relative_to(ROOT).as_posix()
        report["support_sha256"] = sha256(receipt.read_bytes()).hexdigest()
        with GATE.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(report,indent=2,default=str)+"\n")
    summary = {kind:{key:float(v) for key,v in improvement[kind].items()}
               for kind in ("old","scaled","intersection","old_over_scaled")}
    print(json.dumps({"status":report["status"],"checks":len(report["checks"]),"moment_keys":len(keys),
                      "metrics":summary},allow_nan=False),flush=True)


if __name__ == "__main__":
    run()
