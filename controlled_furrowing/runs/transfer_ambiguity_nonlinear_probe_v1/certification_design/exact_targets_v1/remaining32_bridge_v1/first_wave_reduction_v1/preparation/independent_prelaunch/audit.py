"""Independent synthetic audit of the first-wave reducer and its signed kernel."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "PLAN.md").is_file())
sys.path.insert(0, str(ROOT / "scripts"))
import transfer_ambiguity_remaining32_first_wave_reduce as C

SOURCE = "58de9f5e96b43005b0ba83dd0f4d8cd1b35f0690eef58290a75ef5abce692d92"
PROTOCOL = "4b9ecc9794f84e00d62f1fbc090b3a96739af044a2a426194d8fba07967fe487"
PREFLIGHT = "2f1c69917e32d25d22a28577f80dc6a2c8c7fab025d93967e73fcf48cffeb79c"


def exists(doc):
    return doc["strict_contraction_and_selfmap"] is True


def parse(bounds):
    return tuple(map(F, bounds))


def envelope(rows):
    return tuple((min(x), max(x)) for x in zip(*rows, strict=True))


def row_norm(rows):
    return max(sum(map(abs, row)) for row in rows)


def affine(r, v):
    return tuple(sum(a * b for a, b in zip(row, v, strict=True)) for row in r)


def bounds_contain(out, ref):
    return all(a <= c <= d <= b for (a, b), (c, d) in zip(map(parse, out), ref, strict=True))


def component(slot, ids, force, H):
    route = {"kind": "bulk", "slot": slot, "global_indices": ids}
    payload = {"point_force": {"scope": C.S.POINT, "F_I_minus_F_cache_intervals_exact": force},
               "whole_box_hessian": {"scope": C.S.BOX, "radius_exact": str(C.S.RADIUS),
                                    "H_I_minus_H_cache_upper_exact": H}}
    return route, {**route, **payload}


def run():
    a, worker = C.static_phase()
    p, f = (C.read(C.PREP / fn, h) for fn, h in (("protocol.json", PROTOCOL), ("preflight.json", PREFLIGHT)))
    assert p["source_sha256"] == f["source_sha256"] == C.sha(C.SOURCE) == SOURCE
    assert f["protocol_sha256"] == PROTOCOL and C.sha(C.PREP / "reducer_frozen.py") == SOURCE
    assert p["worker_runtime_binding"] == worker and p["worker_target_authority"] == a["targets"]
    assert C.V.yes(f["checks"])
    assert all(C.V.zero(p[key]) and C.V.zero(f[key]) for key in C.ZERO)
    assert not C.ACTUAL.exists()
    K = C.K
    K.ctx.prec = 640
    route0, p0 = component(0, [2, 0], [["2","3"],["-1","1"]],
        [[0,0,"-1","2"],[0,1,"3","4"],[1,1,"2","3"]])
    route1, p1 = component(1, [1, 0], [["0","2"],["-2","0"]],
        [[0,0,"0","1"],[0,1,"-2","-1"],[1,1,"-3","-2"]])
    force, H = K.assemble([p0,p1], [route0,route1], 3)
    reference_F = tuple(map(parse, [("-3","1"),("0","2"),("2","3")]))
    reference_H = {k: parse(v) for k, v in {
        (0,0):("-1","1"), (1,1):("0","1"), (2,2):("-1","2"),
        (0,1):("-2","-1"), (0,2):("3","4"), (1,2):("0","0")}.items()}
    assert tuple(force) == reference_F
    assert all(H.get(k, (F(0), F(0))) == v for k, v in reference_H.items())
    r = ((F(5,7),F(-3,11),F(1,19)), (F(-13,17),F(0),F(2,3)), (F(1,2),F(-1,7),F(-1,13)))
    anchor = ((F(-3,19),F(1,7)), (F(2,5),F(3,7)), (F(3,11),F(4,13)))
    zq = F(9,64)
    reduced = K.contract(r, force, H, anchor, zq, F(1))
    vectors = list(product(*reference_F))
    got_defect = envelope([affine(r, v) for v in vectors])
    got_full = envelope([tuple(x + y for x, y in zip(affine(r, v), u, strict=True))
                         for v in vectors for u in product(*anchor)])
    assert bounds_contain(reduced["signed_RF_defect_rows_exact"], got_defect)
    assert bounds_contain(reduced["signed_RF_I_rows_exact"], got_full)
    H_keys = tuple(reference_H)
    H_vertices = [dict(zip(H_keys, vertex, strict=True)) for vertex in product(*reference_H.values())]
    def RH(m):
        return tuple(tuple(sum(r[i][k] * m[tuple(sorted((k,j)))] for k in range(3))
                     for j in range(3)) for i in range(3))
    rh_norm = max(map(lambda m: row_norm(RH(m)), H_vertices))
    assert F(reduced["z_upper_exact"]) == F(reduced["b_upper_exact"]) + zq
    assert F(reduced["b_upper_exact"]) >= rh_norm
    assert F(reduced["b_upper_exact"]) - rh_norm <= F(1, 2**24), "outward generic balls should retain useful precision"
    rv = C.S.RADIUS
    strict = K.classify([(rv/16, rv/8)], F(1,2), rv)
    tie = K.classify([(rv/2, rv/2)], F(1,2), rv)
    away = K.classify([(-rv*2,-rv*F(7,4))], F(1,2), rv)
    assert exists(strict) and F(strict["posterior_error_upper_exact"]) == rv/4
    assert not exists(tie) and not tie["strict_exclusion_of_this_cube"]
    assert away["strict_exclusion_of_this_cube"] and not exists(away)
    sigma, xi = rv/16, rv/64
    target = {"archived_center_shift_exact":str(sigma),"original_addition_gap_exact":str(xi),
              "observer_goal_exact":str(sigma + rv/4 + rv/1024)}
    premise = {"inert_test_for_same_Sym_R_SPD": True}
    converted = C.transfer(strict, target, premise, F(4))
    assert converted["strict_shift_inclusive_observer_budget_candidate"] is True
    assert converted["independent_observer_transfer_still_required"] is True
    assert F(converted["physical_nodal_error_um_upper_exact"]) == F(67,2)*(sigma+rv/4+xi)
    assert F(converted["candidate_whole_box_H_euclidean_lower_exact"]) == F(1,8)
    assert not C.transfer(strict, target, None, F(4))["candidate_strict_minimum_requires_independent_final_review"]
    assert not C.transfer(strict, {**target, "observer_goal_exact":str(sigma + rv/4)}, premise, F(4))[
        "strict_shift_inclusive_observer_budget_candidate"]
    no_min = C.transfer(away, target, premise, F(4))
    assert no_min["physical_nodal_error_um_upper_exact"] is None and no_min["candidate_whole_box_H_euclidean_lower_exact"] is None
    checks = {
        "source_current_frozen_and_worker_authority_match":True,
        "permuted_sparse_signed_scatter_off_diagonal_not_doubled":True,
        "signed_nondyadic_three_dimensional_RF_encloses_fraction_corner_oracle":True,
        "signed_preconditioned_cached_anchor_added_exactly_once":True,
        "RH_upper_from_symmetric_H_encloses_fraction_corner_norm":True,
        "strict_banach_boundary_and_local_only_exclusion":True,
        "comparison_goal_already_reserves_xi_and_physical_transfer_adds_it_once":True,
        "minimum_requires_same_R_SPD_premise_and_strict_root":True,
        "producer_fixture_hull_guards_repeat":C.V.yes(C.fixtures()),
        "all_physical_inputs_remain_externally_pending":True,
    }
    assert C.V.yes(checks)
    print(C.save(HERE / "exact_arithmetic_fixtures.json", {
        "status": "PASS","review_fixture_source_sha256": C.sha(Path(__file__)),
        "source_sha256":SOURCE, "protocol_sha256":PROTOCOL, "preflight_sha256":PREFLIGHT,
        "checks":checks, "actual_raw_or_hull_packets_opened":0, "actual_spatial_or_global_arithmetic_jobs":0,
        "synthetic_corners": {"force_plus_anchor":len(vectors)*2**len(anchor), "symmetric_H":len(H_vertices)},
    }))


if __name__ == "__main__":
    run()
