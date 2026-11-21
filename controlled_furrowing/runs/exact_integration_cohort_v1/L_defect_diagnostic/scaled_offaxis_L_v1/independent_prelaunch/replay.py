"""Independent 132-chart preparation replay, without force/H or cache sums."""
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
from unittest.mock import patch

OUT = Path(__file__).resolve().parent
STAGE = OUT.parent
PREP = STAGE / "preparation"
ROOT = STAGE.parents[3]
SOURCE = "scripts/exact_scaled_offaxis_cohort.py"
EXPECTED = {
    "source_sha256": "41a3745cd85c52ca7cc3d11701617e7d2735d6721704e0a6ce86209ddb6678af",
    "protocol_sha256": "1c07ccbebcbed226c3f8eeb29c707ea0d8524c516de74f86f22dd8267db5db61",
    "preflight_sha256": "54b8776d29337d872287f35aa73f57ae6ccd40505cab848ebf9c04fb79341be3",
    "selectors_sha256": "f95c210036c3fbc13e1acb25477c58fa1b785f77e90ec2b93f3b6fee3e60a2be",
}


def q(x):
    return Q.from_float(float(x)) if isinstance(x, float) else Q(x)


def canonical(value):
    return sha256(json.dumps(value, default=str, separators=(",",":"), sort_keys=True).encode()).hexdigest()


def reject(action, wanted):
    try:
        action()
    except ValueError as error:
        assert wanted in str(error), str(error)
    else:
        raise AssertionError("malformed selector or endpoint admitted")


def independent_kinematics(A, g, field, local, bary):
    R,Z = [sum(p[k]*w for p,w in zip(g.corners,bary)) for k in (0,1)]
    ((r0,z0),(r1,z1),(r2,z2)) = g.corners
    a,b,c,d = r1-r0,r2-r0,z1-z0,z2-z0
    det = a*d-b*c
    er,ez = [Q(c-d,det),Q(d,det),Q(-c,det)],[Q(b-a,det),Q(-b,det),Q(a,det)]
    edges = ((0,1),(1,2),(2,0))
    N = [v*(2*v-1) for v in bary]+[4*bary[i]*bary[j] for i,j in edges]
    grad = lambda v:[(4*l-1)*u for l,u in zip(bary,v)]+[4*(bary[i]*v[j]+bary[j]*v[i]) for i,j in edges]
    DR,DZ = grad(er),grad(ez)
    r,z = zip(*local)
    apply = lambda u,v:sum(x*y for x,y in zip(u,v))
    want_q = [apply(DR,r),apply(DZ,r),apply(DR,z),apply(DZ,z),apply(N,r)/R]
    want_B = []
    for k in field["free_local"]:
        j,axis = divmod(k,2)
        want_B.append([DR[j],DZ[j],0,0,N[j]/R] if axis == 0 else [0,0,DR[j],DZ[j],0])
    u,v = (R-g.R_center)/g.R_scale,(Z-g.Z_center)/g.Z_scale
    Rhat = R/g.R_center
    val = lambda p,k:A.T.evaluate(p,u,v)/(Rhat if k == 4 else 1)
    question,B = A.S.numerators(g,local,field["free_local"])
    assert [val(p,k) for k,p in enumerate(question)] == want_q
    assert [[val(p,k) for p in row] for k,row in enumerate(B)] == list(map(list,zip(*want_B)))


def verify_local(A,ctx,nodes,rule,body,selected,proof,old,center,i):
    source = selected[i]
    tri = [int(k) for k in ctx.target.a["triangles"][i]]
    free_index = {int(dof):pos for pos,dof in enumerate(ctx.free)}
    dofs = [2*k+axis for k in tri for axis in (0,1)]
    free_local = [j for j,k in enumerate(dofs) if k in free_index]
    global_indices = [free_index[dofs[j]] for j in free_local]
    assert source["element"] == old["element"] == proof["element"] == i
    assert source["chart"] == old["chart"] == proof["chart"] == "off_axis"
    assert tuple(map(int,proof["local_nodes"])) == tuple(tri)
    assert proof["global_indices"] == old["global_indices"] == source["global_indices"] == global_indices
    assert proof["free_local"] == old["free_local"] == source["free_local"] == free_local
    assert set(proof["fixed_local_positions"]) == set(range(12))-set(free_local)
    assert proof["record_id"] == old["record_id"] == source["record_id"] == ctx.record_id
    assert proof["dose_exact"] == old["dose_exact"] == source["dose_exact"] == str(ctx.dose)
    assert proof["parameter_cell_sha256"] == old["parameter_cell_sha256"] == ctx.row["parameter_cell_sha256"]
    assert proof["model_observation_contract_sha256"] == old["model_observation_contract_sha256"] == ctx.row["model_observation_contract_sha256"]
    assert proof["point_cache_sha256"] == old["point_cache_sha256"] == ctx.row["point_cache_sha256"]
    assert proof["predictor_sha256"] == old["predictor_sha256"] == ctx.row["predictor_sha256"]
    assert proof["precision_bits"] == old["precision_bits"] == A.precision.prec == 640
    assert proof["active"] == old["active"] == str(int(ctx.target.a["active"][i,0]))
    assert proof["nodal_hulls_sha256"] == old["nodal_hulls_sha256"] == source["nodal_hulls_sha256"]
    assert proof["auxiliary_element_rule_sha256"] == old["auxiliary_element_rule_sha256"] == source["auxiliary_element_rule_sha256"]
    assert proof["unchanged_degree"] == old["degree"] == source["degree"]
    assert proof["unchanged_rho_exact"] == old["rho_exact"] == source["rho_exact"]
    assert proof["physical_selector_sha256"] == canonical(source)
    assert proof["body_element_sha256"] == canonical(body[i])
    assert old["status"] == "UNIFORM_COHORT_BULK_COMPONENT_UNREVIEWED"
    expected_weights = [q(w)/ctx.target.energy_scale for w in ctx.target.a["volume_weights"][i]]
    assert rule["volume_weights"][i] == expected_weights and min(expected_weights) > 0
    assert len(expected_weights) == 64

    fresh = A.prove_selector(ctx,nodes,rule,body[i],source,i,center)
    assert canonical(fresh) == canonical(proof)
    n = len(free_local)
    expected_identities = 5*(n+1) + 1+n+n*(n+1)//2
    assert expected_identities == proof["polynomial_identity_count"]

    g = A.S.chart([ctx.target.a["reference"][j] for j in tri[:3]])
    assert g.Rhat_min == min(r for r,z in g.corners)/g.R_center > 0
    assert g.R_scale*g.Z_scale*g.det_scaled == g.det_physical > 0
    mapped = []
    for cached_N, lp in zip(ctx.target.a["shape"],rule["volume_barycentric"]):
        _,v1,v2,e01,e12,e20 = map(q,cached_N)
        la,lb = v1+(e01+e12)/2,v2+(e12+e20)/2
        assert lp == [1-la-lb,la,lb] and min(lp) > 0
        R,Z = [sum(p[k]*w for p,w in zip(g.corners,lp)) for k in (0,1)]
        x, y, radial = (R-g.R_center)/g.R_scale, (Z-g.Z_center)/g.Z_scale, R/g.R_center
        assert abs(x) <= 1 and abs(y) <= 1 and radial >= g.Rhat_min
        assert A.S.chart_point(g,lp) == (x,y,radial)
        mapped.append((x,y,radial))
    assert len(mapped) == proof["invariant_Qstar_node_count"] == 64
    assert canonical(mapped) == proof["scaled_Qstar_coordinates_sha256"]
    independent_kinematics(A,g,proof,[center[j] for j in tri],rule["volume_barycentric"][i%64])
    f_budget,h_budget = Q(source["force_tail_budget_exact"]),Q(source["hessian_tail_budget_exact"])
    assert 0 <= Q(proof["force_tail_max_exact"]) <= f_budget == Q(proof["force_tail_budget_exact"])
    assert 0 <= Q(proof["H_entry_tail_max_exact"]) <= h_budget == Q(proof["H_entry_tail_budget_exact"])
    return expected_identities


def test_adversarial(A,ctx,nodes,rule,body,selected,proofs,old_docs,indices):
    i,j = indices[:2]
    bad = deepcopy(selected[i])
    bad["rho_exact"] = "0"
    reject(lambda:A.prove_selector(ctx,nodes,rule,body[i],bad,i),"metadata changed")
    reject(lambda:A.prove_selector(ctx,nodes,rule,body[i],selected[j],i),"same immutable off-axis selector")
    reject(lambda:A.scaled_bulk_component(ctx,nodes,rule,body[i],selected[i],proofs[i],old_docs[j]),"original off-axis packet binding mismatch")
    pack = lambda a,b:{"lower_exact":str(a),"upper_exact":str(b)}
    reject(lambda:A._endpoint_intersection(pack(1,3),pack(-2,0)),"disjoint")
    assert A._endpoint_intersection(pack(-4,1),pack(-1,2)) == pack(-1,1)
    for mutated in ("scalar_force","scalar_hessian"):
        case = {"max_DJ_norm_upper_exact":"2","max_D2J_norm_upper_exact":"3",
                "phi_prime_tail_upper_exact":"1","phi_second_tail_upper_exact":"1",
                "force_tail_budget_exact":"2","hessian_tail_budget_exact":"7"}
        f,h = A._tail_bounds(case,Q(1),[Q(1),Q(5)],{(0,0):Q(4),(0,1):Q(1),(1,1):Q(8)})
        assert f == [Q(1),Q(2)] and h == {(0,0):Q(4),(0,1):Q(3),(1,1):Q(7)}
        key = "force_tail_budget_exact" if mutated == "scalar_force" else "hessian_tail_budget_exact"
        case[key] = "0"
        reject(lambda:A._tail_bounds(case,Q(1),[Q(1),Q(5)],{(0,0):Q(4),(0,1):Q(1),(1,1):Q(8)}),"original local")


def run():
    if not __debug__:
        raise RuntimeError("Replay requires enabled assertions")
    bindings = {}

    def raw(name, expected=None):
        file = ROOT/name
        assert file.resolve().is_relative_to(ROOT)
        data = file.read_bytes()
        value = sha256(data).hexdigest()
        assert expected is None or value == expected, name
        assert bindings.setdefault(name,value) == value, name
        return data

    prefix = PREP.relative_to(ROOT).as_posix()
    source = raw(SOURCE,EXPECTED["source_sha256"])
    assert source == raw(f"{prefix}/exact_scaled_offaxis_cohort_frozen.py",EXPECTED["source_sha256"])
    proto,pre,prior = [json.loads(raw(f"{prefix}/{key}.json",EXPECTED[f"{key}_sha256"]))
                       for key in ("protocol","preflight","selectors")]
    assert proto["source_sha256"] == pre["source_sha256"] == EXPECTED["source_sha256"]
    assert pre["protocol_sha256"] == EXPECTED["protocol_sha256"]
    assert pre["selectors_sha256"] == proto["selectors_sha256"] == EXPECTED["selectors_sha256"]
    assert pre["status"] == "PASS_132_OFFAXIS_PREPARATION_ONLY"
    assert pre["checks"] and all(v is True for v in pre["checks"].values())
    assert pre["actual_force_or_Hessian_run"] is proto["actual_force_or_Hessian_run"] is False
    for path,h in proto["input_hashes"].items():
        raw(path,h)
    sys.path.insert(0,str(ROOT/"scripts"))
    import exact_scaled_offaxis_cohort as A
    forbidden = AssertionError("prelaunch attempted actual moments, contractions or local cache sites")
    with patch.object(A.S.Functional,"__init__",side_effect=forbidden), \
            patch.object(A.S,"_eval_parts",side_effect=forbidden), \
            patch.object(A.M,"triangle_moment",side_effect=forbidden), \
            patch.object(A.U,"uniform_cache",side_effect=forbidden), \
            patch.object(A.U.C,"exact_bulk_maps",side_effect=forbidden):
        hashes,ctx,nodes,rule,body,domains,originals,indices,axis = A._pinned_L_inputs()
        assert hashes == proto["input_hashes"] and originals == proto["original_bulk_packet_hashes"]
        assert ctx.record_id == proto["record_id"] == "optoGEF2_small_col_L"
        assert ctx.row["parameter_cell_sha256"] == proto["parameter_cell_sha256"]
        assert ctx.old_radius == Q(proto["original_cube_radius_exact"])
        assert ctx.dose == Q(proto["dose_exact"])
        excluded = [*range(6),*range(72,78)]
        assert axis == proto["original_axis_element_indices"] == excluded
        assert indices == proto["offaxis_element_indices"] == [i for i in range(144) if i not in excluded]
        assert [v["element"] for v in prior] == indices and len(prior) == pre["elements"] == 132
        output = proto["actual_job_output"]
        assert proto["actual_expected_files"] == [f"{output}/elements/element_{i:03d}.json" for i in indices]
        old_docs = {}
        for i in range(144):
            path = f"{A.ORIGINAL}/bulk/element_{i:03d}.json"
            old_docs[i] = json.loads(raw(path,originals[path]))
        for i in excluded:
            actual = old_docs[i]
            assert actual["chart"] in ("axis_edge_physical","axis_vertex")
            assert actual["element"] == i and actual["record_id"] == ctx.record_id
        admission = json.loads(raw(A.S.F.REDUCTION+A.RECORD+"_input_admission.json"))
        surface = admission["expected_packets"]["surface"]
        assert len(surface) == 24
        assert [entry["file"] for entry in surface] == [
            f"{A.ORIGINAL}/surface/face{face}_edge{edge:02d}.json" for face in (0,1) for edge in range(12)]
        for k,entry in enumerate(surface):
            prev = json.loads(raw(entry["file"],entry["sha256"]))
            face,edge = divmod(k,12)
            assert (prev["face"],prev["edge"]) == (face,edge)
            assert prev["record_id"] == ctx.record_id and prev["status"] == "UNIFORM_COHORT_SURFACE_COMPONENT_UNREVIEWED"
            assert prev["parameter_cell_sha256"] == ctx.row["parameter_cell_sha256"]
        center = ctx.target.nodes()
        identity_count = sum(verify_local(A,ctx,nodes,rule,body["triangles"],domains["bulk_selectors"],p,
                                          old_docs[i],center,i) for i,p in zip(indices,prior))
        assert identity_count == pre["exact_polynomial_identities"] == 19824
        assert sum(v["invariant_Qstar_node_count"] for v in prior) == 8448
        test_adversarial(A,ctx,nodes,rule,body["triangles"],domains["bulk_selectors"],
                         {p["element"]:p for p in prior},old_docs,indices)
    assert A.S is sys.modules["scaled_element21"] and A.S.U is A.U
    assert A.S.F.ELEMENT == 21 and A.S.F.RECORD == "optoGEF2_small_col_L"
    assert A.precision.prec == 640
    script = Path(__file__).relative_to(ROOT).as_posix()
    raw(script)
    for path,h in tuple(bindings.items()):
        raw(path,h)
    report = {
        **EXPECTED,"status":"PASS",
        "scope":"Independent prelaunch for the 132 off-axis scaled L components on the frozen original cube/C only. Reviewed e21 moment/AD/cache kernels are shared premises. No new integration, cache sums, global reduction or root assertion.",
        "checks":{
            "frozen_current_source_protocol_preflight_complete_selector_set_and_inputs":True,
            "132_off_axis_only_exactly_excludes_0_to_5_and_72_to_77":True,
            "original_144_bulk_and_24_surface_digests_retained":True,
            "all_local_context_indices_nodes_cache_rule_J_and_rho_match_original":True,
            "19824_exact_chart_identities_replayed_and_independent_P2_sample_for_every_element":True,
            "8448_exact_sites_and_original_positive_binary_normalized_weights_replayed":True,
            "per_entry_min_of_independent_derivative_majorants_retains_original_tail_budgets":True,
            "negative_cross_record_selector_packet_rho_tail_and_intersection_fixtures_rejected":True,
            "explicit_sequential_context_no_cross_element_moment_cache_or_global_rebinding":True,
            "actual_moment_contraction_and_cache_entrypoints_blocked_for_replay":True,
        },
        "record_id":ctx.record_id,"offaxis":132,"axis_untouched":12,"surface_untouched":24,
        "exact_polynomial_identities":identity_count,"exact_transformed_Qstar_sites":8448,
        "actual_force_or_Hessian_run":False,"actual_cache_sums_run":False,
        "input_hashes":dict(sorted(bindings.items())),
        "review_source_file":script,"review_source_sha256":bindings[script],
    }
    if "--freeze" in sys.argv:
        receipt = OUT/"replay.json"
        gate = PREP/"prelaunch_independent_review.json"
        assert not receipt.exists() and not gate.exists()
        with receipt.open("x",encoding="utf8",newline="\n") as stream:
            stream.write(json.dumps(report,indent=2)+"\n")
        report["support_file"] = receipt.relative_to(ROOT).as_posix()
        report["support_sha256"] = sha256(receipt.read_bytes()).hexdigest()
        with gate.open("x",encoding="utf8",newline="\n") as stream:
            stream.write(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":report["status"],"checks":len(report["checks"]),"offaxis":132,
                      "identities":identity_count,"Qstar_sites":8448}),flush=True)


if __name__ == "__main__":
    run()
