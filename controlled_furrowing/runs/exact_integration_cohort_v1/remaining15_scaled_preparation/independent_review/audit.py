"""Independent module review. No record domain, cache sum, or moment calculation."""

import ast
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import importlib
import json
from pathlib import Path
import sys
from unittest.mock import Mock, patch

import numpy as np
from flint import arb, ctx

OUT = Path(__file__).resolve().parent
BATCH = OUT.parent
ROOT = OUT.parents[3]
BASE = BATCH.parent.relative_to(ROOT).as_posix()
SOURCE = "scripts/exact_scaled_cohort_bulk.py"
SFILE = BASE + "/L_defect_diagnostic/scaled_chart_v1/scaled_element21.py"
TRANSITIVE = BASE + "/L_defect_diagnostic/element21_force.py"
TOYFILE = BASE + "/integrator_preparation/preflight_integrator.py"
MANIFEST = BASE + "/remaining15_preparation/manifest.json"
DOMAIN_BATCH_REVIEW = BASE + "/remaining15_preparation/independent_prelaunch/batch_independent_review.json"
PINS = {
    SOURCE: "6639bb2fae3c3892ec44115ff90e27242586744a2900daa16177cc5a091f1526",
    f"{BASE}/remaining15_scaled_preparation/protocol.json": "5103212ad5870d370a9fa51000c9bfcdd1496f99652b065622c19c5a262fb656",
    f"{BASE}/remaining15_scaled_preparation/preflight.json": "8cd2dcd90a06b029a4b454b036889943e51649c160a8a0bf762fb8fb2ae37641",
    SFILE: "a173b3759fb66fe66e5f3e442c1bcce9a5a185e67515649a94499f82d1d67f60",
    TRANSITIVE: "9a5f93bb3d4a39020e28bee165ba606d1d17f63e44a17ee758959aa05109d329",
    MANIFEST: "93a2f5858a450ca9a8db9f9bc944b918e2a7446abba614f97de5e47fe434495d",
    DOMAIN_BATCH_REVIEW: "34bc25ad801e579ff49d12b992b1dff53a1f2d19d7796e060d08d23ef4a59e26",
}


def require(predicate, message):
    if not predicate:
        raise ValueError(message)


def hash_bytes(raw):
    return sha256(raw).hexdigest()


def bound(path, pin, parse=False):
    p = (ROOT / path).resolve()
    require(p.is_relative_to(ROOT), f"unexpected source path {path}")
    raw = p.read_bytes()
    require(hash_bytes(raw) == pin, f"changed source or input {path}")
    return json.loads(raw) if parse else raw


def rel(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def all_true(v):
    return isinstance(v, dict) and bool(v) and all(x is True for x in v.values())


def rejects(call):
    try:
        call()
    except ValueError:
        return True
    return False


def structure(A):
    tree = ast.parse((ROOT/SOURCE).read_text(encoding="utf8"))
    top_methods = {n.name for n in tree.body if isinstance(n, ast.FunctionDef)}
    require(top_methods == {"_setup","_tail_bounds","_proof","prepare_component","_scaled_component","bulk_component"},
            "unexpected entry point")
    for node in ast.walk(tree):
        require(not isinstance(node, (ast.Global, ast.Nonlocal)), "implicit adapter rebinding")
        if isinstance(node, ast.Call):
            name = node.func.attr if isinstance(node.func, ast.Attribute) else getattr(node.func, "id", "")
            require(name not in ("open","read_text","read_bytes","load_context","_geometry","_certificate","_asset","prepare","run"),
                    "unexpected file or pilot execution dependency")
        if isinstance(node, ast.Attribute) and isinstance(node.ctx, ast.Store):
            raise ValueError("module/attribute assignment")
    require(rel(A.__file__) == SOURCE and rel(A.S.__file__) == SFILE and rel(A.S.F.__file__) == TRANSITIVE,
            "imported source path")
    require(all(getattr(A, name) is getattr(A.S, name) for name in ("D","E","U","M","T","JET","core","A")),
            "pilot transitive numerical aliases agree with reviewed shared engine")
    require(hasattr(A.S, "material_moment") and hasattr(A.S, "rational_identities"), "reviewed shared scaled primitives")


def profiles(A, C, p, manifest, domain_batch):
    require(domain_batch["status"] == "PASS_REMAINING15_DOMAIN_PRELAUNCH_ONLY", "prepared domain-only review")
    ids = [f"optoGEF2_small_col_{x}" for x in "CDEFGHIJKMNOPQR"]
    profiles = p["record_profiles"]
    require([r["record_id"] for r in profiles] == [r["record_id"] for r in manifest["records"]]
            == [r["record_id"] for r in domain_batch["approved_gates"]] == ids, "selected 15 records")
    old = C.load_context("optoGEF2_small_col_L")
    baseline = {k: (v.dtype.str,list(v.shape), hash_bytes(v.tobytes())) for k,v in old.target.a.items()
                if k not in ("u","preconditioner")}
    require("u" not in baseline and "preconditioner" not in baseline, "comparison is static geometry/rule only")
    report = {}
    seen_anchors = set()
    checks = {}
    chart_cache = {}
    for profile, row in zip(profiles, manifest["records"]):
        record = profile["record_id"]
        own = C.load_context(record)
        require(profile["selected_context"] == row["selected_context"] == own.row
                and Q(profile["exact_dose"]) == own.dose and Q(profile["exact_drive"]) == own.dose*own.base["common_gain"]
                and Q(profile["old_radius_exact"]) == own.old_radius
                and Q(profile["existing_original_force_tail_budget_exact"]) == Q(row["force_tail_budget_exact"])
                and Q(profile["existing_original_H_entry_tail_budget_exact"]) == Q(row["H_entry_tail_budget_exact"]),
                "preserved record context and exact budgets")
        require({k: hash_bytes(getattr(own,k).tobytes()) for k in ("u","R","S")} == profile["own_u_R_S_binary_sha256"],
                "own binary centers and arrays")
        seen_anchors.add(profile["own_u_R_S_binary_sha256"]["u"])
        cache = own.target.a
        static = {k:(v.dtype.str,list(v.shape),hash_bytes(v.tobytes())) for k,v in cache.items()
                  if k not in ("u","preconditioner")}
        require(static == baseline, "same exact mesh and cached rule bytes as static witness")
        rule = own.target.auxiliary_rule()
        route = {"off_axis": [], "axis_vertex": [], "axis_edge_physical":[]}
        for j, elm in enumerate(cache["triangles"]):
            coords = tuple(tuple(Q.from_float(float(x)) for x in cache["reference"][i]) for i in elm[:3])
            radios = tuple(c[0] for c in coords)
            require(min(radios) >= 0 and max(radios) > 0, "reference triangle radial geometry")
            nz = radios.count(Q(0))
            name = ("off_axis","axis_vertex","axis_edge_physical")[nz]
            route[name].append(j)
            weights = rule["volume_weights"][j]
            require(len(weights) == 64 and min(weights) > 0
                    and all(w == Q.from_float(float(b))/own.target.energy_scale
                            for w,b in zip(weights, cache["volume_weights"][j])),
                    "same per-element normalized 64 auxiliary weights")
            if name == "off_axis":
                if coords not in chart_cache:
                    geom = A.S.chart(coords)
                    dr1,dz1 = [coords[1][k] - coords[0][k] for k in (0,1)]
                    dr2,dz2 = [coords[2][k] - coords[0][k] for k in (0,1)]
                    require(geom.det_physical == dr1*dz2-dr2*dz1 > 0
                            and geom.Rhat_min == min(radios)/geom.R_center > 0
                            and geom.R_scale*geom.Z_scale*geom.det_scaled == geom.det_physical,
                            "exact radius positivity and affine measure")
                    rcoords = []
                    for lam in rule["volume_barycentric"]:
                        R, Z = [sum(a*c[k] for a,c in zip(lam,coords)) for k in (0,1)]
                        x,y,r = A.S.chart_point(geom, lam)
                        require(R == geom.R_center+geom.R_scale*x and Z == geom.Z_center+geom.Z_scale*y
                                and r == R/geom.R_center > 0 and abs(x) <= 1 and abs(y) <= 1,
                                "unchanged physical auxiliary node")
                        rcoords.append((x,y,r))
                    chart_cache[coords] = rcoords
                require(len(chart_cache[coords]) == len(rule["volume_barycentric"]) == 64, "all Qstar nodes")
        require(route == profile["choice_of_chart"] and route["axis_vertex"] == list(range(6))
                and route["axis_edge_physical"] == list(range(72,78)) and len(route["off_axis"]) == 132,
                "static 132/6/6 route dispatch")
        report[record] = {k:len(v) for k,v in route.items()}
    require(len(seen_anchors) == 15 and len(chart_cache) == 132, "distinct anchors, 132 shared charts")
    checks["replayed_all_2160_static_routes_and_normalized_weights"] = True
    checks["all_8448_distinct_static_offaxis_sites_preserve_physical_coordinates"] = True
    return report, checks


def toy_proofs(A, toy):
    E = A.E
    opaque = AssertionError("actual cache or moment calculation during prep review")
    checks = {}
    cases = [
        ("off_axis", [(2,-1),(4,0),(3,2)], 1),
        ("axis_vertex", [(0,0),(2,0),(2,3)], 0),
        ("axis_edge_physical", [(0,0),(2,0),(0,3)], 1),
    ]
    for name, vertices, active in cases:
        context, nodes, rule, body = toy(vertices, active)
        # Remove a non-axis axial DOF to exercise a strict projected set on every chart.
        free = [i for i in context.target.index if i != 3]
        context.target.index = {i: pos for pos,i in enumerate(free)}
        context.target.a["free"] = np.asarray(free, dtype=np.int64)
        context.target.a["u"] = np.zeros(len(free))
        selector = E.prepare_bulk_selector(context, nodes, rule, body, 0, Q(1,10**11), Q(1,10**8))
        proof = A.prepare_component(context, nodes, rule, body, selector, 0)
        check_route = "scaled_positive_R_affine_material" if name == "off_axis" else name
        require(proof["chart"] == name and proof["integration_chart"] == check_route
                and proof["free_local"] == free and proof["global_indices"] == list(range(len(free)))
                and proof["unchanged_degree"] == selector["degree"] and proof["unchanged_rho_exact"] == selector["rho_exact"],
                f"projected route and selector for {name}")
        if name == "off_axis":
            n = len(free)
            require(proof["scaled_polynomial_identity_count"] == 5*(1+n)+1+n+n*(n+1)//2,
                    "q, B, J, DJ and all mixed D2J symbolic checks")
        identity = E._bulk_setup(context, nodes, rule, 0)[0]
        with (patch.object(E,"bulk_component", return_value=identity) as axis,
              patch.object(A,"_scaled_component", return_value=identity) as scaled):
            result = A.bulk_component(context,nodes,rule,body,selector,0,proof)
            require(result["integration_chart"] == check_route and result["record_component_proof_sha256"] == A.D._digest(proof),
                    "packet route proof digest")
            require((axis.call_count,scaled.call_count) == ((0,1) if name == "off_axis" else (1,0)),
                    "only intended route called")
            actual = scaled.call_args.args if name == "off_axis" else axis.call_args.args
            require(all(a is b for a,b in zip(actual,(context,nodes,rule))), "same explicit object arguments")
            foreign = deepcopy(context)
            foreign.record_id = foreign.row["record_id"] = "synthetic_other"
            bad_rule = {**rule, "volume_weights": [[2*rule["volume_weights"][0][0],*rule["volume_weights"][0][1:]]]}
            bad_nodes = deepcopy(nodes)
            bad_nodes[0][1] += Q(1,10**4)
            require(all(rejects(c) for c in (
                lambda:A.prepare_component(foreign,nodes,rule,body,selector,0),
                lambda:A.prepare_component(context,nodes,bad_rule,body,selector,0),
                lambda:A.prepare_component(context,bad_nodes,rule,body,selector,0),
                lambda:A.prepare_component(context,nodes,rule,body,{**selector,"degree":selector["degree"]+1},0),
                lambda:A.prepare_component(context,nodes,rule,{**body,"positive_J_for_entire_box":False},selector,0),
                lambda:A.bulk_component(context,nodes,rule,{**body,"perturb_proof":1},selector,0,proof),
                lambda:A.bulk_component(context,nodes,rule,body,selector,0,{**proof,"unchanged_rho_exact":"0"}),
            )), f"record, node, rule, degree, body and proof bindings {name}")
        checks[f"{name}_masked_map_mock_routing_and_rejections"] = True

    # Signed series/kernel plumbing without any moment or same-Y cache job.
    context, nodes, rule, body = toy(cases[0][1], 1)
    n = len(context.target.index)
    selector = E.prepare_bulk_selector(context,nodes,rule,body,0,Q(1,10**11),Q(1,10**8))
    setup = A._setup(context,nodes,rule,body,selector,0)
    geometry = A.S.chart(setup[1]["corners"])
    keys = [(i,j) for i in range(n) for j in range(i,n)]
    force = [arb(i+1) for i in range(n)]
    hess = {ij:arb(100+k) for k,ij in enumerate(keys)}
    g = [Q(2) for _ in range(n)]
    h = dict.fromkeys(keys,Q(3))
    corrections_f = [-v for v in force]
    corrections_h = {ij:-v for ij,v in hess.items()}
    expected_mass = geometry.det_physical * sum(r for r,_ in geometry.corners)/6
    fun = Mock()
    fun.integral.return_value = context.target.normalization()*A.A(expected_mass)
    fun.weight_mass = sum(rule["volume_weights"][0],Q(0))
    fun.defect.cache_info.return_value.misses = 1
    cache_metadata = {"all64_convex_strain_lines_positive":True, "synthetic_signed_cache":True}
    with (patch.object(A.S,"Functional",return_value=fun) as f_factory,
          patch.object(A.S,"_eval_parts",return_value=(force,hess,g,h)) as eval_parts,
          patch.object(A.U,"uniform_cache",return_value=(corrections_f,corrections_h,cache_metadata)) as correction):
        result = A._scaled_component(context,nodes,rule,0,selector,setup,geometry)
        require(f_factory.call_args.args == (geometry,context.target.energy_scale,rule["volume_barycentric"],
                                             rule["volume_weights"][0]) and eval_parts.call_args.args[3:] == (fun,setup[2],selector["degree"]),
                "physical normalization, rule, degree and full gain cell forwarded")
        args = correction.call_args.args
        require(args[1] is setup[1] and args[2] == 0 and args[3] is nodes and args[4] is setup[2],
                "shared cache correction sees the same Y, element and gain cell")
        require(result["cache_consistency"] is cache_metadata and len(result["force_integral_minus_cache"]) == n
                and [(v["i"],v["j"]) for v in result["H_integral_minus_cache_upper"]] == keys,
                "all projected force and upper Hessian entries")
        for one in result["force_integral_minus_cache"] + [v["enclosure"] for v in result["H_integral_minus_cache_upper"]]:
            require(Q(one["lower_exact"]) <= 0 <= Q(one["upper_exact"]), "signed correction before remainder inflation")
        require(all(m.cache_clear.call_count == 1 for m in (fun.integral,fun.site_power,fun.defect)),
                "per-component caches cleared")
        correction.return_value = corrections_f, corrections_h, {"all64_convex_strain_lines_positive":False}
        require(rejects(lambda:A._scaled_component(context,nodes,rule,0,selector,setup,geometry)),
                "invalid cache strain path rejected")
    checks["scaled_packet_plumbing_signed_inflation_and_cache_cleanup_mocked"] = True
    # Finite nonzero-width coefficients touching/crossing zero reach the moment contraction.
    T = A.T
    c = arb(0, "1/8")
    poly = T.mul({(1,0):c},{(0,1):Q(2)})
    dual = A.JET.moment_dual(T.constant(c), {(1,1)}, lambda i,j: arb(4))
    value = A.JET.contract_factor(poly, dual)
    require((1,1) in poly and not value.is_zero() and 0 in value, "zero-crossing uncertainty preserved")
    checks["zero_crossing_coefficients_retained_to_signed_factor_contraction"] = True
    selector = {"max_DJ_norm_upper_exact":"3","max_D2J_norm_upper_exact":"5",
                "phi_prime_tail_upper_exact":"1/11","phi_second_tail_upper_exact":"1/7",
                "force_tail_budget_exact":"4","hessian_tail_budget_exact":"8"}
    f,H = A._tail_bounds(selector,Q(2),[Q(7),Q(2)],{(0,0):Q(8),(0,1):Q(4),(1,1):Q(1)})
    require(f == [Q(6,11),Q(4,11)] and H[(0,1)] == 2*(Q(6,7)+Q(4,11))
            and rejects(lambda:A._tail_bounds({**selector,"hessian_tail_budget_exact":"1"},Q(2),[Q(7)],
                                              {(0,0):Q(8)})), "mixed two-mass bound and overflow rejection")
    checks["mixed_H_tail_minimum_uses_two_valid_majorants_and_retains_budget"] = True
    return checks


def run():
    for p,h in PINS.items():
        bound(p,h)
    p = bound(rel(BATCH/"protocol.json"),PINS[rel(BATCH/"protocol.json")],True)
    pre = bound(rel(BATCH/"preflight.json"),PINS[rel(BATCH/"preflight.json")],True)
    require(p["source_sha256"] == pre["source_sha256"] == PINS[SOURCE]
            and pre["protocol_sha256"] == PINS[rel(BATCH/"protocol.json")]
            and p["actual_component_runs_started"] is pre["actual_component_runs_started"] is False
            and p["precision_bits"] == 640 and all_true(pre["checks"]) and pre["check_count"] == len(pre["checks"]) == 26,
            "sole prep protocol and preflight")
    bound(rel(BATCH/"exact_scaled_cohort_bulk_frozen.py"),PINS[SOURCE])
    for path,dig in {**p["input_hashes"], **p["failed_attempts"],
                     p["preflight_source_file"]:p["preflight_source_sha256"]}.items():
        bound(path,dig)
    absent_from_producer = TRANSITIVE not in p["input_hashes"]
    require(absent_from_producer and p["input_hashes"][SFILE] == PINS[SFILE], "transitive binding condition")
    sys.path.insert(0,str(ROOT/"scripts"))
    sys.path.insert(0,str((ROOT/TOYFILE).parent))
    A = importlib.import_module("exact_scaled_cohort_bulk")
    C = importlib.import_module("exact_cohort_context")
    toy = importlib.import_module("preflight_integrator")
    require(rel(toy.__file__) == TOYFILE and TOYFILE in p["input_hashes"], "bound synthetic fixture")
    structure(A)
    ctx.prec = 640
    A.core.A.cache_clear()
    forbidden = AssertionError("pilot data or actual moment/cache job during preparation")
    with (patch.object(A.S,"_certificate",side_effect=forbidden),
          patch.object(A.S.F,"_geometry",side_effect=forbidden),
          patch.object(A.S.F,"_asset",side_effect=forbidden),
          patch.object(A.S,"material_moment",side_effect=forbidden),
          patch.object(A.core,"local_terms",side_effect=forbidden)):
        report, checks = profiles(A,C,p,bound(MANIFEST,PINS[MANIFEST],True),
                                  bound(DOMAIN_BATCH_REVIEW,PINS[DOMAIN_BATCH_REVIEW],True))
        checks.update(toy_proofs(A,toy.toy))
    checks.update(dict.fromkeys(("frozen_input_sources_and_prior_development_history",
        "import_alias_identity_and_explicit_runtime_record","transitive_pilot_import_independently_bound",
        "memoized_moment_keys_and_physical_measure_algebra_reviewed",
        "actual_L_packet_absent_from_generic_evaluation_call_path","preserved_existing_axis_engines"), True))
    receipt = {
        "status":"PASS",
        "classification":"PASS_RECORD_NEUTRAL_SCALED_BULK_KERNEL_PREPARATION_ONLY",
        "source_sha256":PINS[SOURCE],
        "protocol_sha256":PINS[rel(BATCH/"protocol.json")],
        "preflight_sha256":PINS[rel(BATCH/"preflight.json")],
        "review_source_sha256":hash_bytes(Path(__file__).read_bytes()),
        "checks":checks, "chart_counts_per_record":report,
        "actual_record_domains_reviewed_here":False,
        "actual_moments_cache_sums_force_or_H_computed":False,
        "authorizes_actual_jobs_alone":False,
        "scope":"Module arithmetic and explicit dataflow for each record's original domain. The record's terminal geometry/domain proof, reviewed selector and exact derived route proof still gate every actual component.",
        "supplemental_runtime_source_hashes":{TRANSITIVE:PINS[TRANSITIVE]},
        "supplemental_source_gate":"The production launcher must verify the supplemental transitive source together with frozen input hashes and runtime import paths before evaluating any component. The preparer's protocol omitted this direct hash.",
    }
    data = (json.dumps(receipt,indent=2,allow_nan=False)+"\n").encode()
    with (BATCH/"independent_review.json").open("xb") as f:
        f.write(data)
    print(json.dumps({"status":receipt["classification"],"sha256":hash_bytes(data),"checks":len(checks),"routes":sum(sum(r.values()) for r in report.values())}))


if __name__ == "__main__":
    run()

