"""Frozen direct reducer prelaunch readback; preparation and terminal metadata only."""
import ast
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
from random import Random
import sys
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
ROOT = OUT.parents[2]
SOURCE = "scripts/exact_cohort_record_reduction.py"
PINS = {
    "source_sha256": "8e3177434299fbffeee2fa5161dd4025fc0e9cd0b7f46be85fceaeccae204ed2",
    "protocol_sha256": "772f47398fd7a376fab274bb1d34c9cf7e3fe0bdd661fdd0b1f8b790eccdcc28",
    "preflight_sha256": "a598a545502f09287f909cf8d2753757da1c02863fe2619d06841b8ce1e79496",
}
CD = "runs/exact_integration_cohort_v1/remaining15_components_observed/CD_terminal_observation.json"
CD_SHA = "95697af788d6f746f236f634c06a103c0c5ecc0e00428b52e80817e6e804e94e"


def insist(condition, detail):
    if not condition:
        raise ValueError(detail)


def sha(raw):
    return sha256(raw).hexdigest()


def read(path, expected=None, parse=True):
    safe = (ROOT/path).resolve()
    insist(safe.is_relative_to(ROOT) and str(path) == safe.relative_to(ROOT).as_posix(),
           f"foreign file {path}")
    insist("/components/" not in str(path), "actual component read during prelaunch")
    raw = safe.read_bytes()
    insist(expected is None or expected == sha(raw), f"changed input {path}")
    return json.loads(raw) if parse else raw


def no(fun):
    try:
        fun()
    except ValueError:
        return True
    return False


def yes(v):
    return isinstance(v, dict) and bool(v) and all(x is True for x in v.values())


def pack(a,b):
    return {"lower_exact":str(a),"upper_exact":str(b)}


def arithmetic(D):
    rng = Random(611)
    n = 9
    f, H = [D.ZERO]*n, {}
    expected_f = [[Q(0),Q(0)] for _ in range(n)]
    full = [[[Q(0),Q(0)] for _ in range(n)] for _ in range(n)]
    seen = set()
    for _ in range(24):
        ids = rng.sample(range(n),rng.randint(2,6))
        g, h = [], []
        for i,a in enumerate(ids):
            l,u = sorted((Q(rng.randint(-81,81),2**rng.randint(0,11)) for _ in range(2)))
            g.append(pack(l,u))
            expected_f[a][0] += l
            expected_f[a][1] += u
            for j,b in enumerate(ids[i:],i):
                ml,mu = sorted((Q(rng.randint(-81,81),2**rng.randint(0,11)) for _ in range(2)))
                h.append({"i":i,"j":j,"enclosure":pack(ml,mu)})
                for x,y in {(a,b),(b,a)}:
                    full[x][y][0] += ml
                    full[x][y][1] += mu
                seen.add(tuple(sorted((a,b))))
        D._scatter(f,H,ids,g,h)
    exp_f = list(map(tuple,expected_f))
    insist(f == exp_f and H == {(i,j):tuple(full[i][j]) for i,j in seen},
           "signed force and reordered upper-H scatter")
    h_abs = [sum(max(abs(a),abs(b)) for a,b in row) for row in full]
    insist(D._h_rows(n,H) == h_abs,"mirror each off-diagonal exactly once")
    R = [[float(Q(rng.randint(-128,128),2**rng.randint(0,15))) for _ in range(n)] for _ in range(n)]
    expect_rf, expect_b, norms = [], [], []
    for row in R:
        exact = list(map(Q.from_float,row))
        endpoint_products = [sorted((a*l,a*u)) for a,(l,u) in zip(exact,expected_f)]
        expect_rf.append(tuple(sum(x) for x in zip(*endpoint_products)))
        expect_b.append(sum(abs(a)*y for a,y in zip(exact,h_abs)))
        norms.append(sum(map(abs,exact)))
    insist(D._dyadic_reduce(R,f,h_abs) == (expect_rf, expect_b, max(norms)), "dyadic vs independent Fraction endpoints")
    insist(D._strict_test(Q(1,8),Q(1,4),Q(1,8),Q(1,8),Q(1,2))[-2:] == (True,Q(2,5))
           and D._strict_test(Q(1,8),Q(1,4),Q(1,8),Q(1,4),Q(1,2))[-2:] == (False,None),
           "strict root self-map including equality")
    return {"synthetic_signed_scatter_and_matrix_endpoints":True,"strict_gate_equality_unresolved":True}


def source_checks(A,D):
    ast_tree = ast.parse(read(SOURCE,PINS["source_sha256"],False))
    src = {n.name:ast.unparse(n) for n in ast_tree.body if isinstance(n,ast.FunctionDef)}
    insist(src["admit"].index("_terminal_part(") < src["admit"].index("D._write(")
           and "/components/" not in src["admit"] and "_reduce(" not in src["admit"]
           and src["reduce"].index("_terminal_part(") < src["reduce"].index("_reduce(")
           and 'prefix.mkdir(exist_ok=False)' in src["reduce"],
           "separate immutable pre-payload admission before any scatter")
    terminal = src["_terminal_part"]
    insist(".os_terminal.json" in terminal and "D._terminal(" in terminal
           and "os_exit['pid'] == terminal['pid'] == started['pid']" in terminal
           and "type(os_exit['exit_code']) is int" in terminal and ".returncode" not in terminal
           and "results.json" in terminal and "_manifest(" in terminal, "external per-part OS, PID and manifest assumptions")
    src_reduce = src["_reduce"]
    insist(src_reduce.index("D._scatter(") < src_reduce.index("D._h_rows(") < src_reduce.index("D._dyadic_reduce(")
           < src_reduce.index("D._strict_test(") < src_reduce.index("D._minimum(")
           and "if passes else None" in src_reduce and "'root_assertion_independently_reviewed': False" in src_reduce,
           "signed arithmetic and conditional SPD only")
    return {"admission_and_actual_read_separated":True,"exact_preconditioned_endpoint_algebra_and_conditional_spd":True,
            "worker_manifest_and_wait_receipt_enforced_without_whole_supervisor":True}


def spd_binding(D,ctx,z,r,provider,summaries,reviews,bindings):
    rec = ctx.record_id
    summary, = [a for a in summaries["records"] if a["record_id"] == rec]
    wrapper, = [a for a in reviews["record_audits"] if a["record_id"] == rec]
    c = D._bound(summary["file"], summary["sha256"], bindings)
    a = D._bound(wrapper["review_file"],wrapper["review_sha256"],bindings)
    b = c["root_binding"]
    row = ctx.row
    insist(a == wrapper["audit"] and a["producer_result_sha256"] == summary["sha256"]
           and a["status"] == "PASS" and yes(a["checks"])
           and c["status"] == "CERTIFIED_UNIFORM_STRICT_LOCAL_ENERGY_MINIMUM"
           and c["whole_tube_Hessian_positive"] is True
           and sha(ctx.R.tobytes(order="C")) == b["preconditioner_binary_sha256"]
           and Q(a["normal_root_radius_exact"]) == r
           and Q(a["independent_fulltube_Z_upper_exact"]) == z
           and b["independent_Z_source"]["file"] == provider
           and all(b[k] == row[k] for k in ("record_id","root_certificate_sha256","parameter_cell_sha256",
               "model_observation_contract_sha256","branch_domain_id","point_cache_file","point_cache_sha256","frozen_binary_cache_sha256"))
           and Q(a["uniform_Hessian_lower_bound_exact"]) ==
               Q(c["uniform_free_coordinate_Hessian_lower_bound_exact"]) > 0,
           "same-record SPD hypotheses are available conditionally")
    return {"cached_spd_audit_file":wrapper["review_file"],"cached_spd_audit_sha256":wrapper["review_sha256"],
            "exact_radius":str(r),"cached_Z_upper_exact":str(z),"cached_provider":provider,
            "center":row["proof_center_state_file"]}


def format_tests(A,D,ctx,p,domains,body,routes):
    # These payloads are invented in memory using the frozen case maps.
    identity = {"record_id":ctx.record_id,"dose_exact":str(ctx.dose),"precision_bits":640,
        **{k:ctx.row[k] for k in ("parameter_cell_sha256","point_cache_sha256","model_observation_contract_sha256","predictor_sha256")}}
    for i in (0,6,72):
        case, proof = domains["bulk_selectors"][i], routes["bulk"][i]
        ids = case["global_indices"]
        n = len(ids)
        packet = {**identity, **{k:case[k] for k in ("element","chart","free_local","active","global_indices",
                      "nodal_hulls_sha256","auxiliary_element_rule_sha256","rho_exact","degree")},
                  "status":"UNIFORM_COHORT_BULK_COMPONENT_UNREVIEWED",
                  "integration_chart":proof["integration_chart"],"record_component_proof_sha256":D._digest(proof),
                  "cache_consistency":{"all64_convex_strain_lines_positive":True,
                      "minimum_cache_convex_line_J_lower":pack(Q(1,4),Q(1,2)),
                      "max_coefficientwise_delta_q_upper":pack(0,0)}}
        fmax = proof.get("scaled_force_tail_max_exact","0")
        hmax = proof.get("scaled_H_entry_tail_max_exact","0")
        def contain(x):
            x = Q(x)
            return pack(0,0) if x == 0 else pack(0,1)
        packet["force_two_mass_tail_upper_exact"] = [fmax]*n
        packet["force_two_mass_tail_upper"] = [contain(fmax)]*n
        packet["H_integral_minus_cache_upper"] = [{"i":j,"j":k,"enclosure":pack(0,0),
            "two_mass_tail_upper_exact":hmax,"two_mass_tail_upper":contain(hmax)}
            for j in range(n) for k in range(j,n)]
        if i == 6:
            validator = lambda pkt:A._scaled_packet(D,ctx,i,case,body["triangles"][i],proof,pkt,p)
        else:
            validator = lambda pkt:D._packet_meta(ctx,"bulk",i,case,pkt,body,p)
        insist(validator(packet) == ids and all(no(lambda:validator({**packet,key:value}))
                    for key,value in (("record_id","foreign"),("nodal_hulls_sha256","foreign"),("degree",case["degree"]+1))),
               f"synthetic original/chart and malformed bulk packet {i}")
    case = domains["surface_cases"][0]
    ids = case["global_indices"]
    n = len(ids)
    packet = {**identity,"status":"UNIFORM_COHORT_SURFACE_COMPONENT_UNREVIEWED","global_indices":ids,
        "face":0,"edge":0,"order":case["order"],"rho_exact":case["domain"]["rho_exact"],
        "q0_exact":case["domain"]["q0_exact"],"material_mass":pack(1,1),"auxiliary_mass":pack(1,1),
        "force_two_mass_tail_upper":[pack(0,0)]*n,
        "H_integral_minus_cache_upper":[{"two_mass_tail_upper":pack(0,0)}]*(n*(n+1)//2),
        "actual_cache_sites":[{"quadrature_index":i,"weight_applied_exactly_once":s["auxiliary_normalized_weight_exact"],
            "unweighted_F_max":pack(0,0),"unweighted_H_max":pack(0,0)}
            for i,s in enumerate(case["cache_site_domains"])]}
    check = lambda pkt:D._packet_meta(ctx,"surface",0,case,pkt,body,p)
    insist(check(packet) == ids and no(lambda:check({**packet,"actual_cache_sites":packet["actual_cache_sites"][:-1]}))
           and no(lambda:check({**packet,"q0_exact":"0"})), "synthetic surface sites, trigger domain and order")
    return {"synthetic_packet_metadata_all_four_routes":True}


def run():
    folder = OUT.relative_to(ROOT).as_posix()
    read(SOURCE,PINS["source_sha256"],False)
    read(f"{folder}/exact_cohort_record_reduction_frozen.py",PINS["source_sha256"],False)
    proto = read(f"{folder}/protocol.json",PINS["protocol_sha256"])
    pre = read(f"{folder}/preflight.json",PINS["preflight_sha256"])
    insist(proto["status"] == "FROZEN_REMAINING15_DIRECT_REDUCTION_PREPARATION_REVIEW_PENDING"
           and pre["status"] == "PASS_GENERIC_DIRECT_REDUCTION_PREPARATION_ONLY" and yes(pre["checks"])
           and proto["source_sha256"] == pre["source_sha256"] == PINS["source_sha256"]
           and pre["protocol_sha256"] == PINS["protocol_sha256"]
           and pre["records"] == len(proto["records"]) == 15
           and all(d["actual_component_files_opened"] is d["actual_global_sums_run"] is False for d in (proto,pre)),
           "frozen preparation state")
    for key, val in {**proto["input_hashes"],**proto["failed_attempts"]}.items():
        read(key,val,False)
    sys.path.insert(0,str(ROOT/"scripts"))
    import exact_cohort_record_reduction as A
    rt,plan,jobs,expect,hashes = A._boot()
    D=rt.D
    insist(A.RECORDS == proto["records"] == ["optoGEF2_small_col_"+c for c in "CDEFGHIJKMNOPQR"],
           "all remaining IDs")
    checks = {**source_checks(A,D),**arithmetic(D)}
    summary = D._bound(D.ENERGY_RESULT,D.PREPARATION_PINS[D.ENERGY_RESULT],hashes)
    review = D._bound(D.ENERGY_REVIEW,D.PREPARATION_PINS[D.ENERGY_REVIEW],hashes)
    insist(review["status"] == "PASS" and yes(review["checks"]), "reviewed cached-energy cohort")
    premises, observed = {}, {}
    obs = read(CD,CD_SHA)
    insist(obs["status"] == "TWO_COMPLETE_RECORD_PARTS_OBSERVED_NOT_ROOTS"
           and [v["record_id"] for v in obs["records"]] == A.RECORDS[:2], "bound root's explicit terminal observation")
    blocked = AssertionError("prelaunch must not open any local packet or compute a force/cache")
    checked = D._bound
    def check_source(path,*args,**kwargs):
        if "/components/" in path:
            raise blocked
        return checked(path,*args,**kwargs)
    with (patch.object(D,"_bound",side_effect=check_source),
          patch.object(A,"_reduce",side_effect=blocked),
          patch.object(D,"_minimum",side_effect=blocked),
          patch.object(rt.E.core,"local_terms",side_effect=blocked),
          patch.object(rt.I,"bulk_component",side_effect=blocked),
          patch.object(rt.E,"surface_component",side_effect=blocked)):
        for record in A.RECORDS:
            ctx,p,domains,body,routes,links,stage = A._load_record(rt,plan,jobs,expect,record,hashes)
            eta,z,r,provider = D._old_premises(ctx,p,hashes)
            insist(eta+z*r < r and r == ctx.old_radius, "same selected cached premise")
            premises[record] = spd_binding(D,ctx,z,r,provider,summary,review,hashes)
            if record in A.RECORDS[:2]:
                w, = [v for v in obs["records"] if v["record_id"] == record]
                parts = {}
                for part,total in (("bulk",144),("surface",24)):
                    pins = w["parts"][part]
                    declared = A._terminal_part(D,stage,links,record,part,pins,hashes)
                    insist(len(declared) == pins["components"] == total, "external component metadata staged")
                    parts[part] = total
                    insist(no(lambda:A._terminal_part(D,stage,links,record,part,
                                     {**pins,"os_terminal_sha256":"0"*64},{})), "foreign OS hash rejected")
                observed[record] = parts
            if record == A.RECORDS[0]:
                checks.update(format_tests(A,D,ctx,p,domains,body,routes))
    insist(all(premises["optoGEF2_small_col_"+c]["center"].endswith("/proof_center_state.json") for c in "FK")
           and premises["optoGEF2_small_col_M"]["cached_provider"] == "evidence/parametric_actual_M_independent_review.json",
           "exact F/K recentered and M nested cached branches")
    checks.update(dict.fromkeys(("frozen_inputs","all_15_exact_same_record_C_S_R_cache_and_original_radius",
         "all_15_conditional_same_record_spd_bindings","all_168_canonical_metadata_each_C_D_from_external_terminals",
         "failed_wrong_OS_digest_rejected","no_component_payload_opened","F_K_recentered_and_M_provider"),True))
    result = {"status":"PASS",**PINS,"checks":checks,
              "scope":"Single-record direct reduction only after both explicit successful part terminals. All computed results require independent actual review before any root or minimum claim.",
              "component_payloads_opened":False,"global_record_reductions_run":False,
              "actual_terminal_metadata_observation":{"file":CD,"sha256":CD_SHA,"packet_count_by_record":observed},
              "old_proof_and_spd_bindings":premises,
              "review_source_sha256":sha(Path(__file__).read_bytes())}
    raw = (json.dumps(result,indent=2,allow_nan=False)+"\n").encode()
    gate_path=OUT/"prelaunch_independent_review.json"
    with gate_path.open("xb") as writer:
        writer.write(raw)
    gate_sha = sha(raw)
    A._prelaunch(D,gate_sha,{})
    print(json.dumps({"status":"PASS","gate_sha256":gate_sha,"checks":len(checks),"records":len(premises)}))


if __name__ == "__main__":
    run()

