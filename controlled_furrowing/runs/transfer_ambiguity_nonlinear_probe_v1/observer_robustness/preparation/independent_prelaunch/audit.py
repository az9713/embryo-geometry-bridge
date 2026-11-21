"""Independent no-field preparation audit of four fixed stored-field observers."""
import ast
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
OUT=HERE.parent
ROOT=OUT.parents[3]
SOURCE="scripts/transfer_ambiguity_observer_sweep.py"
ID={"source_sha256":"9eacb1de680cefb84d90d603cc1aca4dd6ef4892372a91cffa3480651c5874a8",
    "protocol_sha256":"dc48cee2cdf5f0da3b49f6972dc37191cf482a238f37b3185f8ceaf1f93bbeea",
    "preflight_sha256":"1e9e318c7f6ae06d99ddb1a12eb64c78f101812dc16631b56cc582e474f5247b",
    "primary_protocol_sha256":"9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94"}
NATIVE="e1ff63fddd6981f9ef8f05183c884537546cb66d59b9ec063911b22a1c50d27d"
V=("current_disk__spatial_max","current_disk__regional_mean",
   "current_radial_line__spatial_max","current_radial_line__regional_mean")


def need(v, why):
    if not v: raise ValueError(why)


def sha(b):
    return sha256(b).hexdigest()


def bound(path,h,parse=True):
    p=(ROOT/path).resolve()
    need(p.is_relative_to(ROOT),"local artifact")
    raw=p.read_bytes()
    need(sha(raw)==h,f"changed {path}")
    return json.loads(raw) if parse else raw


def rel(p):
    return p.resolve().relative_to(ROOT).as_posix()


def rejects(f):
    try:
        f()
    except ValueError:
        return True
    return False


def yes(d):
    return isinstance(d,dict) and bool(d) and all(v is True for v in d.values())


def exact_graph(A):
    # Total shortening 1+6r-4r² peaks at 3/4; component optima separated.
    top=[[(0,9),(Q(1,2),Q(37,4)),(1,10)]]
    floor=[[(0,0),(Q(1,2),Q(9,4)),(1,3)]]
    vals=A.extract_variants(A.observer.observe_graphs(top,floor,roi_um=(0,1),reference_apical_um=10))
    keys=A.contract.OBSERVABLES
    expected=((Q(13,4),Q(3),Q(1),Q(3)),(Q(13,4),Q(3),Q(1,2),Q(5,2)),
              (Q(13,4),Q(8,3),Q(1),Q(3)),(Q(13,4),Q(8,3),Q(2,3),Q(2)))
    need(tuple(vals)==V and all(list(map(Q,vals[v][k]))==[x,x] for v,row in zip(V,expected) for k,x in zip(keys,row))
         and Q(13,4)<1+3,"independent line/disk integrals and true combined height maximum")
    negative=A.extract_variants(A.observer.observe_graphs(
        [[(0,11),(Q(1,2),11),(1,11)]],[[(0,-2),(Q(1,2),-2),(1,-2)]],
        roi_um=(0,1),reference_apical_um=10))
    need(all([list(map(Q,row[k])) for k in keys]==[[Q(-3)]*2,[Q(-3)]*2,[Q(-1)]*2,[Q(-2)]*2]
             for row in negative.values()),"signed inward/outward graph components retained")


def native_fixture(A, protocol):
    primary=bound(A.PRIMARY_PROT,ID["primary_protocol_sha256"])
    native=bound(A.REPLAY+"/replay_frozen.py",NATIVE,parse=False)
    syntax={f.name:ast.unparse(f) for f in ast.parse(native).body if isinstance(f,ast.FunctionDef)}
    rp=syntax["replay"]
    need(all(repr(k) in rp for k in ("actual_results_sha256","checkpoints","independently_recomputed_small_residuals",
                                     "all_hashes","endpoints","new_equilibrium_solves","heldout_response_tables_read"))
         and "'PASS_INDEPENDENT_ACTUAL_NONLINEAR_STORED_FIELD_REPLAY'" in rp
         and "'field_manifest_sha256'" not in rp
         and "'receipt_sha256': sha(output / 'receipt.json')" in rp,
         "observer consumes exact native independent-review and terminal schema")
    native_p=bound(A.REPLAY+"/protocol.json",A.PINS[A.REPLAY+"/protocol.json"])
    archive=dict(native_p["prerequisite_sha256"])
    prepmap={A.REPLAY+"/protocol.json":A.PINS[A.REPLAY+"/protocol.json"],
             A.REPLAY+"/preflight.json":A.PINS[A.REPLAY+"/preflight.json"],A.REPLAY+"/replay_frozen.py":NATIVE}
    need(all(({**archive,**prepmap}).get(k)==v for k,v in A.PINS.items()),"every required primary input retained by real replay")
    flat=[]
    for e in primary["endpoints"]:
        eid=e["endpoint_id"]
        groups=[("small_zero","small_zero",e["small_zero"])]
        groups.extend((f"small_{i+1:02d}","small_record",r) for i,r in enumerate(e["small_records"]))
        groups.extend((f"large_{i:02d}",k,x) for i,(k,x) in enumerate(zip(
            ("large_zero","large_checkpoint","large_checkpoint","large_final"),e["large_checkpoints"],strict=True)))
        for label,role,step in groups:
            folder=f"{A.PRIMARY_ACTUAL}/{eid}/{label}"
            flat.append({"endpoint_id":eid,"checkpoint":label,"role":role,"step":step,
                         "field":folder+"/field.json","verification":folder+"/verification.json"})
    need(len(flat)==44==len({p["field"] for p in flat}) and sum(x["role"]=="small_record" for x in flat)==34
         and flat==protocol["frozen_checkpoint_plan"]==A.plan_checkpoints(primary),"independent exact frozen primary order and plan")
    rows=[{**x,"field_sha256":sha(x["field"].encode()),"verification_sha256":sha(x["verification"].encode())} for x in flat]
    rs,ts="a"*64,"b"*64
    result={"status":"COMPLETE_CONDITIONAL_STORED_FIELD_COMPARISON","verified_outcome_fields":36,"small_observer_intervals":136,
            "endpoints":[{"endpoint_id":eid,"status":"COMPLETE_CONDITIONAL_STORED_FIELDS","verified_final_outcome_fields":18,
              "checkpoints":[{k:v for k,v in x.items() if k!="endpoint_id"} for x in rows if x["endpoint_id"]==eid],
              "small_residual":{"all68_residuals":{"constraints":68}}} for eid in A.ENDPOINTS]}
    manifest={"status":"COMPLETE_44_PRIMARY_CHECKPOINTS","primary_results_sha256":rs,"primary_terminal_sha256":ts,"fields":rows}
    hashes={**A.PINS,A.PRIMARY_RESULT:rs,A.PRIMARY_TERMINAL:ts,
            **{x[k]:x[k+"_sha256"] for x in rows for k in ("field","verification")}}
    review={"status":"PASS_INDEPENDENT_ACTUAL_NONLINEAR_STORED_FIELD_REPLAY","actual_results_sha256":rs,"checkpoints":44,
            "independently_recomputed_small_residuals":136,"new_equilibrium_solves":False,"heldout_response_tables_read":False,
            "all_hashes":hashes,"endpoints":[{"endpoint_id":eid,"checkpoints":[{"checkpoint":x["checkpoint"],
                 "field_sha256":x["field_sha256"]} for x in rows if x["endpoint_id"]==eid]} for eid in A.ENDPOINTS]}
    need(A.planned_manifest(result,flat)==rows,"frozen primary result parsing")
    A.require_primary_review(review,manifest,rs,ts,rows)
    actual_step=deepcopy(result)
    actual_step["endpoints"][1]["checkpoints"][5]["step"]={"record_id":"foreign"}
    need(rejects(lambda:A.planned_manifest(actual_step,flat)),"mixed record/step in primary manifest rejected")
    for k in ("field","verification"):
        bad={**review,"all_hashes":{**hashes,rows[5][k]:"5"*64}}
        need(rejects(lambda:A.require_primary_review(bad,manifest,rs,ts,rows)),"stale field and verification links both rejected")
    for bad in ({**manifest,"fields":rows[:-1]},{**manifest,"fields":[rows[1],rows[0],*rows[2:]]},
                {**manifest,"primary_terminal_sha256":"4"*64}):
        need(rejects(lambda:A.require_primary_review(review,bad,rs,ts,rows)),"partial, reordered and stale terminal manifest rejected")
    for changes in ({"status":"FAIL"},{"checkpoints":43},{"new_equilibrium_solves":"false"},
                    {"actual_results_sha256":"3"*64},{"independently_recomputed_small_residuals":135}):
        need(rejects(lambda:A.require_primary_review({**review,**changes},manifest,rs,ts,rows)),"non-native/failed actual receipt rejected")
    return primary,rows,result


def synthetic_run(A,protocol,primary,rows,result):
    # One synthetic interval shift per observer, no input file parses or writes.
    delta=Q(1)
    shifts=((Q(-1),Q(-1)),(Q(2),Q(5,2)),(Q(-3,2),Q(0)),(Q(1,10),Q(1,2)))
    obs={r["record_id"]:r["observed_exact"] for r in primary["endpoints"][0]["small_records"]}
    def measurements(state,params,step):
        rid=step["record_id"]
        return {"variant_observables_um_exact":{
            mode:{k:[str(Q(v)+a),str(Q(v)+b)] for k,v in obs[rid].items()}
            for mode,(a,b) in zip(V,shifts)}}
    data={r:measurements({}, {}, {"record_id":r})["variant_observables_um_exact"][V[0]] for r in obs}
    for stage,case in zip(result["endpoints"],primary["endpoints"]):
        stage["small_residual"]=A.threshold(case["small_records"],data,delta)
    def fake_read(path,expected):
        return primary if rel(path)==A.PRIMARY_PROT else {}
    with (patch.object(A,"pinned",side_effect=fake_read) as opened,
          patch.object(A,"measure_field",side_effect=measurements) as measured,
          patch.object(A,"save",return_value="4"*64)):
        p={**protocol,"conditional_source_allowance_um_exact":"1"}
        out=A.run(p,result,rows,{})
        need(opened.call_count==35 and measured.call_count==34,"exactly one observer parse of each of 34 small and none of large/zero cases")
        need(out["primary_replay_exactly_matches"] is True and out["measured_small_fields"]==34
             and out["small_constraint_evaluations"]==544 and out["additional_equilibria"]==out["additional_force_evaluations"]==0
             and out["retuning"] is False,"same fields over four independently predeclared finite conventions")
        expected=[("ALLOWANCE_SATISFIED",[Q(1),Q(1)]),("ALLOWANCE_EXCEEDED",[Q(2),Q(5,2)]),
                  ("UNRESOLVED",[Q(0),Q(3,2)]),("ALLOWANCE_SATISFIED",[Q(1,10),Q(1,2)])]
        for v,(status,rng) in zip(V,expected):
            q=out["common_pair_requirement_by_variant"][v]
            need(q["status_at_frozen_allowance"]==status and list(map(Q,q["required_allowance_um_exact"]))==rng,
                 "correct max over both endpoints and 68 signed residuals")
        total=out["common_pair_requirement_across_all_four_conventions"]
        need(total["status_at_frozen_allowance"]=="ALLOWANCE_EXCEEDED" and total["required_allowance_um_exact"]==["2","5/2"]
             and set(out["primary_threshold_outcome_by_endpoint"].values())=={"ALLOWANCE_SATISFIED"},
             "overall failed alternative does not change primary decision")
        stale=deepcopy(result)
        stale["endpoints"][0]["small_residual"]["declared_allowance_um_exact"]="2"
        need(rejects(lambda:A.run(p,stale,rows,{})),"primary equality is exact dictionary replay, no numerical or policy relaxation")


def verify_source_order(src):
    f={node.name:ast.unparse(node) for node in ast.parse(src).body if isinstance(node,ast.FunctionDef)}
    s=f["admit_actual"]
    need(s.index("review = pinned(")<s.index("manifest = pinned(")<s.index("require_primary_review(")
         <s.index("verify_primary_supervisor(")<s.index("sha(ROOT / row[kind])"),
         "all 44 logical links and both process exits must admit before any field or verification byte hashed")
    need(f["launch"].index("prelaunch_bindings(")<f["launch"].index("admit_actual(")<f["launch"].index("run(protocol"),
         "independent gate and admission before any observer loop")
    k=f["verify_primary_supervisor"]
    need("type(worker['observed_exit_code']) is int" in k and "worker['observed_exit_code'] == 0" in k
         and "child['pid'] == worker['pid']" in k and "terminal_sha256" in k
         and "worker and supervisor results disagree" in k and "receipt_sha256" in s,
         "native worker/supervisor terminal success and running process boundaries")
    loop=f["run"]
    need("primary_receipts[eid]" in loop and "2 * len(VARIANTS) * 68" in loop
         and "if row['endpoint_id'] != eid or row['role'] != 'small_record'" in loop,
         "fixed outcomes no post-selection, all comparisons from fixed existing fields")


def run():
    folder=rel(OUT)
    raw=bound(SOURCE,ID["source_sha256"],False)
    bound(folder+"/transfer_ambiguity_observer_sweep_frozen.py",ID["source_sha256"],False)
    p=bound(folder+"/protocol.json",ID["protocol_sha256"])
    pre=bound(folder+"/preflight.json",ID["preflight_sha256"])
    need(pre["status"]=="PASS" and yes(pre["checks"]) and all(pre[k]==ID[k] for k in ("source_sha256","protocol_sha256","primary_protocol_sha256"))
         and p["source_sha256"]==ID["source_sha256"] and tuple(p["variant_order"])==V
         and p["primary_kept"]==V[0] and p["status"]=="PREDECLARED_SECONDARY_OBSERVER_SWEEP"
         and all(pre[k]==0 for k in ("primary_actual_fields_opened","actual_observer_loops_run","additional_force_evaluations","additional_equilibria"))
         and all(p[k] is False for k in ("additional_force_evaluations","additional_equilibrium_solves","retuning","heldout_large_response_read")),
         "frozen no-primary-field prep and universal static scope")
    for path,h in {**p["source_hashes"],**p["input_hashes"]}.items():
        need("/actual/" not in path,"no actual primary field in preparation hash set")
        bound(path,h,False)
    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_observer_sweep as A
    original_open=Path.open
    def forbid_primary_actual(path,*a,**kw):
        need(not path.resolve().is_relative_to(ROOT/A.PRIMARY_ACTUAL),"no primary actual bytes accessed in independent preparation")
        return original_open(path,*a,**kw)
    with patch.object(Path,"open",forbid_primary_actual):
        need(p==A.prepare_protocol(),"independent source and selection binding readback")
        verify_source_order(raw)
        exact_graph(A)
        primary,rows,result=native_fixture(A,p)
        synthetic_run(A,p,primary,rows,result)
    variants=[{"status":"PASS","checks":{"gate":True}},{"status":"PASS","checks":{}},
              {"status":"PASS","checks":{"gate":1}},{"status":"FAIL","checks":{"gate":True}},
              {"status":"PASS","checks":{"gate":"true"}}]
    need(list(map(A.good_review,variants))==[True,False,False,False,False],"literal review check flags required")
    checks=dict.fromkeys(("all_current_source_and_prep_bindings","exact_independent_four_graph_conventions_and_noncoincident_maxima",
        "signed_negative_components","native_primary_actual_receipt_and_true_OS_input_schema",
        "all_44_plan_rows_and_88_field_verification_hash_links","separate_whole_admission_before_field_bytes",
        "mixed_stale_partial_failed_and_nonliteral_reviews_rejected","all_544_synthetic_residuals_and_all_four_quantifiers",
        "exact_primary_replay_and_separate_secondary_decisions","no_actual_field_open_force_or_solves"),True)
    doc={"status":"PASS",**ID,"checks":checks,"review_source_sha256":sha(Path(__file__).read_bytes()),
         "scope":"Preparation-only four fixed static observer variants. No primary appearance, fit, source claim, response retuning, or equilibrium solve. A future sweep requires complete independently reviewed externally pinned primary output and all 44 ordered field/verification pairs before any field read.",
         "actual_field_payloads_opened":0,"actual_observer_sweeps_run":0,
         "expected_future_small_observer_fields":34,"expected_future_constraint_comparisons":544}
    b=(json.dumps(doc,indent=2,allow_nan=False)+"\n").encode()
    receipt=OUT/"prelaunch_independent_review.json"
    with receipt.open("xb") as f:f.write(b)
    digest=sha(b)
    need(all(A.prelaunch_bindings(receipt,digest)[k]==v for k,v in ID.items()),"real gate readback")
    print(json.dumps({"status":"PASS","review_sha256":digest,"checks":len(checks),"actual_fields_opened":0}))


if __name__=="__main__":
    run()

