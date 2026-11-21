"""Independent preparation-only centered reduction audit; no packet or cache evaluations."""
import ast
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from pathlib import Path
from random import Random
import sys
from types import ModuleType
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
OUT=HERE.parent
ROOT=OUT.parents[2]
SOURCE="scripts/exact_centered_cohort_reduction.py"
ARCHIVE=OUT/"development_revision_0"
PINS={"source_sha256":"e5e11f7a048af122e9fd0d41629fcd345bf47f0f721991b464f897d2f119170a",
      "protocol_sha256":"e708a46f4c85afc3bf4d3229e57ac376c74777de223155ec7a104e6af67da66e",
      "preflight_sha256":"7be7c4fb84b2a23b61a594404a16c205be49d93bb16279f4a4941564746ee6f5"}
CALLER_REVIEW="31ed5e5a15e73242414a44ed3b9412679d01d0e590eedbf4d94973dcbb556ca2"
RECORD="optoGEF2_small_col_C"


def need(v,msg):
    if not v:
        raise ValueError(msg)


def sha(a):
    return sha256(a).hexdigest()


def read(p,h,parse=True):
    file=(ROOT/p).resolve()
    need(file.is_relative_to(ROOT) and file.relative_to(ROOT).as_posix()==p,"canonical workspace path")
    need("/components/" not in p,"no component opens")
    raw=file.read_bytes()
    need(sha(raw)==h,f"changed {p}")
    return json.loads(raw) if parse else raw


def rel(p):
    return Path(p).resolve().relative_to(ROOT).as_posix()


def yes(d):
    return isinstance(d,dict) and bool(d) and all(x is True for x in d.values())


def rejects(f):
    try:
        f()
    except ValueError:
        return True
    return False


def arith(D,A):
    rng=Random(197)
    q=lambda a,b:Q(a,2**b)
    rnd=lambda:q(rng.randrange(-17,18),rng.randrange(0,6))
    R=[[float(rnd()) for _ in range(4)] for _ in range(4)]
    S=[[float(rnd()) for _ in range(5)] for _ in range(4)]
    H={k:tuple(sorted((rnd(),rnd()))) for k in combinations_with_replacement(range(4),2)}
    H[(0,3)]=(Q(0),Q(0))
    widths=[Q(17,113),Q(2,15),Q(3,4),Q(7,31)]

    def add(seq):
        seq=list(seq)
        return tuple(sum(a) for a in zip(*seq)) if seq else (Q(0),Q(0))

    def mul(x,y):
        vals=[a*b for a in x for b in y]
        return min(vals),max(vals)
    Rh=[[add(mul((Q.from_float(R[i][k]),)*2,H[tuple(sorted((k,j)))]) for k in range(4))
         for j in range(4)] for i in range(4)]
    rhs=[[add(mul(row[k],(Q.from_float(S[k][j]),)*2) for k in range(4)) for j in range(5)] for row in Rh]
    bound=[sum(max(map(abs,x))*r for x,r in zip(row,widths)) for row in Rh]
    require=A._bulk_RH_effects(R,H,widths,S)
    need(require==(bound,rhs),"nontrivial interval RH before abs and signed RH*S on five columns")
    p=[tuple(sorted((rnd(),rnd()))) for _ in range(4)]
    g=[tuple(sorted((rnd(),rnd()))) for _ in range(4)]
    s=[tuple(sorted((rnd(),rnd()))) for _ in range(4)]
    d=[(Q(-1,19),Q(1,20)),(Q(-3,23),Q(1,20)),(Q(0),Q(3,289)),(Q(-1,4),Q(3,20)),(Q(-6,11),Q(18,41))]
    expected=[add([p[i],s[i],*(mul(add([a,g[i]]) if j==4 else a,d[j]) for j,a in enumerate(row))])
              for i,row in enumerate(rhs)]
    need(A._coupled_predictor_force(p,rhs,g,s,d,4)==expected,"all five signed shared parameters applied once")
    need(A._coupled_predictor_force([(Q(0),Q(0))],
         [[(Q(0),Q(0))]*4+[(Q(-2),Q(-2))]],[(Q(2),Q(2))],[(Q(1),Q(1))],d,4)==[(Q(1),Q(1))],
         "gain column cancels exactly before interval gain shift")
    f=[D.ZERO]*4
    A._force_scatter(f,[3,1],[(Q(-3),Q(1)),(Q(1),Q(2))])
    A._force_scatter(f,[3,1],[(Q(1),Q(4)),(Q(-5),Q(-2))])
    need(f==[D.ZERO,(Q(-4),Q(0)),D.ZERO,(Q(-2),Q(5))],"signed local mapped accumulation")
    need(D._strict_test(Q(1,8),Q(1,4),Q(1,8),Q(1,4),Q(1,2))[-2:]==(False,None),"strict equality")
    return {"signed_RH_before_absolute_row_bound":True,"shared_five_offsets_and_gain_cancellation":True,
            "signed_scatter_strict_radius_equality_unresolved":True}


def structural(M):
    code={x.name:ast.unparse(x) for x in ast.parse((ARCHIVE/"exact_centered_cohort_reduction_frozen.py").read_bytes()).body if isinstance(x,ast.FunctionDef)}
    c=code["_calculate"]
    need(all(z not in c for z in ("_source_inputs(", "_geometry(", "load_context(", "bulk_component("))
         and "W._factory(" in c and "W._element_guard(" in c,
         "only pure explicit arithmetic helpers reused, no L assets or force work")
    need(c.index("Htotal.update(Hbulk)")<c.index('packets[\'surface\']')<c.index("D._dyadic_reduce(Rmat, f")
         and "D._digest(entries) == D._digest(prior" in c and "norm == Q(prior" in c and "D._minimum(" in c
         and "if ok else None" in c and "defect, b, r" in c and "D._old_premises(ctx, p, hashes)" in c,
         "original 24 surfaces, original full H and z(tube), conditional same-record SPD")
    need(c.index("A._bulk_RH_effects(Rmat, Hbulk")<c.index("A._coupled_predictor_force(")
         <c.index("D._strict_test(") and "ctx.S" in c and "ctx.cell[k][0] - ctx.base[k]" in c,
         "predictor-only RHbulkS; five same-cell offsets, valid only for residual")
    part=code["_terminal_inputs"]
    need("os_exit['pid'] == child['pid'] == terminal['pid']" in part and "type(os_exit['exit_code']) is int" in part
         and "os_exit['exit_code'] == 0" in part and "D._terminal(terminal, pins['manifest_sha256'])" in part
         and "_complete_files(" in part and "/results.json" in part,
         "no manifest-only approval; observed OS exit and worker PID, digest")
    need(code["admit"].index("_terminal_inputs(")<code["admit"].index("D._write(")
         and "_calculate(" not in code["admit"] and code["reduce"].index("_packet_review(")<code["reduce"].index("_calculate("),
         "admission and independent actual review before global reduction")
    return {"manifest_requires_canonical_coverage_and_matching_observed_OS_exit":True,
            "prior_H_b_R_surface_domain_preserved":True,"pure_helpers_explicit_record_scope":True,
            "independent_actual_gate_precedes_calculation_in_source":True}


def fixtures(M,rt,D,A,ready):
    paths=[{"element":i,"file":f"synthetic/components/element_{i:03d}.json","sha256":"3"*64} for i in range(144)]
    links={"record_id":RECORD,"factory_nodeball_sha256":"f"*64}
    manifest={**links,"status":"COMPLETE_144_CENTERED_FORCE_PACKETS_UNREVIEWED","files":paths,
              "new_Hessian_integrations":False,"root_asserted":False}
    need(M._complete_files(D,manifest,"synthetic",links)==paths,"canonical full fake manifest")
    bad={"partial":paths[:-1],"duplicate":paths[:-1]+paths[:1],
         "wrongpath":[{**paths[0],"file":"foreign"+paths[0]["file"]},*paths[1:]]}
    need(all(rejects(lambda:M._complete_files(D,{**manifest,"files":q},"synthetic",links)) for q in bad.values())
         and rejects(lambda:M._complete_files(D,{**manifest,"root_asserted":True},"synthetic",links))
         and rejects(lambda:D._terminal({"status":"FINISHED","pid":True,"exit_code":0,"results_sha256":"a"*64})),
         "invalid completeness, root flag and terminal type rejected")
    pins={"manifest_sha256":"a"*64,"terminal_sha256":"b"*64}
    good={"status":"PASS_CENTERED_PACKET_REPLAY","checks":{"fixture":True},"record_id":RECORD,"components":144,
          "admission_sha256":"c"*64,"centered_manifest_sha256":"a"*64,"centered_terminal_sha256":"b"*64,**PINS}
    with patch.object(D,"_bound",return_value=good):
        M._packet_review(D,RECORD,"d"*64,"c"*64,pins,PINS,{})
    for k,v in (("centered_terminal_sha256","x"),("admission_sha256","x"),("components",143),
                ("status","MANIFEST_ONLY"),("checks",{}),("record_id","other")):
        with patch.object(D,"_bound",return_value={**good,k:v}):
            need(rejects(lambda:M._packet_review(D,RECORD,"d"*64,"c"*64,pins,PINS,{})),
                 "mixed independent actual packet review rejected")
    blob={}
    def capture(file,obj):
        blob.update(obj)
        return "a"*64
    _,_,runtime,_,plan,jobs,expect,_,_=rt
    with (patch.object(M,"_boot",return_value=rt),patch.object(M,"_own_gate",return_value=PINS),
          patch.object(M,"_declared_record",return_value=({},{})),
          patch.object(M,"_terminal_inputs",return_value=(paths,links)),
          patch.object(D,"_write",side_effect=capture),patch.object(D,"_bound",return_value=blob),
          patch.object(M,"_packet_review",side_effect=AssertionError("expected failure before actual packet gate"))):
        M.admit(RECORD,"d"*64,pins)
        need(blob["local_centered_payloads_opened"] is False and blob["expected_packets"]==paths,"actual admit output keys")
        try:
            M.reduce(RECORD,"c"*64,"e"*64)
        except KeyError as defect:
            need(defect.args==("expected_worker_links",) and blob["expected_current_worker_links"]==links,
                 "precise admit/reduce schema defect")
        else:
            raise AssertionError("expected current schema to fail closed before composition")
    return {"canonical_worker_manifest_and_negative_packet_reviews":True,
            "admit_reduce_roundtrip":False}, sorted(blob)


def main():
    d=rel(ARCHIVE)
    archived=read(f"{d}/exact_centered_cohort_reduction_frozen.py",PINS["source_sha256"],False)
    protocol,pre=[read(f"{d}/{n}.json",PINS[n+"_sha256"]) for n in ("protocol","preflight")]
    need(protocol["status"]=="FROZEN_CENTERED_COHORT_REDUCTION_PREPARATION_REVIEW_PENDING"
         and pre["status"]=="PASS_GENERIC_CENTERED_REDUCTION_PREPARATION_ONLY" and yes(pre["checks"])
         and protocol["source_sha256"]==pre["source_sha256"]==PINS["source_sha256"]
         and pre["protocol_sha256"]==PINS["protocol_sha256"]
         and protocol["actual_component_payloads_opened_or_global_sums"] is pre["actual_component_payloads_opened_or_global_sums"] is False,
         "bound prep")
    for p,h in {**protocol["input_hashes"],**protocol["failed_attempts"]}.items():
        read(p,h,False)
    sys.path.insert(0,str(ROOT/"scripts"))
    M=ModuleType("archived_centered_reduction")
    M.__file__=str(ROOT/SOURCE)
    exec(compile(archived,M.__file__,"exec"),M.__dict__)
    rt=M._boot()
    W,H,runtime,A,plan,jobs,expect,hashes,ready=rt
    D=runtime[1].D
    need(rel(M.__file__)==SOURCE and set(W.RECORDS)==set(ready["records"])==set(protocol["currently_eligible_fully_reviewed_records"]),
         "only current C/D source scope")
    for r in W.RECORDS:
        M._worker_gate(D,CALLER_REVIEW,r,hashes)
        M._declared_record(D,ready,r,hashes)
    a=arith(D,A)
    c=structural(M)
    f,keys=fixtures(M,rt,D,A,ready)
    receipt={"status":"BLOCKED_ADMISSION_LINK_SCHEMA_MISMATCH",**PINS,
             "checks":{**a,**c,**f},"no_prelaunch_PASS_gate_issued":True,
             "review_source_sha256":sha(Path(__file__).read_bytes()),"synthetic_admission_keys":keys,
             "finding":"admit writes expected_current_worker_links. reduce and _calculate read expected_worker_links. The synthetic admit->reduce transition raises KeyError before independent packet review.",
             "actual_payloads_opened_or_global_sums":False}
    raw=(json.dumps(receipt,indent=2,allow_nan=False)+"\n").encode()
    with (HERE/"blocked_e5e11f7.json").open("xb") as p:
        p.write(raw)
    print(json.dumps({"status":receipt["status"],"receipt_sha256":sha(raw)}))


if __name__=="__main__":
    main()

