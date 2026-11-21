"""Independent centered supervisor prelaunch review and mocked OS tests."""
import ast
from dataclasses import replace
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
from unittest.mock import Mock,patch

from flint import arb,fmpq,ctx as precision

OUT=Path(__file__).resolve().parent
BATCH=OUT.parent
ROOT=OUT.parents[3]
BASE=BATCH.parent.relative_to(ROOT).as_posix()
PREFIX=BATCH.relative_to(ROOT).as_posix()
SOURCE="scripts/remaining_cohort_centered.py"
PINS={"source_sha256":"bd384270f40fd08962d1c82c33ce5203af7c3f40c79654b8275801ffa22e397c",
      "protocol_sha256":"1a4ffdd07db4b1073efc2715eb6a548f2011da9d44e286d115713fe7f1d2b845",
      "preflight_sha256":"3b10c934c20f6c802355257b7cc710663b873c34f39201816bce8ca4ac88d77f",
      "manifest_sha256":"8cf34d5a9fded6cc6c5a761d57e157da292c315d871bf34c2ec352d707b9ce65"}
HELPER_REVIEW="ffcebc86827a4e12cdbec8e0a7fd2b5be84e884eb83f88078d9e61356e7b3769"
RECORDS=tuple("optoGEF2_small_col_"+r for r in "CD")


def check(test,why):
    if not test:
        raise ValueError(why)


def sha(raw):
    return sha256(raw).hexdigest()


def stable(x):
    return sha(json.dumps(x,default=str,sort_keys=True,separators=(",",":")).encode())


def rel(p):
    return Path(p).resolve().relative_to(ROOT).as_posix()


def bound(p,h,parse=True):
    local=(ROOT/p).resolve()
    check(local.is_relative_to(ROOT) and rel(local)==p,"canonical study location")
    raw=local.read_bytes()
    check(sha(raw)==h,f"changed input {p}")
    return json.loads(raw) if parse else raw


def yes(x):
    return isinstance(x,dict) and bool(x) and all(v is True for v in x.values())


def reject(func):
    try:
        func()
    except ValueError:
        return True
    return False


def box_independent(ctx, box):
    # Build each original coordinate using the binary arrays and exact affine extrema.
    r=ctx.old_radius
    order=("apical_area_modulus","basal_area_modulus","apical_surface_tension","basal_surface_tension","common_gain")
    delta=[tuple(x-ctx.base[k] for x in ctx.cell[k]) for k in order]
    free={int(i):j for j,i in enumerate(ctx.free)}
    v=[]
    for j in range(325):
        row=[]
        for k in (0,1):
            parent=Q.from_float(float(ctx.target.a["reference"][j,k]))
            pos=free.get(2*j+k)
            if pos is None:
                y,s=parent,Q(0)
            else:
                y=parent+Q.from_float(float(ctx.u[pos]))
                prods=[sorted(Q.from_float(float(t))*v for v in d) for t,d in zip(ctx.S[pos],delta)]
                lo,hi=[sum(c) for c in zip(*prods)]
                s=r+max(-lo,hi)
            check(y == Q(box["current_centers_exact"][j][k]) and s == Q(box["coordinate_radii_exact"][j][k])
                  and parent==Q(box["reference_nodes_exact"][j][k]),"exact reconstructed 650-coordinate box")
            c,d=y-s,y+s
            qlo=arb(fmpq(c.numerator,c.denominator))
            qhi=arb(fmpq(d.numerator,d.denominator))
            ball=qlo.union(qhi)
            lower,upper=hull(ball)
            check(Q(lower)<=c<=d<=Q(upper) and (k!=0 or parent!=0 or lower==upper=="0"),
                  "original force/H domain, outward 640-bit endpoints and exact axis")
            row.append([lower,upper])
        v.append(row)
    return v


def hull(b):
    return str(b.lower().fmpq()),str(b.upper().fmpq())


def runtime_hulls(nodes):
    return [[list(hull(x)) for x in n] for n in nodes]


def source_order():
    source=bound(SOURCE,PINS["source_sha256"],False)
    f={x.name:ast.unparse(x) for x in ast.parse(source).body if isinstance(x,ast.FunctionDef)}
    factory=f["_factory"]
    check(factory.index("precision.prec = 640") < factory.index("op.E.core.A.cache_clear()")
          < factory.index("H.bind_reviewed_record(") < factory.index("common_nodes("),
          "precision initialization before domain import and convex holder")
    leaf=f["_element_guard"]
    check(leaf.index("_node_guard(") < leaf.index("D._bound(raw['file']") < leaf.index("H.prepare_component("),
          "per-element independent route identity before prior packet opens")
    worker=f["worker"]
    check(worker.index("_gate(")<worker.index("_factory(")<worker.index("exist_ok=False")
          <worker.index("_element_guard(")<worker.index("H.bulk_component(")
          and 'receipt == expected' in worker and "'root_asserted': False" in worker,
          "nonduplicate execution, source receipt matches, no root inference")
    launch=f["launch"]
    check(launch.index("_gate(")<launch.index("subprocess.Popen(")<launch.index("child.wait()")
          <launch.index(".os_terminal.json")<launch.index('require(code == 0')
          and "terminal['pid'] == child.pid" in launch and "D._terminal(terminal, rsha)" in launch,
          "explicit gate and authentic process exit")
    for k in ("_boot","_factory","_element_guard","worker","launch"):
        check(all(marker not in f[k] for marker in ("_geometry(", "_source_inputs(", "load_context(")),
              "no alternative L/imported data route")


def replay_all(A,H,runtime,plan,jobs,expect,hashes,manifest):
    _,op,_=runtime
    D=op.D
    protected=AssertionError("new centered work in preparation-only review")
    stats={}
    with (patch.object(H,"bulk_component",side_effect=protected),
          patch.object(op.F,"_force_parts",side_effect=protected),
          patch.object(op.E.U.C,"exact_bulk_maps",side_effect=protected),
          patch.object(op.E.M,"triangle_moment",side_effect=protected),
          patch.object(A.subprocess,"Popen",side_effect=protected)):
        for record in RECORDS:
            m=manifest["records"][record]
            p=m["helper_proof_file"]
            check(p==f"{A.PREPARED}/{record}_proofs.json" and m["record_id"]==record and A.PIN[p]==m["helper_proof_sha256"],
                  "same helper proof pin")
            proof=bound(p,A.PIN[p])
            calls=[]
            bind=H.bind_reviewed_record

            def captured(*args,**kw):
                check(precision.prec==640 and op.E.core.A.cache_info().currsize==0,"precision was set and A cleared before binder")
                result=bind(*args,**kw)
                calls.append(result[0])
                return result

            precision.prec=192
            op.E.core.A(Q(123457,1217))
            with patch.object(H,"bind_reviewed_record",side_effect=captured):
                factory=A._factory(H,runtime,plan,jobs,expect,record,proof,hashes)
            dom,actual,routes,pack,canonical,key=factory
            ref=box_independent(dom.context,dom.box)
            check(precision.prec==640 and calls==[dom] and calls[0] is dom
                  and ref==runtime_hulls(dom.nodes)==runtime_hulls(canonical) and stable(ref)==key==m["factory_nodeball_sha256"],
                  "actual binder reference preserved and all independent rational endpoints match")
            check(proof["retained_surface_manifest"]==pack["surface"]==m["retained_surface_manifest"]
                  and proof["retained_upper_Hessian_entries"]==m["unchanged_full_H_upper_entries"]==6676
                  and proof["retained_fullbox_b_exact"]==m["unchanged_fullbox_b_exact"]
                  and m["operator_review_sha256"]==proof["actual_direct_review_sha256"]
                  and len(m["elements"])==len(pack["bulk"])==144 and len(pack["surface"])==24,
                  "same actual terminal leaf and retained operators")
            for i,(entry,prior,route,leaf) in enumerate(zip(m["elements"],pack["bulk"],routes["bulk"],proof["bulk_proofs"],strict=True)):
                raw=bound(prior["file"],prior["sha256"])
                nodes=list(map(int,dom.context.target.a["triangles"][i]))
                check(raw["record_id"]==dom.context.record_id==record==entry["record_id"]
                      and nodes==route["local_nodes"]==leaf["local_nodes"]
                      and raw["nodal_hulls_sha256"]==route["nodal_hulls_sha256"]==stable([ref[j] for j in nodes])
                      and entry["nodeball_sha256"]==raw["nodal_hulls_sha256"]
                      and entry["prior_uniform_file"]==prior["file"] and entry["prior_uniform_sha256"]==prior["sha256"]
                      and entry["centered_proof_sha256"]==stable(leaf) and entry["original_route_sha256"]==stable(route)
                      and entry["chart"]==leaf["chart"] and entry["point_force_tail_upper_exact"]==leaf["point_force_tail_max_exact"],
                      "local canonical nodeball and prior scoped manifest hashes across all cases")
                recheck,_,_,_,_=A._element_guard(H,runtime,dom,canonical,actual,routes,pack,proof,i,hashes)
                check(recheck==entry,"worker guard replay agrees with current frozen all-case manifest")
            # Test independent failures with actual, previously reviewed protocol and proof.
            ii=6
            first=proof["bulk_proofs"][ii]["local_nodes"][0]
            enlarged=[r.copy() for r in dom.nodes]
            enlarged[first][0]=enlarged[first][0]+arb(0,"1/1024")
            changed=replace(dom,nodes=enlarged)
            original=D._bound
            def deny_payload(file,*args,**kwargs):
                if "/components/" in file:
                    raise protected
                return original(file,*args,**kwargs)
            with patch.object(D,"_bound",side_effect=deny_payload):
                check(reject(lambda:A._element_guard(H,runtime,changed,canonical,actual,routes,pack,proof,ii,{})),
                      "same midpoint enlarged box rejected before packet opens")
                foreign=replace(dom,context=replace(dom.context,record_id=RECORDS[1] if record==RECORDS[0] else RECORDS[0]))
                check(reject(lambda:A._element_guard(H,runtime,foreign,canonical,actual,routes,pack,proof,ii,{})),
                      "foreign record rejected before packet opens")
                precision.prec=320
                try:
                    check(reject(lambda:A._element_guard(H,runtime,dom,canonical,actual,routes,pack,proof,ii,{})),
                          "wrong execution precision rejected before proof use")
                finally:
                    precision.prec=640
            loose=(changed,actual,routes,pack,{"upper_Hessian_entries":6676,"b_upper_exact":m["unchanged_fullbox_b_exact"]})
            with patch.object(H,"bind_reviewed_record",return_value=loose):
                check(reject(lambda:A._factory(H,runtime,plan,jobs,expect,record,proof,{})),
                      "full factory enlargement rejected despite identical exact center")
            stats[record]={"nodal_coordinates":650,"local_ball_and_route_proofs":144,
                           "retained_prior_surface_receipts":24,"realized_640bit_factory_sha256":key}
            print(record,"checked",flush=True)
    return stats


def fake_children(A,D,runtime):
    loc=OUT/"supervisor_fixtures"
    loc.mkdir(exist_ok=False)
    for tag,returncode,mismatch in (("success",0,False),("failure",7,False),("wrong_pid",0,True)):
        folder=loc/tag
        folder.mkdir()
        def spawned(args,**kw):
            check(gate.called and args[4]=="worker" and args[-1]==RECORDS[0] and kw["cwd"]==ROOT,
                  "gate observes one explicit worker before fake spawn")
            p=Mock(pid=55)
            p.wait.return_value=returncode
            if returncode==0:
                d=folder/"actual"/RECORDS[0]
                d.mkdir(parents=True)
                (d/"results.json").write_bytes(b"{}")
                (d/"terminal.json").write_text(json.dumps({"status":"FINISHED","exit_code":0,"pid":55+int(mismatch),
                       "results_sha256":sha(b"{}")}),encoding="utf8")
            children.append(p)
            return p
        children=[]
        with (patch.object(A,"OUT",folder),patch.object(A,"_boot",return_value=(None,runtime,None,None,None,None)),
              patch.object(A,"_gate",return_value=({},{})) as gate,patch.object(A.subprocess,"Popen",side_effect=spawned)):
            failed=reject(lambda:A.launch("a"*64,"b"*64,[RECORDS[0]]))
        success=returncode==0 and not mismatch
        check(len(children)==1 and children[0].wait.call_count==1 and failed is not success,
              "mocked wait and terminal/PID outcomes")
        work,=(folder/"launches").iterdir()
        os_exit=json.loads((work/f"{RECORDS[0]}.os_terminal.json").read_bytes())
        check(os_exit["pid"]==55 and os_exit["exit_code"]==returncode
              and (work/"failure.json").exists() is not success,"true OS evidence and failed artifacts preserved")


def main():
    bound(f"{PREFIX}/remaining_cohort_centered_frozen.py",PINS["source_sha256"],False)
    protocol,pre,manifest=[bound(f"{PREFIX}/{n}.json",PINS[f"{n}_sha256"]) for n in ("protocol","preflight","manifest")]
    check(protocol["status"]=="FROZEN_CENTERED_SUPERVISOR_PREPARATION_REVIEW_PENDING"
          and pre["status"]=="PASS_CENTERED_ACTUAL_CALLER_PREPARATION_ONLY" and yes(pre["checks"])
          and manifest["status"]=="READY_FOR_CENTERED_LAUNCHER_INDEPENDENT_PRELAUNCH"
          and protocol["records"]==list(RECORDS)==list(manifest["records"])
          and all(all(v[k]==PINS[k] for k in v if k in PINS)
                  and v["actual_centered_force_or_child_launch"] is False for v in (protocol,pre,manifest)),
          "bound preparation only, no centered executions")
    for name,h in protocol["input_hashes"].items():
        if "/components/" in name:
            check(any(f"/records/{r}/scaled_integrals_v1/actual/" in name for r in RECORDS),
                  "prior read only two terminal independently reviewed records")
        bound(name,h,False)
    source_order()
    sys.path.insert(0,str(ROOT/"scripts"))
    import remaining_cohort_centered as A
    H,rt,plan,jobs,expect,hashes=A._boot()
    D=rt[1].D
    check(rel(A.__file__)==SOURCE and A.RECORDS==RECORDS and A.PIN[A.PREPARED+"/independent_review.json"]==HELPER_REVIEW,
          "current supervisor, dependency and allowed scope")
    stats=replay_all(A,H,rt,plan,jobs,expect,hashes,manifest)
    fake_children(A,D,rt)
    check(not (BATCH/"actual").exists(),"no real centered processes, outputs or packets")
    checks=dict.fromkeys(("frozen_sources_and_manifest","precision_640_and_cache_clear_before_authentic_factory",
                           "independent_1300_material_coordinate_hulls_and_288_local_ball_replays",
                           "authentic_factory_object_and_bound_all_centered_proofs",
                           "only_reviewed_C_and_D_prior_uniform_packets","retained_surface_full_H_and_b",
                           "oversized_wrong_record_and_precision_rejected_before_prior_packet_open",
                           "worker_gate_before_child_wait_and_exact_OS_exit_and_PID",
                           "mocked_nonzero_exit_and_mismatch_fail_closed"),True)
    result={"status":"PASS",**PINS,"helper_review_sha256":HELPER_REVIEW,"records":list(RECORDS),
            "checks":checks,"source_review_sha256":sha(Path(__file__).read_bytes()),"coverage":stats,
            "actual_centered_force_or_child_launch_by_review":False,
            "scope":"Prelaunch review only for C/D original centered bulk point force jobs. Bind current supervisor and actual independent reviewed prior domains; all actual outputs remain unreviewed until admitted against externally pinned worker OS exit, matching terminal and complete 144 packet manifest.",
            "execution":"One child per supervisor at a time. Respect active process cap through the caller. Existing per-record outputs are exclusive. Full original 24 surfaces and 6676 Hessian entries stay in the subsequent signed reduction."}
    gate=BATCH/"prelaunch_independent_review.json"
    data=(json.dumps(result,indent=2,allow_nan=False)+"\n").encode()
    digest=sha(data)
    read=D._bound

    def staged(p,h,use,parse=True):
        if p==rel(gate):
            check(h==digest and parse,"exact staged gate pin")
            return result
        return read(p,h,use,parse=parse)
    with patch.object(D,"_bound",side_effect=staged):
        m,links=A._gate(D,PINS["manifest_sha256"],digest,list(RECORDS))
        check(m==manifest and links["prelaunch_sha256"]==digest
              and all(reject(lambda:A._gate(D,PINS["manifest_sha256"],digest,rows))
                      for rows in ([RECORDS[0]]*2,["optoGEF2_small_col_E"],[])),
              "explicit record inclusion/duplicate/exclusion gate readback")
    with gate.open("xb") as dst:
        dst.write(data)
    A._gate(D,PINS["manifest_sha256"],digest,list(RECORDS))
    print(json.dumps({"status":"PASS","gate_sha256":digest,"records":len(stats),"local_proofs":288}))


if __name__=="__main__":
    main()

