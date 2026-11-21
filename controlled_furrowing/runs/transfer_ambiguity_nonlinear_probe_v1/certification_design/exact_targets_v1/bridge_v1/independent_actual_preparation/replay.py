"""Authenticate and independently replay one completed 168-component singleton."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib
import json
import os
from pathlib import Path
import re
import sys
import time
import traceback

SOURCE=Path(__file__).resolve()
HERE=SOURCE.parent
BRIDGE=HERE.parent
DESIGN=BRIDGE.parents[1]
ROOT=DESIGN.parents[2]
PRODUCER=ROOT/"scripts/transfer_ambiguity_exact_bridge.py"
MATH=HERE/"replay_math.py"
R=Q(1,2**34)
ORDER=tuple((k,i) for k,n in (("bulk",144),("surface",24)) for i in range(n))
DX_SOURCE="scripts/transfer_ambiguity_cached_derivative_extension.py"
DERIV_ADMISSION=DESIGN/"exact_targets_v1/cached_derivative_extension_v1/independent_actual_preparation/bridge_admission.json"
DX_SHA="7ad0b7d4d42efe634123989c6092ccf768cedb92a79244cc7aaa5c10a62f370c"
DERIV_SHA="269229312793ec101cc301086d634dd4de89c62a5050ccfc60f185cf761dbc7d"
SCALED="runs/exact_integration_cohort_v1/L_defect_diagnostic/scaled_chart_v1/scaled_element21.py"
OTHER_PINS={
    DX_SOURCE:DX_SHA,
    "scripts/exact_cohort_integrator.py":"092dc51161e708cdba3791ec47a7f8e843dd0e1afc866b6f65f286b8657d1745",
    "scripts/review_algebraic_root.py":"6af2311d0667a28af263b6fb6237063481aa8e65209ae9722c2b0246a8cff27a",
    "scripts/review_centered_cohort_packets.py":"4f7e6450f57f0c368f3334646df32a9411e7e333402d5c8ec4b05d17397117b3",
    "scripts/exact_cache_consistency.py":"d363567cf2efbf34b2b3b39fb5282dc7dcabced1b3b96a9717fbcf382d5f17d6",
    "scripts/exact_surface_cache_consistency.py":"02c39fa9dfa459a1cc174211a4731dc54966e6aeb12a0133fd9129b100644d2e",
    "scripts/exact_uniform_surface_bridge_B_safe.py":"61cf47d2b1a9e7ee46f9fa680c24b4fdf7b4e31e4889379c56951ab6a7455ee9",
    "scripts/exact_uniform_bulk_defect_B_safe.py":"7c43896c220090025e3bc98969340f98c2b145099905649d12b16bb57d6fac9c",
    "scripts/exact_scaled_cohort_bulk.py":"6639bb2fae3c3892ec44115ff90e27242586744a2900daa16177cc5a091f1526",
    SCALED:"a173b3759fb66fe66e5f3e442c1bcce9a5a185e67515649a94499f82d1d67f60",
}
SCOPE=("Local point F_I-F_Q and separate entire-box H_I-H_Q bridge enclosures. "
       "No root existence, cached derivative extension, actual reduction or energy minimum. "
       "No continuum or source-data premise.")


def name(p):
    return p.relative_to(ROOT).as_posix()


def require(ok, msg):
    if not ok:
        raise ValueError(msg)


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def pinned(p,h):
    require(type(h) is str and re.fullmatch("[0-9a-f]{64}",h),f"missing valid digest {p}")
    raw=p.read_bytes()
    require(sha256(raw).hexdigest()==h,f"wrong digest {p}")
    return json.loads(raw)


def save(p,x):
    with p.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(x,f,allow_nan=False,indent=2)
        f.write("\n")
    return sha(p)


def yes(x):
    return type(x) is dict and bool(x) and all(y is True for y in x.values())


def same_order(components,order):
    require(all(type(x["slot"]) is int for x in components)
            and tuple((x["kind"],x["slot"]) for x in components)==tuple(order), "unexpected component coverage")


def plan():
    pins=OTHER_PINS|{name(SOURCE):sha(SOURCE),name(MATH):sha(MATH),name(DERIV_ADMISSION):DERIV_SHA}
    require(all(sha(ROOT/f)==h for f,h in pins.items()),"cold mathematical primitive pin changed")
    return {"status":"FROZEN_INDEPENDENT_SINGLETON_COMPONENT_REPLAY",
            "source_sha256":sha(SOURCE),"math_source_sha256":sha(MATH),"source_and_premise_hashes":pins,
            "radius_exact":str(R),"scope":SCOPE,"precision_bits":640,
            "counts":{"bulk":144,"surface":24},
            "independence":"Separate exact barycentric and trace cache maps, direct point strain gradients, "
                           "separate weighted cache and variational moment composition, exact containment tests. "
                           "Shared reviewed lower moments/polynomial J/G jets, axis and scaled coordinate transforms, "
                           "third-tensor local-unweighted corrections and outward Arb arithmetic. "
                           "No producer component functions or legacy input admission loaders."}


def prepare():
    protocol=plan()
    sys.path.insert(0,str(ROOT/"scripts"))
    import replay_math as M
    M.E.core.ctx.prec=640
    M.E.core.A.cache_clear()
    M.ad.A.cache_clear()
    q=Q(1,10)
    # A nonlinear surface site with known direct signed first variation:
    # W=K/2 (a*b-1)^2+T*(a*b), c=0. Canonical radial column doubled.
    y=[Q(1),Q(2),Q(0)]
    new=[[Q(2),Q(0),Q(0)],[Q(0),Q(1),Q(0)],[Q(0),Q(0),Q(1)]]
    old=[[Q(1),Q(0),Q(0)],[Q(0),Q(1),Q(0)],[Q(0),Q(0),Q(1)]]
    site=M.project_direct(new,old,y,[0,1,2],(Q(1),q),True)
    expected=[4*(3+q)-2*(1+q),2*(3+q)-(1+q),Q(0)]
    require(all(M.contain(M.ends(b),(x,x)) for b,x in zip(site,expected,strict=True)),"direct point oracle")
    require(M.contain((-2,3),(-1,2)) and not M.contain((-2,1),(-1,2)),"containment not overlap")
    same_order([{"kind":k,"slot":i} for k,i in ORDER],ORDER)
    try:
        same_order([{"kind":"bulk","slot":False}],[("bulk",0)])
    except ValueError:
        pass
    else:
        raise ValueError("Boolean component scope index")
    ph=save(HERE/"protocol.json",protocol)
    fh=save(HERE/"preflight.json",{"status":"PASS","source_sha256":sha(SOURCE),"protocol_sha256":ph,
                    "synthetic_point_projection_and_nonstrict_overlap_rejection":True,
                    "actual_packet_payloads_opened":0,"actual_moment_force_or_H_work":0})
    print(ph,fh,sha(SOURCE),flush=True)


def admitted(path,h,pilot):
    a=pinned(path,h)
    require(a["status"]=="PASS_ROOT_OBSERVED_SINGLETON_BRIDGE_COMPONENT_TERMINALS" and a["pilot_id"]==pilot
            and type(a["observed_exit_code"]) is int and a["observed_exit_code"]==0
            and type(a["unified_exec_session"]) is int and a["unified_exec_session"]>0,
            "external actual successful supervisor observation")
    pdir=BRIDGE/"plans"/pilot
    run=BRIDGE/"actual"/pilot
    paths=a["sha256"]
    required=[run/f"{s}.json" for s in ("started","terminal","manifest")]
    required += [run/k/f"{s}.json" for k in ("bulk","surface")
                  for s in ("started","supervisor_child_started","terminal","observed_os_exit","manifest")]
    required += [pdir/"plan.json",pdir/"started.json",pdir/"terminal.json"]
    required += [ROOT/a[label] for label in ("source_review_file","plan_review_file","producer_protocol_file",
                                               "producer_preflight_file","producer_frozen_source_file")]
    require(all(p.resolve().is_relative_to(BRIDGE.resolve()) and name(p.resolve())==name(p) for p in required),
            "all processes and proofs inside managed bridge archive")
    require(all(name(p) in paths for p in [*required,PRODUCER]),"incomplete external parent pins")
    require(all((ROOT/f).resolve().is_relative_to(ROOT) and name((ROOT/f).resolve())==f
                and sha(ROOT/f)==d for f,d in paths.items()),"producer source or metadata fingerprint changed")
    read=lambda p:pinned(p,paths[name(p)])
    start,end,manifest=(read(run/f"{kind}.json") for kind in ("started","terminal","manifest"))
    plan,pstart,pend=(read(pdir/f"{kind}.json") for kind in ("plan","started","terminal"))
    sr,pr,prod_protocol,prod_preflight=[read(ROOT/a[label]) for label in (
              "source_review_file","plan_review_file","producer_protocol_file","producer_preflight_file")]
    links=start["bindings"]
    require(start["pid"]==end["pid"]==a["supervisor_pid"] and end["status"]=="FINISHED"
            and all(x["pilot_id"]==pilot for x in (start,end,manifest))
            and end["manifest_sha256"]==paths[name(run/"manifest.json")]
            and start["bindings"]==end["bindings"]==manifest["bindings"],"same complete supervisor phase")
    require(manifest["status"]=="COMPLETE_168_SINGLETON_BRIDGE_PACKETS_UNREVIEWED"
            and manifest["independent_actual_component_review_required"] is True and manifest["root_claimed"] is False,
            "complete component-only packet manifest required")
    require(paths[name(PRODUCER)]==paths[a["producer_frozen_source_file"]]==links["source_sha256"]
            and paths[a["producer_protocol_file"]]==links["protocol_sha256"]
            and paths[a["producer_preflight_file"]]==links["preflight_sha256"]
            and paths[a["source_review_file"]]==links["source_prelaunch_review_sha256"]
            and paths[a["plan_review_file"]]==links["plan_review_sha256"],
            "actual producer source/protocol/preflight and both independent gates")
    require(sr["status"]=="PASS" and yes(sr["checks"])
            and all(sr[k]==links[k] for k in ("source_sha256","protocol_sha256","preflight_sha256"))
            and prod_protocol["source_sha256"]==links["source_sha256"]
            and prod_preflight["status"]=="PASS" and yes(prod_preflight["checks"]),"valid source preparation and gate")
    require(pr["status"]=="PASS_INDEPENDENT_ACTUAL_SINGLETON_BRIDGE_PLAN" and yes(pr["checks"])
            and pr["pilot_id"]==pilot and type(pr["plan_observed_os_exit_code"]) is int
            and pr["plan_observed_os_exit_code"]==0
            and pr["plan_observed_pid"]==pstart["pid"]==pend["pid"]
            and pr["source_sha256"]==links["source_sha256"]
            and pr["plan_terminal_sha256"]==paths[name(pdir/"terminal.json")]
            and pr["plan_sha256"]==paths[name(pdir/"plan.json")]==links["plan_sha256"]
            and pr["derivative_admission_sha256"]==links["derivative_admission_sha256"]==DERIV_SHA,
            "same current independent selector plan")
    require(all(p[k]==plan[k] for p in (pstart,pend) for k in ("pilot_id",))
            and pend["status"]=="FINISHED" and type(pend["exit_code"]) is int and pend["exit_code"]==0
            and pend["plan_sha256"]==links["plan_sha256"]
            and pstart["bindings"]==plan["bindings"]
            and all(links[k]==v for k,v in plan["bindings"].items())
            and plan["status"]=="COMPLETE_POINT_FORCE_WHOLE_H_SELECTOR_PLAN_UNREVIEWED",
            "same completed point/whole selector job")
    same_order(manifest["files"],ORDER)
    same_order(plan["components"],ORDER)
    out=[]
    require([x["part"] for x in end["workers"]]==["bulk","surface"]
            and len({x["pid"] for x in end["workers"]})==2,"two serial child processes")
    for w in end["workers"]:
        part=w["part"]
        d=run/part
        e,child,childstart,parentstart,partmanifest=[read(d/f"{s}.json") for s in (
              "observed_os_exit","terminal","started","supervisor_child_started","manifest")]
        require(type(e["observed_exit_code"]) is int and e["observed_exit_code"]==0 and child["status"]=="FINISHED"
                and e["terminal_sha256"]==paths[name(d/"terminal.json")]
                and w==e|{"os_exit_receipt_sha256":paths[name(d/"observed_os_exit.json")]}
                and child["manifest_sha256"]==paths[name(d/"manifest.json")],"true bound child exit and manifest")
        identity={"pid":e["pid"],"supervisor_pid":start["pid"],"pilot_id":pilot,"bindings":links,"part":part}
        require(all(all(doc[k]==v for k,v in identity.items()) for doc in (child,childstart,partmanifest))
                and parentstart["pid"]==e["pid"] and parentstart["supervisor_pid"]==start["pid"],
                "same worker process identity")
        same_order(partmanifest["files"],[(k,i) for k,i in ORDER if k==part])
        out+=partmanifest["files"]
    require(out==manifest["files"],"complete union of the bulk and surface terminal manifests")
    # Authenticate all actual and planned packet bytes before any individual payload is parsed.
    packet_hashes={}
    proof_hashes={}
    for f,p in zip(manifest["files"],plan["components"],strict=True):
        k,i=f["kind"],f["slot"]
        rel=f"{k}_{i:03d}.json"
        require(p["file"]==rel and f["file"]==name(run/k/rel),"bounded stage-specific packet path")
        require(sha(ROOT/f["file"])==f["sha256"] and sha(pdir/rel)==p["sha256"],"packet or plan byte changed")
        packet_hashes[f["file"]]=f["sha256"]
        proof_hashes[name(pdir/rel)]=p["sha256"]
    return a,plan,links,manifest,packet_hashes,proof_hashes


def run(args):
    pre=pinned(HERE/"preflight.json",args.preflight_sha256)
    ph=pinned(HERE/"protocol.json",args.protocol_sha256)
    require(ph==plan() and pre["status"]=="PASS" and pre["source_sha256"]==sha(SOURCE)
            and pre["protocol_sha256"]==args.protocol_sha256
            and pre["actual_packet_payloads_opened"]==pre["actual_moment_force_or_H_work"]==0,"unchanged reviewer phase")
    a,p,links,manifest,packet_hashes,proof_hashes=admitted(args.admission.resolve(),args.admission_sha256,args.pilot_id)
    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_cached_derivative_extension as DX
    import transfer_ambiguity_exact_targets as EX
    import replay_math as M
    scaled=importlib.import_module("exact_scaled_cohort_bulk").S
    require(name(Path(scaled.__file__).resolve())==SCALED and M.E is EX.E,"same explicit math aliases")
    out=HERE/"actual_replay"/args.pilot_id
    out.mkdir(parents=True,exist_ok=False)
    mark={"pid":os.getpid(),"source_sha256":sha(SOURCE),"protocol_sha256":args.protocol_sha256,
          "preflight_sha256":args.preflight_sha256,"root_observation_sha256":args.admission_sha256}
    save(out/"started.json",mark|{"status":"RUNNING"})
    terminal=mark|{"status":"FAILED","exit_code":1}
    t0=time.perf_counter()
    try:
        _,cases,man,data=DX.authority()
        case,=[x for x in cases if x["pilot_id"]==args.pilot_id]
        ctx,rule,entry=EX.open_context(args.pilot_id,man,data)
        box=pinned(ROOT/case["box"],case["box_sha256"])
        require(box==EX.exact_node_box(ctx,R) and p["context_box_sha256"]==case["box_sha256"]==entry["box_1_sha256"]
                and ctx.record_id==p["record_id"],"exact current factory box and physical law")
        proof=pinned(ROOT/case["verified_point"],case["verified_point_sha256"])
        require(proof["status"]=="PASS_INDEPENDENT_CACHED_LOCAL_ROOT" and proof["exact_strength"]==str(ctx.drive)
                and case["preconditioner_binary_sha256"]==sha256(ctx.R.tobytes(order="C")).hexdigest(),
                "same exact-requested-drive anchor and binary R")
        deriv=pinned(DERIV_ADMISSION,DERIV_SHA)
        dr,=[x for x in deriv["reviewed_cases"] if x["pilot_id"]==args.pilot_id]
        zn=pinned(ROOT/dr["native_case_file"],dr["native_case_sha256"])
        z=max(Q(zn["reviewed_ZQ_upper_exact"]),Q(zn["producer_ZQ_upper_exact"]))
        require(Q(p["derivative"]["ZQ_upper_exact"])==z and z<1
                and dr["box_sha256"]==case["box_sha256"] and dr["cache_packet_sha256"]==case["cache_packet_sha256"]
                and Q(dr["radius_exact"])==R,"same admitted cached whole-box ZQ")
        tails=tuple(map(Q,p["tail_caps_exact"]))
        require(tails==tuple(Q(p["derivative"]["planning"][k]) for k in (
           "per_component_force_entry_tail_cap_exact","per_component_H_entry_tail_cap_exact")),"same selector tail capsule")
        geometry=pinned(ROOT/case["geometry"],case["geometry_sha256"])
        body=pinned((ROOT/case["geometry"]).with_name("body_geometry.json"),geometry["body_geometry_sha256"])
        require(body["status"]=="CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING","same actual spatial geometry")
        route=pinned(ROOT/entry["routes"],entry["routes_sha256"])
        ctxnodes=ctx.target.nodes()
        nodes=M.SUR.common_nodes(ctx.target,box)
        # Mark the center with zero radii independently and match the planned hash.
        zeros=[["0","0"] for _ in ctxnodes]
        exactpoint={**box,"coordinate_radii_exact":zeros,"coordinate_hulls_exact":[
                [[str(x),str(x)] for x in row] for row in ctxnodes],"corrector_radius_exact":"0","scope":M.POINT}
        require(M.E._digest(exactpoint)==p["point_box_sha256"]
                and Q(p["whole_box_radius_exact"])==R,"point is the same exact X+u and derivative has full radius")
        rows=[]
        any_failure=False
        for x,y in zip(manifest["files"],p["components"],strict=True):
            k,i=x["kind"],x["slot"]
            pf=BRIDGE/"plans"/args.pilot_id/y["file"]
            pfdoc=pinned(pf,y["sha256"])
            actual=pinned(ROOT/x["file"],x["sha256"])
            require(actual["bindings"]==links and actual["plan_file_sha256"]==y["sha256"]
                    and actual["record_id"]==ctx.record_id and actual["precision_bits"]==640
                    and actual["proof_sha256"]==M.E._digest(pfdoc),"same packet/source/plan element proof identity")
            require(route[k][i]["global_indices"]==pfdoc["global_indices"],"independently admitted geometry map")
            M.ad.A.cache_clear()
            M.E.core.A.cache_clear()
            local=M.replay_component(ctx,rule,ctxnodes,nodes,body,i,k,pfdoc,actual,tails,scaled)
            paths={"kind":k,"slot":i,"global_indices":local["global_indices"],"point_force_scope":M.POINT,
                   "whole_box_H_scope":M.BOX,"plan_file_sha256":y["sha256"],"packet_file":x["file"],"packet_sha256":x["sha256"]}
            local|=paths
            report=out/x["file"].split("/")[-1]
            digest=save(report,local)
            rows.append(paths|{"replay_file":name(report),"replay_sha256":digest})
            any_failure |= any(local["containment_failures"].values())
            print(json.dumps({"kind":k,"slot":i,"containment":not any(local["containment_failures"].values()),
                              "seconds":local["seconds"]}),flush=True)
        require(len(rows)==168 and all(sha(ROOT/f)==h for f,h in
                    (packet_hashes|proof_hashes|a["sha256"]).items()),"all168 complete and unmodified")
        rec={"status":"HOLD_INDEPENDENT_COMPONENT_CONTAINMENT" if any_failure else
             "PASS_INDEPENDENT_ACTUAL_SINGLETON_BRIDGE_COMPONENT_REPLAY",
             "checks":{"true_os_complete_producer_and_both_child_terminals":True,
                       "same_exact_point_and_whole_box_all_168_bindings":True,
                       "all_point_cache_force_reconstructed_and_contained":not any_failure,
                       "all_point_force_and_whole_H_contained":not any_failure},
             "scope":SCOPE,"pilot_id":args.pilot_id,"record_id":ctx.record_id,**mark,
             "producer_source_sha256":links["source_sha256"],
             "producer_manifest_sha256":a["sha256"][name(BRIDGE/"actual"/args.pilot_id/"manifest.json")],
             "producer_terminal_sha256":a["sha256"][name(BRIDGE/"actual"/args.pilot_id/"terminal.json")],
             "producer_observed_os_exit_code":a["observed_exit_code"],"producer_observed_pid":a["supervisor_pid"],
             "producer_observed_exec_session":a["unified_exec_session"],
             "plan_sha256":links["plan_sha256"],"plan_review_sha256":links["plan_review_sha256"],
             "derivative_admission_sha256":DERIV_SHA,"box_sha256":case["box_sha256"],
             "parameter_cell_sha256":ctx.row["parameter_cell_sha256"],"cache_packet_sha256":case["cache_packet_sha256"],
             "binary_preconditioner_sha256":case["preconditioner_binary_sha256"],
             "signed_cached_anchor_file":case["verified_point"],"signed_cached_anchor_sha256":case["verified_point_sha256"],
             "reviewed_cached_ZQ_upper_exact":str(z),"radius_exact":str(R),"checked_packet_sha256":packet_hashes,
             "plan_packet_sha256":proof_hashes,"components":rows,"seconds":time.perf_counter()-t0}
        digest=save(out/"receipt.json",rec)
        terminal|={"status":"FINISHED_WITH_HOLD" if any_failure else "FINISHED","exit_code":2 if any_failure else 0,
                   "receipt_sha256":digest}
        if any_failure:
            return 2
        return 0
    except BaseException as error:
        save(out/"failure.json",{"error":str(error),"traceback":traceback.format_exc()})
        raise
    finally:
        save(out/"terminal.json",terminal)
        print(json.dumps(terminal),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    sub=parser.add_subparsers(dest="action",required=True)
    sub.add_parser("prepare")
    act=sub.add_parser("run")
    act.add_argument("--admission",type=Path,required=True)
    act.add_argument("--admission-sha256",required=True)
    act.add_argument("--protocol-sha256",required=True)
    act.add_argument("--preflight-sha256",required=True)
    act.add_argument("--pilot-id",required=True)
    a=parser.parse_args()
    if a.action=="prepare":
        prepare()
    else:
        require(re.fullmatch("[a-z0-9_]+",a.pilot_id),"known singleton ID form")
        sys.exit(run(a))
