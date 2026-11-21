"""First-wave raw spatial packets to independently certified hulls."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
import os
from pathlib import Path
import sys
from time import monotonic
import traceback
from types import SimpleNamespace
from unittest.mock import patch

import hull_math as hm

HERE=Path(__file__).resolve().parent
BR=HERE.parent
EX=BR.parent
ROOT=next(p for p in EX.parents if (p/"PLAN.md").is_file())
SRC=Path(__file__).resolve()
WORKER=ROOT/"scripts/transfer_ambiguity_remaining32_spatial_worker.py"
BRIDGE=ROOT/"scripts/transfer_ambiguity_remaining32_exact_bridge.py"
OLD=EX/"bridge_v1/independent_actual_preparation/replay_math.py"
RAW=BR/"actual_first_wave"
WP=BR/"worker_preparation"
GATE=WP/"independent_prelaunch/review.json"
ADMISSION=BR/"admission/first_wave_worker_admission.json"
PINS={
    WORKER:"61b0082a68a297f3ae287232704a9b72ed4fbf53e96090db990ba0d8a8458f2d",
    BRIDGE:"4d52e8fb8dbc01ba58c38052bbff537443693234499a9a3b9b3e138c82d9b51b",
    OLD:"e7bb66066d54e131fd009d2bf48e94f3f1fb7bc1955023f66bc96cbe3c6d78e9",
    WP/"protocol.json":"434c3786447d92c75bd13880d57ad090cbc36c2aff5b795515ceed526d110e42",
    WP/"preflight.json":"2f9c6d5e1f3c88f50389774e6afc67b1e3ed18030c6b3f12af7804cb4f41d918",
    GATE:"5463da7c9c36a30007513dce4788f7114878978bb7274d2e484b470c6916e85b",
    ADMISSION:"8468da6cd555797298df182d6213a93a4d054029bdfca9b332aea80ca39f6855",
    BR/"admission/root_hull_design_review.json":"6442253147ebc1a73dff1358f8a44ebba8425ddb3e0af18ba01fb3b6a8fefc03",
    BR/"admission/INDEPENDENT-HULL-DESIGN.md":"bfd06186f1c6637d004e420ba595f6b61ceb3397ab8914b7155e7958efc40e59",
}
FIRST=("max_contrast_change_small_01","max_contrast_change_small_13",
       "max_negative_contrast_change_small_01","max_negative_contrast_change_small_16")
COUNTS={"bulk":144,"surface":24}
JOBS=tuple((case,part) for case in FIRST for part in COUNTS)
STAGE="remaining32_spatial_first_wave"
POINT="EXACT_CENTER_POINT_ONLY"
BOX="ENTIRE_CLOSED_FIXED_PARAMETER_CUBE"
RADIUS=Q(1,2**34)
SCOPE=("Same-target independently enclose exact-integration minus frozen cached force at the exact point "
       "and Hessian over the matched closed cube. Exact hulls are valid through the independent enclosures. "
       "Disjointness anywhere is retained as HOLD, nesting misses are retained diagnostics. No root claim.")


def require(ok,why):
    if not ok:
        raise ValueError(why)


def sha(raw):
    return sha256(raw).hexdigest()


def name(path):
    return path.resolve().relative_to(ROOT).as_posix()


def read(path,want,checked,parse=True):
    p=(ROOT/path).resolve()
    k=name(p)
    require(not isinstance(path,str) or k==path,"wrong canonical file path")
    raw=p.read_bytes()
    h=sha(raw)
    require(h==want==checked.get(k,h),f"compiled evidence changed: {k}")
    checked[k]=h
    return json.loads(raw) if parse else raw


def save(path,doc):
    with path.open("x",encoding="utf-8",newline="\n") as w:
        json.dump(doc,w,indent=2)
        w.write("\n")
    return sha(path.read_bytes())


def yes(x):
    return type(x) is dict and bool(x) and all(v is True for v in x.values())


def zero(x):
    return type(x) is int and x==0


def pid(x):
    return type(x) is int and x>0


def cold(checked):
    p={k:read(k,h,checked,k.suffix==".json") for k,h in PINS.items()}
    w=p[WP/"protocol.json"]
    a=p[ADMISSION]
    gate=p[GATE]
    require(w["status"]=="FROZEN_RAW_EXACT_SPATIAL_FIRST_WAVE_WORKER"
            and a["status"]=="PASS_INDEPENDENT_REMAINING32_FIRST_WAVE_SPATIAL_ADMISSION"
            and gate["status"]=="PASS_INDEPENDENT_REMAINING32_FIRST_WAVE_SPATIAL_WORKER_COLD"
            and yes(gate["checks"]) and yes(a["checks"]) and a["no_later_waves_authorized"] is True
            and a["independent_packet_replay_required"] is True and tuple(w["first_wave"])==FIRST
            and tuple(x["case_id"] for x in a["plans"])==FIRST and a["targets"]==w["targets"]
            and w["jobs"]==[list(j) for j in JOBS],"wrong exact-cube first-wave authority")
    return a


def full_actual(outer_hash,worker_bindings,admission,checked):
    ext=read(RAW/"root_terminal_admission.json",outer_hash,checked)
    require(ext["status"]=="PASS_ROOT_OBSERVED_REMAINING32_FIRST_WAVE_SPATIAL_TERMINALS"
            and pid(ext["pid"]) and pid(ext["unified_exec_session"])
            and zero(ext["observed_exit_code"]),"external real complete producer process required")
    md=("started","terminal","manifest")
    sub=("launched","started","terminal","observed_os_exit","manifest")
    files=ext["sha256"]
    expected={name(RAW/(k+".json")) for k in md}
    for case,part in JOBS:
        folder=RAW/"cases"/case/part
        expected|={name(folder/(k+".json")) for k in sub}
        expected|={name(folder/f"{part}_{i:03d}.json") for i in range(COUNTS[part])}
    require(set(files)==expected and len(expected)==715,"all metadata and all 672 raw packet SHA pins required")
    def seal(path,parse=True):
        return read(path,files[name(path)],checked,parse)
    start,end,batch=(seal(RAW/(k+".json")) for k in md)
    require(start["status"]=="RUNNING" and end["status"]=="FINISHED" and zero(end["exit_code"])
            and end["manifest_sha256"]==files[name(RAW/"manifest.json")]
            and batch["status"]=="COMPLETE_FOUR_RAW_SPATIAL_BRIDGES_AWAITING_INDEPENDENT_PACKET_REPLAY"
            and batch["raw_components"]==672 and batch["root_asserted"] is False
            and batch["later_waves_launched"] is False and batch["independent_packet_replay_required"] is True
            and all(x["pid"]==ext["pid"] and pid(x["pid"]) and x["stage"]==STAGE
                    and x["jobs"]==[list(j) for j in JOBS] and x["bindings"]==worker_bindings for x in (start,end,batch)),
            "raw supervisor identity and real observed phase mismatched")
    require(tuple(row["case_id"] for row in batch["cases"])==FIRST,"all four same routed records")
    rows={}
    worker_pids=set()
    for e in batch["cases"]:
        case=e["case_id"]
        require(e["target"]==admission["targets"][case] and tuple(x["part"] for x in e["parts"])==tuple(COUNTS),
                "mixed first-wave target or parts")
        for index in e["parts"]:
            part=index["part"]
            folder=RAW/"cases"/case/part
            launched,started,term,obs,m=(seal(folder/(k+".json")) for k in sub)
            worker=index["pid"]
            require(pid(worker) and worker not in worker_pids,"distinct literal owned child PID required")
            worker_pids.add(worker)
            common={"pid":worker,"supervisor_pid":ext["pid"],"stage":STAGE,"case_id":case,"part":part,
                    "bindings":worker_bindings}
            require(all(all(x[k]==v for k,v in common.items()) and pid(x["pid"]) and pid(x["supervisor_pid"])
                        for x in (index,launched,started,term,obs,m))
                    and zero(index["observed_exit_code"]) and zero(obs["observed_exit_code"])
                    and started["status"]=="RUNNING" and term["status"]=="FINISHED" and zero(term["exit_code"])
                    and obs["terminal_sha256"]==index["terminal_sha256"]==files[name(folder/"terminal.json")]
                    and term["results_sha256"]==index["manifest_sha256"]==files[name(folder/"manifest.json")]
                    and index["observed_os_exit_sha256"]==files[name(folder/"observed_os_exit.json")]
                    and m["status"]=="COMPLETE_RAW_REMAINING32_SPATIAL_PART_AWAITS_INDEPENDENT_REPLAY"
                    and m["target"]==e["target"] and m["root_asserted"] is False and m["plan"] in admission["plans"]
                    and m["plan"]["case_id"]==case and (m["point_force_scope"],m["whole_H_scope"])==(POINT,BOX)
                    and m["radius_exact"]==str(RADIUS) and m["packet_independent_replay_required"] is True,
                    "foreign, partial, nonliteral or unsuccessful owned raw child proof")
            cat=m["components"]
            require(len(cat)==m["counts"]["components"]==COUNTS[part],"complete actual raw part expected")
            for i,k in enumerate(cat):
                c=folder/f"{part}_{i:03d}.json"
                require(k["kind"]==part and k["slot"]==i and type(k["slot"]) is int
                        and k["sha256"]==files[name(c)] and k["file"]==c.name,"noncanonical actual leaf")
                seal(c,False)
            rows[case,part]=(folder,m)
    require(len(rows)==8,"eight distinct observed parts required")
    return rows,ext


def load_math():
    spec=importlib.util.spec_from_file_location("_firstwave_independent_math",OLD)
    M=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(M)
    return M


def prepare():
    checked={}
    a=cold(checked)
    # Strict one-sided containment diagnostic, exact boundary intersection,
    # and truly disjoint counterparts. Neither the code nor the samples use
    # producer intervals to produce the independent enclosure.
    def run(A,B):
        return hm.hull_and_compare(A,B)
    checks={
        "identity_only_first_wave_protocol":list(a["targets"])==list(FIRST),
        "independent_nonnested_overlap_does_not_narrow_truth":run((0,2),(1,3))==((Q(0),Q(3)),False,False),
        "nested_producer_contains_own":run((1,2),(0,3))==((Q(0),Q(3)),True,False),
        "exact_endpoint_touch_is_intersection":run((0,1),(1,2))==((Q(0),Q(2)),False,False),
        "disjoint_is_separate_hold_even_when_hull_exists":run((0,1),(2,3))==((Q(0),Q(3)),False,True),
        "no_actual_component_payloads_opened":not any(name(RAW) in k for k in checked),
    }
    require(yes(checks),"cold independent hull scalar fixtures failed")
    math_file=Path(hm.__file__).resolve()
    closure={"source_sha256":sha(SRC.read_bytes()),"local_bound_and_hull_source_sha256":sha(math_file.read_bytes())}
    p=save(HERE/"protocol.json",{"status":"FROZEN_FIRST_WAVE_INDEPENDENT_COMPONENT_HULL_REPLAY",**closure,
                                  "scope":SCOPE,"cold_inputs_sha256":checked,"first_wave":list(FIRST),
                                  "raw_full_closure_files":715,"case_components":168,"sufficient_tail_policy":"own tails against caps, separate raw tail checks",
                                  "disjoint_policy":"HOLD_ALL_REDUCTION_FOR_CASE; diagnostics retained",
                                  "actual_fresh_moment_force_H_or_AD_jobs":0})
    f=save(HERE/"preflight.json",{"status":"PASS",**closure,"protocol_sha256":p,"checks":checks,
                                  "actual_fresh_moment_force_H_or_AD_jobs":0,"actual_packet_payloads_read":0})
    print(json.dumps({**closure,"protocol_sha256":p,"preflight_sha256":f}))


def run(args):
    checked={}
    p=read(HERE/"protocol.json",args.protocol_sha256,checked)
    f=read(HERE/"preflight.json",args.preflight_sha256,checked)
    g=read(HERE/"root_source_review.json",args.source_review_sha256,checked)
    source={"source_sha256":sha(SRC.read_bytes()),
            "local_bound_and_hull_source_sha256":sha(Path(hm.__file__).read_bytes())}
    phase={**source,"protocol_sha256":args.protocol_sha256,"preflight_sha256":args.preflight_sha256}
    require(all(p[k]==f[k]==g[k]==v for k,v in source.items()) and yes(g["checks"]) and yes(f["checks"])
            and f["status"]==g["status"]=="PASS" and g["protocol_sha256"]==f["protocol_sha256"]==args.protocol_sha256
            and g["preflight_sha256"]==args.preflight_sha256,"current independent mathematical source review required")
    pins={}
    a=cold(pins)
    require(pins==p["cold_inputs_sha256"],"cold proof inputs drifted")
    checked.update(pins)
    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_remaining32_spatial_worker as W
    require(Path(W.__file__).resolve()==WORKER,"mixed actual admitted worker source")
    wargs=SimpleNamespace(protocol_sha256=PINS[WP/"protocol.json"],preflight_sha256=PINS[WP/"preflight.json"],
                          review=GATE,review_sha256=PINS[GATE],admission=ADMISSION,admission_sha256=PINS[ADMISSION])
    auth,wphase=W.gate(wargs)
    require(auth==a,"runtime reviewed producer admission differs")
    raw,os_ob=full_actual(args.admission_sha256,wphase,a,checked)
    M=load_math()
    row,d,tail,plans=W.reviewed_context(a,args.case_id)
    _,scaled=W.S.kernels()
    require(M.E is W.S.E and M.RADIUS==RADIUS and M.E.core.ctx.prec==640 and tuple(plans)==W.S.OPS,
            "reviewed saved finite-element math alias lost exact scope or precision")
    M.ad.ctx.prec=640
    case=args.case_id
    folder=HERE/"actual_replay"/case
    folder.mkdir(parents=True,exist_ok=False)
    phase|={"pid":os.getpid(),"case_id":case,"root_source_review_sha256":args.source_review_sha256,
            "root_producer_admission_sha256":args.admission_sha256}
    end={"status":"FAILED","exit_code":1,**phase}
    save(folder/"started.json",{"status":"RUNNING",**phase})
    t0=monotonic()
    leaves=[]
    disagreements={"point_force":0,"whole_box_H":0,"point_same_Y_cache_correction":0}
    disjoint=dict.fromkeys(disagreements,0)
    try:
        with patch.object(W.B,"evaluate",side_effect=AssertionError("producer mathematical path prohibited")), \
             patch.object(M,"replay_component",side_effect=AssertionError("old containment-gated wrapper prohibited")):
            for kind,cnt in COUNTS.items():
                path,rawpart=raw[case,kind]
                for i in range(cnt):
                    s=plans[kind,i]
                    meta=rawpart["components"][i]
                    actual=read(path/meta["file"],meta["sha256"],checked)
                    W.validate_packet(actual,s,row["target"],i,kind,tail)
                    require(actual["trial_box_sha256"]==d["entry"]["trial"]["trial_box_sha256"]
                            ==row["target"]["trial_box_sha256"] and actual["selector_plan_sha256"]==row["plan_sha256"]
                            and all(actual["worker_identity"][k]==rawpart[k] for k in
                                    ("pid","supervisor_pid","case_id","part","stage","bindings")),
                            "same root-observed producer packet and reviewed selector")
                    M.E.core.A.cache_clear()
                    M.ad.A.cache_clear()
                    own=hm.local(M,d,kind,i,s,actual,tail,scaled.S)
                    output=folder/meta["file"]
                    h=save(output,{"phase":phase,"packet_sha256":meta["sha256"],
                                   "selector_leaf_content_sha256":W.B.digest(s),**own})
                    leaves.append({"kind":kind,"slot":i,"native_file":name(output),"native_sha256":h,
                                   "packet_file":name(path/meta["file"]),"packet_sha256":meta["sha256"],
                                   "status":own["status"],"global_indices":own["global_indices"]})
                    for k,v in own["signed_channels"].items():
                        disagreements[k]+=len(v["producer_containment_misses"])
                        disjoint[k]+=len(v["disjoint_indices"])
                    print(json.dumps({"case_id":case,"kind":kind,"slot":i,"status":own["status"]}),flush=True)
        hold=any(disjoint.values())
        require(len(leaves)==168,"all 168 independently evaluated components required")
        admission_doc={"status":"HOLD_TRUE_DISJOINT_COMPONENT_DIAGNOSTIC" if hold else
                              "PASS_INDEPENDENT_ACTUAL_REMAINING32_COMPONENT_HULLS",
                       **phase,"scope":SCOPE,"source_and_actual_sha256":checked,"case_components":leaves,
                       "producer_observed_pid":os_ob["pid"],"producer_observed_unified_exec_session":os_ob["unified_exec_session"],
                       "record_id":d["entry"]["record_id"],"trial_box_sha256":row["target"]["trial_box_sha256"],
                       "model_contract_sha256":row["target"]["model_contract_sha256"],
                       "exact_drive":row["target"]["exact_drive"],"original_archived_shift_exact":row["target"]["archived_center_shift_exact"],
                       "cache_packet_sha256":row["target"]["point_cache_sha256"],"R_binary_sha256":row["target"]["R_binary_sha256"],
                       "center_binary_sha256":row["target"]["center_binary_sha256"],
                       "radius_exact":str(RADIUS),"plan_sha256":row["plan_sha256"],
                       "plan_review_sha256":row["independent_plan_review_sha256"],
                       "producer_containment_miss_counts":disagreements,"disjoint_counts":disjoint,
                       "full_independent_hulls_computed":True,"admit_hull_only_global_reduction_after_separate_review":not hold,
                       "overlap_or_producer_containment_used_as_proof":False,"root_asserted":False,
                       "seconds":monotonic()-t0}
        h=save(folder/"receipt.json",admission_doc)
        end.update(status="FINISHED_WITH_HOLD" if hold else "FINISHED",exit_code=2 if hold else 0,receipt_sha256=h)
        print(json.dumps({"receipt_sha256":h,"disjoint_counts":disjoint,"containment_misses":disagreements}),flush=True)
        return end["exit_code"]
    except BaseException as error:
        save(folder/"failure.json",{"error":str(error),"traceback":traceback.format_exc(),"complete_leaves":leaves,**phase})
        raise
    finally:
        save(folder/"terminal.json",end)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    modes=cli.add_subparsers(dest="stage",required=True)
    modes.add_parser("prepare")
    act=modes.add_parser("run")
    act.add_argument("--case-id",required=True,choices=FIRST)
    for n in ("admission-sha256","source-review-sha256","protocol-sha256","preflight-sha256"):
        act.add_argument("--"+n,required=True)
    args=cli.parse_args()
    if args.stage=="prepare":
        prepare()
    else:
        sys.exit(run(args))

