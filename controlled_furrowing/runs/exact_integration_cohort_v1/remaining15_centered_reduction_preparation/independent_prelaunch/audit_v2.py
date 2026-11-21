"""Independent follow-up for centered reduction admission schema v2."""
import ast
from copy import deepcopy
import json
from pathlib import Path
import sys
from unittest.mock import patch

import audit as T

HERE,OUT,ROOT,SOURCE=T.HERE,T.OUT,T.ROOT,T.SOURCE
PINS={"source_sha256":"af64f05e165a119d14be761ea79f258d10ac371f6f97de653b85c1c5f79f18a2",
      "protocol_sha256":"ec5419359e9ecfe21ba73a01ae1a1b72fed8561942920eaea43e27921236999d",
      "preflight_sha256":"984954626ade96dbb6ad802f7e465d674aed9087555ee3b419970adc85c22b82"}
OLD="51ed30eb7f7e9016f5f77853716d67392ff25440c7d39b02b603af5c1de0d583"
ARCHIVE="d597da6cd49ef8194878d253601bf133fd21270cbb8a0e3dddece15f84ab3dff"
RECORD=T.RECORD


def functions(raw):
    return {x.name:ast.unparse(x) for x in ast.parse(raw).body if isinstance(x,ast.FunctionDef)}


def checked_delta(current, archive):
    a,b=functions(archive),functions(current)
    same=("_boot","_declared_record","_own_gate","_worker_gate","_complete_files","_terminal_inputs",
          "_packet_review","_packet_metadata","_bounds","_packet_force")
    T.need(all(a[k]==b[k] for k in same),"proven caller/domain and terminal admission helpers changed")
    T.need(a["_calculate"].replace("admission['expected_worker_links']","admission[ADMITTED_WORKER_LINKS]")==b["_calculate"],
           "proven predictor-only and full derivative arithmetic changed beyond schema key")
    T.need("ADMITTED_WORKER_LINKS: links" in b["admit"] and "admission_schema_version" in b["admit"]
           and "prep['admission_schema_version'] == ADMISSION_VERSION" in b["reduce"]
           and "links == prep[ADMITTED_WORKER_LINKS]" in b["reduce"]
           and b["reduce"].index("_packet_review(")<b["reduce"].index("_calculate("),"consistent versioned admission and gate")


def independent_roundtrip(M,rt,D):
    folder=HERE/"mocked_schema_v2"
    admission=T.rel(folder/f"{RECORD}_input_admission.json")
    review_file=T.rel(folder/"independent_actual"/f"{RECORD}_packet_review.json")
    pins={"manifest_sha256":"a"*64,"terminal_sha256":"b"*64}
    links={"record_id":RECORD,"factory_nodeball_sha256":"c"*64}
    paths=[{"element":i,"file":f"fixtures/components/element_{i:03d}.json","sha256":"d"*64} for i in range(144)]
    docs, events = {}, []

    def writer(path,obj):
        key=T.rel(path)
        T.need(key not in docs,"synthetic write-once output")
        docs[key]=deepcopy(obj)
        events.append(("write",key))
        return D._digest(obj)

    def reader(key,expected,bindings,parse=True):
        T.need("/components/" not in key,"no component read")
        T.need(D._digest(docs[key])==expected,"mock packet-review/admission digest")
        events.append(("read",key))
        bindings[key]=expected
        return deepcopy(docs[key]) if parse else None

    def mkdir(path,**kwargs):
        T.need(path==folder/RECORD and kwargs=={"exist_ok":False},"exclusive post-review new record directory")
        events.append(("mkdir",T.rel(path)))

    def compute(*args):
        a=args[-2]
        T.need(args[7]==RECORD and a["expected_worker_links"]==links and a["admission_schema_version"]==2
               and a["expected_packets"]==paths and review_file in args[-1],"current exact schema and independently reviewed gate")
        events.append(("calculate",RECORD))
        return {"status":"TEST_ONLY","record_id":RECORD}

    with (patch.object(M,"OUT",folder),patch.object(M,"_boot",side_effect=lambda:(*rt[:7],{},rt[-1])),
          patch.object(M,"_own_gate",return_value=PINS),patch.object(M,"_terminal_inputs",return_value=(paths,links)),
          patch.object(M,"_declared_record",return_value=({},{})),patch.object(M,"_calculate",side_effect=compute),
          patch.object(D,"_bound",side_effect=reader),patch.object(D,"_write",side_effect=writer),patch.object(Path,"mkdir",mkdir)):
        M.admit(RECORD,"e"*64,pins)
        a=docs[admission];a_sha=D._digest(a)
        approved={"status":"PASS_CENTERED_PACKET_REPLAY","checks":{"fixture":True},"record_id":RECORD,
                  "components":144,"admission_sha256":a_sha,"centered_manifest_sha256":pins["manifest_sha256"],
                  "centered_terminal_sha256":pins["terminal_sha256"],**PINS}
        docs[review_file]={**approved,"admission_sha256":"9"*64}
        T.need(T.rejects(lambda:M.reduce(RECORD,a_sha,D._digest(docs[review_file]))),"wrong packet replay was accepted")
        T.need(not any(x[0] in ("mkdir","calculate") for x in events),"review must succeed before output and scatter")
        docs[review_file]=approved
        bad={**a,"expected_worker_links":{**links,"record_id":"foreign"}}
        docs[admission]=bad
        n=len(events)
        T.need(T.rejects(lambda:M.reduce(RECORD,D._digest(bad),D._digest(approved)))
               and ("read",review_file) not in events[n:],"altered worker links rejected before review")
        stale={**a,"admission_schema_version":1}
        docs[admission]=stale
        T.need(T.rejects(lambda:M.reduce(RECORD,D._digest(stale),D._digest(approved))),"old schema claimed")
        docs[admission]=a
        M.reduce(RECORD,a_sha,D._digest(approved))
    result=docs[T.rel(folder/RECORD/"results.json")]
    finished=docs[T.rel(folder/RECORD/"terminal.json")]
    T.need(result["admission_sha256"]==a_sha and result["actual_packet_review_sha256"]==D._digest(approved)
           and finished["exit_code"]==0 and finished["results_sha256"]==D._digest(result)
           and events.count(("calculate",RECORD))==1,"successful synthetic exact schema read/write roundtrip")
    return {"independent_v2_admit_reduce_mock_roundtrip":True,
            "stale_schema_and_foreign_worker_links_rejected_before_reduction":True,
            "wrong_actual_review_pin_rejected_before_exclusive_directory_and_assembly":True}


def run():
    var=T.rel(OUT)
    source=T.read(SOURCE,PINS["source_sha256"],False)
    old=T.read(T.rel(T.ARCHIVE/"exact_centered_cohort_reduction_frozen.py"),T.PINS["source_sha256"],False)
    T.read(T.rel(T.ARCHIVE/"preservation_manifest.json"),ARCHIVE,False)
    T.read(T.rel(HERE/"blocked_e5e11f7.json"),OLD,False)
    T.read(f"{var}/exact_centered_cohort_reduction_frozen.py",PINS["source_sha256"],False)
    p,pre=[T.read(f"{var}/{n}.json",PINS[n+"_sha256"]) for n in ("protocol","preflight")]
    T.need(p["admission_schema_version"]==2 and p["admitted_worker_links_field"]=="expected_worker_links"
          and p["status"]=="FROZEN_CENTERED_COHORT_REDUCTION_PREPARATION_REVIEW_PENDING"
          and pre["status"]=="PASS_GENERIC_CENTERED_REDUCTION_PREPARATION_ONLY" and T.yes(pre["checks"])
          and p["source_sha256"]==pre["source_sha256"]==PINS["source_sha256"]
          and pre["protocol_sha256"]==PINS["protocol_sha256"]
          and pre["actual_component_payloads_opened_or_global_sums"] is p["actual_component_payloads_opened_or_global_sums"] is False,
          "new frozen prep, no real packet read")
    for n,h in {**p["input_hashes"],**p["failed_attempts"]}.items():
        T.read(n,h,False)
    sys.path.insert(0,str(ROOT/"scripts"))
    import exact_centered_cohort_reduction as M
    rt=M._boot()
    W,H,runtime,A,_,_,_,b,ready=rt
    D=runtime[1].D
    T.need(T.rel(M.__file__)==SOURCE and M.ADMISSION_VERSION==2 and M.ADMITTED_WORKER_LINKS=="expected_worker_links",
           "current imported v2")
    checked_delta(source,old)
    for r in W.RECORDS:
        M._worker_gate(D,T.CALLER_REVIEW,r,b)
        M._declared_record(D,ready,r,b)
    checks={**T.arith(D,A),**T.structural(M),**independent_roundtrip(M,rt,D)}
    checks.update(dict.fromkeys(("bound_new_source_and_preserved_failed_revision","only_current_C_D_prior_C_radius_H_and_surface_route",
                                "actual_worker_OS_and_manifest_required","no_real_centered_payloads_or_global_sums"),True))
    receipt={"status":"PASS",**PINS,"checks":checks,"source_review_sha256":T.sha(Path(__file__).read_bytes()),
             "preserved_failed_review_sha256":OLD,"corrected_admission_schema_version":2,"record_scope":list(W.RECORDS),
             "actual_component_payloads_opened_or_global_sums":False,
             "authorization_scope":"Prelaunch proof review for the current C/D reduction. Require independently observed successful OS/worker and full manifest admission, then a separately frozen and independently reviewed actual packet replay before any global composition.",
             "proof_scope":"Use shared signed RH_bulk and common five parameter offsets at predictors only. The derivative bound remains the full original retained H and b, including all 24 surfaces."}
    raw=(json.dumps(receipt,indent=2,allow_nan=False)+"\n").encode()
    sha=T.sha(raw)
    with (OUT/"prelaunch_independent_review.json").open("xb") as f:
        f.write(raw)
    T.need(M._own_gate(D,sha,{})==PINS,"real self-gate readback")
    print(json.dumps({"status":"PASS","gate_sha256":sha}))


if __name__=="__main__":
    run()

