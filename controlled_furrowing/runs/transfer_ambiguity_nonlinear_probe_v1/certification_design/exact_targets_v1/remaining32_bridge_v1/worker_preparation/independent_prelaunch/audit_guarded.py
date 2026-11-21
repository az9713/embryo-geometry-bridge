"""Independent review of the first-wave worker scheduler and scientific admission."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import inspect
import json
from pathlib import Path
import sys
from threading import Event, Lock
from time import monotonic, sleep
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
P=HERE.parent
BR=P.parent
ROOT=next(p for p in BR.parents if (p/"PLAN.md").is_file())
SRC=ROOT/"scripts/transfer_ambiguity_remaining32_spatial_worker.py"
ARCHIVED_HOLDS={
    "hold_300923.json":"1b2b81878bfa1af5663461ac4a51be0f20c02976151a21357379d605b2658993",
    "hold_301fc5.json":"c8c9237b02a1dba3ed2226e2cbb7eec385ab5fcddc691a13ee084ffa4eeddc5e",
}
BINDING={"bridge_manifest_sha256":"d42b7260283ba02ccc34b9e9ba69ce76d58507115c22615c90d28fa38020cb41",
         "bridge_cold_review_sha256":"89f9b5cacf89b4055fc404c00d3e3c42cf10fed20337fe455e996c2944ef73b2"}


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def exact_hashtable(path,h):
    assert sha(path)==h, path
    return json.loads(path.read_bytes())


def failed_second_thread_adjust(W):
    """Actual existing thread picks a queued wrapper whose submit never returned."""
    allow_first=Event()
    completed_first=Event()
    illicit_task=Event()
    entered_unreturned_wrapper=Event()
    held=Lock()
    pools=[]

    class InjectWorkItemNotReturned(ThreadPoolExecutor):
        def _adjust_thread_count(self):
            if not self._threads:
                pools.append(self)
                return super()._adjust_thread_count()
            allow_first.set()
            assert completed_first.wait(5)
            deadline=monotonic()+5
            while not self._work_queue.empty() and monotonic()<deadline:
                sleep(.001)
            assert self._work_queue.empty()
            entered_unreturned_wrapper.set()
            raise RuntimeError("independent injected thread adjustment failure")

    def task(name):
        if name=="armed":
            assert allow_first.wait(5)
            completed_first.set()
        else:
            illicit_task.set()
        return name

    try:
        with held,patch.object(W,"ThreadPoolExecutor",InjectWorkItemNotReturned):
            try:
                W.bounded(("armed","unreturned"),task,2)
            except RuntimeError as error:
                assert str(error)=="independent injected thread adjustment failure"
            else:
                raise AssertionError("thread adjustment fault swallowed")
        assert completed_first.is_set() and entered_unreturned_wrapper.is_set()
        assert held.acquire(blocking=False)
        held.release()
    finally:
        allow_first.set()
        for p in pools:
            p.shutdown(wait=True,cancel_futures=True)
    return not illicit_task.is_set()


def fake_plan_shape(W):
    cases, binding=W.bridge_manifest()
    expect=set(W.FIRST)
    assert expect=={"max_contrast_change_small_01","max_contrast_change_small_13",
                    "max_negative_contrast_change_small_01","max_negative_contrast_change_small_16"}
    assert binding["source_sha256"]==W.BRIDGE_SOURCE and binding["cold_review_sha256"]==BINDING["bridge_cold_review_sha256"]
    args={"source_sha256":"b"*64,"protocol_sha256":"c"*64,"preflight_sha256":"d"*64,**BINDING}
    target={k:W.case_payload(cases[k]) for k in W.FIRST}
    r={"status":W.REVIEW,**args,"checks":{"cold_only":True},"permitted_first_wave_case_ids":list(W.FIRST),
       "maximum_owned_workers":4,"no_later_waves_authorized":True}
    a={"status":W.ADMISSION,**args,"checks":{"actual_observations":True},
       "plans":[{"case_id":k} for k in W.FIRST],"targets":target,"no_later_waves_authorized":True,
       "independent_packet_replay_required":True}
    W.worker_gate_data(r,a,args,target)
    def rejects(b):
        try:
            W.worker_gate_data(r,b,args,target)
        except ValueError:
            return True
        return False
    return all(rejects(b) for b in ({**a,"plans":list(reversed(a["plans"]))},
                                    {**a,"independent_packet_replay_required":False},
                                    {**a,"targets":{**target,W.FIRST[1]:target[W.FIRST[0]]}},
                                    {**a,"bridge_cold_review_sha256":"f"*64}))


def main(src_hash,pro_hash,pf_hash):
    assert sha(SRC)==sha(P/"worker_frozen.py")==src_hash
    p=exact_hashtable(P/"protocol.json",pro_hash)
    f=exact_hashtable(P/"preflight.json",pf_hash)
    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_remaining32_spatial_worker as W
    assert Path(W.__file__).resolve()==SRC
    assert src_hash==p["source_sha256"]==f["source_sha256"] and f["protocol_sha256"]==pro_hash
    assert p["status"]=="FROZEN_RAW_EXACT_SPATIAL_FIRST_WAVE_WORKER" and f["status"]=="PASS" and W.V.yes(f["checks"])
    assert p["first_wave"]==list(W.FIRST) and p["jobs"]==[list(x) for x in W.JOBS]
    assert p["later_waves_authorized"] is p["root_asserted"] is False
    assert type(p["max_owned_workers"]) is int and p["max_owned_workers"]==4
    assert not W.ACTUAL.exists()
    held_paths={}
    for name,h in ARCHIVED_HOLDS.items():
        matching=[x for x in BR.glob("worker_development_revision_*/independent_prelaunch/"+name) if sha(x)==h]
        item,=matching
        held_paths[W.rel(item)]=h
    submit_source=inspect.getsource(ThreadPoolExecutor.submit)
    assert submit_source.index("self._work_queue.put(w)")<submit_source.index("self._adjust_thread_count()")<submit_source.index("return f")
    tests={"prior_two_cold_HOLDs_retained_unchanged":True,
           "new_repeated_console_and_drain_interrupts_are_retained":W.fixture_repeated_interrupts(),
           "unreturned_queued_future_cannot_start_under_injected_native_thread_adjust_failure":failed_second_thread_adjust(W),
           "pool_never_exceeds_four_and_returns_canonical_order":W.fixture_parallel(),
           "typed_raw_signed_point_and_closed_cube_hessian_and_same_Y_correction":W.fixture_raw(),
           "independent_target_and_plan_admission_rejects_wrong_case_phase_and_later_waves":fake_plan_shape(W),
           "source_matches_current_frozen_worker_gates_and_no_mechanics_run":True}
    assert W.V.yes(tests),tests
    output={"status":W.REVIEW,"source_sha256":src_hash,"protocol_sha256":pro_hash,"preflight_sha256":pf_hash,**BINDING,
            "checks":tests,"permitted_first_wave_case_ids":list(W.FIRST),"maximum_owned_workers":4,
            "no_later_waves_authorized":True,"prior_hold_hashes":held_paths,"cold_review_source_sha256":sha(Path(__file__)),
            "scope":"Four cases, eight possible parts, raw first wave only. Repeated console interrupts and failed unreturned thread submissions tested. Forced process termination does not confer a clean proof. Separate external actual admission required.",
            "actual_spatial_workers_started":0}
    dest=HERE/"review.json"
    with dest.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(output,f,indent=2)
        f.write("\n")
    print(sha(dest))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for x in ("source-sha256","protocol-sha256","preflight-sha256"):
        p.add_argument("--"+x,required=True)
    a=p.parse_args()
    main(a.source_sha256,a.protocol_sha256,a.preflight_sha256)

