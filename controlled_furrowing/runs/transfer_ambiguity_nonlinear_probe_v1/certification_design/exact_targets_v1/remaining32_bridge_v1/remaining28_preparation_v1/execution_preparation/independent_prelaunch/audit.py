"""Independent cold audit of the late-wave scheduling and receipt adapter."""
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
STUDY = next(p for p in HERE.parents if (p/"PLAN.md").is_file())
sys.path.insert(0,str(HERE.parents[1]))
import wave_execution as E

SOURCE = "44c4fea65168a77926f93a08c3817edc0645c94ec7ef625fa34b4656beade8b7"
PROTOCOL = "5057fecf76205bd0eee84cde3aff730b0011299d87c36fdebaad1be53cf10c72"
PREFLIGHT = "aa03be5fd4d2341326a830bea1081e421b3739711d51b50d9398dbe6c764e214"


def reject(fn):
    try:
        fn()
    except (ValueError, KeyError, FileNotFoundError):
        return True
    return False


def phase_assertions(p,f):
    b={"source_sha256":SOURCE,"protocol_sha256":PROTOCOL,"preflight_sha256":E.sha(E.PREP/"preflight.json"),
       "profiles_sha256":E.PROFILE_SHA,"joined_bridge_input_admission_sha256":E.JOINED_SHA}
    cold={"status":E.SOURCE_REVIEW,"checks":{"mock":True},**b,"maximum_owned_scientific_processes":4,
          "permitted_waves":list(E.WAVES),"actual_release_required":True}
    E.source_gate(cold,p,f,b)
    assert all(reject(lambda: E.source_gate(x,p,f,b)) for x in (
        {**cold,"checks":{"mock":1}}, {**cold,"actual_release_required":False}, {**cold,"source_sha256":"0"*64},
        {**cold,"permitted_waves":list(E.WAVES[1:])},{**cold,"maximum_owned_scientific_processes":5}))
    rel={**b,"status":E.RELEASE,"checks":{"release":True},"case_order":list(E.W.FIRST),"following_waves":list(E.WAVES),
         "first_wave_independent_receipt_sha256":"a"*64,"selector_invocations_only_under_this_release":True,
         "per_wave_plan_review_and_raw_admission_required":True}
    obs={"status":E.FIRST_OS,"pid":2098,"unified_exec_session":1796,"observed_exit_code":0,"receipt_sha256":"a"*64}
    term={"status":"FINISHED","pid":2098,"exit_code":0,"receipt_sha256":"a"*64}
    native={"status":"COMPLETE_INDEPENDENT_FIRST_WAVE_SIGNED_ENDPOINT_REPLAY","pid":2098,
            "case_results":[{"case_id":case} for case in E.W.FIRST],"full_five_gate_complete":False}
    E.release_summary(rel,obs,term,native,b)
    assert all(reject(lambda:E.release_summary(rel,x,term,native,b)) for x in (
        {**obs,"unified_exec_session":True},{**obs,"observed_exit_code":False},
        {**obs,"pid":2098.0},{**obs,"receipt_sha256":"0"*64}))
    assert all(reject(lambda:E.release_summary(rel,obs,term,x,b)) for x in (
        {**native,"case_results":native["case_results"][::-1]},{**native,"case_results":native["case_results"][:2]},
        {**native,"full_five_gate_complete":True}))


def native_lifecycle(tmp):
    V=E.V
    worker_dir=tmp/"inert_child"
    worker_dir.mkdir()
    handles=[]
    real_popen=subprocess.Popen
    def popen(*args,**kwargs):
        p=real_popen(*args,**kwargs)
        handles.append(p)
        return p
    def fail_after_popen(path,obj):
        if path.name=="launched.json":
            assert handles and handles[-1].poll() is None
            assert reject(lambda: acquired(tmp/"shared")), "competing wave cannot reacquire active lock"
            raise OSError("inject a failure after native process creation, before launch proof")
        return E.save(path,obj)
    def inert():
        return V.run_owned_child([sys.executable,"-B","-c","import time;time.sleep(30)"],worker_dir,
            {"supervisor_pid":1,"wave_id":"inert","case_id":"only","part":"inert","bindings":{}},[])
    with patch.object(V.subprocess,"Popen",side_effect=popen), patch.object(V,"save",side_effect=fail_after_popen):
        with V.stage_lock(tmp/"shared"):
            try:
                E.W.bounded(["one"],lambda _:inert(),1)
            except OSError as er:
                assert "inject" in str(er)
            else:
                raise RuntimeError("injected native launch receipt failure did not propagate")
            assert len(handles)==1 and type(handles[0].poll()) is int
    row=E.read(worker_dir/"abnormal_child_reap.json",E.sha(worker_dir/"abnormal_child_reap.json"))
    assert row["status"]=="REAPED_OWNED_CHILD_AFTER_SUPERVISOR_ERROR" and row["pid"]==handles[0].pid
    assert row["observed_exit_code"]==handles[0].returncode and row["scientific_candidate_admitted"] is False
    acquired(tmp/"shared")
    return {"observed_inert_child_exit":row["observed_exit_code"],
            "owned_child_reaped_before_competing_lock_released":True}


def acquired(folder):
    with E.V.stage_lock(folder):
        return True


def completed_dummy_raw(tmp,p):
    wave=p["wave_id"]
    bind={"inert_only":True}
    args=SimpleNamespace(wave=wave,wave_admission_sha256="a"*64)
    root=tmp/wave/"raw"
    root.mkdir(parents=True)
    a={"plans":[{"case_id":case,"target":{"case_id":case}} for case in p["case_order"]]}
    def part_mark(case,part,pid):
        return {"pid":pid,"supervisor_pid":71,"wave_id":wave,"case_id":case,"part":part,"bindings":bind}
    with patch.object(E,"RUNS",tmp),patch.object(E,"raw_gate",return_value=(p,a,bind)):
        cases=[]
        links={}
        def put(path,data):
            sha=E.save(path,data)
            links[E.rel(path)]=sha
            return sha
        for z,parent in enumerate(a["plans"]):
            case=parent["case_id"]
            parts=[]
            for j,(kind,n) in enumerate(E.COUNTS.items()):
                pid=120+2*z+j
                ident=part_mark(case,kind,pid)
                folder=E.part_folder(wave,case,kind)
                folder.mkdir(parents=True)
                components=[]
                for i in range(n):
                    name=f"{kind}_{i:03d}.json"
                    path=folder/name
                    # No scientific packet bytes exist in this controlled fixture.
                    links[E.rel(path)]="c"*64
                    components.append({"kind":kind,"slot":i,"file":name,"sha256":"c"*64,"force_entries":1,"upper_H_entries":1})
                m=put(folder/"manifest.json",{"status":E.RAW_STATUS,**ident,"target":parent["target"],"plan":parent,
                      "packet_independent_replay_required":True,"point_force_scope":E.S.POINT,"whole_H_scope":E.S.BOX,
                      "radius_exact":str(E.S.RADIUS),"components":components,
                      "counts":{"components":n,"force_entries":n,"upper_H_entries":n},"root_asserted":False})
                t=put(folder/"terminal.json",{"status":"FINISHED",**ident,"exit_code":0,"results_sha256":m})
                obs={**ident,"observed_exit_code":0,"terminal_sha256":t}
                exit_sha=put(folder/"observed_os_exit.json",obs)
                put(folder/"launched.json",ident)
                put(folder/"started.json",{**ident,"status":"RUNNING"})
                parts.append({**obs,"manifest_sha256":m,"observed_os_exit_sha256":exit_sha})
            cases.append({"case_id":case,"target":parent["target"],"parts":parts})
        parent_ident={"pid":71,"wave_id":wave,"bindings":bind,"jobs":p["part_jobs"],"case_order":p["case_order"]}
        m=put(root/"manifest.json",{"status":"COMPLETE_RAW_LATER_WAVE_SPATIAL_BRIDGES_AWAITING_INDEPENDENT_REPLAY",
            **parent_ident,"cases":cases,"raw_components":672,"root_asserted":False,"independent_packet_replay_required":True})
        put(root/"terminal.json",{"status":"FINISHED",**parent_ident,"exit_code":0,"manifest_sha256":m})
        put(root/"started.json",{"status":"RUNNING",**parent_ident})
        outer={"status":"PASS_ROOT_OBSERVED_LATER_WAVE_RAW_SPATIAL_TERMINALS","wave_id":wave,"pid":71,"unified_exec_session":222,
               "observed_exit_code":0,"raw_wave_admission_sha256":"a"*64,"sha256":links}
        def run(name,test):
            scratch=root/name
            h=E.save(scratch,test)
            with patch.object(E,"read",side_effect=lambda p,s:E.S.read(scratch if p.name=="root_terminal_admission.json" else p,s)):
                return E.open_raw_phase(args,h)
        p0,a0,refs,x=run("dummy_valid.json",outer)
        assert p0==p and a0==a and len(refs)==8 and len(x["sha256"])==715
        fewer=dict(links)
        fewer.pop(next(k for k in fewer if k.endswith("bulk_000.json")))
        assert reject(lambda:run("dummy_underfilled.json",{**outer,"sha256":fewer}))
        assert reject(lambda:run("dummy_alias.json",{**outer,"observed_exit_code":False}))
        assert reject(lambda:run("dummy_scope.json",{**outer,"wave_id":"wave_08" if wave!="wave_08" else "wave_07"}))
        # A hard-stop reservation lacks a terminal. Normal or failed typed closes are accepted.
        prev=E.raw_folder(E.WAVES[-1])
        prev.mkdir(parents=True)
        E.save(prev/"started.json",{"pid":877,"wave_id":E.WAVES[-1],"status":"RUNNING"})
        assert reject(E.check_closed_prior_waves)
        E.save(prev/"terminal.json",{"pid":877,"wave_id":E.WAVES[-1],"status":"FAILED","exit_code":1})
        E.check_closed_prior_waves()
        return True


def main():
    p=E.read(E.PREP/"protocol.json",PROTOCOL)
    ph=PREFLIGHT
    f=E.read(E.PREP/"preflight.json",ph)
    assert E.sha(E.SOURCE)==E.sha(E.PREP/"execution_adapter_frozen.py")==SOURCE and E.V.yes(f["checks"])
    scopes,bind=E.cold()
    assert p["bridge_binding"]==bind and p["wave_ids"]==list(scopes)==list(E.WAVES)
    flat=[k for wave in scopes.values() for k in wave["case_order"]]
    frozen=E.P.get(E.P.BRIDGE/"preparation/manifest.json")["cases"]
    assert len(flat)==len(set(flat))==28 and set(flat)=={z["case_id"] for z in frozen}-set(E.W.FIRST)
    assert all(p["per_wave_jobs"][k]==v["part_jobs"] for k,v in scopes.items())
    phase_assertions(p,f)
    assert E.V.yes(E.fixtures(scopes))
    with TemporaryDirectory(dir=HERE,prefix="scratch_") as tmp:
        l=native_lifecycle(Path(tmp))
        assert completed_dummy_raw(Path(tmp)/"mock_waves",scopes[E.WAVES[0]])
    # Recheck the same already-observed first-wave plan replay metadata only.
    first=E.read(E.W.WAVE_GATE,"8468da6cd555797298df182d6213a93a4d054029bdfca9b332aea80ca39f6855")
    assert all(E.independent_plan_exit(z,z["case_id"]) is None for z in first["plans"])
    assert reject(lambda:E.independent_plan_exit(first["plans"][0],first["plans"][1]["case_id"]))
    print(E.save(HERE/"final_fixtures.json",{"status":"PASS","source_sha256":SOURCE,"protocol_sha256":PROTOCOL,"preflight_sha256":ph,
          "fixture_source_sha256":E.sha(Path(__file__)),"checks":{"all_28_targets_exact_seven_profiles":True,
          "root_release_type_scope_and_sibling_plan_spoofs_rejected":True,"same_715_typed_references_and_eight_os_children":True,
          "raw_reader_partial_or_wrong_wave_refused_and_stranded_parent_blocks_all_waves":True,
          "four_real_inert_owners_in_hardened_primitive_fixture_replayed":True,**{
          "owned_child_reaped_before_competing_lock_released":l["owned_child_reaped_before_competing_lock_released"]}},
          "inert_child_reap_exit_code":l["observed_inert_child_exit"],"actual_scientific_packets_opened":0,
          "actual_selector_spatial_hull_or_matrix_jobs":0}))


if __name__=="__main__":
    main()
