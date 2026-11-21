"""Inert synthetic test of the later selector's genuine-observation gate."""
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent))
import replay as R


def reject(fn,why=None):
    try:
        fn()
    except (ValueError,KeyError) as e:
        return why is None or why in str(e)
    return False


def fixture(M,case):
    with TemporaryDirectory(dir=R.HERE) as root:
        base=Path(root)
        p=base/case
        p.mkdir()
        bind={"fixture":True}
        meta={"case_id":case,"pid":27,"bindings":bind}
        leaves=[]
        for k,i in M.ORDER:
            file=p/f"{k}_{i:03d}.json"
            h=M.save(file,{"inert_not_scientific_data":True})
            leaves.append({"kind":k,"slot":i,"file":file.name,"sha256":h})
        st=M.save(p/"started.json",{"status":"RUNNING_SELECTOR_PLAN",**meta})
        ph=M.save(p/"plan.json",{"status":"COMPLETE_REMAINING32_CASE_SELECTOR_PLAN_AWAITS_INDEPENDENT_REVIEW",
              **meta,"components":leaves,"spatial_workers_authorized":False,"new_spatial_integrations":0,
              "whole_box_radius_exact":str(M.R),"point_force_scope":M.POINT,"whole_H_scope":M.BOX})
        th=M.save(p/"terminal.json",{"status":"FINISHED","exit_code":0,"plan_sha256":ph,**meta})
        native={"status":"PASS_OBSERVED_REMAINING32_SELECTOR_TRUE_OS_EXIT","case_id":case,"pid":27,
                    "unified_exec_session":8307,"observed_exit_code":0,"plan_sha256":ph,"terminal_sha256":th,"started_sha256":st}
        outer=p/"root_os_exit_observation.json"
        def run(body):
            if outer.is_file():
                outer.unlink()
            h=M.save(outer,body)
            links={}
            with patch.object(R.E.S,"PLAN",base):
                directory,obs,doc=R.admitted(M,case,h,links)
            return directory==p and obs==body and len(doc["components"])==168 and len(links)==172
        assert run(native)
        assert reject(lambda:run({**native,"observed_exit_code":False}),"literal observed")
        assert reject(lambda:run({**native,"pid":27.0}),"literal observed")
        assert reject(lambda:run({**native,"case_id":"alien"}),"literal observed")
        assert reject(lambda:run({k:v for k,v in native.items() if k!="started_sha256"}),"started_sha256")
        assert reject(lambda:run({**native,"started_sha256":"e"*64}),"changed checked bytes")
        path=p/leaves[-1]["file"]
        path.write_text('{"mutation": true}',encoding="utf-8")
        assert reject(lambda:run(native),"changed checked bytes")


def main():
    M,_,_,profiles,_=R.inputs()
    first=profiles[R.E.WAVES[0]]["case_order"][0]
    assert reject(lambda:R.route(profiles,R.E.WAVES[1],first))
    fixture(M,first)
    print(R.save(R.HERE/"synthetic_observation_readback.json",{"status":"PASS","source_sha256":R.sha(R.HERE/"replay.py"),
        "synthetic_source_sha256":R.sha(Path(__file__)),"checks":{"full_168_inert_leaves_and_exact_process_metadata":True,
        "missing_or_changed_start_hash_rejected":True,"wrong_case_foreign_wave_boolean_exit_and_float_pid_rejected":True,
        "corrupted_plan_leaf_bytes_rejected":True},"actual_new_selectors_or_mechanics":0,"actual_plan_leaf_reads":0}))


if __name__=="__main__":
    main()
