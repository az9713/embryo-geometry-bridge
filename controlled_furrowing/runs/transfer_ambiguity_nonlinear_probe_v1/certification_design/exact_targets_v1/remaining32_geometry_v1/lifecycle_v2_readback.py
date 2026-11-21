"""Real synthetic process check with write failure. No geometry or AD entrypoint."""
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
sys.path.insert(0,str(ROOT/"scripts"))
import transfer_ambiguity_remaining32_geometry_variation as m

assert m.sha(m.SOURCE)=="68bca869ee486d7dc25b34cae3ceee1ddc2adc53500a20301c0038c259a2b26a"
m.bindings()
real_save=m.save
real_popen=m.subprocess.Popen
owned=[]

def fixture_popen(*args,**kwargs):
    proc=real_popen(*args,**kwargs)
    owned.append(proc)
    return proc

def fail_launch(path,value):
    if path.name=="launched.json":
        raise OSError("cold synthetic launch-start record write failure")
    return real_save(path,value)

with TemporaryDirectory(dir=HERE) as tmp:
    phase=Path(tmp)
    cases=phase/"cases"
    child=cases/"fixture"
    child.mkdir(parents=True)
    ident={"case_id":"fixture","stage":"geometry","supervisor_pid":m.os.getpid(),"launch_group":"a"*64,"bindings":{}}
    observed=[]
    with m.stage_lock(phase),patch.object(m,"save",side_effect=fail_launch), \
         patch.object(m.subprocess,"Popen",side_effect=fixture_popen):
        try:
            m.run_owned_child([sys.executable,"-B","-c","import time; time.sleep(60)"],child,ident,observed)
        except OSError as e:
            assert str(e)=="cold synthetic launch-start record write failure"
        else:
            raise AssertionError("expected cold injected write error")
        assert len(owned)==len(observed)==1
        proc=owned[0]
        assert type(proc.poll()) is int and proc.returncode==observed[0]["observed_exit_code"]
        assert m.positive_int(proc.pid) and proc.pid==observed[0]["pid"]
        assert proc.returncode!=0 and observed[0]["scientific_candidate_admitted"] is False
        m.guard_previous_children(cases,"geometry")
    with m.stage_lock(phase):
        pass

result={"status":"PASS_REAL_SYNTHETIC_OWNED_CHILD_REAP","source_sha256":m.sha(m.SOURCE),
        "fixture_sha256":sha256(Path(__file__).read_bytes()).hexdigest(),
        "dummy_child_observed_exit_code":proc.returncode,"dummy_child_pid":proc.pid,
        "real_reap_confirmed_before_stage_unlock":True,"actual_geometry_or_AD_or_solve_jobs":0}
print(m.save(HERE/"lifecycle_v2_readback.json",result))
