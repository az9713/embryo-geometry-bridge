"""Reproduce the old shutdown hazard with actual dummy child, injected interrupts."""
import importlib.util
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from threading import Event
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
sys.path.insert(0,str(ROOT/"scripts"))
archive=HERE/"worker_development_revision_0/worker_frozen.py"
spec=importlib.util.spec_from_file_location("old_spatial_supervisor",archive)
old=importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
assert old.sha(archive)=="3009234a126b5b1b6997e4795889ee2a8ba9c3cc595882b247e21903fc7e27b6"

live=[]
retained=[]
child_started=Event()
calls=[]
base=old.ThreadPoolExecutor

class DoubleInterruptedExecutor(base):
    def shutdown(self,wait=True,*,cancel_futures=False):
        calls.append(("shutdown",wait))
        if wait and len(calls)==2:
            raise KeyboardInterrupt("second interrupt during blocking join")
        return super().shutdown(wait=wait,cancel_futures=cancel_futures)

def first_wait(*args,**kwargs):
    assert child_started.wait(10)
    calls.append(("wait",True))
    raise KeyboardInterrupt("first interrupt during normal wait")

def dummy(_):
    proc=subprocess.Popen([sys.executable,"-c","import time;time.sleep(60)"])
    live.append(proc)
    child_started.set()
    try:
        retained.append(proc.wait())
    finally:
        retained.append(None)

with TemporaryDirectory(dir=HERE) as tmp:
    probe=Path(tmp)
    try:
        with old.V.stage_lock(probe),patch.object(old,"ThreadPoolExecutor",DoubleInterruptedExecutor), \
             patch.object(old,"wait",side_effect=first_wait):
            old.bounded(["owned"],dummy,1)
    except KeyboardInterrupt as e:
        assert str(e)=="second interrupt during blocking join"
    try:
        assert len(live)==1 and live[0].poll() is None
        with old.V.stage_lock(probe):
            assert live[0].poll() is None
        defect=True
    finally:
        live[0].terminate()
        code=live[0].wait()
        assert type(code) is int
    assert calls==[("wait",True),("shutdown",True)]
print(old.save(HERE/"worker_development_revision_0/double_interrupt_reproduction.json",
    {"status":"REPRODUCED_OLD_LOCK_RELEASE_DURING_SECOND_INTERRUPT","source_sha256":old.sha(archive),
    "lock_reacquired_while_owned_child_alive":defect,"dummy_child_pid":live[0].pid,
    "genuinely_reaped_dummy_child_exit_code":code,"new_spatial_workers_or_mechanics":0}))
