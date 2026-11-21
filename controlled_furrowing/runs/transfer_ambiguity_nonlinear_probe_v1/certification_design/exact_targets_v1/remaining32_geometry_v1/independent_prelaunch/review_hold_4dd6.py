"""Reproduce two launch-admission gaps against the immutable 4dd6 candidate."""
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch


HERE=Path(__file__).resolve().parent
EXACT=HERE.parents[1]
ROOT=EXACT.parents[3]
SOURCE=HERE/'reviewed_launcher_4dd6.py'
EXPECTED='4dd6256d2eec12158f022f80c5b4c38aedc7d540ff8999a9c5f913e9751290c4'


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


class FailedLaunchRecord(Exception):
    """An enrolled subprocess has started but its launch record could not be written."""


class SyntheticLiveChild:
    pid=9457

    def __init__(self):
        self.events=[]
        self.alive=True

    def wait(self,*args,**kwargs):
        self.events.append('wait')
        self.alive=False
        return 0

    def poll(self):
        self.events.append('poll')
        return None if self.alive else 0

    def terminate(self):
        self.events.append('terminate')
        self.alive=False

    kill=terminate


def witness(launcher):
    child=SyntheticLiveChild()
    write=launcher.save
    def failing_save(path,value):
        if path.name=='launched.json':
            raise FailedLaunchRecord('synthetic error after actual Popen seam')
        return write(path,value)
    with TemporaryDirectory() as scratch:
        root=Path(scratch); actual=root/'actual'
        args=SimpleNamespace(stage='geometry',case_id=['one'],review=HERE/'fixture.json',review_sha256='a'*64,
                             target_review=HERE/'fixture.json',target_review_sha256='b'*64)
        with patch.object(launcher,'OUTS',{'geometry':root}), \
             patch.object(launcher.subprocess,'Popen',return_value=child) as spawned, \
             patch.object(launcher,'save',side_effect=failing_save):
            try:
                with launcher.stage_lock(actual):
                    launcher.run_locked(args,{'case_order':['one','two']},{})
            except FailedLaunchRecord:
                pass
            else:
                raise AssertionError('fixture did not reach failing launch-record write')
            assert spawned.call_count==1
            terminal=next((actual/'launches').glob('*/terminal.json'))
            end=json.loads(terminal.read_bytes())
            assert child.alive and child.events==[]
            assert end['status']=='FAILED' and end['exit_code']==1 and end['workers']==[]
            assert not (actual/'cases/one/observed_os_exit.json').exists()
            with launcher.stage_lock(actual):
                assert child.alive and child.events==[]
                second_case_unreserved=not (actual/'cases/two').exists()
            assert second_case_unreserved
    literal={'pid':7,'supervisor_pid':9,'stage':'geometry'}
    end={**literal,'status':'FINISHED','exit_code':0,'results_sha256':'d'*64}
    malformed={**literal,'pid':7.0,'supervisor_pid':9.0,'root_asserted':False}
    launcher.child_admit(0,malformed,end,literal,'d'*64)
    return {'launch_record_failure_releases_stage_lock_with_unreaped_child':True,
            'subsequent_disjoint_case_remains_unreserved':True,
            'supervisor_writes_failed_native_terminal_but_no_child_os_record':True,
            'foreign_JSON_floating_child_result_pid_and_supervisor_pid_compare_equal':True,
            'mock_child_exit_cleanup_events':child.events}


def main():
    assert sha(SOURCE)==EXPECTED
    sys.path.insert(0,str(ROOT/'scripts'))
    spec=importlib.util.spec_from_file_location('launcher_4dd6_hold_snapshot',SOURCE)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    actual=witness(module)
    receipt={'status':'HOLD_INDEPENDENT_REMAINING32_GEOMETRY_VARIATION_PRELAUNCH',
             'source_sha256':EXPECTED,'review_source_sha256':sha(Path(__file__)),
             'geometry_protocol_sha256':'cb0b2f1054110d3b08b0cbc54399d66bfadf7231118afb0cbbcc249cc8ae31db',
             'geometry_preflight_sha256':'a26e7ec9bf59d894ed21f2ed57a4b119e480cf918aa30caa1fee4e65f53d380b',
             'variation_protocol_sha256':'d570e6d4dfc90a5dfd306900709d8863e6105866c350da526c389fbaf9ddf57a',
             'variation_preflight_sha256':'fd161e448464fa233f4078e3b5f335e215522e372b2b5af28ce714d6ce3d5f7f',
             'target_context_review_sha256':'67d33bbc3171839b1c578b02b5cfbf9ea6ba5d6cb351d0655dd14aaf99df75f9',
             'observed_witnesses':actual,'permitted_case_ids':{'geometry':[],'variation':[]},
             'stage_launch_authorized':False,'exact_spatial_bridge_authorized':False,
             'real_child_processes_created':0,'new_geometry_certificates':0,'new_AD_derivatives':0,
             'new_spatial_integrations':0,'new_solves':0,
             'scope':'Failure injection at the existing subprocess/write boundary and synthetic JSON type-alias test. No mechanics called.'}
    path=HERE/'hold_4dd6.json'
    with path.open('x',encoding='utf8') as handle:
        json.dump(receipt,handle,indent=2)
        handle.write('\n')
    print(json.dumps({'status':receipt['status'],'SHA256':sha(path)}))


if __name__=='__main__':
    main()
