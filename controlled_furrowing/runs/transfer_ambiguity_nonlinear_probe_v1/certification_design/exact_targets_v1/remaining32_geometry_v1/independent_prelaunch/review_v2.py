"""Cold independent review of the revised remaining-state stage launcher."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

from review_stage_inputs import HERE,ROOT,TREC,TRSHA,TSHA,STAGES,check,geometry_variation_inputs,sha


SOURCE=ROOT/'scripts/transfer_ambiguity_remaining32_geometry_variation.py'
HOLD='8859b9f22e4781f4533e641cbd3270eb5d7723020fd16ad43a86c0869cad9c5e'
BINDING={'source_sha256':'68bca869ee486d7dc25b34cae3ceee1ddc2adc53500a20301c0038c259a2b26a',
         'target_manifest_sha256':TSHA,'target_context_review_sha256':TRSHA,
         'geometry_protocol_sha256':'754da64bb96cd8a19707181315e6c3e3f014e8d37a898b83c1693342c82c9850',
         'geometry_preflight_sha256':'d70a7217c4dbbf68f789bb59193321252ea4192576c3f41790cc7daea38fb15a',
         'variation_protocol_sha256':'6d04fc82779430e66c2d27b5ab7aa508f522f7ffac03b87f8f2e8e2587e1d128',
         'variation_preflight_sha256':'abbbeb7c0c303db2d2675398f8ea09155a19fd969b121e7936d7f4a1a468e264'}
STATUS='PASS_INDEPENDENT_REMAINING32_GEOMETRY_VARIATION_PRELAUNCH'


def rejected(f):
    try:
        f()
    except ValueError:
        return True
    return False


def gate_roundtrip(m,orders,receipt):
    review_path=HERE/'review.json'; review_h='a'*64
    read=m.read
    fixture=[receipt]
    def proxy(path,h,*a,**k):
        if path==review_path:
            check(h==review_h,'mock reviewer digest drift')
            return fixture[0]
        return read(path,h,*a,**k)
    with patch.object(m,'read',side_effect=proxy):
        for stage in STAGES:
            args=SimpleNamespace(stage=stage,review=review_path,review_sha256=review_h,
                                 target_review=TREC,target_review_sha256=TRSHA)
            p,b=m.launch_gate(args)
            check(p['case_order']==orders[stage] and b=={**BINDING,'launcher_review_sha256':review_h},
                  'actual frozen launch_gate phase roundtrip')
        for altered in ({**receipt,'exact_spatial_bridge_authorized':True},
                        {**receipt,'checks':{'integer_truthy':1}},
                        {**receipt,'variation_preflight_sha256':'b'*64},
                        {**receipt,'permitted_case_ids':{**orders,'variation':orders['variation'][::-1]}}):
            fixture[0]=altered
            check(rejected(lambda:m.launch_gate(args)),'wrong stage set or external binding accepted')


def typed_os(m):
    ident={'pid':9781,'supervisor_pid':8451,'case_id':'cold','stage':'geometry'}
    good={**ident,'root_asserted':False}
    terminal={**ident,'status':'FINISHED','exit_code':0,'results_sha256':'9'*64}
    m.child_admit(0,good,terminal,ident,'9'*64)
    for k in ('pid','supervisor_pid'):
        for value in (float(ident[k]),True,0):
            check(rejected(lambda:m.child_admit(0,{**good,k:value},terminal,ident,'9'*64))
                  and rejected(lambda:m.child_admit(0,good,{**terminal,k:value},ident,'9'*64))
                  and rejected(lambda:m.child_admit(0,good,terminal,{**ident,k:value},'9'*64)),
                  'opaque JSON numeric equals typed PID')
    for code in (0.0,False,1):
        check(rejected(lambda:m.child_admit(code,good,terminal,ident,'9'*64))
              and rejected(lambda:m.child_admit(0,good,{**terminal,'exit_code':code},ident,'9'*64)),
              'opaque or failed child exit accepted')
    check(rejected(lambda:m.child_admit(0,good,terminal,ident,'8'*64)),'foreign result accepted')


class WaitInterrupt:
    """Wrap a real OS process and inject an interrupt into exactly its first normal wait."""
    def __init__(self,proc):
        self.proc=proc
        self.injected=False

    def __getattr__(self,key):
        return getattr(self.proc,key)

    def wait(self,timeout=None):
        if timeout is None and not self.injected:
            self.injected=True
            raise KeyboardInterrupt('synthetic supervisor wait interruption')
        return self.proc.wait(timeout=timeout)


def lifecycle(m):
    write=m.save
    real_popen=subprocess.Popen
    census=[]
    with TemporaryDirectory() as tmp:
        scratch=Path(tmp)
        for kind in ('launched_record','wait_interrupt','clean'):
            stage=scratch/kind; cases=stage/'cases'; cases.mkdir(parents=True)
            case=cases/'one'; case.mkdir()
            observed=[]; owned=[]
            holder={'stage':'geometry','launch_group':'9'*64,'case_id':'one','supervisor_pid':9081,'bindings':{}}
            def popen(script,**opts):
                check(script==['test-process-seam'],'would start non-fixture process')
                code='import sys; sys.exit(0)' if kind=='clean' else 'import time; time.sleep(20)'
                proc=real_popen([sys.executable,'-B','-c',code],**opts)
                owned.append(proc)
                return WaitInterrupt(proc) if kind=='wait_interrupt' else proc
            def save(path,value):
                if kind=='launched_record' and path.name=='launched.json':
                    raise FileExistsError('synthetic write error after OS process spawn')
                if path.name=='abnormal_child_reap.json':
                    check(len(owned)==1 and type(owned[0].poll()) is int
                          and rejected(lambda: reacquire(m,stage)),
                          'must observe reaped OS child while unchanged exclusive stage lock remains held')
                return write(path,value)
            try:
                with m.stage_lock(stage),patch.object(m.subprocess,'Popen',side_effect=popen),patch.object(m,'save',side_effect=save):
                    try:
                        ident,code=m.run_owned_child(['test-process-seam'],case,holder,observed)
                    except (KeyboardInterrupt,FileExistsError) as error:
                        check((kind=='wait_interrupt' and isinstance(error,KeyboardInterrupt))
                              or (kind=='launched_record' and isinstance(error,FileExistsError)),
                              'changed original supervisor failure')
                    else:
                        check(kind=='clean' and code==0 and ident['pid']==owned[0].pid,
                              'successful native process code differed')
                    check(len(owned)==1 and type(owned[0].poll()) is int,'owned subprocess exit unobserved')
                    code=owned[0].returncode
                    if kind=='clean':
                        check(not observed and not (case/'abnormal_child_reap.json').exists(), 'clean exit marked abnormal')
                        write(case/'observed_os_exit.json',{**holder,'pid':owned[0].pid,'observed_exit_code':code})
                    else:
                        note,=observed
                        check(note['observed_exit_code']==code and type(code) is int and note['pid']==owned[0].pid
                              and note['scientific_candidate_admitted'] is False
                              and note==json.loads((case/'abnormal_child_reap.json').read_bytes()),
                              'owner OS exit must remain retained even with no admitted result')
                with m.stage_lock(stage):
                    m.guard_previous_children(cases,'geometry')
                    (cases/'lost_observation').mkdir()
                    check(rejected(lambda:m.guard_previous_children(cases,'geometry')),
                          'unclosed previous case must block disjoint restart')
                census.append({'kind':kind,'test_child_pid':owned[0].pid,'literal_observed_child_exit_code':code})
            finally:
                for proc in owned:
                    if proc.poll() is None:
                        proc.kill()
                        proc.wait(timeout=5)
    return census


def reacquire(m,stage):
    with m.stage_lock(stage):
        return


def scalar_math(m):
    m.A.X.bootstrap()
    arb=m.core
    # Distinct signed axes for |R| row multiplication, with exact fractions.
    matrix=arb.arb_mat(3,3,[2,0,1,1,3,0,0,1,5])
    v=arb.arb_mat(3,1,[1,2,3])
    check(m.majorize(matrix,v)==[5,7,17],'positive R row majorant changed')
    bad=arb.arb_mat(3,1,[1,-2,3])
    check(rejected(lambda:m.majorize(matrix,bad)),'signed negative Hessian majorant accepted')


def main():
    check(sha(SOURCE)==BINDING['source_sha256'] and sha(HERE/'hold_4dd6.json')==HOLD,
          'v2 frozen source or retained independent HOLD changed')
    check(all(not (p/'actual').exists() for p in STAGES.values()),'stage has already run')
    matched=geometry_variation_inputs(BINDING)
    sys.path.insert(0,str(ROOT/'scripts'))
    spec=importlib.util.spec_from_file_location('independent_launcher_v2',SOURCE)
    launcher=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(launcher)
    typed_os(launcher)
    witnesses=lifecycle(launcher)
    scalar_math(launcher)
    checks={'all_original_32_cases_same_field_center_R_gain_times_dose_and_operator':True,
            'all_32_trial_650_coordinate_hulls_and_144_bulk_24_surface_free_maps':True,
            '22_native_34_matched_7_fresh_smaller_3_fresh_larger_selector_scope':True,
            'both_recentered_archived_field_xi_and_sigma_goals_remain':True,
            'point_and_whole_cube_separate_cached_derivative_enclosures':True,
            'actual_OS_PIDs_booleans_float_alias_and_mixed_digest_rejections':True,
            'real_dummy_OS_child_reaped_on_record_write_failure_and_wait_interrupt_while_lock_held':True,
            'unclosed_previous_child_reservation_blocks_disjoint_restart_without_replay':True,
            'positive_stage_gate_roundtrip_with_authentic_both_phase_authorities':True,
            'forbid_unapproved_spatial_phase_or_mechanics_result_promotion':True}
    result={'status':STATUS,**BINDING,'checks':checks,'permitted_case_ids':matched['orders'],
            'exact_spatial_bridge_authorized':False,'independent_review_source_sha256':sha(Path(__file__)),
            'independent_metadata_source_sha256':sha(HERE/'review_stage_inputs.py'),
            'predecessor_hold_sha256':HOLD,'exact_metadata_reconstruction':matched,
            'operating_system_lifecycle_fixtures':witnesses,
            'new_geometry_certificates':0,'new_AD_derivatives':0,'new_spatial_integrations':0,'new_solves':0,
            'scope':'Whole-box geometry and finite cached derivative stage admission only. Retain true outer and native process exits, then independently replay numerical candidates. No exact-spatial force, derivative, minimum, root or observer promotion.'}
    gate_roundtrip(launcher,matched['orders'],result)
    check(all(not (p/'actual').exists() for p in STAGES.values()),'cold review created launcher actual folder')
    out=HERE/'review.json'
    with out.open('x',encoding='utf8') as f:
        json.dump(result,f,indent=2)
        f.write('\n')
    print(json.dumps({'status':STATUS,'receipt_sha256':sha(out),'cold_actual_geometry_or_variation_work':0,
                      'geometry':len(matched['orders']['geometry']),'variation':len(matched['orders']['variation'])}))


if __name__=='__main__':
    main()
