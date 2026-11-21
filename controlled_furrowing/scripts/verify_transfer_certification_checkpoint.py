"""Read back newly reviewed cached components; not a mathematical proof replay."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'runs/transfer_ambiguity_nonlinear_probe_v1/certification_design'
DERIV = BASE/'exact_targets_v1/cached_derivative_extension_v1'
PINS = {
    DERIV/'independent_actual_preparation/bridge_admission.json':'269229312793ec101cc301086d634dd4de89c62a5050ccfc60f185cf761dbc7d',
    BASE/'cached_pilots_v1/observer_application_v1/independent_readback/receipt.json':'939212f0690f4ccdf27ede99a5b7e7be6eff4a080adfd2dda78c0352c2e9c2e5',
    BASE/'exact_targets_v1/bridge_reduction_v1/independent_review/final_binding.json':'d6e71bcc68e6293fda36e04c3f2b506a7c34c8516492beeda699ad97c433d4e2',
}


def main():
    checked = {}
    def load(path, expected):
        path = Path(path).resolve()
        path.relative_to(ROOT)
        actual = sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f'changed bound file: {path}')
        checked[path.relative_to(ROOT).as_posix()] = actual
        return json.loads(path.read_text(encoding='utf-8'))
    wrapper, observer, reduction = [load(p,h) for p,h in PINS.items()]
    native = load(ROOT/wrapper['native_independent_review_file'],wrapper['native_independent_review_sha256'])
    terminal = load((ROOT/wrapper['native_independent_review_file']).with_name('terminal.json'),
                    wrapper['native_independent_terminal_sha256'])
    os_exit = load((ROOT/wrapper['native_independent_review_file']).with_name('os_exit_observation.json'),
                   wrapper['native_independent_os_exit_sha256'])
    produced = load(DERIV/'actual/results.json',wrapper['derivative_results_sha256'])
    producer_end = load(DERIV/'actual/terminal.json',wrapper['derivative_terminal_sha256'])
    for path,h in native['producer_and_input_sha256'].items():
        # All files, including binary archives, are hashed before optional JSON reads.
        target = (ROOT/path).resolve()
        target.relative_to(ROOT)
        if sha256(target.read_bytes()).hexdigest() != h:
            raise ValueError(f'changed derivative dependency: {path}')
        checked[path] = h
    for row in native['reviewed_cases']:
        load(ROOT/row['native_case_file'],row['native_case_sha256'])
    checks = {
        'native_actual_derivative_status_and_literal_checks':native['status']=='PASS_INDEPENDENT_FOUR_LARGER_CACHED_DERIVATIVES'
            and all(type(v) is bool and v for v in native['checks'].values()),
        'native_terminal_and_true_os_exit':terminal['status']=='FINISHED' and terminal['exit_code']==0
            and type(os_exit['observed_exit_code']) is int and os_exit['observed_exit_code']==0
            and terminal['pid']==os_exit['pid']==native['source_and_phase']['pid']
            and terminal['receipt_sha256']==os_exit['receipt_sha256']==wrapper['native_independent_review_sha256'],
        'same_producer_results_terminal_four_child_exits':producer_end['status']=='FINISHED'
            and producer_end['results_sha256']==wrapper['derivative_results_sha256']
            and len(producer_end['workers'])==4
            and all(type(w['observed_exit_code']) is int and w['observed_exit_code']==0 for w in producer_end['workers'])
            and produced['status']=='COMPLETE_FOUR_CACHED_DERIVATIVE_ATTEMPTS_UNREVIEWED',
        'wrapper_rows_identical_to_native_proof_rows':wrapper['reviewed_cases']==native['reviewed_cases'],
        'four_positive_whole_radius_headrooms':len(native['reviewed_cases'])==4
            and Q(wrapper['radius_exact'])==Q(1,2**34)
            and all(Q(r['radius_exact'])==Q(1,2**34) and 0<=Q(r['reviewed_ZQ_upper_exact'])<1 for r in native['reviewed_cases']),
        'independent_cached_observer_application_scope':observer['status']=='PASS_INDEPENDENT_CACHED_ROOT_OBSERVER_APPLICATION_READBACK'
            and observer['exact_spatial_roots_certified']==0 and observer['full_goal_claimed'] is False,
    }
    if not all(type(v) is bool and v for v in checks.values()):
        raise ValueError(checks)
    out = BASE/'exact_targets_v1/certification_checkpoint_v1'
    out.mkdir(exist_ok=False)
    result = {'status':'PASS_BOUND_COMPONENT_READBACK_ONLY','checks':checks,'checked_file_sha256':dict(sorted(checked.items())),
              'verifier_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
              'reviewed_ZQ_upper_exact_by_pilot':{r['pilot_id']:r['reviewed_ZQ_upper_exact'] for r in native['reviewed_cases']},
              'cached_observer_comparison_counts':observer['comparison_counts'],
              'cached_observer_large_gap_um_exact':observer['exact_large_gap_um'],
              'pure_reduction_review_binding':reduction,
              'scope':'Hash and scope readback of existing independent mathematical reviews; no new AD/integral replay.',
              'new_exact_spatial_integral_roots_certified':0,'full_five_gate_goal_complete':False}
    path = out/'receipt.json'
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'receipt_sha256':sha256(path.read_bytes()).hexdigest(),'checked_files':len(checked),'checks':checks}))


if __name__ == '__main__':
    main()
