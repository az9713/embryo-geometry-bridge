"""Guarded entry point for reusing the retained material-mean diagnostic.

Default: read-only provenance check plus simulated stale-input probes.
--replay: check first, retain previous output, run existing mean code, check
again and retain the new output. Frozen mechanical modules are not edited.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import material_tangent_inference as m


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def manifest():
    protocol=json.loads((m.OUT/'mechanics_protocol.json').read_text())
    expected={m.ROOT/'scripts/material_tangent_inference.py':protocol['script_sha256']}
    expected.update({m.ROOT/'scripts'/name:value for name,value in protocol['physical_source_hashes'].items()})
    expected[m.ROOT/'scripts/quadratic_mean_derivatives.py']='681b69335a80ae0468dd98856c5562b24524827809ece465747b95454d152bff'
    expected[m.ROOT/'scripts/quadratic_endpoint_observer.py']='075dac2148ed26a3bfe83fadf1b97602e4b08536f2be38cb95c2054776c0b926'
    expected[m.ROOT/'scripts/annulus_observer.py']='77bcd1bc10e557de4a4be9c304041f9f0dfe0fcbac47871cabfbb11ffb7346b9'
    for key,path in m.CASES.items():
        row=json.loads((m.OUT/f'strength_{key}'/'verification.json').read_text())
        if row['input_sha256']!=protocol['input_hashes'][key]:
            raise RuntimeError('mechanical row does not bind the protocol input')
        expected[path]=protocol['input_hashes'][key]
        expected[m.OUT/f'strength_{key}'/'tangents.npz']=row['matrix_archive_sha256']
    return {str(path.relative_to(m.ROOT)):value for path,value in expected.items()}


def mismatches(expected,actual):
    return sorted(key for key in set(expected)|set(actual) if expected.get(key)!=actual.get(key))


def check(replay=False):
    expected=manifest();actual={name:digest(m.ROOT/name) for name in expected}
    bad=mismatches(expected,actual)
    if bad:raise RuntimeError('refusing stale diagnostic inputs: '+', '.join(bad))
    # Inject mismatches only into copied in-memory manifests; never alter the
    # retained source, coefficients or state/NPZ files to test the guard.
    probes={}
    for name in expected:
        altered=dict(actual);altered[name]='0'*64
        probes[name]=mismatches(expected,altered)==[name]
    checks={'all_actual_source_input_archive_hashes_match':not bad,
            'every_dependency_stale_probe_rejected':all(probes.values()),
            'missing_dependency_rejected':bool(mismatches(expected,{}))}
    archive=None
    if replay:
        archive=m.OUT/'guarded_replays'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        archive.mkdir(parents=True,exist_ok=False)
        result=m.ROOT/'evidence/material_tangent_mean_diagnostic.json'
        if result.exists():(archive/'previous_mean_diagnostic.json').write_bytes(result.read_bytes())
        m.means()
        (archive/'new_mean_diagnostic.json').write_bytes(result.read_bytes())
        after={name:digest(m.ROOT/name) for name in expected}
        checks['source_inputs_archives_unchanged_during_replay']=not mismatches(expected,after)
        if not checks['source_inputs_archives_unchanged_during_replay']:
            raise RuntimeError('input change during replay; new output is not accepted')
    receipt={'status':'PASS_GUARDED_MEAN_REUSE' if all(checks.values()) else 'FAIL',
             'checks':checks,'dependency_count':len(expected),'expected_sha256':expected,
             'simulated_stale_probes':probes,'replay_performed':replay,
             'archive':None if archive is None else archive.relative_to(m.ROOT).as_posix(),
             'guard_source_sha256':digest(__file__),
             'scope':'Source/input/archive preconditions for the existing floating diagnostic; no extra physical, algebraic-root or parameter-identifiability certificate.'}
    target=m.ROOT/'evidence/material_tangent_guard_verification.json'
    target.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'checks':checks,'dependencies':len(expected),'replay':replay}),flush=True)
    if not all(checks.values()):raise RuntimeError('mean reuse guard verification failed')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--replay',action='store_true');args=parser.parse_args()
    check(args.replay)
