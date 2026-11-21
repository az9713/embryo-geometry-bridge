"""Sequential supervisor for separately prelaunch-reviewed domain jobs only."""
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent / 'domains_job'
MANIFEST = Path(__file__).with_name('manifest.json')
MANIFEST_SHA = '93a2f5858a450ca9a8db9f9bc944b918e2a7446abba614f97de5e47fe434495d'
DRIVER_SHA = 'a21b680117f236cc38efc46e086be1d9aff80cfb09ae613691b1642db8bed25e'


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def write(name, value):
    with (OUT / name).open('x', encoding='utf8') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def utc():
    return datetime.now(timezone.utc).isoformat()


def main():
    assert sha(MANIFEST) == MANIFEST_SHA
    plan = read(MANIFEST)
    assert len(plan['records']) == 15
    assert [row['record_id'] for row in plan['records']] == plan['expected_record_ids']
    assert len(set(plan['expected_record_ids'])) == 15
    assert set(plan['excluded_completed_records']) == {'optoGEF2_small_col_B', 'optoGEF2_small_col_L'}
    driver = ROOT / 'scripts/cohort_exact_pilot.py'
    assert sha(driver) == DRIVER_SHA
    gates = []
    # Validate every gate and every incoming layout before launching any child.
    for row in plan['records']:
        files = row['files']
        for item in files.values():
            assert sha(ROOT / item['file']) == item['sha256']
        folder = (ROOT / files['protocol.json']['file']).parent
        assert folder.name == 'original_cube_v1' and folder.parent.name == row['record_id']
        gate_path = folder / 'domains_prelaunch_independent_review.json'
        gate = read(gate_path)
        assert gate['status'] == 'PASS' and gate['checks'] and all(v is True for v in gate['checks'].values())
        assert gate['source_sha256'] == DRIVER_SHA
        assert gate['protocol_sha256'] == files['protocol.json']['sha256']
        assert gate['preflight_sha256'] == files['preflight.json']['sha256']
        assert not (folder / 'domains_started.json').exists(), 'never duplicate a started domain job'
        gates.append({'record_id': row['record_id'], 'file': gate_path.relative_to(ROOT).as_posix(),
                      'sha256': sha(gate_path)})
    OUT.mkdir(exist_ok=False)
    started = time.monotonic()
    write('started.json', {'pid': os.getpid(), 'utc': utc(), 'manifest_sha256': MANIFEST_SHA,
                           'driver_sha256': DRIVER_SHA, 'supervisor_sha256': sha(Path(__file__)),
                           'separately_verified_domains_only_gates': gates})
    completed = []
    try:
        for index, row in enumerate(plan['records']):
            record = row['record_id']
            folder = (ROOT / row['files']['protocol.json']['file']).parent
            gate = gates[index]
            assert sha(ROOT / gate['file']) == gate['sha256']
            args = [sys.executable, '-B', '-u', str(driver), 'domains', record]
            with (OUT / f'{record}.stdout.log').open('wb') as stdout, \
                    (OUT / f'{record}.stderr.log').open('wb') as stderr:
                process = subprocess.Popen(args, cwd=ROOT.parent, stdout=stdout, stderr=stderr)
                write(f'{record}.started.json', {'pid': process.pid, 'utc': utc(), 'command': args})
                print(json.dumps({'record': record, 'started_pid': process.pid}), flush=True)
                code = process.wait()
            # Preserve an OS exit even when the child failed before its receipt.
            result = {'record_id': record, 'pid': process.pid, 'exit_code': code, 'utc': utc()}
            write(f'{record}.os_terminal.json', result)
            assert code == 0, f'{record}: child failed; inspect retained logs'
            terminal = read(folder / 'domains_terminal.json')
            assert terminal['exit_code'] == 0 and terminal['pid'] == process.pid
            assert terminal['results_sha256'] == sha(folder / 'domains.json')
            result.update({'domains_sha256': terminal['results_sha256'],
                           'terminal_sha256': sha(folder / 'domains_terminal.json')})
            completed.append(result)
            print(json.dumps({'record': record, 'complete': len(completed), 'total': 15}), flush=True)
        write('results.json', {'status': 'COMPLETE_15_DOMAINS_ONLY_UNREVIEWED', 'records': completed,
                              'force_H_integrals_run': False, 'seconds': time.monotonic() - started})
        write('terminal.json', {'status': 'FINISHED', 'exit_code': 0, 'pid': os.getpid(),
                               'results_sha256': sha(OUT / 'results.json'), 'utc': utc()})
    except BaseException as error:
        write('failure.json', {'error': str(error), 'completed': completed})
        write('terminal.json', {'status': 'FAILED', 'exit_code': 1, 'pid': os.getpid(), 'utc': utc()})
        raise


if __name__ == '__main__':
    main()
