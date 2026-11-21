"""Root observation of completed per-record centered worker and OS receipts."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CALLER = 'runs/exact_integration_cohort_v1/remaining15_centered_actual_preparation'
LAUNCH = CALLER+'/launches/77076fd9deee4f63df93c3a0411b05820c1b0769101d85619d73daf7ac7c7442'


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def main(letter):
    assert letter in ('C', 'D')
    record = 'optoGEF2_small_col_'+letter
    paths = {
        'terminal': f'{CALLER}/actual/{record}/terminal.json',
        'manifest': f'{CALLER}/actual/{record}/results.json',
        'os_terminal': f'{LAUNCH}/{record}.os_terminal.json',
        'child_started': f'{LAUNCH}/{record}.started.json',
    }
    docs = {k: json.loads((ROOT/v).read_bytes()) for k,v in paths.items()}
    t, m, o, s = (docs[k] for k in paths)
    assert t['status'] == 'FINISHED' and t['exit_code'] == o['exit_code'] == 0
    assert t['pid'] == o['pid'] == s['pid']
    assert t['results_sha256'] == sha(ROOT/paths['manifest'])
    assert m['record_id'] == record and len(m['files']) == 144 and not m['root_asserted']
    pins = {k+'_sha256': sha(ROOT/v) for k,v in paths.items()}
    pins['worker_review_sha256'] = '31ed5e5a15e73242414a44ed3b9412679d01d0e590eedbf4d94973dcbb556ca2'
    pins['supervisor_dir'] = LAUNCH
    output = {'status': 'ROOT_OBSERVED_COMPLETE_CENTERED_RECORD_BEFORE_PACKET_READ',
              'record_id': record, 'pid': t['pid'], 'unified_exec_session': 30794,
              'source_sha256': sha(Path(__file__)), 'paths': paths, 'pins': pins,
              'local_component_payloads_opened': False}
    target = Path(__file__).parent/f'{letter}_centered_terminal_observation.json'
    with target.open('x', encoding='utf-8') as f:
        json.dump(output, f, indent=2)
        f.write('\n')
    print(json.dumps({'observation_sha256': sha(target), 'pins': pins}))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('letter', choices=['C','D'])
    main(p.parse_args().letter)
