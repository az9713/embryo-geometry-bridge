"""Bind root-observed completed nonlinear workers before independent field reads."""
from hashlib import sha256
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
WORK = BASE/'actual'
ENDPOINTS = ('max_contrast_change', 'max_negative_contrast_change')


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def main():
    names = [WORK/x for x in ('started.json','terminal.json','results.json')]
    for name in ENDPOINTS:
        names.extend(WORK/name/x for x in ('supervisor_child_started.json','worker_started.json',
                                          'supervisor_os_exit.json','terminal.json','results.json'))
    assert len(names) == 13
    top = read(WORK/'terminal.json')
    assert top['status'] == 'FINISHED' and top['pid'] == 35436
    assert top['results_sha256'] == sha(WORK/'results.json')
    for name in ENDPOINTS:
        folder = WORK/name
        term, os_exit, start = (read(folder/x) for x in ('terminal.json','supervisor_os_exit.json','worker_started.json'))
        assert term['status'] == 'FINISHED'
        assert term['pid'] == os_exit['pid'] == start['pid']
        assert os_exit['observed_exit_code'] == 0 and os_exit['supervisor_pid'] == 35436
        assert os_exit['terminal_sha256'] == sha(folder/'terminal.json')
        assert term['results_sha256'] == sha(folder/'results.json')
    receipt = {'status':'PASS_ROOT_OBSERVED_NONLINEAR_TERMINAL', 'supervisor_pid':35436,
               'observed_exit_code':0, 'unified_exec_session':90479,
               'sha256':{x.relative_to(ROOT).as_posix():sha(x) for x in names},
               'observation_source_sha256':sha(Path(__file__)),
               'field_payloads_opened':False,
               'scope':'Root observed unified-exec exit0 and both true subprocess OS exits; all field validation remains to run.'}
    path = BASE/'terminal_admission.json'
    with path.open('x',encoding='utf-8') as f:
        json.dump(receipt,f,indent=2)
        f.write('\n')
    print(json.dumps({'admission_sha256':sha(path),'files':len(names)}))


if __name__=='__main__':
    main()
