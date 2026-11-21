"""Pin all four completed cached-pilot workers before independent cache reads."""
from hashlib import sha256
import json
from pathlib import Path

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[3]
WORK=BASE/'actual'


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_bytes())


def main():
    start=read(WORK/'started.json'); done=read(WORK/'terminal.json')
    assert start['pid']==done['pid']==27828 and done['status']=='FINISHED'
    assert done['results_sha256']==sha(WORK/'results.json')
    ids=start['ordered_pilot_ids']
    assert len(ids)==4 and len(set(ids))==4
    paths=[WORK/x for x in ('started.json','terminal.json','results.json')]
    for name in ids:
        folder=WORK/name
        paths.extend(folder/x for x in ('child_started.json','worker_started.json','observed_os_exit.json','terminal.json','results.json'))
        t,o,s=(read(folder/x) for x in ('terminal.json','observed_os_exit.json','worker_started.json'))
        assert t['status']=='FINISHED' and t['pid']==o['pid']==s['pid']
        assert o['observed_exit_code']==0 and o['supervisor_pid']==27828
        assert o['terminal_sha256']==sha(folder/'terminal.json') and t['results_sha256']==sha(folder/'results.json')
    output={'status':'PASS_ROOT_OBSERVED_CACHED_PILOT_TERMINALS','supervisor_pid':27828,
            'unified_exec_session':46941,'observed_exit_code':0,
            'sha256':{x.relative_to(ROOT).as_posix():sha(x) for x in paths},
            'observation_source_sha256':sha(Path(__file__)),
            'cache_payloads_opened':False,'scope':'Terminal provenance admission only; numerical proof replay remains separate.'}
    assert len(output['sha256'])==23
    p=BASE/'terminal_admission.json'
    with p.open('x',encoding='utf-8') as f:
        json.dump(output,f,indent=2); f.write('\n')
    print(json.dumps({'admission_sha256':sha(p),'metadata_files':len(paths)}))


if __name__=='__main__':
    main()
