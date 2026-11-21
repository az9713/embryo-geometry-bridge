"""Record a root-observed real OS exit; authenticate metadata, not mathematics.

Use only after the root has actually polled the stated original exec handle and
observed its terminal exit. This script cannot observe or infer that tool event.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('runs/transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1/bridge_v1')
PILOTS = ('max_contrast_change_small_17','max_contrast_change_large_03',
          'max_negative_contrast_change_small_17','max_negative_contrast_change_large_03')
PRODUCER = Path('scripts/transfer_ambiguity_exact_bridge.py')
PRODUCER_SHA = '58ae92f849414f185c4e28bbb2ec3b9d05361391c4d6e3fa71f519162e30e5e7'


def require(ok,message):
    if not ok:
        raise ValueError(message)


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def read(path):
    return json.loads((ROOT/path).read_text(encoding='utf-8'))


def record(args):
    require(args.observed_exit_code==0 and args.session>0 and args.pid>0,
            'caller must supply an actually observed successful root exec event')
    b,p=BASE/'actual'/args.pilot_id,BASE/'plans'/args.pilot_id
    start,end,manifest=[read(b/(n+'.json')) for n in ('started','terminal','manifest')]
    require(end['status']=='FINISHED' and type(end['pid']) is int and end['pid']==start['pid']==args.pid,
            'same terminal supervisor and explicit root-observed PID required')
    require(end['manifest_sha256']==digest(b/'manifest.json')
            and start['bindings']==end['bindings']==manifest['bindings'], 'terminal manifest/bindings changed')
    require(len(end['workers'])==2 and [w['part'] for w in end['workers']]==['bulk','surface']
            and all(type(w['observed_exit_code']) is int and w['observed_exit_code']==0
                    and type(w['pid']) is int and w['pid']>0 for w in end['workers']),
            'two real successful typed child exits required')
    order=[(k,i) for k,n in [('bulk',144),('surface',24)] for i in range(n)]
    require([(r['kind'],r['slot']) for r in manifest['files']]==order
            and all(type(r['slot']) is int and digest(Path(r['file']))==r['sha256']
                    for r in manifest['files']), 'complete unchanged 168-component payloads required')
    pre=BASE/'preparation_v2'
    review=(BASE/'independent_plan_plus_L/review.json' if args.pilot_id==PILOTS[0] else
            BASE/'independent_plan_remaining_three'/args.pilot_id/'review.json')
    gates={'source_review_file':pre/'prelaunch_independent_review.json','plan_review_file':review,
           'producer_protocol_file':pre/'protocol.json','producer_preflight_file':pre/'preflight.json',
           'producer_frozen_source_file':pre/'transfer_ambiguity_exact_bridge_frozen.py'}
    require(digest(PRODUCER)==digest(gates['producer_frozen_source_file'])==PRODUCER_SHA,
            'same approved producer source required')
    files=[b/(n+'.json') for n in ('started','terminal','manifest')]
    files += [p/(n+'.json') for n in ('plan','started','terminal')]+list(gates.values())+[PRODUCER]
    for part in ('bulk','surface'):
        files += [b/part/(n+'.json') for n in ('started','supervisor_child_started','terminal','observed_os_exit','manifest')]
    hashes={f.as_posix():digest(f) for f in files}
    doc={'status':'PASS_ROOT_OBSERVED_SINGLETON_BRIDGE_COMPONENT_TERMINALS','pilot_id':args.pilot_id,
         'supervisor_pid':args.pid,'unified_exec_session':args.session,'observed_exit_code':0,
         **{k:v.as_posix() for k,v in gates.items()},'sha256':hashes,
         'root_observation_writer_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
         'observation':'Root supplied a successful terminal observation of the original exec handle. '
                       'Metadata and complete payload bytes verified here; mathematical replay remains separate.'}
    target=b/'root_terminal_admission.json'
    with (ROOT/target).open('x',encoding='utf-8') as output:
        json.dump(doc,output,indent=2);output.write('\n')
    print(json.dumps({'admission_file':target.as_posix(),'admission_sha256':digest(target),
                      'terminal_sha256':digest(b/'terminal.json'),'metadata_files':len(files),'packets':168}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pilot-id',choices=PILOTS,required=True)
    parser.add_argument('--session',type=int,required=True)
    parser.add_argument('--pid',type=int,required=True)
    parser.add_argument('--observed-exit-code',type=int,required=True)
    record(parser.parse_args())
