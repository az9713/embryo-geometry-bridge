"""Retain a single observer-composition process; no numerical method changes."""
from pathlib import Path
import argparse,datetime,hashlib,json,os,sys,traceback

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import cohort_static_observers as producer

def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic(p,data):
    p.parent.mkdir(parents=True,exist_ok=True)
    temp=p.with_suffix('.tmp');temp.write_text(json.dumps(data,indent=2)+'\n',encoding='utf8');temp.replace(p)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--selection',required=True);args=parser.parse_args()
    review='runs/cohort_static_observers_v1/prelaunch_independent_review.json'
    out=ROOT/'runs/cohort_static_observers_v1/composition_final'
    jobpath=out/'job.json'
    if jobpath.exists():raise RuntimeError('existing job retained; do not restart this launch')
    job={'status':'RUNNING','pid':os.getpid(),'started_utc':now(),
        'launcher_sha256':digest(__file__),'producer_source_sha256':digest(producer.__file__),
        'observer_protocol_sha256':digest(producer.OUT/'protocol.json'),
        'review_file':review,'review_sha256':digest(ROOT/review),
        'selection_file':args.selection,'selection_sha256':digest(ROOT/args.selection),
        'no_equilibrium_optimization':True}
    atomic(jobpath,job)
    try:
        answer=producer.run(review,args.selection,'composition_final')
        job.update(status='FINISHED',finished_utc=now(),results_sha256=digest(out/'results.json'),producer_status=answer['status'])
        atomic(jobpath,job)
    except BaseException:
        job.update(status='FAILED',failed_utc=now(),traceback=traceback.format_exc());atomic(jobpath,job);raise
