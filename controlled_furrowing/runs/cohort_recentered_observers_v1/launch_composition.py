"""Retain one F/K center-offset observer process; no numerical method changes."""
from pathlib import Path
import datetime,hashlib,json,os,sys,traceback
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import recentered_static_observers as producer

def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic(p,data):
    p.parent.mkdir(parents=True,exist_ok=True);temp=p.with_suffix('.tmp');temp.write_text(json.dumps(data,indent=2)+'\n',encoding='utf8');temp.replace(p)

if __name__=='__main__':
    review='runs/cohort_recentered_observers_v1/prelaunch_independent_review.json'
    jobpath=producer.OUT/'job.json'
    if jobpath.exists():raise RuntimeError('existing observer job retained; do not restart')
    job={'status':'RUNNING','pid':os.getpid(),'started_utc':now(),'launcher_sha256':sha(__file__),
        'source_sha256':sha(producer.__file__),'protocol_sha256':sha(producer.OUT/'protocol.json'),
        'review_file':review,'review_sha256':sha(ROOT/review),'no_root_computation_or_optimization':True}
    atomic(jobpath,job)
    try:
        answer=producer.run(review);job.update(status='FINISHED',finished_utc=now(),results_sha256=sha(producer.OUT/'results.json'),producer_status=answer['status']);atomic(jobpath,job)
    except BaseException:
        job.update(status='FAILED',failed_utc=now(),traceback=traceback.format_exc());atomic(jobpath,job);raise
