"""Single approved parameter-polynomial job launcher; no physical edits."""
from pathlib import Path
import hashlib,json,subprocess,sys,os,datetime,time
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent
SOURCE=ROOT/'scripts/exact_parameter_defect_B_safe.py'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def write(p,d):
 if p.exists():raise ValueError('write-once launcher artifact exists '+str(p))
 p.write_text(json.dumps(d,indent=2)+'\n')
review=read(OUT/'prelaunch_independent_review.json');checks=review.get('checks')
if sha(OUT/'prelaunch_independent_review.json')!='8a36c79100da725db3f77252f7446cdd1eb56bdf96d0cb0ac2975355bdfc5412'or review.get('status')!='PASS'or not isinstance(checks,dict)or not checks or not all(v is True for v in checks.values()):raise ValueError('actual bound root approval required')
for field,path in [('source_sha256',SOURCE),('protocol_sha256',OUT/'protocol.json'),('preflight_sha256',OUT/'preflight.json')]:
 if review[field]!=sha(path):raise ValueError('changed approved '+field)
if any((OUT/n).exists()for n in ['job.json','started.json','results.json','terminal.json','stdout.log','execution.json']):raise ValueError('sole approved run already started')
cmd=[sys.executable,'-u','-B',str(SOURCE),'run'];start=time.monotonic()
with (OUT/'stdout.log').open('wb')as log:
 child=subprocess.Popen(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
 job={'status':'RUNNING','pid':child.pid,'launcher_pid':os.getpid(),'launch_utc':now(),'command':cmd,'cwd':str(ROOT),'source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json'),'prelaunch_review_sha256':sha(OUT/'prelaunch_independent_review.json'),'launcher_source_sha256':sha(__file__)}
 write(OUT/'job.json',job);print(json.dumps(job),flush=True);exitcode=child.wait()
terminal={'status':'FINISHED'if exitcode==0 else'FAILED','exit_code':exitcode,'pid':child.pid,'launcher_pid':os.getpid(),'finished_utc':now(),'seconds':time.monotonic()-start,'stdout_sha256':sha(OUT/'stdout.log'),'job_sha256':sha(OUT/'job.json'),'results_sha256':sha(OUT/'results.json')if(OUT/'results.json').exists()else None,'producer_terminal_sha256':sha(OUT/'terminal.json')if(OUT/'terminal.json').exists()else None,'all_approved_source_protocol_preflight_unchanged':all(review[k]==sha(v)for k,v in [('source_sha256',SOURCE),('protocol_sha256',OUT/'protocol.json'),('preflight_sha256',OUT/'preflight.json')])}
write(OUT/'execution.json',terminal);print(json.dumps(terminal),flush=True);sys.exit(exitcode)
