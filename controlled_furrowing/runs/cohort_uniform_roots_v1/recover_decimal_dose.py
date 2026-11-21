"""Versioned metadata-only recovery of the frozen cohort driver.

The original files stay unchanged. The sole adapter normalizes the spelling of
the reused M dose AFTER proving exact Fraction equality and checking its hash.
No physical function, parameter, predictor, radius policy or kernel is patched.
Requires a separate hash-bound independent recovery review before execution.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,copy,hashlib,json,os,sys

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
import cohort_uniform_roots as frozen

SOURCE='889182fc142ba34abda7f35758809f0e84e6b7f3d923e42237e247c16f3e7d12'
PROTOCOL='e1a6efd2a3a09611f60419b36845db296f8b75e875a150b8ce2384cf623e0982'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))


def preflight():
    p=read(OUT/'protocol.json');job=read(OUT/'job.json');oldpath=ROOT/p['actualM_reuse']['root_certificate'];old=read(oldpath)
    row=next(x for x in p['records'] if x['mode']=='reuse_immutable_actual_M')
    if sha(Path(frozen.__file__))!=SOURCE or sha(OUT/'protocol.json')!=PROTOCOL:raise RuntimeError('original frozen cohort context changed')
    frozen.assert_frozen(p)
    if job['status']!='FAILED' or job['source_sha256']!=SOURCE or job['protocol_sha256']!=PROTOCOL:raise RuntimeError('specific failed original job required')
    current=frozen.process_identity(job['pid'])
    if current and current['alive'] and current['creation_ticks']==job['creation_ticks']:raise RuntimeError('original job still live')
    if Q(old['dose_exact'])!=Q(row['dose_exact']) or old['dose_exact']==row['dose_exact']:raise RuntimeError('only proven-equal alternative rational spellings may be normalized')
    if sha(oldpath)!=p['actualM_reuse']['root_certificate_sha256']:raise RuntimeError('reused actualM root changed')
    progress=read(OUT/'progress.json')
    if len(progress['completed'])!=11 or any(sha(ROOT/x['file'])!=x['sha256'] for x in progress['completed']):raise RuntimeError('preserved eleven-record prefix changed')
    return p,row,oldpath,{'status':'PASS','checks':{'frozen_original_source_protocol_and_inputs_unchanged':True,'original_job_failed_and_dead':True,
        'alternative_spellings_exactly_equal_as_Fractions':True,'only_reused_M_dose_metadata_is_adapted':True,'eleven_completed_results_preserved_and_hash_checked':True,
        'F_K_unresolved_not_retried':all(x['status']=='UNRESOLVED' for x in progress['completed'] if x['record_id'] in ('optoGEF2_small_col_F','optoGEF2_small_col_K'))},
        'source_sha256':sha(__file__),'protocol_sha256':PROTOCOL,'failed_job_sha256':sha(OUT/'job.json'),
        'preserved_progress_sha256':sha(OUT/'progress.json'),'old_dose_spelling':old['dose_exact'],'canonical_dose_spelling':row['dose_exact'],
        'new_model_or_hessian_calculations':0,'patch_scope':'read-only returned dictionary for ONE immutable actualM root certificate; original bytes and all mathematical targets untouched'}


def run():
    p,row,oldpath,check=preflight();reviewpath=OUT/'recovery_independent_review.json';review=read(reviewpath)
    if review.get('status')!='PASS' or not review.get('checks') or not all(v is True for v in review['checks'].values()):raise RuntimeError('independent recovery review required')
    if review.get('source_sha256')!=sha(__file__) or review.get('protocol_sha256')!=PROTOCOL or review.get('preflight_sha256')!=sha(OUT/'recovery_preflight.json'):raise RuntimeError('recovery review source/protocol/preflight not bound')
    if check!=read(OUT/'recovery_preflight.json'):raise RuntimeError('failed prefix/preflight changed')
    frozen.immutable(OUT/'recovery_launch.json',{'status':'STARTED','pid':os.getpid(),'source_sha256':sha(__file__),
        'independent_review_sha256':sha(reviewpath),'preflight_sha256':sha(OUT/'recovery_preflight.json'),'original_protocol_sha256':PROTOCOL})
    original=frozen.read
    def normalized_read(path):
        value=original(path)
        if Path(path).resolve()==oldpath.resolve():
            if sha(path)!=p['actualM_reuse']['root_certificate_sha256'] or Q(value['dose_exact'])!=Q(row['dose_exact']):raise RuntimeError('M equality/hash guard changed')
            value=copy.deepcopy(value);value['dose_exact']=row['dose_exact']
        return value
    frozen.read=normalized_read
    try:frozen.run(resume=True)
    finally:frozen.read=original
    frozen.immutable(OUT/'recovery_completion.json',{'status':'FINISHED','source_sha256':sha(__file__),'independent_review_sha256':sha(reviewpath),
        'original_driver_sha256':sha(Path(frozen.__file__)),'original_protocol_sha256':sha(OUT/'protocol.json'),
        'results_sha256':sha(OUT/'results.json'),'job_sha256':sha(OUT/'job.json'),'patch_scope':check['patch_scope']})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['preflight','run']);args=parser.parse_args()
    if args.stage=='preflight':
        _,_,_,value=preflight();frozen.immutable(OUT/'recovery_preflight.json',value);print(json.dumps(value,indent=2))
    else:run()
