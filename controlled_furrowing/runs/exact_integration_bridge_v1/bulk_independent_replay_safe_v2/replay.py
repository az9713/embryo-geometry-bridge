"""One frozen six-element SAFE640bit point replay; no optimization."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,datetime,time
from flint import arb,ctx
import kernel as I
from sample import sample_element
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PRODUCER=ROOT/'runs/exact_integration_bridge_v1/bulk_assembly_B_safe_v2'
SURFACE=ROOT/'runs/exact_integration_bridge_v1/surface_assembly_B_safe_v2'
INDICES=[0,6,72,66,48,143]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):
    if p.exists():raise ValueError('preserve '+str(p))
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n',encoding='utf8')
def selfchecks():
    oldprec=ctx.prec;ctx.prec=128;I.C.A.cache_clear();checks=I.selftests();x=arb(0,1)
    checks['safe_T889b']=sha(Path(I.T.__file__))=='889b61cfa794d432021fbd4670b781169033d4b8451554b18e7f7a569d9d3954'
    checks['safe_helper1439']=sha(Path(I.DUAL.__file__))=='1439e20ec127d39565467d0788173b1148b5ab47f659b3b39486b81aaf1ee373'
    checks['uncertain_constant_retained']=(0,0)in I.T.constant(x)
    checks['uncertain_product_retained']=(1,0)in I.T.mul({(1,0):x},I.T.ONE)
    checks['uncertain_sum_retained']=(0,0)in I.T.add({(0,0):arb(1)},{(0,0):arb(-1,1)})
    exact=F(5001,2**1000);aa,bb=I.series({(0,0):F(1,5001)+F(1,2**1000)},1)
    checks['independent_signed_convolution_near_cancelled_constant_retained']=(0,0)in aa and aa[(0,0)].contains(I.C.A(exact))and not aa[(0,0)].is_zero()
    checks['own_kernel_has_no_unsafe_neq_cleanup']='if v!=0'not in (HERE/'kernel.py').read_text()
    assert all(v is True for v in checks.values()),checks
    ctx.prec=oldprec;I.C.A.cache_clear();return checks
def prepare(expected):
    if not expected or sha(PRODUCER/'results.json')!=expected:raise ValueError('parent-authoritative TRUE terminal producer SHA required')
    terminal=read(PRODUCER/'results.json');assert terminal['status']=='SAFE_V2_ALL144_POINT_BULK_ASSEMBLED_UNREVIEWED'and len(terminal['elements'])==144
    ctx.prec=640;tests=selfchecks();p=read(PRODUCER/'protocol_v2.json');bindings=dict(p['source_and_input_hashes']);selected={x['element']:x for x in terminal['elements']}
    for index in INDICES:
        item=selected[index];assert sha(ROOT/item['file'])==item['sha256'];bindings[item['file'].replace('\\','/')]=item['sha256']
    for path in [Path(__file__),HERE/'kernel.py',HERE/'sample.py',HERE/'code_lineage.json',PRODUCER/'results.json',PRODUCER/'protocol_v2.json',PRODUCER/'preflight_v2.json',PRODUCER/'prelaunch_independent_review.json',SURFACE/'results.json',ROOT/'evidence/exact_safe_polynomial_kernel_independent_review.json']:
        bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    assert all(sha(ROOT/k)==v for k,v in bindings.items())
    proto={'version':2,'status':'FROZEN_BEFORE_SIX_REPLAY','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':sha(__file__),'input_hashes':bindings,'elements':INDICES,'order':18,'precision_bits':640,'rho_exact':'21/500','producer_results_sha256':expected,'producer_protocol_sha256':sha(PRODUCER/'protocol_v2.json'),'selfchecks':tests,'scope':'Only six POINT bulk local integrations plus exact signed reductions later. Reused safe sparse arithmetic/DUAL, exact moments and constitutive core; independent cubicJ, coefficient convolution, direct canonical maps/manual lifts and selected direct moment contractions. No parameter tube/PDE/root transfer.'}
    write(HERE/'protocol.json',proto);print(json.dumps({'status':proto['status'],'protocol_sha256':sha(HERE/'protocol.json'),'source_sha256':sha(__file__)}),flush=True)
def run():
    p=read(HERE/'protocol.json');assert sha(__file__)==p['source_sha256']and all(sha(ROOT/k)==v for k,v in p['input_hashes'].items());ctx.prec=640;I.C.A.cache_clear();selfchecks()
    write(HERE/'started.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'protocol_sha256':sha(HERE/'protocol.json')})
    rec=read(ROOT/'runs/cohort_uniform_roots_v2/records/00_optoGEF2_small_col_B/result.json');target=I.T.Target(read(ROOT/rec['point_cache_file']));rule=target.auxiliary_rule();rows=[];checks={};start=time.monotonic()
    for index in p['elements']:
        x=sample_element(I,target,index,rule,rec);prod=read(PRODUCER/'elements'/f'{index:03d}.json');checks[f'element{index}_same_free_map']=x['global_indices']==prod['global_indices']
        checks[f'element{index}_all_force_intervals_overlap']=len(x['force'])==len(prod['force_integral_minus_cache'])and all(I.overlaps(I.C.A(F(a['lower_exact'])).union(I.C.A(F(a['upper_exact']))),b)for a,b in zip(x['force'],prod['force_integral_minus_cache']))
        oldH={(z['i'],z['j']):z['enclosure']for z in prod['H_integral_minus_cache_upper_triangle']};newH={(z['i'],z['j']):z['interval']for z in x['H_upper']};checks[f'element{index}_all_H_intervals_overlap']=set(oldH)==set(newH)and all(I.overlaps(I.C.A(F(v['lower_exact'])).union(I.C.A(F(v['upper_exact']))),oldH[k])for k,v in newH.items())
        for k,v in x['checks'].items():checks[f'element{index}_{k}']=v
        dest=HERE/f'element_{index:03d}.json';write(dest,x);rows.append({'element':index,'file':dest.name,'sha256':sha(dest),'force_entries':len(x['force']),'H_entries':len(x['H_upper']),'canonical_basis_entry_checks':x['canonical_basis_entry_checks']});print(json.dumps({'element':index,'seconds':x['seconds'],'all_checks_pass':all(checks.values())}),flush=True)
    checks['all_frozen_inputs_current']=all(sha(ROOT/k)==v for k,v in p['input_hashes'].items());checks['exactly_six_predeclared_indices']=[x['element']for x in rows]==INDICES
    result={'status':'PASS_SAFE_SIX_POINT_BULK_REPLAYS'if all(checks.values())else'FAIL_SAFE_SIX_POINT_BULK_REPLAYS','checks':checks,'protocol_sha256':sha(HERE/'protocol.json'),'source_sha256':sha(__file__),'elements':rows,'seconds':time.monotonic()-start,'producer_results_sha256':p['producer_results_sha256'],'scope':p['scope'],'old_v1_integration_certificates_restored':False,'root_transfer_proved':False}
    write(HERE/'results.json',result);write(ROOT/'evidence/exact_bulk_independent_replay_safe_v2.json',result);print(json.dumps({'status':result['status'],'seconds':result['seconds'],'receipt_sha256':sha(ROOT/'evidence/exact_bulk_independent_replay_safe_v2.json')}),flush=True)
    assert all(v is True for v in checks.values()),checks
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('command',choices=['selfcheck','prepare','run']);a.add_argument('--producer-results-sha');args=a.parse_args()
    if args.command=='selfcheck':write(HERE/'initial_selfchecks.json',{'status':'PASS','checks':selfchecks(),'source_sha256':sha(__file__),'actual_integrals_run':False});print('SELF_CHECKS_PASS_ONLY')
    elif args.command=='prepare':prepare(args.producer_results_sha)
    else:run()
