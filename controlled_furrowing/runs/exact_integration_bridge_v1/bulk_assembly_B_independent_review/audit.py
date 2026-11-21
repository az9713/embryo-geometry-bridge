"""Exact144-packet reduction; only preselected48/143 local reintegrations."""
from pathlib import Path
from fractions import Fraction as F
import argparse,base64,hashlib,importlib.util,io,json,sys,time
import numpy as np
from flint import arb,ctx
HERE=Path(__file__).resolve().parent;CF=HERE.parents[2];PRODUCER=CF/'runs/exact_integration_bridge_v1/bulk_assembly_B'
HELPER=CF/'runs/exact_integration_bridge_v1/bulk_pilot_independent_replay/replay.py'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):
    p=Path(p)
    if p.exists():raise ValueError('write-once output exists: '+str(p))
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n')
def kernel():
    if sha(HELPER)!='82f19f5aaf7a2ae986768fb2450c882e6b625cd6a03f5e0a6846ff215b927686':raise ValueError('frozen independent kernel changed')
    spec=importlib.util.spec_from_file_location('independent_bulk_point_helper',HELPER);m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
def check_inputs(p):
    for name,d in p['input_hashes'].items():
        if sha(CF/name)!=d:raise ValueError('changed input '+name)
def prepare():
    I=kernel();ctx.prec=640;I.C.A.cache_clear();checks=I.selftests();p=read(PRODUCER/'protocol_v2.json');inputs={}
    if sha(CF/'scripts/exact_bulk_bridge_assembly.py')!='41fc1735566f72ec42e249aa00fbc61f9e26b9ef55ea96f3e790bb82adac83ff'or sha(PRODUCER/'protocol_v2.json')!='8b14c671ca877d96256225014dcc42fef97df6ea6f7cabc3af860b857398cb26':raise ValueError('producer context changed')
    for name,d in p['source_and_input_hashes'].items():
        if sha(CF/name)!=d:raise ValueError('producer input changed')
        key=(CF/name).resolve().relative_to(CF.resolve()).as_posix()
        if key in inputs and inputs[key]!=d:raise ValueError('conflicting input hashes')
        inputs[key]=d
    for f in[Path(__file__),HELPER,PRODUCER/'protocol_v2.json',PRODUCER/'preflight_v2.json',CF/'evidence/exact_bulk_pilot_independent_review.json']:
        inputs[f.resolve().relative_to(CF.resolve()).as_posix()]=sha(f)
    protocol={'source_sha256':sha(Path(__file__)),'input_hashes':inputs,'additional_elements_predeclared':[48,143],'additional_order':18,'additional_precision_bits':640,'additional_rho_exact':'21/500','selfchecks':checks,'coverage':'Exact all144 packet reduction after producerterminal, but only additional48/143 local integrations. Prior four independent local replays are reused transparently.','producer_results_not_read_before_freeze':True,'full_tube_or_root_transfer_claim':False}
    write(HERE/'protocol.json',protocol);write(HERE/'selfchecks.json',{'status':'PASS','checks':checks,'protocol_sha256':sha(HERE/'protocol.json')});print(json.dumps({'status':'PREPARED_ONLY','source_sha256':sha(Path(__file__)),'protocol_sha256':sha(HERE/'protocol.json')}),flush=True)
def sample_element(I,target,index,rule,record):
    start=time.monotonic();f=target.element(index);jet=I.independent_jet(f['q'],f['B']);n=jet['local_dofs'];fun=I.Functional(target,f,index,rule);a,b=I.series(I.T.sub(jet['value'],I.T.ONE),18);vol=I.DUAL.contract_volumetric_defect(jet,a,b,fun.defect);drive=F(record['base_strength_exact'])*I.rational(target.a['active'][index,0]);diagonal=[F(1),1+drive,F(1),1+drive,F(1)];t1,t2=I.M.volumetric_derivative_tail_bounds(18,F(21,500),1,5000);mass=I.C.upper(fun.massI+fun.massQ);jg=[I.majorant(x,f)for x in jet['gradient']]
    force=[];hessian={}
    for i in range(n):
        p={}
        for k in range(5):p=I.T.add(p,I.T.scale(I.T.mul(f['q'][k],f['B'][k][i]),diagonal[k]))
        tail=mass*I.C.A(t1*jg[i]);force.append(fun.apply(p)+vol['gradient'][i]+(-tail).union(tail))
    for i,j in jet['hessian_upper']:
        p={}
        for k in range(5):p=I.T.add(p,I.T.scale(I.T.mul(f['B'][k][i],f['B'][k][j]),diagonal[k]))
        tail=mass*I.C.A(t2*jg[i]*jg[j]+t1*I.majorant(jet['hessian_upper'][i,j],f));hessian[i,j]=fun.apply(p)+vol['hessian_upper'][i,j]+(-tail).union(tail)
    cg,ch,count=I.correction(target,f,index,fun,drive);force=[v+cg[i]for i,v in enumerate(force)];hessian={k:v+ch[k]for k,v in hessian.items()}
    return {'element':index,'order':18,'precision_bits':640,'global_indices':f['global_indices'],'chart':f['chart']['kind'],'canonical_basis_entry_checks':count,'force':[I.C.pack(v)for v in force],'H_upper':[{'i':i,'j':j,'interval':I.C.pack(v)}for(i,j),v in hessian.items()],'seconds':time.monotonic()-start,'scope':'Only preselected POINT element; no global/tube/root claim.'}
def samples():
    p=read(HERE/'protocol.json');check_inputs(p);I=kernel();ctx.prec=640;I.C.A.cache_clear();I.selftests();record=read(CF/'runs/cohort_uniform_roots_v2/records/00_optoGEF2_small_col_B/result.json');target=I.T.Target(read(CF/record['point_cache_file']));rule=target.auxiliary_rule();rows=[]
    write(HERE/'samples_started.json',{'protocol_sha256':sha(HERE/'protocol.json'),'indices':[48,143]})
    for index in p['additional_elements_predeclared']:
        result=sample_element(I,target,index,rule,record);file=HERE/f'additional_{index:03d}.json';write(file,result);rows.append({'element':index,'file':file.name,'sha256':sha(file)});print(json.dumps({'element':index,'seconds':result['seconds']}),flush=True)
    check_inputs(p);write(HERE/'samples_results.json',{'status':'ADDITIONAL_TWO_POINT_REPLAYS_COMPLETE','elements':rows,'all_frozen_inputs_unchanged':True,'protocol_sha256':sha(HERE/'protocol.json')})
def interval(p):
    lo,hi=F(p['lower_exact']),F(p['upper_exact'])
    if lo>hi:raise ValueError('reversed interval')
    return lo,hi
def add(a,b):return a[0]+b[0],a[1]+b[1]
def scalar(x,a):return(x*a[0],x*a[1])if x>=0 else(x*a[1],x*a[0])
def maxabs(a):return max(map(abs,a))
def contains(p,a):
    lo,hi=interval(p);return lo<=a[0]and a[1]<=hi
def overlap(a,b):return a[0]<=b[1]and b[0]<=a[1]
def finalize(expected):
    if len(expected)!=64 or sha(PRODUCER/'results.json')!=expected:raise ValueError('authoritative producer terminal SHA required')
    protocol=read(HERE/'protocol.json');check_inputs(protocol);root=read(PRODUCER/'results.json');p=read(PRODUCER/'protocol_v2.json');sample=read(HERE/'samples_results.json');checks={};record=read(CF/'runs/cohort_uniform_roots_v2/records/00_optoGEF2_small_col_B/result.json');cachefile=CF/record['point_cache_file'];cache=read(cachefile);blob=base64.b64decode(cache['frozen_binary_cache']['data'],validate=True)
    if hashlib.sha256(blob).hexdigest()!=cache['frozen_binary_cache']['sha256']:raise ValueError('R cache corrupt')
    with np.load(io.BytesIO(blob),allow_pickle=False)as z:R=z['preconditioner'].copy()
    if R.shape!=(611,611)or not np.isfinite(R).all():raise ValueError('actual611preconditioner required')
    checks['producer_true_complete_all144']=root['status']=='ALL144_POINT_BULK_ASSEMBLED_UNREVIEWED'and len(root['elements'])==144 and root['all_frozen_bindings_unchanged']is True
    selectors={x['element']:x for x in p['selectors']};Fsum=[(F(0),F(0))for _ in range(611)];Hsum={};files={};tails=[]
    for item in root['elements']:
        file=CF/item['file']
        if sha(file)!=item['sha256']:raise ValueError('producer element changed')
        x=read(file);i=x['element'];sel=selectors[i]
        if i in files:raise ValueError('duplicate producer element')
        files[i]=item;ids=x['global_indices'];n=len(ids)
        if ids!=sel['global_indices']or len(set(ids))!=n or any(not 0<=k<611 for k in ids)or len(x['force_integral_minus_cache'])!=n:raise ValueError('local free map/force mismatch')
        if x['truncation_order']!=sel['degree']or F(x['rho_exact'])!=F(sel['point_rho_exact']):raise ValueError('local degree/rho mismatch')
        ht={(a['i'],a['j']):a['enclosure']for a in x['H_integral_minus_cache_upper_triangle']}
        if len(ht)!=len(x['H_integral_minus_cache_upper_triangle'])or set(ht)!={(a,b)for a in range(n)for b in range(a,n)}:raise ValueError('missing/duplicate local Hessian slot')
        for k,value in zip(ids,x['force_integral_minus_cache']):Fsum[k]=add(Fsum[k],interval(value))
        for(a,b),v in ht.items():
            key=tuple(sorted((ids[a],ids[b])));Hsum[key]=add(Hsum.get(key,(F(0),F(0))),interval(v))
        good=interval(x['force_max_tail_upper'])[1]<=F(p['force_tail_target_exact'])and interval(x['H_max_tail_upper'])[1]<=F(p['H_entry_tail_target_exact']);tails.append({'element':i,'PASS':good});checks[f'element{i:03d}_tails_degree_rho_hash_free_map']=good
    checks['all144_element_IDs_exactly_once']=set(files)==set(range(144))
    checks['all611_force_intervals_contain_exact_endpoint_sums']=len(root['force_integral_minus_cache'])==611 and all(contains(v,s)for v,s in zip(root['force_integral_minus_cache'],Fsum))
    producerH={(a['i'],a['j']):a['enclosure']for a in root['bulk_H_upper']};checks['all6676_Hslots_complete_unique_and_contain_exact_sums']=len(producerH)==len(root['bulk_H_upper'])==len(Hsum)==6676 and set(producerH)==set(Hsum)and all(contains(producerH[k],v)for k,v in Hsum.items())
    rf=[];Rnorm=F(0)
    for row in R:
        exact=[F.from_float(float(x))for x in row];v=(F(0),F(0))
        for a,b in zip(exact,Fsum):v=add(v,scalar(a,b))
        rf.append(v);Rnorm=max(Rnorm,sum(map(abs,exact),F(0)))
    checks['all611_signed_RF_intervals_contain_exact_matrix_endpoint_sums']=len(root['preconditioned_bulk_force_defect'])==611 and all(contains(p,v)for p,v in zip(root['preconditioned_bulk_force_defect'],rf))
    hrows=[F(0)]*611
    for(i,j),v in Hsum.items():
        hrows[i]+=maxabs(v)
        if i!=j:hrows[j]+=maxabs(v)
    fn=max(map(maxabs,Fsum));rfn=max(map(maxabs,rf));hn=max(hrows)
    checks['reported_F_norm_bounds_exact_reduction']=interval(root['bulk_force_norm_upper'])[1]>=fn
    checks['reported_RF_norm_bounds_exact_reduction']=interval(root['preconditioned_bulk_force_norm_upper'])[1]>=rfn
    checks['reported_symmetric_H_row_norm_bounds_exact_reduction']=interval(root['bulk_H_row_norm_upper'])[1]>=hn
    checks['reported_conservative_RH_norm_bounds_exact_Rnorm_times_Hnorm']=interval(root['conservative_preconditioned_bulk_H_norm_upper'])[1]>=Rnorm*hn
    sample_rows=[]
    for item in sample['elements']:
        file=HERE/item['file']
        if sha(file)!=item['sha256']:raise ValueError('independent additional sample changed')
        x=read(file);old=read(CF/files[x['element']]['file']);ph={(v['i'],v['j']):interval(v['enclosure'])for v in old['H_integral_minus_cache_upper_triangle']};nh={(v['i'],v['j']):interval(v['interval'])for v in x['H_upper']}
        good=x['global_indices']==old['global_indices']and len(x['force'])==len(old['force_integral_minus_cache'])and all(overlap(interval(a),interval(b))for a,b in zip(x['force'],old['force_integral_minus_cache']))and set(ph)==set(nh)and all(overlap(ph[k],v)for k,v in nh.items());checks[f'additional{x["element"]}_all_force_H_overlap']=good
        sample_rows.append({'element':x['element'],'independent_file':item['file'],'independent_sha256':item['sha256'],'force_entries':len(x['force']),'H_entries':len(nh),'exact_basis_entry_checks':x['canonical_basis_entry_checks'],'PASS':good})
    checks['only_predeclared_additional_indices48_143']=sorted(x['element']for x in sample_rows)==[48,143]
    checks['all_inputs_still_frozen']=all(sha(CF/k)==v for k,v in protocol['input_hashes'].items())
    exact={'force_endpoint_sums':[[str(a),str(b)]for a,b in Fsum],'H_upper_endpoint_sums':[{'i':i,'j':j,'interval':[str(a),str(b)]}for(i,j),(a,b)in sorted(Hsum.items())],'RF_endpoint_sums':[[str(a),str(b)]for a,b in rf],'R_infinity_norm_exact':str(Rnorm),'force_norm_exact_upper':str(fn),'RF_norm_exact_upper':str(rfn),'H_row_norm_exact_upper':str(hn),'RH_conservative_exact_upper':str(Rnorm*hn)}
    write(HERE/'exact_reduction.json',exact)
    result={'status':'PASS_INDEPENDENT_POINT_BULK_ASSEMBLY'if all(checks.values())else'FAIL_INDEPENDENT_POINT_BULK_ASSEMBLY','checks':checks,'protocol_sha256':sha(HERE/'protocol.json'),'review_source_sha256':sha(Path(__file__)),'producer_results_sha256':expected,'producer_protocol_sha256':sha(PRODUCER/'protocol_v2.json'),'exact_reduction_sha256':sha(HERE/'exact_reduction.json'),'all144_local_packets':list(files.values()),'additional_samples':sample_rows,'tail_checks':tails,'exact_norms':{k:v for k,v in exact.items()if k.endswith('_exact')or k.endswith('_upper')},'local_reintegration_coverage':{'new_indices':[48,143],'prior_indices':[0,6,72,66],'total_distinct_independently_reintegrated':6,'other138':'Use mathematically reviewed generic localengine plus complete localpacket/tail/hash/free-map and exactglobalreduction checks; not individually independently reintegrated.'},'scope':'ALL144Bcenter BULK assembly only. No surface/tube/parameteruniformity/equilibriumtransfer/PDE/sourcevalidation.','root_transfer_proved':False}
    write(HERE/'results.json',result);write(CF/'evidence/exact_bulk_assembly_B_independent_review.json',result);print(json.dumps({'status':result['status'],'checks':len(checks),'failed':[k for k,v in checks.items()if not v],'receipt_sha256':sha(CF/'evidence/exact_bulk_assembly_B_independent_review.json')}),flush=True)
    if not all(checks.values()):raise SystemExit(1)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','samples','finalize']);p.add_argument('--producer-results-sha');a=p.parse_args()
    prepare()if a.command=='prepare'else samples()if a.command=='samples'else finalize(a.producer_results_sha)
