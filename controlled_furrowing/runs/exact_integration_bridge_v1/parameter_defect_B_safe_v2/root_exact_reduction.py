"""Independent integer-dyadic reduction of all48 local parameter packets.

No nonlinear/polynomial integrals are replayed. Local reviewed enclosure
engines and the audited basepoint packets remain explicit premises.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,sys,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import exact_fe_target_safe as T
def read(p):return json.loads(p.read_text(encoding='utf8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def interval(p):return Q(p['lower_exact']),Q(p['upper_exact'])
def add(a,b):return a[0]+b[0],a[1]+b[1]
def scale(a,c):return tuple(sorted([a[0]*c,a[1]*c]))
def sub(a,b):return a[0]-b[1],a[1]-b[0]
def included(a,b):return b[0]<=a[0]<=a[1]<=b[1]
checks={}
def ck(k,v):
    checks[k]=bool(v)
    assert v,k
start=time.monotonic();p=read(HERE/'protocol.json');d=read(HERE/'results.json');terminal=read(HERE/'terminal.json')
ck('literal_producer_result',sha(HERE/'results.json')=='0de3472a8b9293bbae397d854d0ec2b6dfafcc7b20f8e54b0e8d236f19f86164')
ck('actual_terminal_bound',terminal['status']=='FINISHED'and terminal['exit_code']==0 and terminal['results_sha256']==sha(HERE/'results.json'))
ck('bound_source_protocol',d['source_sha256']==sha(ROOT/'scripts/exact_parameter_defect_B_safe.py')and d['protocol_sha256']==sha(HERE/'protocol.json'))
for name,h in p['source_and_input_hashes'].items():ck('binding_'+name,sha(ROOT/name)==h)
target=T.Target(read(ROOT/p['point_cache_file']));order=p['parameter_order'];columns=[[(Q(0),Q(0))for _ in range(611)]for _ in range(5)]
surfaces=set();active=set()
for item in d['local_packets']:
    path=ROOT/item['file'];ck('local_hash_'+item['file'],sha(path)==item['sha256']);packet=read(path);ids=packet['global_indices']
    if item['kind']=='surface':
        face,edge=packet['face'],packet['edge'];surfaces.add((face,edge));fields=target.surface(face,edge)
        case=next(v for v in p['surface_edges']if(v['face'],v['edge'])==(face,edge));base=read(ROOT/case['file'])
        ck(f'audited_base_hash_{face}_{edge}',packet['source_base_packet_sha256']==sha(ROOT/case['file']))
        K,D=Q(case['base_K_exact']),Q(case['T_minus_K_exact'])
        for j,(F,P)in enumerate(zip(base['force_integral_minus_cache'],packet['poly_force_defect'])):
            q=scale(sub(interval(F),scale(interval(P),K)),1/D)
            ck(f'local_T_algebra_{face}_{edge}_{j}',included(q,interval(packet['partial_T'][j])))
            ck(f'local_K_algebra_{face}_{edge}_{j}',included(sub(interval(P),q),interval(packet['partial_K'][j])))
        column_specs=[(1 if face==0 else 0,packet['partial_K']),(3 if face==0 else 2,packet['partial_T'])]
    elif item['kind']=='bulk_gain':
        element=packet['element'];active.add(element);fields=target.element(element)
        ck(f'active_mask_{element}',fields['active']==1)
        column_specs=[(4,packet['gain_partial_force_defect'])]
    else:raise AssertionError('unknown local packet')
    ck('free_map_'+item['file'],ids==fields['global_indices'])
    for col,values in column_specs:
        ck('local_count_'+item['file']+'_'+str(col),len(values)==len(ids))
        for i,v in zip(ids,values):columns[col][i]=add(columns[col][i],interval(v))
ck('all24_unique_surface_packets',surfaces=={(a,b)for a in(0,1)for b in range(12)})
ck('all24_active_bulk_packets',active==set(p['bulk_active_elements'])and len(active)==24)
R=[[Q.from_float(float(v))for v in row]for row in target.a['preconditioner']]
def exponent(q):
    den=q.denominator
    if den&(den-1):raise ValueError('exact dyadic endpoint required')
    return den.bit_length()-1
kr=max(exponent(v)for row in R for v in row)
RI=[[v.numerator<<(kr-exponent(v))for v in row]for row in R]
partials=[];norms=[]
for j,(key,values)in enumerate(zip(order,columns)):
    produced=d['parameter_partials'][j]
    ck('parameter_order_'+key,produced['parameter']==key)
    ck('all611_force_sums_'+key,all(included(v,interval(w))for v,w in zip(values,produced['force_partial'])))
    kx=max(exponent(v)for pair in values for v in pair)
    XI=[(a.numerator<<(kx-exponent(a)),b.numerator<<(kx-exponent(b)))for a,b in values]
    denom=1<<(kr+kx);result=[]
    for row in RI:
        lo=sum(c*(v[0]if c>=0 else v[1])for c,v in zip(row,XI))
        hi=sum(c*(v[1]if c>=0 else v[0])for c,v in zip(row,XI))
        result.append((Q(lo,denom),Q(hi,denom)))
    ck('all611_signed_R_products_'+key,all(included(v,interval(w))for v,w in zip(result,produced['preconditioned_force_partial'])))
    norm=max(max(abs(a),abs(b))for a,b in result);norms.append(norm)
    ck('producer_norm_encloses_exact_reduce_'+key,norm<=Q(produced['infinity_norm_upper']['upper_exact']))
    partials.append({'parameter':key,'exact_force_intervals':[[str(a),str(b)]for a,b in values],
                     'exact_preconditioned_intervals':[[str(a),str(b)]for a,b in result],
                     'norm_upper_exact':str(norm),'R_common_denominator_exponent':kr,'force_common_denominator_exponent':kx})
delta=[max(abs(Q(a)),abs(Q(b)))for a,b in p['parameter_deltas_exact']]
parameter_bound=sum(v*r for v,r in zip(norms,delta))
ck('parameter_bound_enclosed_by_producer',parameter_bound<=Q(d['parameter_only_variation_upper_exact']))
ck('no_V_or_root_claim',d['V_established']is False and d['exact_integral_root_certified']is False)
out={'status':'PASS_EXACT_PARAMETER_COMPONENT_REDUCTION','checks':checks,'producer_results_sha256':sha(HERE/'results.json'),
     'protocol_sha256':sha(HERE/'protocol.json'),'review_source_sha256':sha(Path(__file__)),
     'parameter_partials':partials,'parameter_only_variation_upper_exact':str(parameter_bound),
     'predictor_displacement_norm_upper_exact':p['predictor_displacement_norm_upper_exact'],
     'local_nonlinear_or_polynomial_reintegrations_by_reviewer':False,'V_established':False,'root_certified':False,
     'seconds':time.monotonic()-start,
     'scope':'Independent exact local algebra, all48 signed local/global sums and integer-dyadic R products, bound inputs and exact parameter-weighted norm reduction. Reviewed local polynomial moments/cache maps reused; no full-b/V/root/PDE assertion.'}
dest=ROOT/'evidence/exact_parameter_defect_B_independent_reduction.json';assert not dest.exists()
dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps({'status':out['status'],'checks':len(checks),'parameter_bound':float(parameter_bound),'seconds':out['seconds'],'review_sha256':sha(dest)}),flush=True)
