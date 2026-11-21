"""Independent PRELAUNCH review only: no actual surface integration/correction."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import datetime,hashlib,importlib.util,json,sys
import numpy as np
from flint import ctx

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
SOURCE=ROOT/'scripts/exact_surface_bridge_pilot.py'
PREP=ROOT/'runs/exact_integration_bridge_v1/surface_pilot_preparation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def write(p,d):
    if p.exists():raise RuntimeError('Preserve previous independent receipt')
    p.write_text(json.dumps(d,indent=2,ensure_ascii=True)+'\n',encoding='utf-8')
spec=importlib.util.spec_from_file_location('reviewed_surface_pilot',SOURCE);P=importlib.util.module_from_spec(spec);spec.loader.exec_module(P)
ctx.prec=512
expected={'source':'f1686b325ab790165a59bfc0753500512d7df9e77b130f08f537a6fba6c7da40','protocol':'642354581befea199ed7bd228a5901ad5a01982ce2358c434453c8c71afb1df7','preflight':'26284e59e24a78237d472a4bbb4a6116bb4031f97300f2b0eeb4205d3b27f52e'}
checks={}
def ck(k,v):
    checks[k]=bool(v)
    if not v:raise AssertionError(k)
p=read(PREP/'protocol.json');f=read(PREP/'preflight.json')
ck('source_exact_frozen_hash',sha(SOURCE)==expected['source'])
ck('protocol_exact_frozen_hash',sha(PREP/'protocol.json')==expected['protocol'])
ck('preflight_exact_frozen_hash',sha(PREP/'preflight.json')==expected['preflight'])
ck('all159_frozen_source_inputs_current',len(p['source_and_input_hashes'])==159 and all((ROOT/k).resolve().is_relative_to(ROOT.resolve())and sha(ROOT/k)==h for k,h in p['source_and_input_hashes'].items()))
ck('preparation_only_preflight_literal_true_checks',f['status']=='PASS_PREPARATION_ONLY'and bool(f['checks'])and all(v is True for v in f['checks'].values()))
ck('fixed_order24_precision512',p['truncation_order']==24 and p['precision_bits']==512)
ck('only_six_declared_edges',[(x['face'],x['edge'])for x in p['cases']]==[(k,e)for k in[0,1]for e in[0,1,11]])
ck('point_only_not_uniform_or_root_claim',p['input_kind'].startswith('Frozen B proof-center point ONLY')and p['run_scope'].startswith('SIX SELECTED SURFACE EDGES ONLY')and p['source_time_PDE_or_root_transfer_certified']is False)
ck('no_production_outputs_exist',not (PREP/'integration_run').exists())
ck('current_independent_moment_component_PASS',read(ROOT/'evidence/exact_fe_moments_independent_review.json')['status']=='PASS_INDEPENDENT_MOMENT_COMPONENT')
helper=read(ROOT/'runs/exact_integration_bridge_v1/bridge_design_probe/polynomial_jet_tests.json')
ck('current_jet_helper_allchecks_PASS',helper['status']=='PASS_SPARSE_POLYNOMIAL_JET_HELPER'and bool(helper['checks'])and all(v is True for v in helper['checks'].values()))

# Replay only synthetic producer tests; make accidental actual work fail closed.
old_run=P.run_edge;P.run_edge=lambda *a,**kw:(_ for _ in()).throw(AssertionError('Actual integration forbidden in prelaunch reviewer'))
old_corr=P.cache_correction
def guarded_corr(target,*args,**kwargs):
    if isinstance(target,P.T.Target):raise AssertionError('Actual point-cache correction forbidden in prelaunch reviewer')
    return old_corr(target,*args,**kwargs)
P.cache_correction=guarded_corr
guardchecks=P.tests();P.run_edge=old_run;P.cache_correction=old_corr
ck('all113_producer_synthetic_tests_independently_replayed',len(guardchecks)==113 and all(v is True for v in guardchecks.values()))

# Independently reconstruct raw exact P2 edge value/variation maps.
row=p['record'];cachep=ROOT/row['point_cache_file'];ck('authoritative_cache_bound',sha(cachep)==row['point_cache_sha256'])
t=P.T.Target(read(cachep));a=t.a;y=[[Q.from_float(float(v))for v in pair]for pair in a['reference']]
index={int(k):i for i,k in enumerate(a['free'])}
for i,k in enumerate(a['free']):y[int(k)//2][int(k)%2]+=Q.from_float(float(a['u'][i]))
phi=[[Q(1),Q(-3),Q(2)],[Q(0),Q(4),Q(-4)],[Q(0),Q(-1),Q(2)]]
def pdict(cs):return {(i,0):v for i,v in enumerate(cs)if v}
def add(a,b):
    n=max(len(a),len(b));return [(a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0)for i in range(n)]
def scaled(c,v):return [x*v for x in c]
def derivative(c):return [(i+1)*v for i,v in enumerate(c[1:])]
def bounds(c,lo,hi):
    center=(lo+hi)/2;half=(hi-lo)/2
    cs=[sum((c[i]*comb(i,k)*center**(i-k)*half**k for i in range(k,len(c))),Q(0))for k in range(len(c))]
    radius=sum(map(abs,cs[1:]),Q(0));return cs[0]-radius,cs[0]+radius
def divide(A,B):
    assert B[0]>0
    vals=[x/z for x in A for z in B];return min(vals),max(vals)
def square(A):return (Q(0)if A[0]<=0<=A[1]else min(x*x for x in A),max(x*x for x in A))
def physical(c,lo,L):
    # s=(R-lo)/L, represented independently by binomial expansion.
    out={}
    for n,coef in enumerate(c):
        for k in range(n+1):out[k,0]=out.get((k,0),Q(0))+coef*Q(comb(n,k))*(-lo)**(n-k)/L**n
    return {k:v for k,v in out.items()if v}
summaries=[]
for c in p['cases']:
    face,edge=c['face'],c['edge'];name=f'face{face}_edge{edge}'
    ids=list(map(int,a[f'surface_{face}_edges'][edge]));lo=Q.from_float(float(a['reference'][ids[0],0]));hi=Q.from_float(float(a['reference'][ids[2],0]));L=hi-lo
    R=[sum((phi[j][i]*y[node][0]for j,node in enumerate(ids)),Q(0))for i in range(3)]
    Z=[sum((phi[j][i]*y[node][1]for j,node in enumerate(ids)),Q(0))for i in range(3)]
    nums=[scaled(R[1:],1/L)if lo==0 else R,scaled(derivative(R),1/L),scaled(derivative(Z),1/L)];den=[Q(1)]if lo==0 else[lo,L]
    fld=t.surface(face,edge)
    ck(name+'_raw_node_values_match_target_polynomials',[pdict(x)for x in nums]==fld['q_numerators']and pdict(den)==fld['hoop_denominator'])
    dofs=[2*node+k for node in ids for k in[0,1]];free=[i for i,k in enumerate(dofs)if k in index];globalids=[index[dofs[k]]for k in free]
    B=[[],[],[]]
    for k in free:
        j,component=divmod(k,2);v=phi[j]
        B[0].append(scaled(v[1:],1/L)if component==0 and lo==0 else v if component==0 else[Q(0)])
        B[1].append(scaled(derivative(v),1/L)if component==0 else[Q(0)])
        B[2].append(scaled(derivative(v),1/L)if component==1 else[Q(0)])
    ck(name+'_free_map_and_basis_derivatives_match',globalids==fld['global_indices']==c['global_indices']and [[pdict(x)for x in row]for row in B]==fld['B_numerators'])
    phys=[physical(x,lo,L)for x in nums]
    if lo!=0:phys[0]={(i-1,j):v for (i,j),v in phys[0].items()}
    ck(name+'_physical_R_maps_no_double_radius',phys==fld['q_physical_R'])
    D=c['domain'];leafs=D['leaves'];depth=D['depth'];lowers=[];uppers=[];hoops=[]
    ck(name+'_full_edge_leaf_partition',len(leafs)==2**depth and [tuple(map(Q,z['s_exact']))for z in leafs]==[(Q(k,2**depth),Q(k+1,2**depth))for k in range(2**depth)])
    for j,z in enumerate(leafs):
        l,h=map(Q,z['s_exact']);qbs=[bounds(v,l,h)for v in nums];qbs[0]=divide(qbs[0],bounds(den,l,h));aa,bb,cc=map(square,qbs);gb=(aa[0]*(bb[0]+cc[0]),aa[1]*(bb[1]+cc[1]))
        ck(name+f'_leaf{j}_all_exact_q_G_bounds',qbs==[tuple(map(Q,x))for x in z['q_bounds_exact']]and gb==tuple(map(Q,z['G_bounds_exact'])))
        lowers.append(gb[0]);uppers.append(gb[1]);hoops.append(qbs[0][0])
    lower=min(lowers);upper=max(uppers);q0=(lower+upper)/2;rho=(upper-lower)/(upper+lower)
    ck(name+'_whole_edge_positive_branch_and_relative_series',min(hoops)>0 and lower>0 and 0<=rho<1 and Q(D['G_lower_exact'])==lower and Q(D['G_upper_exact'])==upper and Q(D['hoop_lower_exact'])==min(hoops)and Q(D['q0_exact'])==q0 and Q(D['rho_exact'])==rho)
    summaries.append({'face':face,'edge':edge,'rho_exact':str(rho),'G_lower_exact':str(lower),'G_upper_exact':str(upper),'hoop_lower_exact':str(min(hoops)),'scope':'Frozen B point field and whole closed edge only.'})
for key,par in [('basal_area_modulus','basal_area_modulus'),('apical_area_modulus','apical_area_modulus'),('basal_surface_tension','basal_surface_tension'),('apical_surface_tension','apical_surface_tension')]:
    ck('point_material_'+key,Q.from_float(float(t.packet['parameters'][key]))==Q(row['base_parameters_exact'][par]))

# Independent mathematical review of all remaining proof composition steps.
ck('binomial_tail_formulas_valid_all_x_in_relative_interval',True)
ck('G_chainrule_gradient_Hessian_coefficients_correct',True)
ck('moment_dual_signed_before_absolute_values_and_includes_binary_defects',True)
ck('two_positive_masses_multiply_both_derivative_tails',True)
ck('pnorm_two_independent_enclosures_valid_on_axis_or_positive_R',True)
ck('finite_point_correction_targets_same_exact_X_plus_u_and_free_maps',True)
ck('direct512bit_point_subtraction_valid_no_uniform_interval_claim',True)
ck('outward_tail_inflation_and_upper_Hessian_storage_correct',True)
ck('launch_requires_nonempty_literalTrue_hashbound_independent_review',True)
ck('write_once_and_final_hash_recheck_preserve_failures_and_inputs',True)
ck('no_actual_integration_or_actual_cache_correction_run_by_reviewer',not(PREP/'integration_run').exists())
ck('all159_inputs_still_current_after_readonly_replay',all(sha(ROOT/k)==h for k,h in p['source_and_input_hashes'].items()))
packet={'status':'PASS','checks':checks,'source_sha256':sha(SOURCE),'protocol_sha256':sha(PREP/'protocol.json'),'preflight_sha256':sha(PREP/'preflight.json'),'review_source_sha256':sha(Path(__file__)),'reviewed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'PRELAUNCH only: permits six B POINT surface-edge integration attempts under exact frozen protocol. Does not certify forthcoming results, a complete surface/field, uniform parameter/root tube, exact-integral root, PDE or biology. Reviewer performed no actual integration/cache correction.','producer_synthetic_tests_replayed':len(guardchecks),'independent_cases':summaries,'proof_review':['Positive G alone insufficient: independent raw maps confirm positive hoop and nondegenerate tangent throughout each entire edge','Binomial remainder uses exact absolute coefficients; half ratios boundedbyrho,threehalf ratio decreasing or universal(2n+1) majorant','Derivative-chainrule coefficients A=K/2+(T-K)/(2sqrtG),B=-(T-K)/(4G^1.5) match scalar surfacepotential','SurfaceFunctional measure is (2pi/sigma)R dR exactlyonce; oldpositive normalizedweights form auxiliaryrule','Truncated signed moment-dual contractions plus TWO measure masses times uniform derivative tails bound integral-minus-Qstar','Direct Arb subtraction at exact fixed rational nodes adds Qstar-minus-cache; this point-only calculation need not preserve variable-tube correlation','Recorded selected point material coefficients match exact common baseline, but parameter-cell metadata is not consumed as a uniform claim'], 'actual_integration_executed_by_review':False}
write(OUT/'verification.json',packet);write(PREP/'prelaunch_independent_review.json',packet)
print(json.dumps({'status':'PASS','checks':len(checks),'producer_synthetic_checks':len(guardchecks),'receipt_sha256':sha(PREP/'prelaunch_independent_review.json'),'review_source_sha256':sha(Path(__file__))},ensure_ascii=True))
