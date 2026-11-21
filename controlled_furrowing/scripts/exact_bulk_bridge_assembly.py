"""Whole B-center bulk integration discrepancy using the reviewed local kernel.

No surface/full-tube/root claim. Orders depend on rigorous POINT chart bounds.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,argparse,hashlib,datetime,time
from flint import arb,arb_mat,ctx
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import exact_bulk_bridge_pilot as P
T=P.T;M=P.M;JET=P.JET;core=P.core
OUT=ROOT/'runs/exact_integration_bridge_v1/bulk_assembly_B'
SOURCE=Path(__file__)
PROTOCOL=OUT/'protocol_v2.json'
PREFLIGHT=OUT/'preflight_v2.json'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def write(p,d):P.write(p,d)
def pack(v):return core.pack(v)
def local_rho(fields,poly):
    """An exact rational bound on |J-1| on this entire POINT element."""
    if fields['chart']['kind']=='off_axis':
        rep=JET.shifted_numerator(poly,fields['corners'])
        z=max(abs(v[1]-fields['zcenter']) for v in fields['corners'])
        bound=sum((abs(v)*z**j for (i,j),v in rep['numerator'].items()),F(0))/rep['Rlo']**rep['denominator_power']
    else:
        shifted=T.shift(poly,F(1,2),F(1,2))
        bound=sum((abs(v)*F(1,2)**(i+j) for (i,j),v in shifted.items()),F(0))
    return min(F(21,500),bound),bound
def select(target,index,rule):
    f=target.element(index);jet=JET.bulk_jacobian_jet(f['q'],f['B']);rho,bound=local_rho(f,T.sub(jet['value'],T.ONE))
    fun=P.Functional(target,f,index,rule);mass=F(str(core.upper(fun.mass_integral+fun.mass_rule).fmpq()))
    ju=[P.pnorm(p,f) for p in jet['gradient']]
    maxju=max(ju);maxj2=max(P.pnorm(p,f) for p in jet['hessian_upper'].values())
    for n in range(4,25):
        t1,t2=M.volumetric_derivative_tail_bounds(n,rho,1,5000)
        fbound=mass*t1*maxju;hbound=mass*(t2*maxju**2+t1*maxj2)
        if fbound<=F(1,10**22) and hbound<=F(1,10**17):
            return {'element':index,'chart':f['chart']['kind'],'point_rho_exact':str(rho),'independent_coefficient_majorant_exact':str(bound),'degree':n,'force_tail_budget_upper_exact':str(fbound),'H_entry_tail_budget_upper_exact':str(hbound),'global_indices':f['global_indices']}
    raise ValueError('No allowed degree meets declared POINT tail budget for element '+str(index))
def prepare():
    ctx.prec=512
    pilot=read(P.OUT/'protocol_v2.json');review=read(ROOT/'evidence/exact_bulk_pilot_prelaunch_independent_review.json')
    if not review['status'].startswith('PASS') or not review['checks'] or not all(review['checks'].values()):raise ValueError('reviewed local kernel required')
    bindings={}
    for collection in ('source_and_input_hashes','target_preparation_input_hashes'):
        for p,h in pilot[collection].items():
            if sha(ROOT/p)!=h:raise ValueError('changed reviewed input '+p)
            bindings[p]=h
    for p in (SOURCE,P.OUT/'protocol_v2.json',ROOT/'evidence/exact_bulk_pilot_prelaunch_independent_review.json'):
        bindings[str(p.relative_to(ROOT))]=sha(p)
    record=read(ROOT/P.RECORD);target=T.Target(read(ROOT/record['point_cache_file']));rule=target.auxiliary_rule()
    guards={
        'axis_centered_majorant':local_rho({'chart':{'kind':'axis_vertex'}},{(1,0):F(1,100)})==(F(1,100),F(1,100)),
        'offaxis_inverse_majorant':local_rho({'chart':{'kind':'off_axis'},'corners':[(F(1),F(0)),(F(2),F(0)),(F(1),F(1))],'zcenter':F(0)},{(-1,0):F(1,100)})==(F(1,100),F(1,100)),
        'global_rho_fallback':local_rho({'chart':{'kind':'axis_vertex'}},{(0,0):F(1,10)})==(F(21,500),F(1,10))}
    if not all(guards.values()):raise AssertionError(guards)
    selectors=[select(target,i,rule) for i in range(144)]
    protocol={'utc':P.now(),'source_and_input_hashes':bindings,'record_id':record['record_id'],'precision_bits':512,'selectors':selectors,'force_tail_target_exact':'1/10000000000000000000000','H_entry_tail_target_exact':'1/100000000000000000','policy':'Smallest degree4..24 meeting exact POINT majorant budgets; rho=min(reviewed global21/500,independent exact shifted-polynomial coefficient bound).','scope':'ALL144 BULK elements at ONE saved center. Surface, tube uncertainty and equilibrium transfer excluded.','source_data_fit':False,'root_transfer_proved':False}
    write(PROTOCOL,protocol);write(PREFLIGHT,{'status':'PASS_PREPARATION_ONLY','guards':guards,'all144_orders_selected':True,'order_histogram':{str(n):sum(s['degree']==n for s in selectors) for n in sorted({s['degree'] for s in selectors})},'protocol_sha256':sha(PROTOCOL),'source_sha256':sha(SOURCE),'actual_integrals_run':False})
    print(json.dumps({'status':'PASS_PREPARATION_ONLY','source_sha256':sha(SOURCE),'protocol_sha256':sha(PROTOCOL),'degrees':sorted({s['degree'] for s in selectors})}),flush=True)
def run():
    p=read(PROTOCOL);ctx.prec=p['precision_bits']
    review=read(OUT/'prelaunch_independent_review.json')
    if review.get('status')!='PASS' or not isinstance(review.get('checks'),dict) or not review['checks'] or not all(v is True for v in review['checks'].values()) or review['source_sha256']!=sha(SOURCE) or review['protocol_sha256']!=sha(PROTOCOL) or review['preflight_sha256']!=sha(PREFLIGHT):raise ValueError('bound nonempty independent launch review required')
    for path,h in p['source_and_input_hashes'].items():
        if sha(ROOT/path)!=h:raise ValueError('changed frozen input '+path)
    record=read(ROOT/P.RECORD);target=T.Target(read(ROOT/record['point_cache_file']));rule=target.auxiliary_rule();started=time.monotonic()
    write(OUT/'started.json',{'utc':P.now(),'protocol_sha256':sha(PROTOCOL),'source_sha256':sha(SOURCE)})
    g=[arb(0) for _ in range(611)];H={};files=[]
    for selector in p['selectors']:
        i=selector['element'];drive=F(record['base_strength_exact'])*F.from_float(float(target.a['active'][i,0]))
        result=P.run_element(target,i,rule,selector['degree'],F(selector['point_rho_exact']),drive)
        indices=result['global_indices']
        for k,value in zip(indices,result['force_integral_minus_cache']):g[k]+=core.A(F(value['lower_exact'])).union(core.A(F(value['upper_exact'])))
        for item in result['H_integral_minus_cache_upper_triangle']:
            key=tuple(sorted((indices[item['i']],indices[item['j']])));v=item['enclosure'];H[key]=H.get(key,arb(0))+core.A(F(v['lower_exact'])).union(core.A(F(v['upper_exact'])))
        dest=OUT/'elements'/f'{i:03d}.json';write(dest,result);files.append({'element':i,'file':str(dest.relative_to(ROOT)),'sha256':sha(dest)})
        print(json.dumps({'phase':'element_complete','element':i,'degree':selector['degree'],'seconds':result['seconds'],'completed':len(files)}),flush=True)
    R=arb_mat(611,611,[core.A(float(v)) for v in target.a['preconditioner'].ravel()]);RF=R*arb_mat(611,1,g)
    rows=[arb(0) for _ in g]
    for (i,j),v in H.items():
        rows[i]+=core.upper(v)
        if i!=j:rows[j]+=core.upper(v)
    rnorm=max(core.upper(sum((core.upper(R[i,j]) for j in range(611)),arb(0))) for i in range(611))
    hnorm=max(core.upper(v) for v in rows)
    for path,h in p['source_and_input_hashes'].items():
        if sha(ROOT/path)!=h:raise ValueError('input changed during assembly '+path)
    result={'status':'ALL144_POINT_BULK_ASSEMBLED_UNREVIEWED','utc':P.now(),'protocol_sha256':sha(PROTOCOL),'elements':files,'free_dofs':611,'H_upper_slots':len(H),'force_integral_minus_cache':[pack(v) for v in g],'preconditioned_bulk_force_defect':[pack(RF[i,0]) for i in range(611)],'bulk_force_norm_upper':pack(max(core.upper(v) for v in g)),'preconditioned_bulk_force_norm_upper':pack(max(core.upper(RF[i,0]) for i in range(611))),'bulk_H_row_norm_upper':pack(hnorm),'conservative_preconditioned_bulk_H_norm_upper':pack(rnorm*hnorm),'bulk_H_upper':[{'i':i,'j':j,'enclosure':pack(v)} for (i,j),v in sorted(H.items())],'seconds':time.monotonic()-started,'all_frozen_bindings_unchanged':True,'root_transfer_proved':False,'scope':p['scope']}
    write(OUT/'results.json',result);print(json.dumps({'status':result['status'],'seconds':result['seconds'],'RF_bulk_upper':result['preconditioned_bulk_force_norm_upper']['ball']}),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--prepare',action='store_true');ap.add_argument('--run',action='store_true');args=ap.parse_args()
    if args.prepare==args.run:raise ValueError('choose prepare OR run')
    prepare() if args.prepare else run()
