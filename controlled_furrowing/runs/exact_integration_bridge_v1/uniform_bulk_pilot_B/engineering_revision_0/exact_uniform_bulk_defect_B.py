"""Preparation-first uniform bulk discrepancy on the NEW common B domain.
No root is assumed. Pilot covers seven predeclared bulk elements, not all144.
"""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,argparse,time
from functools import lru_cache
from flint import arb,arb_mat,ctx
import exact_bulk_bridge_pilot as P
import exact_cache_consistency as C
T=P.T;M=P.M;JET=P.JET;core=P.core
ROOT=Path(__file__).resolve().parents[1]
SOURCE=Path(__file__)
OUT=ROOT/'runs/exact_integration_bridge_v1/uniform_bulk_pilot_B'
GEOM=ROOT/'runs/exact_integration_bridge_v1/retube_B_geometry'
REVIEW=ROOT/'evidence/exact_retube_B_geometry_independent_review.json'
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):P.write(p,d)
def A(v):return v if isinstance(v,arb) else core.A(F(v))
def coeff_abs(v):return F(str(core.upper(v).fmpq())) if isinstance(v,arb) else abs(F(v))
def pnorm(poly,fields):
    if fields['chart']['kind']=='off_axis':
        lo=min(v[0] for v in fields['corners']);hi=max(v[0] for v in fields['corners']);z=max(abs(v[1]-fields['zcenter']) for v in fields['corners'])
        return sum((coeff_abs(v)*(hi**i if i>=0 else lo**i)*z**j for (i,j),v in poly.items()),F(0))
    return sum((coeff_abs(v) for v in poly.values()),F(0))
def series(poly,n):
    c1,c2=M.volumetric_derivative_coefficients(n,1,5000);power=T.ONE;AA={};BB={};poly={k:A(v) for k,v in poly.items()}
    for i in range(n+1):
        AA=T.add(AA,T.scale(power,A(c1[i])));BB=T.add(BB,T.scale(power,A(c2[i])))
        if i<n:power=T.mul(power,poly)
    return AA,BB
class RuleTarget:
    """Identical frozen target with memoized auxiliary-rule construction only."""
    def __init__(self,target,rule):self.target=target;self.rule=rule
    def __getattr__(self,k):return getattr(self.target,k)
    def auxiliary_rule(self):return self.rule
def field_nodes(box):
    return [[A(F(c)-F(r)).union(A(F(c)+F(r))) for c,r in zip(cr,rr)] for cr,rr in zip(box['current_centers_exact'],box['coordinate_radii_exact'])]
def gain_drive(protocol,active):
    dose=F(protocol['dose_exact']);lo,hi=map(F,protocol['parameter_cell']['common_gain'])
    return A(F(0)) if active==0 else A(dose*lo).union(A(dose*hi))
def uniform_cache(target,fields,element,Y,drive):
    """Signed mean-value force and third-tensor Hessian discrepancies."""
    n=len(fields['global_indices']);Fout=[arb(0) for _ in range(n)];Hout=arb_mat(n,n);maxdq=arb(0);minJ=None
    for iq in range(64):
        maps=C.exact_bulk_maps(target,element,iq)
        if maps['global_indices']!=fields['global_indices'] or maps['selected_local']!=fields['free_local']:raise ValueError('identical free projection required')
        nodes=[Y[j][k] for j in maps['nodes'] for k in (0,1)]
        result=C.consistency_bound(maps['star_full'],maps['cache_full'],nodes,maps['selected_local'],1,5000,drive)
        _,Hline,_=core.bulk_local(result['convex_line_box'],A(1),A(5000),drive)
        gc,_,_=core.bulk_local(result['cache_q'],A(1),A(5000),drive)
        delta=arb_mat(5,n,[A(row[i]) for row in result['delta_B_full'] for i in maps['selected_local']])
        dg=arb_mat(5,5,[v for row in Hline for v in row])*arb_mat(5,1,result['delta_q'])
        df=delta.transpose()*arb_mat(5,1,gc)+result['star_B_free'].transpose()*dg
        w=A(maps['auxiliary_normalized_weight_exact'])
        for i in range(n):
            Fout[i]+=w*df[i,0]
            for j in range(n):Hout[i,j]+=w*result['difference'][i,j]
        maxdq=max(maxdq,max(core.upper(v) for v in result['delta_q']))
        lower=result['domain']['J'].lower();minJ=lower if minJ is None else min(minJ,lower)
    return Fout,Hout,{'max_coefficientwise_delta_q_upper':core.pack(maxdq),'all64_convex_strain_lines_positive':True,'minimum_cache_convex_line_J_lower':core.pack(minJ)}
def choose(target,fields,index,rule,rho):
    jet=JET.bulk_jacobian_jet(fields['q'],fields['B']);fun=P.Functional(target,fields,index,rule)
    mass=F(str(core.upper(fun.mass_integral+fun.mass_rule).fmpq()));g=max(pnorm(v,fields) for v in jet['gradient']);h=max(pnorm(v,fields) for v in jet['hessian_upper'].values())
    for n in range(2,19):
        t1,t2=M.volumetric_derivative_tail_bounds(n,rho,1,5000);fb=mass*t1*g;hb=mass*(t2*g*g+t1*h)
        if fb<=F(1,10**22) and hb<=F(1,10**16):return {'element':index,'rho_exact':str(rho),'degree':n,'force_tail_upper_exact':str(fb),'H_entry_tail_upper_exact':str(hb),'active':str(fields['active'])}
    raise ValueError('declared uniform pilot order budget unresolved')
def evaluate(target,index,Y,rule,selector,protocol):
    start=time.monotonic();fields=target.element(index,Y);jet=JET.bulk_jacobian_jet(fields['q'],fields['B']);fun=P.Functional(target,fields,index,rule);rho=F(selector['rho_exact']);n=selector['degree'];drive=gain_drive(protocol,fields['active'])
    AA,BB=series(T.sub(jet['value'],T.ONE),n);poly=JET.contract_volumetric_defect(jet,AA,BB,fun.defect)
    mass=core.upper(fun.mass_integral+fun.mass_rule);t1,t2=M.volumetric_derivative_tail_bounds(n,rho,1,5000);norms=[pnorm(v,fields) for v in jet['gradient']]
    diag=[A(1),A(1)+drive,A(1),A(1)+drive,A(1)];d=len(fields['global_indices']);fg=[];HH={};ft=[];ht={}
    for i in range(d):
        shear={}
        for k in range(5):shear=T.add(shear,T.scale(T.mul(fields['q'][k],fields['B'][k][i]),diag[k]))
        fg.append(poly['gradient'][i]+fun.apply(shear));ft.append(mass*A(t1*norms[i]))
    for i,j in jet['hessian_upper']:
        shear={}
        for k in range(5):shear=T.add(shear,T.scale(T.mul(fields['B'][k][i],fields['B'][k][j]),diag[k]))
        HH[i,j]=poly['hessian_upper'][i,j]+fun.apply(shear);ht[i,j]=mass*A(t2*norms[i]*norms[j]+t1*pnorm(jet['hessian_upper'][i,j],fields))
    cf,ch,cache_info=uniform_cache(RuleTarget(target,rule),fields,index,Y,drive)
    inflate=lambda v,t:v+(-core.upper(t)).union(core.upper(t))
    FF=[inflate(v,ft[i])+cf[i] for i,v in enumerate(fg)];HH={k:inflate(v,ht[k])+ch[k[0],k[1]] for k,v in HH.items()}
    # A known exact base-point field belongs to this common geometry/parameter
    # box. Independent local point packets must overlap every uniform entry.
    point=read(ROOT/f'runs/exact_integration_bridge_v1/bulk_assembly_B/elements/{index:03d}.json')
    overlap=lambda a,b:not(F(a['upper_exact'])<F(b['lower_exact']) or F(b['upper_exact'])<F(a['lower_exact']))
    Fpacked=[core.pack(v) for v in FF];Hpacked=[{'i':i,'j':j,'enclosure':core.pack(v),'two_mass_tail_upper':core.pack(ht[i,j])} for (i,j),v in HH.items()]
    pH={(v['i'],v['j']):v['enclosure'] for v in point['H_integral_minus_cache_upper_triangle']}
    checks={'base_point_free_map_matches':point['global_indices']==fields['global_indices'],'all_base_point_force_entries_overlap':all(overlap(a,b) for a,b in zip(Fpacked,point['force_integral_minus_cache'])),'all_base_point_H_entries_overlap':all(overlap(v['enclosure'],pH[v['i'],v['j']]) for v in Hpacked),'all64_cache_convex_lines_positive':cache_info['all64_convex_strain_lines_positive']}
    if not all(checks.values()):raise AssertionError(checks)
    return {'status':'UNIFORM_COMMON_BOX_BULK_ELEMENT_UNREVIEWED','element':index,'active':str(fields['active']),'degree':n,'rho_exact':str(rho),'global_indices':fields['global_indices'],'force_integral_minus_cache':Fpacked,'H_integral_minus_cache_upper':Hpacked,'force_two_mass_tail_upper':[core.pack(v) for v in ft],'cache_consistency':cache_info,'checks':checks,'seconds':time.monotonic()-start,'scope':'Every canonical field in NEW common nodal box and every shared gain in SAME cell; ONE bulk element. Not assembled611force/H, surface discrepancy or root certificate.'}
def tests():
    toy={'chart':{'kind':'off_axis'},'corners':[(F(1),F(0)),(F(2),F(0)),(F(1),F(1))],'zcenter':F(0)}
    zero=arb(0);small=A(F(-1,3)).union(A(F(1,2)))
    checks={'interval_coefficient_majorant':pnorm({(-2,1):small},toy)>=F(1,2),'zero_coeff_retained_correctly':pnorm({(0,0):zero},toy)==0}
    aa,bb=series({(1,0):A(F(1,10))},2)
    checks['uniform_series_exact_point_fixture']=aa[(1,0)].overlaps(A(F(5001,10))) and bb[(2,0)].overlaps(A(F(27503,100)))
    # Mean-value gradient difference formula is a direct consequence of the
    # fundamental theorem on the same analytic segment; matched maps vanish.
    eye=[[F(int(i==j)) for j in range(5)] for i in range(5)];Y=[A(1),A(F(1,20)),A(0),A(1),A(1)]
    c=C.consistency_bound(eye,eye,Y,list(range(5)),1,5000,A(1).union(A(2)))
    checks['matched_maps_zero_interval_H_difference']=all(v.is_zero() for v in c['difference'].entries()) and all(v.is_zero() for v in c['delta_q'])
    if not all(checks.values()):raise AssertionError(checks)
    return checks
def prepare():
    ctx.prec=512;core.A.cache_clear();checks=tests();geomreview=read(REVIEW)
    if not geomreview['status'].startswith('PASS') or not all(geomreview['checks'].values()):raise ValueError('independent NEW common geometry required')
    geom=read(GEOM/'certification_run/results.json');body=read(GEOM/'certification_run/whole_body_geometry.json');box=read(GEOM/'coordinate_boxes.json');record=read(ROOT/P.RECORD);target=T.Target(read(ROOT/record['point_cache_file']));rule=target.auxiliary_rule();Y=field_nodes(box)
    bind={}
    for p in [SOURCE,Path(P.__file__),Path(C.__file__),Path(T.__file__),Path(M.__file__),Path(JET.__file__),Path(core.__file__),REVIEW,GEOM/'certification_run/results.json',GEOM/'certification_run/whole_body_geometry.json',GEOM/'coordinate_boxes.json',ROOT/P.RECORD,ROOT/record['point_cache_file'],ROOT/'evidence/exact_cache_consistency_independent_review.json']:
        bind[p.relative_to(ROOT).as_posix()]=sha(p)
    # Transitive geometry protocol binds the candidate, factory and all original
    # physical arrays; adding its hashes makes later source drift explicit.
    gp=read(GEOM/'protocol.json')
    for path,h in gp['source_and_input_hashes'].items():
        if sha(ROOT/path)!=h:raise ValueError('geometry source/input changed '+path)
        if path in bind and bind[path]!=h:raise ValueError('conflicting binding '+path)
        bind[path]=h
    selectors=[]
    for i in (0,2,6,8,72,74,66):
        jb=body['triangles'][i]['J_bound'];rho=max(abs(F(jb['lower_exact'])-1),abs(F(jb['upper_exact'])-1))
        if not 0<=rho<1 or F(jb['lower_exact'])<=F(1,10**6):raise ValueError('fresh element analytic domain required')
        selectors.append(choose(target,target.element(i,Y),i,rule,rho))
        pp=ROOT/f'runs/exact_integration_bridge_v1/bulk_assembly_B/elements/{i:03d}.json';bind[pp.relative_to(ROOT).as_posix()]=sha(pp)
    p={'utc':P.now(),'source_and_input_hashes':bind,'geometry_domain_id':box['domain_id'],'geometry_boxes_sha256':sha(GEOM/'coordinate_boxes.json'),'parameter_cell':record['parameter_cell'],'parameter_cell_sha256':record['parameter_cell_sha256'],'dose_exact':record['dose_exact'],'precision_bits':512,'selectors':selectors,'force_tail_target_exact':'1/10000000000000000000000','H_entry_tail_target_exact':'1/10000000000000000','scope':'SEVEN bulk elements on NEW common nodal box and SAME parameter cell; both active and inactive cases. No611field/surface/root certificate.','root_certified':False}
    write(OUT/'protocol.json',p);write(OUT/'preflight.json',{'status':'PASS_UNIFORM_BULK_PREPARATION_ONLY','checks':checks,'source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json'),'selected_elements':[v['element'] for v in selectors],'selected_degrees':[v['degree'] for v in selectors],'actual_integrals_or_64site_cache_sums_run':False})
    print(json.dumps({'status':'PASS_UNIFORM_BULK_PREPARATION_ONLY','source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json'),'degrees':[v['degree'] for v in selectors]}),flush=True)
def run():
    p=read(OUT/'protocol.json');ctx.prec=p['precision_bits'];core.A.cache_clear()
    review=read(OUT/'prelaunch_independent_review.json')
    if review.get('status')!='PASS' or not isinstance(review.get('checks'),dict) or not review['checks'] or not all(v is True for v in review['checks'].values()) or review['source_sha256']!=sha(SOURCE) or review['protocol_sha256']!=sha(OUT/'protocol.json') or review['preflight_sha256']!=sha(OUT/'preflight.json'):raise ValueError('bound nonempty review required')
    for path,h in p['source_and_input_hashes'].items():
        if sha(ROOT/path)!=h:raise ValueError('changed frozen source/input '+path)
    if sha(GEOM/'coordinate_boxes.json')!=p['geometry_boxes_sha256']:raise ValueError('changed geometry boxes')
    record=read(ROOT/P.RECORD);target=T.Target(read(ROOT/record['point_cache_file']));rule=target.auxiliary_rule();Y=field_nodes(read(GEOM/'coordinate_boxes.json'));outputs=[]
    write(OUT/'started.json',{'protocol_sha256':sha(OUT/'protocol.json'),'utc':P.now()})
    for selector in p['selectors']:
        result=evaluate(target,selector['element'],Y,rule,selector,p);dest=OUT/f'element_{selector["element"]:03d}.json';write(dest,result);outputs.append({'element':selector['element'],'file':dest.name,'sha256':sha(dest)});print(json.dumps({'phase':'uniform_element_complete','element':selector['element'],'seconds':result['seconds']}),flush=True)
    for path,h in p['source_and_input_hashes'].items():
        if sha(ROOT/path)!=h:raise ValueError('input changed during uniform pilot '+path)
    write(OUT/'results.json',{'status':'SEVEN_UNIFORM_BULK_ELEMENTS_UNREVIEWED','protocol_sha256':sha(OUT/'protocol.json'),'elements':outputs,'all_frozen_bindings_current':True,'scope':p['scope'],'assembled611_bound_or_root_certified':False})
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['prepare','run']);a=ap.parse_args();prepare() if a.command=='prepare' else run()
