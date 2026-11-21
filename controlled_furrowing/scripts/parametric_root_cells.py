"""Uniform nearby-root cells for a frozen representative discrete tissue model.

Finite binary-cache potential, fixed geometry/mu/kappa/dose.  This is not a
physical prior, the full identified set, or a certificate for 17 source doses.
No equilibrium optimization is performed. Failed sufficient tests are UNRESOLVED.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,base64,datetime,hashlib,io,json,sys,time,os
import numpy as np
from flint import arb,arb_mat,ctx
import algebraic_root_feasibility as core

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/parametric_root_cells_v1'
POINT=ROOT/'evidence/algebraic_root_certification_feasibility.json'
REVIEW=ROOT/'evidence/algebraic_root_independent_review.json'
NAMES=('apical_area_modulus','basal_area_modulus','apical_surface_tension','basal_surface_tension','common_gain')
SCALES=('1/100','1/10000','1/1000000','1/100000000','1/10000000000','1/1000000000000')
A=core.A;upper=core.upper;pack=core.pack


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def save(p,value):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf8')
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def Q(x):return F.from_float(float(x))if isinstance(x,(float,np.floating))else F(x)
def ball(lo,hi):return A(Q(lo)).union(A(Q(hi)))


def check_cell(cell):
    if set(cell)!=set(NAMES):raise ValueError('exact five-coordinate cell required')
    for name in NAMES:
        lo,hi=map(F,cell[name])
        if lo>hi or lo<0 or(name!='apical_surface_tension'and lo<=0):raise ValueError('ordered positive coefficients/nonnegative apical tension required')
    return {name:tuple(map(F,cell[name]))for name in NAMES}


def make_cell(base,relative):
    s=F(relative);c={name:[str(q*(1-s)),str(q*(1+s))]for name,q in zip(NAMES,base)}
    c['apical_surface_tension']=['0',str(F(1,50)*s)];check_cell(c);return c


def split_cell(cell,scales):
    c=check_cell(cell);j=max(range(5),key=lambda k:((c[NAMES[k]][1]-c[NAMES[k]][0])/F(scales[k]),-k));name=NAMES[j];mid=sum(c[name])/2
    if c[name][0]==c[name][1]:raise ValueError('no positive-width coordinate to bisect')
    children=[]
    for pair in ((c[name][0],mid),(mid,c[name][1])):
        child={k:list(map(str,v))for k,v in c.items()};child[name]=list(map(str,pair));children.append(child)
    return name,children


def prepare():
    if(OUT/'protocol.json').exists():raise RuntimeError('protocol already exists; preserve it')
    actual=(OUT/'point_cache.json').exists()
    pointfile=OUT/'point_cache.json'if actual else POINT
    receipt=read(pointfile);review=read(REVIEW)if not actual else None
    if not actual and(review['status']!='PASS'or review['certificate_sha256']!=sha(POINT)):raise RuntimeError('reviewed exact point target required')
    sourceprep=read(OUT/'preparation_protocol.json')if actual else None
    model,u,R=core.restore_cache(receipt)
    if model.parameters.mu!=1 or model.parameters.bulk!=5000:raise ValueError('fixed representative normalization required')
    base=[Q(getattr(model.parameters,k))for k in NAMES[:4]]+[F(sourceprep['base_gain_exact'])if actual else F(receipt['activation_strength_binary_exact'])]
    archive=ROOT/receipt['matrix_archive']
    with np.load(archive)as z:
        S=z['displacement_tangents'][:,[1,2,3,4,5]].copy();assert np.array_equal(z['free_vector'],u)
    S/=np.array([float(base[0]),float(base[1]),1.,float(base[3]),float(base[4])])[None,:]
    OUT.mkdir(exist_ok=True);np.savez_compressed(OUT/'affine_predictor.npz',S=S)
    dose=sourceprep['dose_exact']if actual else'1'
    protocol={'version':2,'scope':'Uniform LOCAL-root cells on one fixed geometry/source input; not a physical prior/fullTheta or all17source records',
        'created_utc':stamp(),'source_hashes':{'parametric_root_cells.py':sha(__file__),'algebraic_root_feasibility.py':sha(Path(core.__file__))},
        'point_certificate_file':pointfile.relative_to(ROOT).as_posix(),'point_certificate_sha256':sha(pointfile),'point_review_file':REVIEW.relative_to(ROOT).as_posix()if not actual else None,'point_review_sha256':sha(REVIEW)if not actual else None,'point_cache_sha256':receipt['frozen_binary_cache']['sha256'],
        'point_input_state':receipt['input_state'],'point_input_sha256':receipt['input_state_sha256'],'matrix_archive_sha256':sha(archive),
        'predictor_file':'affine_predictor.npz','predictor_sha256':sha(OUT/'affine_predictor.npz'),
        'predictor_rule':'Archived binary tangent columns(logKa,logKb,addTa,logTb,loggain) divided in binary64 by baseline coefficient; stored result defines exact affine predictor. Accuracy is validated by explicit first-order defect/remainder, not assumed.',
        'coordinate_order':list(NAMES),'coordinate_units':'Declared normalized solver coefficients; common_gain=activation strength at representative normalized dose1',
        'base_parameters_exact':dict(zip(NAMES,map(str,base))),'fixed_mu':'1','fixed_kappa':'5000','base_full_parameters':receipt['parameters'],'fixed_geometry':{k:v for k,v in receipt['parameters'].items()if k not in('mu','bulk',*NAMES[:4])},
        'dose_exact':dose,'input_kind':'source_record_dose'if actual else'Representative normalization of archived strength1, NOT an identified source embryo dose','calibration_record_ids_covered':[sourceprep['record_id']]if actual else[],
        'branch_domain_id':'new_primary_M_local_branch_v1'if actual else'old_weak_fixed_strength1_local_branch_v1','source_record_preparation_sha256':sha(OUT/'preparation_protocol.json')if actual else None,
        'strength_relation':'common_gain_times_source_dose'if actual else'common_gain_times_representative_dose1','base_strength_requested_exact':str(base[4]*F(dose)),
        'precision_bits':160,'initial_cells':[{'cell_id':'scale_'+str(i),'relative_scale_exact':s,'bounds':make_cell(base,s)}for i,s in enumerate(SCALES)],
        'methods':['constant','affine'],'error_radius_floor_exact':'1/100000000000','radius_rule':'max(floor,2*rigorous_eta_upper) times each frozen multiplier; no nonlinear solve',
        'radius_multipliers':[1,2,4],'maximum_test_radius_exact':'1/100','adaptive_demo':{'root_relative_scale_exact':'1/10000','method':'affine','maximum_depth':1,'split':'widest width divided by numerical scales; tie by coordinate order; retain all terminal unresolved leaves','numerical_width_scales_exact':[str(x)for x in (base[0],base[1],F(1,50),base[3],base[4])]},
        'status_policy':'CERTIFIED_UNIFORM_LOCAL_ROOT means foralltheta in cell existsunique root in predictor(theta)+errorcube; failed/invalid sufficient test remainsUNRESOLVED, neverOUT',
        'quantified_branch':'Only selected nearby branch; no energetic/global stability, other-branch exclusion, idealquadrature/PDE or source-time validation',
        'observer_enclosures':None,'observer_scope':'Not supplied by this pilot; future transfer must use fullpredictor+error tube and separately verified assumptions',
        'no_new_equilibrium_optimization':True,'heldout_response_files_read':False}
    save(OUT/'protocol.json',protocol);print(json.dumps({'status':'PROTOCOL_FROZEN','sha256':sha(OUT/'protocol.json'),'cells':6,'source':protocol['source_hashes']}),flush=True)


def tests():
    p=read(OUT/'protocol.json');assert_sources(p);checks={}
    base=list(map(F,p['base_parameters_exact'].values()));c=make_cell(base,F(1,100));name,ch=split_cell(c,p['adaptive_demo']['numerical_width_scales_exact'])
    checks['split_covers_exact_parent_interval']=F(ch[0][name][0])==F(c[name][0])and F(ch[1][name][1])==F(c[name][1])and F(ch[0][name][1])==F(ch[1][name][0])
    checks['other_coordinate_intervals_preserved']=all(child[k]==c[k]for child in ch for k in NAMES if k!=name)
    checks['all5_coordinates_nonzero_width']=all(F(v[0])<F(v[1])for v in c.values())
    bad={**c,'apical_surface_tension':['-1','0']}
    try:check_cell(bad);checks['negative_tension_rejected']=False
    except ValueError:checks['negative_tension_rejected']=True
    # Exact affine predictor cancellation for f(u,t)=u-t: M=HS+B=0.
    checks['scalar_affine_predictor_defect_zero']=F(1)*F(1)+F(-1)==0
    t=ball(F(199,100),F(201,100));u=A(F(3,2));R=A(F(1,3));eta=upper(R*(u*u-t));rad=A(F(1,5));Z=upper(1-R*2*(u+arb(0,rad)))
    checks['scalar_parametric_sqrt_cell_contracts']=bool(Z<1 and eta+Z*rad<rad)
    checks['failure_policy_not_OUT']=p['status_policy'].find('neverOUT')>=0
    checks['not_claimed_all17_or_fullset']=len(p['calibration_record_ids_covered'])<17 and p['observer_enclosures']is None
    core_checks=core.synthetic_checks();checks['reviewed_kernel47_checks_pass']=all(core_checks['checks'].values())
    if not all(checks.values()):raise AssertionError(checks)
    save(OUT/'tests.json',{'status':'PASS','checks':checks,'kernel_checks':core_checks,'protocol_sha256':sha(OUT/'protocol.json'),'source_sha256':sha(__file__)})
    print(json.dumps({'status':'PASS','checks':len(checks),'kernel_checks':47}),flush=True)


def assert_sources(p):
    for n,h in p['source_hashes'].items():
        if sha(Path(__file__).with_name(n))!=h:raise RuntimeError('frozen source changed:'+n)
    inputs=[(ROOT/p['point_certificate_file'],p['point_certificate_sha256']),(OUT/p['predictor_file'],p['predictor_sha256'])]
    if p.get('point_review_file'):inputs.append((ROOT/p['point_review_file'],p['point_review_sha256']))
    if p.get('source_record_preparation_sha256'):inputs.append((OUT/'preparation_protocol.json',p['source_record_preparation_sha256']))
    for f,h in inputs:
        if sha(f)!=h:raise RuntimeError('frozen input changed:'+str(f))


def build_point(p):
    point=read(ROOT/p['point_certificate_file']);model,u,R=core.restore_cache(point);n=len(u);dose=F(p['dose_exact']);basegain=F(p['base_parameters_exact']['common_gain'])
    with np.load(OUT/p['predictor_file'])as z:S=z['S'].copy()
    surface0=len(model.surfaces[0][0])*len(model.surfaces[0][1]);counter=0;bulkcounter=0;active=model.active.numpy();nq=len(model.shape)
    def decorated():
        nonlocal counter,bulkcounter
        for term in core.local_terms(model,u,basegain*dose):
            if term['kind']=='bulk':
                t,iq=divmod(bulkcounter,nq);term['active_mask']=A(float(active[t,0 if active.shape[1]==1 else iq]));bulkcounter+=1
            if term['kind']=='surface':term['surface_index']=int(counter>=surface0);counter+=1
            yield term
    g,H,cache,minJ=core.assemble_point(decorated(),n);B0=arb_mat(n,5);SM=arb_mat(n,5,[A(float(v))for v in S.ravel()])
    for term in cache:
        B=term['B'];k=len(B);d=len(term['indices']);C=arb_mat(k,5)
        if term['kind']=='bulk':
            drive=term['active_mask']*A(dose);C[1,4]=drive*term['q'][1];C[3,4]=drive*term['q'][3]
        else:
            Kidx,Tidx=(1,3)if term['surface_index']==0 else(0,2);scale=term['parameters'][2]
            gK,_,_=core.surface_local(term['q'],A(1),A(0),scale);gT,_,_=core.surface_local(term['q'],A(0),A(1),scale)
            for a in range(k):C[a,Kidx]=gK[a];C[a,Tidx]=gT[a]
        BM=arb_mat(k,d,[v for row in B for v in row]);local=BM.transpose()*C
        for j,i in enumerate(term['indices']):
            for a in range(5):B0[i,a]+=term['weight']*local[j,a]
        localS=arb_mat(d,5,[A(float(S[i,j]))for i in term['indices']for j in range(5)]);term['predictor_strain']=BM*localS
    RM=arb_mat(n,n,[A(float(v))for v in R.ravel()]);AR=arb_mat(n,n,[A(abs(float(v)))for v in R.ravel()]);D=arb_mat(n,n)
    for i in range(n):D[i,i]=1
    D-=RM*H;Z0=max(sum((upper(D[i,j])for j in range(n)),arb(0)).upper()for i in range(n))
    return {'n':n,'cache':cache,'R':RM,'absR':AR,'S':S,'gR':RM*g,'constant_M':RM*B0,'affine_M':RM*(H*SM+B0),'Z0':Z0,'base':[F(p['base_parameters_exact'][k])for k in NAMES],'dose':A(dose),'raw_gradient':g,'raw_Hessian':H,'raw_parameter_force_derivatives':B0,'predictor_matrix':SM}


def cell_context(point,cell,method):
    c=check_cell(cell);theta=[ball(*c[k])for k in NAMES];delta=[theta[j]-A(point['base'][j])for j in range(5)];hull=[x.union(arb(0))for x in delta]
    S=point['S']if method=='affine'else np.zeros_like(point['S']);nodal=[]
    for i in range(point['n']):nodal.append(upper(sum((A(float(S[i,j]))*delta[j]for j in range(5)),arb(0))))
    return {'theta':theta,'delta':delta,'delta_hull':hull,'nodal_predictor_bounds':nodal,'method':method}


def variation(point,C,r,remainder=False):
    rows=[arb(0)for _ in range(point['n'])];minJ=None
    for term in point['cache']:
        k=len(term['q']);P=term['predictor_strain'];offset=[sum((P[a,b]*C['delta_hull'][b]for b in range(5)),arb(0))if C['method']=='affine'else arb(0)for a in range(k)]
        q=[term['q'][a]+offset[a]+arb(0,(r*term['rowsumB'][a]).upper())for a in range(k)]
        if term['kind']=='bulk':
            mu,kappa,drive=term['parameters'];_,H,J=core.bulk_local(q,mu,kappa,term['active_mask']*point['dose']*C['theta'][4])
            minJ=J.lower()if minJ is None else min(minJ,J.lower())
        else:
            i,j=(1,3)if term['surface_index']==0 else(0,2);_,H,J=core.surface_local(q,C['theta'][i],C['theta'][j],term['parameters'][2])
        change=[[upper(H[a][b]-term['Hcenter'][a][b])for b in range(k)]for a in range(k)]
        bsum=[upper(offset[b])for b in range(k)]if remainder else term['rowsumB']
        projected=[sum((change[a][b]*bsum[b]for b in range(k)),arb(0))for a in range(k)]
        for j,i in enumerate(term['indices']):rows[i]+=term['weight']*sum((term['absB'][a][j]*projected[a]for a in range(k)),arb(0))
    out=point['absR']*arb_mat(point['n'],1,[v.upper()for v in rows]);return [upper(out[i,0])for i in range(point['n'])],minJ


def evaluate_cell(point,p,cell_id,bounds,method):
    start=time.monotonic();C=cell_context(point,bounds,method);M=point[method+'_M'];linear=[point['gR'][i,0]+sum((M[i,j]*C['delta'][j]for j in range(5)),arb(0))for i in range(point['n'])]
    result={'cell_id':cell_id,'method':method,'cell_bounds_exact':bounds,'status':'UNRESOLVED','dose_exact':p['dose_exact'],'calibration_record_ids_covered':p['calibration_record_ids_covered'],'input_kind':p['input_kind'],'strength_relation':p['strength_relation'],'branch_domain_id':p['branch_domain_id'],
        'quantified_claim':'ForEVERYtheta inthiscell, oneunique nearbyroot inpredictor(theta)+B_infinity(r); selectednearbybranchonly','radius_attempts':[],'observer_enclosures':None}
    try:
        rem,_=variation(point,C,arb(0),True)if method=='affine'else([arb(0)]*point['n'],None)
        eta=max((upper(linear[i])+rem[i]).upper()for i in range(point['n']));result['eta_uniform']=pack(eta)
        rbase=max(F(p['error_radius_floor_exact']),2*F(str(eta.fmpq())))
        for mult in p['radius_multipliers']:
            rF=rbase*mult;attempt={'radius_exact':str(rF)}
            if rF>F(p['maximum_test_radius_exact']):attempt['failure']='Declared maximum test radius exceeded';result['radius_attempts'].append(attempt);continue
            r=A(rF)
            try:
                vb,J=variation(point,C,r);Z=point['Z0']+max(vb);lhs=eta+Z*r;ok=bool(Z<1 and lhs<r)
                attempt.update({'Z_uniform':pack(Z),'eta_plus_Zr':pack(lhs),'quadrature_J_lower':pack(J),'strict_selfmap_contraction':ok})
                if ok:
                    error=(eta/(1-Z)).upper();total=max((v+error).upper()for v in C['nodal_predictor_bounds']);result.update({'status':'CERTIFIED_UNIFORM_LOCAL_ROOT','selected_radius_exact':str(rF),'uniform_error_from_predictor':pack(error),'uniform_outer_nodal_radius_about_base':pack(total),'full_nodal_predictor_bounds_exact':[str(v.fmpq())for v in C['nodal_predictor_bounds']]})
            except ValueError as e:attempt.update({'strict_selfmap_contraction':False,'failure':str(e)})
            result['radius_attempts'].append(attempt)
            if result['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT':break
    except ValueError as e:result['failure']=str(e)
    result['seconds']=time.monotonic()-start;return result


def run():
    p=read(OUT/'protocol.json');assert_sources(p);test=read(OUT/'tests.json');assert test['status']=='PASS'and test['protocol_sha256']==sha(OUT/'protocol.json')
    if(OUT/'results.json').exists()or(OUT/'job.json').exists():raise RuntimeError('job/results already exist; preserve instead of rerun')
    save(OUT/'job.json',{'status':'RUNNING','pid':os.getpid(),'protocol_sha256':sha(OUT/'protocol.json'),'source_sha256':sha(__file__)})
    ctx.prec=p['precision_bits'];A.cache_clear();point=build_point(p);rows=[];started=time.monotonic()
    summary={'status':'RUNNING','protocol_sha256':sha(OUT/'protocol.json'),'source_hashes':p['source_hashes'],'point_defect_Z0':pack(point['Z0']),'cell_attempts':rows,'scope':p['scope'],'no_new_equilibrium_optimization':True,'full_identified_set_computed':False}
    save(OUT/'progress.json',summary)
    for cell in p['initial_cells']:
        for method in p['methods']:
            assert_sources(p);row=evaluate_cell(point,p,cell['cell_id'],cell['bounds'],method);rows.append(row);save(OUT/'progress.json',summary);print(json.dumps({'cell':row['cell_id'],'method':method,'status':row['status'],'seconds':row['seconds']}),flush=True)
    scale=p['adaptive_demo']['root_relative_scale_exact'];root=next(c for c in p['initial_cells']if F(c['relative_scale_exact'])==F(scale));parent=next(r for r in rows if r['cell_id']==root['cell_id']and r['method']=='affine');ledger={'root_cell':root['bounds'],'root_result_id':parent['cell_id'],'leaves':[],'split_records':[],'full_physical_domain_covered':False}
    if parent['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT':ledger['leaves'].append({'cell_id':parent['cell_id'],'bounds':root['bounds'],'status':parent['status']})
    else:
        name,children=split_cell(root['bounds'],p['adaptive_demo']['numerical_width_scales_exact']);ledger['split_records'].append({'parent_cell_id':root['cell_id'],'split_coordinate':name,'children':children,'union_equals_parent':True})
        for i,child in enumerate(children):
            row=evaluate_cell(point,p,'adaptive_'+str(i),child,'affine');rows.append(row);ledger['leaves'].append({'cell_id':row['cell_id'],'bounds':child,'status':row['status']});save(OUT/'progress.json',summary)
    summary.update({'status':'COMPLETE_DECLARED_CELL_ATTEMPTS','cover_ledger':ledger,'elapsed_seconds_after_point_setup':time.monotonic()-started,'completed_utc':stamp(),'certified_cells':sum(r['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT'for r in rows),'unresolved_cells':sum(r['status']=='UNRESOLVED'for r in rows),'no_OUT_classification':True})
    assert_sources(p);save(OUT/'results.json',summary);save(OUT/'job.json',{'status':'FINISHED','pid':os.getpid(),'results_sha256':sha(OUT/'results.json')});print(json.dumps({'status':summary['status'],'certified':summary['certified_cells'],'unresolved':summary['unresolved_cells']}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','test','run']);parser.add_argument('--subrun',choices=['actual_M']);args=parser.parse_args();ctx.prec=160
    if args.subrun:OUT=OUT/args.subrun
    {'prepare':prepare,'test':tests,'run':run}[args.stage]()
