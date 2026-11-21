"""Fresh finite-cache derivative bound on the frozen common B comparison tube.
Preparation never assembles the611-coordinate Hessian. A reviewed run proves
only a cached-Q derivative bound, not a root or exact-integral transfer.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,ast,base64,datetime,hashlib,io,json,os,time,traceback
import numpy as np
from flint import arb,ctx
import parametric_root_cells as P
core=P.core
ROOT=Path(__file__).resolve().parents[1];SOURCE=Path(__file__)
OUT=ROOT/'runs/exact_integration_bridge_v1/cached_Q_common_B_v1'
RETUBE=ROOT/'runs/exact_integration_bridge_v1/retube_B_preparation'
GEOM=ROOT/'runs/exact_integration_bridge_v1/retube_B_geometry'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def q(v):return Q.from_float(float(v))if isinstance(v,(float,np.floating))else Q(v)
def write(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():raise ValueError('preserve '+str(p))
    p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n',encoding='utf8')
def alltrue(c):return isinstance(c,dict)and bool(c)and all(v is True for v in c.values())
def validate_review(review):
    if not isinstance(review,dict)or review.get('status')!='PASS'or not alltrue(review.get('checks'))or review.get('source_sha256')!=sha(SOURCE)or review.get('protocol_sha256')!=sha(OUT/'protocol.json')or review.get('preflight_sha256')!=sha(OUT/'preflight.json'):raise ValueError('exact-bound nonempty literalTrue parent prelaunch review required')
def assert_inputs(p):
    for name,h in p['source_and_input_hashes'].items():
        if sha(ROOT/name)!=h:raise ValueError('frozen input changed '+name)
def decode(packet):
    raw=base64.b64decode(packet['frozen_binary_cache']['data'],validate=True)
    if hashlib.sha256(raw).hexdigest()!=packet['frozen_binary_cache']['sha256']:raise ValueError('embedded binary cache changed')
    with np.load(io.BytesIO(raw),allow_pickle=False)as z:a={k:z[k].copy()for k in z.files}
    if not all(np.isfinite(v).all()for v in a.values()):raise ValueError('finite cache required')
    return a
def synthetic():
    checks={};r=Q(1,10);R=Q(1,3)
    checks['scalar_Z0_plus_H_variation_bound']=max(abs(1-R*(3+v))for v in [-r,r])==r/3
    c=2;d=Q(7,10);rr=Q(1,9);common=rr+abs(d)
    checks['candidate_and_swept_scalar_contained']=c-common<=c-d-rr<=c-d+rr<=c+common and c-common<=c+min(0,-d)-rr and c+max(0,-d)+rr<=c+common
    try:core.bulk_local([arb(0,1),arb(0),arb(0),arb(1),arb(1)],arb(1),arb(5000),arb(1));checks['cached_J_domain_guard']=False
    except ValueError:checks['cached_J_domain_guard']=True
    try:core.surface_local([arb(0,1),arb(1),arb(0)],arb(1),arb(0),arb(1));checks['cached_surface_domain_guard']=False
    except ValueError:checks['cached_surface_domain_guard']=True
    # Matrix norm majorization with signs, prior to preconditioning.
    B=[[Q(1),Q(-2)],[Q(3),Q(1)]];D=[[Q(2),Q(-1)],[Q(-1),Q(4)]]
    exact=[[sum(B[a][i]*D[a][b]*B[b][j]for a in range(2)for b in range(2))for j in range(2)]for i in range(2)]
    major=[sum(abs(B[a][i])*sum(abs(D[a][b])*sum(map(abs,B[b]))for b in range(2))for a in range(2))for i in range(2)]
    checks['local_projected_row_majorant_with_signs']=all(sum(map(abs,exact[i]))<=major[i]for i in range(2))
    source=ast.parse(SOURCE.read_text());calls={name:[]for name in ['build_point','variation']}
    for fn in source.body:
        if isinstance(fn,ast.FunctionDef):
            for call in ast.walk(fn):
                if isinstance(call,ast.Call)and isinstance(call.func,ast.Attribute)and isinstance(call.func.value,ast.Name)and call.func.value.id=='P'and call.func.attr in calls:calls[call.func.attr].append(fn.name)
    checks['all611_build_and_variation_only_in_gated_run']=calls=={'build_point':['run'],'variation':['run']}
    checks['no_unsafe_sparse_polynomial_import']=not any(isinstance(x,ast.Import)and any(y.name.startswith(('exact_fe_target','polynomial_jet'))for y in x.names)for x in source.body)
    assert alltrue(checks),checks
    return checks
def prepare():
    if (OUT/'protocol.json').exists():raise ValueError('preparation already frozen')
    ctx.prec=192;core.A.cache_clear();checks=synthetic();row=read(RETUBE/'protocol.json')['selected_record'];gp=read(GEOM/'protocol.json');box=read(GEOM/'coordinate_boxes.json')
    cert=read(ROOT/row['point_cache_file']);a=decode(cert)
    with np.load(ROOT/row['predictor_file'],allow_pickle=False)as z:S=z['S'].copy()
    with np.load(RETUBE/'candidate_vectors.npz',allow_pickle=False)as z:u0=z['u0'].copy();u1=z['u1'].copy()
    r=Q(gp['old_radius_exact']);d=[q(x)-q(y)for x,y in zip(u0,u1)];dn=max(map(abs,d));common=Q(gp['new_common_enclosing_radius_exact'])
    checks['same_original_611_u0_R_S']=bool(a['u'].shape==(611,)and np.array_equal(a['u'],u0)and a['preconditioner'].shape==(611,611)and S.shape==(611,5)and np.isfinite(S).all())
    checks['same_exact_cell_base_dose']=row['parameter_cell']==gp['parameter_cell']and row['parameter_cell_sha256']==gp['parameter_cell_sha256']and row['base_parameters_exact']==gp['base_parameters_exact']and row['record_id']=='optoGEF2_small_col_B'
    checks['same_exact_common_radius_r_plus_applied_decrement']=common==r+dn and dn==Q(gp['actual_applied_decrement_norm_exact'])
    cell=P.check_cell(row['parameter_cell']);base=[Q(row['base_parameters_exact'][k])for k in P.NAMES];deltas=[(cell[k][0]-b,cell[k][1]-b)for k,b in zip(P.NAMES,base)];free={int(k):i for i,k in enumerate(a['free'])};rows=[]
    def contains(o,i):return o[0]<=i[0]<=i[1]<=o[1]
    for dof in range(650):
        node,c=divmod(dof,2);x=q(a['reference'][node,c]);lo=hi=Q(0);shift=Q(0);rr=Q(0);big=Q(0)
        if dof in free:
            j=free[dof];x+=q(a['u'][j]);pairs=[(q(v)*p,q(v)*m)for v,(p,m)in zip(S[j],deltas)];lo=sum(min(z)for z in pairs);hi=sum(max(z)for z in pairs);shift=d[j];rr=r;big=common
        cm=(x+lo-big,x+hi+big);old=(x+lo-rr,x+hi+rr);candidate=(x-shift+lo-rr,x-shift+hi+rr);path=(x+lo+min(0,-shift),x+hi+max(0,-shift));swept=(path[0]-rr,path[1]+rr);gc=Q(box['current_centers_exact'][node][c]);gr=Q(box['coordinate_radii_exact'][node][c]);gb=(gc-gr,gc+gr)
        good={'old_cube':contains(cm,old),'candidate_cube':contains(cm,candidate),'translation_path':contains(cm,path),'swept_cube':contains(cm,swept),'canonical_geometry_box':contains(gb,cm)}
        if not alltrue(good):raise AssertionError((dof,good))
        rows.append({'global_dof':dof,'free':dof in free,'checks':good,'common_exact_interval':list(map(str,cm))})
    checks['all650_exact_containment_rows']=len(rows)==650 and all(alltrue(x['checks'])for x in rows)
    checks['base_in_closed_C']=all(lo<=b<=hi for(lo,hi),b in zip([cell[k]for k in P.NAMES],base))
    checks['fixed_mu1_kappa5000']=q(cert['parameters']['mu'])==1 and q(cert['parameters']['bulk'])==5000
    oldroot=read(ROOT/row['root_certificate_file']);checks['exact_gain_times_source_dose_not_solver_rounding']=base[4]*Q(row['dose_exact'])==Q(oldroot['base_strength_exact'])
    inputs={}
    paths=[SOURCE,Path(P.__file__),Path(core.__file__),RETUBE/'protocol.json',RETUBE/'candidate_vectors.npz',GEOM/'protocol.json',GEOM/'coordinate_boxes.json',GEOM/'certification_run/results.json',ROOT/'evidence/exact_retube_B_geometry_independent_review.json',ROOT/row['point_cache_file'],ROOT/row['predictor_file'],ROOT/row['root_certificate_file'],ROOT/cert['input_state']]
    for path in paths:inputs[path.relative_to(ROOT).as_posix()]=sha(path)
    checks['authoritative_cache_predictor_hashes']=sha(ROOT/row['point_cache_file'])==row['point_cache_sha256']==oldroot['artifact_hashes'][row['point_cache_file']]and sha(ROOT/row['predictor_file'])==row['predictor_sha256']==oldroot['artifact_hashes'][row['predictor_file']]
    state=read(ROOT/cert['input_state']);checks['original_state_hash_and_free_vector']=sha(ROOT/cert['input_state'])==cert['input_state_sha256']and np.array_equal(a['u'],np.asarray(state['free_displacement']))and np.array_equal(a['reference'],np.asarray(state['reference_nodes']))
    gr=read(ROOT/'evidence/exact_retube_B_geometry_independent_review.json');checks['fresh_geometry_review_bound']=gr['status'].startswith('PASS')and alltrue(gr['checks'])and gr['protocol_sha256']==sha(GEOM/'protocol.json')and gr['producer_results_sha256']==sha(GEOM/'certification_run/results.json')
    assert alltrue(checks),checks
    write(OUT/'exact_containment.json',{'status':'PASS_EXACT_COMMON_TUBE_CONTAINMENT','rows':rows,'r_exact':str(r),'decrement_norm_exact':str(dn),'common_radius_exact':str(common),'no_integral_Newton_step_interpretation_required':True})
    build={'point_certificate_file':row['point_cache_file'],'predictor_file':row['predictor_file'],'dose_exact':row['dose_exact'],'base_parameters_exact':row['base_parameters_exact']}
    p={'version':1,'utc':now(),'source_and_input_hashes':inputs,'record_id':row['record_id'],'parameter_cell':row['parameter_cell'],'parameter_cell_sha256':row['parameter_cell_sha256'],'base_parameters_exact':row['base_parameters_exact'],'dose_exact':row['dose_exact'],'base_strength_exact':str(base[4]*Q(row['dose_exact'])),'build_point_parameters':build,'predictor_path_policy':'Resolve ROOT/predictor_file to absolute path in memory before pure P.build_point; P.OUT is never modified.','precision_bits':192,'old_radius_exact':str(r),'actual_applied_decrement_norm_exact':str(dn),'common_radius_exact':str(common),'geometry_domain_id':box['domain_id'],'geometry_boxes_sha256':sha(GEOM/'coordinate_boxes.json'),'exact_containment_sha256':sha(OUT/'exact_containment.json'),'target':'The unchanged finite binary-cache Q potential on u=p0(theta)+e, theta in SAME C, ||e||infinity<=common_radius; includes closed straight segments to the base.','bound':'Fresh Z0>=||I-R H_Q(u0,theta0)||infinity plus fresh preconditioned Hessian-variation row majorants on the full common tube.','no_sparse_interval_polynomial_used':True,'arbitrary_candidate_interpretation':'Candidate vector and radius are fixed geometric inputs. Its historical origin does not supply a validated integral-residual premise.','actual_root_or_selfmap_claim':False,'integral_H_or_PDE_bound':False,'full_identified_set_or_all17_computed':False,'prelaunch_review_required':True}
    write(OUT/'protocol.json',p);write(OUT/'preflight.json',{'status':'PASS_PREPARATION_ONLY','checks':checks,'source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json'),'containment_rows':650,'exact_containments':3250,'actual_611_point_assembly_or_contraction_run':False});print(json.dumps({'status':'PASS_PREPARATION_ONLY','checks':len(checks),'source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json')}),flush=True)
def cache_domains(point,C,r):
    rows=[];minJ=None;minR=None;minN=None;bulk=surface=0
    for i,term in enumerate(point['cache']):
        k=len(term['q']);Pstrain=term['predictor_strain'];offset=[sum((Pstrain[a,b]*C['delta_hull'][b]for b in range(5)),arb(0))for a in range(k)];z=[term['q'][a]+offset[a]+arb(0,(r*term['rowsumB'][a]).upper())for a in range(k)]
        if term['kind']=='bulk':
            det=z[0]*z[3]-z[1]*z[2];J=det*z[4]
            if not det>0 or not z[4]>0 or not J>core.A(Q(1,10**6)):raise ValueError('fresh cached bulk convex-domain failure')
            value=Q(str(J.lower().fmpq()));minJ=value if minJ is None else min(minJ,value);bulk+=1;detail={'det_lower_exact':str(det.lower().fmpq()),'hoop_lower_exact':str(z[4].lower().fmpq()),'J_lower_exact':str(value)}
        else:
            n=z[1]*z[1]+z[2]*z[2]
            if not z[0]>0 or not n>0 or not term['parameters'][2]>0:raise ValueError('fresh cached surface convex-domain failure')
            rr=Q(str(z[0].lower().fmpq()));nn=Q(str(n.lower().fmpq()));minR=rr if minR is None else min(minR,rr);minN=nn if minN is None else min(minN,nn);surface+=1;detail={'cached_radial_lower_exact':str(rr),'cached_tangent_norm2_lower_exact':str(nn)}
        if not term['weight']>0:raise ValueError('positive normalized weights required')
        rows.append({'term':i,'kind':term['kind'],**detail})
    return {'status':'PASS_FRESH_CACHED_Q_CONVEX_DOMAINS','terms':len(rows),'bulk_terms':bulk,'surface_terms':surface,'minimum_bulk_J_lower_exact':str(minJ),'minimum_surface_radial_lower_exact':str(minR),'minimum_surface_tangent_norm2_lower_exact':str(minN),'rows':rows,'scope':'Cached Q analytic domains, separate from canonical whole-body geometry; not an integral-root claim.'}
def run():
    p=read(OUT/'protocol.json');validate_review(read(OUT/'prelaunch_independent_review.json'));assert_inputs(p)
    if (OUT/'job.json').exists()or(OUT/'results.json').exists():raise ValueError('one versioned job only')
    import torch
    torch.set_num_threads(1);ctx.prec=p['precision_bits'];core.A.cache_clear();start=time.monotonic();write(OUT/'job.json',{'status':'RUNNING','pid':os.getpid(),'utc':now(),'protocol_sha256':sha(OUT/'protocol.json'),'source_sha256':sha(SOURCE)})
    try:
        build=dict(p['build_point_parameters']);build['predictor_file']=str((ROOT/build['predictor_file']).resolve());point=P.build_point(build)
        write(OUT/'fresh_point_summary.json',{'Z0':core.pack(point['Z0']),'free_dofs':point['n'],'terms':len(point['cache']),'base_strength_exact':p['base_strength_exact'],'seconds_after_point_setup':time.monotonic()-start,'root_certified':False})
        C=P.cell_context(point,p['parameter_cell'],'affine');radius=core.A(Q(p['common_radius_exact']));domain=cache_domains(point,C,radius);write(OUT/'cached_domain_bounds.json',domain)
        variation,minJ=P.variation(point,C,radius);z0=Q(str(point['Z0'].upper().fmpq()));v=[Q(str(x.upper().fmpq()))for x in variation];uniform=[z0+x for x in v];ZQ=max(uniform)
        write(OUT/'fresh_row_bounds.json',{'Z0_upper_exact':str(z0),'variation_row_upper_exact':list(map(str,v)),'uniform_row_upper_exact':list(map(str,uniform)),'ZQ_upper_exact':str(ZQ),'common_radius_exact':p['common_radius_exact'],'scope':'Cached-Q derivative rows only; not a self-map or integral-root proof.'})
        checks={'same611_free_coordinates':point['n']==611 and len(v)==611,'all9504_cached_domains_positive':domain['terms']==9504 and domain['bulk_terms']==9216 and domain['surface_terms']==288,'variation_domain_J_agrees':Q(str(minJ.lower().fmpq()))==Q(domain['minimum_bulk_J_lower_exact']),'nonnegative_finite_row_bounds':z0>=0 and all(x>=0 for x in v),'ZQ_exact_Z0_plus_max_row_variation':ZQ==z0+max(v),'all_frozen_inputs_current':all(sha(ROOT/n)==h for n,h in p['source_and_input_hashes'].items())}
        if not alltrue(checks):raise AssertionError(checks)
        out={'status':'CERTIFIED_FRESH_CACHED_Q_CONTRACTION_BOUND'if ZQ<1 else'VALIDATED_CACHED_Q_DERIVATIVE_BOUND_NOT_CONTRACTIVE','checks':checks,'utc':now(),'source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json'),'prelaunch_review_sha256':sha(OUT/'prelaunch_independent_review.json'),'record_id':p['record_id'],'parameter_cell':p['parameter_cell'],'parameter_cell_sha256':p['parameter_cell_sha256'],'dose_exact':p['dose_exact'],'base_strength_exact':p['base_strength_exact'],'common_radius_exact':p['common_radius_exact'],'Z0_upper_exact':str(z0),'variation_row_upper_exact':list(map(str,v)),'uniform_row_upper_exact':list(map(str,uniform)),'ZQ_upper_exact':str(ZQ),'contraction_less_than_one':ZQ<1,'cached_domain_bounds_sha256':sha(OUT/'cached_domain_bounds.json'),'point_summary_sha256':sha(OUT/'fresh_point_summary.json'),'seconds':time.monotonic()-start,'quantified_claim':'For every theta in SAME closed C and every u=p0(theta)+e with ||e||infinity<=common_radius, ||I-R H_Q(u,theta)||infinity<=ZQ.','root_selfmap_or_integral_transfer_certified':False,'integral_H_error_bound_computed':False,'full_17_or_Theta_or_PDE_computed':False}
        assert_inputs(p);write(OUT/'results.json',out);write(OUT/'terminal.json',{'status':'FINISHED','pid':os.getpid(),'utc':now(),'results_sha256':sha(OUT/'results.json'),'seconds':out['seconds']});print(json.dumps({'status':out['status'],'Z0':float(z0),'max_variation':float(max(v)),'ZQ':float(ZQ),'seconds':out['seconds'],'results_sha256':sha(OUT/'results.json')}),flush=True)
    except BaseException as error:
        write(OUT/'failure.json',{'status':'UNRESOLVED_FRESH_Q_BOUND','error':str(error),'traceback':traceback.format_exc(),'source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json'),'seconds':time.monotonic()-start});write(OUT/'terminal.json',{'status':'FAILED','pid':os.getpid(),'utc':now(),'failure_sha256':sha(OUT/'failure.json')});raise
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['prepare','run']);x=ap.parse_args();prepare()if x.command=='prepare'else run()
