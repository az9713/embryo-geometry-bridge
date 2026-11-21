"""Joint integration/mesh check for a preselected actual strong-load case.

The 16x8 default-rule state is only an initial guess under Q8. No target data,
gain fits, or 24x12 jobs. Sparse Newton must be independently verified first.
"""
from pathlib import Path
from fractions import Fraction as Q
from datetime import datetime,timezone
from dataclasses import asdict
import argparse,hashlib,json,math,os,time
import numpy as np
import torch
from quadratic_quadrature_study import Parameters,IntegrationPolicy,evaluate,exact_geometry
from robust_quadratic import SolvePolicy
from run_activated_mesh_study import QuadraticGraph,compare_graphs
from annulus_observer import serial

ROOT=Path(__file__).resolve().parents[1]
FINE=ROOT/'runs/activated_mesh/lambda_5000.0_target_10.0/final.json'
COARSE=ROOT/'runs/activated_quadrature/lambda_5000.0_strength_10.0/reequilibrated.json'
OUT=ROOT/'runs/joint_refinement'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,value):Path(path).write_text(json.dumps(serial(value),indent=2)+'\n')

def reference_curves(result,surface):
    n=max(i for t in result['triangles'] for i in t[:3])+1
    ref=result['reference_nodes'][:n];z=max(p[1] for p in ref) if surface=='top' else min(p[1] for p in ref)
    radius=sorted(set(Q(p[0]) for p in ref if p[1]==z))
    controls=result['boundary_'+surface+'_quadratic_controls']
    if len(radius)!=len(controls)+1:raise ValueError('reference edge/curve count mismatch')
    return [(radius[i],radius[i+1],[tuple(Q(v) for v in p) for p in edge]) for i,edge in enumerate(controls)]

def material_curve_difference(a,b,surface):
    aa=reference_curves(a,surface);bb=reference_curves(b,surface)
    if aa[0][0]!=bb[0][0] or aa[-1][1]!=bb[-1][1]:raise ValueError('common reference curve domain required')
    knots=sorted(set([r for curve in (aa,bb) for edge in curve for r in edge[:2]]))
    def at(curves,r):
        edge=next(v for v in curves if v[0]<=r<=v[1]);t=(r-edge[0])/(edge[1]-edge[0]);p=edge[2]
        return tuple((1-t)*(1-2*t)*p[0][j]+4*t*(1-t)*p[1][j]+t*(2*t-1)*p[2][j] for j in (0,1))
    upper=Q(0);lower=Q(0);records=[]
    for left,right in zip(knots[:-1],knots[1:]):
        delta=[]
        for r in (left,(left+right)/2,right):
            av=at(aa,r);bv=at(bb,r);delta.append(tuple(av[j]-bv[j] for j in (0,1)))
        bezier=[delta[0],tuple(2*delta[1][j]-(delta[0][j]+delta[2][j])/2 for j in (0,1)),delta[2]]
        up=max(sum(v*v for v in c) for c in bezier);lo=max(sum(v*v for v in c) for c in delta)
        upper=max(upper,up);lower=max(lower,lo)
        records.append({'reference_interval':[left,right],'difference_Bezier_controls':bezier,'squared_norm_upper':up})
    return {'squared_same_parameter_displacement_lower_exact':str(lower),'squared_same_parameter_displacement_upper_exact':str(upper),
            'lower_display_um':33.5*math.sqrt(float(lower)),'upper_display_um':33.5*math.sqrt(float(upper)),
            'records':records,'scope':'Exact common-material-parameter vector difference enclosure; upper bound also bounds geometric Hausdorff distance. Lower bound is for that parameter pairing only, NOT a Hausdorff lower bound or Eulerian vertical discrepancy.'}

def compare_results(fine,coarse):
    pf=dict(fine['parameters']);pc=dict(coarse['parameters'])
    for k in ('nr','nz'):pf.pop(k);pc.pop(k)
    if pf!=pc or fine['activation_strength']!=coarse['activation_strength']:raise ValueError('unmatched physical parameters or load')
    if fine['integration_policy']!=coarse['integration_policy']:raise ValueError('unmatched actual integration policies')
    result={'fine_mesh':[fine['parameters']['nr'],fine['parameters']['nz']],'coarse_mesh':[coarse['parameters']['nr'],coarse['parameters']['nz']],
            'shortening_difference_um':33.5*(fine['maximum_activated_column_shortening']-coarse['maximum_activated_column_shortening']),
            'central_apical_difference_um':33.5*(fine['central_apical_displacement']-coarse['central_apical_displacement']),
            'central_basal_difference_um':33.5*(fine['central_basal_displacement']-coarse['central_basal_displacement'])}
    for surface in ('top','bottom'):
        key='boundary_'+surface+'_quadratic_controls'
        result[surface+'_Eulerian_vertical']=compare_graphs(QuadraticGraph(fine[key]),QuadraticGraph(coarse[key]),tolerance_um=Q('0.01'))
        result[surface+'_material_curve']=material_curve_difference(fine,coarse,surface)
    return result

def run(verification):
    from sparse_quadratic import AxisymmetricModel,SparsePolicy
    OUT.mkdir(parents=True,exist_ok=True);live=OUT/'active_job.json'
    if live.exists() and json.loads(live.read_text()).get('status')=='RUNNING':raise RuntimeError('Joint16 job already active; retain its handle rather than duplicate')
    verified=json.loads(Path(verification).read_text())
    checks=verified.get('checks',{})
    if not checks or not all(checks.values()):raise RuntimeError('Sparse solver verification has not passed')
    sparse_path=ROOT/'scripts/sparse_quadratic.py'
    if verified.get('source_sha256')!=sha(sparse_path):raise RuntimeError('Sparse source differs from verification receipt')
    independent_path=ROOT/'evidence/joint_refinement_sparse_review.json'
    independent=json.loads(independent_path.read_text())
    if independent.get('status')!='PASS_INDEPENDENT_SPARSE_REVIEW' or not all(independent['checks'].values()) or independent.get('source_sha256')!=sha(sparse_path):raise RuntimeError('Current sparse source must pass independent review before the joint solve')
    old=json.loads(FINE.read_text());coarse=json.loads(COARSE.read_text())
    if not old['converged'] or not coarse['converged']:raise RuntimeError('Input fields must be verified in their own integration rules')
    if (old['parameters']['nr'],old['parameters']['nz'],old['parameters']['bulk'],old['activation_strength'])!=(16,8,5000.,10.):raise RuntimeError('Preselected fine-case inputs changed')
    if (coarse['parameters']['nr'],coarse['parameters']['nz'],coarse['parameters']['bulk'],coarse['activation_strength'])!=(12,6,5000.,10.):raise RuntimeError('Preselected coarse-case inputs changed')
    params=Parameters(**old['parameters']);policy=IntegrationPolicy(8,12)
    sources={name:sha(ROOT/'scripts'/name) for name in ('run_joint_refinement.py','sparse_quadratic.py','quadratic_quadrature_study.py','robust_quadratic.py','axisymmetric_quadratic.py','axisymmetric_mixed.py','axisymmetric_prestress.py','axisymmetric_mechanics.py','run_activated_mesh_study.py','annulus_observer.py','mixed_geometry_certificate.py','mixed_geometry_certificate_v2.py')}
    protocol={'scope':'Target-free joint integration+mesh check, no gain refit; strongest default-discrepancy case preselected',
              'lambda':5000.,'strength':10.,'mesh':[16,8],'integration_policy':asdict(policy),'cross_policy':asdict(IntegrationPolicy(12,16)),
              'old_fine_input':FINE.relative_to(ROOT).as_posix(),'old_fine_sha256':sha(FINE),
              'coarse_highQ_input':COARSE.relative_to(ROOT).as_posix(),'coarse_highQ_sha256':sha(COARSE),
              'sparse_verification_sha256':sha(verification),'independent_sparse_review_sha256':sha(independent_path),'sources':sources,
              'initialization':'Old Q0 fine state is an initial guess only; fixed-load Q8 relaxation, no unverified previous-state admission',
              'comparison_width_um':'0.01 is reporting resolution, not biological/PDE acceptance','no_24x12_job':True}
    write(ROOT/'evidence/joint_refinement_protocol.json',protocol)
    write(live,{'status':'RUNNING','pid':os.getpid(),'UTC_start':datetime.now(timezone.utc).isoformat(),'sources':sources})
    try:
        model=AxisymmetricModel(params,policy=SolvePolicy(),integration_policy=policy,sparse_policy=SparsePolicy(seed_chunk=4))
        before=evaluate(model,old);started=time.monotonic()
        result=model._attempt(10.,np.array(old['free_displacement']))
        result['numerical_initial_state']={'path':FINE.relative_to(ROOT).as_posix(),'sha256':sha(FINE),'scope':protocol['initialization']}
        result['joint_refinement_source_hashes']=sources
        write(OUT/'final.json',result)
        crossmodel=AxisymmetricModel(params,policy=SolvePolicy(),integration_policy=IntegrationPolicy(12,16),sparse_policy=SparsePolicy(seed_chunk=4))
        cross=evaluate(crossmodel,result);geom=exact_geometry(result)
        checks={'actual_fixed_load_reached':result['converged'] and result['activation_strength']==10.,
                'higher_rule_cross_force':cross['scaled_force_linf']<1e-7,
                'exact_local_J':geom['all_exact_J_positive'],'exact_global_embedding':geom['global_embedding_certified'],
                'source_hashes_unchanged':all(sha(ROOT/'scripts'/name)==value for name,value in sources.items())}
        receipt={'status':'PASS_JOINT_DISCRETE_EQUILIBRIUM' if all(checks.values()) else 'FAIL','checks':checks,
                 'seconds':time.monotonic()-started,'protocol':protocol,'initial_guess_Q8_evaluation':before,
                 'final_Q8_force_linf':result['scaled_free_force_linf'],'cross_Q12_evaluation':cross,
                 'geometry':geom,'final_sha256':sha(OUT/'final.json'),
                 'scope':'Verified discrete integration/mesh experiment; no PDE-error or novelty claim.'}
        write(ROOT/'evidence/joint_refinement_solve.json',receipt)
        if not all(checks.values()):raise RuntimeError('Retained joint solve failed verification')
        comparison=compare_results(result,coarse)
        write(ROOT/'evidence/joint_refinement_comparison.json',{'status':'COMPLETE','comparison':comparison,
              'fine_sha256':sha(OUT/'final.json'),'coarse_sha256':sha(COARSE),'sources':sources,
              'scope':'Same higher integration policy, different P2 meshes. All comparison intervals concern discrete curves, not continuum/PDE uncertainty.'})
        write(live,{'status':'FINISHED','pid':os.getpid(),'sources':sources})
        print(json.dumps(serial({'status':receipt['status'],'seconds':receipt['seconds'],'force8':receipt['final_Q8_force_linf'],'force12':cross['scaled_force_linf'],
              'apical_vertical_sup_um':comparison['top_Eulerian_vertical']['supremum_absolute_discrete_height_difference_um'],
              'shortening_difference_um':comparison['shortening_difference_um']})),flush=True)
    except BaseException as error:
        write(live,{'status':'STOPPED_WITH_ERROR','pid':os.getpid(),'error':str(error),'sources':sources});raise

if __name__=='__main__':
    torch.set_num_threads(1);parser=argparse.ArgumentParser();parser.add_argument('--verification',type=Path,required=True);args=parser.parse_args();run(args.verification)
