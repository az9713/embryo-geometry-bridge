"""Matched16x8 activated-body study with frozen energies and numerical policy.

No response data or fitting. Physical Eulerian graph comparisons use rational
quadratic inversion; their enclosures concern two discrete curves, not PDE error.
"""
from pathlib import Path
from fractions import Fraction as Q
from dataclasses import asdict
from datetime import datetime,timezone
import argparse,hashlib,heapq,json,os,platform,time
import numpy as np
import torch
from robust_quadratic import Parameters,SolvePolicy,build_model
from annulus_observer import _quadratic_controls,_root_bracket,_poly,Interval,intersection,serial

ROOT=Path(__file__).resolve().parents[1]
FROZEN={
 'axisymmetric_mechanics.py':'ca431e30a2574a02f4a21c98935c50217aa79bc5512cf7279143f22ce1a580e0',
 'axisymmetric_prestress.py':'141b4b2ea92f4e04f4076b18edc9bafd6218f75a7c4c9633e66c5d2a9fc03c41',
 'axisymmetric_mixed.py':'daf9b4f52181b18b1a75cc69909444bf2a85725c3896b32b125b2f20f10b1f14',
 'axisymmetric_quadratic.py':'5bb64d6b441f6db8ada9da0f1bcb44716c86ed5d816e44f75bd4e9e021fdf410',
 'robust_quadratic.py':'be6f17e620ddf04533d456fff9d48cd48517f4a0709d921ab510616926b4176b',
 'annulus_observer.py':'77bcd1bc10e557de4a4be9c304041f9f0dfe0fcbac47871cabfbb11ffb7346b9'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,value):Path(path).write_text(json.dumps(serial(value),indent=2)+'\n')


class QuadraticGraph:
    def __init__(self,controls,height_um=Q('33.5')):
        physical=[[(Q(r)*height_um,Q(z)*height_um) for r,z in segment] for segment in controls]
        self.curves=_quadratic_controls(physical)
        self.domain=(self.curves[0][0][0][0],self.curves[-1][0][2][0])

    def height_range(self,lo,hi):
        radius=Interval(lo,hi)
        if lo<self.domain[0] or hi>self.domain[1]:raise ValueError('common physical radius interval not covered')
        values=[]
        for q,c0,c1,c2 in self.curves:
            cut=intersection(radius,Interval(q[0][0],q[2][0]))
            if cut is None:continue
            left=_root_bracket(c0[0],c1[0],c2[0],cut.lo,64)[0]
            right=_root_bracket(c0[0],c1[0],c2[0],cut.hi,64)[1]
            values += [_poly(c0[1],c1[1],c2[1],left),_poly(c0[1],c1[1],c2[1],right)]
            if c2[1]!=0:
                critical=-c1[1]/(2*c2[1])
                if left<=critical<=right:values.append(_poly(c0[1],c1[1],c2[1],critical))
        return min(values),max(values)

    def knots(self):
        points=[]
        for q,c0,c1,c2 in self.curves:
            points += [q[0][0],q[2][0]]
            if c2[1]!=0:
                t=-c1[1]/(2*c2[1])
                if 0<t<1:points.append(_poly(c0[0],c1[0],c2[0],t))
        return points


def point_difference(a,b,radius):
    alo,ahi=a.height_range(radius,radius);blo,bhi=b.height_range(radius,radius)
    return alo-bhi,ahi-blo


def compare_graphs(a,b,tolerance_um=Q('0.01'),maximum_cells=4096):
    """Rational discrete-curve supremum enclosure using radial interval cells."""
    lo=max(a.domain[0],b.domain[0]);hi=min(a.domain[1],b.domain[1])
    if lo>=hi:raise ValueError('nonempty common physical radial domain required')
    knots=sorted(set([lo,hi]+[x for x in a.knots()+b.knots() if lo<=x<=hi]+[lo+(hi-lo)*Q(i,128) for i in range(129)]))
    lower=Q(0);counter=0;heap=[]
    def cell(left,right):
        nonlocal lower,counter
        al,ah=a.height_range(left,right);bl,bh=b.height_range(left,right)
        upper=max(abs(al-bh),abs(ah-bl))
        for r in [left,(left+right)/2,right]:
            dl,dh=point_difference(a,b,r)
            lower=max(lower,dl if dl>0 else (-dh if dh<0 else Q(0)))
        counter+=1;return (-upper,counter,left,right)
    for left,right in zip(knots[:-1],knots[1:]):heapq.heappush(heap,cell(left,right))
    while -heap[0][0]-lower>tolerance_um and len(heap)<maximum_cells:
        _,_,left,right=heapq.heappop(heap);mid=(left+right)/2
        heapq.heappush(heap,cell(left,mid));heapq.heappush(heap,cell(mid,right))
    upper=-heap[0][0];sampled=[]
    for i in range(257):
        r=lo+(hi-lo)*Q(i,256);dl,dh=point_difference(a,b,r)
        al,ah=a.height_range(r,r);bl,bh=b.height_range(r,r)
        sampled.append({'radius_um':float(r),'fine_z_um':float((al+ah)/2),'coarse_z_um':float((bl+bh)/2),'difference_lower_um':float(dl),'difference_upper_um':float(dh)})
    return {'common_physical_radius_domain_um':(lo,hi),'supremum_absolute_discrete_height_difference_um':Interval(lower,upper),
            'requested_enclosure_width_um':tolerance_um,'enclosure_width_um':upper-lower,'requested_width_achieved':upper-lower<=tolerance_um,
            'active_radial_cells':len(heap),'maximum_cells':maximum_cells,'root_bisection_steps':64,'sampled_eulerian_curves':sampled,
            'interpretation':'Exact rational enclosure for supplied stored discrete P2 curves after declared33.5um conversion; NOT a continuum PDE error bound or physical uncertainty interval.'}


def _geometry_test():
    same1=QuadraticGraph([[(0,0),(Q(1,2),Q(1,2)),(1,1)]],Q(1))
    same2=QuadraticGraph([[(0,0),(Q(1,4),Q(1,4)),(1,1)]],Q(1))
    lo,hi=point_difference(same1,same2,Q(1,4))
    constant1=QuadraticGraph([[(0,1),(Q(1,2),1),(1,1)]],Q(1))
    constant2=QuadraticGraph([[(0,Q(11,10)),(Q(1,2),Q(11,10)),(1,Q(11,10))]],Q(1))
    result=compare_graphs(constant2,constant1)
    interval=result['supremum_absolute_discrete_height_difference_um']
    checks={'same_Eulerian_graph_different_material_parameterization':lo<=0<=hi and hi-lo<Q(1,10**15),
            'exact_constant_offset_supremum':interval.lo==interval.hi==Q(1,10),
            'full_common_physical_domain':result['common_physical_radius_domain_um']==(Q(0),Q(1))}
    receipt={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'source_sha256':sha(__file__),'scope':'Synthetic exactgraphcomparison only; no mechanical solve or target data.'}
    write(ROOT/'evidence/activated_mesh_graph_tests.json',receipt);print(json.dumps(receipt));return all(checks.values())


def compare_saved():
    comparisons=[]
    for lam in [500.,5000.]:
        for target in [1.,5.,10.]:
            fine_path=ROOT/f'runs/activated_mesh/lambda_{lam}_target_{target}/final.json'
            coarse_path=ROOT/f'runs/quadratic_robust/lambda_{lam}_target_{target}/final.json'
            if not fine_path.exists():continue
            fine=json.loads(fine_path.read_text());coarse=json.loads(coarse_path.read_text())
            if not fine['converged'] or not fine['reached_target']:continue
            if fine['activation_strength']!=target or coarse['activation_strength']!=target:raise ValueError('gain mismatch; no refitting allowed')
            pc=dict(coarse['parameters']);pf=dict(fine['parameters'])
            for key in ['nr','nz']:pc.pop(key);pf.pop(key)
            if pc!=pf:raise ValueError('constitutive/geometry/quadrature mismatch beyond mesh dimensions')
            row={'lambda':lam,'strength':target,'fine_field_sha256':sha(fine_path),'coarse_field_sha256':sha(coarse_path),'coarse_mesh':[12,6],'fine_mesh':[16,8],
                 'central_apical_displacement_difference_um':(fine['central_apical_displacement']-coarse['central_apical_displacement'])*33.5,
                 'central_basal_displacement_difference_um':(fine['central_basal_displacement']-coarse['central_basal_displacement'])*33.5,
                 'material_column_max_shortening_coarse_um':coarse['maximum_activated_column_shortening']*33.5,
                 'material_column_max_shortening_fine_um':fine['maximum_activated_column_shortening']*33.5,
                 'material_column_max_shortening_difference_um':(fine['maximum_activated_column_shortening']-coarse['maximum_activated_column_shortening'])*33.5}
            for surface in ['top','bottom']:
                controls='boundary_'+surface+'_quadratic_controls'
                row[surface+'_Eulerian_comparison']=compare_graphs(QuadraticGraph(fine[controls]),QuadraticGraph(coarse[controls]))
            comparisons.append(row);write(ROOT/'evidence/activated_mesh_comparison.json',{'status':'COMPLETE' if len(comparisons)==6 else 'PARTIAL','comparisons':comparisons,'source_sha256':sha(__file__),'scope':'Matched discrete mesh changes without refit; quadrature unchanged; no PDE error bound.'})
            print(json.dumps({'comparison_done':{'lambda':lam,'strength':target,'shortening_difference_um':row['material_column_max_shortening_difference_um'],
                             'apical_curve_difference_um':[float(row['top_Eulerian_comparison']['supremum_absolute_discrete_height_difference_um'].lo),float(row['top_Eulerian_comparison']['supremum_absolute_discrete_height_difference_um'].hi)]}}),flush=True)
    return comparisons


def run(resume=False):
    out=ROOT/'runs/activated_mesh';out.mkdir(parents=True,exist_ok=True)
    live=out/'active_job.json'
    if live.exists() and json.loads(live.read_text()).get('status')=='RUNNING':raise RuntimeError('An active job is already recorded. Use its live handle; do not launch a duplicate.')
    hashes={name:sha(ROOT/'scripts'/name) for name in FROZEN}
    if hashes!=FROZEN:raise RuntimeError('A frozen dependency changed; stop before solving')
    parameters={'nr':16,'nz':8,'height':1.,'outer_radius':3.,'patch_radius':9.35512037/33.5,'activation_z_min':8.5/33.5,'activation_z_max':28.5/33.5,'mu':1.,'basal_area_modulus':2.,'apical_area_modulus':.2,'basal_surface_tension':2.,'apical_surface_tension':0.}
    protocol={'scope':'Target-free matched activated mesh refinement, no response reads or gain refitting','physical_height_um':33.5,'fine_parameters_except_lambda':parameters,'coarse_states':'runs/quadratic_robust/lambda_*_target_*/final.json','lambda':[500.,5000.],'targets':[1.,5.,10.],'policy':asdict(SolvePolicy()),'integration':'Frozen6pointdegree4Dunavantvolume/surfaceGauss8; unchanged between meshes, separate from quadraturestudy','comparison':'CommonEulerianphysicalradii; exactrationalquadraticinverse/enclosures; tolerance0.01um is reporting resolution, not scientific acceptance or PDEerror','source_sha256':sha(__file__),'frozen_sources':hashes}
    write(ROOT/'evidence/activated_mesh_protocol.json',protocol)
    write(out/'environment.json',{'python':platform.python_version(),'platform':platform.platform(),'numpy':np.__version__,'torch':torch.__version__,'pid':os.getpid(),'UTC_start':datetime.now(timezone.utc).isoformat()})
    write(live,{'status':'RUNNING','pid':os.getpid(),'source_sha256':sha(__file__)})
    rows=[]
    try:
        for lam in protocol['lambda']:
            model=build_model(Parameters(bulk=lam,**parameters));previous=None
            for target in protocol['targets']:
                directory=out/f'lambda_{lam}_target_{target}';cached=directory/'final.json'
                if cached.exists():
                    if not resume:raise RuntimeError('Existing target output; use --resume only after confirming prior job stopped')
                    result=json.loads(cached.read_text());model._verify_state(result)
                    if not result.get('reached_target') or result['activation_strength']!=target:raise RuntimeError('Cached target failed; preserve and audit before retrying')
                else:
                    result=model.advance(target,previous,output_directory=directory)
                row={'lambda':lam,'strength':target,'reached':result['reached_target'],'verified_strength':result['activation_strength'],'force_linf':result['scaled_free_force_linf'],
                     'sampled_J_min':result['sampled_minimum_J'],'sampled_J_max':result['sampled_maximum_J'],'Bernstein_J_lower':result['minimum_J'],'Bernstein_J_upper':result['maximum_J'],
                     'geometry_certified':result['orientation_certificate']['certified_positive'],'relative_volume_change':result['relative_volume_change'],
                     'central_apical_displacement_um':result['central_apical_displacement']*33.5,'central_basal_displacement_um':result['central_basal_displacement']*33.5,
                     'material_column_max_shortening_um':result['maximum_activated_column_shortening']*33.5,'attempts':len(result['robust_continuation']['attempts']),
                     'rejected_attempts':sum(not a['accepted'] for a in result['robust_continuation']['attempts']),'seconds':result['robust_continuation']['seconds'],'field_sha256':sha(cached)}
                rows.append(row);write(ROOT/'evidence/activated_mesh_results.json',{'status':'COMPLETE_VERIFIED_TARGETS' if len(rows)==6 and all(r['reached'] for r in rows) else 'PARTIAL_OR_RUNNING','rows':rows,'source_sha256':sha(__file__),'protocol':protocol})
                print(json.dumps(row),flush=True)
                if not result['reached_target']:break
                previous=result
        compare_saved()
        write(live,{'status':'FINISHED','pid':os.getpid(),'targets_completed':len(rows),'all_reached':len(rows)==6 and all(r['reached'] for r in rows),'source_sha256':sha(__file__)})
    except BaseException as error:
        write(live,{'status':'STOPPED_WITH_ERROR','pid':os.getpid(),'error':str(error),'source_sha256':sha(__file__)})
        raise


if __name__=='__main__':
    torch.set_num_threads(1);parser=argparse.ArgumentParser();parser.add_argument('--test-graphs',action='store_true');parser.add_argument('--compare-only',action='store_true');parser.add_argument('--resume',action='store_true');args=parser.parse_args()
    if args.test_graphs:raise SystemExit(0 if _geometry_test() else 1)
    elif args.compare_only:compare_saved()
    else:run(args.resume)
