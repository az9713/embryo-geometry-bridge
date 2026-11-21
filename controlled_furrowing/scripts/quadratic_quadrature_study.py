"""Declared integration-policy study on the unchanged P2 finite element space.

No global monkeypatches or edits to frozen model/controller modules. Rebuilds
quadrature caches and evaluates the same full pointwise constitutive energy.
"""
from dataclasses import dataclass,asdict
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
import torch
from robust_quadratic import AxisymmetricModel as RobustModel,SolvePolicy
from axisymmetric_quadratic import Parameters,AxisymmetricModel as DefaultModel
from axisymmetric_mixed import _shape
from mixed_geometry_certificate import certify_flat_rectangular_boundary
from mixed_geometry_certificate_v2 import correlated_J_certificate,strong_monotonicity_certificate

ROOT=Path(__file__).resolve().parents[1]

@dataclass(frozen=True)
class IntegrationPolicy:
    volume_duffy_gauss_order:int=8
    surface_gauss_order:int=12
    volume_chunk_triangles:int=12
    def __post_init__(self):
        if self.volume_duffy_gauss_order not in (4,6,8,12,16):raise ValueError('declared volume orders are4,6,8,12,16')
        if self.surface_gauss_order not in (8,12,16,20):raise ValueError('declared surface orders are8,12,16,20')
        if self.volume_chunk_triangles<1:raise ValueError('positive chunk size required')

def duffy_gauss(order):
    q,w=np.polynomial.legendre.leggauss(order);q=(q+1)/2;w=w/2
    bary=np.array([[1-u,u*(1-v),u*v] for u in q for v in q])
    weights=np.array([2*w[i]*w[j]*u for i,u in enumerate(q) for j,v in enumerate(q)])
    return bary,weights

class AxisymmetricModel(RobustModel):
    def __init__(self,parameters=Parameters(),policy=SolvePolicy(),integration_policy=IntegrationPolicy()):
        self.integration_policy=integration_policy
        super().__init__(parameters,policy)
        p=parameters;tri=self.pressure_triangles;vertices=self.pressure_vertices
        X=vertices[tri];dm=np.stack((X[:,1]-X[:,0],X[:,2]-X[:,0]),axis=2)
        grad_lambda=np.array([[-1.,-1.],[1.,0.],[0.,1.]])[None]@np.linalg.inv(dm)
        self.bary,self.weights=duffy_gauss(integration_policy.volume_duffy_gauss_order)
        N,grad=_shape(self.bary,grad_lambda)
        self.shape=torch.tensor(N,dtype=torch.float64);self.shape_grad=torch.tensor(grad,dtype=torch.float64)
        qr=np.einsum('qa,tad->tqd',self.bary,X)[:,:,0]
        self.qr=torch.tensor(qr,dtype=torch.float64)
        weight=2*np.pi*self.area[:,None]*qr*self.weights[None,:]
        self.integration=torch.tensor(weight,dtype=torch.float64)
        self.reference_volume=float(weight.sum());self.energy_scale=p.mu*self.reference_volume
        self.pressure_shape=torch.tensor(self.bary,dtype=torch.float64)
        mass=np.zeros((len(vertices),len(vertices)))
        for i in range(3):
            for j in range(3):np.add.at(mass,(tri[:,i],tri[:,j]),np.sum(weight*self.bary[None,:,i]*self.bary[None,:,j],axis=1))
        self.mass=torch.tensor(mass,dtype=torch.float64);self.mass_inverse=torch.tensor(np.linalg.inv(mass),dtype=torch.float64)
        self.surfaces=[]
        q,w=np.polynomial.legendre.leggauss(integration_policy.surface_gauss_order);q=(q+1)/2;w=w/2
        N=np.column_stack(((1-q)*(1-2*q),4*q*(1-q),q*(2*q-1)))
        dN=np.column_stack((-3+4*q,4-8*q,-1+4*q))
        for edges,modulus,tension in ((self.bottom_edges,p.basal_area_modulus,p.basal_surface_tension),(self.top_edges,p.apical_area_modulus,p.apical_surface_tension)):
            xy=self.reference[edges];r0=np.einsum('qa,tad->tqd',N,xy)[:,:,0];length=np.linalg.norm(xy[:,2]-xy[:,0],axis=1)
            self.surfaces.append(tuple(torch.tensor(a,dtype=torch.long if k==0 else torch.float64) for k,a in enumerate((edges,N,dN,r0,length,2*np.pi*r0*length[:,None]*w)))+(modulus,tension))

    def energy_components(self,current,activation_strength=0.):
        # Small element chunks limit vectorized Hessian peak memory. Summing
        # these chunks changes no finite-dimensional mechanical potential.
        p=self.parameters;zero=torch.zeros((),dtype=torch.float64)
        shear=zero;full=zero;active=zero;invalid=zero;F=[];JJ=[];hh=[];ll=[];dd=[]
        chunk=self.integration_policy.volume_chunk_triangles
        for start in range(0,len(self.triangles),chunk):
            sl=slice(start,start+chunk);y=current[self.tri[sl]]
            fm=torch.einsum('tad,tqaj->tqdj',y,self.shape_grad[sl]);det=torch.linalg.det(fm)
            radius=torch.einsum('qa,tad->tqd',self.shape,y)[:,:,0];hoop=radius/self.qr[sl]
            J=det*hoop;logJ=torch.log(torch.clamp(J,min=1e-8));w=self.integration[sl]
            tr=(fm*fm).sum(dim=(2,3))+hoop**2
            shear=shear+(w*p.mu/2*(tr-3-2*logJ)).sum()
            full=full+(w*p.bulk/2*logJ**2).sum()
            active=active+(w*self.active[sl]*activation_strength/2*(fm[:,:,:,1]**2).sum(dim=2)).sum()
            invalid=invalid+(1e8*w*torch.relu(1e-6-J)**2).sum()
            F.append(fm);JJ.append(J);hh.append(hoop);ll.append(logJ);dd.append(det)
        logJ=torch.cat(ll);J=torch.cat(JJ)
        # The projected pressure is diagnostic only in the full-penalty law.
        # Detaching this auxiliary calculation avoids traversing a cancelled
        # projection graph when forming Hessians of the physical potential.
        local_b=torch.einsum('tq,tq,qi->ti',self.integration,logJ.detach(),self.pressure_shape)
        b=torch.zeros(self.pressure_count,dtype=torch.float64).index_add(0,self.pressure_tri_t.ravel(),local_b.ravel())
        pressure=p.bulk*(self.mass_inverse@b);projected=.5*(b@pressure)
        surface=zero;tension_energy=zero
        for edges,N,dN,r0,l0,w,modulus,tension in self.surfaces:
            y=current[edges];R=torch.einsum('qa,tad->tqd',N,y)[:,:,0]
            tangent=torch.einsum('qa,tad->tqd',dN,y);Js=R/r0*torch.linalg.vector_norm(tangent,dim=2)/l0[:,None]
            surface=surface+(w*modulus/2*(Js-1)**2).sum();tension_energy=tension_energy+(w*tension*Js).sum()
        return {'bulk':shear+full,'shear':shear,'pointwise_volumetric':full,'projected_volumetric':projected,
                'unresolved_volumetric':full-projected,'surface':surface,'tension':tension_energy,'active':active,
                'trial_inversion_penalty':invalid,'total':shear+full+surface+tension_energy+active+invalid,
                'J':J,'logJ':logJ,'F_meridional':torch.cat(F),'meridional_determinant':torch.cat(dd),
                'hoop':torch.cat(hh),'pressure':pressure,'pressure_rhs':b,'pointwise_pressure':p.bulk*logJ}

    def diagnostics(self,vector,strength,history=None):
        answer=super().diagnostics(vector,strength,history)
        answer['integration_policy']=asdict(self.integration_policy)
        answer['nominal_parent_quadrature_fields']='Parent parameter quadrature=4/surface_quadrature=8 describe frozen construction only; integration_policy above is the actual rule.'
        answer['scope']='Same P2 finite-dimensional field space and full pointwise constitutive law, independently declared higher numerical integration policy; no continuum error certificate or target data.'
        return answer

def source_hashes():
    names=('quadratic_quadrature_study.py','robust_quadratic.py','axisymmetric_quadratic.py','axisymmetric_mixed.py','axisymmetric_prestress.py','axisymmetric_mechanics.py')
    return {name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}

def evaluate(model,state):
    e,g=model.objective(state['free_displacement'],state['activation_strength'])
    parts=model.energy_components(model.coordinates(torch.tensor(state['free_displacement'],dtype=torch.float64)),state['activation_strength'])
    return {'scaled_energy':e,'energy':e*model.energy_scale,'scaled_force_linf':float(np.max(abs(g))),
            'energy_components':{k:float(parts[k]) for k in ('bulk','surface','tension','active','pointwise_volumetric','total')},
            'sampled_J_min':float(parts['J'].min()),'sampled_J_max':float(parts['J'].max())}

def run_evaluation():
    output=ROOT/'runs/activated_quadrature';output.mkdir(parents=True,exist_ok=True)
    protocol={'scope':'Target-free same-space integration sensitivity on actual activated small geometry','mesh':[12,6],
              'inputs':'runs/quadratic_robust/lambda_{500.0,5000.0}_target_{1.0,5.0,10.0}/final.json',
              'policies':[asdict(IntegrationPolicy(8,12)),asdict(IntegrationPolicy(12,16))],
              'force_threshold':1e-7,'source_hashes':source_hashes(),
              're_equilibration_if_needed':'Warm fixed-strength solve under the new rule, retaining the original field; not a physical activation trajectory.'}
    (ROOT/'evidence/activated_quadrature_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    rows=[]
    for lam in (500.,5000.):
        for strength in (1.,5.,10.):
            path=ROOT/'runs/quadratic_robust'/f'lambda_{lam}_target_{strength}'/'final.json';state=json.loads(path.read_text())
            params=Parameters(**state['parameters']);default=evaluate(DefaultModel(params),state)
            policies=[]
            for policy in (IntegrationPolicy(8,12),IntegrationPolicy(12,16)):
                model=AxisymmetricModel(params,integration_policy=policy);result=evaluate(model,state)
                policies.append({'policy':asdict(policy),**result,'force_gate':result['scaled_force_linf']<1e-7,'energy_difference_from_default':result['energy']-default['energy']})
            row={'lambda':lam,'strength':strength,'input':path.relative_to(ROOT).as_posix(),
                 'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'default_evaluation':default,'higher_evaluations':policies}
            rows.append(row);print(json.dumps(row),flush=True)
            (ROOT/'evidence/activated_quadrature_evaluation.json').write_text(json.dumps({'status':'COMPLETE' if len(rows)==6 else 'RUNNING','rows':rows,'protocol':protocol},indent=2)+'\n')

def self_test():
    p=Parameters(nr=6,nz=3,basal_surface_tension=2.,apical_surface_tension=0.)
    model=AxisymmetricModel(p);checks={}
    _,gradient=model.objective(np.zeros(len(model.free)),0.)
    checks['zero_reference_force']=float(np.max(abs(gradient)))<1e-12
    lam=1.1;y=model.X.clone();y[:,0]*=lam;y[:,1]/=lam**2
    actual=float(model.energy_components(y,0.)['total'])
    exact=model.reference_volume*p.mu/2*(2*lam**2+lam**-4-3)+np.pi*p.outer_radius**2*((p.basal_area_modulus+p.apical_area_modulus)/2*(lam**2-1)**2+(p.basal_surface_tension+p.apical_surface_tension)*lam**2)
    checks['affine_isochoric_energy']=abs(actual-exact)<1e-10
    checks['positive_quadrature_weights']=bool(np.all(model.integration.numpy()>0))
    checks['reference_volume']=abs(model.reference_volume-np.pi*p.outer_radius**2*p.height)<1e-11
    checks['pressure_mass_positive_definite']=float(torch.linalg.eigvalsh(model.mass).min())>0
    rng=np.random.default_rng(91719);u=rng.normal(size=len(model.free))*1e-4;direction=rng.normal(size=len(u));direction/=np.linalg.norm(direction)
    e,g=model.objective(u,5.);h=1e-6;fd=(model.objective(u+h*direction,5.,False)-model.objective(u-h*direction,5.,False))/(2*h)
    rel=abs(float(g@direction)-fd)/max(abs(fd),abs(float(g@direction)),1e-12)
    checks['AD_vs_independent_directional_FD']=rel<2e-5
    other=AxisymmetricModel(p,integration_policy=IntegrationPolicy(8,12,1));ee,gg=other.objective(u,5.)
    checks['chunk_size_potential_invariant']=abs(e-ee)<1e-12
    checks['chunk_size_gradient_invariant']=float(np.max(abs(g-gg)))<1e-11
    receipt={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'AD_FD_relative_error':rel,'source_hashes':source_hashes(),
             'scope':'Synthetic cache/integration/derivative checks; no continuum error certificate or target data.'}
    (ROOT/'evidence/activated_quadrature_synthetic.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)
    if not all(checks.values()):raise RuntimeError('quadrature component synthetic verification failed')

def exact_geometry(result):
    from fractions import Fraction as Q
    n=max(i for tri in result['triangles'] for i in tri[:3])+1
    certs=[correlated_J_certificate([result['reference_nodes'][i] for i in tri[:3]],[result['current_nodes'][i] for i in tri]) for tri in result['triangles']]
    lower=min(Q(c['J_bound']['lower']['exact']) for c in certs);upper=max(Q(c['J_bound']['upper']['exact']) for c in certs)
    boundary=certify_flat_rectangular_boundary(result['reference_nodes'][:n],result['current_nodes'],result['triangles'])
    monotone=strong_monotonicity_certificate(result['reference_nodes'][:n],result['current_nodes'],result['triangles'])
    return {'all_exact_J_positive':all(c['positive_J_certified'] for c in certs),'J_lower_exact':str(lower),'J_upper_exact':str(upper),
            'J_lower_display':float(lower),'J_upper_display':float(upper),'triangle_certificates':certs,'boundary':boundary,'strong_monotonicity':monotone,
            'global_embedding_certified':boundary['status']=='CERTIFIED_JORDAN_BOUNDARY' or monotone['global_injectivity_certified']}

def run_reequilibration():
    self_test();run_evaluation();rows=[]
    for lam in (500.,5000.):
        for strength in (1.,5.,10.):
            path=ROOT/'runs/quadratic_robust'/f'lambda_{lam}_target_{strength}'/'final.json';old=json.loads(path.read_text());p=Parameters(**old['parameters'])
            output=ROOT/'runs/activated_quadrature'/f'lambda_{lam}_strength_{strength}';output.mkdir(parents=True,exist_ok=True)
            policy=IntegrationPolicy(8,12);model=AxisymmetricModel(p,integration_policy=policy)
            started=time.monotonic();result=model._attempt(strength,np.array(old['free_displacement']))
            result['numerical_initial_state']={'path':path.relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                'scope':'Same-strength relaxation of old-rule field under new quadrature; not a physical activation path.'}
            result['study_source_hashes']=source_hashes();resultpath=output/'reequilibrated.json';resultpath.write_text(json.dumps(result,indent=2)+'\n')
            high=AxisymmetricModel(p,integration_policy=IntegrationPolicy(12,16));cross=evaluate(high,result)
            oldvec=np.array(old['free_displacement']);newvec=np.array(result['free_displacement'])
            # Exact P2 boundary curves have the same parameter intervals on
            # this same mesh. Bernstein control differences bound the vector
            # displacement difference over every common material interval.
            curve_bounds={}
            for name,edges in [('apical',model.top_edges),('basal',model.bottom_edges)]:
                delta=np.array(result['current_nodes'])[edges]-np.array(old['current_nodes'])[edges]
                bezier=np.stack((delta[:,0],2*delta[:,1]-.5*(delta[:,0]+delta[:,2]),delta[:,2]),axis=1)
                curve_bounds[name]={'material_parameter_curve_displacement_bound':float(np.linalg.norm(bezier,axis=2).max()),
                    'note':'Same-material-reference parameter difference; not an Eulerian graph difference or rigorous floating-point enclosure.'}
            geom=exact_geometry(result)
            row={'lambda':lam,'strength':strength,'seconds':time.monotonic()-started,'file':resultpath.relative_to(ROOT).as_posix(),
                'file_sha256':hashlib.sha256(resultpath.read_bytes()).hexdigest(),'integration_policy':asdict(policy),
                'converged':result['converged'],'re_equilibrated_force_linf':result['scaled_free_force_linf'],
                'cross_evaluation_Duffy12_surface16':cross,'cross_force_gate':cross['scaled_force_linf']<1e-7,
                'shortening_before':old['maximum_activated_column_shortening'],'shortening_after':result['maximum_activated_column_shortening'],
                'shortening_change':result['maximum_activated_column_shortening']-old['maximum_activated_column_shortening'],
                'free_vector_linf_change':float(np.max(abs(newvec-oldvec))),'curve_differences':curve_bounds,
                'exact_geometry':geom}
            rows.append(row)
            (ROOT/'evidence/activated_quadrature_reequilibration.json').write_text(json.dumps({'status':'COMPLETE' if len(rows)==6 else 'RUNNING',
                'rows':rows,'source_hashes':source_hashes(),'scope':'Same-space quadrature sensitivity; no mesh/PDE error certificate or target data.'},indent=2)+'\n')
            print(json.dumps({k:v for k,v in row.items() if k not in ('exact_geometry','curve_differences')}),flush=True)
            if not result['converged'] or not row['cross_force_gate'] or not geom['all_exact_J_positive'] or not geom['global_embedding_certified']:
                raise RuntimeError('retained reequilibrated state failed a declared readiness check')

if __name__=='__main__':
    torch.set_num_threads(1);parser=argparse.ArgumentParser();parser.add_argument('--evaluate',action='store_true');parser.add_argument('--solve',action='store_true');parser.add_argument('--self-test',action='store_true');args=parser.parse_args()
    if args.evaluate:run_evaluation()
    if args.self_test:self_test()
    if args.solve:run_reequilibration()
