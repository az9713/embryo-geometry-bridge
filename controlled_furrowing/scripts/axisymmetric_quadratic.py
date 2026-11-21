"""P2 axisymmetric full pointwise logarithmic-volume penalty.

Restores the original pointwise constitutive energy on the higher-order
geometry. Also supplies a target-independent forced manufactured equilibrium.
Neither frozen mixed nor original P1/prestress sources are modified.
"""
from dataclasses import dataclass,asdict
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
import torch
from axisymmetric_mixed import Parameters,AxisymmetricModel as MixedModel


class AxisymmetricModel(MixedModel):
    def energy_components(self,current,activation_strength=0.):
        result=super().energy_components(current,activation_strength)
        full=(self.integration*self.parameters.bulk/2*result['logJ']**2).sum()
        difference=full-result['projected_volumetric']
        result['total']=result['total']+difference
        result['bulk']=result['shear']+full
        result['pointwise_volumetric']=full
        result['unresolved_volumetric']=difference
        result['pointwise_pressure']=self.parameters.bulk*result['logJ']
        return result

    def diagnostics(self,vector,strength,history=None):
        result=super().diagnostics(vector,strength,history)
        parts=self.energy_components(self.coordinates(torch.tensor(vector,dtype=torch.float64)),strength)
        result['energy'].update({key:float(parts[key]) for key in ['pointwise_volumetric','unresolved_volumetric']})
        result['pointwise_pressure_at_quadrature']=parts['pointwise_pressure'].detach().numpy().tolist()
        result['pressure_interpretation']='Pointwise volumetric stress coefficient lambda*logJ. Continuous P1 projected coefficients/constraint are auxiliary diagnostics, not the governing constitutive law.'
        result['scope']='P2 full pointwise logarithmic-volume penalty, not a P1-pressure incompressibility cure; no target data or continuum error certificate.'
        return result


def build_model(parameters=None,**kwargs):
    return AxisymmetricModel(Parameters(**kwargs) if parameters is None else parameters)


class ManufacturedModel(AxisymmetricModel):
    """Exact nonlinear isochoric map with independently derived body force.

    All boundary displacement nodes use the exact map. P2 interpolation on
    those edges approximates the exact nonpolynomial Dirichlet boundary.
    """
    def __init__(self,parameters=Parameters(),epsilon=.2):
        if epsilon==0 or not np.isfinite(epsilon):raise ValueError('finite nonzero epsilon required')
        if any(getattr(parameters,key)!=0 for key in ['basal_area_modulus','apical_area_modulus','basal_surface_tension','apical_surface_tension']):raise ValueError('manufactured bulk benchmark requires zero surface energies')
        super().__init__(parameters);self.epsilon=epsilon
        reference=self.reference;p=parameters
        self.fixed=np.repeat(((reference[:,0]==0)|(reference[:,0]==p.outer_radius)|(reference[:,1]==0)|(reference[:,1]==p.height))[:,None],2,axis=1)
        self.free=np.flatnonzero(~self.fixed.ravel());self.free_t=torch.tensor(self.free,dtype=torch.long)
        self.bounds=[(-(1-1e-7)*reference.ravel()[i],None) if i%2==0 else (None,None) for i in self.free]
        exact=self.exact_coordinates(reference)
        prescribed=reference.copy();prescribed[self.fixed]=exact[self.fixed]
        self.prescribed=torch.tensor(prescribed,dtype=torch.float64)
        self.initial_interpolant=(exact-reference).ravel()[self.free]
        qx=np.einsum('qa,tad->tqd',self.bary,self.pressure_vertices[self.pressure_triangles])
        self.reference_quadrature=qx
        f=np.stack((-p.mu*epsilon**2*qx[:,:,0]*np.exp(epsilon*qx[:,:,1]),2*p.mu*epsilon*np.exp(-2*epsilon*qx[:,:,1])),axis=2)
        self.body_force=torch.tensor(f,dtype=torch.float64)

    def exact_coordinates(self,reference):
        y=np.asarray(reference).copy();e=self.epsilon
        y[...,0]=reference[...,0]*np.exp(e*reference[...,1])
        y[...,1]=-np.expm1(-2*e*reference[...,1])/(2*e)
        return y

    def coordinates(self,free_displacement):
        return self.prescribed+torch.zeros(self.X.numel(),dtype=torch.float64).scatter(0,self.free_t,free_displacement).reshape_as(self.X)

    def energy_components(self,current,activation_strength=0.):
        if activation_strength!=0:raise ValueError('manufactured body-force benchmark has no activation')
        result=super().energy_components(current,0.)
        yq=torch.einsum('qa,tad->tqd',self.shape,current[self.tri])
        force=-(self.integration*(self.body_force*yq).sum(dim=2)).sum()
        result['body_force_potential']=force;result['total']=result['total']+force
        return result

    def errors(self,vector):
        current=self.coordinates(torch.tensor(vector,dtype=torch.float64));parts=self.energy_components(current,0.)
        yq=torch.einsum('qa,tad->tqd',self.shape,current[self.tri]).detach().numpy()
        exact=self.exact_coordinates(self.reference_quadrature);q=self.reference_quadrature;e=self.epsilon
        fm=np.zeros_like(parts['F_meridional'].detach().numpy())
        fm[:,:,0,0]=np.exp(e*q[:,:,1]);fm[:,:,0,1]=e*q[:,:,0]*np.exp(e*q[:,:,1]);fm[:,:,1,1]=np.exp(-2*e*q[:,:,1])
        dF=parts['F_meridional'].detach().numpy()-fm
        dhoop=parts['hoop'].detach().numpy()-np.exp(e*q[:,:,1])
        weight=self.integration.numpy();J=parts['J'].detach().numpy()
        return {'displacement_L2_RMS':float(np.sqrt(np.sum(weight*np.sum((yq-exact)**2,axis=2))/self.reference_volume)),
                'gradient_H1_RMS':float(np.sqrt(np.sum(weight*(np.sum(dF*dF,axis=(2,3))+dhoop*dhoop))/self.reference_volume)),
                'J_L2_RMS':float(np.sqrt(np.sum(weight*(J-1)**2)/self.reference_volume)),
                'J_max_sampled_error':float(np.max(abs(J-1)))}

    def diagnostics(self,vector,strength,history=None):
        result=super().diagnostics(vector,strength,history)
        parts=self.energy_components(self.coordinates(torch.tensor(vector,dtype=torch.float64)),0.)
        result['energy']['body_force_potential']=float(parts['body_force_potential'])
        result['scope']='Target-independent forced manufactured nonlinear equilibrium; analytic nodal Dirichlet boundary, no activation/surface energy/biological data. Sampled errors are not rigorous PDE error bounds.'
        result['manufactured_epsilon']=self.epsilon
        return result


def _self_test():
    root=Path(__file__).resolve().parents[1];out=root/'runs/quadratic_checks';out.mkdir(parents=True,exist_ok=True)
    protocol={'scope':'Synthetic P2 pointwise penalty consistency','thresholds':{'projection_residual_identity':1e-11,'AD_relative':2e-5,'force':1e-7},'source_parents_frozen':True}
    (root/'evidence/quadratic_protocol.json').write_text(json.dumps(protocol,indent=2))
    m=build_model(Parameters(nr=6,nz=3,basal_surface_tension=2.,apical_surface_tension=.2))
    rng=np.random.default_rng(1911);u=rng.normal(size=len(m.free))*1e-4;d=rng.normal(size=len(u));d/=np.linalg.norm(d)
    current=m.coordinates(torch.tensor(u,dtype=torch.float64));c=m.energy_components(current,.2)
    coeff=c['pressure']/m.parameters.bulk;projected=torch.einsum('qi,ti->tq',m.pressure_shape,coeff[m.pressure_tri_t])
    residual=.5*m.parameters.bulk*(m.integration*(c['logJ']-projected)**2).sum()
    identity_error=abs(float(residual-c['unresolved_volumetric']))
    _,g=m.objective(u,.2);h=1e-6;fd=(m.objective(u+h*d,.2,False)-m.objective(u-h*d,.2,False))/(2*h);ad=float(g@d);derror=abs(fd-ad)/max(abs(fd),abs(ad),1e-10)
    zero=m.relax(0.);active=m.relax(.2)
    checks={'projection_kernel_energy_identity':identity_error<1e-11,'AD_vs_FD':derror<2e-5,'zero_reference':zero['converged'],'activated_stationary_positive':active['converged']}
    for name,result in [('zero',zero),('activated',active)]:
        (out/(name+'.json')).write_text(json.dumps(result,indent=2))
    receipt={'status':'PASS_SYNTHETIC_COMPONENT' if all(checks.values()) else 'FAIL','checks':{k:bool(v) for k,v in checks.items()},'identity_error':identity_error,'AD_relative_error':derror,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'parents':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ['axisymmetric_mixed.py','axisymmetric_prestress.py','axisymmetric_mechanics.py']}}
    (root/'evidence/quadratic_verification.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2));return all(checks.values())


def _manufactured_study():
    root=Path(__file__).resolve().parents[1];out=root/'runs/quadratic_checks/manufactured';out.mkdir(parents=True,exist_ok=True)
    protocol={'scope':'Forced nonlinear continuum equilibrium with lambda-independent exact map; sampled errors are not rigorous PDE error bounds','epsilon':.2,'mu':1.,'outer_radius':3.,'height':1.,'grids':[[6,3],[8,4],[12,6]],'lambda':[5.,50.,500.],
              'forcing':['f_r=-mu epsilon^2 r exp(epsilon z)','f_z=2mu epsilon exp(-2epsilon z)'],'boundary':'All-boundary displacement nodal values set to exact map; quadratic edge interpolation','initial':'Analytic nodal interpolant','gain':0.,'surface_energies':0.,'solver_maxiter':200,'thresholds':{'force_linf':1e-7,'positive_orientation':True},'predeclared_evaluation':['Allcase stationarity+geometry','L2/H1/J errors by mesh and lambda','Do not equate sampled refinement with a rigorous continuum error certificate']}
    (root/'evidence/quadratic_manufactured_protocol.json').write_text(json.dumps(protocol,indent=2));rows=[]
    for nr,nz in protocol['grids']:
        for lam in protocol['lambda']:
            start=time.monotonic();m=ManufacturedModel(Parameters(nr=nr,nz=nz,bulk=lam,basal_area_modulus=0.,apical_area_modulus=0.))
            before=m.errors(m.initial_interpolant);result=m.relax(0.,m.initial_interpolant,maxiter=200);after=m.errors(result['free_displacement'])
            path=out/f'mesh{nr}x{nz}_lambda{lam}.json';path.write_text(json.dumps(result,indent=2))
            row={'nr':nr,'nz':nz,'lambda':lam,'seconds':time.monotonic()-start,'file':str(path),'converged':result['converged'],'force_linf':result['scaled_free_force_linf'],'sampled_min_J':result['sampled_minimum_J'],'sampled_max_J':result['sampled_maximum_J'],'Bernstein_J_lower':result['minimum_J'],'Bernstein_J_upper':result['maximum_J'],'interpolation_errors':before,'solved_errors':after}
            rows.append(row);print(json.dumps(row),flush=True)
            receipt={'status':'COMPLETED_FORCED_NUMERICAL_BENCHMARK' if len(rows)==9 else 'RUNNING','complete':len(rows)==9,'all_stationary_positive_cert':all(q['converged'] for q in rows),'rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'protocol':protocol}
            (root/'evidence/quadratic_manufactured.json').write_text(json.dumps(receipt,indent=2))


if __name__=='__main__':
    torch.set_num_threads(1);parser=argparse.ArgumentParser();parser.add_argument('--manufactured',action='store_true');args=parser.parse_args()
    if args.manufactured:_manufactured_study()
    else:raise SystemExit(0 if _self_test() else 1)
