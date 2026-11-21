"""Flat-reference axisymmetric solid with physical surface area tensions.

Extends the frozen axisymmetric_mechanics implementation without modifying it.
Uniform tension on the planar surfaces is balanced by the clamped rim. No
experimental outcomes or specimen data are read by this module.
"""
from dataclasses import dataclass, asdict, replace
from pathlib import Path
import hashlib, json, time
import numpy as np
import torch
from axisymmetric_mechanics import Parameters as BaseParameters
from axisymmetric_mechanics import AxisymmetricModel as BaseModel

FROZEN_BASE_SHA256='ca431e30a2574a02f4a21c98935c50217aa79bc5512cf7279143f22ce1a580e0'


@dataclass(frozen=True)
class Parameters(BaseParameters):
    basal_surface_tension: float = 0.
    apical_surface_tension: float = 0.

    def __post_init__(self):
        super().__post_init__()
        if self.basal_surface_tension<0 or self.apical_surface_tension<0:
            raise ValueError('surface tensions must be finite and nonnegative')
        if self.reference_curvature!=0:
            raise ValueError('prestress extension requires flat reference; curved pressure/load balance is not implemented')


class AxisymmetricModel(BaseModel):
    def __init__(self, parameters=Parameters()):
        if not isinstance(parameters,Parameters):
            raise TypeError('use axisymmetric_prestress.Parameters')
        super().__init__(parameters)

    def tension_components(self,current):
        energies=[]
        tensions=(self.parameters.basal_surface_tension,self.parameters.apical_surface_tension)
        for (edges,l0,r0,q,weight,_),tension in zip(self.surfaces,tensions):
            yc=current[edges]
            length=torch.linalg.vector_norm(yc[:,1]-yc[:,0],dim=1)
            rq=(1-q[None,:])*yc[:,0,0,None]+q[None,:]*yc[:,1,0,None]
            Js=(rq/r0)*(length/l0)[:,None]
            energies.append((weight*tension*Js).sum())
        return {'basal_tension':energies[0],'apical_tension':energies[1]}

    def energy_components(self,current,activation_strength=0.):
        result=super().energy_components(current,activation_strength)
        tension=self.tension_components(current)
        result.update(tension)
        result['total']=result['total']+tension['basal_tension']+tension['apical_tension']
        return result

    def diagnostics(self,vector,strength,history=None):
        result=super().diagnostics(vector,strength,history)
        tensions=self.tension_components(self.coordinates(torch.tensor(vector,dtype=torch.float64)))
        result['energy'].update({name:float(value) for name,value in tensions.items()})
        result['reference_tension_energy']=float(np.pi*self.parameters.outer_radius**2*(self.parameters.basal_surface_tension+self.parameters.apical_surface_tension))
        result['scope']='Flat idealized finite-height axisymmetric quasistatic equilibrium with physical surface tension; no target data, no transient biological prediction.'
        return result


def build_model(parameters=None,**kwargs):
    return AxisymmetricModel(Parameters(**kwargs) if parameters is None else parameters)


def _self_test():
    root=Path(__file__).resolve().parents[1];out=root/'runs/prestress_checks';out.mkdir(parents=True,exist_ok=True)
    evidence=root/'evidence';evidence.mkdir(exist_ok=True)
    source=Path(__file__);base=source.with_name('axisymmetric_mechanics.py')
    protocol={'scope':'Synthetic flat-reference surface-tension extension, no outcome reads',
              'parameters':asdict(Parameters(nr=8,nz=4,basal_surface_tension=2.,apical_surface_tension=.2)),
              'thresholds':{'force_linf':1e-7,'zero_displacement':1e-8,'reference_force':1e-12,'AD_relative_error':2e-5,'scale_displacement':2e-5,'warm_cold_displacement':2e-5,'work_balance':1e-10,'area_energy_relative':1e-12},
              'test_gains':[.2,1.,4.],'basal_tension_scenarios':[.2,2.,20.],
              'base_source_sha256_required':FROZEN_BASE_SHA256,
              'limitations':['No continuum mesh error certificate','No experimentally identified tension range','No prediction of transient peak without kinetics']}
    (evidence/'prestress_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    checks={'original_source_unchanged':hashlib.sha256(base.read_bytes()).hexdigest()==FROZEN_BASE_SHA256}
    p=Parameters(**protocol['parameters']);m=build_model(p);zero=m.relax(0.)
    checks['flat_zero_dose_reference']=zero['converged'] and np.max(abs(np.array(zero['displacement'])))<1e-8
    _,reference_gradient=m.objective(np.zeros(len(m.free)),0.)
    checks['flat_reference_first_variation']=np.max(abs(reference_gradient))<1e-12
    reference_tension=float(sum(m.tension_components(m.X).values()))
    area_exact=np.pi*p.outer_radius**2*(p.basal_surface_tension+p.apical_surface_tension)
    area_error=abs(reference_tension-area_exact)/area_exact
    checks['reference_true_area_energy']=area_error<1e-12
    # Independent exact axisymmetric area of each straight, revolved P1 edge.
    rng=np.random.default_rng(2011);u=rng.normal(size=len(m.free))*1e-4
    y=m.coordinates(torch.tensor(u,dtype=torch.float64));yn=y.detach().numpy()
    exact=0.
    for nodes,tension in [(m.bottom,p.basal_surface_tension),(m.top,p.apical_surface_tension)]:
        edges=yn[nodes];exact+=tension*float(np.sum(np.pi*(edges[:-1,0]+edges[1:,0])*np.linalg.norm(np.diff(edges,axis=0),axis=1)))
    computed=float(sum(m.tension_components(y).values()))
    checks['deformed_true_area_crosscheck']=abs(computed-exact)<1e-11
    direction=rng.normal(size=len(u));direction/=np.linalg.norm(direction)
    _,gradient=m.objective(u,.3);h=1e-6
    fd=(m.objective(u+h*direction,.3,False)-m.objective(u-h*direction,.3,False))/(2*h)
    ad=float(gradient@direction);derivative_error=abs(fd-ad)/max(abs(fd),abs(ad),1e-10)
    checks['AD_vs_finite_difference']=derivative_error<2e-5
    rejected=False
    try:Parameters(reference_curvature=.01)
    except ValueError:rejected=True
    checks['unbalanced_curvature_rejected']=rejected
    start=time.monotonic();initial=m.relax(.2);warm=m.relax(1.,initial['free_displacement']);cold=m.solve(1.,maxiter=100)
    checks['warm_cold_stationarity']=all(q['converged'] for q in (initial,warm,cold))
    difference=float(np.max(abs(np.array(warm['displacement'])-np.array(cold['displacement']))))
    checks['warm_cold_agreement']=difference<2e-5
    ramp=m.solve(1.,initial=initial['free_displacement'],continuation_start_strength=.2,continuation_steps=1,maxiter=100)
    checks['work_ledger']=ramp['converged'] and abs(ramp['quasistatic_work_ledger']['balance_residual'])<1e-10 and ramp['quasistatic_work_ledger']['all_activation_increments_nonnegative'] and ramp['quasistatic_work_ledger']['all_relaxation_drops_nonnegative']
    scaled=build_model(replace(p,mu=7*p.mu,bulk=7*p.bulk,basal_area_modulus=7*p.basal_area_modulus,apical_area_modulus=7*p.apical_area_modulus,basal_surface_tension=7*p.basal_surface_tension,apical_surface_tension=7*p.apical_surface_tension)).relax(7.,np.array(warm['free_displacement']))
    scaling_difference=float(np.max(abs(np.array(scaled['displacement'])-np.array(warm['displacement']))))
    checks['material_and_tension_scaling']=scaled['converged'] and scaling_difference<2e-5
    scenarios={}
    for tension in protocol['basal_tension_scenarios']:
        model=build_model(replace(p,basal_surface_tension=tension,patch_radius=9.35512037/33.5,activation_z_min=(33.5-25)/33.5,activation_z_max=(33.5-5)/33.5))
        last=None;results=[]
        for gain in protocol['test_gains']:
            result=model.relax(gain,last);last=result['free_displacement']
            results.append({k:result[k] for k in ['activation_strength','converged','scaled_free_force_linf','minimum_J','central_apical_displacement','central_basal_displacement','maximum_activated_column_shortening','single_valued_surface_graphs']})
            (out/f'Tb_{tension}_gain_{gain}.json').write_text(json.dumps(result,indent=2)+'\n')
        scenarios[str(tension)]=results
    checks['all_scenarios_stationary_positiveJ']=all(q['converged'] and q['minimum_J']>0 for rows in scenarios.values() for q in rows)
    for name,result in [('zero',zero),('warm',warm),('cold',cold),('ramp',ramp),('scaled',scaled)]:
        (out/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    checks={name:bool(value) for name,value in checks.items()}
    receipt={'status':'PASS_SYNTHETIC_COMPONENT' if all(checks.values()) else 'FAIL_SYNTHETIC_CHECK','checks':checks,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'base_source_sha256':hashlib.sha256(base.read_bytes()).hexdigest(),'derivative_relative_error':derivative_error,'reference_area_relative_error':area_error,'reference_force_linf':float(np.max(abs(reference_gradient))),'warm_cold_displacement_difference':difference,'scaling_displacement_difference':scaling_difference,'scenarios':scenarios,'seconds':time.monotonic()-start,'limits':protocol['limitations']}
    (evidence/'prestress_verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
    return all(checks.values())


if __name__=='__main__':
    torch.set_num_threads(1)
    raise SystemExit(0 if _self_test() else 1)
