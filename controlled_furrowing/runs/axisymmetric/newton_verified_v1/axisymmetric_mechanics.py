"""Finite-height axisymmetric solid mechanics, with no experimental-data reads.

Reference/current meridians are revolved about the z axis. Units are declared
by the caller; defaults use thickness H=1 and shear modulus mu=1. This is a
conditional idealized equilibrium model, not a specimen or transient fit.
"""
from dataclasses import dataclass, asdict, replace
from pathlib import Path
import argparse, hashlib, json, math, time
import numpy as np
import torch
from scipy.optimize import minimize

DTYPE = torch.float64


@dataclass(frozen=True)
class Parameters:
    nr: int = 12
    nz: int = 4
    height: float = 1.
    outer_radius: float = 3.
    patch_radius: float = .28
    activation_z_min: float = .22
    activation_z_max: float = .84
    mu: float = 1.
    bulk: float = 50.
    basal_area_modulus: float = 2.
    apical_area_modulus: float = .2
    reference_curvature: float = 0.
    quadrature: int = 2
    radial_inner_cells: int = 0

    def __post_init__(self):
        if self.nr < 6 or self.nz < 3 or int(self.nr) != self.nr or int(self.nz) != self.nz:
            raise ValueError('nr>=6 and nz>=3 integer cells required')
        if not 0 < self.patch_radius < self.outer_radius or self.height <= 0:
            raise ValueError('positive height and patch radius below outer radius required')
        if not 0 < self.activation_z_min < self.activation_z_max < 1:
            raise ValueError('strictly interior activation depth fractions required')
        if self.mu <= 0 or self.bulk <= 0 or min(self.basal_area_modulus,self.apical_area_modulus) < 0:
            raise ValueError('positive bulk/shear and nonnegative surface moduli required')
        if self.quadrature not in (1,2,4):
            raise ValueError('triangle quadrature must be1,2,4')


def _layer_grid(n, breaks):
    lengths=np.diff(breaks);counts=np.ones(len(lengths),int)
    for _ in range(n-len(lengths)):
        counts[np.argmax(lengths/counts)]+=1
    return np.concatenate([np.linspace(a,b,c+1)[:-1] for a,b,c in zip(breaks[:-1],breaks[1:],counts)]+[[breaks[-1]]])


def _quadrature(order):
    if order == 1: return np.array([[1/3,1/3,1/3]]),np.array([1.])
    if order == 2: return np.array([[2/3,1/6,1/6],[1/6,2/3,1/6],[1/6,1/6,2/3]]),np.full(3,1/3)
    a,b=.445948490915965,.108103018168070
    c,d=.091576213509771,.816847572980459
    return np.array([[a,a,b],[a,b,a],[b,a,a],[c,c,d],[c,d,c],[d,c,c]]),np.array([.223381589678011]*3+[.109951743655322]*3)


class AxisymmetricModel:
    def __init__(self, parameters=Parameters()):
        self.parameters=p=parameters
        ni=p.radial_inner_cells or max(3,p.nr//3)
        if not 1 <= ni < p.nr: raise ValueError('inner cells must lie below nr')
        self.radial=np.r_[np.linspace(0,p.patch_radius,ni+1),np.linspace(p.patch_radius,p.outer_radius,p.nr-ni+1)[1:]]
        self.layers=_layer_grid(p.nz,np.array([0,p.activation_z_min,p.activation_z_max,1.]))
        self.logical=np.array([(r,s) for r in self.radial for s in self.layers])
        self.reference=self.logical.copy()
        self.reference[:,1]=p.height*self.logical[:,1]+.5*p.reference_curvature*self.logical[:,0]**2
        idx=lambda i,j:i*(p.nz+1)+j
        self.triangles=np.array([(idx(i,j),idx(i+1,j),idx(i+1,j+1)) for i in range(p.nr) for j in range(p.nz)]+[(idx(i,j),idx(i+1,j+1),idx(i,j+1)) for i in range(p.nr) for j in range(p.nz)],int)
        self.bottom=np.array([idx(i,0) for i in range(p.nr+1)])
        self.top=np.array([idx(i,p.nz) for i in range(p.nr+1)])
        fixed=np.zeros(self.reference.shape,bool);fixed[self.logical[:,0]==0,0]=True
        fixed[self.logical[:,0]==p.outer_radius,:]=True
        self.free=np.flatnonzero(~fixed.ravel())
        self.fixed=fixed
        self.bounds=[(-(1-1e-7)*self.reference.ravel()[i],None) if i%2==0 else (None,None) for i in self.free]
        X=self.reference[self.triangles]
        dm=np.stack((X[:,1]-X[:,0],X[:,2]-X[:,0]),axis=2)
        self.area=np.linalg.det(dm)/2
        assert self.area.min()>0
        self.bary,self.weights=_quadrature(p.quadrature)
        qx=np.einsum('qa,tad->tqd',self.bary,X)
        qlogical=np.einsum('qa,tad->tqd',self.bary,self.logical[self.triangles])
        centre=self.logical[self.triangles].mean(axis=1)
        active=(centre[:,0] < p.patch_radius)&(centre[:,1] > p.activation_z_min)&(centre[:,1] < p.activation_z_max)
        self.activation_elements=active
        self.dm_inverse=torch.tensor(np.linalg.inv(dm),dtype=DTYPE)
        self.qr=torch.tensor(qx[:,:,0],dtype=DTYPE)
        self.integration=torch.tensor(2*np.pi*self.area[:,None]*qx[:,:,0]*self.weights[None,:],dtype=DTYPE)
        self.n0=torch.tensor(np.stack((-p.reference_curvature*qlogical[:,:,0],np.ones_like(qlogical[:,:,0])),axis=2),dtype=DTYPE)
        self.n0=self.n0/torch.linalg.vector_norm(self.n0,dim=2,keepdim=True)
        self.active=torch.tensor(active[:,None],dtype=DTYPE)
        self.X=torch.tensor(self.reference,dtype=DTYPE)
        self.tri=torch.tensor(self.triangles,dtype=torch.long)
        self.bary_t=torch.tensor(self.bary,dtype=DTYPE)
        self.free_t=torch.tensor(self.free,dtype=torch.long)
        self.surfaces=[]
        # Three-point Gauss is used for every surface energy.
        ab,aw=np.polynomial.legendre.leggauss(3);ab=(ab+1)/2;aw=aw/2
        for nodes,modulus in [(self.bottom,p.basal_area_modulus),(self.top,p.apical_area_modulus)]:
            edges=np.column_stack((nodes[:-1],nodes[1:]));xx=self.reference[edges]
            length=np.linalg.norm(xx[:,1]-xx[:,0],axis=1)
            rr=(1-ab[None,:])*xx[:,0,0,None]+ab[None,:]*xx[:,1,0,None]
            self.surfaces.append((torch.tensor(edges,dtype=torch.long),torch.tensor(length,dtype=DTYPE),torch.tensor(rr,dtype=DTYPE),torch.tensor(ab,dtype=DTYPE),torch.tensor(2*np.pi*length[:,None]*rr*aw[None,:],dtype=DTYPE),modulus))
        self.reference_volume=float(self.integration.sum())
        self.energy_scale=p.mu*self.reference_volume

    def coordinates(self, free_displacement):
        vector=torch.zeros(self.X.numel(),dtype=DTYPE).scatter(0,self.free_t,free_displacement)
        return self.X+vector.reshape_as(self.X)

    def energy_components(self, current, activation_strength=0.):
        p=self.parameters;y=current[self.tri]
        ds=torch.stack((y[:,1]-y[:,0],y[:,2]-y[:,0]),dim=2)
        fm=ds@self.dm_inverse;detm=torch.linalg.det(fm)
        R=torch.einsum('qa,tad->tqd',self.bary_t,y)[:,:,0]
        hoop=R/self.qr;J=detm[:,None]*hoop
        safeJ=torch.clamp(J,min=1e-8);logJ=torch.log(safeJ)
        trC=(fm*fm).sum(dim=(1,2))[:,None]+hoop**2
        bulk_density=p.mu/2*(trC-3-2*logJ)+p.bulk/2*logJ**2
        bulk=(self.integration*bulk_density).sum()
        fibre=torch.einsum('tij,tqj->tqi',fm,self.n0)
        active=(self.integration*self.active*activation_strength/2*(fibre*fibre).sum(dim=2)).sum()
        surface=torch.zeros((),dtype=DTYPE)
        for edges,l0,r0,q,weight,modulus in self.surfaces:
            yc=current[edges];length=torch.linalg.vector_norm(yc[:,1]-yc[:,0],dim=1)
            rq=(1-q[None,:])*yc[:,0,0,None]+q[None,:]*yc[:,1,0,None]
            Js=(rq/r0)*(length/l0)[:,None]
            surface=surface+(weight*modulus/2*(Js-1)**2).sum()
        # Extension used only to reject inverted optimizer trial states. It is
        # identically zero on every accepted positive-J physical state.
        invalid=(1e8*self.integration*torch.relu(1e-6-J)**2).sum()
        return {'bulk':bulk,'surface':surface,'active':active,'trial_inversion_penalty':invalid,'total':bulk+surface+active+invalid,'J':J,'meridional_determinant':detm,'F_meridional':fm,'hoop':hoop}

    def objective(self, vector, activation_strength, gradient=True):
        u=torch.tensor(np.asarray(vector),dtype=DTYPE,requires_grad=gradient)
        energy=self.energy_components(self.coordinates(u),activation_strength)['total']/self.energy_scale
        if gradient:
            derivative=torch.autograd.grad(energy,u)[0]
            return float(energy.detach()),derivative.detach().numpy()
        return float(energy.detach())

    def _newton_polish(self, vector, strength, iterations=5):
        history=[]
        for _ in range(iterations):
            old,gradient=self.objective(vector,strength)
            if np.max(abs(gradient))<2e-10:break
            value=torch.tensor(vector,dtype=DTYPE,requires_grad=True)
            hessian=torch.autograd.functional.hessian(lambda q:self.energy_components(self.coordinates(q),strength)['total']/self.energy_scale,value,vectorize=True).detach().numpy()
            hessian=(hessian+hessian.T)/2
            eigen_min=float(np.linalg.eigvalsh(hessian)[0])
            regularization=max(0.,1e-10-eigen_min)
            direction=np.linalg.solve(hessian+regularization*np.eye(len(vector)),-gradient)
            accepted=False
            for exponent in range(25):
                fraction=2.**(-exponent);trial=vector+fraction*direction
                if any(low is not None and trial[i]<low for i,(low,_) in enumerate(self.bounds)):continue
                trial_value=self.objective(trial,strength,False)
                if trial_value <= old+1e-4*fraction*(gradient@direction)+1e-16:
                    vector=trial;accepted=True;break
            history.append({'newton_eigenvalue_min':eigen_min,'regularization':regularization,'line_fraction':fraction,'accepted':accepted})
            if not accepted:break
        return vector,history

    def solve(self, activation_strength, continuation_steps=4, initial=None, maxiter=1200):
        if not np.isfinite(activation_strength) or activation_strength<0:raise ValueError('nonnegative finite activation required')
        vector=np.zeros(len(self.free)) if initial is None else np.asarray(initial,float).copy()
        history=[];previous_gain=0.;initial_energy=self.objective(vector,0,False)*self.energy_scale
        for gain in np.linspace(0,activation_strength,max(2,continuation_steps+1)):
            energy_at_previous_gain=self.objective(vector,previous_gain,False)*self.energy_scale
            energy_before_relaxation=self.objective(vector,float(gain),False)*self.energy_scale
            result=minimize(lambda u:self.objective(u,float(gain)),vector,jac=True,method='L-BFGS-B',bounds=self.bounds,options={'maxiter':maxiter,'ftol':1e-14,'gtol':2e-9,'maxls':40,'maxcor':30})
            vector,polish=self._newton_polish(result.x,float(gain))
            final_energy,final_gradient=self.objective(vector,float(gain))
            history.append({'activation_strength':float(gain),'optimizer_success':bool(result.success),'message':str(result.message),'iterations':int(result.nit),'objective':final_energy,'scaled_gradient_linf':float(np.max(abs(final_gradient))),'newton_polish':polish,'energy_before_relaxation':energy_before_relaxation,'energy_after_relaxation':final_energy*self.energy_scale,'external_activation_work':energy_before_relaxation-energy_at_previous_gain,'relaxation_energy_drop':energy_before_relaxation-final_energy*self.energy_scale})
            previous_gain=float(gain)
        answer=self.diagnostics(vector,float(activation_strength),history)
        work=sum(r['external_activation_work'] for r in history);drop=sum(r['relaxation_energy_drop'] for r in history)
        answer['quasistatic_work_ledger']={'initial_energy':initial_energy,'external_activation_work':work,'relaxation_energy_drop':drop,'energy_change':answer['energy']['total']-initial_energy,'balance_residual':answer['energy']['total']-initial_energy-work+drop,'all_activation_increments_nonnegative':all(r['external_activation_work']>=-1e-12 for r in history),'all_relaxation_drops_nonnegative':all(r['relaxation_energy_drop']>=-1e-12 for r in history),'interpretation':'Mechanical potential bookkeeping for externally ramped activation and idealized relaxation; no chemical/ATP/photon budget or biological timescale.'}
        return answer

    def diagnostics(self, vector, strength, history=None):
        p=self.parameters
        coordinates=self.coordinates(torch.tensor(vector,dtype=DTYPE))
        parts=self.energy_components(coordinates,strength)
        fm=parts['F_meridional'].detach().numpy();detm=np.linalg.det(fm)
        current=coordinates.detach().numpy()
        # R/r is a weighted average of non-axis nodal R_i/r_i. These vertex
        # extrema, times constant det(F_mer), give exact P1-element J extrema
        # in real arithmetic; no quadrature-only positivity assertion.
        ratios=np.divide(current[:,0],self.reference[:,0],out=np.full(len(current),np.nan),where=self.reference[:,0]>0)
        element_ratios=ratios[self.triangles]
        minJ=float(np.min(detm*np.nanmin(element_ratios,axis=1)))
        maxJ=float(np.max(detm*np.nanmax(element_ratios,axis=1)))
        volume=float((self.integration*parts['J']).sum())
        energy={k:float(parts[k]) for k in ['bulk','surface','active','trial_inversion_penalty','total']}
        _,gradient=self.objective(vector,strength)
        residual=float(np.max(abs(gradient)))
        top=current[self.top];bottom=current[self.bottom]
        top0=self.reference[self.top];bottom0=self.reference[self.bottom]
        monotone=bool(np.all(np.diff(top[:,0])>0) and np.all(np.diff(bottom[:,0])>0))
        graph={}
        if monotone:
            grid=np.unique(np.r_[top[:,0],bottom[:,0]])
            a=np.interp(grid,top[:,0],top[:,1]);b=np.interp(grid,bottom[:,0],bottom[:,1])
            graph={'graph_radial_grid':grid.tolist(),'graph_apical_z':a.tolist(),'graph_basal_z':b.tolist(),'graph_initial_apical_z':np.interp(grid,top0[:,0],top0[:,1]).tolist(),'graph_initial_basal_z':np.interp(grid,bottom0[:,0],bottom0[:,1]).tolist(),'minimum_graph_thickness':float(np.min(a-b))}
        vertical_height=top[:,1]-bottom[:,1]
        inside=self.radial<=p.patch_radius+1e-12
        physical=bool(minJ>1e-6 and np.all(current[:,0]>=0) and energy['trial_inversion_penalty']==0)
        converged=bool(history and all(r['optimizer_success'] for r in history) and residual<1e-7 and physical)
        return {'parameters':asdict(p),'activation_strength':strength,'converged':converged,'physical_local_orientation':physical,'scaled_free_force_linf':residual,'minimum_J':minJ,'maximum_J':maxJ,'volume':volume,'reference_volume':self.reference_volume,'relative_volume_change':volume/self.reference_volume-1,'energy':energy,'continuation':history or [],'radial_material_grid':self.radial.tolist(),'top_profile':top.tolist(),'bottom_profile':bottom.tolist(),'initial_top_profile':top0.tolist(),'initial_bottom_profile':bottom0.tolist(),'displacement':(current-self.reference).tolist(),'current_nodes':current.tolist(),'reference_nodes':self.reference.tolist(),'triangles':self.triangles.tolist(),'vertical_material_column_heights':vertical_height.tolist(),'material_fibre_lengths':np.linalg.norm(top-bottom,axis=1).tolist(),'maximum_activated_column_shortening':float(np.max(p.height-vertical_height[inside])),'central_apical_displacement':float(top[0,1]-top0[0,1]),'central_basal_displacement':float(bottom[0,1]-bottom0[0,1]),'apical_radial_monotone':bool(np.all(np.diff(top[:,0])>0)),'basal_radial_monotone':bool(np.all(np.diff(bottom[:,0])>0)),'single_valued_surface_graphs':monotone,'units':'Caller declared; defaults thickness H=1 and stress mu=1; +z apically/outward. No dimensional force or time calibration.', 'scope':'Idealized finite-height axisymmetric quasistatic equilibrium; no target data, no ellipse, no transient peak-to-recovery prediction.',**graph}


def build_model(parameters=None, **kwargs):
    return AxisymmetricModel(Parameters(**kwargs) if parameters is None else parameters)


def _self_test():
    root=Path(__file__).resolve().parents[1]
    out=root/'runs/axisymmetric';out.mkdir(parents=True,exist_ok=True)
    evidence=root/'evidence';evidence.mkdir(exist_ok=True)
    protocol={'scope':'Synthetic idealized axisymmetric mechanics, no biological fit','base':asdict(Parameters(nr=8,nz=4)),'activation_strength':.2,'thresholds':{'zero_displacement':1e-8,'affine_J_and_volume_error':1e-12,'affine_energy_relative_error':1e-11,'directional_derivative_relative_error':2e-5,'scaled_gradient_linf':1e-7,'material_scale_displacement':2e-5},'sensitivity_rule':'Report mesh/quadrature/domain/surface/reference-curvature changes; no empirical differences relabeled rigorous error bounds.'}
    (evidence/'axisymmetric_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    model=build_model(Parameters(nr=8,nz=4));zero=model.solve(0,continuation_steps=1)
    checks={'zero_dose_reference':max(abs(np.array(zero['displacement'])).ravel())<1e-8 and zero['converged']}
    lam=1.1;affine=model.X.clone();affine[:,0]*=lam;affine[:,1]/=lam**2
    comp=model.energy_components(affine,0)
    exact=model.reference_volume*model.parameters.mu/2*(2*lam**2+lam**-4-3)+(model.parameters.basal_area_modulus+model.parameters.apical_area_modulus)/2*(lam**2-1)**2*np.pi*model.parameters.outer_radius**2
    checks['isochoric_affine_J']=float(torch.max(abs(comp['J']-1)))<1e-12
    checks['isochoric_affine_volume']=abs(float((model.integration*comp['J']).sum())/model.reference_volume-1)<1e-12
    affine_error=abs(float(comp['total'])-exact)/max(abs(exact),1e-20)
    checks['affine_bulk_surface_energy']=affine_error<1e-11
    rng=np.random.default_rng(1911);u=rng.normal(size=len(model.free))*1e-5;direction=rng.normal(size=len(u));direction/=np.linalg.norm(direction)
    _,grad=model.objective(u,.2);step=1e-6
    finite=(model.objective(u+step*direction,.2,False)-model.objective(u-step*direction,.2,False))/(2*step)
    ad=float(grad@direction);derivative_error=abs(finite-ad)/max(abs(finite),abs(ad),1e-10)
    checks['automatic_derivative_vs_central_difference']=derivative_error<2e-5
    start=time.monotonic();base=model.solve(.2)
    (out/'base.json').write_text(json.dumps(base,indent=2)+'\n')
    checks['activated_converged_positive_J']=base['converged'] and base['minimum_J']>0
    checks['external_activation_work_balance']=abs(base['quasistatic_work_ledger']['balance_residual'])<1e-10 and base['quasistatic_work_ledger']['all_activation_increments_nonnegative'] and base['quasistatic_work_ledger']['all_relaxation_drops_nonnegative']
    scaled_parameters=replace(model.parameters,mu=7,bulk=350,basal_area_modulus=14,apical_area_modulus=1.4)
    scaled=build_model(scaled_parameters).solve(1.4)
    scale_difference=float(np.max(abs(np.array(base['displacement'])-np.array(scaled['displacement']))))
    checks['material_scale_invariance']=scaled['converged'] and scale_difference<2e-5
    cases={}
    for name,parameter in [('mesh_12x6',replace(model.parameters,nr=12,nz=6)),('mesh_16x8',replace(model.parameters,nr=16,nz=8)),('quadrature_4',replace(model.parameters,quadrature=4)),('outer_radius_4',replace(model.parameters,outer_radius=4)),('basal_modulus_0',replace(model.parameters,basal_area_modulus=0)),('curvature_positive',replace(model.parameters,reference_curvature=.05)),('curvature_negative',replace(model.parameters,reference_curvature=-.05))]:
        result=build_model(parameter).solve(.2)
        (out/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
        cases[name]={k:result[k] for k in ['converged','minimum_J','relative_volume_change','central_apical_displacement','central_basal_displacement','maximum_activated_column_shortening','single_valued_surface_graphs','scaled_free_force_linf']}
    checks['all_sensitivity_solves_converged']=all(v['converged'] for v in cases.values())
    receipt={'verification':'PASS_SYNTHETIC_COMPONENT' if all(checks.values()) else 'FAIL_SYNTHETIC_CHECK','checks':checks,'affine_energy_relative_error':affine_error,'directional_derivative_relative_error':derivative_error,'scale_displacement_difference':scale_difference,'base':{k:base[k] for k in ['converged','minimum_J','relative_volume_change','central_apical_displacement','central_basal_displacement','maximum_activated_column_shortening','scaled_free_force_linf']},'sensitivity':cases,'seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'limits':['Quasistatic local equilibrium only; no transient recovery/peak mapping.','Geometry and material ratios are declared idealizations.','P1 near-incompressibility and boundary sensitivity are not rigorous continuum error bounds.','Local positiveJ does not prove global injectivity or absence of self-contact.','Activation is externally controlled mechanical potential; no ATP/photon budget or novelty claim.']}
    (evidence/'axisymmetric_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
    return all(checks.values())


if __name__=='__main__':
    torch.set_num_threads(1)
    parser=argparse.ArgumentParser();parser.add_argument('--self-test',action='store_true');parser.add_argument('--gain',type=float,default=.2);args=parser.parse_args()
    if args.self_test:raise SystemExit(0 if _self_test() else 1)
    print(json.dumps(build_model().solve(args.gain),indent=2))
