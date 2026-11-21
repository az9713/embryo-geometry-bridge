"""Flat axisymmetric P2 displacement / continuous P1 log-pressure model.

The projected volumetric energy is a finite-lambda mixed formulation, not the
pointwise penalty energy of the two frozen P1 implementations. No data reads.
"""
from dataclasses import dataclass,asdict,replace
from pathlib import Path
import hashlib,json,time
import numpy as np
import torch
from axisymmetric_prestress import Parameters as ParentParameters
from axisymmetric_prestress import AxisymmetricModel as ParentModel
from axisymmetric_mechanics import _quadrature


@dataclass(frozen=True)
class Parameters(ParentParameters):
    nr:int=8
    nz:int=4
    quadrature:int=4
    surface_quadrature:int=8

    def __post_init__(self):
        super().__post_init__()
        if self.quadrature!=4:raise ValueError('mixed model requires degree4 volume quadrature')
        if self.surface_quadrature<5:raise ValueError('mixed model requires at least5 Gauss surface points')


def _shape(bary,grad_lambda):
    """P2 Lagrange basis in [v0,v1,v2,e01,e12,e20] order."""
    q=np.asarray(bary);pairs=((0,1),(1,2),(2,0))
    N=np.column_stack([q[:,i]*(2*q[:,i]-1) for i in range(3)]+[4*q[:,i]*q[:,j] for i,j in pairs])
    grad=np.stack([(4*q[None,:,i,None]-1)*grad_lambda[:,None,i,:] for i in range(3)]+[4*(q[None,:,i,None]*grad_lambda[:,None,j,:]+q[None,:,j,None]*grad_lambda[:,None,i,:]) for i,j in pairs],axis=2)
    return N,grad


def _bernstein(values):
    """Degree2 nodal values at vertices and edge midpoints -> Bernstein."""
    c=np.array(values,copy=True)
    for edge,(i,j) in enumerate(((0,1),(1,2),(2,0))):
        c[...,3+edge]=2*values[...,3+edge]-.5*(values[...,i]+values[...,j])
    return c


class AxisymmetricModel(ParentModel):
    def __init__(self,parameters=Parameters()):
        if not isinstance(parameters,Parameters):raise TypeError('use axisymmetric_mixed.Parameters')
        super().__init__(parameters)
        p=parameters;vertices=self.reference.copy();tri=self.triangles.copy();self.pressure_vertices=vertices
        self.pressure_triangles=tri;self.pressure_count=len(vertices)
        midpoint={};reference=vertices.tolist();tri6=[]
        for triangle in tri:
            mids=[]
            for a,b in ((triangle[0],triangle[1]),(triangle[1],triangle[2]),(triangle[2],triangle[0])):
                key=tuple(sorted((int(a),int(b))))
                if key not in midpoint:
                    midpoint[key]=len(reference);reference.append(((vertices[a]+vertices[b])/2).tolist())
                mids.append(midpoint[key])
            tri6.append([*triangle,*mids])
        self.reference=np.array(reference);self.triangles=np.array(tri6,int)
        def enrich(nodes):
            edges=np.array([[int(a),midpoint[tuple(sorted((int(a),int(b))))],int(b)] for a,b in zip(nodes[:-1],nodes[1:])])
            ordered=np.r_[edges[:,:2].ravel(),edges[-1,-1]]
            return ordered,edges
        self.bottom,self.bottom_edges=enrich(self.bottom);self.top,self.top_edges=enrich(self.top)
        fixed=np.zeros_like(self.reference,dtype=bool)
        fixed[self.reference[:,0]==0,0]=True;fixed[self.reference[:,0]==p.outer_radius,:]=True
        self.fixed=fixed;self.free=np.flatnonzero(~fixed.ravel())
        self.bounds=[(-(1-1e-7)*self.reference.ravel()[i],None) if i%2==0 else (None,None) for i in self.free]
        self.X=torch.tensor(self.reference,dtype=torch.float64);self.free_t=torch.tensor(self.free,dtype=torch.long)
        self.tri=torch.tensor(self.triangles,dtype=torch.long)
        dX=vertices[tri];dm=np.stack((dX[:,1]-dX[:,0],dX[:,2]-dX[:,0]),axis=2)
        inv=np.linalg.inv(dm);grad_lambda=np.array([[-1.,-1.],[1.,0.],[0.,1.]])[None,:,:]@inv
        self.bary,self.weights=_quadrature(4);N,grad=_shape(self.bary,grad_lambda)
        self.shape=torch.tensor(N,dtype=torch.float64);self.shape_grad=torch.tensor(grad,dtype=torch.float64)
        qx=np.einsum('qa,tad->tqd',self.bary,dX)
        self.qr=torch.tensor(qx[:,:,0],dtype=torch.float64)
        self.integration=torch.tensor(2*np.pi*self.area[:,None]*qx[:,:,0]*self.weights[None,:],dtype=torch.float64)
        self.reference_volume=float(self.integration.sum());self.energy_scale=p.mu*self.reference_volume
        self.pressure_tri_t=torch.tensor(tri,dtype=torch.long)
        local_mass=np.einsum('tq,qi,qj->tij',self.integration.numpy(),self.bary,self.bary)
        mass=np.zeros((len(vertices),len(vertices)))
        for i in range(3):
            for j in range(3):np.add.at(mass,(tri[:,i],tri[:,j]),local_mass[:,i,j])
        self.mass=torch.tensor(mass,dtype=torch.float64);self.mass_inverse=torch.tensor(np.linalg.inv(mass),dtype=torch.float64)
        self.pressure_shape=torch.tensor(self.bary,dtype=torch.float64)
        # Gradients at degree2 interpolation nodes certify det(Fm) polynomials.
        check=np.array([[1,0,0],[0,1,0],[0,0,1],[.5,.5,0],[0,.5,.5],[.5,0,.5]])
        _,cert_grad=_shape(check,grad_lambda);self.cert_grad=torch.tensor(cert_grad,dtype=torch.float64)
        self.reference_r_bernstein=_bernstein(self.reference[self.triangles,0])
        self.surfaces=[]
        q,w=np.polynomial.legendre.leggauss(p.surface_quadrature);q=(q+1)/2;w=w/2
        shape=np.column_stack(((1-q)*(1-2*q),4*q*(1-q),q*(2*q-1)))
        deriv=np.column_stack((-3+4*q,4-8*q,-1+4*q))
        for edges,modulus,tension in [(self.bottom_edges,p.basal_area_modulus,p.basal_surface_tension),(self.top_edges,p.apical_area_modulus,p.apical_surface_tension)]:
            xy=self.reference[edges];r0=np.einsum('qa,tad->tqd',shape,xy)[:,:,0]
            length0=np.linalg.norm(xy[:,2]-xy[:,0],axis=1)
            self.surfaces.append((torch.tensor(edges,dtype=torch.long),torch.tensor(shape,dtype=torch.float64),torch.tensor(deriv,dtype=torch.float64),torch.tensor(r0,dtype=torch.float64),torch.tensor(length0,dtype=torch.float64),torch.tensor(2*np.pi*r0*length0[:,None]*w,dtype=torch.float64),modulus,tension))

    def energy_components(self,current,activation_strength=0.):
        p=self.parameters;y=current[self.tri]
        fm=torch.einsum('tad,tqaj->tqdj',y,self.shape_grad)
        detm=torch.linalg.det(fm);R=torch.einsum('qa,tad->tqd',self.shape,y)[:,:,0]
        hoop=R/self.qr;J=detm*hoop;logJ=torch.log(torch.clamp(J,min=1e-8))
        trC=(fm*fm).sum(dim=(2,3))+hoop**2
        shear=(self.integration*p.mu/2*(trC-3-2*logJ)).sum()
        local_b=torch.einsum('tq,tq,qi->ti',self.integration,logJ,self.pressure_shape)
        b=torch.zeros(self.pressure_count,dtype=torch.float64).index_add(0,self.pressure_tri_t.ravel(),local_b.ravel())
        pressure=p.bulk*(self.mass_inverse@b)
        volumetric=.5*(b@pressure)
        # Flat prescribed fibers are vertical, so their current vector is Fm e_z.
        active=(self.integration*self.active*activation_strength/2*(fm[:,:,:,1]**2).sum(dim=2)).sum()
        surface=torch.zeros((),dtype=torch.float64);tension_energy=torch.zeros_like(surface)
        for edges,N,dN,r0,l0,weight,modulus,tension in self.surfaces:
            yy=current[edges];radius=torch.einsum('qa,tad->tqd',N,yy)[:,:,0]
            tangent=torch.einsum('qa,tad->tqd',dN,yy)
            stretch=torch.linalg.vector_norm(tangent,dim=2)/l0[:,None]
            Js=radius/r0*stretch
            surface=surface+(weight*modulus/2*(Js-1)**2).sum()
            tension_energy=tension_energy+(weight*tension*Js).sum()
        invalid=(1e8*self.integration*torch.relu(1e-6-J)**2).sum()
        return {'bulk':shear+volumetric,'shear':shear,'projected_volumetric':volumetric,'surface':surface,'tension':tension_energy,'active':active,'trial_inversion_penalty':invalid,'total':shear+volumetric+surface+tension_energy+active+invalid,'J':J,'logJ':logJ,'F_meridional':fm,'meridional_determinant':detm,'hoop':hoop,'pressure':pressure,'pressure_rhs':b}

    def orientation_certificate(self,current):
        y=current[self.tri]
        fm=torch.einsum('tad,tqaj->tqdj',y,self.cert_grad).detach().numpy()
        det_coefficient=_bernstein(np.linalg.det(fm))
        radius_coefficient=_bernstein(current.detach().numpy()[self.triangles,0])
        reference=self.reference_r_bernstein;positive=reference>1e-14
        ratios=np.divide(radius_coefficient,reference,out=np.full_like(reference,np.nan),where=positive)
        axis_zero=bool(np.all(abs(radius_coefficient[~positive])<1e-12))
        dmin=det_coefficient.min(axis=1);dmax=det_coefficient.max(axis=1)
        hmin=np.nanmin(ratios,axis=1);hmax=np.nanmax(ratios,axis=1)
        certified=bool(axis_zero and np.all(dmin>0) and np.all(hmin>0))
        return {'certified_positive':certified,'axis_zero_coefficients_verified':axis_zero,
                'determinant_Bernstein_min':float(dmin.min()),'hoop_ratio_Bernstein_min':float(hmin.min()),
                'J_lower_bound':float(np.min(dmin*hmin)) if certified else None,
                'J_upper_bound':float(np.max(dmax*hmax)) if certified else None,
                'method':'Sufficient degree2 Bernstein coefficient bounds in floating point; not interval arithmetic; failure is uncertified, not a proof of inversion.'}

    def diagnostics(self,vector,strength,history=None):
        p=self.parameters;current=self.coordinates(torch.tensor(vector,dtype=torch.float64));parts=self.energy_components(current,strength)
        y=current.detach().numpy();certificate=self.orientation_certificate(current)
        _,gradient=self.objective(vector,strength);residual=float(np.max(abs(gradient)))
        stationary=bool(history and all(q['stationarity_verified'] for q in history) and residual<1e-7)
        physical=certificate['certified_positive'] and float(parts['trial_inversion_penalty'])==0.
        volume=float((self.integration*parts['J']).sum())
        b=parts['pressure_rhs'];pressure=parts['pressure']
        pressure_residual=float(torch.max(abs(self.mass@pressure/p.bulk-b)))
        top=y[self.top];bottom=y[self.bottom];top0=self.reference[self.top];bottom0=self.reference[self.bottom]
        surface_cert=self._surface_certificates(y)
        heights=top[:,1]-bottom[:,1];inside=top0[:,0]<=p.patch_radius+1e-12
        energy={name:float(parts[name]) for name in ['bulk','shear','projected_volumetric','surface','tension','active','trial_inversion_penalty','total']}
        fm=parts['F_meridional'].detach().numpy();singular=np.linalg.svd(fm,compute_uv=False)
        return {'parameters':asdict(p),'activation_strength':float(strength),'converged':bool(stationary and physical),'stationary_discrete_equations':stationary,'physical_local_orientation':bool(physical),'orientation_certificate':certificate,'scaled_free_force_linf':residual,'minimum_J':certificate['J_lower_bound'],'maximum_J':certificate['J_upper_bound'],'sampled_minimum_J':float(parts['J'].min()),'sampled_maximum_J':float(parts['J'].max()),'volume':volume,'reference_volume':self.reference_volume,'relative_volume_change':volume/self.reference_volume-1,'energy':energy,'continuation':history or [],'free_displacement':np.asarray(vector).tolist(),'displacement':(y-self.reference).tolist(),'current_nodes':y.tolist(),'reference_nodes':self.reference.tolist(),'triangles':self.triangles.tolist(),'activation_elements':self.activation_elements.tolist(),'top_profile':top.tolist(),'bottom_profile':bottom.tolist(),'initial_top_profile':top0.tolist(),'initial_bottom_profile':bottom0.tolist(),'vertical_material_column_heights':heights.tolist(),'maximum_activated_column_shortening':float(np.max(p.height-heights[inside])),'central_apical_displacement':float(top[0,1]-top0[0,1]),'central_basal_displacement':float(bottom[0,1]-bottom0[0,1]),'pressure_coefficients':pressure.detach().numpy().tolist(),'pressure_constraint_linf':pressure_residual,'mean_logJ':float((self.integration*parts['logJ']).sum()/self.reference_volume),'sampled_principal_stretch_min':float(min(singular.min(),parts['hoop'].min())),'sampled_principal_stretch_max':float(max(singular.max(),parts['hoop'].max())),'surface_geometry_representation':'Piecewise quadratic, not a node polyline; use exact edge controls for extrema and observer','boundary_top_quadratic_controls':y[self.top_edges].tolist(),'boundary_bottom_quadratic_controls':y[self.bottom_edges].tolist(),'surface_certificates':surface_cert,'single_valued_surface_graphs':surface_cert['radially_monotone'],'units':'Caller declared; defaults H=1,mu=1; +Z apical. No time or ATP calibration.','scope':'Flat idealized P2 displacement/P1 projected-log-pressure quasistatic mechanics; stationarity and sufficient positivity certificate distinct; no target data or continuum error certificate.'}

    def _surface_certificates(self,y):
        # In t in[0,1], y=c0+c1*t+c2*t^2 for each exact boundary edge.
        def controls(edges):
            q=y[edges];c0=q[:,0];c1=-3*q[:,0]+4*q[:,1]-q[:,2];c2=2*q[:,0]-4*q[:,1]+2*q[:,2]
            radial=np.minimum(c1[:,0],c1[:,0]+2*c2[:,0])
            values=[q[:,0,1],q[:,2,1]]
            roots=np.divide(-c1[:,1],2*c2[:,1],out=np.full(len(q),np.nan),where=abs(c2[:,1])>1e-15)
            valid=(roots>0)&(roots<1)
            if valid.any():values.append((c0[valid]+roots[valid,None]*c1[valid]+roots[valid,None]**2*c2[valid])[:,1])
            z=np.concatenate(values)
            return {'radial_derivative_min':float(radial.min()),'Z_min':float(z.min()),'Z_max':float(z.max())}
        a=controls(self.top_edges);b=controls(self.bottom_edges)
        axis=y[np.flatnonzero(self.reference[:,0]==0)]
        axis=axis[np.argsort(self.reference[self.reference[:,0]==0,1])]
        # Boundary axis uses vertices and edge mids; derivative endpoint bounds.
        axis_mono=True
        for z0,zm,z1 in zip(axis[:-2:2,1],axis[1:-1:2,1],axis[2::2,1]):
            d0=-3*z0+4*zm-z1;d1=z0-4*zm+3*z1
            axis_mono &= min(d0,d1)>0
        radial=a['radial_derivative_min']>0 and b['radial_derivative_min']>0
        simple=bool(radial and axis_mono and a['Z_min']>b['Z_max'])
        return {'apical':a,'basal':b,'radially_monotone':bool(radial),'axis_monotone':bool(axis_mono),'strict_vertical_band_separation':bool(a['Z_min']>b['Z_max']),'sufficient_simple_meridian_boundary':simple,'scope':'Sufficient exact-quadratic algebra evaluated in floating point. No interval rounding or general global-injectivity theorem claimed.'}


def build_model(parameters=None,**kwargs):
    return AxisymmetricModel(Parameters(**kwargs) if parameters is None else parameters)


def _self_test():
    root=Path(__file__).resolve().parents[1];out=root/'runs/mixed_checks';out.mkdir(parents=True,exist_ok=True)
    evidence=root/'evidence';evidence.mkdir(exist_ok=True)
    protocol={'scope':'Synthetic P2/P1 projected-log-pressure mechanics; no outcome reads','parameters':asdict(Parameters(nr=6,nz=3,basal_surface_tension=2.,apical_surface_tension=.2)),
              'thresholds':{'reference_force':1e-12,'zero_displacement':1e-8,'affine_J_energy_relative':1e-11,'AD_relative':2e-5,'pressure_constraint':1e-11,'force_linf':1e-7,'scaled_displacement':2e-5},
              'interpretation':'Projection is a changed finite-lambda mixed approximation, not pointwise incompressibility or an automatic cure.'}
    (evidence/'mixed_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    m=build_model(Parameters(**protocol['parameters']));checks={}
    zero=m.relax(0.);checks['zero_reference_stationary']=zero['converged'] and np.max(abs(np.array(zero['displacement'])))<1e-8
    checks['reference_positive_Bernstein']=zero['orientation_certificate']['certified_positive']
    checks['pressure_mass_positive_definite']=float(torch.linalg.eigvalsh(m.mass).min())>0
    # Exact homogeneous isochoric deformation, independent of the fixed-rim solve.
    lam=1.1;affine=m.X.clone();affine[:,0]*=lam;affine[:,1]/=lam**2
    parts=m.energy_components(affine,0.);p=m.parameters
    exact=m.reference_volume*p.mu/2*(2*lam**2+lam**-4-3)+np.pi*p.outer_radius**2*((p.basal_area_modulus+p.apical_area_modulus)/2*(lam**2-1)**2+(p.basal_surface_tension+p.apical_surface_tension)*lam**2)
    energy_error=abs(float(parts['total'])-exact)/abs(exact)
    checks['isochoric_affine_J']=float(torch.max(abs(parts['J']-1)))<1e-11
    checks['isochoric_affine_volume']=abs(float((m.integration*parts['J']).sum())/m.reference_volume-1)<1e-11
    checks['isochoric_affine_energy']=energy_error<1e-11
    checks['affine_orientation_certificate']=m.orientation_certificate(affine)['certified_positive']
    rng=np.random.default_rng(1811);u=rng.normal(size=len(m.free))*1e-5;d=rng.normal(size=len(u));d/=np.linalg.norm(d)
    _,g=m.objective(u,.2);h=1e-6;fd=(m.objective(u+h*d,.2,False)-m.objective(u-h*d,.2,False))/(2*h);ad=float(g@d)
    derivative_error=abs(fd-ad)/max(abs(fd),abs(ad),1e-10);checks['AD_vs_FD']=derivative_error<2e-5
    started=time.monotonic();activated=m.relax(.2,maxiter=150)
    checks['activated_stationarity_positive_cert']=activated['converged']
    checks['projected_pressure_constraint']=activated['pressure_constraint_linf']<1e-11
    scaled=build_model(replace(p,mu=7*p.mu,bulk=7*p.bulk,basal_area_modulus=7*p.basal_area_modulus,apical_area_modulus=7*p.apical_area_modulus,basal_surface_tension=7*p.basal_surface_tension,apical_surface_tension=7*p.apical_surface_tension)).relax(1.4,activated['free_displacement'])
    scale_difference=float(np.max(abs(np.array(scaled['displacement'])-np.array(activated['displacement']))));checks['material_scaling']=scaled['converged'] and scale_difference<2e-5
    for name,result in [('zero',zero),('activated',activated),('scaled',scaled)]:
        (out/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    checks={name:bool(value) for name,value in checks.items()}
    receipt={'status':'PASS_SYNTHETIC_COMPONENT' if all(checks.values()) else 'FAIL_SYNTHETIC_CHECK','checks':checks,'affine_energy_relative_error':energy_error,'AD_relative_error':derivative_error,'scaling_displacement_difference':scale_difference,'active_summary':{k:activated[k] for k in ['converged','scaled_free_force_linf','minimum_J','maximum_J','sampled_minimum_J','sampled_maximum_J','pressure_constraint_linf','relative_volume_change']},'seconds':time.monotonic()-started,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'parent_sources':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ['axisymmetric_mechanics.py','axisymmetric_prestress.py']},'limitations':['Projected logJ does not imply pointwise J close to1','Bernstein floating-point sufficiency, not validated interval bound','No mesh-convergence or biology claim','Quadratic surfaces require an appropriate observer']}
    (evidence/'mixed_verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));return all(checks.values())


if __name__=='__main__':
    torch.set_num_threads(1)
    raise SystemExit(0 if _self_test() else 1)
