from pathlib import Path
from dataclasses import replace
import json,hashlib,sys,time
import numpy as np
import torch
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from axisymmetric_mixed import Parameters,build_model
root=Path('controlled_furrowing');out=root/'runs/mixed_checks';torch.set_num_threads(1)
protocol={'map':'R=r exp(e z), Z=(1-exp(-2e z))/(2e), exact J=1','epsilon':.2,'grids':[[6,3],[8,4],[12,6],[16,8]],'scope':'Interpolation/energy consistency, not forced PDE convergence or validation; map need not satisfy clamped-rim boundary condition','nonisochoric_affine_tolerance':1e-10}
(root/'evidence/mixed_manufactured_protocol.json').write_text(json.dumps(protocol,indent=2))
rows=[]
for nr,nz in protocol['grids']:
 m=build_model(Parameters(nr=nr,nz=nz,basal_surface_tension=0,apical_surface_tension=0))
 eps=protocol['epsilon'];X=m.reference.copy();y=X.copy();y[:,0]=X[:,0]*np.exp(eps*X[:,1]);y[:,1]=-np.expm1(-2*eps*X[:,1])/(2*eps)
 c=m.energy_components(torch.tensor(y,dtype=torch.float64),0.)
 H=m.parameters.height;R=m.parameters.outer_radius;I2=np.expm1(2*eps*H)/(2*eps);Im4=-np.expm1(-4*eps*H)/(4*eps)
 shear_exact=np.pi*m.parameters.mu*(R**2/2*(2*I2+Im4-3*H)+eps**2*R**4/4*I2)
 row={'nr':nr,'nz':nz,'maximum_sampled_J_error':float(torch.max(abs(c['J']-1))),'shear_energy':float(c['shear']),'shear_energy_exact_continuum':float(shear_exact),'relative_shear_energy_error':abs(float(c['shear'])/shear_exact-1),'projected_volumetric_energy':float(c['projected_volumetric']),'volume_error':float((m.integration*c['J']).sum()/m.reference_volume-1),'orientation':m.orientation_certificate(torch.tensor(y,dtype=torch.float64))}
 rows.append(row);print({k:row[k] for k in ['nr','nz','maximum_sampled_J_error','relative_shear_energy_error']},flush=True)
# Nonisochoric homogeneous affine map tests nonzero pressure exactly.
m=build_model(Parameters(nr=6,nz=3));y=m.X.clone();y[:,0]*=1.02;y[:,1]*=.98;c=m.energy_components(y,0.);J=1.02**2*.98;expected=m.parameters.bulk*np.log(J)
p_error=float(torch.max(abs(c['pressure']-expected)));vexpected=.5*m.parameters.bulk*m.reference_volume*np.log(J)**2;energy_error=abs(float(c['projected_volumetric'])-vexpected)
receipt={'status':'PASS_CONSISTENCY_COMPONENT' if p_error<1e-10 and energy_error<1e-10 and all(q['orientation']['certified_positive'] for q in rows) else 'FAIL','rows':rows,'nonisochoric_affine_pressure_error':p_error,'nonisochoric_affine_volume_energy_error':energy_error,'source_sha256':hashlib.sha256(Path('controlled_furrowing/scripts/axisymmetric_mixed.py').read_bytes()).hexdigest(),'scope':protocol['scope']}
(root/'evidence/mixed_manufactured.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
