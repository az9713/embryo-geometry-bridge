from pathlib import Path
from dataclasses import fields
import sys,json,hashlib
import numpy as np
import torch
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from axisymmetric_prestress import build_model,Parameters
from axisymmetric_mechanics import build_model as build_base,Parameters as BaseParameters
torch.set_num_threads(1)
p=Parameters(nr=8,nz=4,basal_surface_tension=2.,apical_surface_tension=.2);m=build_model(p)
v=np.zeros_like(m.reference);r=m.reference[:,0];s=m.logical[:,1];v[:,1]=(1-(r/p.outer_radius)**2)**2*(1-s)
t=torch.tensor(0.,dtype=torch.float64,requires_grad=True);value=sum(m.tension_components(m.X+t*torch.tensor(v,dtype=torch.float64)).values());first=torch.autograd.grad(value,t,create_graph=True)[0];second=float(torch.autograd.grad(first,t)[0])
expected=0.
for nodes,tension in [(m.bottom,p.basal_surface_tension),(m.top,p.apical_surface_tension)]:
 rr=m.reference[nodes,0];vv=v[nodes,1];expected+=float(np.sum(np.pi*tension*(rr[1:]**2-rr[:-1]**2)*(np.diff(vv)/np.diff(rr))**2))
p0=Parameters(nr=8,nz=4);new=build_model(p0);base=build_base(BaseParameters(**{f.name:getattr(p0,f.name) for f in fields(BaseParameters)}));u=np.random.default_rng(19).normal(size=len(new.free))*1e-5
e,g=new.objective(u,.3);eb,gb=base.objective(u,.3)
checks={'true_area_normal_second_variation':abs(second-expected)<1e-11,'zero_tension_exact_original_objective':e==eb and np.array_equal(g,gb)}
receipt={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'second_variation_autodiff':second,'second_variation_independent_edge_formula':expected,'absolute_error':abs(second-expected),'zero_tension_objective_difference':e-eb,'zero_tension_gradient_linf_difference':float(np.max(abs(g-gb))),'source_sha256':hashlib.sha256(Path('controlled_furrowing/scripts/axisymmetric_prestress.py').read_bytes()).hexdigest(),'base_source_sha256':hashlib.sha256(Path('controlled_furrowing/scripts/axisymmetric_mechanics.py').read_bytes()).hexdigest(),'scope':'Synthetic mathematical verification, not experimental evidence'}
Path('controlled_furrowing/evidence/prestress_second_variation.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
