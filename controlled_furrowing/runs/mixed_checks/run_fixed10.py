from pathlib import Path
from dataclasses import asdict
import sys,time,json,hashlib
import torch
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from axisymmetric_mixed import build_model,Parameters
root=Path('controlled_furrowing');out=root/'runs/mixed_checks/fixed10';out.mkdir(parents=True,exist_ok=True);torch.set_num_threads(1)
protocol={'scope':'Synthetic fixedactivation10 comparison; same material geometry as original no-tension P1 study, changed finite-element/projected-volume formulation','grids':[[8,4],[12,6],[16,8]],'lambda':[5.,50.,500.],'gain':10.,'continuation_steps':10,'maxiter_lbfgs':100,'parameters':{'height':1.,'outer_radius':3.,'patch_radius':9.35512037/33.5,'activation_z_min':1-25/33.5,'activation_z_max':1-5/33.5,'mu':1.,'basal_area_modulus':2.,'apical_area_modulus':.2,'basal_surface_tension':0.,'apical_surface_tension':0.},'report_rules':['No scientific PASS from solver stationarity alone','Bernstein uncertified not inversion','Pointwise quadrature compression distinct from conservative coefficient bound','No locking isolation or continuum accuracy without matched manufactured forced problem']}
(root/'evidence/mixed_fixed10_protocol.json').write_text(json.dumps(protocol,indent=2));rows=[]
for nr,nz in protocol['grids']:
 for lam in protocol['lambda']:
  t=time.monotonic();m=build_model(Parameters(nr=nr,nz=nz,bulk=lam,**protocol['parameters']))
  r=m.solve(10.,continuation_steps=10,maxiter=100)
  path=out/f'mesh{nr}x{nz}_lambda{lam}.json';path.write_text(json.dumps(r,indent=2))
  row={'nr':nr,'nz':nz,'lambda':lam,'seconds':time.monotonic()-t,'file':str(path),**{k:r[k] for k in ['converged','stationary_discrete_equations','scaled_free_force_linf','minimum_J','maximum_J','sampled_minimum_J','sampled_maximum_J','relative_volume_change','pressure_constraint_linf','maximum_activated_column_shortening','central_apical_displacement','central_basal_displacement']},'geometry_certified':r['orientation_certificate']['certified_positive'],'simple_boundary':r['surface_certificates']['sufficient_simple_meridian_boundary']}
  rows.append(row);print(json.dumps(row),flush=True)
  receipt={'status':'COMPLETED_NUMERICAL_SENSITIVITY_NOT_CONTINUUM_CERTIFICATE','rows':rows,'source_sha256':hashlib.sha256(Path('controlled_furrowing/scripts/axisymmetric_mixed.py').read_bytes()).hexdigest(),'complete':len(rows)==9,'protocol':protocol}
  (root/'evidence/mixed_fixed10.json').write_text(json.dumps(receipt,indent=2))
