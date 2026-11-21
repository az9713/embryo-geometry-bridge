from pathlib import Path
import sys,time,json,hashlib
import numpy as np
import torch
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from axisymmetric_mechanics import build_model,Parameters
torch.set_num_threads(1)
m=build_model(Parameters(nr=8,nz=4,patch_radius=9.35512037/33.5,activation_z_min=(33.5-25)/33.5,activation_z_max=(33.5-5)/33.5,basal_area_modulus=.2))
checks=[];last=None
for gain in [.2,.5,1.,2.,4.]:
 t=time.monotonic();r=m.relax(gain,last,maxiter=100);last=r['free_displacement'];checks.append({'gain':gain,'seconds':time.monotonic()-t,'force':r['scaled_free_force_linf'],'converged':r['converged'],'minJ':r['minimum_J'],'stages':len(r['continuation'])})
previous=json.loads(Path('controlled_furrowing/runs/axisymmetric/warm_gain_4.0.json').read_text())
error=float(np.max(abs(np.asarray(previous['displacement'])-np.asarray(r['displacement']))))
receipt={'source_sha256':hashlib.sha256(Path('controlled_furrowing/scripts/axisymmetric_mechanics.py').read_bytes()).hexdigest(),'checks':checks,'previous_default_result_difference':error,'status':'PASS' if all(q['converged'] and q['stages']==1 for q in checks) and error<1e-10 else 'FAIL'}
Path('controlled_furrowing/evidence/axisymmetric_warm_start_final.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
