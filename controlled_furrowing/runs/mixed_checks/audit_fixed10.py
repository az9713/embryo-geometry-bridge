from pathlib import Path
import sys,json,hashlib
import numpy as np
import torch
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from axisymmetric_mixed import Parameters,build_model
torch.set_num_threads(1);root=Path('controlled_furrowing');rows=[]
for path in sorted((root/'runs/mixed_checks/fixed10').glob('*.json')):
 r=json.loads(path.read_text());m=build_model(Parameters(**r['parameters']));u=np.array(r['free_displacement']);c=m.energy_components(m.coordinates(torch.tensor(u,dtype=torch.float64)),10.)
 full=float((m.integration*m.parameters.bulk/2*c['logJ']**2).sum());projected=float(c['projected_volumetric']);kernel=full-projected
 _,g=m.objective(u,10.);finalforce=float(np.max(abs(g)))
 failed=[{'gain':s['activation_strength'],'force':s['scaled_gradient_linf']} for s in r['continuation'] if not s['stationarity_verified']]
 rows.append({'file':str(path),'lambda':m.parameters.bulk,'nr':m.parameters.nr,'nz':m.parameters.nz,'final_force_replayed':finalforce,'stored_force_difference':abs(finalforce-r['scaled_free_force_linf']),'strict_failed_stages':failed,'pointwise_logJ_RMS':np.sqrt(2*full/(m.parameters.bulk*m.reference_volume)),'projected_logJ_RMS':np.sqrt(2*projected/(m.parameters.bulk*m.reference_volume)),'unresolved_logJ_RMS':np.sqrt(max(0,2*kernel/(m.parameters.bulk*m.reference_volume))),'unpenalized_kernel_energy_at_same_lambda':kernel})
receipt={'status':'REPLAYED_ALL9_FINAL_FIELDS','all_final_forces_below_1e_minus7':all(q['final_force_replayed']<1e-7 for q in rows),'all_strict_continuation_gates_pass':all(not q['strict_failed_stages'] for q in rows),'rows':rows,'source_sha256':hashlib.sha256(Path('controlled_furrowing/scripts/axisymmetric_mixed.py').read_bytes()).hexdigest(),'scope':'Retains failed prior stages even if final state stationary; projection kernel not governed by lambda penalty'}
(root/'evidence/mixed_kernel_audit.json').write_text(json.dumps(receipt,indent=2));print(json.dumps({'all_final_stationary':receipt['all_final_forces_below_1e_minus7'],'all_strict_continuation_pass':receipt['all_strict_continuation_gates_pass'],'last':rows[-1]},indent=2))

