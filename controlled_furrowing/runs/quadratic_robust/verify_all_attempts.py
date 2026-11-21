from pathlib import Path
import json,hashlib,sys
import numpy as np
import torch
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from robust_quadratic import Parameters,build_model
from axisymmetric_quadratic import AxisymmetricModel as Frozen
root=Path('controlled_furrowing');torch.set_num_threads(1);rows=[]
for folder in sorted((root/'runs/quadratic_robust').glob('lambda_*')):
 ledger=json.loads((folder/'attempt_ledger.json').read_text())
 for a in ledger['attempts']:
  path=Path(a['result_path']);r=json.loads(path.read_text());m=build_model(Parameters(**r['parameters']));v=np.array(r['free_displacement']);gain=r['activation_strength'];e,g=m.objective(v,gain);ef,gf=Frozen(Parameters(**r['parameters'])).objective(v,gain);cert=m.orientation_certificate(m.coordinates(torch.tensor(v,dtype=torch.float64)));force=float(np.max(abs(g)))
  rows.append({'target_folder':folder.name,'attempt':a['attempt'],'gain':gain,'accepted':a['accepted'],'force_replayed':force,'stored_force_error':abs(force-r['scaled_free_force_linf']),'geometry_certified':cert['certified_positive'],'energy_gradient_bit_identical':bool(e==ef and np.array_equal(g,gf))})
checks={'all18_actual_attempts_present':len(rows)==18,'all12_accepted_force_geometry_pass':sum(q['accepted'] for q in rows)==12 and all(q['force_replayed']<1e-7 and q['geometry_certified'] for q in rows if q['accepted']),'all6_failed_attempts_retained':sum(not q['accepted'] for q in rows)==6,'stored_force_matches_replay':all(q['stored_force_error']<1e-12 for q in rows),'all_attempt_energies_gradients_bit_identical':all(q['energy_gradient_bit_identical'] for q in rows)}
receipt={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'rows':rows,'controller_sha256':hashlib.sha256(Path('controlled_furrowing/scripts/robust_quadratic.py').read_bytes()).hexdigest(),'scope':'Actualscientificnumericalattemptsonly; excludesclearlylabeledcontrollerinjectiontests; inputhashchainverifiedseparately.'}
(root/'evidence/quadratic_robust_attempt_replay.json').write_text(json.dumps(receipt,indent=2));print(json.dumps({'status':receipt['status'],'checks':checks},indent=2))
