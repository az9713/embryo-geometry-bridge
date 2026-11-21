from pathlib import Path
import json,hashlib,sys,platform
import numpy as np,torch,scipy
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from robust_quadratic import Parameters,build_model,_vector_hash
from run_activated_mesh_study import FROZEN
root=Path('controlled_furrowing');torch.set_num_threads(1);rows=[];chains=[];models={}
for lam in [500.,5000.]:
 expected=None
 for strength in [1.,5.,10.]:
  folder=root/f'runs/activated_mesh/lambda_{lam}_target_{strength}';final=json.loads((folder/'final.json').read_text());m=models.setdefault(lam,build_model(Parameters(**final['parameters'])))
  if expected is None:expected=_vector_hash(np.zeros(len(m.free)))
  good=True
  for attempt in final['robust_continuation']['attempts']:
   good &= expected==attempt['initial_verified_vector_sha256'];p=Path(attempt['result_path']);r=json.loads(p.read_text());vector=np.array(r['free_displacement']);_,g=m.objective(vector,r['activation_strength']);force=float(np.max(abs(g)));cert=m.orientation_certificate(m.coordinates(torch.tensor(vector,dtype=torch.float64)))
   rows.append({'lambda':lam,'target':strength,'attempt':attempt['attempt'],'gain':r['activation_strength'],'accepted':attempt['accepted'],'replayed_force':force,'stored_force_difference':abs(force-r['scaled_free_force_linf']),'positive_geometry':cert['certified_positive'],'field_sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
   if attempt['accepted']:expected=attempt['final_vector_sha256']
  good &= expected==_vector_hash(final['free_displacement']);chains.append({'lambda':lam,'target':strength,'verified_chain':bool(good)})
actual={n:hashlib.sha256((root/'scripts'/n).read_bytes()).hexdigest() for n in FROZEN};comp=json.loads((root/'evidence/activated_mesh_comparison.json').read_text())
checks={'six_targets_reached':len(chains)==6 and all(c['verified_chain'] for c in chains),'all22_attempts_retained':len(rows)==22,'all14_accepted_force_geometry_pass':sum(r['accepted'] for r in rows)==14 and all(r['replayed_force']<1e-7 and r['positive_geometry'] for r in rows if r['accepted']),'all8_rejections_retained':sum(not r['accepted'] for r in rows)==8,'stored_forces_replay':all(r['stored_force_difference']<1e-12 for r in rows),'frozen_sources_unchanged':actual==FROZEN,'all12_curve_widths_achieved':len(comp['comparisons'])==6 and all(r[s+'_Eulerian_comparison']['requested_width_achieved'] for r in comp['comparisons'] for s in ['top','bottom'])}
receipt={'status':'PASS_NUMERICAL_EXECUTION_NOT_MESH_ADEQUACY' if all(checks.values()) else 'FAIL','checks':checks,'attempts':rows,'chains':chains,'frozen_sources':actual,'runner_source_sha256':hashlib.sha256((root/'scripts/run_activated_mesh_study.py').read_bytes()).hexdigest(),'software':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'torch':torch.__version__},'scope':'FixedQ0 integration. Replays discrete force/geometry and lastverifiedhashchains; nonnegligiblemeshchangesremain. No PDEerrororbiologicalaccuracycertificate.'}
(root/'evidence/activated_mesh_verification.json').write_text(json.dumps(receipt,indent=2));print(json.dumps({'status':receipt['status'],'checks':checks},indent=2))
