from pathlib import Path
import sys,json,hashlib
import numpy as np,torch
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from sparse_quadratic import AxisymmetricModel,Parameters,HighModel,source_hashes
from robust_quadratic import _vector_hash
root=Path('controlled_furrowing');torch.set_num_threads(1);out=root/'runs/sparse_quadratic/24x12_lambda5000';initial=json.loads((out/'gain_1.0/final.json').read_text());m=AxisymmetricModel(Parameters(**initial['parameters']));parent=HighModel(Parameters(**initial['parameters']));expected=_vector_hash(np.zeros(len(m.free)));rows=[];chains=[];policy={'nonzero_shift_trials':0,'accepted_newton_steps':0,'maximum_linear_residual':0.}
for gain in [1.,5.,10.]:
 final=json.loads((out/f'gain_{gain}/final.json').read_text());chain=True
 for a in final['robust_continuation']['attempts']:
  chain &= expected==a['initial_verified_vector_sha256'];p=Path(a['result_path']);r=json.loads(p.read_text());v=np.array(r['free_displacement']);g=r['activation_strength'];e,gradient=m.objective(v,g);ee,gg=parent.objective(v,g);cert=m.orientation_certificate(m.coordinates(torch.tensor(v,dtype=torch.float64)));force=float(np.max(abs(gradient)))
  rows.append({'requested_target':gain,'attempt':a['attempt'],'gain':g,'accepted':a['accepted'],'force_replayed':force,'stored_force_difference':abs(force-r['scaled_free_force_linf']),'positive_geometry':cert['certified_positive'],'energy_gradient_bit_identical':bool(e==ee and np.array_equal(gradient,gg)),'field_sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
  for stage in r['continuation']:
   for n in stage['newton_polish']:
    policy['accepted_newton_steps']+=int(n['accepted'])
    for trial in n['sparse_linear_trials']:
     policy['nonzero_shift_trials']+=int(trial['diagonal_shift']>0)
     policy['maximum_linear_residual']=max(policy['maximum_linear_residual'],trial.get('linear_residual_linf',0.))
  if a['accepted']:expected=a['final_vector_sha256']
 chain &= expected==_vector_hash(final['free_displacement']);chains.append({'target':gain,'reached':final['reached_target'],'verified_chain':bool(chain),'work_balance_residual':final['accepted_path_work_ledger']['balance_residual']})
verified=json.loads((root/'evidence/sparse_quadratic_verification.json').read_text());checks={'all_three_targets_reached':all(x['reached'] for x in chains),'all11_attempts_retained':len(rows)==11,'all7_accepted_force_geometry_pass':sum(x['accepted'] for x in rows)==7 and all(x['force_replayed']<1e-7 and x['positive_geometry'] for x in rows if x['accepted']),'all4_rejected_attempts_retained':sum(not x['accepted'] for x in rows)==4,'all_energy_gradient_bit_identical':all(x['energy_gradient_bit_identical'] for x in rows),'stored_forces_match_replay':all(x['stored_force_difference']<1e-12 for x in rows),'last_verified_hash_chains':all(x['verified_chain'] for x in chains),'work_balance':all(abs(x['work_balance_residual'])<1e-9 for x in chains),'frozen_parent_sources_unchanged':source_hashes()==verified['parents']}
receipt={'status':'PASS_COMPLETED_DISCRETE_24_PATH' if all(checks.values()) else 'FAIL','checks':checks,'attempts':rows,'chains':chains,'recorded_Newton_policy':policy,'source_sha256':hashlib.sha256((root/'scripts/sparse_quadratic.py').read_bytes()).hexdigest(),'parents':source_hashes(),'scope':'Allactualaccepted/rejected24fields replayed; identicalphysicalenergy/gradient, differentnumericalpolicy. No PDEerror, stabilityorbiologicalpredictioncertificate.'}
(root/'evidence/sparse_quadratic_24_verification.json').write_text(json.dumps(receipt,indent=2));print(json.dumps({'status':receipt['status'],'checks':checks,'policy':policy},indent=2))
