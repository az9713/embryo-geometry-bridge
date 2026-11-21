from pathlib import Path
from dataclasses import replace
from copy import deepcopy
import sys,json,hashlib
import numpy as np
import torch
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from robust_quadratic import Parameters,SolvePolicy,build_model,AxisymmetricModel,_vector_hash
from axisymmetric_quadratic import AxisymmetricModel as Frozen
root=Path('controlled_furrowing');out=root/'runs/quadratic_robust';torch.set_num_threads(1)
checks={};m=build_model(Parameters(nr=6,nz=3,basal_surface_tension=2.));zero=m._attempt(0.,np.zeros(len(m.free)))
for tag,change in [('nonfinite_gain',{'activation_strength':float('nan')}),('failed_previous',{'converged':False}),('nonfinite_vector',{'free_displacement':[float('nan')]*len(m.free)})]:
 r={**zero,**change}
 try:m._verify_state(r);checks[tag+'_rejected']=False
 except ValueError:checks[tag+'_rejected']=True
original=m.objective;m.objective=lambda vector,strength,gradient=True:(0.,np.full(len(m.free),np.nan))
try:m._verify_state(zero);checks['nonfinite_gradient_rejected']=False
except ValueError:checks['nonfinite_gradient_rejected']=True
m.objective=original
class FailOnce(AxisymmetricModel):
 def __init__(self,p,policy):super().__init__(p,policy);self.injected=False
 def _attempt(self,strength,initial):
  r=super()._attempt(strength,initial)
  if strength>0 and not self.injected:
   self.injected=True;r['converged']=False;r['free_displacement']=(np.asarray(r['free_displacement'])+99).tolist();r['injected_failure_for_controller_unit_test']=True
  return r
p=Parameters(nr=6,nz=3,basal_surface_tension=2.)
s=FailOnce(p,SolvePolicy());r=s.advance(.2,output_directory=out/'guard_tests/bisection');a=r['robust_continuation']['attempts']
checks['injected_failure_bisected']=r['reached_target'] and len(a)==3 and [q['accepted'] for q in a]==[False,True,True]
checks['failed_vector_never_reused']=a[1]['initial_verified_vector_sha256']==a[0]['initial_verified_vector_sha256'] and a[1]['initial_verified_vector_sha256']!=a[0]['final_vector_sha256'] and a[2]['initial_verified_vector_sha256']==a[1]['final_vector_sha256']
s=FailOnce(p,SolvePolicy(minimum_step=.2));r=s.advance(.2,output_directory=out/'guard_tests/minstep');checks['minimum_step_stops_without_target']=not r['reached_target'] and not r['converged'] and r['activation_strength']==0 and len(r['robust_continuation']['attempts'])==1
s=FailOnce(p,SolvePolicy(maximum_attempts=1));r=s.advance(.2,output_directory=out/'guard_tests/maxattempt');checks['attempt_limit_stops_without_target']=not r['reached_target'] and not r['converged'] and r['activation_strength']==0 and len(r['robust_continuation']['attempts'])==1
rows=[]
for lam in [500.,5000.]:
 previous_hash=None
 for target in [1.,5.,10.]:
  path=out/f'lambda_{lam}_target_{target}'/'final.json';r=json.loads(path.read_text());model=build_model(Parameters(**r['parameters']));v,gain=model._verify_state(r);frozen=Frozen(Parameters(**r['parameters']));e,gradient=model.objective(v,gain);ef,gf=frozen.objective(v,gain)
  attempts=r['robust_continuation']['attempts'];expected=previous_hash or _vector_hash(np.zeros(len(v)));chain=True
  for a in attempts:
   chain &= a['initial_verified_vector_sha256']==expected
   if a['accepted']:expected=a['final_vector_sha256']
  actual=_vector_hash(v);chain &= expected==actual;previous_hash=actual
  rows.append({'lambda':lam,'target':target,'state_reverified':True,'input_hash_chain_verified':bool(chain),'energy_gradient_bit_identical':bool(e==ef and np.array_equal(gradient,gf)),'work_balance':r['accepted_path_work_ledger']['balance_residual'],'failed_attempts':sum(not a['accepted'] for a in attempts)})
checks['all6_scientific_states_reverified']=len(rows)==6 and all(q['state_reverified'] and q['input_hash_chain_verified'] and q['energy_gradient_bit_identical'] and abs(q['work_balance'])<1e-9 for q in rows)
receipt={'status':'PASS' if all(checks.values()) else 'FAIL','checks':{k:bool(v) for k,v in checks.items()},'rows':rows,'final_controller_sha256':hashlib.sha256(Path('controlled_furrowing/scripts/robust_quadratic.py').read_bytes()).hexdigest(),'numerical_batch_controller_sha256':hashlib.sha256((out/'pre_guard_hardening/robust_quadratic.py').read_bytes()).hexdigest(),'scope':'Final input guards and failure-controller unit tests; all6 prior numerical states replayed. Failure injections explicitly synthetic, not scientific runs. Guards/terminal ledger persistence do not alter finite-input energy/optimizer arithmetic.'}
(root/'evidence/quadratic_robust_guard_review.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
