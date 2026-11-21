"""Freeze an actual source-M saved-state cache; no equilibrium optimization."""
from pathlib import Path
from fractions import Fraction as F
import datetime,hashlib,json,sys
import numpy as np
from scipy.sparse.linalg import splu
import torch

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
import algebraic_root_feasibility as core
from material_tangent_inference import mechanical_tangent,check_sources

sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
save=lambda p,v:Path(p).write_text(json.dumps(v,indent=2)+'\n',encoding='utf8')

if(OUT/'point_cache.json').exists():raise RuntimeError('preserve existing cache; no duplicate regeneration')
protocol=read(ROOT/'runs/material_cohort_direct_v1/protocol.json')
result=read(ROOT/'runs/material_cohort_direct_v1/results.json')
assert result['status']=='COMPLETE_PRIMARY_MATERIAL_STEP_DIRECT_REPLAY'
source=next(x for x in protocol['records']if x['record_id']=='optoGEF2_small_col_M')
rec=next(x for x in result['records']if x['record_id']==source['record_id'])
assert sha(ROOT/rec['file'])==rec['sha256']
state=read(ROOT/rec['file']);assert state['parameters']==protocol['new_parameters']
preparation={'status':'FROZEN_BEFORE_LOCAL_HESSIAN_EVALUATION','source_sha256':sha(__file__),
 'direct_protocol_sha256':sha(ROOT/'runs/material_cohort_direct_v1/protocol.json'),
 'direct_results_sha256':sha(ROOT/'runs/material_cohort_direct_v1/results.json'),
 'direct_review_sha256':sha(ROOT/'evidence/material_cohort_direct_independent_review.json'),
 'record_id':source['record_id'],'dose_exact':source['dose_exact'],'base_gain_exact':protocol['new_gain_exact'],
 'requested_strength_exact':source['requested_strength_exact'],'solver_strength_binary_exact':source['solver_input_strength_exact'],
 'input_state':rec['file'],'input_state_sha256':rec['sha256'],
 'physical_sources':check_sources(),'tangent_source_sha256':sha(ROOT/'scripts/material_tangent_inference.py'),
 'scope':'Stored-state Hessian/preconditioner/predictor preparation only; no new optimization, root certificate or fullparameter fit'}
save(OUT/'preparation_protocol.json',preparation)
torch.set_num_threads(1)
model,u,H,B,S,verification=mechanical_tangent(state)
R=splu(H.tocsc()).solve(np.eye(len(u)))
np.savez_compressed(OUT/'tangents.npz',free_vector=u,free_indices=model.free,
 H_data=H.data,H_indices=H.indices,H_indptr=H.indptr,H_shape=np.array(H.shape),
 parameter_force_derivatives=B,displacement_tangents=S)
save(OUT/'tangent_verification.json',verification)
packet={'status':'FROZEN_SAVED_STATE_CACHE_NO_ROOT_CLAIM','parameters':state['parameters'],
 'activation_strength_binary_exact':str(F.from_float(float(state['activation_strength']))),
 'input_state':rec['file'],'input_state_sha256':rec['sha256'],
 'matrix_archive':(OUT/'tangents.npz').relative_to(ROOT).as_posix(),'matrix_archive_sha256':sha(OUT/'tangents.npz'),
 'frozen_binary_cache':core.freeze_cache(model,u,R),'preparation_protocol_sha256':sha(OUT/'preparation_protocol.json'),
 'no_new_equilibrium_optimization':True,
 'uniform_target_warning':'Actual cell must use exact common_gain*dose, not inherit binarysolverstrength as its parameterlaw'}
assert all(verification['checks'].values());save(OUT/'point_cache.json',packet)
print(json.dumps({'status':packet['status'],'sha256':sha(OUT/'point_cache.json'),
 'record':source['record_id'],'dose':source['dose_exact'],'strength':state['activation_strength'],
 'dof':len(u),'no_new_equilibrium_optimization':True}),flush=True)
