from pathlib import Path
import sys,json,hashlib,time
import numpy as np,torch
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from quadratic_quadrature_study import AxisymmetricModel,Parameters,IntegrationPolicy,evaluate,exact_geometry,source_hashes
root=Path('controlled_furrowing');torch.set_num_threads(1);rows=[]
for gain in (5.,10.):
 path=root/f'runs/sparse_quadratic/24x12_lambda5000/gain_{gain}/final.json';state=json.loads(path.read_text());start=time.perf_counter()
 model=AxisymmetricModel(Parameters(**state['parameters']),integration_policy=IntegrationPolicy(12,16));cross=evaluate(model,state)
 row={'gain':gain,'field':path.relative_to(root).as_posix(),'field_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'Duffy12_surface16':cross,'cross_force_gate':cross['scaled_force_linf']<1e-7,'seconds':time.perf_counter()-start};rows.append(row)
 print(json.dumps(row),flush=True)
 if gain==10.:
  geom=exact_geometry(state);(root/'evidence/sparse_quadratic_24_final_geometry.json').write_text(json.dumps({'field_sha256':row['field_sha256'],'exact_geometry':geom,'scope':'Exact rational certificate of stored binary P2 coordinates, not PDE or biological error.'},indent=2),encoding='utf-8')
  row['geometry_summary']={k:v for k,v in geom.items() if k not in ('triangle_certificates','boundary','strong_monotonicity')}
receipt={'status':'PASS' if all(r['cross_force_gate'] for r in rows) and geom['all_exact_J_positive'] and geom['global_embedding_certified'] else 'FAIL','rows':rows,'source_hashes':source_hashes(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Evaluation-only final24 gains5/10, no new equilibrium or biological data.'}
(root/'evidence/sparse_quadratic_24_crosscheck.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(json.dumps(receipt,indent=2),flush=True)
