from pathlib import Path
import sys,json,hashlib,time
import numpy as np,torch
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from sparse_quadratic import AxisymmetricModel,Parameters
root=Path('controlled_furrowing');path=root/'runs/sparse_quadratic/24x12_lambda5000/gain_1.0/attempt_001_gain_0.5.json';r=json.loads(path.read_text());torch.set_num_threads(1);m=AxisymmetricModel(Parameters(**r['parameters']));v=np.array(r['free_displacement']);m.objective(v,.5)
with torch.profiler.profile(activities=[torch.profiler.ProfilerActivity.CPU],record_shapes=False,profile_memory=False) as prof:
 m.objective(v,.5)
items=sorted(prof.key_averages(),key=lambda x:x.self_cpu_time_total,reverse=True);total=sum(x.self_cpu_time_total for x in items);rows=[{'operator':x.key,'self_CPU_us':x.self_cpu_time_total,'CPU_total_us':x.cpu_time_total,'count':x.count,'fraction_of_profiled_self_CPU':x.self_cpu_time_total/total} for x in items[:20]]
receipt={'scope':'One warmedCPUenergy+gradient evaluation at existingverified24gain0.5; profiling only, no optimization/modelchange; overheadnotwholejobtiming','self_CPU_us_total':total,'top_operators':rows,'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256((root/'scripts/sparse_quadratic.py').read_bytes()).hexdigest()}
(root/'evidence/sparse_quadratic_gradient_profile.json').write_text(json.dumps(receipt,indent=2));print(json.dumps({'total_profiled_CPU_ms':total/1000,'top5':rows[:5]},indent=2))
