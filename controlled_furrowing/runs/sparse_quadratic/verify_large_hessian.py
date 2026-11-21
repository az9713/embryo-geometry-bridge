from pathlib import Path
import sys,json,hashlib
import numpy as np,torch
from scipy.sparse import save_npz
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from sparse_quadratic import AxisymmetricModel,Parameters,IntegrationPolicy
root=Path('controlled_furrowing');path=root/'runs/sparse_quadratic/24x12_lambda5000/gain_1.0/attempt_001_gain_0.5.json';r=json.loads(path.read_text());assert r['converged'];torch.set_num_threads(1);m=AxisymmetricModel(Parameters(**r['parameters']),integration_policy=IntegrationPolicy(8,12));v=np.array(r['free_displacement']);H,meta=m.recover_hessian(v,.5);rng=np.random.default_rng(191124);rows=[]
for k in range(3):
 d=rng.normal(size=len(v));d/=np.linalg.norm(d);u=torch.tensor(v,dtype=torch.float64,requires_grad=True);e=m.energy_components(m.coordinates(u),.5)['total']/m.energy_scale;g=torch.autograd.grad(e,u,create_graph=True)[0];ad=torch.autograd.grad(g,u,grad_outputs=torch.tensor(d,dtype=torch.float64))[0].detach().numpy();s=H@d;h=1e-6;fd=(m.objective(v+h*d,.5)[1]-m.objective(v-h*d,.5)[1])/(2*h);scale=max(1.,float(np.max(abs(ad))))
 rows.append({'direction':k,'AD_relative_linf':float(np.max(abs(s-ad)))/scale,'FD_relative_linf':float(np.max(abs(s-fd)))/max(1.,float(np.max(abs(fd))))})
checks={'large_pattern_row_colors_unique':m.pattern_receipt['row_color_collision_free'],'three_AD_Hv_checks':all(q['AD_relative_linf']<1e-11 for q in rows),'three_FD_Hv_checks':all(q['FD_relative_linf']<2e-5 for q in rows)}
save_npz(root/'runs/sparse_quadratic/hessian24_at_verified_gain0p5.npz',H)
receipt={'status':'PASS_LARGE_DERIVATIVE_SPOTCHECK' if all(checks.values()) else 'FAIL','checks':checks,'rows':rows,'pattern':m.pattern_receipt,'input':str(path.relative_to(root)),'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256((root/'scripts/sparse_quadratic.py').read_bytes()).hexdigest(),'scope':'Derivativeevaluationonlyattheexistingverified24x12gain0.5field; nooptimizationjob, noPDEerrororbiologyclaim.'}
(root/'evidence/sparse_quadratic_24_hvp_check.json').write_text(json.dumps(receipt,indent=2));print(json.dumps({'status':receipt['status'],'rows':rows},indent=2))
