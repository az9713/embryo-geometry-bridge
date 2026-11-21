from pathlib import Path
import sys,json,hashlib,time
import numpy as np,torch
from scipy.linalg import ldl
from scipy.sparse.linalg import eigsh
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from sparse_quadratic import AxisymmetricModel,Parameters
root=Path('controlled_furrowing');torch.set_num_threads(1);rows=[]
protocol={'scope':'Floating arithmetic Hessian curvature/inertia diagnostic only, not a rigorous true Hessian or PDE stability certificate; no optimization or target fit.','threshold':'100*machine_epsilon*matrix_infinity_norm for inertia zero classification; eigenpair residual<1e-8 relative.','inputs':['runs/joint_refinement/final.json','runs/sparse_quadratic/24x12_lambda5000/gain_10.0/final.json'],'eigensolver':'Sparse shift-invert near zero; numerical dense symmetric LDL inertia independently determines sign counts.','force_tolerance':1e-7,'source_sha256':hashlib.sha256((root/'scripts/sparse_quadratic.py').read_bytes()).hexdigest()}
(root/'evidence/sparse_quadratic_stability_protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
for tag,rel in zip(['16x8','24x12'],protocol['inputs']):
 start=time.perf_counter();path=root/rel;state=json.loads(path.read_text());model=AxisymmetricModel(Parameters(**state['parameters']));u=np.array(state['free_displacement']);gain=state['activation_strength'];H,meta=model.recover_hessian(u,gain);A=H.toarray();scale=np.linalg.norm(A,np.inf);tolerance=100*np.finfo(float).eps*scale
 lu,D,perm=ldl(A,lower=True,hermitian=True,check_finite=True);blocks=[];i=0
 while i<len(D):
  if i+1<len(D) and D[i,i+1]!=0.:
   vals=np.linalg.eigvalsh(D[i:i+2,i:i+2]);i+=2
  else:vals=np.array([D[i,i]]);i+=1
  blocks.extend(vals.tolist())
 blocks=np.array(blocks);inertia={'negative':int(sum(blocks < -tolerance)),'near_zero':int(sum(abs(blocks)<=tolerance)),'positive':int(sum(blocks>tolerance))}
 reconstruction=float(np.linalg.norm(A-lu@D@lu.T,np.inf)/scale)
 vals,vecs=eigsh(H,k=3,sigma=0.,which='LM',tol=1e-11,maxiter=5000,v0=np.random.default_rng(1911).normal(size=len(u)));order=np.argsort(vals);vals=vals[order];vecs=vecs[:,order]
 residuals=[float(np.linalg.norm(H@vecs[:,i]-vals[i]*vecs[:,i])/max(np.linalg.norm(H@vecs[:,i]),abs(vals[i]),1e-30)) for i in range(3)]
 # Independent AD product of actual frozen physical potential at smallest reported pair.
 x=torch.tensor(u,dtype=torch.float64,requires_grad=True);e=model.energy_components(model.coordinates(x),gain)['total']/model.energy_scale;g=torch.autograd.grad(e,x,create_graph=True)[0];v=vecs[:,0];hv=torch.autograd.grad(g,x,grad_outputs=torch.tensor(v,dtype=torch.float64))[0].detach().numpy();ad_res=float(np.linalg.norm(hv-vals[0]*v)/max(np.linalg.norm(hv),abs(vals[0]),1e-30))
 energy_checks=[]
 if inertia['negative'] or vals[0]<-tolerance:
  # If near-zero eigenvectors are positive but inertia is negative, explicitly seek algebraic smallest.
  if vals[0]>=-tolerance:
   negvalues,negvectors=eigsh(H,k=1,which='SA',tol=1e-9,maxiter=20000,v0=v);direction=negvectors[:,0];negative_value=float(negvalues[0])
  else:direction=v;negative_value=float(vals[0])
  initial=model.objective(u,gain,False)
  for h in (1e-3,1e-4,1e-5):
   for sign in (-1,1):
    candidate=u+sign*h*direction;cert=model.orientation_certificate(model.coordinates(torch.tensor(candidate,dtype=torch.float64)));value=model.objective(candidate,gain,False)
    energy_checks.append({'step':sign*h,'energy_difference':value-initial,'positive_geometry':cert['certified_positive'],'rayleigh':negative_value})
 else:negative_value=None
 np.savez_compressed(root/f'runs/sparse_quadratic/stability_{tag}.npz',data=H.data,indices=H.indices,indptr=H.indptr,shape=H.shape,eigenvalues=vals,eigenvectors=vecs,LDL_D=D,LDL_permutation=perm,LDL_block_eigenvalues=blocks)
 row={'mesh':tag,'input':rel,'field_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'seconds':time.perf_counter()-start,'matrix_infinity_norm':float(scale),'LDL_sign_tolerance':float(tolerance),'LDL_inertia':inertia,'LDL_relative_reconstruction_residual':reconstruction,'near_zero_eigenvalues':vals.tolist(),'eigenpair_relative_residuals':residuals,'independent_AD_smallest_pair_relative_residual':ad_res,'raw_hessian_asymmetry_linf':meta['raw_asymmetry_linf'],'force_linf':float(np.max(abs(meta['gradient']))),'negative_direction_energy_replays':energy_checks,'interpretation':'Numerically positive definite discrete Hessian at this saved field.' if inertia['positive']==len(u) and vals[0]>tolerance else 'Nonpositive or numerically unresolved curvature; see inertia and energy replays.'};rows.append(row)
 print(json.dumps(row),flush=True)
 receipt={'protocol':protocol,'rows':rows,'status':'COMPLETE' if len(rows)==2 else 'PARTIAL','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Numerical discrete equilibrium curvature diagnostic only. No rigorous stability, continuum error, global minimum or branch uniqueness theorem.'};(root/'evidence/sparse_equilibrium_stability_diagnostic.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
