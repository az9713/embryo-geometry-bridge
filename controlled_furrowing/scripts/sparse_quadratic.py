"""Sparse AD Hessian recovery and damped Newton for the fixed high-Q energy.

No energy/gradient override. Structural element cliques + column coloring recover
the Hessian from chunked batched AD Hessian-vector products. Solver policy differs
from dense spectral regularization; no identical-trajectory claim is made.
"""
from dataclasses import dataclass,asdict
from pathlib import Path
import argparse,hashlib,json,os,platform,subprocess,sys,threading,time
import numpy as np
import psutil
import scipy
from scipy.sparse import csr_matrix,eye
from scipy.sparse.linalg import splu
import torch
from quadratic_quadrature_study import AxisymmetricModel as HighModel,IntegrationPolicy,source_hashes
from robust_quadratic import Parameters,SolvePolicy

ROOT=Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class SparsePolicy:
    seed_chunk:int=4
    maximum_shift_trials:int=12
    line_search_steps:int=30
    armijo:float=1e-4
    relative_shift_floor:float=1e-9
    maximum_newton_displacement:float=.25

    def __post_init__(self):
        if self.seed_chunk<1 or self.maximum_shift_trials<1 or self.line_search_steps<1:raise ValueError('positive sparse-policy budgets required')
        if not 0<self.armijo<1 or self.relative_shift_floor<=0 or self.maximum_newton_displacement<=0:raise ValueError('valid positive descent-policy constants required')


def measured(call):
    process=psutil.Process();before=process.memory_info().rss;peak=[before];stop=threading.Event()
    def sample():
        while not stop.wait(.005):peak[0]=max(peak[0],process.memory_info().rss)
    thread=threading.Thread(target=sample,daemon=True);thread.start();start=time.perf_counter()
    try:value=call()
    finally:
        peak[0]=max(peak[0],process.memory_info().rss);stop.set();thread.join()
    return value,{'seconds':time.perf_counter()-start,'RSS_before_bytes':before,'sampled_peak_RSS_bytes':peak[0],'sampled_incremental_RSS_bytes':peak[0]-before,'measurement':'ProcessRSSsampled5ms; includesnativeTorchallocations; maymissbriefpeaks; notGPUmemory.'}


class AxisymmetricModel(HighModel):
    def __init__(self,parameters=Parameters(),policy=SolvePolicy(),integration_policy=IntegrationPolicy(8,12),sparse_policy=SparsePolicy()):
        self.sparse_policy=sparse_policy
        super().__init__(parameters,policy,integration_policy)
        n=len(self.free);mapping=np.full(self.X.numel(),-1,int);mapping[self.free]=np.arange(n)
        rowsets=[{i} for i in range(n)]
        elements=list(self.triangles)+list(self.top_edges)+list(self.bottom_edges)
        for nodes in elements:
            local=mapping[np.array([[2*int(k),2*int(k)+1] for k in nodes]).ravel()];local=local[local>=0]
            for i in local:rowsets[i].update(local.tolist())
        self.pattern_rows=np.concatenate([np.full(len(row),i,int) for i,row in enumerate(rowsets)])
        self.pattern_cols=np.concatenate([np.array(sorted(row),int) for row in rowsets])
        self.pattern_indptr=np.r_[0,np.cumsum([len(row) for row in rowsets])]
        conflicts=[set() for _ in range(n)]
        for row in rowsets:
            for j in row:conflicts[j].update(row)
        colors=np.full(n,-1,int)
        for j in sorted(range(n),key=lambda i:(-len(conflicts[i]),i)):
            occupied={int(colors[k]) for k in conflicts[j] if colors[k]>=0};color=0
            while color in occupied:color+=1
            colors[j]=color
        self.colors=colors;self.color_count=int(colors.max()+1)
        self.color_positions=[np.flatnonzero(colors[self.pattern_cols]==c) for c in range(self.color_count)]
        unique=all(len({int(colors[j]) for j in row})==len(row) for row in rowsets)
        symmetric=all(i in rowsets[j] for i,row in enumerate(rowsets) for j in row)
        if not unique or not symmetric:raise RuntimeError('invalid symmetric Hessian pattern/coloring')
        pattern_hash=hashlib.sha256(self.pattern_rows.astype('<i8').tobytes()+self.pattern_cols.astype('<i8').tobytes()+colors.astype('<i8').tobytes()).hexdigest()
        self.pattern_receipt={'free_DOF':n,'structural_nonzeros':len(self.pattern_cols),'colors':self.color_count,'seed_chunk':sparse_policy.seed_chunk,'row_color_collision_free':unique,'pattern_symmetric':symmetric,'pattern_color_sha256':pattern_hash,'dense_matrix_bytes':8*n*n,'sparse_values_and_int64_indices_bytes':16*len(self.pattern_cols)+8*(n+1),'scope':'Structural12DOFtriangle/6DOFsurfacecliques oflocalhighQpointwiseenergy; nohand-codedconstitutiveHessian.'}

    def recover_hessian(self,vector,strength):
        start=time.perf_counter();u=torch.tensor(np.asarray(vector),dtype=torch.float64,requires_grad=True)
        energy=self.energy_components(self.coordinates(u),strength)['total']/self.energy_scale
        gradient=torch.autograd.grad(energy,u,create_graph=True)[0]
        g=gradient.detach().numpy().copy();data=np.empty(len(self.pattern_cols));n=len(vector)
        for first in range(0,self.color_count,self.sparse_policy.seed_chunk):
            selected=list(range(first,min(first+self.sparse_policy.seed_chunk,self.color_count)))
            seeds=torch.zeros((len(selected),n),dtype=torch.float64)
            for k,color in enumerate(selected):seeds[k,torch.tensor(self.colors==color)]=1.
            hv=torch.autograd.grad(gradient,u,grad_outputs=seeds,is_grads_batched=True,retain_graph=selected[-1]+1<self.color_count)[0].detach().numpy()
            for k,color in enumerate(selected):
                positions=self.color_positions[color];data[positions]=hv[k,self.pattern_rows[positions]]
        raw=csr_matrix((data,self.pattern_cols,self.pattern_indptr),shape=(n,n));skew=raw-raw.T
        asymmetry=float(np.max(abs(skew.data))) if skew.nnz else 0.
        H=(raw+raw.T)*.5
        return H,{'seconds':time.perf_counter()-start,'raw_asymmetry_linf':asymmetry,'energy':float(energy.detach()),'gradient':g,'pattern':self.pattern_receipt}

    def _newton_polish(self,vector,strength,iterations=None):
        history=[];count=self.policy.newton_iterations if iterations is None else iterations
        for _ in range(count):
            old,g=self.objective(vector,strength)
            if np.max(abs(g))<2e-10:break
            H,meta=self.recover_hessian(vector,strength);diag_scale=max(1.,float(np.max(abs(H.diagonal()))));attempts=[];accepted=False
            for trial in range(self.sparse_policy.maximum_shift_trials):
                shift=0. if trial==0 else diag_scale*self.sparse_policy.relative_shift_floor*10.**(trial-1)
                A=H if shift==0 else H+shift*eye(len(vector),format='csr')
                record={'diagonal_shift':shift}
                try:
                    direction=splu(A.tocsc(),permc_spec='COLAMD').solve(-g)
                    dot=float(g@direction);record['gradient_dot_direction']=dot
                    record['linear_residual_linf']=float(np.max(abs(A@direction+g)))
                    if not np.all(np.isfinite(direction)) or dot>=-1e-14*np.linalg.norm(g)*np.linalg.norm(direction):
                        record['reason']='not a verified descent direction';attempts.append(record);continue
                    cap=min(1.,self.sparse_policy.maximum_newton_displacement*self.parameters.height/max(float(np.max(abs(direction))),1e-300))
                    for exponent in range(self.sparse_policy.line_search_steps):
                        fraction=cap*2.**(-exponent);candidate=vector+fraction*direction
                        if any(lo is not None and candidate[i]<lo for i,(lo,_) in enumerate(self.bounds)):continue
                        value=self.objective(candidate,strength,False)
                        if value>old+self.sparse_policy.armijo*fraction*dot+1e-16:continue
                        current=self.coordinates(torch.tensor(candidate,dtype=torch.float64))
                        if not self.orientation_certificate(current)['certified_positive']:continue
                        vector=candidate;accepted=True;record.update({'accepted':True,'line_fraction':fraction,'objective_after':value});break
                    if not accepted:record['reason']='Armijo/geometry backtracking exhausted'
                except (RuntimeError,ValueError,np.linalg.LinAlgError) as error:record['reason']=str(error)
                attempts.append(record)
                if accepted:break
            history.append({'accepted':accepted,'gradient_before_linf':float(np.max(abs(g))),'hessian_seconds':meta['seconds'],'raw_hessian_asymmetry':meta['raw_asymmetry_linf'],'sparse_linear_trials':attempts,'policy_difference':'SparseLUdescent/diagonalshifts/stepcapandpositivegeometrylineacceptance, notdenseeigenvalueregularization; trajectoriesneednotmatch.'})
            if not accepted:break
        return vector,history

    def diagnostics(self,vector,strength,history=None):
        answer=super().diagnostics(vector,strength,history);answer['sparse_solver_policy']=asdict(self.sparse_policy);answer['sparse_hessian_pattern']=self.pattern_receipt
        return answer


def build_model(parameters=None,policy=SolvePolicy(),integration_policy=IntegrationPolicy(8,12),sparse_policy=SparsePolicy(),**kwargs):
    return AxisymmetricModel(Parameters(**kwargs) if parameters is None else parameters,policy,integration_policy,sparse_policy)


def _benchmark(kind):
    m=build_model(Parameters(nr=6,nz=3,bulk=500.,basal_surface_tension=2.));vector=np.random.default_rng(1911).normal(size=len(m.free))*1e-5
    if kind=='sparse':_,metrics=measured(lambda:m.recover_hessian(vector,1.))
    else:
        u=torch.tensor(vector,dtype=torch.float64,requires_grad=True)
        _,metrics=measured(lambda:torch.autograd.functional.hessian(lambda x:m.energy_components(m.coordinates(x),1.)['total']/m.energy_scale,u,vectorize=True))
    result={'kind':kind,**metrics,'pattern':m.pattern_receipt,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Smallmeshseparateprocessbenchmark; not24meshmemorypredictionorphysicalaccuracy.'}
    (ROOT/'evidence'/f'sparse_quadratic_benchmark_{kind}.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


def self_test():
    out=ROOT/'runs/sparse_quadratic';out.mkdir(parents=True,exist_ok=True)
    protocol={'scope':'SparseADHessian/numericalpolicyverification beforeany24x12run','mesh':[6,3],'volume':'Duffy8','surface':'Gauss12','seed_chunk':4,'thresholds':{'entire_H_relative':1e-11,'offpattern_absolute':1e-10,'asymmetry_absolute':1e-10,'Hv_AD_relative':1e-11,'Hv_FD_relative':2e-5,'energy_gradient_identity':'bitforbit'},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (ROOT/'evidence/sparse_quadratic_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    p=Parameters(nr=6,nz=3,bulk=500.,basal_surface_tension=2.);m=build_model(p);parent=HighModel(p,integration_policy=IntegrationPolicy(8,12));rng=np.random.default_rng(1911);vector=rng.normal(size=len(m.free))*1e-5
    H,meta=m.recover_hessian(vector,1.);u=torch.tensor(vector,dtype=torch.float64,requires_grad=True)
    dense=torch.autograd.functional.hessian(lambda x:m.energy_components(m.coordinates(x),1.)['total']/m.energy_scale,u,vectorize=True).detach().numpy()
    delta=float(np.max(abs(H.toarray()-dense)));scale=max(1.,float(np.max(abs(dense))));mask=np.ones(dense.shape,bool);mask[m.pattern_rows,m.pattern_cols]=False
    off=float(np.max(abs(dense[mask]))) if mask.any() else 0.
    checks={'row_colors_collision_free':m.pattern_receipt['row_color_collision_free'],'symmetric_pattern':m.pattern_receipt['pattern_symmetric'],'entire_H_matches_dense_AD':delta/scale<1e-11,'dense_offpattern_zero':off<1e-10,'raw_recovery_symmetric':meta['raw_asymmetry_linf']<1e-10};hvs=[]
    for k in range(3):
        direction=rng.normal(size=len(vector));direction/=np.linalg.norm(direction);uu=torch.tensor(vector,dtype=torch.float64,requires_grad=True)
        energy=m.energy_components(m.coordinates(uu),1.)['total']/m.energy_scale;g=torch.autograd.grad(energy,uu,create_graph=True)[0];ad=torch.autograd.grad(g,uu,grad_outputs=torch.tensor(direction,dtype=torch.float64))[0].detach().numpy();sparse=H@direction
        h=1e-6;fd=(m.objective(vector+h*direction,1.)[1]-m.objective(vector-h*direction,1.)[1])/(2*h)
        err=float(np.max(abs(sparse-ad)))/max(1.,float(np.max(abs(ad))));fd_err=float(np.max(abs(sparse-fd)))/max(1.,float(np.max(abs(fd))))
        checks[f'Hv_AD_{k}']=err<1e-11;checks[f'Hv_FD_{k}']=fd_err<2e-5;hvs.append({'AD_relative':err,'FD_relative':fd_err})
    for gain in [0.,1.,5.]:
        e,g=m.objective(vector,gain);ee,gg=parent.objective(vector,gain);checks[f'energy_gradient_bitidentity_gain{gain}']=e==ee and np.array_equal(g,gg)
    solved=m._attempt(.2,np.zeros(len(m.free)));checks['sparse_policy_activated_force_geometry']=solved['converged']
    (out/'small_activated.json').write_text(json.dumps(solved,indent=2)+'\n');np.savez_compressed(out/'small_hessian_comparison.npz',dense=dense,rows=m.pattern_rows,cols=m.pattern_cols,sparse_data=H[m.pattern_rows,m.pattern_cols].A.ravel(),vector=vector,colors=m.colors)
    benchmarks=[]
    for kind in ['dense','sparse']:
        completed=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--benchmark',kind],capture_output=True,text=True,check=True)
        benchmarks.append(json.loads(completed.stdout))
    receipt={'status':'PASS_SPARSE_COMPONENT' if all(checks.values()) else 'FAIL','checks':{k:bool(v) for k,v in checks.items()},'entire_H_absolute_difference':delta,'entire_H_relative_difference':delta/scale,'offpattern_absolute':off,'raw_asymmetry':meta['raw_asymmetry_linf'],'Hv_checks':hvs,'pattern':m.pattern_receipt,'benchmarks':benchmarks,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'parents':source_hashes(),'scope':'Samephysicalenergy/gradient, differentverifiednumericalpolicy; nophysical/PDEerrorclaim.'}
    (ROOT/'evidence/sparse_quadratic_verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2),flush=True)
    if not all(checks.values()):raise RuntimeError('sparse component not verified; no24run permitted')


def run24():
    source=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();verified=json.loads((ROOT/'evidence/sparse_quadratic_verification.json').read_text())
    if verified['status']!='PASS_SPARSE_COMPONENT' or verified['source_sha256']!=source:raise RuntimeError('current sparse source must verify before24x12run')
    out=ROOT/'runs/sparse_quadratic/24x12_lambda5000';out.mkdir(parents=True,exist_ok=True);lock=out/'job.json'
    if lock.exists():raise RuntimeError('24x12job record already exists; do not duplicate it')
    p=Parameters(nr=24,nz=12,bulk=5000.,patch_radius=9.35512037/33.5,activation_z_min=8.5/33.5,activation_z_max=28.5/33.5,basal_area_modulus=2.,apical_area_modulus=.2,basal_surface_tension=2.,apical_surface_tension=0.)
    model=build_model(p);protocol={'scope':'ONEtargetfree24x12highQactualsmallgeometrypath; noresponsefit','parameters':asdict(p),'integration_policy':asdict(model.integration_policy),'solver_policy':asdict(model.policy),'sparse_policy':asdict(model.sparse_policy),'pattern':model.pattern_receipt,'targets':[1.,5.,10.],'source_sha256':source,'parents':source_hashes()}
    (ROOT/'evidence/sparse_quadratic_24_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n');lock.write_text(json.dumps({'status':'RUNNING','pid':os.getpid(),'source_sha256':source},indent=2)+'\n');rows=[];previous=None
    try:
        for gain in protocol['targets']:
            (result,memory)=measured(lambda:model.advance(gain,previous,output_directory=out/f'gain_{gain}'))
            row={'strength':gain,'reached':result['reached_target'],'force_linf':result['scaled_free_force_linf'],'shortening_um':result['maximum_activated_column_shortening']*33.5,'sampled_J_min':result['sampled_minimum_J'],'sampled_J_max':result['sampled_maximum_J'],'conservative_J_lower':result['minimum_J'],'conservative_J_upper':result['maximum_J'],'attempts':len(result['robust_continuation']['attempts']),'rejected_attempts':sum(not a['accepted'] for a in result['robust_continuation']['attempts']),'resource_measurement':memory}
            rows.append(row);(ROOT/'evidence/sparse_quadratic_24_results.json').write_text(json.dumps({'status':'COMPLETE' if len(rows)==3 and all(r['reached'] for r in rows) else 'RUNNING_OR_STOPPED','rows':rows,'protocol':protocol},indent=2)+'\n');print(json.dumps(row),flush=True)
            if not result['reached_target']:break
            previous=result
        lock.write_text(json.dumps({'status':'FINISHED','all_targets_reached':len(rows)==3 and all(r['reached'] for r in rows),'pid':os.getpid(),'source_sha256':source},indent=2)+'\n')
    except BaseException as error:
        lock.write_text(json.dumps({'status':'STOPPED_WITH_ERROR','error':str(error),'pid':os.getpid(),'source_sha256':source},indent=2)+'\n');raise


if __name__=='__main__':
    torch.set_num_threads(1);parser=argparse.ArgumentParser();parser.add_argument('--benchmark',choices=['dense','sparse']);parser.add_argument('--run24',action='store_true');args=parser.parse_args()
    if args.benchmark:_benchmark(args.benchmark)
    elif args.run24:run24()
    else:self_test()
