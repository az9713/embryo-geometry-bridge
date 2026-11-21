"""Read-only F/K bound decomposition; no new center/root/optimizer.

Exact Arb sufficient bounds are separated from floating archived-H diagnostics.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,base64,hashlib,io,json,os,time
import numpy as np
from scipy.linalg import cholesky,eigvalsh
from scipy.sparse import csr_matrix
from flint import arb,ctx
import parametric_root_cells as kernel

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/cohort_bound_diagnostics_v1'
OLD=ROOT/'runs/cohort_uniform_roots_v1'
SUP=ROOT/'runs/cohort_radius_supplement_v1'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
rel=lambda p:Path(p).resolve().relative_to(ROOT.resolve()).as_posix()


def save(p,v):
    if Path(p).exists():raise RuntimeError('retain existing diagnostic artifact')
    Path(p).parent.mkdir(parents=True,exist_ok=True);Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf8')


def freeze():
    p=read(OLD/'protocol.json');supp=read(SUP/'results.json');records=[];inputs={rel(OLD/'protocol.json'):sha(OLD/'protocol.json'),rel(SUP/'results.json'):sha(SUP/'results.json')}
    for row in p['records']:
        if row['record_id'] not in ('optoGEF2_small_col_F','optoGEF2_small_col_K'):continue
        path=ROOT/row['output_folder']/'result.json';prior=read(path);rec=next(r for r in supp['records'] if r['record_id']==row['record_id']);later=read(ROOT/rec['file'])
        if prior['status']!='UNRESOLVED' or later['status']!='UNRESOLVED':raise RuntimeError('retained unresolved inputs required')
        records.append({**row,'prior_result_file':rel(path),'prior_result_sha256':sha(path),'supplement_result_file':rec['file'],'supplement_result_sha256':rec['sha256'],
            'point_cache_file':prior['point_cache_file'],'predictor_file':prior['predictor_file']})
        inputs[rel(path)]=sha(path);inputs[rec['file']]=rec['sha256'];inputs.update(prior['artifact_hashes'])
    result={'source_sha256':sha(__file__),'source_hashes':{**p['source_hashes'],'cohort_bound_diagnostics.py':sha(__file__)},
        'input_hashes':inputs,'records':records,'shared_parameter_cell':p['shared_parameter_cell'],'shared_parameter_cell_sha256':p['shared_parameter_cell_sha256'],
        'base_parameters_exact':p['base_parameters_exact'],'precision_bits':160,'policy':'Reassemble exact cached point, predictor residual decomposition and predictor-only r=0 variation. Read historical radius tests. No fresh radius/root/center/optimizer.',
        'floating_policy':'Archived sparseH eigenvalue extrema/Cholesky, cachedR norm products; diagnostics only, not rigorous stability or condition bounds.',
        'runtime_versions':{'numpy':np.__version__},'no_nonlinear_optimization':True,'no_new_center':True}
    save(OUT/'protocol.json',result);print(json.dumps({'status':'PROTOCOL_FROZEN','source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json')}),flush=True)


def assert_inputs(p):
    if any(sha(ROOT/'scripts'/k)!=v for k,v in p['source_hashes'].items()) or any(sha(ROOT/k)!=v for k,v in p['input_hashes'].items()):raise RuntimeError('frozen diagnostic input/source changed')


def run():
    p=read(OUT/'protocol.json');assert_inputs(p)
    if (OUT/'results.json').exists():raise RuntimeError('retain diagnostic results')
    if os.environ.get('OPENBLAS_NUM_THREADS')!='1':raise RuntimeError('singlethread BLAS required')
    ctx.prec=p['precision_bits'];results=[];started=time.time()
    for row in p['records']:
        assert_inputs(p);kernel.OUT=(ROOT/row['predictor_file']).parent;kernel.A.cache_clear();prior=read(ROOT/row['prior_result_file']);later=read(ROOT/row['supplement_result_file'])
        local={'point_certificate_file':row['point_cache_file'],'predictor_file':Path(row['predictor_file']).name,'dose_exact':row['dose_exact'],'base_parameters_exact':p['base_parameters_exact']}
        point=kernel.build_point(local);C=kernel.cell_context(point,p['shared_parameter_cell'],'affine');n=point['n'];upper=kernel.upper;pack=kernel.pack
        g=[upper(point['gR'][i,0]) for i in range(n)]
        defect=[upper(sum((point['affine_M'][i,j]*C['delta'][j] for j in range(5)),arb(0))) for i in range(n)]
        linear=[upper(point['gR'][i,0]+sum((point['affine_M'][i,j]*C['delta'][j] for j in range(5)),arb(0))) for i in range(n)]
        rem,_=kernel.variation(point,C,arb(0),True);eta=max((linear[i]+rem[i]).upper() for i in range(n));imax=max(range(n),key=lambda i:float((linear[i]+rem[i]).upper()))
        variation0,J=kernel.variation(point,C,arb(0),False);Zpred=point['Z0']+max(variation0)
        pc=read(ROOT/row['point_cache_file'])
        with np.load(ROOT/pc['matrix_archive'],allow_pickle=False) as z:H=csr_matrix((z['H_data'],z['H_indices'],z['H_indptr']),shape=tuple(z['H_shape'])).toarray()
        raw=base64.b64decode(pc['frozen_binary_cache']['data'])
        with np.load(io.BytesIO(raw),allow_pickle=False) as z:R=z['preconditioner'].copy()
        smallest=float(eigvalsh(H,subset_by_index=[0,0])[0]);largest=float(eigvalsh(H,subset_by_index=[n-1,n-1])[0])
        try:cholesky(H,lower=True);chol=True
        except Exception:chol=False
        historical=[{'lineage':'v1','radius_exact':a['radius_exact'],'Z_uniform':a.get('Z_uniform'),'eta_plus_Zr':a.get('eta_plus_Zr'),'strict_selfmap_contraction':a.get('strict_selfmap_contraction')} for a in prior['uniform_cell_result']['radius_attempts']]
        historical +=[{'lineage':'radius_supplement','radius_exact':a['radius_exact'],'Z_uniform':a.get('Z_uniform'),'eta_plus_Zr':a.get('eta_plus_Zr'),'strict_selfmap_contraction':a.get('strict_selfmap_contraction')} for a in later['uniform_cell_result']['radius_attempts']]
        oldeta=Q(prior['uniform_cell_result']['eta_uniform']['upper_exact']);neweta=Q(str(eta.fmpq()))
        result={'record_id':row['record_id'],'free_dofs':n,'local_terms':len(point['cache']),'preconditioned_center_residual_linf':pack(max(g)),
            'affine_linear_defect_over_cell_linf':pack(max(defect)),'combined_center_plus_linear_linf':pack(max(linear)),
            'predictor_Taylor_remainder_linf':pack(max(rem)),'full_eta_uniform':pack(eta),
            'eta_worst_free_index':imax,'eta_worst_row_decomposition':{'center':pack(g[imax]),'linear_defect':pack(defect[imax]),'center_plus_linear':pack(linear[imax]),'remainder':pack(rem[imax])},
            'point_preconditioned_H_defect_Z0':pack(point['Z0']),'predictor_only_H_variation_linf':pack(max(variation0)),'predictor_only_total_Z_r0':pack(Zpred),
            'predictor_only_quadrature_J_lower':pack(J),'historical_tested_radius_bounds':historical,
            'eta_historical_comparison':{'old_exact_upper':str(oldeta),'fresh_exact_upper':str(neweta),'difference_exact':str(neweta-oldeta),'same_160bit_upper_exact':neweta==oldeta},
            'floating_archived_matrix_diagnostics':{'archived_H_smallest_eigenvalue':smallest,'archived_H_largest_eigenvalue':largest,'archived_H_Cholesky_succeeds':chol,
                'spectral_condition_diagnostic':largest/smallest if smallest>0 else None,'cached_R_linf_norm':float(np.abs(R).sum(axis=1).max()),'archived_H_linf_norm':float(np.abs(H).sum(axis=1).max()),
                'norm_product_Harchive_Rcache':float(np.abs(R).sum(axis=1).max()*np.abs(H).sum(axis=1).max()),
                'scope':'Floating archivedH at rounded solverstrength, not intervalH at exact gain*dose or exact unknownroot. No rigorous positive definiteness/stability/conditioning conclusion.'}}
        save(OUT/(row['record_id']+'.json'),result);results.append(result)
        print(json.dumps({'record':row['record_id'],'gR':float(max(g)),'linear_defect':float(max(defect)),'remainder':float(max(rem)),'eta':float(eta),'Z0':float(point['Z0']),'Zpred_r0':float(Zpred),'same_eta_160':neweta==oldeta}),flush=True)
    assert_inputs(p);summary={'status':'COMPLETE_BOUND_DIAGNOSTICS_NO_NEW_ROOT_CLAIM','protocol_sha256':sha(OUT/'protocol.json'),'source_sha256':sha(__file__),
        'records':results,'checks':{'both_F_K_retained':len(results)==2,'all611dof9504terms':all(r['free_dofs']==611 and r['local_terms']==9504 for r in results),
            'same_precision_eta_exactly_reproduced':all(r['eta_historical_comparison']['same_160bit_upper_exact'] for r in results),'all_sources_and_inputs_unchanged':True,
            'no_new_root_center_or_optimization':True},'elapsed_seconds':time.time()-started,'scope':'Decomposition of sufficient bounds, not evidence of root nonexistence or biological inadequacy; matrix spectra are floating diagnostics only.'}
    save(OUT/'results.json',summary);print(json.dumps({'status':summary['status'],'checks':summary['checks'],'elapsed':summary['elapsed_seconds'],'results_sha256':sha(OUT/'results.json')}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['freeze','run']);args=parser.parse_args();freeze() if args.stage=='freeze' else run()
