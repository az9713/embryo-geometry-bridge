"""Implicit material sensitivities of saved high-Q equilibria, no optimization.

The implicit-function theorem is conditional on an exact regular root. Here
matrices are floating diagnostics at stored tolerance-verified states; neither
validated root existence nor global/experimental identifiability is claimed.
"""
from pathlib import Path
from dataclasses import replace
import argparse
import hashlib
import json
import numpy as np
import torch
from scipy.sparse.linalg import splu
from sparse_quadratic import AxisymmetricModel,Parameters

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/material_tangent_inference'
CASES={str(f):ROOT/'runs/activated_quadrature'/f'lambda_5000.0_strength_{float(f)}'/'reequilibrated.json'
       for f in (1,5,10)}
COLUMNS=('log_kappa','log_apical_area','log_basal_area',
         'apical_tension_additive','log_basal_tension','log_common_gain')
ATTRS=('bulk','apical_area_modulus','basal_area_modulus','apical_surface_tension','basal_surface_tension',None)
LOCKED={
 'sparse_quadratic.py':'d4833318151c7b94b4c56d40a96158d55a06258619ef080548d489e2be0394f9',
 'quadratic_quadrature_study.py':'0eac4bdcfd2bb0b725940d15874fc762c16065d2a29f2a10d0da70c5725ac5dc',
 'robust_quadratic.py':'be6f17e620ddf04533d456fff9d48cd48517f4a0709d921ab510616926b4176b',
 'axisymmetric_quadratic.py':'5bb64d6b441f6db8ada9da0f1bcb44716c86ed5d816e44f75bd4e9e021fdf410',
 'axisymmetric_mixed.py':'daf9b4f52181b18b1a75cc69909444bf2a85725c3896b32b125b2f20f10b1f14',
 'axisymmetric_prestress.py':'141b4b2ea92f4e04f4076b18edc9bafd6218f75a7c4c9633e66c5d2a9fc03c41',
 'axisymmetric_mechanics.py':'ca431e30a2574a02f4a21c98935c50217aa79bc5512cf7279143f22ce1a580e0'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,value):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def check_sources():
    actual={name:sha(Path(__file__).with_name(name)) for name in LOCKED}
    if actual!=LOCKED:raise RuntimeError('frozen physical sources changed')
    return actual


def mechanical_tangent(state):
    p=Parameters(**state['parameters']); model=AxisymmetricModel(p)
    u=np.asarray(state['free_displacement'],float);f=float(state['activation_strength'])
    if not state['converged'] or f<=0:raise ValueError('positive-strength verified state required')
    if any(getattr(p,k)<=0 for k in ('bulk','apical_area_modulus','basal_area_modulus','basal_surface_tension')):
        raise ValueError('log-coordinate coefficients must be strictly positive')
    energy,gradient=model.objective(u,f)
    H,meta=model.recover_hessian(u,f)
    if np.max(abs(gradient))>=1e-7:raise ValueError('saved state fails force replay')
    rhs=[];fd_rows=[];h=1e-3
    for label,attr in zip(COLUMNS,ATTRS):
        if attr is None:
            col=model.objective(u,2*f)[1]-gradient
            fd=(model.objective(u,f*np.exp(h))[1]-model.objective(u,f*np.exp(-h))[1])/(2*h)
        else:
            value=getattr(p,attr)
            additive=label=='apical_tension_additive'
            if additive:
                shifted=AxisymmetricModel(replace(p,**{attr:value+1.}))
                col=shifted.objective(u,f)[1]-gradient
                near=AxisymmetricModel(replace(p,**{attr:value+h}))
                fd=(near.objective(u,f)[1]-gradient)/h
            else:
                shifted=AxisymmetricModel(replace(p,**{attr:2*value}))
                col=shifted.objective(u,f)[1]-gradient
                plus=AxisymmetricModel(replace(p,**{attr:value*np.exp(h)}))
                minus=AxisymmetricModel(replace(p,**{attr:value*np.exp(-h)}))
                fd=(plus.objective(u,f)[1]-minus.objective(u,f)[1])/(2*h)
        relative=float(np.linalg.norm(col-fd)/max(np.linalg.norm(col),1e-14))
        fd_rows.append({'coordinate':label,'fixed_state_FD_relative_error':relative,
                        'FD_pass':relative<2e-4,'FD_step':h,
                        'FD_kind':'one-sided additive' if label=='apical_tension_additive' else 'central multiplicative'})
        rhs.append(col)
    B=np.stack(rhs,axis=1);factor=splu(H.tocsc());S=factor.solve(-B)
    residual=np.linalg.norm(H@S+B,axis=0)/np.maximum(np.linalg.norm(B,axis=0),1e-14)
    # Check the known common-force gauge, including the normalization's mu
    # derivative. All physical coefficients scale together, so the normalized
    # force is invariant even away from equilibrium (trial penalty inactive).
    uu=torch.tensor(u,dtype=torch.float64,requires_grad=True)
    parts=model.energy_components(model.coordinates(uu),f)
    shear_gradient=torch.autograd.grad(parts['shear']/model.energy_scale,uu)[0].detach().numpy()
    mu_rhs=shear_gradient-gradient
    gauge=B[:,0]+B[:,1]+B[:,2]+p.apical_surface_tension*B[:,3]+B[:,4]+B[:,5]+mu_rhs
    gauge_error=float(np.max(abs(gauge)))
    checks={'force_gate':float(np.max(abs(gradient)))<1e-7,
            'trial_penalty_inactive':float(parts['trial_inversion_penalty'])==0.,
            'all_fixed_state_material_FD_checks':all(r['FD_pass'] for r in fd_rows),
            'all_implicit_linear_systems_replayed':bool(np.max(residual)<1e-9),
            'common_force_gauge_cancels':gauge_error<1e-10}
    return model,u,H,B,S,{'checks':checks,'force_linf':float(np.max(abs(gradient))),
        'columns':COLUMNS,'fixed_state_derivative_checks':fd_rows,'linear_solve_relative_residuals':residual.tolist(),
        'gauge_rhs_linf':gauge_error,'Hessian_pattern':model.pattern_receipt,
        'interpretation':'Floating linearized response at a stored state. IFT derivative requires an exact regular equilibrium; no validated root or global identified set.'}


def mechanics():
    OUT.mkdir(parents=True,exist_ok=True)
    protocol_path=OUT/'mechanics_protocol.json'
    if protocol_path.exists():raise RuntimeError('existing tangent protocol; retain rather than overwrite or duplicate')
    protocol={'scope':'No-optimization material tangent audit on three saved12x6 high-Q fields',
              'columns':COLUMNS,'input_hashes':{k:sha(v) for k,v in CASES.items()},
              'physical_source_hashes':check_sources(),'script_sha256':sha(__file__),
              'FD_step':1e-3,'FD_relative_tolerance':2e-4,'linear_solve_relative_tolerance':1e-9,
              'caveats':['conditional implicit-function theorem, not a validated root',
                 'additive apical tension at zero is a one-sided admissible coordinate',
                 'finite compressibility kappa sensitivity is not a measured physiological modulus',
                 'no source observations, fitting, PDE certificate or full parameter-set inference']}
    save(protocol_path,protocol);rows=[]
    for key,path in CASES.items():
        state=json.loads(path.read_text());model,u,H,B,S,row=mechanical_tangent(state)
        folder=OUT/f'strength_{key}';folder.mkdir(exist_ok=True)
        np.savez_compressed(folder/'tangents.npz',free_vector=u,free_indices=model.free,
            H_data=H.data,H_indices=H.indices,H_indptr=H.indptr,H_shape=np.array(H.shape),
            parameter_force_derivatives=B,displacement_tangents=S)
        row.update({'strength':float(key),'input':path.relative_to(ROOT).as_posix(),'input_sha256':sha(path),
                    'matrix_archive_sha256':sha(folder/'tangents.npz')})
        save(folder/'verification.json',row);rows.append(row)
        print(json.dumps({'strength':key,'checks':row['checks'],
                          'max_FD_relative':max(r['fixed_state_FD_relative_error'] for r in row['fixed_state_derivative_checks']),
                          'max_solve_residual':max(row['linear_solve_relative_residuals'])}),flush=True)
        if not all(row['checks'].values()):raise RuntimeError('tangent checks failed; do not infer parameter combinations')
    save(ROOT/'evidence/material_tangent_mechanics_verification.json',{
        'status':'PASS_FLOATING_TANGENT_AUDIT','rows':rows,'protocol_sha256':sha(protocol_path),
        'sources_unchanged':check_sources()==protocol['physical_source_hashes'],'scope':protocol['scope']})


def means():
    from quadratic_mean_derivatives import mean_state_gradients
    source=Path(__file__).with_name('quadratic_mean_derivatives.py')
    protocol=json.loads((OUT/'mechanics_protocol.json').read_text())
    if protocol['script_sha256']!=sha(__file__) or check_sources()!=protocol['physical_source_hashes']:
        raise RuntimeError('tangent mechanical source changed')
    rows=[];jacobians={w:[] for w in ('current_disk','current_radial_line')}
    # The smooth apical/basal means are independent channels. The height mean
    # is their sum and must not be counted as a third independent observable.
    for key,path in CASES.items():
        state=json.loads(path.read_text());saved=np.load(OUT/f'strength_{key}'/'tangents.npz')
        free=saved['free_indices'];S=saved['displacement_tangents']
        for weighting in jacobians:
            result=mean_state_gradients(state,weighting=weighting)
            channel_rows=[]
            for name in ('apical_regional_mean_downward','basal_regional_mean_upward'):
                intervals=result['node_gradients'][name]
                gradient=np.array([[float((v.lo+v.hi)/2)*33.5 for v in node] for node in intervals]).ravel()[free]
                channel_rows.append(gradient@S)
            J=np.stack(channel_rows);jacobians[weighting].extend(channel_rows)
            rows.append({'strength':float(key),'weighting':weighting,'mean_material_jacobian':J.tolist(),
                         'gradient_scope':'Midpoints of exact physical-node derivative intervals, scaled by fixedL0=33.5 and selected freeDOFs; floating implicit tangents.'})
    summaries={}
    for weighting,data in jacobians.items():
        J=np.asarray(data);U,s,V=np.linalg.svd(J,full_matrices=False)
        summaries[weighting]={'matrix':J.tolist(),'singular_values':s.tolist(),
            'right_singular_vectors':V.tolist(),'columns':COLUMNS,
            'rank_warning':'Coordinates, chosen strengths and observation weights define this diagnostic; no noise scale, uncertainty interval or global identifiability theorem.',
            'feasible_direction_warning':'Apical tension zero permits only its nonnegative side physically; unconstrained SVD ignores that tangent-cone restriction.'}
    save(ROOT/'evidence/material_tangent_mean_diagnostic.json',{
        'status':'COMPLETE_FLOATING_LOCAL_MEAN_DIAGNOSTIC','observer_derivative_source_sha256':sha(source),
        'rows':rows,'summaries':summaries,
        'scope':'Three fixed diagnostic strengths, smoothmeans only. No exact maximum derivatives, data fit, full cohort rank, physical posterior or identified set.'})
    print(json.dumps({'status':'COMPLETE_FLOATING_LOCAL_MEAN_DIAGNOSTIC',
                      'singular_values':{k:v['singular_values'] for k,v in summaries.items()}}),flush=True)


if __name__=='__main__':
    torch.set_num_threads(1)
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['mechanics','means']);args=parser.parse_args()
    {'mechanics':mechanics,'means':means}[args.stage]()
