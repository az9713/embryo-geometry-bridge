"""Compose an audited energy-gain tube with exact whole-body geometry.

The full metric correction radius is converted to outward rational coordinate
radii via row norms of the fixed W. No new roots or fitted values are computed.
"""
from pathlib import Path
from fractions import Fraction as Q
import base64,hashlib,io,json
import numpy as np
from review_energy_gain_cell import sqrt_upper
from root_body_geometry import certify_box

ROOT=Path(__file__).resolve().parents[1]


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))


def main():
    out=ROOT/'evidence/energy_gain_geometry.json'
    if out.exists():raise RuntimeError('existing composition receipt; preserve instead of overwrite')
    ep=ROOT/'evidence/energy_gain_cell.json';er=ROOT/'evidence/energy_gain_cell_independent_review.json'
    bp=ROOT/'evidence/root_body_geometry_independent_review.json';energy=read(ep);review=read(er);body=read(bp)
    if review['status']!='PASS' or not review['checks'] or not all(v is True for v in review['checks'].values()):raise RuntimeError('independent energy review required')
    if review['producer_certificate_sha256']!=sha(ep) or review['producer_source_sha256']!=energy['script_sha256'] or sha(ROOT/'scripts/energy_gain_cell.py')!=energy['script_sha256']:raise RuntimeError('energy review/source binding differs')
    if not body['status'].startswith('PASS') or not body.get('checks') or not all(body['checks'].values()):raise RuntimeError('independent interval geometry review required')
    for path,value in energy['input_hashes'].items():
        if sha(ROOT/path)!=value:raise RuntimeError('energy input changed')
    point=read(ROOT/'evidence/algebraic_root_certification_feasibility.json');norm=read(ROOT/'evidence/parametric_norm_feasibility.json')
    raw=base64.b64decode(norm['validated_nodal_energy']['W_cache']['data'])
    if hashlib.sha256(raw).hexdigest()!=energy['metric_cache_sha256']:raise RuntimeError('metric digest differs')
    with np.load(io.BytesIO(raw),allow_pickle=False) as z:W=z['W'].copy()
    raw=base64.b64decode(point['frozen_binary_cache']['data'])
    if hashlib.sha256(raw).hexdigest()!=energy['model_cache_sha256']:raise RuntimeError('model digest differs')
    with np.load(io.BytesIO(raw),allow_pickle=False) as z:reference=z['reference'].copy();free=z['free'].copy();u=z['u'].copy()
    with np.load(ROOT/point['matrix_archive'],allow_pickle=False) as z:S=z['displacement_tangents'][:,5].copy()
    h=Q(energy['gain_relative_halfwidth_exact']);rho=Q(energy['metric_error_radius_exact'])
    errors=[]
    for i,row in enumerate(W):
        norm_squared=sum((Q.from_float(float(v))**2 for v in row),Q(0))
        errors.append(h*abs(Q.from_float(float(S[i])))+rho*sqrt_upper(norm_squared))
    index={int(k):i for i,k in enumerate(free)};centers=[];radii=[]
    for node,p in enumerate(reference):
        centers.append([Q.from_float(float(v))+(Q.from_float(float(u[index[2*node+j]])) if 2*node+j in index else Q(0)) for j,v in enumerate(p)])
        radii.append([errors[index[2*node+j]] if 2*node+j in index else Q(0) for j in (0,1)])
    state=read(ROOT/point['input_state']);geometry=certify_box(state['reference_nodes'],centers,state['triangles'],radii)
    result={'status':'CERTIFIED_ENERGY_GAIN_FAMILY_GEOMETRY' if geometry['global_injectivity_for_entire_box'] else 'UNRESOLVED_ENERGY_GAIN_GEOMETRY',
        'energy_certificate_sha256':sha(ep),'energy_review_sha256':sha(er),'body_review_sha256':sha(bp),
        'parameter_interval':'Relative active-gain increment [-1/1000,1/1000] at fixed old weak-load material',
        'full_metric_corrector_radius_exact':str(rho),'coordinate_radius_formula':'h*abs(S_i)+rho*sqrt(sum_j W_ij^2), exact binary inputs and outward rational sqrt',
        'free_coordinate_radii_exact':list(map(str,errors)),'max_free_coordinate_radius_exact':str(max(errors)),
        'geometry':geometry,'source_hashes':{name:sha(Path(__file__).with_name(name)) for name in ('compose_energy_geometry.py','root_body_geometry.py','review_energy_gain_cell.py')},
        'scope':'Every stationary root supplied by the audited one-parameter energy-gain theorem has an axis-preserving, positive-J, globally injective P2 body. Not actualM/all17/fullmaterialdomain/continuum/biology.'}
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'triangles':len(geometry['triangles']),'max_nodal_radius':float(max(errors)),
                     'J':{k:float(Q(v)) for k,v in geometry['J_bound'].items()} if geometry['J_bound'] else None}))


if __name__=='__main__':main()
