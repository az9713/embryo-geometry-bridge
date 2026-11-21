"""Apply independently reviewed fixed-quadrature root error to its endpoints.

This consumes a real certified radius, includes exact addition-rounding gaps,
and requires the independent root and transfer reviews. It certifies only the
algebraic-error contribution at the one archived field, never PDE accuracy.
"""
from pathlib import Path
from fractions import Fraction as Q
import base64,hashlib,io,json,math
import numpy as np
from quadratic_endpoint_observer import _state_controls
from quadratic_mean_derivatives import _boundary_node_ids
from nodal_observer_error_transfer import endpoint_bounds

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'evidence/algebraic_root_certification_feasibility.json'
REVIEW=ROOT/'evidence/algebraic_root_independent_review.json'
TRANSFER_REVIEW=ROOT/'evidence/nodal_observer_error_independent_review.json'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def outward_upper(value):
    exact=Q(value);number=float(exact)
    if exact<0 or not math.isfinite(number):raise ValueError('finite nonnegative display upper required')
    if Q.from_float(number)<exact:number=math.nextafter(number,math.inf)
    return number


def main():
    if not REVIEW.exists():raise RuntimeError('independent root review not yet available; no applied physical error claim')
    c=json.loads(CERT.read_text(encoding='utf-8'));review=json.loads(REVIEW.read_text(encoding='utf-8'));tr=json.loads(TRANSFER_REVIEW.read_text(encoding='utf-8'))
    if c['status']!='CERTIFIED_LOCAL_ROOT_BINARY_QUADRATURE':raise RuntimeError('actual root certificate required')
    if review.get('status')!='PASS' or not review.get('checks') or not all(v is True for v in review['checks'].values()):raise RuntimeError('independent root audit has not passed')
    if review.get('certificate_sha256')!=sha(CERT) or review.get('source_sha256')!=c['script_sha256']:raise RuntimeError('root review does not bind current certificate/source')
    if not tr['status'].startswith('PASS') or not tr['checks'] or not all(tr['checks'].values()):raise RuntimeError('independent conditional error-transfer review required')
    for path,value in tr['reviewed_hashes'].items():
        if sha(ROOT/path)!=value:raise RuntimeError('reviewed transfer artifact changed: '+path)
    if sha(ROOT/'scripts/algebraic_root_feasibility.py')!=c['script_sha256']:raise RuntimeError('root producer source changed')
    state_path=ROOT/c['input_state']
    if sha(state_path)!=c['input_state_sha256']:raise RuntimeError('root input state changed')
    state=json.loads(state_path.read_text(encoding='utf-8'));blob=base64.b64decode(c['frozen_binary_cache']['data'])
    if hashlib.sha256(blob).hexdigest()!=c['frozen_binary_cache']['sha256']:raise RuntimeError('frozen binary cache changed')
    with np.load(io.BytesIO(blob),allow_pickle=False) as z:
        reference=z['reference'];free=z['free'];u=z['u']
        if not np.array_equal(reference,np.asarray(state['reference_nodes'])) or not np.array_equal(u,np.asarray(state['free_displacement'])):
            raise RuntimeError('cache/reference/vector mismatch')
        selected={int(index):Q.from_float(float(v)) for index,v in zip(free,u)}
        exact=[Q.from_float(float(v))+selected.get(i,Q(0)) for i,v in enumerate(reference.ravel())]
    current=[Q.from_float(float(v)) for pair in state['current_nodes'] for v in pair]
    if len(current)!=len(exact):raise RuntimeError('nodal dimensions mismatch')
    gap=max(abs(a-b) for a,b in zip(current,exact))
    attempts=[r for r in c['radius_attempts'] if r.get('strict_selfmap_and_contraction')]
    if not attempts:raise RuntimeError('no successful actual root ball')
    producer_radius=min(Q(a['root_distance_upper']['upper_exact']) for a in attempts)
    independent_radius=Q(review['root_distance_upper_exact'])
    # Use the looser of two independently retained valid certificates, making
    # the downstream number directly supported by either proof calculation.
    radius=max(producer_radius,independent_radius)
    L=Q('33.5');epsilon=L*(radius+gap)
    ids=_boundary_node_ids(state);free_set=set(map(int,free));endpoint_dofs=[]
    for name in ('top','bottom'):
        endpoint_dofs.extend([2*ids[name][0][0],2*ids[name][-1][-1]])
    if any(index in free_set for index in endpoint_dofs):raise RuntimeError('fixed radial endpoint assumption not verified by actual free set')
    if any(exact[index]!=current[index] for index in endpoint_dofs):raise RuntimeError('fixed radial endpoints differ between exact center and observer')
    controls,H=_state_controls(state,L);roi=(Q(0),Q.from_float(float(state['parameters']['patch_radius']))*L)
    transferred=endpoint_bounds(controls['top'],controls['bottom'],epsilon,fixed_radial_endpoints=True,roi_um=roi)
    if transferred['status']!='CONDITIONAL_EXACT_ERROR_TRANSFER':raise RuntimeError('actual radius fails transfer prerequisites')
    result={'status':'APPLIED_REVIEWED_ALGEBRAIC_ERROR_TO_STATIC_ENDPOINTS','checks':{
        'independent_root_review_bound':True,'independent_transfer_review_bound':True,
        'exact_cache_and_state_match':True,'actual_fixed_radial_endpoints_verified':True,
        'positive_width_shared_ROI_and_strict_radial_margin':True,'addition_rounding_gap_included':True},
        'certificate_sha256':sha(CERT),'root_review_sha256':sha(REVIEW),'transfer_review_sha256':sha(TRANSFER_REVIEW),
        'input_state_sha256':sha(state_path),'script_sha256':sha(__file__),
        'root_distance_dimensionless_upper_exact':str(radius),
        'producer_root_distance_upper_exact':str(producer_radius),
        'independent_root_distance_upper_exact':str(independent_radius),
        'root_bound_selection':'Larger of producer and independent exact upper bounds',
        'serialized_node_addition_gap_exact':str(gap),
        'physical_nodal_coordinate_error_um_exact':str(epsilon),'transferred':transferred,
        'endpoint_algebraic_error_upper_um':{k:outward_upper(v) for k,v in transferred['endpoint_error_um_exact'].items()},
        'source_hashes':{name:sha(Path(__file__).with_name(name)) for name in ('nodal_observer_error_transfer.py',
            'quadratic_endpoint_observer.py','quadratic_mean_derivatives.py','annulus_observer.py')},
        'scope':'One independently reviewed finite-binary-quadrature root. Only algebraic root error relative to the stored curves; no quadrature/PDE, energetic stability, time/ROI/reference uncertainty, biological validity or full parameter-set certificate.'}
    path=ROOT/'evidence/algebraic_root_observer_error.json'
    if path.exists():raise RuntimeError('application receipt exists; preserve it rather than overwrite')
    path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'algebraic_endpoint_errors_um':result['endpoint_algebraic_error_upper_um']}))


if __name__=='__main__':main()
