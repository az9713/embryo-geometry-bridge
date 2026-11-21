"""Independent proof audit plus isolated replay of the uniform energy cell.

The scalar acceptance arithmetic is independently reconstructed with exact
Fractions and an outward rational square root. Full Arb bound production is
replayed from its unchanged source into a separate output, never canonical
receipts. Replaying the producer is labelled as such, not independent AD.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import isqrt
import argparse,base64,hashlib,io,json,time,traceback
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs/energy_gain_cell_independent_review'
PRODUCER=ROOT/'evidence/energy_gain_cell.json'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def upper(v):return Q(v['upper_exact'])
def lower(v):return Q(v['lower_exact'])


def sqrt_upper(value,digits=80):
    value=Q(value)
    if value<0:raise ValueError('nonnegative radical required')
    a=isqrt(value.numerator);b=isqrt(value.denominator)
    if a*a==value.numerator and b*b==value.denominator:return Q(a,b)
    scale=10**digits;n=isqrt(value.numerator*scale*scale//value.denominator)
    if Q(n,scale)**2<value:n+=1
    return Q(n,scale)


def exact_review(c):
    checks={};rows=[];pieces=c['predictor_integral_pieces'];h=Q(c['gain_relative_halfwidth_exact']);rho=Q(c['metric_error_radius_exact'])
    plans=c['completed_plan_bounds']
    checks['valid_positive_gain_interval_and_metric_radius']=0<h<1 and rho>0
    checks['all_input_and_source_checks_true']=bool(c.get('checks')) and all(v is True for v in c['checks'].values())
    checks['all_9504_terms']=c['terms']==9504
    checks['correct_path_partition_count']=len(plans)==pieces+1 and pieces>0
    intervals=[]
    for i,p in enumerate(plans):
        if i<pieces:
            interval=tuple(map(Q,p['time_interval']));intervals.append(interval)
            checks[f'path_piece_{i}_exact_domain']=interval==(Q(i,pieces),Q(i+1,pieces)) and Q(p['shape_halfwidth_exact'])==h*interval[1] and Q(p['metric_error_radius_exact'])==0
        else:checks['whole_error_tube_exact_domain']=Q(p['shape_halfwidth_exact'])==h and Q(p['metric_error_radius_exact'])==rho
        a=upper(p['rank_center_squared_norm_upper']);d=upper(p['rank_change_squared_norm_upper'])
        rank=2*sqrt_upper(a*d)+d
        total=rank+upper(p['bulk_logQ_variation_bound'])+upper(p['surface_variation_bound'])+upper(p['gain_parameter_variation_bound'])
        checks[f'plan_{i}_rank_arithmetic_upper']=rank<=upper(p['volumetric_rank_variation_bound'])
        checks[f'plan_{i}_assembled_arithmetic_upper']=total<=upper(p['total_transformed_Hessian_variation_bound'])
        margins=p['analytic_domain_strict_lower_margins']
        checks[f'plan_{i}_all_analytic_margins_positive']=all(lower(v)>0 for v in margins.values()) and lower(margins['J'])>Q(1,10**6)
        rows.append({'name':p['name'],'independent_rank_upper_exact':str(rank),'independent_total_upper_exact':str(total)})
    integral=sum((upper(p['total_transformed_Hessian_variation_bound'])/pieces for p in plans[:-1]),Q(0))
    eta=upper(c['center_residual_energy_norm'])+h*upper(c['predictor_linear_defect_energy_norm'])+h*upper(c['predictor_metric_norm_upper'])*integral
    Z=upper(c['center_metric_defect'])+upper(plans[-1]['total_transformed_Hessian_variation_bound'])
    checks['path_integral_rational_upper']=integral<=upper(c['path_integral_H_variation_upper'])
    checks['eta_rational_upper']=eta<=upper(c['uniform_predictor_eta'])
    checks['Z_rational_upper']=Z<=upper(c['uniform_contraction_Z'])
    checks['strict_exact_contraction']=Z<1
    checks['strict_exact_self_mapping']=eta+Z*rho<rho
    # Stored author error uses unrounded Arb operands rather than the looser
    # separately serialized scalar uppers; do not demand a smaller bound here.
    error=eta/(1-Z)
    checks['exact_radius_error_below_declared_tube']=error<rho
    return checks,{'per_plan_arithmetic':rows,'eta_exact_upper':str(eta),'Z_exact_upper':str(Z),
                   'eta_plus_Zrho_exact_upper':str(eta+Z*rho),'root_metric_error_exact_upper':str(error)}


def run():
    if RUN.exists():raise RuntimeError('existing independent-review directory; inspect actual state, never duplicate')
    RUN.mkdir(parents=True);c=read(PRODUCER)
    if c['status']!='CERTIFIED_UNIFORM_ENERGY_GAIN_CELL':raise RuntimeError('successful original receipt required for this audit')
    checks,arithmetic=exact_review(c)
    for name,value in c['source_hashes'].items():checks['source_current_'+name]=sha(ROOT/'scripts'/name)==value
    for name,value in c['input_hashes'].items():checks['input_current_'+name]=sha(ROOT/name)==value
    checks['producer_source_current']=sha(ROOT/'scripts/energy_gain_cell.py')==c['script_sha256']
    norm=read(ROOT/'evidence/parametric_norm_feasibility.json');point=read(ROOT/'evidence/algebraic_root_certification_feasibility.json')
    cache=base64.b64decode(point['frozen_binary_cache']['data'])
    checks['binary_cache_hash_matches']=hashlib.sha256(cache).hexdigest()==c['model_cache_sha256']
    with np.load(io.BytesIO(cache),allow_pickle=False) as z:
        scale=Q.from_float(float(z['energy_scale'][0]));weights=list(z['volume_weights'].ravel())
        for k in (0,1):weights+=list(z[f'surface_{k}_weights'].ravel())
        checks['all_9504_cached_weights_positive']=len(weights)==9504 and scale>0 and all(Q.from_float(float(v))>0 for v in weights)
    metric=base64.b64decode(norm['validated_nodal_energy']['W_cache']['data'])
    checks['binary_metric_hash_matches']=hashlib.sha256(metric).hexdigest()==c['metric_cache_sha256']
    with np.load(io.BytesIO(metric),allow_pickle=False) as z:
        W=z['W'];checks['exact_binary_metric_triangular_nonsingular']=W.shape==(611,611) and np.array_equal(W,np.triu(W)) and bool(np.all(np.diag(W)!=0))
    # Exact independent identities used by the proof, checked at their
    # nontrivial boundary cases; general inequalities are justified in the doc.
    checks['sqrt_upper_exact_square']=sqrt_upper(Q(4,9))==Q(2,3)
    checks['sqrt_upper_outward_irrational']=sqrt_upper(Q(2))**2>=2
    checks['sqrt_upper_tight_irrational']=sqrt_upper(Q(2))**2-Q(2)<Q(1,10**70)
    initial={'status':'PASS_EXACT_REVIEW_PENDING_FULL_REPLAY' if all(checks.values()) else 'FAIL_EXACT_REVIEW',
             'checks':checks,'producer_sha256':sha(PRODUCER),'review_source_sha256':sha(__file__),
             'independent_arithmetic':arithmetic,'scope':'Scalar proof arithmetic and input/metric/domain checks; full bound reproduction pending.'}
    save(RUN/'preflight.json',initial)
    if not all(checks.values()):raise RuntimeError('independent exact preflight failed: '+str([k for k,v in checks.items() if not v]))
    frozen={'producer_sha256':sha(PRODUCER),'producer_source_sha256':sha(ROOT/'scripts/energy_gain_cell.py'),
            'source_hashes':c['source_hashes'],'input_hashes':c['input_hashes'],'reviewer_source_sha256':sha(__file__),
            'output_override_only':'Producer.OUT points to separate producer_replay.json; numerical function source and inputs unmodified.'}
    save(RUN/'protocol.json',frozen)
    # A deliberately labelled unchanged-source reproduction. The independent
    # certificate combination and mathematical audit above do not call producer.
    import energy_gain_cell as producer
    producer.OUT=RUN/'producer_replay.json'
    save(RUN/'job.json',{'status':'RUNNING_READONLY_REPLAY','pid':__import__('os').getpid(),'protocol_sha256':sha(RUN/'protocol.json')})
    started=time.monotonic()
    try:
        producer.attempt(Q(c['gain_relative_halfwidth_exact']),Q(c['metric_error_radius_exact']),c['predictor_integral_pieces'])
        replay=read(producer.OUT)
        selected=('center_metric_defect','center_residual_energy_norm','predictor_linear_defect_energy_norm','predictor_metric_norm_upper',
                  'path_integral_H_variation_upper','uniform_predictor_eta','uniform_contraction_Z','eta_plus_Zrho','selfmap_margin')
        checks['full_9504_term_replay_pass']=replay['status']=='CERTIFIED_UNIFORM_ENERGY_GAIN_CELL' and all(replay['checks'].values())
        checks['all_serialized_final_bounds_identical']=all(replay[k]==c[k] for k in selected)
        checks['all_plan_bounds_identical']=replay['completed_plan_bounds']==c['completed_plan_bounds']
        checks['original_receipt_and_source_unchanged']=sha(PRODUCER)==frozen['producer_sha256'] and sha(ROOT/'scripts/energy_gain_cell.py')==frozen['producer_source_sha256']
        checks['reviewer_source_unchanged']=sha(__file__)==frozen['reviewer_source_sha256']
        result={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,
            'producer_certificate_sha256':sha(PRODUCER),'producer_source_sha256':frozen['producer_source_sha256'],
            'reviewer_source_sha256':sha(__file__),'protocol_sha256':sha(RUN/'protocol.json'),
            'replay_receipt_sha256':sha(producer.OUT),'independent_arithmetic':arithmetic,
            'mathematical_review':['Relative gain changes only the linear active coefficient; F(u0;delta)=F0+Bdelta exactly.',
                'Predictor-path integral uses the full parameter variation on every t piece and only reduces the positional excursion.',
                'Rank split kappa vv^T: cross terms bounded by2||A0W||||DeltaA W|| plus||DeltaA W||^2.',
                'Symmetric row-diagonal majorants and positive weights imply Loewner bounds before exact metric congruence.',
                'Frozen W invertible; all spectral norm uses are for exactly symmetric targets.',
                'Uniform selfmap and contraction apply on fixed z ball with parameter-dependent affine predictor.',
                'Strict transformed Hessian positivity at the interior fixed point implies a local energy minimum, not global or dynamical stability.'],
            'known_guard_limitation':'Original producer status does not itself combine source/input booleans. This audit explicitly requires every acceptance and metadata check true, and verifies frozen bytes before/after. Actual run passed; a drifted future packet is not accepted.',
            'elapsed_seconds':time.monotonic()-started,
            'scope':'One fixed-material original weak-load gain interval only. Independent exact scalar audit plus unchanged-source full Arb reproduction, not independent reimplementation of all majorants. No actualM/all17/fullTheta, continuum, physical time or biological claim.'}
        save(ROOT/'evidence/energy_gain_cell_independent_review.json',result)
        save(RUN/'job.json',{'status':'FINISHED','pid':__import__('os').getpid(),'review_status':result['status'],
             'result_sha256':sha(ROOT/'evidence/energy_gain_cell_independent_review.json')})
        print(json.dumps({'status':result['status'],'checks':len(checks),'failed':[k for k,v in checks.items() if not v],'seconds':result['elapsed_seconds']}),flush=True)
        if result['status']!='PASS':raise RuntimeError('independent replay/binding checks failed')
    except BaseException as error:
        save(RUN/'job.json',{'status':'STOPPED_WITH_ERROR','pid':__import__('os').getpid(),'error':repr(error),'traceback':traceback.format_exc()});raise


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',action='store_true');args=parser.parse_args()
    if not args.run:raise SystemExit('Explicit --run required; no equilibrium optimization')
    run()
