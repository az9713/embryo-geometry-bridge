"""Relink the unchanged candidate to the NEW independently reviewed point D."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib, json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = ROOT / 'runs/exact_integration_bridge_v1'
def read(p): return json.loads(p.read_text(encoding='utf8'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p, d): p.write_text(json.dumps(d, indent=2) + '\n', encoding='utf8')

pointfile = BASE / 'complete_point_B_safe_v2/results.json'
reviewfile = ROOT / 'evidence/exact_complete_point_B_independent_review_safe_v2.json'
assert sha(pointfile) == 'bc4984893680ee945b09797847a8926ece5ab482ff75c932bef73a6a20fbe8e0'
assert sha(reviewfile) == '433125c73e1ee6856f425ab7dc1d62658d9505a7c7ff8968c5ad2d4f7363dfc2'
point, review = read(pointfile), read(reviewfile)
assert review['status'] == 'PASS_SAFE_V2_COMPLETE_POINT_REDUCTION'
assert review['checks'] and all(v is True for v in review['checks'].values())
ledgerfile = BASE / 'retube_B_preparation/rounding_ledger.json'
vectorfile = BASE / 'retube_B_preparation/candidate_vectors.npz'
geomfile = BASE / 'retube_B_geometry/protocol.json'
geometry_review = ROOT / 'evidence/exact_retube_B_geometry_independent_review.json'
assert sha(geometry_review) == '8114999e0dd022fb00cb8aa95040e8bf6096c636fea4978e4bb7cc40c6f40002'
geometry = read(geomfile)
for p, h in geometry['source_and_input_hashes'].items():
    assert sha(ROOT / p) == h, p
ledger = read(ledgerfile)
with np.load(vectorfile, allow_pickle=False) as z:
    u0, u1 = z['u0'].copy(), z['u1'].copy()
d = [Q.from_float(float(a)) - Q.from_float(float(b)) for a, b in zip(u0, u1)]
checks = {}
def ck(k, v):
    checks[k] = bool(v)
    assert v, k
ck('candidate_actual_binary_decrement_all611_matches_ledger', len(d) == 611 and list(map(str, d)) == ledger['applied_decrement_exact'])
D = point['preconditioned_force_defect']
ck('safe_point_exact611_vector', len(D) == len(d))
residuals = []
for interval, value in zip(D, d):
    lo, hi = Q(interval['lower_exact']) - value, Q(interval['upper_exact']) - value
    ck('ordered_interval_' + str(len(residuals)), lo <= hi)
    residuals.append({'lower_exact': str(lo), 'upper_exact': str(hi)})
bound = max(max(abs(Q(v['lower_exact'])), abs(Q(v['upper_exact']))) for v in residuals)
dnorm = max(map(abs, d))
radius = Q(geometry['new_common_enclosing_radius_exact'])
ck('unchanged_common_radius_exact', radius == Q(geometry['old_radius_exact']) + dnorm)
ck('same_C_between_geometry_and_safe_point', geometry['parameter_cell'] == read(BASE / 'retube_B_preparation/protocol.json')['selected_record']['parameter_cell'])
ck('positive_safe_mismatch_below_one_millionth_old_radius', 0 < bound < Q(geometry['old_radius_exact']) / 10**6)
inputs = [pointfile, reviewfile, ledgerfile, vectorfile, geomfile, geometry_review, Path(__file__)]
result = {'status': 'PASS_SAFE_POINT_TO_UNCHANGED_CANDIDATE_RELINK', 'checks': checks,
          'input_hashes': {p.relative_to(ROOT).as_posix(): sha(p) for p in inputs},
          'safe_D_minus_actual_applied_d_norm_upper_exact': str(bound),
          'safe_D_minus_actual_applied_d_norm_upper_display': float(bound),
          'actual_applied_d_norm_exact': str(dnorm), 'common_radius_exact': str(radius),
          'coordinate_residuals': residuals,
          'old_unsafe_point_approval_used_as_premise': False,
          'new_integral_root_certified': False,
          'remaining_lemma_terms': ['fresh uniform cached-Q contraction ZQ', 'uniform full bulk+surface derivative discrepancy b', 'parameter/predictor variation V'],
          'scope': 'Exact rational subtraction of NEW SAFE independently reviewed D enclosures and unchanged actual binary candidate decrement. Existing candidate/common geometry unchanged. No new integration, residual evaluation, root transfer or PDE claim.'}
write(HERE / 'results.json', result)
print(json.dumps({'status': result['status'], 'checks': len(checks), 'safe_mismatch_upper': float(bound), 'results_sha256': sha(HERE / 'results.json')}))
