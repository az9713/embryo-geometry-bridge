"""Exact all-four-observable constraints for the 17 small-region records.

This checks bounds supplied by a forward model. It does NOT establish their
physical validity, observation semantics, PDE error, or a full identified set.
All four discrepancy allowances must be explicit; no noise law is invented.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import csv
import hashlib
import itertools
import json
from numbers import Integral

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/derived/calibration_small_18um.csv'
OBSERVABLES = ('maximum_shortening_um', 'region_averaged_shortening_um',
               'apical_max_shortening_um', 'basal_max_shortening_um')


def rational(value):
    if isinstance(value, bool):
        raise ValueError('Boolean is not a physical value')
    try:
        return value if isinstance(value, Q) else Q(str(value))
    except (ValueError, ZeroDivisionError, TypeError) as error:
        raise ValueError('finite rational numeric value required') from error


def load_observations():
    with DATA.open(newline='', encoding='utf-8-sig') as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 17 or len({r['record_id'] for r in rows}) != 17:
        raise ValueError('expected exactly 17 distinct source-column records')
    result = {r['record_id']:{k:rational(r[k]) for k in OBSERVABLES} for r in rows}
    if any(v < 0 for row in result.values() for v in row.values()):
        raise ValueError('this source cohort has nonnegative recorded responses')
    return result


def interval_residual(observed, interval):
    """Sharp min/max absolute residual on one closed prediction interval."""
    if not isinstance(interval,(tuple,list)) or len(interval) != 2:
        raise ValueError('two ordered closed interval endpoints required')
    y = rational(observed)
    lo, hi = map(rational, interval)
    if lo > hi:
        raise ValueError('prediction interval endpoints out of order')
    return max(lo-y, y-hi, Q(0)), max(abs(lo-y), abs(hi-y))


def evaluate(observations, predictions, discrepancies):
    """Conservative IN/OUT/UNRESOLVED predicate for a model's enclosed output.

    IN means every value in every supplied interval satisfies its constraint.
    OUT means some interval cannot intersect the corresponding allowed band.
    Otherwise the bounds do not decide membership. In particular interval
    overlap never proves a physically realizable compatible parameter exists.
    """
    if not observations or set(predictions) != set(observations):
        raise ValueError('exact same nonempty record IDs required; no skipped embryos')
    if set(discrepancies) != set(OBSERVABLES):
        raise ValueError('explicit discrepancy allowance for each of four observables required')
    delta = {k:rational(discrepancies[k]) for k in OBSERVABLES}
    if min(delta.values()) < 0:
        raise ValueError('nonnegative discrepancy allowances required')
    rows = []
    lower = {k:Q(0) for k in OBSERVABLES}
    upper = {k:Q(0) for k in OBSERVABLES}
    for rid, observed in observations.items():
        if set(observed) != set(OBSERVABLES) or set(predictions[rid]) != set(OBSERVABLES):
            raise ValueError('every record must supply all four observables')
        for key in OBSERVABLES:
            infimum, supremum = interval_residual(observed[key], predictions[rid][key])
            lower[key] = max(lower[key], infimum)
            upper[key] = max(upper[key], supremum)
            status = 'OUT' if infimum > delta[key] else ('IN' if supremum <= delta[key] else 'UNRESOLVED')
            rows.append({'record_id':rid, 'observable':key, 'observed_exact':str(observed[key]),
                         'prediction_interval_exact':list(map(lambda x:str(rational(x)), predictions[rid][key])),
                         'minimum_interval_residual_exact':str(infimum),
                         'maximum_interval_residual_exact':str(supremum),
                         'discrepancy_allowance_exact':str(delta[key]), 'status':status})
    status = 'OUT' if any(r['status'] == 'OUT' for r in rows) else (
        'IN' if all(r['status'] == 'IN' for r in rows) else 'UNRESOLVED')
    return {'status':status, 'records':len(observations), 'constraints':len(rows),
            'discrepancy_lower_requirements_exact':{k:str(v) for k,v in lower.items()},
            'discrepancy_sufficient_requirements_exact':{k:str(v) for k,v in upper.items()},
            'constraint_rows':rows,
            'scope':'Membership conditional on valid forward intervals and matching observation contract; no automatic physical/continuum/identified-set certificate.'}


def matching_certificate(adjacency, right_count):
    """Maximum matching and, if deficient, a directly checkable Hall witness."""
    if isinstance(right_count,bool) or not isinstance(right_count,Integral) or right_count < 0:
        raise ValueError('nonnegative non-Boolean integer right vertex count required')
    if not isinstance(adjacency,(tuple,list)) or any(not isinstance(row,(tuple,list,set,range)) for row in adjacency):
        raise ValueError('adjacency must contain collections of integer right vertices')
    if any(isinstance(j,bool) or not isinstance(j,Integral) for row in adjacency for j in row):
        raise ValueError('non-Boolean integer right vertices required')
    graph = [tuple(sorted(set(map(int,row)))) for row in adjacency]
    if any(j < 0 or j >= right_count for row in graph for j in row):
        raise ValueError('right vertex outside graph')
    right_to_left = [-1]*right_count

    def augment(i, visited):
        for j in graph[i]:
            if j in visited:
                continue
            visited.add(j)
            if right_to_left[j] == -1 or augment(right_to_left[j], visited):
                right_to_left[j] = i
                return True
        return False

    for i in range(len(graph)):
        augment(i,set())
    left_to_right = [-1]*len(graph)
    for j,i in enumerate(right_to_left):
        if i >= 0:
            left_to_right[i] = j
    result = {'matched_count':sum(j >= 0 for j in left_to_right),
              'left_to_right':left_to_right, 'all_left_matched':all(j >= 0 for j in left_to_right)}
    if not result['all_left_matched']:
        left = {i for i,j in enumerate(left_to_right) if j == -1}
        right = set()
        queue = list(left)
        while queue:
            i = queue.pop()
            for j in graph[i]:
                if left_to_right[i] == j:
                    continue
                right.add(j)
                paired = right_to_left[j]
                if paired >= 0 and paired not in left:
                    left.add(paired); queue.append(paired)
        neighbors = {j for i in left for j in graph[i]}
        if neighbors != right or len(neighbors) >= len(left):
            raise RuntimeError('matching algorithm did not yield a valid Hall obstruction')
        result['Hall_witness'] = {'left_indices':sorted(left), 'neighbor_right_indices':sorted(neighbors),
                                  'deficit':len(left)-len(neighbors)}
    return result


def evaluate_unknown_component_pairing(observations, predictions, discrepancies):
    """Optional single permutation of intact observed (apical,basal) pairs.

    Dose/maximum/mean record links remain fixed. This is a different, explicitly
    weaker contract than canonical source-column pairing, not an assertion that
    the source columns were unordered or biological samples were dependent.
    """
    canonical = evaluate(observations,predictions,discrepancies)
    ids = list(observations)
    fixed = [row for row in canonical['constraint_rows'] if row['observable'] in OBSERVABLES[:2]]
    outer,inner = [],[]
    for rid in ids:
        allowed,guaranteed = [],[]
        for j,source in enumerate(ids):
            bounds = [interval_residual(observations[source][k],predictions[rid][k])
                      for k in OBSERVABLES[2:]]
            if all(pair[0] <= rational(discrepancies[k]) for k,pair in zip(OBSERVABLES[2:],bounds)):
                allowed.append(j)
            if all(pair[1] <= rational(discrepancies[k]) for k,pair in zip(OBSERVABLES[2:],bounds)):
                guaranteed.append(j)
        outer.append(allowed); inner.append(guaranteed)
    oc = matching_certificate(outer,len(ids)); ic = matching_certificate(inner,len(ids))
    if any(row['status'] == 'OUT' for row in fixed) or not oc['all_left_matched']:
        status = 'OUT'
    elif all(row['status'] == 'IN' for row in fixed) and ic['all_left_matched']:
        status = 'IN'
    else:
        status = 'UNRESOLVED'
    return {'status':status, 'record_order':ids,
            'pairing_contract':'One unknown bijection of intact observed apical/basal pairs; dose/max/mean associations fixed',
            'fixed_maximum_mean_constraints':fixed,
            'outer_graph':outer, 'inner_graph':inner,
            'outer_matching_certificate':oc, 'inner_matching_certificate':ic,
            'canonical_status':canonical['status'],
            'scope':'Optional pairing-nuisance sensitivity only; not verified specimen linkage or a full physical identified set.'}


def zero_activation_check():
    """Joint exact residuals of the common gain=0 flat-equilibrium subfamily.

    All passive coefficients may vary within the positive/nonnegative physical
    domain, with flat fixed-rim geometry and the declared zero reference. This
    is one real, uninformative subfamily, not an estimate of biological noise.
    """
    observed = load_observations()
    predictions = {rid:{k:('0','0') for k in OBSERVABLES} for rid in observed}
    thresholds = {k:max(row[k] for row in observed.values()) for k in OBSERVABLES}
    evaluation = evaluate(observed, predictions, thresholds)
    facts = []
    for rid, row in observed.items():
        m, mean, apical, basal = (row[k] for k in OBSERVABLES)
        facts.append({'record_id':rid, 'maximum_minus_component_sum_exact':str(m-apical-basal),
                      'mean_minus_maximum_exact':str(mean-m)})
    result = {'status':'COMPLETE_EXACT_ZERO_ACTIVATION_BASELINE',
              'thresholds_um_exact':{k:str(v) for k,v in thresholds.items()},
              'thresholds_um_display':{k:float(v) for k,v in thresholds.items()},
              'joint_check':evaluation,
              'unknown_component_pairing_check':evaluate_unknown_component_pairing(observed,predictions,thresholds),
              'source_inequality_checks':facts,
              'source_inequalities_all_satisfied':all(
                  rational(r['maximum_minus_component_sum_exact']) <= 0 and
                  rational(r['mean_minus_maximum_exact']) <= 0 for r in facts),
              'source_sha256':hashlib.sha256(DATA.read_bytes()).hexdigest(),
              'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'assumptions':['flat reference and clamped rim', 'common gain exactly zero',
                             'zero reference component displacements', 'nonnegative surface tensions',
                             'positive shear and log-volume coefficients'],
              'limitations':['not a fitted model', 'not a measured discrepancy law',
                  'does not resolve source component aggregation/reference ambiguities',
                  'does not enumerate the full compatible mechanical set',
                  'observed inequalities alone do not identify the source observation algorithm']}
    (ROOT/'evidence/joint_zero_activation_baseline.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def self_test():
    checks = {}
    for i,(y,pair,expected) in enumerate([
            (3,(1,2),(1,2)), (3,(2,5),(0,2)), (3,(3,3),(0,0)),
            (-1,(-4,-2),(1,3)), (0,(-2,2),(0,2))]):
        checks[f'interval_example_{i}'] = interval_residual(y,pair) == tuple(map(Q,expected))
    o = {'x':{k:Q(1) for k in OBSERVABLES}}
    p = {'x':{k:('0','2') for k in OBSERVABLES}}
    checks['IN_inclusive_at_upper_threshold'] = evaluate(o,p,{k:'1' for k in OBSERVABLES})['status'] == 'IN'
    checks['overlap_is_not_membership'] = evaluate(o,p,{k:'0' for k in OBSERVABLES})['status'] == 'UNRESOLVED'
    checks['OUT_requires_disjoint_interval'] = evaluate(o,{'x':{k:('3','4') for k in OBSERVABLES}},
                                                     {k:'1' for k in OBSERVABLES})['status'] == 'OUT'
    # Independent exhaustive enumeration on seven rational-grid intervals:
    # all exact min/max residuals are attained by endpoints or the observed y.
    samples = [(-2,(-3,0)), (0,(1,3)), (5,(2,5)), (1,(-1,4)),
               (3,(0,2)), (-4,(-8,-7)), (2,(2,2))]
    grids = [sorted(set((Q(lo),Q(hi),max(Q(lo),min(Q(hi),Q(y)))))) for y,(lo,hi) in samples]
    all_costs = [max(abs(z-Q(sample[0])) for z,sample in zip(values,samples))
                 for values in itertools.product(*grids)]
    theoretical = [interval_residual(y,pair) for y,pair in samples]
    checks['independent_product_grid_minimum_attained'] = min(all_costs) == max(x[0] for x in theoretical)
    checks['independent_product_grid_maximum_attained'] = max(all_costs) == max(x[1] for x in theoretical)
    # Exhaust every 3x3 bipartite graph. The independent brute-force side tries
    # assignments, including unmatched vertices, instead of augmenting paths.
    graph_checks = []
    for mask in range(1 << 9):
        graph = [[j for j in range(3) if mask & (1 << (3*i+j))] for i in range(3)]
        cert = matching_certificate(graph,3)
        maximum = 0
        for assignment in itertools.product(range(-1,3),repeat=3):
            used = [j for j in assignment if j >= 0]
            if len(set(used)) == len(used) and all(j < 0 or j in graph[i] for i,j in enumerate(assignment)):
                maximum = max(maximum,len(used))
        good = cert['matched_count'] == maximum
        if not cert['all_left_matched']:
            hall = cert['Hall_witness']; subset = hall['left_indices']
            neighbors = sorted({j for i in subset for j in graph[i]})
            good = good and neighbors == hall['neighbor_right_indices'] and len(neighbors) < len(subset)
        graph_checks.append(good)
    checks['all_512_bipartite_graphs_and_Hall_witnesses'] = all(graph_checks)
    pair_observed = {str(i):{k:Q(i if k == OBSERVABLES[2] else 2-i if k == OBSERVABLES[3] else 0)
                                  for k in OBSERVABLES} for i in range(3)}
    pair_predicted = {str(i):{k:(str(pair_observed[str(2-i)][k]),)*2 for k in OBSERVABLES} for i in range(3)}
    pair_delta = {k:'0' for k in OBSERVABLES}
    checks['canonical_rejected_unknown_pairing_certified'] = (
        evaluate(pair_observed,pair_predicted,pair_delta)['status'] == 'OUT' and
        evaluate_unknown_component_pairing(pair_observed,pair_predicted,pair_delta)['status'] == 'IN')
    checks['unknown_pairing_keeps_components_together'] = (
        evaluate_unknown_component_pairing(pair_observed,pair_predicted,pair_delta)
        ['inner_matching_certificate']['left_to_right'] == [2,1,0])
    invalid = [lambda:interval_residual(0,(2,1)), lambda:interval_residual('nan',(0,1)),
               lambda:evaluate(o,{}, {k:'1' for k in OBSERVABLES}),
               lambda:evaluate(o,p,{}), lambda:evaluate(o,p,{k:'-1' for k in OBSERVABLES}),
               lambda:evaluate(o,{'x':{}},{k:'1' for k in OBSERVABLES}),
               lambda:matching_certificate([[0]],True),
               lambda:matching_certificate([[True]],2),
               lambda:matching_certificate([[1.0]],2),
               lambda:matching_certificate([[0.5]],2),
               lambda:interval_residual(0,'12'),
               lambda:interval_residual(0,{'1':0,'2':0})]
    for i,fn in enumerate(invalid):
        try:
            fn(); checks[f'invalid_input_{i}'] = False
        except ValueError:
            checks[f'invalid_input_{i}'] = True
    actual = zero_activation_check()
    checks['all_68_source_constraints_present'] = actual['joint_check']['constraints'] == 68
    checks['zero_family_membership_at_explicit_thresholds'] = actual['joint_check']['status'] == 'IN'
    checks['zero_family_also_passes_unknown_pairing'] = actual['unknown_component_pairing_check']['status'] == 'IN'
    checks['actual_source_inequalities'] = actual['source_inequalities_all_satisfied']
    checks['one_component_tightening_rejects'] = all(
        evaluate(load_observations(), {rid:{k:('0','0') for k in OBSERVABLES} for rid in load_observations()},
                 {k:rational(actual['thresholds_um_exact'][k])-(Q(1,10**8) if k == target else 0)
                  for k in OBSERVABLES})['status'] == 'OUT' for target in OBSERVABLES)
    receipt = {'status':'PASS' if all(checks.values()) else 'FAIL', 'checks':checks,
               'source_sha256':hashlib.sha256(DATA.read_bytes()).hexdigest(),
               'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'scope':'Exact joint interval-constraint logic and real 68-entry zero-activation baseline only, not completed mechanical inference.'}
    (ROOT/'evidence/joint_calibration_contract_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'], 'checks':len(checks),
                      'failures':[k for k,v in checks.items() if not v],
                      'zero_baseline_thresholds':actual['thresholds_um_display']}),flush=True)
    if not all(checks.values()):
        raise RuntimeError('joint constraint verification failed')


if __name__ == '__main__':
    self_test()
