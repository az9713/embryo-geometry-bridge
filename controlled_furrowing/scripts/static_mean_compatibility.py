"""Sharp kinematic consistency of STATIC regional-mean endpoint summaries.

No constitutive equilibrium, common material, temporal peak, noise model or
physical identified set is inferred. Fractions refer to the recorded decimals.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import itertools
import json
from joint_calibration_contract import load_observations,OBSERVABLES,matching_certificate

ROOT=Path(__file__).resolve().parents[1]


def correction(row):
    M,V,A,B=(row[k] for k in OBSERVABLES)
    if V>M:
        raise ValueError('this exact construction requires observed mean <= maximum')
    t=(A+B-V)/3
    aa,bb,vv=A-t,B-t,V+t
    mm=max(M,vv)
    error=max(abs(mm-M),abs(vv-V),abs(aa-A),abs(bb-B))
    if vv != aa+bb or mm < vv or error != abs(t):
        raise RuntimeError('static mean matching construction failed')
    return {'M':mm,'V':vv,'A':aa,'B':bb,'error':error,'t':t}


def polynomial_moments(corrected,weighting):
    """Core r/R=t, A(t)=Aconst, S(t)=M-Ct^2, B=S-A."""
    M,V,A,B=(corrected[k] for k in ('M','V','A','B'))
    factor=Q(3) if weighting=='line' else Q(2) if weighting=='disk' else None
    if factor is None:raise ValueError('line or disk required')
    C=factor*(M-V)
    # Normalized moments of t^2 on[0,1] are1/3(line),1/2(disk).
    meanS=M-C/factor
    return {'S_coefficients':[M,Q(0),-C],'A_coefficients':[A],
            'B_coefficients':[M-A,Q(0),-C], 'mean_S':meanS,'mean_A':A,
            'mean_B':meanS-A,'max_S':M,'min_S':M-C}


def qserial(value):
    if isinstance(value,Q):return str(value)
    if isinstance(value,dict):return {k:qserial(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [qserial(v) for v in value]
    return value


def main():
    observed=load_observations();ids=list(observed);height=Q('33.5')
    rows=[];checks={};canonical=max(abs(row[OBSERVABLES[2]]+row[OBSERVABLES[3]]-row[OBSERVABLES[1]])/3
                                 for row in observed.values())
    for rid,row in observed.items():
        c=correction(row)
        if c['M']>=height:raise RuntimeError('height insufficient for the kinematic construction')
        variants={w:polynomial_moments(c,w) for w in ('line','disk')}
        checks[rid+'_all_four_residuals'] = c['error']<=canonical
        checks[rid+'_positive_thickness'] = height-c['M']>0
        checks[rid+'_both_mean_constructions'] = all(
            v['mean_S']==c['V'] and v['mean_A']==c['A'] and v['mean_B']==c['B']
            and v['max_S']==c['M'] for v in variants.values())
        rows.append({'record_id':rid,'original':row,'correction':c,'constructions':variants,
                     'core_J_lower_bound_exact':(height-c['M'])/height})
    sums=sorted((row[OBSERVABLES[2]]+row[OBSERVABLES[3]],rid) for rid,row in observed.items())
    means=sorted((row[OBSERVABLES[1]],rid) for rid,row in observed.items())
    permutation=max(abs(a[0]-b[0])/3 for a,b in zip(sums,means))
    pairing=[{'predicted_record_id':b[1],'component_pair_source_id':a[1],
              'required_discrepancy':abs(a[0]-b[0])/3} for a,b in zip(sums,means)]
    costs=[[abs(observed[j][OBSERVABLES[2]]+observed[j][OBSERVABLES[3]]-observed[i][OBSERVABLES[1]])/3
            for j in ids] for i in ids]
    thresholds=sorted({Q(0)}|{cost for row in costs for cost in row})
    index=thresholds.index(permutation)
    def graph(delta):return [[j for j,c in enumerate(row) if c<=delta] for row in costs]
    witness=matching_certificate(graph(permutation),len(ids))
    predecessor=thresholds[index-1] if index else None
    obstruction=matching_certificate(graph(predecessor),len(ids)) if predecessor is not None else None
    checks['unknown_pairing_matching_at_sorted_optimum']=witness['all_left_matched']
    checks['unknown_pairing_predecessor_Hall_obstruction']=predecessor is not None and not obstruction['all_left_matched']
    # Check exact sharpness/constructive moments under the selected pairing too.
    for pair in pairing:
        fixed=observed[pair['predicted_record_id']];comp=observed[pair['component_pair_source_id']]
        row={**fixed,OBSERVABLES[2]:comp[OBSERVABLES[2]],OBSERVABLES[3]:comp[OBSERVABLES[3]]}
        c=correction(row)
        checks[pair['predicted_record_id']+'_unknown_pairing_primal']=c['error']<=permutation and c['M']<height
        pair['correction']=c
    # An arbitrary3point exact fixture independently enumerates all bijections.
    fixture_s=[Q(7),Q(1),Q(5)];fixture_v=[Q(0),Q(3),Q(6)]
    brute=min(max(abs(fixture_s[j]-fixture_v[i])/3 for i,j in enumerate(p))
              for p in itertools.permutations(range(3)))
    sorted_cost=max(abs(a-b)/3 for a,b in zip(sorted(fixture_s),sorted(fixture_v)))
    checks['independent_three_pair_permutation_enumeration']=brute==sorted_cost
    result={'status':'PASS_SHARP_STATIC_KINEMATIC_COMPATIBILITY' if all(checks.values()) else 'FAIL',
            'checks':checks,'source_csv_sha256':hashlib.sha256((ROOT/'data/derived/calibration_small_18um.csv').read_bytes()).hexdigest(),
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'static_regional_mean_canonical_floor_exact':canonical,
            'static_regional_mean_canonical_floor_display':float(canonical),
            'canonical_rows':rows,'witness_records':[r['record_id'] for r in rows if r['correction']['error']==canonical],
            'unknown_intact_pair_permutation_floor_exact':permutation,
            'unknown_intact_pair_permutation_floor_display':float(permutation),
            'unknown_pairing_construction':pairing,'matching_at_optimum':witness,
            'previous_critical_threshold':predecessor,'previous_threshold_Hall_obstruction':obstruction,
            'conditional_positive_height_um':height,
            'scope':['STATIC shared-ROI/shared-reference regional means only',
                'sharp over independent kinematically admissible fields, not common mechanical equilibria',
                'no statement that source component means were measured this way',
                'not valid for independently maximized temporal component means',
                'not a measurement-noise estimate or a full material identified set',
                'standard additive linear constraint and one-dimensional bottleneck matching, no novelty claim']}
    target=ROOT/'evidence/static_mean_compatibility.json'
    target.write_text(json.dumps(qserial(result),indent=2)+'\n')
    print(json.dumps({'status':result['status'],'checks':len(checks),
                      'canonical_floor_um':float(canonical),'unknown_pairing_floor_um':float(permutation),
                      'failures':[k for k,v in checks.items() if not v]}),flush=True)
    if not all(checks.values()):raise RuntimeError('static compatibility checks failed')


if __name__=='__main__':main()
