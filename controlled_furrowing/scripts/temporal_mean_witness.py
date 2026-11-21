"""Exact two-peak kinematic paths matching the four recorded summaries.

This is a counterexample to equating separately maximized means with one
static state. No common material, force history, rate or biology is inferred.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json
from joint_calibration_contract import load_observations,OBSERVABLES
from static_isochoric_construction import construct

ROOT=Path(__file__).resolve().parents[1]


def main():
    rows=[];checks={}
    for rid,observed in load_observations().items():
        M,V,A,B=(observed[k] for k in OBSERVABLES)
        if min(M,V,A,B)<0 or V>M or V>A+B:
            raise ValueError('nonnegative summaries with V<=M and V<=A+B required')
        # State0 is the undeformed reference. Atstate1/2, regional component
        # maxima occur separately while total mean shortening stays V.
        first={'M':M,'V':V,'A':A,'B':V-A}
        second={'M':M,'V':V,'A':V-B,'B':B}
        peak={'M':max(Q(0),first['M'],second['M']),
              'V':max(Q(0),first['V'],second['V']),
              'A':max(Q(0),first['A'],second['A']),
              'B':max(Q(0),first['B'],second['B'])}
        checks[rid+'_all_four_temporal_peaks_exact']=tuple(peak[k] for k in ('M','V','A','B'))==(M,V,A,B)
        checks[rid+'_static_identity_each_state']=all(s['V']==s['A']+s['B'] for s in (first,second))
        constructions={}
        for weight in ('line','disk'):
            c1=construct(first,weight);c2=construct(second,weight)
            # S and therefore the compensated height/radial map are exactly
            # equal during the transition between these two component peaks.
            checks[rid+'_'+weight+'_constant_height_during_peak_transfer']=(
                c1['core_S_coefficients_normalized_radius']==c2['core_S_coefficients_normalized_radius'] and
                c1['annular_S_coefficients_normalized_annulus_radius']==c2['annular_S_coefficients_normalized_annulus_radius'])
            checks[rid+'_'+weight+'_volume_and_height_each_state']=all(
                Q(c['total_integral_rS_dr_exact'])==0 and Q(c['global_minimum_height_exact'])>0 for c in (c1,c2))
            constructions[weight]={'first':c1,'second':c2}
        rows.append({'record_id':rid,'recorded':{k:str(v) for k,v in zip(('M','V','A','B'),(M,V,A,B))},
                     'state1':{k:str(v) for k,v in first.items()},'state2':{k:str(v) for k,v in second.items()},
                     'temporal_peaks':{k:str(v) for k,v in peak.items()},'constructions':constructions})
    result={'status':'PASS_EXACT_TEMPORAL_KINEMATIC_WITNESSES' if all(checks.values()) else 'FAIL',
            'checks':checks,'records':rows,
            'source_hashes':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                             for name in ('temporal_mean_witness.py','static_isochoric_construction.py','joint_calibration_contract.py')},
            'source_csv_sha256':hashlib.sha256((ROOT/'data/derived/calibration_small_18um.csv').read_bytes()).hexdigest(),
            'continuous_time_construction':[
                '0<=t<=1: S(t,r)=t*S_star(r), a_mean=t*A, b_mean=t*(V-A); annulus compensation also scales linearly.',
                '1<=t<=2: S=S_star, a_mean=(2-t)*A+(t-1)*(V-B), b_mean=V-a_mean.',
                'Construct the isochoric radial-transport map at each time. Positiveheight follows from convex interpolation and M<H.',
                'Each mean is piecewise affine in t, so its supremum occurs at an endpoint of one time interval.'
            ],
            'scope':['separately maximized regional means and peak spatial height/mean summaries only',
                     'distinct kinematic path per record and weighting',
                     'volume-preserving embedded configurations, not force-balanced solutions',
                     'time is an arbitrary parameter, not seconds or a measured reaction/relaxation rate',
                     'no common constitutive law, gain, activation history or biological explanation',
                     'does not establish the original source used this particular observation algorithm',
                     'no noise estimate, identified set or novelty claim']}
    target=ROOT/'evidence/temporal_mean_witness.json';target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'checks':len(checks),'records':len(rows),
                      'failures':[k for k,v in checks.items() if not v]}))
    if not all(checks.values()):raise RuntimeError('temporal witness checks failed')


if __name__=='__main__':main()
