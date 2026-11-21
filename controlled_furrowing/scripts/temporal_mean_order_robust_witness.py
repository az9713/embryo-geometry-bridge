"""Nonnegative-ROI temporal witnesses matching BOTH mean/max orderings.

Exact kinematics only, no shared mechanical law or force history. Earlier
parabolic witnesses remain unchanged and retain their narrower interpretation.
"""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import hashlib,json
from joint_calibration_contract import load_observations,OBSERVABLES

ROOT=Path(__file__).resolve().parents[1]


def profile(M,V,weight,H=Q('33.5'),radius=Q('9.35512037'),outer=Q('100.5')):
    if not 0<V<=M<H or not 0<radius<outer:
        raise ValueError('positive mean<=maximum<height and radii required')
    p=2
    def moment(n):
        if weight=='disk':return Q(1,n+1)
        if weight=='line':return sum((Q((-1)**k*comb(n,k),2*k+1) for k in range(n+1)),Q(0))
        raise ValueError('line or disk required')
    F=moment(p)
    while F>V/M:
        p+=1
        if p>10000:raise RuntimeError('computational polynomial degree budget exceeded')
        F=moment(p)
    C=(V-M*F)/(1-F)
    coefficients=[Q(0)]*(2*p+1);coefficients[0]=C
    for k in range(p+1):coefficients[2*k]+=(M-C)*(-1)**k*comb(p,k)
    actual_mean=sum(c*Q(1,j+1) for j,c in enumerate(coefficients)) if weight=='line' else sum(c*Q(2,j+2) for j,c in enumerate(coefficients))
    if actual_mean!=V or not 0<=C<=M:raise RuntimeError('nonnegative exact-mean profile failed')
    core_integral=radius**2/2*(C+(M-C)/Q(p+1))
    d=outer-radius;IL=d*(radius/2+d/6);IB=d*(radius/6+d/12)
    K=(core_integral+C*IL)/IB
    # Sann=C(1-t)-K*t(1-t)<=C<=M because K>=0. Thus height>=H-M.
    if K<0 or core_integral+C*IL-K*IB!=0:raise RuntimeError('exact volume compensation failed')
    return {'degree_parameter':p,'bump_mean_exact':str(F),'nonnegative_baseline_exact':str(C),
            'maximum_exact':str(M),'mean_exact':str(V),'core_polynomial_coefficients':list(map(str,coefficients)),
            'core_integral_rS_dr_exact':str(core_integral),'annulus_coefficients':list(map(str,(C,-C-K,K))),
            'total_integral_rS_dr_exact':'0','minimum_height_lower_exact':str(H-M),
            'factorization':'S=C+(M-C)*(1-(r/R)^2)^p >=0 onROI; J=1a.e. through compensated radial transport'}


def main():
    rows=[];checks={};boundary_rows=[]
    for rid,row in load_observations().items():
        M,V,A,B=(row[k] for k in OBSERVABLES)
        if min(M,V,A,B)<0 or V>A+B:raise ValueError('source component-peak compatibility required')
        first=(A,V-A);second=(V-B,B)
        checks[rid+'_exact_component_peaks']=max(Q(0),first[0],second[0])==A and max(Q(0),first[1],second[1])==B
        checks[rid+'_instantaneous_mean_identity']=sum(first)==sum(second)==V
        variants={}
        if V==0 and M>0:
            # Never remove the exceptional row or replace zero by an invented
            # measured epsilon. Exact local-peak averaging cannot attain this
            # vector under pointwise-zero initial reference and continuity.
            sequence={}
            for weight in ('line','disk'):
                terms=[]
                for p in (8,32,128,512):
                    F=Q(1,p+1)
                    if weight=='line':
                        F=Q(1)
                        for n in range(1,p+1):F*=Q(2*n,2*n+1)
                    epsilon=M*F
                    checks[rid+'_'+weight+'_approximation_'+str(p)]=0<epsilon<=A+B
                    terms.append({'p':p,'mean_error_exact':str(epsilon),'mean_error_display':float(epsilon),
                                  'core_profile':'M*(1-(r/R)^2)^p',
                                  'component_state1':[str(A),str(epsilon-A)],
                                  'component_state2':[str(epsilon-B),str(B)]})
                checks[rid+'_'+weight+'_decreasing_mean_sequence']=all(Q(b['mean_error_exact'])<Q(a['mean_error_exact']) for a,b in zip(terms,terms[1:]))
                sequence[weight]=terms
            boundary_rows.append({'record_id':rid,'maximum_exact':str(M),'reported_mean_exact':'0',
                'status':'EXACT_NONATTAINMENT_UNDER_INITIAL_ZERO_AND_LOCAL_PEAK_MEAN',
                'infimum_mean_error_exact':'0','approximating_sequences':sequence,
                'caveat':'Continuity, nonnegative normalized full-support weights and pointwise-zero initial reference are model assumptions, not verified source implementation.'})
            rows.append({'record_id':rid,'status':'NONATTAINMENT_BOUNDARY_RETAINED','observed_exact':{k:str(v) for k,v in zip(('M','V','A','B'),(M,V,A,B))}})
            continue
        for weight in ('line','disk'):
            result=profile(M,V,weight);C=Q(result['nonnegative_baseline_exact']);p=result['degree_parameter']
            # Independent moment formula by a product recurrence forline,
            # rather than the coefficient-sum formula used to choose p.
            recurrence=Q(1)
            for n in range(1,p+1):recurrence*=Q(2*n,2*n+1)
            independent_mean=C+(M-C)*(recurrence if weight=='line' else Q(1,p+1))
            checks[rid+'_'+weight+'_mean_exact']=independent_mean==V
            checks[rid+'_'+weight+'_ROI_profile_nonnegative']=0<=C<=M
            checks[rid+'_'+weight+'_positive_height']=Q(result['minimum_height_lower_exact'])>0
            checks[rid+'_'+weight+'_volume_exact']=Q(result['total_integral_rS_dr_exact'])==0
            variants[weight]=result
        rows.append({'record_id':rid,'status':'EXACT_BOTH_ORDERING_WITNESS','observed_exact':{k:str(v) for k,v in zip(('M','V','A','B'),(M,V,A,B))},
                     'component_state1':list(map(str,first)),'component_state2':list(map(str,second)),
                     'profiles':variants})
    checks['all_seventeen_records_retained']=len(rows)==17
    checks['actual_zero_mean_positive_max_boundary_is_O']=len(boundary_rows)==1 and boundary_rows[0]['record_id'].endswith('_O')
    result={'status':'COMPLETE_BOTH_ORDERING_CLASSIFICATION_WITH_NONATTAINMENT' if all(checks.values()) else 'FAIL',
            'checks':checks,'rows':rows,'nonattainment_boundary_rows':boundary_rows,
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'source_csv_sha256':hashlib.sha256((ROOT/'data/derived/calibration_small_18um.csv').read_bytes()).hexdigest(),
            'proof':['OnROI Sstar>=0, so ramp0->Sstar thenhold gives sup_t S(t,r)=Sstar(r) pointwise.',
                     'For V>0, both sup_t mean(S) andmean(sup_t S) therefore equalV.',
                     'Component means follow0->(A,V-A)->(V-B,B); their separatepeaks areA,B.',
                     'The maximum spatialheightshortening isM; annuluscompensation andradialtransport giveJ=1a.e.,positiveheight,fixedrim.',
                     'For V=0<M, continuity+pointwisezero initial reference imply mean(sup_tS)>0, so exact matching is impossible in this class.',
                     'M*(1-(r/R)^2)^p has fixed positive maximumM and mean->0 by dominated convergence, so the nonattained error infimum is zero.',
                     'Each record/weight usesitsownpath. No common constitutive/activation law or timescale is realized.'],
            'scope':'Sixteen exact both-ordering kinematic witnesses plus one retained nonattainment boundary; no fully source-aligned dynamics, measured error or common material fit.'}
    target=ROOT/'evidence/temporal_mean_order_robust_witness.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'checks':len(checks),'rows':len(rows),
                      'exact_records':sum(r['status']=='EXACT_BOTH_ORDERING_WITNESS' for r in rows),
                      'boundary_records':len(boundary_rows),
                      'max_degree_parameter':max(v['degree_parameter'] for r in rows if 'profiles' in r for v in r['profiles'].values())}))
    if not all(checks.values()):raise RuntimeError('order-robust witness checks failed')


if __name__=='__main__':main()
