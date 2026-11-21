"""Volume-preserving kinematic witnesses for the static-mean compatibility floor.

Exact polynomial volume/height checks. No equilibrium or common material law
is constructed. The continuum map is defined by monotone radial transport.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]


def maximum_quadratic(coefficients):
    c,b,a=coefficients
    values=[c,c+b+a]
    if a:
        t=-b/(2*a)
        if 0<t<1:values.append(c+b*t+a*t*t)
    return max(values)


def construct(corrected,weighting,H=Q('33.5'),A=Q('9.35512037'),Ro=Q('100.5')):
    M,V,aa,bb=(Q(corrected[k]) for k in ('M','V','A','B'))
    if not 0<A<Ro or H<=M or V>M or V!=aa+bb:
        raise ValueError('positive geometry and compatible corrected static means required')
    factor={'line':Q(3),'disk':Q(2)}[weighting]
    C=factor*(M-V);edge=M-C;d=Ro-A
    core_integral=A*A*(M/2-C/4)
    L_integral=d*(A/2+d/6)
    bump_integral=d*(A/6+d/12)
    K=(core_integral+edge*L_integral)/bump_integral
    annulus=[edge,-edge-K,K]  # t=(r-A)/(Ro-A), S=edge(1-t)-K*t(1-t).
    total=core_integral+edge*L_integral-K*bump_integral
    maximum=max(M,maximum_quadratic(annulus))
    if total!=0 or maximum>=H:
        raise RuntimeError('annulus did not preserve volume with positive height')
    return {'height_um':str(H),'ROI_radius_um':str(A),'outer_radius_um':str(Ro),
            'core_S_coefficients_normalized_radius':list(map(str,(M,Q(0),-C))),
            'core_apical_downward_constant':str(aa),
            'annular_S_coefficients_normalized_annulus_radius':list(map(str,annulus)),
            'annular_apical_downward_coefficients':list(map(str,(aa,-aa))),
            'core_integral_rS_dr_exact':str(core_integral),
            'total_integral_rS_dr_exact':str(total),
            'global_maximum_S_exact':str(maximum),
            'global_minimum_height_exact':str(H-maximum),
            'outer_S_exact':str(sum(annulus)),
            'total_target_volume_over_pi_exact':str(H*Ro*Ro),
            'J_identity':'rho(r)^2=(2/H)*integral_0^r s*(H-S(s))ds; R=rho_inverse; J=(R/rho)*R_prime*(H-S)/H=1',
            'map':'R=R(rho), Z=(S(R)-A_down(R))+(1-S(R)/H)*z',
            'scope':'Piecewise smooth globally embedded isochoric continuum configuration, NOT an equilibrium, P2 representation or common material/gain realization.'}


def main():
    source=ROOT/'evidence/static_mean_compatibility.json'
    data=json.loads(source.read_text());rows=[];checks={}
    inputs=[('canonical',r['record_id'],r['correction']) for r in data['canonical_rows']]
    inputs += [('unknown_pairing',r['predicted_record_id'],r['correction']) for r in data['unknown_pairing_construction']]
    for pairing,rid,correction in inputs:
        for weight in ('line','disk'):
            result=construct(correction,weight)
            name=pairing+'_'+rid+'_'+weight
            checks[name+'_volume_exact']=Q(result['total_integral_rS_dr_exact'])==0
            checks[name+'_height_positive']=Q(result['global_minimum_height_exact'])>0
            checks[name+'_rim_identity']=Q(result['outer_S_exact'])==0
            rows.append({'pairing':pairing,'record_id':rid,'weighting':weight,**result})
    result={'status':'PASS_EXACT_VOLUME_AND_HEIGHT_CONSTRUCTIONS' if all(checks.values()) else 'FAIL',
            'checks':checks,'constructions':rows,'count':len(rows),
            'source_compatibility_receipt_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Exact-volume and positive-height certificates plus analytic radial-transport J=1 proof. Separate fields for each row/weighting, no mechanical equilibrium or common inferred parameters.'}
    (ROOT/'evidence/static_isochoric_construction.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'checks':len(checks),'constructions':len(rows),
                     'minimum_height_um':min(float(Q(r['global_minimum_height_exact'])) for r in rows)}))
    if not all(checks.values()):raise RuntimeError('isochoric construction failed')


if __name__=='__main__':main()
