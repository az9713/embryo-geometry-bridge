"""Independent analytic geometry, gauge and numerical-integral checks.

No optimization, fitting or new mechanics state. Numerical quadrature is an
independent crosscheck, not the basis for the observer's exact enclosures.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib
import json
import math
import numpy as np
from scipy.integrate import quad
from quadratic_endpoint_observer import observe_graphs,observe_state

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = '075dac2148ed26a3bfe83fadf1b97602e4b08536f2be38cb95c2054776c0b926'


def graph(r,z):
    return [[(r(t),z(t)) for t in (Q(0),Q(1,2),Q(1))]]


def contains(v,x):
    return v.lo <= x <= v.hi


def run():
    source = Path(__file__).with_name('quadratic_endpoint_observer.py')
    checks = {'frozen_observer_source':hashlib.sha256(source.read_bytes()).hexdigest() == EXPECTED}
    a = graph(lambda t:t+t*t,lambda t:10-t)
    b = graph(lambda t:2*t,lambda t:1-t)
    # Material-paired S=1. At current R=t+t^2, the two graph values instead
    # give S=1+t/2-t^2/2, with maximum9/8, line mean13/12, disk mean43/40.
    for weight,mean in [('current_radial_line',Q(13,12)),('current_disk',Q(43,40))]:
        out = observe_graphs(a,b,roi_um=(0,2),reference_apical_um=10,
                             weighting=weight,tolerance_um='0.0001')
        end = out['four_endpoints']
        checks[weight+'_Eulerian_maximum_9_over_8'] = contains(end['maximum_Eulerian_height_shortening'],Q(9,8))
        checks[weight+'_mean_analytic'] = contains(end['regional_mean_Eulerian_height_shortening'],mean)
        checks[weight+'_separate_maxima_are_one_each'] = (contains(end['apical_downward'],Q(1)) and contains(end['basal_upward'],Q(1)))
        checks[weight+'_material_pair_proxy_is_different'] = end['maximum_Eulerian_height_shortening'].lo > 1
        shifted = observe_graphs(a,b,roi_um=(0,2),reference_apical_um=Q(21,2),reference_basal_um=Q(1,2),
                                 weighting=weight,tolerance_um='0.0001')['four_endpoints']
        checks[weight+'_reference_gauge_preserves_height'] = (
            shifted['maximum_Eulerian_height_shortening'] == end['maximum_Eulerian_height_shortening'] and
            shifted['regional_mean_Eulerian_height_shortening'] == end['regional_mean_Eulerian_height_shortening'])
        checks[weight+'_reference_gauge_moves_components_oppositely'] = (
            shifted['apical_downward'].lo-end['apical_downward'].lo == Q(1,2) and
            shifted['basal_upward'].hi-end['basal_upward'].hi == -Q(1,2))
        means = observe_graphs(a,b,roi_um=(0,2),reference_apical_um=10,weighting=weight,
                               component_aggregation='regional_mean',tolerance_um='0.0001')['four_endpoints']
        aa,bb = ((Q(7,12),Q(1,2)) if weight == 'current_radial_line' else (Q(89,120),Q(1,3)))
        checks[weight+'_component_mean_analytic'] = contains(means['apical_downward'],aa) and contains(means['basal_upward'],bb)
    path = ROOT/'runs/activated_quadrature/lambda_5000.0_strength_10.0/reequilibrated.json'
    state = json.loads(path.read_text()); L=33.5
    # Reconstruct Lagrange polynomials directly from saved controls, without
    # the observer's Curve, polynomial or inverse routines.
    def make_value(controls):
        segments=[]
        for edge in controls:
            p=np.asarray(edge,float)*L
            segments.append((p[0],-3*p[0]+4*p[1]-p[2],2*p[0]-4*p[1]+2*p[2],p[0,0],p[2,0]))
        def value(radius):
            for c0,c1,c2,lo,hi in segments:
                if lo-1e-10 <= radius <= hi+1e-10:
                    if abs(c2[0]) < 1e-13:
                        t=(radius-c0[0])/c1[0]
                    else:
                        roots=np.roots([c2[0],c1[0],c0[0]-radius])
                        valid=[x.real for x in roots if abs(x.imag)<1e-10 and -1e-9 <= x.real <= 1+1e-9]
                        if len(valid)!=1:raise RuntimeError('independent radial inverse ambiguous')
                        t=valid[0]
                    return float(c0[1]+c1[1]*t+c2[1]*t*t)
            raise ValueError('radius not covered')
        return value,[s[3] for s in segments]+[segments[-1][4]]
    top,kt=make_value(state['boundary_top_quadratic_controls'])
    bottom,kb=make_value(state['boundary_bottom_quadratic_controls'])
    R=float(state['parameters']['patch_radius'])*L
    cuts=sorted(set([0.,R]+[x for x in kt+kb if 0<x<R]))
    rows=[]
    for weighting in ('current_radial_line','current_disk'):
        out=observe_state(state,weighting=weighting,component_aggregation='regional_mean')
        radial_weight=(lambda r:1.) if weighting=='current_radial_line' else (lambda r:r)
        norm=R if weighting=='current_radial_line' else R*R/2
        for key,fn in [('apical_downward',lambda r:L-top(r)),('basal_upward',bottom),
                       ('regional_mean_Eulerian_height_shortening',lambda r:L-top(r)+bottom(r))]:
            estimate=sum(quad(lambda r:radial_weight(r)*fn(r),lo,hi,epsabs=1e-10,epsrel=1e-10)[0]
                         for lo,hi in zip(cuts,cuts[1:]))/norm
            bound=out['four_endpoints'][key]
            error=abs(estimate-float((bound.lo+bound.hi)/2))
            checks['independent_float_integral_'+weighting+'_'+key]=error < 1e-9
            rows.append({'weighting':weighting,'observable':key,'independent_quad':estimate,
                         'exact_interval_midpoint_display':float((bound.lo+bound.hi)/2),'difference':error})
    receipt={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,
             'observer_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
             'actual_field_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
             'actual_integral_crosschecks':rows,
             'analytic_case':{'top':'R=t+t^2,Z=10-t','bottom':'R=2t,Z=1-t',
                  'material_paired_shortening':'1','Eulerian_peak':'9/8','line_mean':'13/12','disk_mean':'43/40'},
             'scope':'Independent exact analytic/gauge cases and numerical-integral replay on one saved field; no fit, continuum-error or biological-observer equivalence claim.'}
    target=ROOT/'evidence/quadratic_endpoint_observer_independent_review.json'
    target.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'checks':len(checks),'failures':[k for k,v in checks.items() if not v],
                      'maximum_independent_integral_difference':max(r['difference'] for r in rows)}),flush=True)
    if not all(checks.values()):raise RuntimeError('independent observer checks failed')


if __name__=='__main__':
    run()
