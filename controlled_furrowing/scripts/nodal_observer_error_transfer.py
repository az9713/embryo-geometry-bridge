"""Conditional exact P2 nodal-coordinate error -> Eulerian endpoint errors.

An actual verified nodal bound must be supplied. This module does not invent
one from a force residual or a mesh comparison. Constants are conservative.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json
from annulus_observer import _quadratic_controls


def curve_bound(controls,nodal_error_um,*,fixed_radial_endpoints=False,roi_um=None):
    epsilon=Q(nodal_error_um)
    if epsilon<0:raise ValueError('nonnegative verified-or-explicitly-assumed nodal error required')
    curves=_quadratic_controls(controls)
    radial=min(min(c1[0],c1[0]+2*c2[0]) for q,c0,c1,c2 in curves)
    if radial<=8*epsilon:
        return {'status':'UNRESOLVED_RADIAL_MONOTONICITY','minimum_radial_derivative_exact':str(radial),
                'required_margin_exact':str(8*epsilon),'scope':'No finite bound accepted when the sufficient strict margin fails.'}
    e=Q(5,4)*epsilon
    domain=(curves[0][0][0][0],curves[-1][0][-1][0])
    if roi_um is None:roi_um=domain
    left,right=map(Q,roi_um)
    if left>right or left<domain[0] or right>domain[1]:raise ValueError('base ROI not covered')
    if not fixed_radial_endpoints and (left<domain[0]+e or right>domain[1]-e):
        return {'status':'UNRESOLVED_COMMON_RADIAL_DOMAIN','scope':'Supply verified fixed radial endpoints or an ROI inside the guaranteed common domain.'}
    slope=max(max(abs(c1[1]),abs(c1[1]+2*c2[1]))/min(c1[0],c1[0]+2*c2[0]) for q,c0,c1,c2 in curves)
    error=e*(1+slope)
    return {'status':'CONDITIONAL_EXACT_ERROR_TRANSFER','nodal_coordinate_error_um_exact':str(epsilon),
            'parametric_coordinate_error_um_exact':str(e),'radial_derivative_margin_exact':str(radial-8*epsilon),
            'base_graph_Lipschitz_upper_exact':str(slope),'Eulerian_vertical_error_um_exact':str(error),
            'fixed_radial_endpoints_assumed':bool(fixed_radial_endpoints),'ROI_um_exact':[str(left),str(right)],
            'preconditions':['Corresponding conforming P2 nodes each differ by at most epsilon in BOTH coordinates.',
                'Only the stated fixed ROI, normalized nonnegative weights and unchanged reference profiles are compared.',
                'Fixed endpoint radii, if selected, require independent verification; they are not inferred from this function.'],
            'scope':'Exact transfer of an input error bound; not proof that any mechanical root/PDE/experiment satisfies that input.'}


def endpoint_bounds(apical,basal,nodal_error_um,*,fixed_radial_endpoints=False,roi_um=None):
    a=curve_bound(apical,nodal_error_um,fixed_radial_endpoints=fixed_radial_endpoints,roi_um=roi_um)
    b=curve_bound(basal,nodal_error_um,fixed_radial_endpoints=fixed_radial_endpoints,roi_um=roi_um)
    if any(v['status']!='CONDITIONAL_EXACT_ERROR_TRANSFER' for v in (a,b)):
        return {'status':'UNRESOLVED_ERROR_TRANSFER','apical':a,'basal':b}
    if a['ROI_um_exact']!=b['ROI_um_exact']:
        return {'status':'UNRESOLVED_COMMON_SURFACE_ROI','apical':a,'basal':b,
                'scope':'Both surfaces must be compared on the same explicit ROI before height-summary errors can be combined.'}
    if Q(a['ROI_um_exact'][0])>=Q(a['ROI_um_exact'][1]):
        raise ValueError('positive-width common ROI required for normalized regional means')
    ea=Q(a['Eulerian_vertical_error_um_exact']);eb=Q(b['Eulerian_vertical_error_um_exact'])
    return {'status':'CONDITIONAL_EXACT_ERROR_TRANSFER','apical':a,'basal':b,
            'endpoint_error_um_exact':{'maximum_shortening_um':str(ea+eb),
                'region_averaged_shortening_um':str(ea+eb),
                'apical_max_shortening_um':str(ea),'basal_max_shortening_um':str(eb)},
            'scope':'Both line/disk means and spatial-max components at ONE static state. No time-supremum or changing-mask/reference bound.'}


def tests():
    graph=lambda slope,intercept=0:[[(Q(t),Q(intercept)+Q(slope)*Q(t)) for t in (Q(0),Q(1,2),Q(1))]]
    checks={};e=Q(1,1000)
    a=curve_bound(graph(4),e,fixed_radial_endpoints=True)
    checks['tilted_graph_exact_bound']=Q(a['Eulerian_vertical_error_um_exact'])==Q(25,4)*e
    flat=curve_bound(graph(0),e,fixed_radial_endpoints=True)
    checks['flat_graph_exact_bound']=Q(flat['Eulerian_vertical_error_um_exact'])==Q(5,4)*e
    out=endpoint_bounds(graph(4,10),graph(-2),e,fixed_radial_endpoints=True)
    checks['height_bound_is_sum_of_components']=Q(out['endpoint_error_um_exact']['maximum_shortening_um'])==10*e
    checks['do_not_assume_fixed_ROI_coverage']=curve_bound(graph(1),e)['status']=='UNRESOLVED_COMMON_RADIAL_DOMAIN'
    checks['sufficient_interior_ROI_coverage']=curve_bound(graph(1),e,roi_um=(Q(1,10),Q(9,10)))['status']=='CONDITIONAL_EXACT_ERROR_TRANSFER'
    checks['radial_margin_failure_unresolved']=curve_bound(graph(1),Q(1,8),fixed_radial_endpoints=True)['status']=='UNRESOLVED_RADIAL_MONOTONICITY'
    checks['zero_error_is_exact']=curve_bound(graph(1),0,fixed_radial_endpoints=True)['Eulerian_vertical_error_um_exact']=='0'
    longer=[[(2*r,z) for r,z in edge] for edge in graph(1)]
    checks['different_default_surface_ROIs_unresolved']=endpoint_bounds(graph(1),longer,e,fixed_radial_endpoints=True)['status']=='UNRESOLVED_COMMON_SURFACE_ROI'
    # Sharp trace constants at independent analytical extremizers.
    t=Q(1,4);N=[(1-t)*(1-2*t),4*t*(1-t),t*(2*t-1)]
    checks['sharp_value_constant']=sum(map(abs,N))==Q(5,4)
    checks['sharp_derivative_constant']=abs(-3)+abs(4)+abs(-1)==8
    for i,call in enumerate([lambda:curve_bound(graph(1),-1),lambda:curve_bound(graph(1),0,roi_um=(-1,1)),
                            lambda:endpoint_bounds(graph(0),graph(0),0,fixed_radial_endpoints=True,roi_um=(Q(1,2),Q(1,2)))]):
        try:call();checks['invalid_input_'+str(i)]=False
        except ValueError:checks['invalid_input_'+str(i)]=True
    receipt={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'tilted_example':a,'two_surface_example':out,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'annulus_dependency_sha256':hashlib.sha256(Path(__file__).with_name('annulus_observer.py').read_bytes()).hexdigest(),
        'scope':'Analytical exact synthetic checks of conditional error-transfer formulas; no actual root/PDE error radius supplied.'}
    path=Path(__file__).resolve().parents[1]/'evidence/nodal_observer_error_transfer_tests.json'
    path.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'checks':len(checks)}))
    if receipt['status']!='PASS':raise RuntimeError('conditional error-transfer checks failed')


if __name__=='__main__':tests()
