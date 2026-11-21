"""Bind/read-only check of the design audit; no exact-integral claim."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
import sympy as S

HERE=Path(__file__).resolve().parent;CF=HERE.parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    scale=json.loads((HERE/'scales.json').read_text()); checks={}
    checks['all_scale_inventory_checks_PASS']=all(v is True for v in scale['checks'].values())
    checks['all_bound_inputs_still_current']=all(sha(CF/k)==v for k,v in scale['bindings'].items())
    checks['scale_script_hash_current']=sha(HERE/'audit_scales.py')==scale['source_sha256']
    checks['all17_error_budgets_have_strict_spare_margin']=True
    checks['all17_transferred_Hessian_floor_positive_under_budget']=True
    for row in scale['records']:
        r=Q(row['radius_exact']);eta=Q(row['eta_exact']);Z=Q(row['Z_exact']);slack=Q(row['remaining_selfmap_slack_exact']);N=Q(row['R_infinity_exact'])
        a=slack/4;b=Q(row['sufficient_preconditioned_Hessian_defect_budget_exact'])
        checks['all17_error_budgets_have_strict_spare_margin'] &= eta+a+(Z+b)*r==r-slack/2 and Z+b<1
        checks['all17_transferred_Hessian_floor_positive_under_budget'] &= (1-Z-b)/N>0
    s,t=S.symbols('s t',positive=True)
    checks['axis_vertex_Duffy_example_direction_dependent_but_smooth']=S.cancel((s+s*t)/s)==1+t
    checks['axis_vertex_example_not_one_corner_limit']=(1+t).subs(t,0)!=(1+t).subs(t,1)
    R,lp,r0,K,T=S.symbols('R lp r0 K T',positive=True)
    # lp includes the physical-reference tangent-length scale.
    js=R*lp/r0
    checks['weighted_surface_energy_cancellation_exact']=S.simplify(r0*(K*(js-1)**2/2+T*js)-(K*(R**2*lp**2/r0-2*R*lp+r0)/2+T*R*lp))==0
    checks['no_exact_integration_certificate_claimed']=scale['exact_integration_bridge_certified'] is False and scale['no_error_bound_computed'] is True
    checks['declared_all17_144_element_inventory']=len(scale['records'])==17 and sum(scale['geometry']['triangle_counts'][k] for k in ['axis_edge_triangles','axis_vertex_only_triangles','off_axis_triangles'])==144
    owned=[HERE/'audit_scales.py',HERE/'scales.json',HERE/'DESIGN-REVIEW.md',Path(__file__).resolve()]
    sources=['scripts/axisymmetric_quadratic.py','scripts/axisymmetric_mixed.py','scripts/axisymmetric_mechanics.py','scripts/quadratic_quadrature_study.py','scripts/algebraic_root_feasibility.py','scripts/root_body_geometry.py','scripts/cohort_energy_stability.py']
    bindings={str(p.relative_to(CF)).replace('\\','/'):sha(p) for p in owned}
    bindings.update({k:sha(CF/k) for k in sources})
    result={'status':'PASS_BOUNDED_INDEPENDENT_DESIGN_REVIEW' if all(checks.values()) else 'FAIL_DESIGN_REVIEW','checks':checks,'artifact_sha256':bindings,'scale_input_bindings':scale['bindings'],'exact_integration_bridge_status':'NOT_IMPLEMENTED_OR_CERTIFIED_BY_THIS_AUDIT','review_scope':'Proof obligations and exact read-only arithmetic/computability audit, not a quadrature error calculation or theorem of successful persistence.','criteria':{'gradient_defect':'a >= sup_theta ||R(FI-FQ)(p_theta)||infinity','Hessian_defect':'b >= sup_theta,u_in_B ||R(HI-HQ)(u)||infinity','acceptance':'Z+b<1 and eta+a+(Z+b)*r<r','minimum':'sym(R)>0 with transferred contraction gives HI >=(1-Z-b)/||R||infinity I'},'quantitative_extrema':scale['extrema'],'author_method_comparison':'Proposed axis-regularized spatial series/exact-moment method agrees in principle. Binary rule moment defects and all cache-map rounding must be retained. No implementation yet reviewed.','next_required_execution':['Freeze exactly integrated FE target and source-bound rule/map consistency decomposition.','Implement validated variation integrals with explicit axis cancellation and nonzero binary moment defects.','Compute actual a,b on representative G/L and other15tubes; retain failures.','If necessary create separately reviewed enlarged/recentered/energy-norm proof domains, then revalidate geometry and all observer error transport.'],'scope_limits':['No new equilibrium optimization, Hessian assembly or integral-bound computation.','No PDE/discretization error, source temporal validation, complete parameter/geometry domain, new held-out prediction or confirmed novelty.','No external contact or production/documentation changes.','Original full five-gate project remains incomplete.']}
    out=CF/'evidence/exact_integration_design_review.json';out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8');(HERE/'verification.json').write_bytes(out.read_bytes())
    print(json.dumps({'status':result['status'],'checks':len(checks),'failed':[k for k,v in checks.items()if not v],'receipt_sha256':sha(out)}))
    if not all(checks.values()):raise SystemExit(1)
if __name__=='__main__':main()
