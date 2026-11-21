"""Bounded counterexamples and dependency audit; NO spatial integration/solve."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib, importlib.util, json, sys, datetime, shutil
from flint import arb, ctx

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import exact_fe_target as T
import exact_bulk_bridge_pilot as P
import exact_uniform_bulk_defect_B as U

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def write(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():raise ValueError('preserve '+str(p))
    p.write_text(json.dumps(d,indent=2)+'\n',encoding='utf8')
def spec(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m

def main():
    ctx.prec=256;P.core.A.cache_clear();checks={};witnesses={}
    x=arb(0,1)
    checks['zero_crossing_ball_is_not_exact_zero']=not x.is_zero()
    checks['zero_crossing_ball_neither_equal_nor_unequal_zero']=(x==0) is False and (x!=0) is False
    checks['bug_drops_uncertain_constant']=T.constant(x)=={}
    checks['bug_drops_uncertain_addition']=T.add({(0,0):arb(1)}, {(0,0):arb(-1,1)})=={}
    checks['bug_drops_uncertain_product']=T.mul({(1,0):x},T.ONE)=={}
    checks['counterexample_actual_value_half_omitted']=T.evaluate(T.constant(x),Q(1),Q(1))==0 and x.contains(P.core.A(Q(1,2)))
    checks['safe_exact_zero_predicate_retains_uncertainty']=not x.is_zero() and not(x==0)
    witnesses['comparison_semantics']={'ball':P.core.pack(x),'eq_zero':x==0,'ne_zero':x!=0,'bool_ball':bool(x),'is_zero':x.is_zero(),'cleanup_output':str(T.constant(x)),'missed_admissible_value_exact':'1/2'}

    # Observe real POINT series arithmetic without computing a single moment.
    original=T.clean;stats={};phase=''
    def observed(poly):
        dropped=[(k,v)for k,v in poly.items()if isinstance(v,arb)and not v.is_zero()and not(v!=0)]
        if dropped:
            st=stats.setdefault(phase,{'calls_with_dropped_uncertainty':0,'coefficients_dropped':0,'example':None})
            st['calls_with_dropped_uncertainty']+=1;st['coefficients_dropped']+=len(dropped)
            if st['example'] is None:st['example']={'exponent':list(dropped[0][0]),'ball':P.core.pack(dropped[0][1])}
        return original(poly)
    T.clean=observed
    try:
        rec=read(ROOT/P.RECORD);target=T.Target(read(ROOT/rec['point_cache_file']))
        for element in [0,6]:
            phase=f'actual_point_element_{element}_series_degree_4'
            f=target.element(element);jet=P.JET.bulk_jacobian_jet(f['q'],f['B'])
            P.series(T.sub(jet['value'],T.ONE),4)
        box=read(U.GEOM/'coordinate_boxes.json');Y=U.field_nodes(box)
        phase='actual_uniform_element_0_fields_and_jet_only'
        f=target.element(0,Y);P.JET.bulk_jacobian_jet(f['q'],f['B'])
    finally:T.clean=original
    checks['real_point_series_encounters_dropped_uncertainty']=any(v['coefficients_dropped']>0 for k,v in stats.items()if 'point' in k)
    checks['real_uniform_fields_encounter_dropped_uncertainty']=stats.get('actual_uniform_element_0_fields_and_jet_only',{}).get('coefficients_dropped',0)>0
    witnesses['actual_arithmetic_only_drop_counts']=stats

    # Independently-owned old bulk-cache AD also has the unsafe predicate.
    old=spec('old_independent_cache_AD_for_cleanup_audit',ROOT/'runs/exact_integration_bridge_v1/cache_consistency_independent/review.py')
    checks['independent_bulk_cache_AD_has_same_failure']=old.clean({old.ZERO:x})=={}
    # Current surface-cache jclean is logically different and safe for this issue.
    surface=spec('surface_cache_cleanup_audit',ROOT/'scripts/exact_surface_cache_consistency.py')
    checks['current_surface_cache_jclean_retains_uncertainty']=surface.ZERO in surface.jclean({surface.ZERO:x})
    checks['current_surface_cache_jclean_removes_exact_zero']=surface.jclean({surface.ZERO:arb(0)})=={}

    paths=[Path(T.__file__),Path(P.__file__),Path(U.__file__),Path(P.JET.__file__),Path(old.__file__),Path(surface.__file__),ROOT/'scripts/exact_fe_moments.py',ROOT/'scripts/root_body_geometry.py',ROOT/'scripts/algebraic_root_feasibility.py']
    sources={str(p.relative_to(ROOT)).replace('\\','/'):sha(p)for p in paths}
    for p in paths:
        dest=OUT/'source_snapshots'/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True)
        if not dest.exists():shutil.copyfile(p,dest)
        assert sha(p)==sha(dest)
    receipt_paths=['exact_bulk_pilot_independent_review','exact_bulk_assembly_B_independent_review','exact_complete_point_B_independent_review','exact_integration_bridge_final_delivery','exact_integration_documentation_review','exact_cache_consistency_independent_review','exact_fe_moments_independent_review','exact_retube_B_geometry_independent_review','exact_retube_B_preparation_independent_review']
    receipts={name:sha(ROOT/'evidence'/f'{name}.json')for name in receipt_paths}
    report={'status':'CONFIRMED_ENCLOSURE_BUG_WITHDRAW_AFFECTED_INTEGRATION_CLAIMS','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'witnesses':witnesses,'source_sha256':sha(__file__),'frozen_sources':sources,'receipt_hashes_at_audit':receipts,'actual_spatial_integrals_or_equilibrium_solves_run':False,'scope':'Counterexamples, two actual degree4 POINT series and one uniform-field polynomial construction only. Historical bytes preserved. Existing independent reductions/hash checks do not validate invalid local enclosures.'}
    assert all(v is True for v in checks.values()),checks
    write(OUT/'results.json',report)
    print(json.dumps({'status':report['status'],'checks':len(checks),'stats':stats}),flush=True)
if __name__=='__main__':main()
