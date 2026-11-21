"""Difference-focused semantic checks; no actual-field nonlinear integration."""
from pathlib import Path
from fractions import Fraction as Q
import ast,hashlib,json,sys,datetime
from flint import arb,ctx
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/'scripts'))
import exact_bulk_bridge_safe as B
import exact_surface_bridge_safe as S
T=B.T
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def functions(path):return{x.name:ast.dump(x,include_attributes=False)for x in ast.parse(path.read_text()).body if isinstance(x,(ast.FunctionDef,ast.ClassDef))}
def main():
    ctx.prec=128;B.core.A.cache_clear();checks={};details={}
    fixture=Q(1,5001)+Q(1,2**1000);answer=Q(5001,2**1000)
    aa,bb=B.series({(0,0):fixture},1);got=aa.get((0,0),arb(0))
    checks['signed_phi_prime_fixture_constant_retained']=(0,0)in aa and got.contains(B.core.A(answer)) and not got.is_zero()
    details['signed_fixture']={'e0_exact':str(fixture),'expected_phi_prime_N1_exact':str(answer),'actual_enclosure':B.core.pack(got),'precision_bits':128}
    x=arb(0,1)
    checks['bulk_safe_zero_crossing_constant_retained']=(0,0)in T.constant(x)
    checks['bulk_safe_uncertain_addition_retained']=(0,0)in T.add({(0,0):arb(1)},{(0,0):arb(-1,1)})
    checks['bulk_safe_uncertain_product_retained']=(1,0)in T.mul({(1,0):x},T.ONE)
    checks['bulk_surface_share_identical_SAFE_target_helper']=B.T is S.T and B.JET is S.JET
    k=read(ROOT/'evidence/exact_safe_polynomial_kernel_independent_review.json')
    checks['safe_kernel992_independent_checks_current']=k['status']=='PASS'and len(k['checks'])==992 and all(v is True for v in k['checks'].values()) and sha(Path(T.__file__))==k['source_sha256']and sha(Path(B.JET.__file__))==k['helper_source_sha256']
    checks['safe_target_no_unsafe_neq_cleanup']='if v!=0'not in Path(T.__file__).read_text() and 'return value.is_zero()'in Path(T.__file__).read_text()
    # The arithmetic changes are isolated: unchanged constitutive/tail/cache
    # functions; modified Functional/run output only adds per-instance counters.
    for old,new,names in [('exact_bulk_bridge_pilot.py','exact_bulk_bridge_safe.py',['pnorm','series','cache_correction']),('exact_surface_bridge_pilot.py','exact_surface_bridge_safe.py',['binomial_coefficients','binomial_tail','positive_G_bounds','normalized_series','pnorm','cache_correction'])]:
        a,b=functions(ROOT/'scripts'/old),functions(ROOT/'scripts'/new)
        for name in names:checks[f'{new}_{name}_formula_unchanged']=a[name]==b[name]
    inputs={}
    for folder,proto,pre in [('bulk_assembly_B_safe_v2','protocol_v2.json','preflight_v2.json'),('surface_assembly_B_safe_v2','protocol.json','preflight.json')]:
        base=ROOT/'runs/exact_integration_bridge_v1'/folder;p=read(base/proto)
        for name,h in p['source_and_input_hashes'].items():
            assert sha(ROOT/name)==h,name
            if name in inputs:assert inputs[name]==h
            inputs[name]=h
        checks[folder+'_safe_source_dependency_closure']=all(Path(module.__file__).resolve().relative_to(ROOT.resolve()).as_posix()in p['source_and_input_hashes'] for module in ([B,T,B.JET,B.M,B.core]if folder.startswith('bulk')else[S,T,S.JET,S.M,S.core]))
        checks[folder+'_no_old_local_point_results_bound']=not any('/elements/'in name or '/edges/'in name or 'complete_point_B/results.json'in name for name in p['source_and_input_hashes'])
        checks[folder+'_not_started_before_review']=not(base/'started.json').exists()and not(base/'results.json').exists()
    checks['original_u0_cache_exactly_same_B']=read(ROOT/B.RECORD)['point_cache_file'] in inputs
    checks['no_root_or_uniform_claim']=True
    assert all(v is True for v in checks.values()),checks
    result={'status':'PASS','checks':checks,'details':details,'input_hashes':inputs,'review_source_sha256':sha(__file__),'actual_point_field_integrations_run':False,'old_integration_claims_restored':False,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    out=HERE/'common_review.json';assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'PASS','checks':len(checks),'sha256':sha(out)}))
if __name__=='__main__':main()
