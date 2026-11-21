"""Write-once version of the bounded safe independent bulk-cache AD replay."""
from pathlib import Path
import hashlib,json,datetime
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():raise ValueError('preserve '+str(p))
    p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf8')

SAFE='''def clean(d):
    out={}
    for key,value in d.items():
        if isinstance(value,arb):
            if not value.is_finite():raise ValueError('finite AD coefficient required')
            exact_zero=value.is_zero()
        else:exact_zero=F(value)==0
        if not exact_zero:out[key]=value
    return out
'''
TESTS='''def adversarial_tests():
    checks={};x=arb(0,1);half=P.A(F(1,2))
    checks['safe_zero_crossing_constant_retained']=ZERO in clean({ZERO:x}) and clean({ZERO:x})[ZERO].contains(half)
    checks['safe_exact_zero_removed']=clean({ZERO:arb(0)})=={}
    checks['safe_fraction_exact_zero_removed']=clean({ZERO:F(0)})=={}
    checks['safe_nonzero_fraction_retained']=clean({ZERO:F(2)})[ZERO]==2
    checks['safe_uncertain_addition_retained']=ZERO in add({ZERO:arb(1)},{ZERO:arb(-1,1)}) and add({ZERO:arb(1)},{ZERO:arb(-1,1)})[ZERO].contains(half)
    checks['safe_exact_cancellation_removed']=add({ZERO:arb(1)},{ZERO:arb(-1)})=={}
    checks['safe_uncertain_scaling_retained']=scale({ZERO:x},2)[ZERO].contains(P.A(1))
    checks['safe_uncertain_multiplication_retained']=mul({ZERO:x},{ZERO:arb(1)})[ZERO].contains(half)
    power={ZERO:arb(1)}
    for n in range(1,7):
        power=mul(power,{ZERO:x});checks['safe_uncertain_power_'+str(n)]=ZERO in power and power[ZERO].contains(P.A(F(1,2)**n))
    key=(1,0,0,0,0)
    checks['safe_uncertain_derivative_coefficient_retained']=derivative(clean({key:x}),[0]).contains(half)
    for label,value in [('nan',arb('nan')),('positive_inf',arb('+inf')),('negative_inf',arb('-inf'))]:
        try:clean({ZERO:value});checks['safe_nonfinite_rejected_'+label]=False
        except ValueError:checks['safe_nonfinite_rejected_'+label]=True
    assert all(v is True for v in checks.values()),checks
    return checks

'''

def main():
    old=ROOT/'runs/exact_integration_bridge_v1/cache_consistency_independent/review.py'
    src=old.read_text()
    needle='def clean(d):return{k:v for k,v in d.items()if v!=0}\n'
    assert needle in src
    src=src.replace(needle,SAFE)
    src=src.replace('def main():\n',TESTS+'''def main():
    safe_protocol=read(HERE/'protocol.json')
    assert sha(Path(__file__))==safe_protocol['review_source_sha256']
    assert all(sha(CF/k)==v for k,v in safe_protocol['input_hashes'].items())
''',1)
    src=src.replace('ctx.prec=384;P.core.A.cache_clear();protocol=read(OUT/\'protocol.json\');pre=read(OUT/\'preflight.json\');checks={};bindings={}',"ctx.prec=384;P.core.A.cache_clear();protocol=read(OUT/'protocol.json');pre=read(OUT/'preflight.json');checks=adversarial_tests();bindings={}")
    src=src.replace("'PASS_UNIFORM_LOCAL_BULK_CACHE_COMPONENT_WITH_GUARD_FINDING'", "'PASS_SAFE_V2_UNIFORM_LOCAL_BULK_CACHE_COMPONENT_WITH_GUARD_FINDING'")
    src=src.replace("'checks':checks,'source_sha256':EXPECTED", "'checks':checks,'review_protocol_sha256':sha(HERE/'protocol.json'),'supersedes_flawed_independent_review_sha256':safe_protocol['old_flawed_receipt_sha256'],'unsafe_old_T_polynomial_operations_used_with_Arb':False,'source_sha256':EXPECTED")
    src=src.replace("(HERE/'results.json').write_text", "assert all(sha(CF/k)==v for k,v in safe_protocol['input_hashes'].items())\n    assert not (HERE/'results.json').exists() and not (CF/'evidence/exact_cache_consistency_independent_review_v2.json').exists()\n    (HERE/'results.json').write_text")
    src=src.replace('evidence/exact_cache_consistency_independent_review.json','evidence/exact_cache_consistency_independent_review_v2.json')
    dest=HERE/'review.py'
    if dest.exists():raise ValueError('preserve new source')
    dest.write_text(src,encoding='utf8')
    original=ROOT/'runs/exact_integration_bridge_v1/cache_consistency_preparation/protocol.json'
    inputs=dict(read(original)['source_and_input_hashes'])
    paths=[old,original,original.with_name('preflight.json'),ROOT/'evidence/exact_cache_consistency_independent_review.json',ROOT/'runs/exact_integration_bridge_v1/bulk_pilot_independent_replay/replay.py',ROOT/'runs/exact_integration_bridge_v1/cache_consistency_preparation/local_element0_q0.json',ROOT/'runs/exact_integration_bridge_v1/cache_consistency_preparation/local_element6_q63.json',Path(__file__)]
    for p in paths:inputs[p.relative_to(ROOT).as_posix()]=sha(p)
    assert all(sha(ROOT/k)==v for k,v in inputs.items())
    p={'status':'FROZEN_BEFORE_BOUNDED_REPLAY','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'review_source_sha256':sha(dest),'input_hashes':inputs,'old_flawed_receipt_sha256':sha(ROOT/'evidence/exact_cache_consistency_independent_review.json'),'producer_sha256':sha(ROOT/'scripts/exact_cache_consistency.py'),'precision_bits':384,'predeclared_replay':'Original125 thirds + exact maps/two local uniform fixtures/16 endpoint comparisons/old boolean guards; extra18 cleanup containment guards. No integrations, whole-cache sum or solves.','dependency_limit':'Frozen T.material_basis uses Fraction only on these paths; I.direct_basis is Fraction only; neither uncertain-coefficient T polynomial arithmetic nor unsafe independent integration algorithms are used. Own AD clean rejects nonfinite and removes exact zeros only.','scope':'Producer dense analytic local bulk-cache component; independent generic third-order AD and four-term grouping. Both actual sites inactive;16 endpoint comparisons contain8 distinct Y values, not16 distinct active physical examples.'}
    write(HERE/'protocol.json',p)
    print(json.dumps({'source_sha256':sha(dest),'protocol_sha256':sha(HERE/'protocol.json'),'inputs':len(inputs)}))
if __name__=='__main__':main()
