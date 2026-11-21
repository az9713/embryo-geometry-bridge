"""Executable exact tests for owned sparse-jet/dual helper; no integration."""
from pathlib import Path
from fractions import Fraction as F
import json,sys
sys.path.insert(0,str(Path(__file__).parent))
import polynomial_jet as J
import verify_derivatives as V
T=J.T;ROOT=J.ROOT;OUT=Path(__file__).parent

def main():
    destination=OUT/'polynomial_jet_tests.json'
    if destination.exists():raise ValueError('write-once result exists')
    protocol=T.read(ROOT/'runs/exact_integration_bridge_v1/target_preparation/protocol.json');row=protocol['records'][0];target=T.Target(T.read(ROOT/row['point_cache_file']));probe=T.read(OUT/'results.json');checks={}
    def ck(k,v):
        checks[k]=bool(v)
        if not v:raise AssertionError(k)
    for case in probe['bulk']:
        f=target.element(case['element']);jets=[V.Jet(q,b)for q,b in zip(f['q'],f['B'])];ad=(jets[0]*jets[3]-jets[1]*jets[2])*jets[4];got=J.bulk_jacobian_jet(f['q'],f['B'])
        ck('independent_polynomial_AD_'+case['kind'],ad.v==got['value']and ad.g==got['gradient']and ad.h==got['hessian_upper'])
        if case['kind'].startswith('off_axis'):
            for scaled in(True,False):
                rep=J.shifted_numerator(got['value'],f['corners'],scaled);ck(f"exact_radial_roundtrip_{case['kind']}_{scaled}",J.recover_physical_R(rep)==got['value'])
        z=T.sub(got['value'],T.ONE);z2=T.mul(z,z);A=T.add(T.add(T.ONE,z),z2);B=T.add(T.scale(T.ONE,2),T.add(T.scale(z,3),T.scale(z2,4)))
        # Any signed exact monomial functional must distribute through products.
        def delta(a,b):return F(2)**a*F(1,3)**b-F(2)*F(3)**a*F(-2,7)**b
        def direct(poly):return T.evaluate(poly,F(2),F(1,3))-F(2)*T.evaluate(poly,F(3),F(-2,7))
        d=J.contract_volumetric_defect(got,A,B,delta)
        ck('dual_equals_direct_gradient_'+case['kind'],d['gradient']==[direct(T.mul(A,p))for p in got['gradient']])
        ck('dual_equals_direct_all_H_'+case['kind'],all(d['hessian_upper'][k]==direct(T.add(T.mul(A,p),T.mul(B,T.mul(got['gradient'][k[0]],got['gradient'][k[1]]))))for k,p in got['hessian_upper'].items()))
    for case in probe['surfaces']:
        f=target.surface(case['surface'],case['edge']);q=f['q_physical_R'];B=f['B_physical_R'];jets=[V.Jet(v,b)for v,b in zip(q,B)];ad=jets[0]*jets[0]*(jets[1]*jets[1]+jets[2]*jets[2]);got=J.surface_squared_jacobian_jet(q,B)
        ck(f"independent_surface_polynomial_AD_{case['surface']}_{case['edge']}",ad.v==got['value']and ad.g==got['gradient']and ad.h==got['hessian_upper'])
    try:J.shifted_numerator(T.X,[[F(0),F(0)],[F(1),F(0)],[F(0),F(1)]]);ck('axis_shift_rejected',False)
    except ValueError:checks['axis_shift_rejected']=True
    try:J.contract_factor(T.X,{(0,0):F(1)});ck('missing_dual_support_rejected',False)
    except ValueError:checks['missing_dual_support_rejected']=True
    ck('frozen_target_unchanged',T.sha(ROOT/'scripts/exact_fe_target.py')==J.TARGET_SHA256)
    result={'status':'PASS_SPARSE_POLYNOMIAL_JET_HELPER','checks':checks,'helper_source_sha256':T.sha(Path(J.__file__)),'test_source_sha256':T.sha(__file__),'generic_AD_verifier_sha256':T.sha(Path(V.__file__)),'probe_results_sha256':T.sha(OUT/'results.json'),'exact_target_source_sha256':J.TARGET_SHA256,'scope':'Exact local polynomial derivatives, reversible shifted rational representation, and moment-dual linearity for all local gradient/H entries of four actual B-center elements and six surfaces. Signed point-functional tests are algebra checks, not physical integrals or error certificates. No high-order series coefficients beyond degree2 in J-1, full integration or solve.'}
    T.write(destination,result);print(json.dumps({'status':result['status'],'checks':len(checks),'helper_source_sha256':result['helper_source_sha256'],'receipt_sha256':T.sha(destination)}))
if __name__=='__main__':main()
