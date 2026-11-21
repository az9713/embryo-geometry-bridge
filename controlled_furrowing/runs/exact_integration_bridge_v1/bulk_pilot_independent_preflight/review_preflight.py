"""Prelaunch math/glue audit: synthetic fixtures only; never run_element/main."""
from pathlib import Path
from fractions import Fraction as F
from types import SimpleNamespace
import hashlib,json,sys
import numpy as np
import sympy as s
from flint import arb,ctx
H=Path(__file__).resolve().parent;C=H.parents[2]
EXPECTED='91a4f9a41b69c7b9e9ae5d5e0021087a7dfa4cf9e80625273bb4945d0ac3d0b3'
PROTOCOL='746b260f5fbe7d049958681533165784e4847a15ca7d302ae9cdf5712cda95f0'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    source=C/'scripts/exact_bulk_bridge_pilot.py';pfile=C/'runs/exact_integration_bridge_v1/bulk_pilot/protocol_v2.json'
    if sha(source)!=EXPECTED or sha(pfile)!=PROTOCOL:raise ValueError('frozen source/protocol changed')
    sys.path.insert(0,str(C/'scripts'));import exact_bulk_bridge_pilot as P
    ctx.prec=512;checks={};protocol=json.loads(pfile.read_text());inputs={}
    for name in ['source_and_input_hashes','target_preparation_input_hashes']:
        for f,digest in protocol[name].items():
            if sha(C/f)!=digest:raise ValueError('changed input '+f)
            canonical=(C/f).resolve().relative_to(C.resolve()).as_posix()
            if canonical in inputs and inputs[canonical]!=digest:raise ValueError('conflicting binding')
            inputs[canonical]=digest
    def ck(k,v):checks[k]=bool(v)
    ck('all_frozen_inputs_verified',bool(inputs));ck('four_predeclared_elements_and_policy',protocol['elements']==[0,6,72,66]and protocol['truncation_order']==16 and protocol['precision_bits']==512 and protocol['rho_exact']=='21/500')
    ck('point_only_scope_explicit','no parameter/tube uncertainty' in protocol['input_kind']and protocol['target_bridge_proved']is False)
    off={'chart':{'kind':'off_axis'},'corners':[(F(1),F(-2)),(F(3),F(0)),(F(1),F(2))],'zcenter':F(0)}
    poly={(-2,3):F(-7,3),(2,0):F(5,7)}
    want=F(7,3)*8+F(5,7)*9
    ck('mixed_Laurent_majorant_independent',P.pnorm(poly,off)==want)
    axis={'chart':{'kind':'axis_vertex'}}
    ck('axis_polynomial_majorant_unit_square',P.pnorm({(2,1):F(-3,7),(0,8):F(1,9)},axis)==F(3,7)+F(1,9))
    # Series coefficient and signed dual applied to a linear synthetic Jacobian.
    Ap,Bp=P.series({(1,0):F(1,10)},2)
    want1=[F(-1),F(5001,10),F(-7501,100)]
    want2=[F(5001),F(-15002,10),F(27503,100)]
    def contains_exact(ball,value):return F(str(ball.lower().fmpq()))<=value<=F(str(ball.upper().fmpq()))
    ck('synthetic_phi_prime_series_coefficients',all(contains_exact(Ap[(i,0)],v)for i,v in enumerate(want1)))
    ck('synthetic_phi_second_series_coefficients',all(contains_exact(Bp[(i,0)],v)for i,v in enumerate(want2)))
    toy=SimpleNamespace(energy_scale=F(1),normalization=lambda:2*arb.pi())
    fields={'corners':[(F(1),F(0)),(F(2),F(0)),(F(1),F(1))],'chart':{'kind':'off_axis'},'zcenter':F(0)}
    rule={'volume_barycentric':[[F(1,3)]*3],'volume_weights':[[F(2,7)]]}
    functional=P.Functional(toy,fields,0,rule)
    independent=2*arb.pi()*P.core.A(F(2,3))-P.core.A(F(2,7))
    ck('weighted_constant_moment_measure_once',functional.defect(0,0).overlaps(independent))
    ck('signed_functional_linear_combination',functional.apply({(0,0):F(-3),(1,0):F(2)}).overlaps(-3*functional.defect(0,0)+2*functional.defect(1,0)))
    for a,b in [(-5,2),(-1,0),(0,4),(7,3)]:
        exact_rule=F(2,7)*F(4,3)**a*F(1,3)**b
        known=P.M.normalized_bulk_moment(fields['corners'],a,b,1,0,512)-P.core.A(exact_rule)
        ck(f'outward_auxiliary_power_cache_{a}_{b}',functional.defect(a,b).overlaps(known))
    ck('two_mass_tail_covers_opposite_measure_errors',P.core.upper(functional.mass_integral+functional.mass_rule)>P.core.upper(functional.mass_integral-functional.mass_rule))
    # POINT cache correction, synthetic two-free-DOF setup. Independent direct
    # scalar-energy differentiation checks fixed-coordinate and lift semantics.
    nodes=[[F(1),F(0)],[F(2),F(0)],[F(1),F(1)],[F(3,2),F(0)],[F(3,2),F(1,2)],[F(1),F(1,2)]]
    grads=np.zeros((1,1,6,2));grads[0,0,:3,0]=[-1,1,0];grads[0,0,:3,1]=[-1,0,1]
    arrays={'triangles':np.array([[0,1,2,3,4,5]]),'shape_grad':grads,'shape':np.array([[.5,.5,0,0,0,0]]),'qr':np.array([[1.5]])}
    delta=F(1,1000);localB=[[-1+delta,0],[-1,0],[0,1],[0,0],[F(1,3),0]]
    fields2={'free_local':[0,3],'B':[[P.T.constant(F(v))for v in row]for row in localB],'q':[P.T.constant(v)for v in [1+delta,0,0,1,1]]}
    target=SimpleNamespace(a=arrays,nodes=lambda:nodes);fun2=SimpleNamespace(points=[((F(4,3),F(1,3)),F(2,7))]);drive=F(2,5)
    g,Hess=P.cache_correction(target,fields2,0,fun2,drive)
    alpha,beta=s.symbols('alpha beta');qa=[1-alpha,-alpha,beta,1,1+alpha/3];qb=[1+s.Rational(1,1000)+(-1+s.Rational(1,1000))*alpha,-alpha,beta,1,1+alpha/3]
    def energy(q):
        J=(q[0]*q[3]-q[1]*q[2])*q[4]
        return sum(v*v for v in q)/2-s.log(J)+2500*s.log(J)**2+s.Rational(1,5)*(q[1]**2+q[3]**2)
    E=s.Rational(2,7)*(energy(qb)-energy(qa))
    def evaluate(expr):
        if expr.is_Rational:return P.core.A(F(int(expr.p),int(expr.q)))
        if expr.func==s.log:return evaluate(expr.args[0]).log()
        if expr.is_Add:return sum((evaluate(x)for x in expr.args),arb(0))
        if expr.is_Mul:
            out=arb(1)
            for x in expr.args:out*=evaluate(x)
            return out
        if expr.is_Pow and expr.args[1].is_Integer:return evaluate(expr.args[0])**int(expr.args[1])
        raise ValueError('unexpected symbolic expression')
    for i,x in enumerate([alpha,beta]):
        ck('synthetic_cache_correction_gradient_'+str(i),g[i].overlaps(evaluate(s.diff(E,x).subs({alpha:0,beta:0}))))
        for j,y in enumerate([alpha,beta]):ck(f'synthetic_cache_correction_Hessian_{i}_{j}',Hess[i,j].overlaps(evaluate(s.diff(E,x,y).subs({alpha:0,beta:0}))))
    ck('source_protocol_unchanged_after_synthetic_checks',sha(source)==EXPECTED and sha(pfile)==PROTOCOL)
    out={'status':'PASS_POINT_PILOT_PRELAUNCH_REVIEW'if all(checks.values())else'FAIL_PRELAUNCH_REVIEW','checks':checks,'driver_sha256':EXPECTED,'protocol_sha256':PROTOCOL,'review_source_sha256':sha(Path(__file__)),'frozen_input_hashes':inputs,'mathematical_review':['Integral-minus-Qstar uses normalized axisymmetric material measure exactly once and signed monomial defects.','Both positive integral and discrete masses multiply valid uniform series remainder bounds.','Shear/active polynomial gradient/Hessian and volumetric chain rule are separated correctly.','FiniteQstar-minus-cache correction is a direct outward-rounded point computation using identical exactX+u and projected local free coordinates; no tube uniformity asserted.','Existing exactP2/P1 wholebody geometry box includes chosen exactcenter and coherent interiorauxiliarynodes; globalrho21/500 justifies logseries atbothintegral andQstarpoints.'],'source_target_scope':'RecordB saved exactcenter only, four predeclared bulk elements, no surfaces or assembled611-vector/matrix or any root transfer.','launch_condition':'Only exact frozen driver/protocol reviewed. Driver does not itself require reviewreceipt; parent launch must externally bind this receipt and these hashes.','actual_integration_executed_by_reviewer':False,'root_transfer_proved':False}
    (H/'preflight.json').write_text(json.dumps(out,indent=2)+'\n');(C/'evidence/exact_bulk_pilot_prelaunch_independent_review.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'checks':len(checks),'failed':[k for k,v in checks.items()if not v],'receipt_sha256':sha(C/'evidence/exact_bulk_pilot_prelaunch_independent_review.json')}))
    if not all(checks.values()):raise SystemExit(1)
if __name__=='__main__':main()
