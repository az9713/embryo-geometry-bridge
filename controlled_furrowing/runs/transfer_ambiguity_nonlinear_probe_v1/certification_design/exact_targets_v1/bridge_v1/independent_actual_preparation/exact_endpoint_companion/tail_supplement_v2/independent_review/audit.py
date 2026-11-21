"""Independent cold arithmetic check of the point-force tail companion."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
from unittest.mock import patch

import flint
from flint import arb, ctx, fmpq

SOURCE = Path(__file__).resolve()
HERE = SOURCE.parents[2]
PINS = {
    "tighten_v2.py":"d3d7b931cbf31400215f5c1ce76c8cd4e759dd1e8b47a846dbfb618db53df79d",
    "protocol_v2.json":"2f3e1866b22a7c2fea817ddd74e02f5f99ac7b44abf5f3518405200251f01e15",
    "preflight_v2.json":"ba4685c386508e94c458ccae304e77e79857435bce403bb30206ff09c31ea6dc",
}

def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def reread_module():
    for name,v in PINS.items():
        assert sha(HERE/name)==v
    spec=importlib.util.spec_from_file_location("_cold_exact_tail_companion",HERE/"tighten_v2.py")
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def fraction(x):
    assert x.is_exact()
    return Q(str(x.fmpq()))


def lookup_matrix_rollbacks(m):
    ab,_,_=m.environment()
    assert Path(ab.__file__).resolve()==m.PREP/"replay_math.py"
    assert ab.A is ab.E.U.A and ab.E.core.ctx is ctx and ab.SUR.S.ctx is ctx
    phase=m.phase()
    pf=json.loads((HERE/"preflight_v2.json").read_bytes())
    assert pf["status"]=="PASS" and pf["force_H_cache_or_moment_work"]==pf["actual_raw_component_packets_opened"]==0
    assert phase==json.loads((HERE/"protocol_v2.json").read_bytes())
    assert phase["precision_bits"]==640 and phase["runtime_python_flint"]==flint.__version__
    p=json.loads((m.PREP/"protocol.json").read_bytes())
    assert p["precision_bits"]==640 and p["math_source_sha256"]==sha(Path(ab.__file__))
    for name,digest in p["source_and_premise_hashes"].items():
        assert sha(m.ROOT/name)==digest
    return ab,p


def same(out,x):
    assert out.mid().is_exact() and out.rad().is_exact()
    return (fraction(out.mid()),fraction(out.rad()))==(fraction(x.mid()),fraction(x.rad()))


def safe_erosion(m,A):
    t=Q(1,10**21)/3
    B=A.A(t)
    S=arb(0,A.E.core.upper(B))
    assert S.mid().is_zero()
    s=fraction(S.rad())
    assert s==m.exact_stored_radius(A,B) and s>t
    raw=Q(1,17)
    x=arb(fmpq(raw.numerator,raw.denominator))
    assert not x.is_exact()
    v1=[(x+arb(0))+S,(x+S)+arb(0)]
    for J in v1:
        outer=A.ends(J)
        strict,s1,t1=m.tighten(A,outer,t,Q(1))
        assert s1==s and t1==t and outer[0]<strict[0]<=raw-t<=raw+t<=strict[1]<outer[1]
    # Width is necessary but cannot replace proof of the actually inserted radius.
    # J contains X+[-s,s] for X={0}; a larger guessed radius can still satisfy
    # the width check and incorrectly move the left endpoint past -t.
    fudge=(s-t)/2
    J=(-s,s+2*fudge)
    wrong=m.add_back(J,s+fudge,t)
    assert wrong[0]>-t and J[1]-J[0]==2*(s+fudge)
    # The independent tube may still extend beyond a producer interval.
    right=m.add_back(J,s,t)
    assert right[0]==-t and right[1]>=t
    assert not A.contain((right[0]+fudge,right[1]),right)


def trace_literal_growth(m,A):
    # Synthetic off-axis, curved P2 reference trace. Only call the zeroth
    # mass primitive, with no force/H moments evaluation.
    X=[(Q(1),Q(0)),(Q(3,2),Q(0)),(Q(2),Q(0))]
    Y=[(Q(1),Q(1,20)),(Q(3,2),Q(1,12)),(Q(2),Q(1,10))]
    class Dummy:
        def normalization(self):
            return 2*arb.pi()/A.A(Q(7,3))
        def surface(self,f,e,nodes):
            assert (f,e)==(0,0) and nodes==Y
            f=A.T.surface_fields(X,Y,[0,1,2,3,4,5])
            return {**f,"global_indices":list(range(6))}
    target=Dummy()
    rule={"surface_s":[[Q(i,13) for i in range(1,13)]],
          "surface_weights":[[[Q(i,23) for i in range(1,13)]]]}
    case={"global_indices":list(range(6)),"domain":{"q0_exact":"7/6","rho_exact":"1/64"},"order":12}
    params=(Q(13,11),Q(2,17))
    fields=target.surface(0,0,Y)
    old=A.SUR.S.SurfaceFunctional.integral
    def zero_only(self,a,b):
        assert a==b==0, "cold trace touched a non-zero moment"
        return old(self,a,b)
    with patch.object(A.SUR.S.SurfaceFunctional,"integral",new=zero_only):
        fun=A.SUR.S.SurfaceFunctional(target,fields,0,0,rule)
        mass=A.E.core.upper(fun.mass_integral+fun.mass_rule)
    K,Tn=map(A.A,params)
    q0=Q(case["domain"]["q0_exact"])
    rho=Q(case["domain"]["rho_exact"])
    D=max(map(abs,A.SUR.rational_hull(Tn-K)))
    ea=A.A(D*A.SUR.S.binomial_tail(Q(-1,2),case["order"],rho))/(2*A.A(q0).sqrt())
    J=A.E.JET.surface_squared_jacobian_jet(fields["q_physical_R"],fields["B_physical_R"])
    grads=[A.SUR.poly_norm_upper(p,fields["Rlo"],fields["Rhi"]) for p in J["gradient"]]
    old=[mass*ea*A.A(x) for x in grads]
    new=m.surface_point_tails(A,target,rule,Y,0,case,params)
    assert len(old)==len(new)==6 and any(A.ends(v)[1]>0 for v in old)
    assert all(same(x,y) for x,y in zip(old,new,strict=True))
    for exact,rebuilt in zip(old,new,strict=True):
        assert m.exact_stored_radius(A,rebuilt)==fraction(arb(0,A.E.core.upper(exact)).rad())


def main():
    m=reread_module()
    A,v1=lookup_matrix_rollbacks(m)
    ctx.prec=640
    A.E.core.A.cache_clear()
    A.ad.A.cache_clear()
    deny=AssertionError("cold review attempted real plan/packet or force/H/cache processing")
    with patch.object(m,"compute",side_effect=deny),patch.object(m,"verify_terminal",side_effect=deny), \
         patch.object(A,"ordinary_bulk",side_effect=deny),patch.object(A,"surface_moments",side_effect=deny), \
         patch.object(A,"replay_component",side_effect=deny),patch.object(A.E.core,"local_terms",side_effect=deny), \
         patch.object(A.E.core,"assemble_point",side_effect=deny):
        safe_erosion(m,A)
        trace_literal_growth(m,A)
    assert ctx.prec==640 and not m.OUTPUT.exists()
    info={"status":"PASS_COLD_INDEPENDENT_EXACT_POINT_TAIL_SUPPLEMENT_REVIEW",
          "source_sha256":PINS["tighten_v2.py"],"protocol_sha256":PINS["protocol_v2.json"],"preflight_sha256":PINS["preflight_v2.json"],
          "audit_source_sha256":sha(SOURCE),"raw_replay_source_sha256":v1["source_sha256"],
          "raw_replay_math_sha256":v1["math_source_sha256"],"checks":{
           "original_and_companion_exact_sources_and_raw_transitive_premises_rehashed":True,
           "independent_outward_Arb_both_groupings_and_exact_radius_checked":True,
           "point_trace_zeroth_measure_and_original_rounded_radius_match_literal_640bit_formula":True,
           "exact_dyadic_inflation_and_fractional_unrounded_tail_erosion_synthetic":True,
           "width_guard_alone_does_not_justify_radius_and_noncontainment_still_detected":True,
           "source_proves_only_F_eroded_with_same_F_scope_and_U_unchanged":True,
           "no_actual_packets_or_mechanics":True},
          "actual_gate_scope":"Only after root verified genuinely completed raw v1 PASS/0 or HOLD/2 and all native component packets and source bindings are admitted. Run under the pinned unchanged library and 640 bit setup, separate actual readback required.",
          "scope":"The radius is read from precisely the reconstructed zero-centered insertion. Directed stored endpoints enclose the Minkowski sum of the truncated force/cache terms and that ball. Preserve producer, direct cache, H and original failed F evidence.",
          "caller_limit":"v2 internally checks pinned raw OS/terminal/receipt agreement. External admission should also require exact canonical raw 144+24 paths, scoped root OS identity, all typed maps and complete local H/force rows before accepting a stage. No old reducer can directly admit this phase."}
    out=SOURCE.with_name("receipt.json")
    with out.open("x",encoding="utf-8") as f:
        json.dump(info,f,indent=2)
        f.write("\n")
    print(sha(out))


if __name__=="__main__":
    main()
