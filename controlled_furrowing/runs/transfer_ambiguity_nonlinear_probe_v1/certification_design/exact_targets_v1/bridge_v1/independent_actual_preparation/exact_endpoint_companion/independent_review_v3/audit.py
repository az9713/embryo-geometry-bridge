"""Cold source and arithmetic readback for endpoint companion v3."""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
from itertools import product
import json
from math import factorial
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch

from flint import ctx

SOURCE=Path(__file__).resolve()
BASE=SOURCE.parents[1]
PINS={
    "tighten_v3.py":"d9ca0122c62a6bd5f94e8bdba078c715adca9056d8b8822e693e2d4a866419de",
    "protocol_v3.json":"a825830947d53daee36a9f75904a85d2097375fadcf070fcd3ee30d1ffb7c9d7",
    "preflight_v3.json":"a62668937b49f9e76e65706d70f2ba015820923150183396815ff26c82a3c01d",
}
V2="ef0f7cc5a0f0649c214d26a7c2f00c9a60b34a9cd7a1e0adc634fbafb6a44ea5"

def sha(p):
    return sha256(p.read_bytes()).hexdigest()

def coeff(p,k):
    v=Q(1)
    for j in range(k):
        v*=p-j
    return v/factorial(k)

def arithmetic(v,A):
    ctx.prec=640
    A.E.core.A.cache_clear()
    A.ad.A.cache_clear()
    for n,rho in product(range(2,65),(Q(0),Q(1,17),Q(3,8),Q(7,8))):
        m=n+1
        actual=v.scalar_tail_sharp(A,n,rho)
        x,y=[abs(coeff(Q(-1,2),k))*rho**k for k in (m,m+1)]
        assert actual==x+y/(1-rho)
        old=A.SUR.S.binomial_tail(Q(-1,2),n,rho)
        assert actual<=old and (rho==0 or actual<old)
    # For -3/2 the weighted ratio exceeds rho.  Source must call the old
    # reviewed -3/2 majorant for H; the point-only sharp rule cannot be reused.
    rho=Q(3,4)
    assert abs(coeff(Q(-3,2),4)/coeff(Q(-3,2),3))*rho>rho
    q0,K,T=Q(4),Q(5),Q(9)
    a,b=A.SUR.uniform_series(A.T.X,q0,64,K,T)
    for n in range(65):
        f=(T-K)*coeff(Q(-1,2),n)/4+(K/2 if n==0 else 0)
        h=-(T-K)*coeff(Q(-3,2),n)/32
        assert A.ends(a[(n,0)])==(f,f) and A.ends(b[(n,0)])==(h,h)
    src=(BASE/"tighten_v3.py").read_text(encoding="utf8")
    assert 'eb=A.A(D*A.SUR.S.binomial_tail(Q(-3,2),order,rho))' in src
    assert 's=exact_stored_radius(A,lift)' in src
    assert 'local=A.SUR.L.consistency_bound(a,b,Y,selected,K,Tn)["H_difference"]' in src
    assert 'cache={k:v+A.A(w)*local[k[0],k[1]]' in src
    # Widen with one Arb group and erode using only the exact inserted
    # radius reconstructed from the old tail, even with a different t.
    old,sharp=Q(1,3),Q(1,5)
    x=A.A(Q(2,17))
    old_ball=A.arb(0,A.E.core.upper(A.A(old)))
    base=A.ends((x+old_ball)+A.arb(0))
    tight,s,t,t0=v.tighten(A,base,old,Q(1),A.A(sharp))
    assert s==Q(str(old_ball.rad().fmpq())) and t0==old
    assert sharp<=t<old and tight==(base[0]+s-t,base[1]-s+t)
    assert A.contain(tight,(Q(2,17)-sharp,Q(2,17)+sharp))

def test_exit(v):
    pilot="mock"
    raw=v.RAW/pilot
    term={"status":"FINISHED_WITH_HOLD","exit_code":2,"pid":4,"receipt_sha256":"r",
          **{k:v.PINS[v.PREP/p] for k,p in (("source_sha256","replay.py"),("protocol_sha256","protocol.json"),
                                          ("preflight_sha256","preflight.json"))}}
    receipt={**term,"status":"HOLD_INDEPENDENT_COMPONENT_CONTAINMENT","pilot_id":pilot,"components":[{}]*168,
             "checks":{"true_os_complete_producer_and_both_child_terminals":True,
                       "same_exact_point_and_whole_box_all_168_bindings":True},
             "checked_packet_sha256":{},"plan_packet_sha256":{}}
    obs={"status":"OBSERVED_RAW_BRIDGE_REPLAY_SHELL_TERMINAL","terminal_status":term["status"],
         "observed_exit_code":1,"native_documented_exit_code":2,"pid":4,"unified_exec_session":5,
         "terminal_sha256":"t","receipt_sha256":"r"}
    o=BASE/"mock_os.json"
    args=SimpleNamespace(pilot_id=pilot,raw_terminal_sha256="t",raw_receipt_sha256="r",
                         raw_os_exit=o,raw_os_exit_sha256="o")
    def trial(change):
        with patch.object(v,"read",side_effect=lambda p,h: {o:change,raw/"terminal.json":term,
                                                           raw/"receipt.json":receipt}[p]):
            return v.verify_terminal(args)
    assert trial(obs) is receipt
    for change in ({"observed_exit_code":0},{"observed_exit_code":True},{"native_documented_exit_code":1},{"pid":11}):
        try:
            trial({**obs,**change})
        except AssertionError:
            continue
        raise AssertionError("invalid native/host pairing accepted")


def main():
    for p,h in PINS.items():
        assert sha(BASE/p)==h
    assert sha(BASE/"tail_supplement_v2/independent_review/receipt.json")==V2
    spec=importlib.util.spec_from_file_location("_cold_endpoint_companion_v3",BASE/"tighten_v3.py")
    v=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(v)
    A,EX,DX=v.environment()
    assert A.E is EX.E and DX.EX is EX and Path(A.__file__).resolve()==v.PREP/"replay_math.py"
    assert v.phase()==json.loads((BASE/"protocol_v3.json").read_bytes())
    raw=json.loads((v.PREP/"protocol.json").read_bytes())
    assert raw["precision_bits"]==640 and raw["math_source_sha256"]==sha(Path(A.__file__))
    for p,h in raw["source_and_premise_hashes"].items():
        assert sha(v.ROOT/p)==h
    stop=AssertionError("actual plan, cache sum or moment entered the cold audit")
    with patch.object(v,"compute",side_effect=stop),patch.object(v,"fresh_surface_H",side_effect=stop), \
            patch.object(A,"ordinary_bulk",side_effect=stop),patch.object(A,"surface_moments",side_effect=stop), \
            patch.object(A.E.core,"local_terms",side_effect=stop),patch.object(EX,"open_context",side_effect=stop):
        arithmetic(v,A)
        test_exit(v)
    assert not v.OUTPUT.exists()
    record={"status":"PASS_INDEPENDENT_COLD_SURFACE_H_ENDPOINT_COMPANION",
            "source_sha256":PINS["tighten_v3.py"],"protocol_sha256":PINS["protocol_v3.json"],
            "preflight_sha256":PINS["preflight_v3.json"],"audit_source_sha256":sha(SOURCE),
            "prior_F_only_erosion_review_sha256":V2,
            "checks":{"raw_and_companion_source_transitive_inputs_current":True,
                      "first_omitted_term_plus_geometric_tightening_only_for_negative_half":True,
                      "old_negative_three_halves_tail_kept_for_H_and_canonical_series_coefficients_exact":True,
                      "exact_constructed_old_radius_erosion_and_same_point_force_scope":True,
                      "fresh_24_surface_H_signed_independent_maps_and_single_positive_weight":True,
                      "native_declared_HOLD2_and_observed_host1_separately_checked":True,
                      "no_actual_packet_cache_sum_or_moment_jobs":True},
            "scope":"Cold only. Actual source requires root-observed complete raw v1, original direct cache and bulk H containment and no previous H reuse on any of the 24 surface edges. Independent actual canonical packet/radius/F/H containment and OS readback remain necessary.",
            "shared_premises":"Reviewed canonical SUR.uniform_series, jets, scalar tails, moments and local third-tensor kernel. Replay independently reconstructs exact surface maps and signed weighted whole-box H corrections. Producer intervals remain unchanged."}
    out=SOURCE.with_name("receipt.json")
    h=v.save(out,record)
    print(h)

if __name__=="__main__":
    main()
