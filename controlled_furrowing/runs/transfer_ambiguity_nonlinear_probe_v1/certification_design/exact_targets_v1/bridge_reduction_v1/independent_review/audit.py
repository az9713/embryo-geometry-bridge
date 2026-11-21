"""Standalone interval-corner, scatter and typed-authority review."""
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
from itertools import combinations_with_replacement, product
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
FROZEN=BASE/"preparation/transfer_ambiguity_exact_bridge_reduction_frozen.py"
SOURCE_SHA="605c421d96f62dcc06ef35a898f16d3c79ff18cc693e057ff99d3aad71ab71de"
PROOF_SHA="438c404403fc82833a10f79f107a43b60969af3f03dee7fd4614de9451b208d0"
PREP_PINS={
    BASE/"preparation/protocol.json":"80ac4e7cc7faf99be684ca9f7036e534d8149f7d1c43765f418d0a65b2d389f8",
    BASE/"preparation/preflight.json":"af0016077db6b9dc9ed0d2327df93cc6a43d5de98d1e52b2267edec44d4acd2d",
}


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def load(p):
    spec=importlib.util.spec_from_file_location("pinned_bridge_reduction",p)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def bounds(strings):
    return tuple(map(Q,strings))


def contains(container,value):
    return container[0]<=value<=container[1]


def dot(a,b):
    return sum(x*y for x,y in zip(a,b,strict=True))


def packet(m, kind, slot, ix, fs, h):
    upper=[[i,j,*map(str,pair)] for (i,j),pair in
           zip(combinations_with_replacement(range(len(ix)),2),h,strict=True)]
    return {"kind":kind,"slot":slot,"global_indices":ix,
        "point_force":{"scope":m.POINT,"F_I_minus_F_cache_intervals_exact":[list(map(str,x)) for x in fs]},
        "whole_box_hessian":{"scope":m.BOX,"radius_exact":str(m.RADIUS),"H_I_minus_H_cache_upper_exact":upper}}


def fixture(m):
    p=packet(m,"bulk",0,[2,0,1],[(2,3),(-2,-1),(-1,5)],
             [(1,2),(-4,5),(1,2),(6,7),(-3,1),(-9,-2)])
    q=packet(m,"surface",1,[1,0],[(2,3),(2,5)],[(-1,4),(2,7),(-7,-6)])
    maps=[{k:x[k] for k in ("kind","slot","global_indices")} for x in (p,q)]
    F,H=m.assemble([p,q],maps,3)
    wantF=[(0,4),(1,8),(2,3)]
    wantH={(0,0):(-1,1),(0,1):(-1,8),(0,2):(-4,5),
           (1,1):(-10,2),(1,2):(1,2),(2,2):(1,2)}
    assert F==wantF and H==wantH
    R=[[Q(3,7),Q(-3,11),Q(2,5)], [Q(-1,3),Q(1),Q(-4,13)], [Q(0),Q(-1,2),Q(1,5)]]
    signed_anchor=[(Q(-3,40),Q(-1,50)),(Q(1,3),Q(7,20)),(Q(-1),Q(0))]
    zq=Q(1,8)
    result=m.contract(R,F,H,signed_anchor,zq)
    own_rows=list(map(bounds,result["signed_RF_I_rows_exact"]))
    defect_rows=list(map(bounds,result["signed_RF_defect_rows_exact"]))
    hrows=list(map(Q,result["RH_absolute_row_sum_bounds_exact"]))
    nforce=0
    for u in product(*wantF):
        for v in product(*signed_anchor):
            for i,r in enumerate(R):
                x=dot(r,u)
                assert contains(defect_rows[i],x)
                assert contains(own_rows[i],x+v[i])
            nforce+=1
    ids=tuple(wantH)
    nh=0
    for p0 in product(*(wantH[k] for k in ids)):
        h=dict(zip(ids,p0,strict=True))
        for i,row in enumerate(R):
            explicit=sum(abs(dot(row,[h[tuple(sorted((k,j)))] for k in range(3)])) for j in range(3))
            assert explicit<=hrows[i]
        nh+=1
    assert Q(result["b_upper_exact"])==max(hrows)
    assert Q(result["z_upper_exact"])==zq+max(hrows)

    assert m.classify([(-1,1)],0,1)["classification"]=="UNRESOLVED"
    assert m.classify([(-2,-2)],1,1)["classification"]=="UNRESOLVED"
    assert m.classify([(-3,-2)],1,Q(3,4))["classification"]=="SUFFICIENT_NO_ZERO_IN_THIS_CUBE"
    assert m.classify([(-3,1),(-2,-1)],0,2)["RF_I_norm_lower_exact"]=="1"
    root=m.classify([(-Q(1,16),Q(1,8))],Q(2,3),1)
    assert root["classification"]=="SUFFICIENT_UNIQUE_LOCAL_ZERO"
    assert Q(root["posterior_error_upper_exact"])==Q(3,8)

    def reject(f):
        try:
            f()
        except ValueError:
            return "REJECTED_VALUE_ERROR"
        except TypeError:
            return "UNCAUGHT_TYPE_ERROR"
        return "ACCEPTED"

    def hindex(item, value):
        other=deepcopy(p)
        other["whole_box_hessian"]["H_I_minus_H_cache_upper_exact"][0][item]=value
        return m.assemble([other,q],maps,3)

    bad_expected=deepcopy(maps)
    bad_expected[0]["global_indices"]=[2,False,1]
    issues={
        "boolean_Hessian_scope_index":reject(lambda:hindex(0,False)),
        "float_Hessian_scope_index":reject(lambda:hindex(1,0.0)),
        "boolean_component_slot":reject(lambda:m.assemble([p,q|{"slot":True}],maps,3)),
        "float_component_slot":reject(lambda:m.assemble([p,q|{"slot":1.0}],maps,3)),
        "boolean_expected_free_map":reject(lambda:m.assemble([p,q],bad_expected,3)),
    }
    for b in (True,False,0.0,1.0):
        assert reject(lambda:m.interval([b,1]))=="REJECTED_VALUE_ERROR"
        assert reject(lambda:m.contract([[b]],[(0,0)],{},[(0,0)],0))=="REJECTED_VALUE_ERROR"
    return {"force_and_anchor_corner_pairs":nforce,"symmetric_Hessian_corner_matrices":nh,
            "all_outward_norms_enclose_exact_corner_reference":True,
            "strict_exclusion_at_z_greater_or_equal_one_checked":True,
            "all_rational_interfaces_rejected_actual_float_bool":True,
            "noninteger_structural_cases":issues}


def main():
    p=Path(sys.argv[1]) if len(sys.argv)>1 else FROZEN
    expected=sys.argv[2] if len(sys.argv)>2 else SOURCE_SHA
    assert sha(p)==expected
    if expected==SOURCE_SHA:
        assert sha(BASE/"PROOF.md")==PROOF_SHA
        assert all(sha(f)==h for f,h in PREP_PINS.items())
    m=load(p)
    old=m.ctx.prec
    m.ctx.prec=640
    try:
        case=fixture(m)
    finally:
        m.ctx.prec=old
    issues=case["noninteger_structural_cases"]
    ok=all(value=="REJECTED_VALUE_ERROR" for value in issues.values())
    action="PASS_REDUCTION_KERNEL_MATH_AND_TYPED_INTERFACE" if ok else "HOLD_REDUCTION_TYPED_INTERFACE"
    output={"status":action, "source_sha256":expected,"proof_sha256":sha(BASE/"PROOF.md"),
            "audit_source_sha256":sha(Path(__file__)),"fixture":case,
            "findings": ([] if ok else [
                "Python equality lets a boolean and a float impersonate an integer component slot.",
                "Boolean upper-Hessian indices pass value-based coverage; float values pass coverage and later raise a TypeError.",
                "Externally supplied expected projection maps can compare equal even with a Boolean index.",
            ]),
            "proof_review":"Signed scatter and preconditioned interval sums are valid. Banach gives an interior unique local zero only under strict z and self-map bounds. The Neumann argument proves invertibility of R under z<1. The exclusion test only needs the same whole-cube bound and segment integration; it does not need R invertible or z<1. Each claim is conditional on separately authenticated same-target operands.",
            "external_conditions":"Packet scope radius and contract radius must be the same. Callers must authenticate binary R against the current cache and fix target provenance, all 168 slots, same signed cached anchor and derivative premise; the pure kernel does not establish that authority.",
            "actual_packets_read":0,"physics_jobs":0}
    dest=HERE/("receipt_v1.json" if expected==SOURCE_SHA else "receipt_v2.json")
    with dest.open("x",encoding="utf-8") as f:
        json.dump(output,f,indent=2)
        f.write("\n")
    print(action,sha(dest),issues)


if __name__=="__main__":
    main()
