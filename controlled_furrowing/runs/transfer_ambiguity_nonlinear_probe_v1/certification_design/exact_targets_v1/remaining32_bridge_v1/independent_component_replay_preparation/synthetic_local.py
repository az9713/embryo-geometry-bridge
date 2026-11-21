"""Exercise normalized signed assembly, tail inflation, and disjoint hold without data."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace as NS

from flint import arb,fmpq,ctx

import hull_math as hm

HERE=Path(__file__).resolve().parent
R=Q(1,2**34)
P,B="EXACT_CENTER_POINT_ONLY","ENTIRE_CLOSED_FIXED_PARAMETER_CUBE"


def A(v):
    v=Q(v)
    return arb(fmpq(v.numerator,v.denominator))


def ends(v):
    return tuple(Q(str(x.fmpq())) for x in (v.lower(),v.upper()))


def setup(which):
    body={"triangles":[{"J_bound":{"lower_exact":"7/8","upper_exact":"9/8"},"positive_J_for_entire_box":True}]}
    a={"active":[[1]],"surface_0_edges":[[0,1,2]]}
    # Production arrays have both two-index access and list conversion on slices.
    import numpy as np
    target=NS(a={k:np.asarray(x) for k,x in a.items()})
    proj={"global_indices":[4],"free_local":[0],"chart":{"kind":which}}
    current=NS(record_id="synthetic",target=target,drive=Q(3,7),
               base={"basal_area_modulus":Q(1),"basal_surface_tension":Q(0)})
    nodes=[[Q(0),Q(0)]]*6
    M=NS(T=None,E=NS(core=NS(upper=lambda x:x.upper())),POINT=P,BOX=B,RADIUS=R,arb=arb,ends=ends,A=A)
    selections={side:{"selector":{"rho_exact":"1/8"}} for side in (P,B)}
    d={"ctx":current,"rule":None,"point":nodes,"nodes":nodes,"body":body}
    i=(0,0)
    def ordinary(t,r,x,slot,s,strength,derivative):
        assert strength==Q(3,7)
        if derivative:
            return {i:A(1)},{i:Q(1,8)},proj,[0,1,2,3,4,5]
        return [A(-1)],[Q(1,8)],proj,[0,1,2,3,4,5]
    def scaled(t,r,x,slot,s,strength,scl):
        assert scl=="pure-scaled-premise"
        return ordinary(t,r,x,slot,s,strength,True)
    def surface(t,r,x,f,e,c,physics,derivative):
        assert physics==(Q(1),Q(0))
        if derivative:
            return {i:A(1)},{i:A(Q(1,8))},proj
        return [A(-1)],[A(Q(1,8))],proj
    def correction(maps,point,whole,J,physics,surface):
        assert len(maps)==(12 if surface else 64) and J==[0]
        return [A(2)],{i:A(Q(1,4))}
    M.ordinary_bulk,M.scaled_bulk_H,M.surface_moments=ordinary,scaled,surface
    M.corrected_cache=correction
    M.bulk_maps=lambda *args:range(64)
    M.surface_maps=lambda *args:range(12)
    return M,d,selections


def main():
    ctx.prec=640
    outcomes=[]
    for which in ("off_axis","axis_vertex","axis_edge_physical","surface"):
        M,d,p=setup(which)
        k="surface" if which=="surface" else "bulk"
        p|={"kind":k,"slot":0,"global_indices":[4]}
        raw={"kind":k,"slot":0,"global_indices":[4],"record_id":"synthetic","root_asserted":False,
             "point_force":{"scope":P,"force_tail_upper_exact":["1/16"],
                            "F_I_minus_F_cache_intervals_exact":[["1","9/8"]],
                            "Qstar_minus_cache_correction_exact":[["2","2"]]},
             "whole_box_hessian":{"scope":B,"radius_exact":str(R),"H_I_minus_H_cache_upper_exact":[[0,0,"1","3/2"]]}}
        a=hm.local(M,d,k,0,p,raw,(Q(1,8),Q(1,8)),"pure-scaled-premise")
        assert a["status"]=="PASS_INDEPENDENT_COMPONENT_WITH_EXACT_HULLS"
        assert a["checked_cache_sites"]==(12 if k=="surface" else 64)
        f=a["signed_channels"]["point_force"]["rows"][0]
        h=a["signed_channels"]["whole_box_H"]["rows"][0]
        assert hm.pair(f["independent_exact"])[0]<=Q(7,8) and hm.pair(f["independent_exact"])[1]>=Q(9,8)
        assert hm.pair(h["independent_exact"])[0]<=Q(9,8) and hm.pair(h["independent_exact"])[1]>=Q(11,8)
        assert f["producer_contains_independent"] is False and f["disjoint"] is False
        assert a["signed_channels"]["point_same_Y_cache_correction"]["rows"][0]["hull_exact"]==["2","2"]
        assert a["independent_point_and_whole_tails_exact"]["F"]==["1/8"]
        raw["point_force"]["F_I_minus_F_cache_intervals_exact"]=[["7","9"]]
        bad=hm.local(M,d,k,0,p,raw,(Q(1,8),Q(1,8)),"pure-scaled-premise")
        assert bad["status"]=="HOLD_DISJOINT_SIGNED_ENCLOSURES"
        assert bad["signed_channels"]["point_force"]["disjoint_indices"]==[0]
        outcomes.append(which)
    out=HERE/"synthetic_local_hull_readback.json"
    payload={"status":"PASS_SYNTHETIC_NO_MECHANICS_OWN_MOMENT_TAIL_AND_HULL_WRAPPER",
             "synthetic_routes":outcomes,"prior_authority_read":False,"actual_jobs":0,
             "local_hull_source_sha256":sha256(Path(hm.__file__).read_bytes()).hexdigest(),
             "fixture_sha256":sha256(Path(__file__).read_bytes()).hexdigest(),
             "checks":{"own_tail_can_exceed_producer_printed_tail_still_inside_own_cap":True,
                       "point_AND_H_and_direct_cache_signed_proof_preserved":True,
                       "truly_disjoint_producer_data_hold":True}}
    with out.open("x",encoding="utf8",newline="\n") as f:
        json.dump(payload,f,indent=2);f.write("\n")
    print(sha256(out.read_bytes()).hexdigest())


if __name__=="__main__":
    main()

