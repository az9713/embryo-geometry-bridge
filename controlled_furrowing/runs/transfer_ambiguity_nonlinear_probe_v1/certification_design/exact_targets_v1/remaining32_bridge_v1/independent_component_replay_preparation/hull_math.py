"""Local independent bound plus exact hull, using reviewed explicit-context primitives."""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import ceil


def require(ok,message):
    if not ok:
        raise ValueError(message)


def pair(value):
    require(type(value) in (list,tuple) and len(value)==2 and
            all(type(x) in (int,str,Q) for x in value),"two exact signed endpoints")
    lo,hi=map(Q,value)
    require(lo<=hi,"reversed signed interval")
    return lo,hi


def hull_and_compare(independent,producer):
    a,b=pair(independent)
    c,d=pair(producer)
    return (min(a,c),max(b,d)), c<=a<=b<=d, max(a,c)>min(b,d)


def packed(items):
    return list(map(str,items))


def labelled_channel(own,raw):
    require(own.keys()==raw.keys(),"signed entries must match exactly")
    rows=[]
    nonnested=[]
    disjoint=[]
    for key in own:
        inner,outer=own[key],raw[key]
        combined,nested,separate=hull_and_compare(inner,outer)
        tag=list(key) if type(key) is tuple else key
        rows.append({"index":tag,"independent_exact":packed(inner),"producer_exact":packed(outer),
                     "hull_exact":packed(combined),"producer_contains_independent":nested,"disjoint":separate})
        if not nested:
            nonnested.append(tag)
        if separate:
            disjoint.append(tag)
    return {"rows":rows,"producer_containment_misses":nonnested,"disjoint_indices":disjoint}


def local(M,d,kind,slot,selector,raw,cap,scale):
    ctx,rule,point,nodes=(d[k] for k in ("ctx","rule","point","nodes"))
    T=M.T
    P,B=M.POINT,M.BOX
    C=M.E.core
    require(ctx.record_id==raw["record_id"] and slot==selector["slot"]==raw["slot"]
            and kind==selector["kind"]==raw["kind"] and type(raw["slot"]) is type(selector["slot"]) is int
            and raw["point_force"]["scope"]==P and raw["whole_box_hessian"]["scope"]==B
            and Q(raw["whole_box_hessian"]["radius_exact"])==M.RADIUS and raw["root_asserted"] is False,
            "separate exact point F and whole matched closed-cube H required")
    target=ctx.target
    if kind=="bulk":
        sp=selector[P]["selector"]
        sb=selector[B]["selector"]
        j=d["body"]["triangles"][slot]
        lo,hi=map(Q,(j["J_bound"]["lower_exact"],j["J_bound"]["upper_exact"]))
        rho=Q(ceil(max(abs(lo-1),abs(hi-1))*2**48),2**48)
        require(j["positive_J_for_entire_box"] is True and Q(1,10**6)<lo<=hi
                and Q(sp["rho_exact"])==Q(sb["rho_exact"])==rho<1,
                "independent analytic volume domain differs")
        drive=ctx.drive*int(target.a["active"][slot,0])
        F,Ft,fp,ids=M.ordinary_bulk(target,rule,point,slot,sp,drive,False)
        if fp["chart"]["kind"]=="off_axis":
            H,Ht,fh,trace=M.scaled_bulk_H(target,rule,nodes,slot,sb,drive,scale)
        else:
            H,Ht,fh,trace=M.ordinary_bulk(target,rule,nodes,slot,sb,drive,True)
        require(ids==trace,"same reference triangle required")
        maps=list(M.bulk_maps(target,rule,slot,ids,fp["free_local"]))
        physics=(drive,)
        count=64
    else:
        require(kind=="surface","typed known component only")
        face,i=divmod(slot,12)
        k=("basal","apical")[face]
        physics=(ctx.base[k+"_area_modulus"],ctx.base[k+"_surface_tension"])
        F,Ft,fp=M.surface_moments(target,rule,point,face,i,selector[P],physics,False)
        H,Ht,fh=M.surface_moments(target,rule,nodes,face,i,selector[B],physics,True)
        ids=target.a[f"surface_{face}_edges"][i].tolist()
        maps=list(M.surface_maps(target,rule,face,i,ids))
        count=12
    J=fp["free_local"]
    ids_free=fp["global_indices"]
    require(J==fh["free_local"] and ids_free==fh["global_indices"]==selector["global_indices"]==raw["global_indices"]
            and len(maps)==count and all(type(x) is int for x in ids_free),"same complete free-coordinate projection and cache sites")
    ci,ch=M.corrected_cache(maps,[point[i][j] for i in ids for j in (0,1)],
                          [nodes[i][j] for i in ids for j in (0,1)],J,physics,kind=="surface")
    def upper(t):
        return t if type(t) is Q else M.ends(t)[1]
    ft=list(map(upper,Ft))
    ht={k:upper(t) for k,t in Ht.items()}
    order=list(combinations_with_replacement(range(len(J)),2))
    require(len(ft)==len(J) and list(H)==list(Ht)==order and list(ch)==order and
            all(0<=v<=cap[0] for v in ft) and all(0<=v<=cap[1] for v in ht.values()),
            "independently computed sufficient two-mass tails exceed reviewed caps")
    f_raw=raw["point_force"]
    printed=list(map(Q,f_raw["force_tail_upper_exact"]))
    require(len(printed)==len(J) and all(0<=t<=cap[0] for t in printed),"producer tail outside policy, separately checked")
    # Inflate each locally valid independent truncated expression by a sufficient
    # analytical tail. This step never uses a producer endpoint as a premise.
    own_F=[x+y+M.arb(0,C.upper(M.A(t))) for x,y,t in zip(F,ci,ft,strict=True)]
    own_H={k:H[k]+ch[k]+M.arb(0,C.upper(M.A(ht[k]))) for k in order}
    read_H=raw["whole_box_hessian"]["H_I_minus_H_cache_upper_exact"]
    require(type(read_H) is list and len(read_H)==len(order)
            and all(type(row) is list and len(row)==4 and type(row[0]) is type(row[1]) is int for row in read_H)
            and [tuple(row[:2]) for row in read_H]==order,"complete typed symmetric upper Hessian required")
    rawH={(i,j):pair((lo,hi)) for i,j,lo,hi in read_H}
    def enumerate_ends(values):
        return {i:M.ends(v) for i,v in enumerate(values)}
    rawF={i:pair(v) for i,v in enumerate(f_raw["F_I_minus_F_cache_intervals_exact"])}
    rawC={i:pair(v) for i,v in enumerate(f_raw["Qstar_minus_cache_correction_exact"])}
    F_proof=labelled_channel(enumerate_ends(own_F),rawF)
    H_proof=labelled_channel({k:M.ends(v) for k,v in own_H.items()},rawH)
    CF_proof=labelled_channel(enumerate_ends(ci),rawC)
    holds=any(ch["disjoint_indices"] for ch in (F_proof,H_proof,CF_proof))
    return {"kind":kind,"slot":slot,"record_id":ctx.record_id,"global_indices":ids_free,
            "status":"HOLD_DISJOINT_SIGNED_ENCLOSURES" if holds else "PASS_INDEPENDENT_COMPONENT_WITH_EXACT_HULLS",
            "checked_cache_sites":count,"point_force_scope":P,"whole_H_scope":B,"radius_exact":str(M.RADIUS),
            "independent_point_and_whole_tails_exact":{"F":packed(ft),"H":[[i,j,str(ht[i,j])] for i,j in order]},
            "independent_same_Y_whole_H_correction_exact":[[i,j,*packed(M.ends(ch[i,j]))] for i,j in order],
            "producer_force_tails_exact":packed(printed),
            "signed_channels":{"point_force":F_proof,"whole_box_H":H_proof,"point_same_Y_cache_correction":CF_proof},
            "only_independently_valid_enclosures_used_for_hull_truth":True,
            "overlap_not_used_as_proof":True,"root_asserted":False}

