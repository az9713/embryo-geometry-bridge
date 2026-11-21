"""Pure local bridge replay using reviewed chart, moment and third-jet primitives."""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import ceil
from pathlib import Path
from time import perf_counter

from flint import arb

import exact_cohort_integrator as E
import review_algebraic_root as ad
from review_centered_cohort_packets import bulk_gradient, line_J, physical_derivatives, p2

T, U, SUR, A = E.T, E.U, E.SUR, E.U.A
POINT = "EXACT_CENTER_POINT_ONLY"
BOX = "ENTIRE_CLOSED_FIXED_PARAMETER_CUBE"
RADIUS = Q(1, 2**34)
PAIRS = ((0,1),(1,2),(2,0))


def require(x, msg):
    if not x:
        raise ValueError(msg)


def q(x):
    require(type(x) in (int, str, Q), "explicit exact scalar")
    return Q(x)


def ends(ball):
    require(ball.is_finite(), "finite Arb enclosure")
    return Q(str(ball.lower().fmpq())), Q(str(ball.upper().fmpq()))


def interval(row):
    require(type(row) in (list,tuple) and len(row)==2, "signed endpoint pair")
    lo,hi=map(q,row)
    require(lo<=hi, "ordered exact interval")
    return lo,hi


def contain(outer,inner):
    return outer[0]<=inner[0]<=inner[1]<=outer[1]


def dot(a,b):
    return sum((s*t for s,t in zip(a,b,strict=True)),Q(0))


def norm(p):
    return sum(max(map(abs,SUR.rational_hull(v))) for v in p.values())


def expand(poly, coefs):
    # Same ascending expansion of powers. The reviewed moment functional
    # remains the only source of transcendental interval moments.
    curr=T.ONE
    out={}
    for i,c in enumerate(coefs):
        out=T.add(out,T.scale(curr,c))
        if i+1<len(coefs):
            curr=T.mul(curr,poly)
    return out


def apply_diagonal(left,right,diagonal,apply):
    p={}
    for a,b,s in zip(left,right,diagonal,strict=True):
        p=T.add(p,T.scale(T.mul(a,b),s))
    return apply(p)


def project_direct(star,binary,y,projection,params,surface):
    qs=[dot(row,y) for row in star]
    qb=[dot(row,y) for row in binary]
    if surface:
        K,Tn=params
        gs=ad.surface_AD(list(map(ad.A,qs)),ad.A(K),ad.A(Tn),ad.A(1))[0]
        gb=ad.surface_AD(list(map(ad.A,qb)),ad.A(K),ad.A(Tn),ad.A(1))[0]
    else:
        line_J(qs,qb)
        drive,=params
        gs=bulk_gradient(qs,drive)
        gb=bulk_gradient(qb,drive)
    return [sum((A(row[j])*h for row,h in zip(star,gs,strict=True)),arb(0))-
            sum((A(row[j])*h for row,h in zip(binary,gb,strict=True)),arb(0))
            for j in projection]


def bulk_maps(target,rule,e,nodes,local):
    a=target.a
    corners=[tuple(T.rational(x) for x in a["reference"][j]) for j in nodes[:3]]
    _,dR,dZ=physical_derivatives(corners)
    for k,(lam,w) in enumerate(zip(rule["volume_barycentric"],rule["volume_weights"][e],strict=True)):
        oldN=[T.rational(x) for x in a["shape"][k]]
        oldD=[tuple(map(T.rational,v)) for v in a["shape_grad"][e,k]]
        xi=oldN[1]+(oldN[3]+oldN[4])/2
        eta=oldN[2]+(oldN[4]+oldN[5])/2
        R=dot(lam,[p[0] for p in corners])
        Rc=T.rational(a["qr"][e,k])
        require(tuple(lam)==(1-xi-eta,xi,eta) and min(lam)>0 and R>0 and Rc>0, "positive exact source sites")
        require(w==T.rational(a["volume_weights"][e,k])/target.energy_scale>0,"positive exact bulk weight")
        newN=[x*(2*x-1) for x in lam]+[4*lam[i]*lam[j] for i,j in PAIRS]
        Bn=[list(z) for z in zip(p2(lam,dR),p2(lam,dZ),strict=True)]
        def full(shape,grad,r):
            M=[[Q(0)]*12 for _ in range(5)]
            for p,(d_r,d_z) in enumerate(grad):
                M[0][2*p]=M[2][2*p+1]=d_r
                M[1][2*p]=M[3][2*p+1]=d_z
                M[4][2*p]=shape[p]/r
            return M
        yield full(newN,Bn,R),full(oldN,oldD,Rc),w


def surface_maps(target,rule,f,e,nodes):
    a=target.a
    lo,hi=[T.rational(a["reference"][j,0]) for j in (nodes[0],nodes[-1])]
    L=hi-lo
    require(L>0 and lo>=0,"positive edge material length")
    for k,(s,w) in enumerate(zip(rule["surface_s"][f],rule["surface_weights"][f][e],strict=True)):
        oldN=list(map(T.rational,a[f"surface_{f}_N"][k]))
        oldD=list(map(T.rational,a[f"surface_{f}_dN"][k]))
        r_cache=T.rational(a[f"surface_{f}_r0"][e,k])
        l_cache=T.rational(a[f"surface_{f}_l0"][e])
        require(0<s<1 and s==oldN[2]+oldN[1]/2 and r_cache>0 and l_cache>0,"positive exact surface sites")
        require(w==T.rational(a[f"surface_{f}_weights"][e,k])/target.energy_scale>0,"positive exact edge weight")
        R=lo+s*L
        N=[(1-s)*(1-2*s),4*s*(1-s),s*(2*s-1)]
        D=[4*s-3,4-8*s,4*s-1]
        def full(shape,deriv,hoop):
            M=[[Q(0)]*6 for _ in range(3)]
            for j in range(3):
                M[0][2*j]=hoop*shape[j]
                M[1][2*j]=M[2][2*j+1]=deriv[j]/L
            return M
        yield full(N,D,1/R),full(oldN,oldD,L/(r_cache*l_cache)),w


def corrected_cache(maps, point, box, sel, params, surface):
    n=len(sel)
    keys=tuple(combinations_with_replacement(range(n),2))
    f=[arb(0) for _ in sel]
    H={k:arb(0) for k in keys}
    for star,cached,weight in maps:
        diff=project_direct(star,cached,point,sel,params,surface)
        if surface:
            Hs=SUR.L.consistency_bound(star,cached,box,sel,*map(A,params))["H_difference"]
        else:
            drive,=params
            Hs=U.C.consistency_bound(star,cached,box,sel,1,5000,A(drive).union(A(drive)))["difference"]
        w=A(weight)
        f=[x+w*y for x,y in zip(f,diff,strict=True)]
        H={k:a+w*Hs[k[0],k[1]] for k,a in H.items()}
    return f,H


def cache_clears(fun):
    for name in ("integral","node_power","site_power","defect"):
        method=getattr(fun,name,None)
        if hasattr(method,"cache_clear"):
            method.cache_clear()
    E.M._triangle_moment.cache_clear()


def ordinary_bulk(target,rule,node_box,i,selector,requested,derivative):
    fields,ids=E._bulk_fields(target,node_box,i)
    local=fields["free_local"]
    require(local==selector["free_local"] and fields["global_indices"]==selector["global_indices"],"bulk field map and selector")
    drive=A(requested).union(A(requested)) if derivative else A(requested)
    diagonal=[A(1),A(1)+drive,A(1),A(1)+drive,A(1)]
    jet=E.JET.bulk_jacobian_jet(fields["q"],fields["B"])
    grad,H=jet["gradient"],jet["hessian_upper"]
    N=selector["degree"]
    require(type(N) is int and 2<=N<=32, "explicit bulk order")
    rho=q(selector["rho_exact"])
    a,b=E.M.volumetric_derivative_tail_bounds(N,rho,1,5000)
    require((a,b)==(q(selector["phi_prime_tail_upper_exact"]),q(selector["phi_second_tail_upper_exact"])),
            "same reviewed bulk scalar tails")
    mass=E._geometry_mass_upper(target,fields,i,rule)
    require(mass==q(selector["geometry_two_mass_upper_exact"]), "same two bulk measures")
    mg=[E._norm(x,fields) for x in grad]
    series=E.M.volumetric_derivative_coefficients(N,1,5000)
    polynomial={k:A(v) for k,v in T.sub(jet["value"],T.ONE).items()}
    f1=expand(polynomial,list(map(A,series[0])))
    columns=list(zip(*fields["B"],strict=True))
    fun=E._functional_type(fields)(target,fields,i,rule)
    try:
        if derivative:
            f2=expand(polynomial,list(map(A,series[1])))
            parts=E.JET.contract_volumetric_defect(jet,f1,f2,fun.defect)["hessian_upper"]
            main={k:v+apply_diagonal(columns[k[0]],columns[k[1]],diagonal,fun.apply) for k,v in parts.items()}
            tails={k:mass*(b*mg[k[0]]*mg[k[1]]+a*E._norm(H[k],fields)) for k in H}
        else:
            d=E.JET.moment_dual(f1,{k for p in grad for k in p},fun.defect)
            main=[E.JET.contract_factor(p,d)+apply_diagonal(fields["q"],col,diagonal,fun.apply)
                  for p,col in zip(grad,columns,strict=True)]
            tails=[mass*a*g for g in mg]
    finally:
        cache_clears(fun)
    return main,tails,fields,ids


def scaled_bulk_H(target,rule,node_box,i,selector,requested,scaled):
    fields,ids=E._bulk_fields(target,node_box,i)
    geo=scaled.chart(fields["corners"])
    local=fields["free_local"]
    require(local==selector["free_local"] and fields["global_indices"]==selector["global_indices"],"scaled field map")
    fun=scaled.Functional(geo,target.energy_scale,rule["volume_barycentric"],rule["volume_weights"][i])
    try:
        q0,B=scaled.numerators(geo,[node_box[j] for j in ids],local)
        jet=E.JET.bulk_jacobian_jet(q0,B)
        grad,H=jet["gradient"],jet["hessian_upper"]
        degree=selector["degree"]
        require(type(degree) is int and 2<=degree<=32,"explicit scaled bulk order")
        A0,B0=scaled.series_numerators(T.sub(jet["value"],geo.Rhat),geo.Rhat,degree)
        pp={k:T.mul(grad[k[0]],grad[k[1]]) for k in H}
        # Include the gradient support in the first dual to match the
        # reviewed functional. Only the Hessian is extracted on this box.
        aa=E.JET.moment_dual(A0,{k for p in [*grad,*H.values()] for k in p},
                             lambda a,b:fun.defect(a,b,degree+1))
        bb=E.JET.moment_dual(B0,{k for p in pp.values() for k in p},
                             lambda a,b:fun.defect(a,b,degree+2))
        columns=list(zip(*B,strict=True))
        dr=A(requested).union(A(requested))
        diag=[A(1),A(1)+dr,A(1),A(1)+dr]
        main={}
        for k,p in H.items():
            i0,j0=k
            s=apply_diagonal(columns[i0][:4],columns[j0][:4],diag,lambda f:fun.apply(f,0))
            s+=fun.apply(T.mul(columns[i0][4],columns[j0][4]),2)
            main[k]=E.JET.contract_factor(p,aa)+E.JET.contract_factor(pp[k],bb)+s
        mass=E._geometry_mass_upper(target,fields,i,rule)
        require(mass==q(selector["geometry_two_mass_upper_exact"]),"same scaled two masses")
        d1,d2=E.M.volumetric_derivative_tail_bounds(degree,q(selector["rho_exact"]),1,5000)
        require((d1,d2)==(q(selector["phi_prime_tail_upper_exact"]),q(selector["phi_second_tail_upper_exact"])),
                "same scaled scalar tails")
        g0,h0=map(q,(selector["max_DJ_norm_upper_exact"],selector["max_D2J_norm_upper_exact"]))
        mg=[min(g0,norm(p)/geo.Rhat_min) for p in grad]
        mh={k:min(h0,norm(p)/geo.Rhat_min) for k,p in H.items()}
        tails={k:mass*(d2*mg[k[0]]*mg[k[1]]+d1*v) for k,v in mh.items()}
    finally:
        cache_clears(fun)
    return main,tails,fields,ids


def surface_moments(target,rule,nodes,f,e,case,params,derivative):
    K,Tn=map(A,params)
    fields=target.surface(f,e,nodes=nodes)
    require(fields["global_indices"]==case["global_indices"],"whole trace projection agrees with plan")
    domain=case["domain"]
    SUR.verify_frozen_domain(fields,domain)
    lo,hi=interval(domain["G_bounds_exact"])
    rho0=(hi-lo)/(hi+lo)
    q0=(lo+hi)/2
    rho=Q(ceil(rho0*2**48),2**48)
    require(lo>0 and rho<1 and (rho0,rho,q0)==tuple(map(q,(domain["raw_rho_exact"],domain["rho_exact"],domain["q0_exact"]))),
            "strict same closed surface G and outward rho")
    order=case["order"]
    require(type(order) is int and 2<=order<=64,"explicit surface order")
    J=E.JET.surface_squared_jacobian_jet(fields["q_physical_R"],fields["B_physical_R"])
    pol=T.scale(T.sub(J["value"],T.constant(q0)),1/q0)
    root=A(q0).sqrt()
    A0=expand(pol,[(Tn-K)/(2*root)*A(v) for v in SUR.S.binomial_coefficients(Q(-1,2),order)])
    A0=T.add(T.constant(K/2),A0)
    B0=expand(pol,[-(Tn-K)/(4*A(q0)*root)*A(v) for v in SUR.S.binomial_coefficients(Q(-3,2),order)])
    fun=SUR.S.SurfaceFunctional(target,fields,f,e,rule)
    try:
        moments=E.JET.contract_volumetric_defect(J,A0,B0,fun.defect)
        mass=E.core.upper(fun.mass_integral+fun.mass_rule)
        D=max(map(abs,SUR.rational_hull(Tn-K)))
        ea=A(D*SUR.S.binomial_tail(Q(-1,2),order,rho))/(2*root)
        eb=A(D*SUR.S.binomial_tail(Q(-3,2),order,rho))/(4*A(q0)*root)
        g=[SUR.poly_norm_upper(p,fields["Rlo"],fields["Rhi"]) for p in J["gradient"]]
        if derivative:
            tails={k:mass*(ea*A(SUR.poly_norm_upper(p,fields["Rlo"],fields["Rhi"]))+eb*A(g[k[0]]*g[k[1]]))
                   for k,p in J["hessian_upper"].items()}
            return moments["hessian_upper"],tails,fields
        return moments["gradient"],[mass*ea*A(v) for v in g],fields
    finally:
        cache_clears(fun)


def replay_component(ctx,rule,point,nodes,body,i,kind,plan,packet,tail_caps,scaled):
    started=perf_counter()
    require(type(i) is int and type(plan["slot"]) is int and type(packet["slot"]) is int
            and i==plan["slot"]==packet["slot"] and kind==plan["kind"]==packet["kind"],"component identity")
    require(packet["point_force"]["scope"]==POINT and packet["whole_box_hessian"]["scope"]==BOX
            and q(packet["whole_box_hessian"]["radius_exact"])==RADIUS,"point/whole derivative scopes")
    require(packet["root_asserted"] is False,"component evidence only")
    target=ctx.target
    if kind=="bulk":
        sp=plan[POINT]["selector"]
        sb=plan[BOX]["selector"]
        require(Q(sp["rho_exact"])==Q(sb["rho_exact"]),"same positive J domain")
        j=body["triangles"][i]
        lo,hi=map(q,(j["J_bound"]["lower_exact"],j["J_bound"]["upper_exact"]))
        rho0=max(abs(lo-1),abs(hi-1))
        require(j["positive_J_for_entire_box"] is True and lo>Q(1,10**6)
                and q(sp["rho_exact"])==Q(ceil(rho0*2**48),2**48)<1,"matching positive bulk domain")
        drive=ctx.drive*int(target.a["active"][i,0])
        F,Ft,fp,ids=ordinary_bulk(target,rule,point,i,sp,drive,False)
        if fp["chart"]["kind"]=="off_axis":
            H,Ht,fh,wids=scaled_bulk_H(target,rule,nodes,i,sb,drive,scaled)
        else:
            H,Ht,fh,wids=ordinary_bulk(target,rule,nodes,i,sb,drive,True)
        require(ids==wids,"point and whole bulk triangle match")
        maps=list(bulk_maps(target,rule,i,ids,fp["free_local"]))
        require(len(maps)==64,"all point/whole bulk sites")
        param=(drive,)
    else:
        require(kind=="surface","known component kind")
        f,e=divmod(i,12)
        k=("basal","apical")[f]
        param=(ctx.base[k+"_area_modulus"],ctx.base[k+"_surface_tension"])
        F,Ft,fp=surface_moments(target,rule,point,f,e,plan[POINT],param,False)
        H,Ht,fh=surface_moments(target,rule,nodes,f,e,plan[BOX],param,True)
        ids=target.a[f"surface_{f}_edges"][e].tolist()
        maps=list(surface_maps(target,rule,f,e,ids))
        require(len(maps)==12,"all point/whole surface sites")
    sel=fp["free_local"]
    glob=fp["global_indices"]
    require(glob==fh["global_indices"]==plan["global_indices"]==packet["global_indices"]
            and sel==fh["free_local"] and all(type(i) is int for i in packet["global_indices"]),"same actual free projection")
    yp=[point[i0][k] for i0 in ids for k in (0,1)]
    yw=[nodes[i0][k] for i0 in ids for k in (0,1)]
    cf,ch=corrected_cache(maps,yp,yw,sel,param,kind=="surface")
    f=packet["point_force"]
    used_tails=list(map(q,f["force_tail_upper_exact"]))
    pft=[t if isinstance(t,Q) else ends(t)[1] for t in Ft]
    pht={k:t if isinstance(t,Q) else ends(t)[1] for k,t in Ht.items()}
    require(len(used_tails)==len(pft)==len(sel) and
            all(0<=a<=b<=tail_caps[0] for a,b in zip(pft,used_tails,strict=True))
            and all(0<=t<=tail_caps[1] for t in pht.values()),"same independently bounded point and whole tails")
    if kind=="bulk":
        force=[a+b+arb(0,E.core.upper(A(t))) for a,b,t in zip(F,cf,Ft,strict=True)]
        hess={k:H[k]+ch[k]+arb(0,E.core.upper(A(Ht[k]))) for k in H}
    else:
        force=[a+arb(0,E.core.upper(t))+b for a,b,t in zip(F,cf,Ft,strict=True)]
        hess={k:H[k]+arb(0,E.core.upper(Ht[k]))+ch[k] for k in H}
    exact_f=list(map(ends,force))
    exact_h={k:ends(v) for k,v in hess.items()}
    supplied_f=list(map(interval,f["F_I_minus_F_cache_intervals_exact"]))
    supplied_c=list(map(interval,f["Qstar_minus_cache_correction_exact"]))
    upper=packet["whole_box_hessian"]["H_I_minus_H_cache_upper_exact"]
    require(all(type(i) is type(j) is int for i,j,*_ in upper),"literal upper H indexes")
    supplied_h={(i,j):interval((l,h)) for i,j,l,h in upper}
    require(list(supplied_h)==list(combinations_with_replacement(range(len(sel)),2))
            and len(supplied_h)==len(upper) and supplied_h.keys()==exact_h.keys(),"exact upper Hessian map")
    local_failure={
        "direct_point_cache":[i for i,(outer,v) in enumerate(zip(supplied_c,cf,strict=True))
                              if not contain(outer,ends(v))],
        "point_force":[i for i,(outer,inner) in enumerate(zip(supplied_f,exact_f,strict=True))
                       if not contain(outer,inner)],
        "whole_box_H":[list(k) for k,inner in exact_h.items() if not contain(supplied_h[k],inner)]}
    return {"kind":kind,"slot":i,"global_indices":glob,"checked_cache_sites":len(maps),
            "checked_force_rows":len(exact_f),"checked_upper_H_entries":len(exact_h),
            "containment_failures":local_failure,
            "independent_point_F_exact":[list(map(str,r)) for r in exact_f],
            "independent_whole_box_H_exact":[[*k,*map(str,v)] for k,v in sorted(exact_h.items())],
            "seconds":perf_counter()-started}
