"""Independent closed-envelope body, trace, and cache-line geometry replay."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
import json
from math import ceil, comb, factorial
import os
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PILOT=HERE.parent/"exact_loading_tube_pilot_preparation_v1"
COLD=HERE.parent/"exact_loading_tube_pilot_geometry_source_review_v1"
BODY_REVIEW=ROOT/"runs/exact_integration_bridge_v1/retube_B_geometry_independent_review/review.py"
PINS={PILOT/k:v for k,v in {
    "pilot.py":"d39ae6c3bfa3d932eaba03ddd59fc275fd5b3b1c094e7f9a544810d8c10dc8bc",
    "protocol.json":"fe64074d1a169a7672039c73dc07ca35c728cc1c299cf2a28025e0a0e7a97e43",
    "preflight.json":"911f6839bd9cc49dafe626e688cc30d88a283ad87b5aaf92fecc5b629a684881",
    "plan.json":"c3b1870e5ad815bc5df5ec94d4731b0228e529ba2e0651f847102b8e0ecaa277",
    "envelope.json":"f50c04b1802176abf4f5494898c8128018d4db8f8d0bb410af7774e5a3fff925",
    "independent_geometry_prelaunch.json":"bfe9374e0f8600510bd197c6bc82111a49a881417f358a176b49465f0be6a90e",
    "root_cold_readback.json":"990155d7cb03a68cda6d028e0facd48228b260c1bb83ffe146af671f4d988b15",
}.items()}
PINS.update({COLD/"audit.py":None,BODY_REVIEW:"14e0a6de53c172b0624842db420ca9dc17b3f567334e4b7b945efa226e2be688"})
CHARTS=("off_axis","axis_vertex","axis_edge_physical")
BODY_PASS="CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING"
TRACE_PASS="CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN"
PID="negative_KL_exact_spatial_loading_secant_v1"
ZERO=(Q(0),Q(0))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def rel(path):
    return path.resolve().relative_to(ROOT).as_posix()


def load(path, expected, links):
    b=path.read_bytes()
    require(sha256(b).hexdigest()==expected and links.setdefault(rel(path),expected)==expected,f"bound path {path}")
    return json.loads(b)


def save(path,data):
    with path.open("x",encoding="utf8",newline="\n") as f:
        json.dump(data,f,indent=2,default=str,allow_nan=False)
        f.write("\n")
    return sha(path)


def module(name,path,expected):
    require(sha(path)==expected,f"current source {path}")
    spec=spec_from_file_location(name,path)
    m=module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def interval(a):
    return Q(a["lower_exact"]),Q(a["upper_exact"])


def raw_pair(a):
    x,y=map(Q,a)
    require(x<=y,"ordered interval")
    return x,y


def hull(pairs):
    return min(x[0] for x in pairs),max(x[1] for x in pairs)


def pack(pair):
    return {"lower_exact":str(pair[0]),"upper_exact":str(pair[1])}


def axis_hulls(box):
    ref=[[Q(v) for v in row] for row in box["reference_nodes_exact"]]
    nodes=[[(Q(y)-Q(r),Q(y)+Q(r)) for y,r in zip(a,b,strict=True)]
           for a,b in zip(box["current_centers_exact"],box["coordinate_radii_exact"],strict=True)]
    require(len(nodes)==650//2 and all(n[0]==ZERO for p,n in zip(ref,nodes,strict=True) if p[0]==0),
            "325 nodes and exact fixed axis radii")
    return ref,nodes


def tri_replay(B,x,Y):
    """Use the separately reviewed shape polynomials and exact interval primitives."""
    a,b=x[1][0]-x[0][0],x[2][0]-x[0][0]
    c,d=x[1][1]-x[0][1],x[2][1]-x[0][1]
    det=a*d-b*c
    require(det>0 and min(x0[0] for x0 in x)>=0,"positive material triangle")
    inv=((d/det,-b/det),(-c/det,a/det))
    grad=[]
    vertices=[]
    for xi,eta in ((0,0),(1,0),(0,1)):
        D=[B.shape_derivatives(i,Q(xi),Q(eta)) for i in range(6)]
        g=[[B.sum_intervals(B.scalar(Y[k][i],D[k][j]) for k in range(6)) for j in range(2)] for i in range(2)]
        grad.append(g)
        F=[[B.sum_intervals(B.scalar(g[i][k],inv[k][j]) for k in (0,1)) for j in (0,1)] for i in (0,1)]
        aa,bb,dd=F[0][0],B.scalar(B.add(F[0][1],F[1][0]),Q(1,2)),F[1][1]
        lower=B.sub(B.mul(aa,dd),B.square(bb))[0]
        trace=aa[1]+dd[1]
        ok=aa[0]>0 and lower>0 and trace>0
        vertices.append({"symmetric_gradient":list(map(pack,(aa,bb,dd))),
          "positive_definite_for_entire_box":ok,"lambda_lower_exact":str(lower/trace) if ok else None})
    controls=[]
    for i,j in ((0,0),(1,1),(2,2),(0,1),(1,2),(0,2)):
        A,C=grad[i],grad[j]
        plus=B.add(B.mul(A[0][0],C[1][1]), B.mul(C[0][0],A[1][1]))
        minus=B.add(B.mul(A[0][1],C[1][0]), B.mul(C[0][1],A[1][0]))
        controls.append(B.scalar(B.sub(plus,minus),1/(2*det)))
    rr=[y[0] for y in Y[:3]]+[B.sub(B.scalar(Y[m][0],Q(2)),B.scalar(B.add(Y[i][0],Y[j][0]),Q(1,2)))
                               for i,j,m in ((0,1,3),(1,2,4),(2,0,5))]
    ind=((2,0,0),(0,2,0),(0,0,2),(1,1,0),(0,1,1),(1,0,1))
    coeff={}
    for a,u in zip(ind,controls,strict=True):
        for b,v in zip(ind,rr,strict=True):
            k=tuple(s+t for s,t in zip(a,b,strict=True))
            weight=Q(B.multinomial(a)*B.multinomial(b),B.multinomial(k))
            coeff[k]=B.add(coeff.get(k,ZERO),B.scalar(B.mul(u,v),weight))
    witnesses=[]
    ratios=[]
    axis_ok=True
    for k,v in sorted(coeff.items()):
        d=sum(Q(k[i],4)*x[i][0] for i in range(3))
        if d:
            r=B.scalar(v,1/d)
            ratios.append(r)
        else:
            r=None
            axis_ok=axis_ok and v==ZERO
        witnesses.append({"index":list(k),"numerator":pack(v),"denominator_exact":str(d),
                          "ratio":None if r is None else pack(r)})
    J=hull(ratios) if axis_ok and ratios else None
    s=all(v["positive_definite_for_entire_box"] for v in vertices)
    return {"positive_J_for_entire_box":J is not None and J[0]>0,
            "axis_zero_coefficients_exact_for_entire_box":axis_ok,"J_bound":pack(J) if J else None,
            "correlated_degree4_coefficients":witnesses,"strong_monotonicity_for_entire_box":s,
            "strong_monotonicity_lower_exact":str(min(Q(v["lambda_lower_exact"]) for v in vertices)) if s else None,
            "vertex_symmetric_gradient_bounds":vertices}


def topology(ref,tri):
    n=1+max(i for t in tri for i in t[:3])
    r,z=(sorted({x[k] for x in ref[:n]}) for k in (0,1))
    grid={tuple(x):i for i,x in enumerate(ref[:n])}
    require(len(r)*len(z)==n==len(grid) and r[0]==0,"convex rectangular material grid")
    cells=set()
    for a,b in zip(r[:-1],r[1:],strict=True):
        for c,d in zip(z[:-1],z[1:],strict=True):
            u,v,w,s=[grid[t] for t in ((a,c),(b,c),(b,d),(a,d))]
            cells.update((frozenset((u,v,w)),frozenset((u,w,s))))
    require(len(cells)==len(tri)==144 and cells=={frozenset(row[:3]) for row in tri},"complete structured rectangle")
    share={}
    for t in tri:
        require(len(t)==len(set(t))==6 and all(type(i) is int and 0<=i<len(ref) for i in t),"strict distinct P2 nodes")
        for a,b,m in ((0,1,3),(1,2,4),(2,0,5)):
            e=tuple(sorted((t[a],t[b])))
            require(t[m]>=n and share.setdefault(e,t[m])==t[m],"conforming trace")
    require(len(share)==len(set(share.values())),"one shared midpoint per material edge")


def body_replay(B,box,producer):
    ref,Y=axis_hulls(box)
    tri=box["triangles"]
    topology(ref,tri)
    require(len(producer["triangles"])==144,"all native triangles")
    rebuilt=[]
    for i,t in enumerate(tri):
        record=tri_replay(B,[ref[j] for j in t[:3]],[Y[j] for j in t])
        require(record==producer["triangles"][i],f"independent body coefficient mismatch at element {i}")
        rebuilt.append(record)
    goodJ=all(v["positive_J_for_entire_box"] for v in rebuilt)
    goodS=all(v["strong_monotonicity_for_entire_box"] for v in rebuilt)
    good=goodJ and goodS
    bJ=hull([interval(v["J_bound"]) for v in rebuilt]) if all(v["J_bound"] is not None for v in rebuilt) else None
    m=min(Q(v["strong_monotonicity_lower_exact"]) for v in rebuilt) if goodS else None
    require(producer["status"]==(BODY_PASS if good else "UNRESOLVED_GEOMETRY_BOX")
            and producer["all_element_J_positive"] is goodJ and producer["uniform_strong_monotonicity"] is goodS
            and producer["axis_preserved_for_entire_box"] is True and producer["global_injectivity_for_entire_box"] is good
            and producer["J_bound"]==(None if bJ is None else pack(bJ))
            and producer["strong_monotonicity_lower_exact"]==(None if m is None else str(m))
            and producer["reference_topology"]["convex_reference_rectangle"] is True
            and producer["reference_topology"]["conforming_P2_traces"] is True,"same independent embedding decision")
    return good,ref,Y,rebuilt,{"body_J_bounds_exact":bJ,"strong_monotonicity_lower_exact":m,
       "witnesses_failed_positive_bound":[i for i,r in enumerate(rebuilt) if not r["positive_J_for_entire_box"]],
       "witnesses_failed_monotonicity_bound":[i for i,r in enumerate(rebuilt) if not r["strong_monotonicity_for_entire_box"]]}




def trace_shape(B,Y,x0,x1,a,b):
    L=x1-x0
    require(L>0 and x0>=0 and 0<=a<b<=1,"physical edge and chart interval")
    N=((Q(1),Q(-3),Q(2)),(Q(0),Q(4),Q(-4)),(Q(0),Q(-1),Q(2)))
    P=lambda x:tuple(sum(v*x**k for k,v in enumerate(n)) for n in N)
    D=lambda x:tuple((n[1]+2*x*n[2])/L for n in N)
    z=P(a); w=P(b)
    if x0:
        # Exact restricted Bernstein coefficient for each original P2 basis function.
        center=P((a+b)/2)
        weights=(z,tuple(2*x-(y+w0)/2 for x,y,w0 in zip(center,z,w,strict=True)),w)
        rad=Y
        den=((x0+L*a,x0+L*a),(x0+L*b,x0+L*b))
    else:
        require(Y[0][0]==ZERO,"exact axis trace")
        weights=tuple(tuple((N[j][1]+N[j][2]*p)/L for j in range(3)) for p in (a,b))
        rad=Y
        den=((Q(1),Q(1)),)
    lin=lambda weights,component:B.sum_intervals(B.scalar(Y[j][component],w0) for j,w0 in enumerate(weights))
    exact_coeff=[[lin(weight,0) for weight in weights],
                 [lin(D(p),0) for p in (a,b)],[lin(D(p),1) for p in (a,b)]]
    denh=hull(den)
    hp=B.mul(hull(exact_coeff[0]),(1/denh[1],1/denh[0]))
    q=(hp,hull(exact_coeff[1]),hull(exact_coeff[2]))
    G=B.mul(B.square(q[0]),B.add(B.square(q[1]),B.square(q[2])))
    return q,G,exact_coeff,den


def contains(outer,inner):
    return outer[0]<=inner[0]<=inner[1]<=outer[1]


def trace_replay(B,out,Y,r0,r1):
    depth=out["depth"]
    require(type(depth) is int and 0<=depth<=8 and out["status"]==TRACE_PASS and len(out["leaves"])==2**depth,
            "exact complete closed trace partition")
    own_G=[];own_hoop=[];orig_G=[];orig_hoop=[]
    for i,w in enumerate(out["leaves"]):
        a,b=Q(i,2**depth),Q(i+1,2**depth)
        require(raw_pair(w["s_exact"])==(a,b),"ordered partition leaf")
        q,G,coeff,den=trace_shape(B,Y,r0,r1,a,b)
        # Zero leading coefficients can be absent from a producer's power polynomial.
        # Compare whole coefficient hulls, which enclose the independently reconstructed trace
        # by its Bernstein controls even when the representations differ in degree.
        n=tuple(hull(tuple(map(raw_pair,r))) for r in w["q_numerator_Bernstein_coefficients_exact"])
        d=hull(tuple(map(raw_pair,w["hoop_denominator_Bernstein_coefficients_exact"])))
        require(d[0]>0 and contains(d,hull(den))
                and all(contains(x,hull(v)) for x,v in zip(n,coeff,strict=True)),
                "producer numerator/denominator control hulls cover all own exact controls")
        q0=(B.mul(n[0],(1/d[1],1/d[0])),n[1],n[2])
        G0=B.mul(B.square(q0[0]),B.add(B.square(q0[1]),B.square(q0[2])))
        require(q0==tuple(map(raw_pair,w["q_bounds_exact"])) and raw_pair(w["G_bounds_exact"])==G0
                and q0[0][0]>0 and G0[0]>0 and contains(G0,G) and q[0][0]>0 and G[0]>0,
                "fresh exact own and producer root trace hulls")
        own_G.append(G);own_hoop.append(q[0][0]);orig_G.append(G0);orig_hoop.append(q0[0][0])
    L,U=hull(orig_G);lhoop=min(orig_hoop)
    rho=(U-L)/(U+L)
    require((L,U)==raw_pair(out["G_bounds_exact"]) and lhoop==Q(out["hoop_lower_exact"])
            and Q(out["q0_exact"])==(L+U)/2 and Q(out["rho_exact"])==rho<1,"rekindled trace cone and rho")
    return {"independent_G_bounds_exact":hull(own_G),"independent_hoop_lower_exact":min(own_hoop),
            "producer_rho_exact":rho,"cells":len(own_G)}


def independent_site(B, maps, Y, iq):
    def projection(rows):
        return [B.sum_intervals(B.scalar(y,v) for y,v in zip(Y,row,strict=True)) for row in rows]
    star,cache=maps["star_full"],maps["cache_normalized_full"]
    dq=[[s-c for s,c in zip(a,b,strict=True)] for a,b in zip(star,cache,strict=True)]
    c=projection(cache);d=projection(dq)
    line=[(x+min(0,u),y+max(0,v)) for (x,y),(u,v) in zip(c,d,strict=True)]
    tangent=B.square(line[1])[0]+B.square(line[2])[0]
    weight=maps["auxiliary_normalized_weight_exact"]
    require(line[0][0]>0 and tangent>0 and weight>0,"ambiguous actual cache-line hoop/tangent/weight")
    def pairs(a): return [list(map(str,p)) for p in a]
    return {"quadrature_index":iq,"nodes":maps["nodes"],"global_indices":maps["global_indices"],
            "selected_local":maps["selected_local"],"normalization_D_exact":list(map(str,maps["normalization_D_exact"])),
            "cache_q_exact":pairs(c),"delta_q_exact":pairs(d),"convex_line_exact":pairs(line),
            "tangent_norm_squared_lower_exact":str(tangent),"auxiliary_normalized_weight_exact":str(weight),
            "status":"CERTIFIED_LOCAL_FULLY_COEFFICIENTWISE_LINE_DOMAIN_NO_CACHE_SUM"}


def exact_Rho(j):
    raw=max(abs(j[0]-1),abs(j[1]-1))
    return raw,Q(ceil(raw*2**48),2**48)


def domains(A,B,ctx,plan,box,result,tris,ref,Y,routes):
    require(len(result["bulk"])==len(routes["bulk"])==144 and len(result["traces"])==len(routes["surface"])==24,
            "full deformed body and trace route census")
    free={index:i for i,index in enumerate(box["free"])}
    def route_only(saved,native,ids):
        dofs=[2*n+c for n in ids for c in (0,1)]
        positions=[i for i,k in enumerate(dofs) if k in free]
        require(all(native[k]==v for k,v in saved.items()) and ids==saved["nodes"]
                and positions==saved["free_local"] and [free[dofs[j]] for j in positions]==saved["global_indices"],
                "same immutable full local projection")
    missing=[];rhos=[];charted=Counter()
    for i,(route,row,tri,ids) in enumerate(zip(routes["bulk"],result["bulk"],tris,box["triangles"],strict=True)):
        route_only(route,row,ids)
        c=CHARTS[sum(ref[j][0]==0 for j in ids[:3])]
        require(route["element"]==i and route["chart"]==c and int(route["active"])==(i in plan["active_slots"]),
                "exact route from reference, mask and no circumferential alias")
        charted[c]+=1
        j=interval(tri["J_bound"]);raw,rho=exact_Rho(j)
        if j[0]<=Q(1,10**6) or not rho<1:
            require("unresolved" in row and "rho_exact" not in row,"retained bulk analytic failure")
            missing.append(["bulk",i])
        else:
            require(raw_pair(row["J_bounds_exact"])==j and Q(row["raw_rho_exact"])==raw and Q(row["rho_exact"])==rho,
                    "new exact bulk analytic interval")
            rhos.append(rho)
    require(charted==dict(zip(CHARTS,(132,6,6),strict=True)),"geometry chart census")
    traces=[]; tangents=[]
    for i,(route,row) in enumerate(zip(routes["surface"],result["traces"],strict=True)):
        face,e=divmod(i,12)
        ids=ctx.target.a[f"surface_{face}_edges"][e].tolist()
        route_only(route,row,ids)
        r0,r1=ref[ids[0]][0],ref[ids[-1]][0]
        require((route["face"],route["edge"])==(face,e) and (r0,r1)==(Q(route["Rlo_exact"]),Q(route["Rhi_exact"])),
                "same exact material boundary coordinates")
        y=[Y[j] for j in ids]
        if "domain" not in row:
            require("unresolved" in row and "surface_cache_line_geometry" not in row,"retain failed trace numeric phase")
            missing.append(["trace",face,e])
            continue
        d=trace_replay(B,row["domain"],y,r0,r1)
        traces.append({"face":face,"edge":e,**d})
        require(len(row["surface_cache_line_geometry"])==12,"all twelve exact surface site lines")
        for iq,site in enumerate(row["surface_cache_line_geometry"]):
            maps=A.E.SUR.L.actual_maps(ctx.target,face,e,iq)
            require(maps["nodes"]==ids and maps["global_indices"]==route["global_indices"],"same K frozen site map")
            ours=independent_site(B,maps,[v for node in y for v in node],iq)
            require(ours==site,"exact finite positive same-Y Qstar-to-cache geometry line")
            tangents.append(Q(ours["tangent_norm_squared_lower_exact"]))
    return {"full_domains_reconstructed":not missing,"bulk_chart_counts":dict(charted),
            "certified_bulk_rho_entries":len(rhos),"certified_surface_edges":len(traces),
            "certified_surface_site_lines":len(tangents),"unresolved_components":missing,
            "unresolved_surface_causes_not_independently_quantified":any(x[0]=="trace" for x in missing),
            "max_bulk_rho_exact":max(rhos) if rhos else None,
            "min_surface_site_tangent_square_lower_exact":min(tangents) if tangents else None,"trace_summary":traces}




def positive_int(x):
    return type(x) is int and x>0


def observed(o):
    return (o["status"]=="PASS_ROOT_OBSERVED_KL_LOADING_GEOMETRY_EXIT" and positive_int(o["pid"])
            and type(o["observed_exit_code"]) is int and o["observed_exit_code"]==0
            and (positive_int(o.get("unified_exec_session"))
                 or type(o.get("unified_exec_chunk_id")) is str and bool(o["unified_exec_chunk_id"].strip())))


def producer(actual, outer, b, links):
    o=load(actual/"os_exit_observation.json",outer,links)
    require(observed(o),"independent genuine root process exit required")
    term=load(actual/"terminal.json",o["terminal_sha256"],links)
    start=actual/"start.json"
    stamp={"bindings":b,"pilot_id":PID,"pid":o["pid"],"phase":"geometry_only",
           "independent_geometry_prelaunch_sha256":PINS[PILOT/"independent_geometry_prelaunch.json"]}
    def same(d):
        return positive_int(d.get("pid")) and all(d.get(k)==v for k,v in stamp.items())
    require(term["status"]=="KL_GEOMETRY_STAGE_FINISHED" and same(term)
            and o["result_sha256"]==term["result_sha256"],"observed complete matching geometry actual phase")
    s=load(start,o.get("start_sha256",o.get("started_sha256",sha(start))),links)
    require(same(s),"same completed start identity")
    d=load(actual/"result.json",o["result_sha256"],links)
    require(same(d) and d["status"]==term["mathematical_status"] and d["root_branch_claimed"] is False
            and all(type(d[k]) is int and d[k]==0 for k in ("new_AD_derivatives","new_integrations","new_solves")),
            "same mathematical status and no mechanical computation")
    body=load(actual/"body.json",d["body_sha256"],links)
    return o,d,body


def pinned_sources():
    links={}
    gate=load(PILOT/"independent_geometry_prelaunch.json",PINS[PILOT/"independent_geometry_prelaunch.json"],links)
    b=gate["bindings"]
    names={"source":"pilot.py","protocol":"protocol.json","preflight":"preflight.json","plan":"plan.json","envelope":"envelope.json"}
    for path,h in {**PINS,COLD/"audit.py":gate["independent_audit_source_sha256"]}.items():
        require(sha(path)==h,"frozen pilot or prior independent source")
        links[rel(path)]=h
    require(b=={k+"_sha256":PINS[PILOT/v] for k,v in names.items()}
            and gate["status"]=="PASS_INDEPENDENT_KL_LOADING_GEOMETRY_PRELAUNCH"
            and gate["pilot_id"]==PID and gate["geometry_only_authorized"] is True
            and gate["mechanical_kernels_authorized"] is False,"exact geometry-only source gate")
    phase=load(PILOT/"protocol.json",b["protocol_sha256"],links)
    require(all(gate["checks"].get(x) is True for x in phase["geometry_source_required_checks"]),
            "all required independently approved geometry checks")
    for k in ("geometry_source_tree_sha256","individual_authority_sha256","read_only_method_source_sha256"):
        for name,h in phase[k].items():
            require(sha(ROOT/name)==h,"same pinned sources and crosswalk")
            links[name]=h
    plan=load(PILOT/"plan.json",b["plan_sha256"],links)
    box=load(PILOT/"envelope.json",b["envelope_sha256"],links)
    C=module("_prior_reviewed_cold_envelope",COLD/"audit.py",gate["independent_audit_source_sha256"])
    margin,routes=C.envelope(None,plan,box)
    for name,h in zip(plan["anchor_box_files"],plan["anchor_box_sha256"],strict=True):
        load(ROOT/name,h,links)
    load(ROOT/plan["K_routes_file"],plan["K_routes_sha256"],links)
    B=module("_prior_reviewed_independent_shape_primitives",BODY_REVIEW,PINS[BODY_REVIEW])
    return b,plan,box,margin,routes,B,links


def fixtures(B):
    x=[[Q(0),Q(0)],[Q(1),Q(0)],[Q(0),Q(1)]]
    y=x+[[Q(1,2),Q(0)],[Q(1,2),Q(1,2)],[Q(0),Q(1,2)]]
    iy=[[(u,u) for u in a] for a in y]
    d=tri_replay(B,x,iy)
    require(interval(d["J_bound"])==(1,1) and Q(d["strong_monotonicity_lower_exact"])==Q(1,2)
            and d["axis_zero_coefficients_exact_for_entire_box"] is True,
            "own identity-triangle axis/cylindrical Jacobian")
    f=tri_replay(B,x,[[v[0],(-v[1][1],-v[1][0])] for v in iy])
    require(f["positive_J_for_entire_box"] is f["strong_monotonicity_for_entire_box"] is False
            and interval(f["J_bound"])==(-1,-1),
            "own independently unresolved or negative example")
    for r0 in (Q(0),Q(1)):
        y=[[(v,v),ZERO] for v in (r0,r0+Q(1,2),r0+1)]
        q,G,_,_=trace_shape(B,y,r0,r0+1,Q(0),Q(1))
        require(q[1]==(1,1) and q[2]==ZERO and G[0]>0,"own unweighted trace exact positive bounds")
    maps={"star_full":[[1,0],[1,0],[0,1]],"cache_normalized_full":[[1,0],[Q(1,2),0],[0,1]],
          "nodes":[0],"global_indices":[0,1],"selected_local":[0,1],
          "auxiliary_normalized_weight_exact":Q(3,4),"normalization_D_exact":[1,1,1]}
    s=independent_site(B,maps,[(Q(1),Q(2)),ZERO],0)
    require(Q(s["tangent_norm_squared_lower_exact"])==Q(1,4),"own full cache-to-auxiliary line")
    o={"status":"PASS_ROOT_OBSERVED_KL_LOADING_GEOMETRY_EXIT","observed_exit_code":0,
       "pid":4,"unified_exec_chunk_id":"inert"}
    require(observed(o) and not observed({**o,"pid":4.0}) and not observed({**o,"observed_exit_code":False})
            and not observed({**o,"unified_exec_chunk_id":""}),"typed observed process metadata")
    return {"exact_P2_axis_identity_and_reflected_box":True,"own_trace_axis_and_offaxis_guards":True,
            "exact_inert_cache_line":True,"typed_producer_OS_observation":True}


def expected_protocol(b,B,links):
    return {"status":"FROZEN_INDEPENDENT_KL_LOADING_GEOMETRY_ACTUAL_READBACK",
            "source_sha256":sha(Path(__file__)),"reviewed_bindings":b,
            "frozen_sources":links.copy(),"checks":fixtures(B),"opened_actual_payloads":0,
            "requires_externally_observed_producer_exit_before_payloads":True,
            "shared_primitives":"Reviewed independent exact P2 basis and interval operations; prior exact K context and same-Y cache maps."}


def prepare():
    b,_,_,_,_,B,links=pinned_sources()
    print(save(HERE/"protocol.json",expected_protocol(b,B,links)),flush=True)


def k_context(plan):
    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_remaining32_exact_targets as A
    require(Path(A.__file__).resolve()==ROOT/"scripts/transfer_ambiguity_remaining32_exact_targets.py",
            "same frozen exact K factory")
    A.X.bootstrap()
    case=plan["case_ids"][0]
    target=plan["K_binding"]
    ctx,_,r=A.open_context(case,"dfea134640f99b45e9ca75d1ce718043d67f2917375506997d2465abda8c1242")
    require(A.a_sha(ctx.R)==target["R_binary_sha256"] and A.a_sha(ctx.u)==target["center_binary_sha256"]
            and ctx.drive==Q(plan["nominal_drive_exact"]) and ctx.target.energy_scale==Q(plan["energy_scale_binary_exact"])
            and r["routes_sha256"]==plan["K_routes_sha256"],"same exact K center cache and routing")
    return A,ctx


def run(args):
    b,plan,box,margin,routes,B,links=pinned_sources()
    p=load(HERE/"protocol.json",args.protocol_sha256,{})
    require(p==expected_protocol(b,B,links),"same frozen verifier authority and inert fixtures")
    o,r,body=producer(PILOT/"actual_geometry",args.producer_os_sha256,b,links)
    out=HERE/"actual"
    out.mkdir(exist_ok=False)
    stamp={"source_sha256":sha(Path(__file__)),"protocol_sha256":args.protocol_sha256,
           "producer_os_sha256":args.producer_os_sha256,"producer_result_sha256":o["result_sha256"],
           "producer_terminal_sha256":o["terminal_sha256"],"producer_body_sha256":r["body_sha256"],
           "pid":os.getpid(),"producer_pid":o["pid"],"pilot_id":PID}
    save(out/"started.json",{"status":"RUNNING",**stamp})
    terminal={"status":"FAILED","exit_code":1,**stamp}
    try:
        good,ref,Y,tris,summary=body_replay(B,box,body)
        require(Q(r["radial_trial_bound_margin_normalized_exact"])==margin,"independent exact solver radial floor")
        if not good or margin<=0:
            require(r["status"]=="UNRESOLVED_KL_WHOLE_ENVELOPE_BODY_OR_RADIAL_BOUND"
                    and r["bulk"]==r["traces"]==[],"unresolved body must stop before all analytic routes")
            status="VERIFIED_UNRESOLVED_KL_WHOLE_ENVELOPE_BODY"
            details={"full_domains_reconstructed":False,"body_geometry_recomputed_in_full":True,
                     "bulk_and_trace_domains_not_reached":True}
        else:
            A,ctx=k_context(plan)
            with A.X.no_integrals():
                details=domains(A,B,ctx,plan,box,r,tris,ref,Y,routes)
            complete=details["full_domains_reconstructed"]
            require(r["status"]==("PASS_KL_FRESH_WHOLE_ENVELOPE_GEOMETRY_CANDIDATE" if complete else
                                  "UNRESOLVED_KL_WHOLE_ENVELOPE_BULK_OR_TRACE_DOMAIN"),"exact route completion standing")
            status="PASS_INDEPENDENT_FRESH_KL_LOADING_GEOMETRY_ONLY" if complete else "VERIFIED_UNRESOLVED_KL_ANALYTIC_ROUTES"
        require(all(sha(ROOT/name)==h for name,h in links.items()),"immutable closure bytes")
        h=save(out/"receipt.json",{"status":status,**stamp,"producer_mathematical_status":r["status"],
             "body_summary":summary,"geometry_scope":details,"independently_recomputed_body_triangles":144,
             "degree4_J_coefficients":2160,"vertex_spd_checks":432,"literal_radial_floor_exact":margin,
             "source_and_actual_sha256":links,"force_H_AD_solve_moment_or_spatial_integration_jobs":0,
             "root_or_loading_branch_proved":False})
        terminal.update(status="FINISHED",exit_code=0,receipt_sha256=h,mathematical_status=status)
        print(json.dumps({"status":status,"receipt_sha256":h}),flush=True)
    finally:
        save(out/"terminal.json",terminal)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    sp=cli.add_subparsers(dest="action",required=True)
    sp.add_parser("prepare")
    cmd=sp.add_parser("run")
    cmd.add_argument("--protocol-sha256",required=True)
    cmd.add_argument("--producer-os-sha256",required=True)
    args=cli.parse_args()
    prepare() if args.action=="prepare" else run(args)

