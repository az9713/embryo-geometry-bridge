"""Replay singleton body and trace geometry. No integration or force evaluation."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import ceil
import os
from pathlib import Path
import sys
import time
import traceback

SOURCE=Path(__file__).resolve()
HERE=SOURCE.parent
BASE=HERE.parent
PRIMARY=BASE.parents[1]
ROOT=PRIMARY.parents[1]
PRE=BASE/"preparation"
ACTUAL=BASE/"geometry"
REPLAY=HERE/"replay"
M_SOURCE=ROOT/"scripts/transfer_ambiguity_exact_targets.py"
BEFORE=PRE/"prelaunch_independent_review.json"
CACHE=PRIMARY/"certification_design/cached_pilots_v1/independent_actual_preparation/actual_replay/receipt.json"
PINS={
    M_SOURCE:"53ff92a87d164871e04323759b27e015fc8e19d039dc0dd15c2d6080acda7cc8",
    PRE/"protocol.json":"28efd2467bbeae4da32797e0405b9a9595cdbb03ded9ce94221f18592305687c",
    PRE/"target_manifest.json":"cb25fbeaef753cea700c187cd8b377b7dcf4aa1fe2b0b246ed17134a07a89bc3",
    PRE/"preflight.json":"a4dacbab0dafb730db12321892d741979fea4654c436996b964cc9cfde1b0808",
    BEFORE:"17c0e8feb5f91898787994ba84bd1afd65e4be5e8ff52f01b0f5a48853f3ce6b",
    CACHE:"4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17",
}
OUTPUT={
    ACTUAL/"results.json":"fabf695e7f83189b2a088570723f1c103cafa7c5e973ac704f9c7adc25b9bdc9",
    ACTUAL/"terminal.json":"060ffa111367f95bbb7eccf17084e888ccf6e0036edae6f8f4246ebd47bccad1",
    ACTUAL/"started.json":"e2cfeee2097d97f695a014d258139a08c05194a86da9060323d1818aa3ab28b7",
}
R=(Q(1,2**36),Q(1,2**34))
RHO_DENOM=2**48
CHART=("off_axis","axis_vertex","axis_edge_physical")
SCOPE=("Exact spatial whole-box body and surface geometry for four fixed model inputs at two "
       "coded radii. This supplies closed domains and geometry bounds. No spatial derivatives, "
       "selector or tail bounds, force or Hessian, stability, root, continuum or source inference.")


def name(p):
    return p.relative_to(ROOT).as_posix()


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def pinned(p,d):
    raw=p.read_bytes()
    assert sha256(raw).hexdigest()==d,p
    return json.loads(raw)


def save(p,obj):
    with p.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(obj,f,indent=2,allow_nan=False)
        f.write("\n")
    return sha(p)


def yes(d):
    return bool(d) and all(x is True for x in d.values())


def digest(obj):
    return sha256(json.dumps(obj,sort_keys=True,separators=(",",":")).encode()).hexdigest()


def pair(p):
    a,b=map(Q,p)
    assert a<=b
    return a,b


def hull(pairs):
    pairs=list(pairs)
    return min(a for a,_ in pairs),max(b for _,b in pairs)


def product(a,b):
    terms=[x*y for x in a for y in b]
    return min(terms),max(terms)


def square(a):
    lo,hi=a
    return (Q(0) if lo<=0<=hi else min(lo*lo,hi*hi)),max(lo*lo,hi*hi)


def add(a,b):
    return a[0]+b[0],a[1]+b[1]


def combination(weights,hulls):
    vals=[product((u,u),v) for u,v in zip(weights,hulls,strict=True)]
    return sum(a for a,_ in vals),sum(b for _,b in vals)


def procedure():
    for p,d in PINS.items():
        assert sha(p)==d
    pre=pinned(PRE/"protocol.json",PINS[PRE/"protocol.json"])
    prelaunch=pinned(BEFORE,PINS[BEFORE])
    assert prelaunch["status"]=="PASS" and yes(prelaunch["checks"])
    assert pre["new_spatial_integrals"]==pre["new_cached_AD_assemblies"]==0
    man=pinned(PRE/"target_manifest.json",PINS[PRE/"target_manifest.json"])
    assert tuple(map(Q,pre["radii_exact"]))==tuple(map(Q,man["radii_exact"]))==R
    assert len(man["targets"])==4 and man["pilot_order"]==[r["pilot_id"] for r in man["targets"]]
    return {"source_sha256":sha(SOURCE),"input_hashes":{name(p):h for p,h in PINS.items()},
            "producer_terminal_and_result_external_pins":{name(p):h for p,h in OUTPUT.items()},
            "producer_pid":21768,"producer_exec_session":52282,"producer_root_observed_exit":0,
            "pilot_order":man["pilot_order"],"radii_exact":list(map(str,R)),"scope":SCOPE,
            "model_authority":"Current independently certified cached points, exact optical drives and explicit fixed material coefficients.",
            "shared_engine_scope":"Re-execute the previously reviewed body Bernstein and trace Bernstein routines from independently rebuilt exact boxes. Verify independent arithmetic summaries and surface cache line intervals. Reuse fixed regularized point chart basis transforms."}


def prepare():
    plan=procedure()
    x=Q(1,32)
    assert product((-x,2*x),(3*x,5*x))==(-5*x*x,10*x*x)
    assert square((-x,2*x))==(0,4*x*x)
    assert combination((Q(1),-Q(2)),((0,x),(x,3*x)))==(-6*x,-x)
    raw=Q(1,3)
    ceil_rho=Q(ceil(raw*RHO_DENOM),RHO_DENOM)
    assert raw<=ceil_rho<raw+Q(1,RHO_DENOM)
    h=save(HERE/"protocol.json",plan)
    f=save(HERE/"preflight.json",{"status":"PASS","protocol_sha256":h,"source_sha256":sha(SOURCE),
            "metadata_from_root_observation":True,"actual_geometry_payloads_opened":0,
            "force_moment_or_AD_jobs":0,"independent_interval_arithmetic_fixtures":True})
    print(h,f,sha(SOURCE),flush=True)


def admit(plan):
    starter,t,r=(pinned(ACTUAL/f"{k}.json",OUTPUT[ACTUAL/f"{k}.json"]) for k in ("started","terminal","results"))
    gate=pinned(BEFORE,PINS[BEFORE])
    allowed=("source_sha256","protocol_sha256","preflight_sha256","target_manifest_sha256","cached_root_independent_review_sha256")
    bindings={k:gate[k] for k in allowed}|{"prelaunch_review_sha256":PINS[BEFORE]}
    assert starter["bindings"]==t["bindings"]==r["bindings"]==bindings
    assert starter["pid"]==t["pid"]==plan["producer_pid"]
    assert t["status"]=="FINISHED" and type(t["exit_code"]) is int and t["exit_code"]==0
    assert t["results_sha256"]==OUTPUT[ACTUAL/"results.json"]
    assert r["status"]=="COMPLETE_FRESH_GEOMETRY_ATTEMPTS_UNREVIEWED"
    assert r["independent_actual_geometry_review_required"] is True and r["force_or_H_jobs_authorized"] is False
    schedule=[(pid,i) for pid in plan["pilot_order"] for i in range(2)]
    assert len(r["geometry_cases"])==len(schedule)==8
    Hash={name(p):h for p,h in OUTPUT.items()}
    out={}
    for item,(pid,i) in zip(r["geometry_cases"],schedule,strict=True):
        d=ACTUAL/pid/f"box_{i}"
        assert item["pilot_id"]==pid and item["file"]==name(d/"geometry.json")
        Hash[item["file"]]=item["sha256"]
        geometry=pinned(d/"geometry.json",item["sha256"])
        assert geometry["status"]==item["status"]=="PASS_FRESH_GEOMETRY_ONLY_UNREVIEWED"
        Hash[name(d/"body_geometry.json")]=geometry["body_geometry_sha256"]
        out[(pid,i)]=(geometry,pinned(d/"body_geometry.json",geometry["body_geometry_sha256"]))
    return out,Hash


def coordinate_box(ctx,stored,pilot,proof,radius):
    a=ctx.target.a
    assert proof["field"]==pilot["field"] and proof["field_sha256"]==pilot["field_sha256"]
    state=pinned(ROOT/proof["field"],proof["field_sha256"])
    ref=[list(map(Q,one)) for one in state["reference_nodes"]]
    assert ref==[list(map(Q,one)) for one in a["reference"].tolist()]
    u=list(map(Q,state["free_displacement"]))
    assert u==list(map(Q,a["u"].tolist()))
    assert len(ref)==325 and len(u)==611
    rim=Q(state["parameters"]["outer_radius"])
    free=[2*i+k for i,(x,z) in enumerate(ref) for k in (0,1) if x!=rim and (x!=0 or k!=0)]
    assert free==a["free"].tolist()==stored["free"]==pilot["free_map"]
    fixed=set(range(650))-set(free)
    x=[row.copy() for row in ref]
    eps=[[Q(0),Q(0)] for _ in x]
    for i,d in zip(free,u,strict=True):
        j,k=divmod(i,2)
        x[j][k]+=d
        eps[j][k]=radius
    assert x==[list(map(Q,one)) for one in stored["current_centers_exact"]]
    assert ref==[list(map(Q,one)) for one in stored["reference_nodes_exact"]]
    assert eps==[list(map(Q,one)) for one in stored["coordinate_radii_exact"]]
    assert all(eps[i//2][i%2]==0 for i in fixed)
    assert all(eps[i][0]==x[i][0]==0 for i,(r,z) in enumerate(ref) if r==0)
    assert [[pair(h) for h in row] for row in stored["coordinate_hulls_exact"]]==[
        [(v-e,v+e) for v,e in zip(row,er,strict=True)] for row,er in zip(x,eps,strict=True)]
    rounded=[list(map(Q,one)) for one in state["current_nodes"]]
    assert max(abs(x[i][j]-rounded[i][j]) for i in range(325) for j in (0,1))==Q(pilot["xi_dimensionless_exact"])
    tri=state["triangles"]
    assert tri==a["triangles"].tolist()==stored["triangles"] and len(tri)==144
    assert stored["root_membership_assumed"] is stored["any_spatial_derivative_bound_assumed"] is False
    assert Q(stored["corrector_radius_exact"])==radius
    assert len(stored["predictor_hulls_exact"])==611 and all(pair(d)==(0,0) for d in stored["predictor_hulls_exact"])
    return ref,x,eps,tri,{j:i for i,j in enumerate(free)}


def check_body(body):
    assert body["status"]=="CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING"
    assert all(body[k] is True for k in ("all_element_J_positive","uniform_strong_monotonicity",
                                        "axis_preserved_for_entire_box","global_injectivity_for_entire_box"))
    assert body["reference_topology"]["convex_reference_rectangle"]
    assert body["reference_topology"]["conforming_P2_traces"]
    j=[]
    lam=[]
    assert len(body["triangles"])==144
    for b in body["triangles"]:
        assert b["positive_J_for_entire_box"] is b["axis_zero_coefficients_exact_for_entire_box"] is True
        assert b["strong_monotonicity_for_entire_box"] is True
        aa=[]
        assert len(b["correlated_degree4_coefficients"])==15
        for t in b["correlated_degree4_coefficients"]:
            assert sum(t["index"])==4
            num=pair((t["numerator"]["lower_exact"],t["numerator"]["upper_exact"]))
            den=Q(t["denominator_exact"])
            assert den>=0
            if den==0:
                assert num==(0,0) and t["ratio"] is None
            else:
                ratio=tuple(x/den for x in num)
                assert ratio==pair((t["ratio"]["lower_exact"],t["ratio"]["upper_exact"]))
                aa.append(ratio)
        dom=hull(aa)
        assert dom==pair((b["J_bound"]["lower_exact"],b["J_bound"]["upper_exact"]))
        assert dom[0]>Q(1,1_000_000)
        j.append(dom)
        mu=[]
        assert len(b["vertex_symmetric_gradient_bounds"])==3
        for v in b["vertex_symmetric_gradient_bounds"]:
            A,B,D=[pair((x["lower_exact"],x["upper_exact"])) for x in v["symmetric_gradient"]]
            det=product(A,D)[0]-square(B)[1]
            tr=A[1]+D[1]
            assert v["positive_definite_for_entire_box"] is True and A[0]>0 and det>0 and tr>0
            assert Q(v["lambda_lower_exact"])==det/tr
            mu.append(det/tr)
        assert min(mu)==Q(b["strong_monotonicity_lower_exact"])
        lam.append(min(mu))
    glob=hull(j)
    assert glob==pair((body["J_bound"]["lower_exact"],body["J_bound"]["upper_exact"]))
    assert min(lam)==Q(body["strong_monotonicity_lower_exact"])>0
    return j,glob,min(lam)


def surface_bounds(d):
    assert d["status"]=="CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN"
    width=2**d["depth"]
    assert 0<=d["depth"]<=8 and len(d["leaves"])==width
    entries=[]
    hoops=[]
    for i,l in enumerate(d["leaves"]):
        assert pair(l["s_exact"])==(Q(i,width),Q(i+1,width))
        num=[hull(map(pair,coeff)) for coeff in l["q_numerator_Bernstein_coefficients_exact"]]
        den=hull(map(pair,l["hoop_denominator_Bernstein_coefficients_exact"]))
        assert den[0]>0
        num[0]=product(num[0],(1/den[1],1/den[0]))
        assert num==list(map(pair,l["q_bounds_exact"]))
        assert num[0][0]>0
        A,B,C=map(square,num)
        value=product(A,add(B,C))
        assert value==pair(l["G_bounds_exact"]) and value[0]>0
        hoops.append(num[0][0])
        entries.append(value)
    lo,hi=hull(entries)
    rho=(hi-lo)/(hi+lo)
    assert pair(d["G_bounds_exact"])==(lo,hi) and Q(d["hoop_lower_exact"])==min(hoops)>0
    assert Q(d["q0_exact"])==(lo+hi)/2 and Q(d["rho_exact"])==rho<1
    return lo,hi,rho,min(hoops),width


def independent_site(S,target,box,f,e,iq,nodes,glob,site):
    maps=S.L.actual_maps(target,f,e,iq)
    assert maps["nodes"]==site["nodes"]==nodes
    assert maps["global_indices"]==site["global_indices"]==glob
    assert maps["selected_local"]==site["selected_local"]
    assert iq==site["quadrature_index"]
    H=[pair(box["coordinate_hulls_exact"][i][k]) for i in nodes for k in (0,1)]
    star,cached=maps["star_full"],maps["cache_normalized_full"]
    delta=[[x-y for x,y in zip(a,b,strict=True)] for a,b in zip(star,cached,strict=True)]
    D=[combination(coef,H) for coef in delta]
    C=[combination(coef,H) for coef in cached]
    whole=[(a+min(c,Q(0)),b+max(d,Q(0))) for (a,b),(c,d) in zip(C,D,strict=True)]
    norm2=square(whole[1])[0]+square(whole[2])[0]
    assert list(map(pair,site["cache_q_exact"]))==C
    assert list(map(pair,site["delta_q_exact"]))==D
    assert list(map(pair,site["convex_line_exact"]))==whole
    assert whole[0][0]>0 and Q(site["tangent_norm_squared_lower_exact"])==norm2>0
    assert list(map(Q,site["normalization_D_exact"]))==list(map(Q,maps["normalization_D_exact"]))
    assert Q(site["auxiliary_normalized_weight_exact"])==maps["auxiliary_normalized_weight_exact"]>0
    assert site["status"]=="CERTIFIED_LOCAL_FULLY_COEFFICIENTWISE_LINE_DOMAIN_NO_CACHE_SUM"
    return norm2


def proj(nodes,free):
    ids=[2*i+k for i in nodes for k in (0,1)]
    return [i for i,j in enumerate(ids) if j in free],[free[j] for j in ids if j in free]


def replay_case(m, ctx, rule, entry, pilot, proof, i, records, hashes):
    d,b=records
    r=R[i]
    for k in ("routes","rule","model_contract","box_0","box_1","predictor","row"):
        hashes[entry[k]]=entry[k+"_sha256"]
    for k in ("routes","rule","model_contract"):
        assert d[k+"_sha256"]==entry[k+"_sha256"]
    assert d["box_sha256"]==entry[f"box_{i}_sha256"]
    box=pinned(ROOT/entry[f"box_{i}"],d["box_sha256"])
    route=pinned(ROOT/entry["routes"],d["routes_sha256"])
    contract=pinned(ROOT/entry["model_contract"],d["model_contract_sha256"])
    assert contract["larger_radius_cached_derivative_asserted"] is False
    assert Q(contract["cached_root_radius_exact"])==Q(proof["radius_exact"])==R[0]
    assert tuple(map(Q,contract["geometry_radii_exact"]))==R
    assert Q(d["radius_exact"])==r and d["pilot_id"]==ctx.pilot_id and d["record_id"]==ctx.record_id
    assert all(d[k] is False for k in ("local_component_selector_frozen","new_force_or_hessian_work",
                                       "new_cached_AD_assembly","root_asserted"))
    ref,x,eps,tri,free=coordinate_box(ctx,box,pilot,proof,r)
    assert m.G.certify_box(ref,x,tri,eps)==b
    jb,J,mono=check_body(b)
    S=m.E.SUR
    nodes=S.common_nodes(ctx.target,box)
    assert len(route["bulk"])==len(d["bulk_geometry_only"])==144
    classes=Counter()
    rho_max=Q(0)
    for k,(source,out,ids,Jexact) in enumerate(zip(route["bulk"],d["bulk_geometry_only"],tri,jb,strict=True)):
        fields,geom_ids=m.E._bulk_fields(ctx.target,nodes,k)
        chart=CHART[sum(ref[j][0]==0 for j in ids[:3])]
        classes[chart]+=1
        assert source["nodes"]==geom_ids==ids
        assert source["chart"]==fields["chart"]["kind"]==chart
        assert source["element"]==k and source["auxiliary_element_rule_sha256"]==m.E._bulk_rule(ctx.target,rule,k)
        assert (source["free_local"],source["global_indices"])==proj(ids,free)==(fields["free_local"],fields["global_indices"])
        assert Q(source["active"])==fields["active"]==int(k in pilot["active_elements"])
        raw=max(abs(j-1) for j in Jexact)
        rounded=Q(ceil(raw*RHO_DENOM),RHO_DENOM)
        assert Q(out["raw_rho_exact"])==raw<=rounded==Q(out["rho_exact"])<1
        assert pair(out["J_bounds_exact"])==Jexact
        assert all(source[c]==out[c] for c in source)
        rho_max=max(rho_max,rounded)
    assert dict(classes)=={"off_axis":132,"axis_vertex":6,"axis_edge_physical":6}
    assert len(route["surface"])==len(d["trace_geometry_only"])==24
    surface_rho, g_min, hoop_min, norm_min=[],[],[],[]
    leaves=0
    for i0,(s,t) in enumerate(zip(route["surface"],d["trace_geometry_only"],strict=True)):
        f,e=divmod(i0,12)
        ids=ctx.target.a[f"surface_{f}_edges"][e].tolist()
        assert s["nodes"]==ids and s["face"]==f and s["edge"]==e
        assert all(s[k]==t[k] for k in s)
        m.E._surface_rule(ctx.target,rule,f,e)
        assert s["trace_rule_sha256"]==digest(m.doc([rule["surface_s"][f],rule["surface_weights"][f][e]]))
        fields=ctx.target.surface(f,e,nodes=nodes)
        assert Q(s["Rlo_exact"])==fields["Rlo"]==ref[ids[0]][0]
        assert Q(s["Rhi_exact"])==fields["Rhi"]==ref[ids[-1]][0]
        assert (s["free_local"],s["global_indices"])==proj(ids,free)==(fields["free_local"],fields["global_indices"])
        assert S.uniform_G_bounds(fields,max_depth=8)==t["domain"]
        Glo,_,rho,hoop,count=surface_bounds(t["domain"])
        g_min.append(Glo)
        surface_rho.append(rho)
        hoop_min.append(hoop)
        leaves+=count
        assert len(t["cache_site_geometry"])==12
        for iq,site in enumerate(t["cache_site_geometry"]):
            norm_min.append(independent_site(S,ctx.target,box,f,e,iq,ids,s["global_indices"],site))
    assert len(norm_min)==288
    return {"status":"PASS_FRESH_SINGLETON_WHOLE_GEOMETRY_ONLY","pilot_id":ctx.pilot_id,"record_id":ctx.record_id,
            "radius_exact":str(r),"box_sha256":d["box_sha256"],"body_geometry_sha256":d["body_geometry_sha256"],
            "proof_center_field_sha256":proof["field_sha256"],"cached_proof_radius_exact":proof["radius_exact"],
            "larger_cached_derivative_proof_reused":False,"all_144_body_and_24_trace_domains":True,
            "raw_surface_cache_geometry_lines":288,"body_J_bounds_exact":list(map(str,J)),
            "strong_monotonicity_lower_exact":str(mono),"max_bulk_rho_outward_exact":str(rho_max),
            "max_surface_rho_raw_exact":str(max(surface_rho)),"min_surface_G_exact":str(min(g_min)),
            "min_trace_hoop_exact":str(min(hoop_min)),"min_site_tangent_squared_exact":str(min(norm_min)),
            "surface_Bernstein_leaves":leaves,"scope":SCOPE}


def run(ph,fh):
    plan=pinned(HERE/"protocol.json",ph)
    assert plan==procedure()
    pf=pinned(HERE/"preflight.json",fh)
    assert pf["status"]=="PASS" and pf["source_sha256"]==sha(SOURCE) and pf["protocol_sha256"]==ph
    assert pf["force_moment_or_AD_jobs"]==pf["actual_geometry_payloads_opened"]==0
    cases, hashes=admit(plan)
    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_exact_targets as m
    assert Path(m.__file__).resolve()==M_SOURCE
    REPLAY.mkdir(exist_ok=False)
    tag={"pid":os.getpid(),"source_sha256":sha(SOURCE),"protocol_sha256":ph,"preflight_sha256":fh}
    save(REPLAY/"started.json",tag|{"status":"RUNNING"})
    terminal=tag|{"status":"FAILED","exit_code":1}
    try:
        started=time.perf_counter()
        _,manifest,data=m.prepared()
        rows=[]
        with m.no_integrals():
            for pid in plan["pilot_order"]:
                ctx,rule,entry=m.open_context(pid,manifest,data)
                pilot,_,r,path=data[pid]
                proof=pinned(ROOT/path,r["receipt_sha256"])
                hashes[path]=r["receipt_sha256"]
                hashes[proof["field"]]=proof["field_sha256"]
                for i in range(2):
                    row=replay_case(m,ctx,rule,entry,pilot,proof,i,cases[pid,i],hashes)
                    f=REPLAY/f"{pid}_box_{i}.json"
                    hashes[name(f)]=save(f,row)
                    summary={k:row[k] for k in ("pilot_id","radius_exact","box_sha256",
                                                   "body_geometry_sha256","scope")}
                    rows.append(summary|{"receipt":name(f),"sha256":hashes[name(f)]})
                    print(json.dumps({"pilot_id":pid,"box":i,"J":[float(Q(x)) for x in row["body_J_bounds_exact"]],
                                      "min_embedding":float(Q(row["strong_monotonicity_lower_exact"])),
                                      "max_rho_bulk":float(Q(row["max_bulk_rho_outward_exact"])),
                                      "max_rho_surface":float(Q(row["max_surface_rho_raw_exact"]))}),flush=True)
        for f,h in hashes.items():
            assert sha(ROOT/f)==h
        assert len(rows)==8
        d=save(REPLAY/"receipt.json",{"status":"PASS_INDEPENDENT_ACTUAL_EIGHT_SINGLETON_GEOMETRY_BOXES",
            "source_and_phase":tag,"scope":SCOPE,
            "observed_producer_os":{"pid":21768,"session":52282,"exit_code":0,"observer":"root"},
            "checks":{"current_producer_and_all_proof_bytes_match":True,
                      "all_8_650_coordinate_hulls_independently_rebuilt":True,
                      "all_1152_triangle_and_full_box_embeddings_recomputed":True,
                      "all_192_surface_hoop_G_domains_and_2304_lines_recomputed":True,
                      "diagram_chart_and_rule_bindings_match":True},
            "geometry_records":rows,"checked_file_hashes":hashes,
            "premise_reuse":plan["shared_engine_scope"],"exact_spatial_root_count":0,
            "spatial_derivative_or_tail_bound_count":0,"force_hessian_or_AD_jobs":0,
            "seconds":time.perf_counter()-started})
        terminal.update(status="FINISHED",exit_code=0,receipt_sha256=d)
    except BaseException as err:
        save(REPLAY/"failure.json",{"error":str(err),"traceback":traceback.format_exc()})
        raise
    finally:
        save(REPLAY/"terminal.json",terminal)
        print(terminal,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="stage",required=True)
    sub.add_parser("prepare")
    r=sub.add_parser("run")
    r.add_argument("--protocol-sha256",required=True)
    r.add_argument("--preflight-sha256",required=True)
    a=p.parse_args()
    if a.stage=="prepare":
        prepare()
    else:
        run(a.protocol_sha256,a.preflight_sha256)
