"""Preparation-only review of the record-neutral centered operator."""
import ast
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import Mock, patch

from flint import arb,ctx

OUT = Path(__file__).resolve().parent
BATCH = OUT.parent
ROOT = OUT.parents[3]
BASE = BATCH.parent.relative_to(ROOT).as_posix()
SOURCE = "scripts/exact_centered_cohort_bulk.py"
PINS = {
    "source_sha256":"524c22266a7e7e8138c06c8e4b444f44bac887d8404296d67d58a3541726d02c",
    "protocol_sha256":"dbe461e20104cf60cee9de78010d2c37d086663c420e3fb47869ce5ad2ea339d",
    "preflight_sha256":"9c4ea9ed77c22fda2a2b36ada0bc7ac3c1101f100f8d07d68a9a3bbb30af790c",
}
CD = {
    "optoGEF2_small_col_C":"9c3ff54d940735bae7db651056dfbf44c4dd893596520c4dfda2e276cb773d87",
    "optoGEF2_small_col_D":"13d51f9626f125eeb34f429a55cea05298ff2494643da679fc048d66ee772269",
}
L_ACTUAL = f"{BASE}/L_defect_diagnostic/centered_bulk_L_v1/actual_independent_review.json"
L_ACTUAL_SHA = "89df6cbaaaae0237b07e39c7cef05fbbeba20117d73a67153a4adcb459b6d6d0"


def require(v,label):
    if not v:
        raise ValueError(label)


def digest(b):
    return sha256(b).hexdigest()


def st(x):
    return digest(json.dumps(x,default=str,sort_keys=True,separators=(",",":")).encode())


def path(p):
    v=(ROOT/p).resolve()
    require(v.is_relative_to(ROOT) and v.relative_to(ROOT).as_posix()==p,"canonical local path")
    return v


def bound(p,h,parse=True):
    v=path(p).read_bytes()
    require(digest(v)==h,f"changed {p}")
    return json.loads(v) if parse else v


def yes(d):
    return isinstance(d,dict) and bool(d) and all(x is True for x in d.values())


def rejected(callback):
    try:
        callback()
    except ValueError:
        return True
    return False


def structural(A,L,F):
    methods={n.name:ast.unparse(n) for n in ast.parse(path(SOURCE).read_bytes()).body if isinstance(n,ast.FunctionDef)}
    require(all(x not in methods["bulk_component"] for x in ("_geometry(", "_source_inputs(", "_element(", "_compute(", "load_context("))
            and "_local(runtime, domain, selector, route, uniform, uniform_sha, element)" in methods["bulk_component"]
            and "_force_parts(fields, jet, fun, base, A(fields['active'] * ctx.drive), maps, selector)" in methods["bulk_component"]
            and "L.affine_gain_derivative(ctx.dose, fields['active'], fields, fun, maps, full_Y)" in methods["bulk_component"]
            and "new_Hessian_integrations': False" in methods["bulk_component"],
            "generic explicit numerical callbacks and scope")
    source=ast.parse(path(SOURCE).read_bytes())
    require(not any(isinstance(n,(ast.Global,ast.Nonlocal)) or
                    (isinstance(n,ast.Attribute) and isinstance(n.ctx,ast.Store)
                     and ast.unparse(n)!="precision.prec") for n in ast.walk(source)),
            "no referenced module assignment")
    require(A.L_CENTRED==Path(L.__file__).resolve().relative_to(ROOT).as_posix()
            and L.F is F and L.A is F.A and L.B.D is L.D
            and A.L_REDUCER==Path(L.B.__file__).resolve().relative_to(ROOT).as_posix(),
            "verified shared source import aliases")


def slopes_at_point(corners, lambdas):
    x,y,z=corners
    r1,Z1=y[0]-x[0],y[1]-x[1]
    r2,Z2=z[0]-x[0],z[1]-x[1]
    det=r1*Z2-r2*Z1
    require(det>0,"positive oriented affine reference map")
    dl1,dl2=-r2/det,r1/det
    dl0=-(dl1+dl2)
    d=(dl0,dl1,dl2)
    return [(4*l-1)*v for l,v in zip(lambdas,d)]+[
        4*(lambdas[i]*d[j]+lambdas[j]*d[i]) for i,j in ((0,1),(1,2),(2,0))]


def spot_axial(E,L,fields):
    corners=fields["corners"]
    lam=(Q(1,3),Q(1,4),Q(5,12))
    real=slopes_at_point(corners,lam)
    if fields["chart"]["kind"]=="axis_vertex":
        off=[i for i,p in enumerate(corners) if p[0] != 0]
        s=lam[off[0]]+lam[off[1]]
        coords=(s,lam[off[1]]/s)
    else:
        R,Z=(sum(a*p[k] for a,p in zip(lam,corners)) for k in (0,1))
        coords=(R,Z-fields["zcenter"])
    dz=L.axial_gradient(fields)
    require([E.T.evaluate(v,*coords) for v in dz]==real,
            "six full nodal d/dZ basis values from independent exact barycentric chain rule")
    for j,k in enumerate(fields["free_local"]):
        require(fields["B"][1 if k%2==0 else 3][j]==dz[k//2],"exact symbolic projection in original chart")


def norm(E,fields,poly):
    if fields["chart"]["kind"]!="off_axis":
        if fields["chart"]["kind"]=="axis_edge_physical":
            a=max(v[0] for v in fields["corners"])
            b=max(abs(v[1]-fields["zcenter"]) for v in fields["corners"])
            return sum(abs(c)*a**i*b**j for (i,j),c in poly.items())
        return sum(map(abs,poly.values()))
    a=min(v[0] for v in fields["corners"])
    b=max(v[0] for v in fields["corners"])
    z=max(abs(v[1]-fields["zcenter"]) for v in fields["corners"])
    return sum(abs(c)*(a if i<0 else b)**i*z**j for (i,j),c in poly.items())


def profile(E,L,dom,cases,routes,manifest,prepared):
    ctxr=dom.context
    center=ctxr.target.nodes()
    require(len(center)==325 and all(list(map(Q,coords))==actual for coords,actual in zip(dom.box["current_centers_exact"],center)),
            "entire exact selected current center")
    result={}
    cached_bary=dom.rule["volume_barycentric"]
    for i,p in enumerate(prepared["bulk_proofs"]):
        se,rt,man=cases[i],routes[i],manifest[i]
        v=bound(man["file"],man["sha256"])
        fields,ids=E._bulk_fields(ctxr.target,center,i)
        spot_axial(E,L,fields)
        require(i==v["element"]==p["element"] and p["record_id"]==v["record_id"]==ctxr.record_id
                and v["global_indices"]==p["global_indices"]==fields["global_indices"]
                and p["order"]==v["degree"]==se["degree"] and p["rho_exact"]==v["rho_exact"]==se["rho_exact"]
                and p["uniform_packet_sha256"]==man["sha256"] and p["original_route_proof_sha256"]==st(rt)
                and p["same_C_sha256"]==ctxr.row["parameter_cell_sha256"]
                and p["same_cache_sha256"]==ctxr.row["point_cache_sha256"]
                and p["same_predictor_sha256"]==ctxr.row["predictor_sha256"]
                and p["chart"]==fields["chart"]["kind"] and p["uniform_source_integration_chart"]==rt["integration_chart"]
                and p["local_nodes"]==ids and p["free_local"]==fields["free_local"],
                "per-element packet, identity, chart and reviewed route")
        gl,gh=ctxr.cell["common_gain"];g0=ctxr.base["common_gain"]
        require(p["common_gain_shift_exact"]==[str(gl-g0),str(gh-g0)] and Q(p["base_drive_exact"])==g0*ctxr.dose
                and Q(p["dose_exact"])==ctxr.dose and Q(p["active"])==fields["active"] in (0,1),"same exact gain path")
        def pp(poly):
            return [[i,j,str(x)] for (i,j),x in sorted(poly.items())]
        require(st({"q":[pp(x) for x in fields["q"]],"B":[[pp(x) for x in row] for row in fields["B"]]})
                ==p["point_polynomials_sha256"],"entire exact q/B field at record center")
        point_map=[(ids[j//2],j%2) for j in p["free_local"]]
        radii={(n,a):Q(dom.box["coordinate_radii_exact"][n][a]) for n in ids for a in (0,1)}
        require(p["original_local_radii_exact"]==[str(radii[k]) for k in point_map],"same local radii")
        for (n,a),r in radii.items():
            low,high=E.SUR.rational_hull(dom.nodes[n][a])
            y=center[n][a]
            require(low<=y-r<=y+r<=high and ((n,a) in point_map or r==0),"full same-Y path includes fixed coordinates")
        f=[[entry[k] for k in ("lower_exact","upper_exact")] for entry in v["force_integral_minus_cache"]]
        h=[[w["i"],w["j"],[w["enclosure"][k] for k in ("lower_exact","upper_exact")]]
           for w in sorted(v["H_integral_minus_cache_upper"],key=lambda w:(w["i"],w["j"]))]
        require(len(h)==len(f)*(len(f)+1)//2 and st(f)==p["uniform_force_bounds_sha256"]
                and st(h)==p["uniform_H_bounds_sha256"],"preserves original signed enclosures")
        weights=dom.rule["volume_weights"][i]
        require(all(w==Q.from_float(float(z))/ctxr.target.energy_scale>0
                    for w,z in zip(weights,ctxr.target.a["volume_weights"][i])) and len(weights)==len(cached_bary)==64
                and p["auxiliary_rule_sha256"]==st([[list(map(str,l)) for l in cached_bary],list(map(str,weights))])
                ==rt["auxiliary_element_rule_sha256"],"same normalized positive rule")
        j=E.JET.bulk_jacobian_jet(fields["q"],fields["B"])
        tail=Q(se["geometry_two_mass_upper_exact"])*Q(se["phi_prime_tail_upper_exact"])*max(norm(E,fields,x) for x in j["gradient"])
        require(tail==Q(p["point_force_tail_max_exact"])<=Q(p["force_tail_budget_exact"])
                ==Q(dom.protocol["force_tail_budget_exact"]),"point force majorant for both masses")
        result.setdefault(p["chart"],[]).append(i)
    require({k:len(v) for k,v in result.items()}=={"off_axis":132,"axis_vertex":6,"axis_edge_physical":6},
            "exact 132+6+6 route coverage")
    return result


def synthetic_gain(L):
    fields={"corners":((Q(1),Q(0)),(Q(3),Q(0)),(Q(1),Q(2))),"chart":{"kind":"off_axis"},
            "free_local":[0,1,4,9]}
    all_y=[Q((-1)**i*(i+1),8) for i in range(12)]
    lam=(Q(1,3),Q(1,4),Q(5,12))
    dz=slopes_at_point(fields["corners"],lam)
    X=Q(3,2);Z=Q(5,6)
    require((X,Z)==tuple(sum(l*p[k] for l,p in zip(lam,fields["corners"])) for k in (0,1)),"toy coordinates")
    z0=sum(p[1] for p in fields["corners"])/3
    # Synthetic moment functional is a signed evaluation at one exact reference point.
    T=L.T
    fake=SimpleNamespace(apply=lambda p:L.A(3*T.evaluate(p,X,Z-z0)))
    star=[[Q(0)]*12 for _ in range(5)]
    cache=[[Q(0)]*12 for _ in range(5)]
    for k in (1,3):
        for j in range((k-1)//2,12,2):
            star[k][j]=Q(j+2,7)
            cache[k][j]=Q(5-j,11)
    maps=[{"star_full":star,"cache_full":cache,"auxiliary_normalized_weight_exact":Q(7,13)}]
    base=[]
    for i in fields["free_local"]:
        k=1 if i%2==0 else 3
        expected=Q(5,17)*sum((3*dz[i//2]*dz[j//2]+Q(7,13)*(star[k][i]*star[k][j]-cache[k][i]*cache[k][j]))
                              *all_y[j] for j in range(i%2,12,2))
        base.append(expected)
    out=L.affine_gain_derivative(Q(5,17),1,fields,fake,maps,list(map(L.A,all_y)))
    require(all((L.A(q)-v).contains(0) for q,v in zip(base,out)), "exact affine dose factor, parity and full 12-coordinate gain")
    require(all(x.is_zero() for x in L.affine_gain_derivative(Q(5,17),0,fields,None,None,None)),
            "inactive gain exact zero before touching fields/moments")
    p=[arb(4),arb(7)];d=[arb(0,1),arb(0,2)]
    h={(0,0):(Q(-2),Q(-1)),(0,1):(Q(-1),Q(1)),(1,1):(Q(2),Q(3))}
    o,state,g=L.centered_force(p,d,h,[arb(1),arb(-2)],arb(0,1))
    require(all(v.contains(L.A(lo)) and v.contains(L.A(hi)) for v,lo,hi in zip(state,(-4,-7),(4,7)))
            and all(v.contains(L.A(lo)) and v.contains(L.A(hi)) for v,lo,hi in zip(o,(-1,-2),(9,16))),
            "signed Hessian symmetry, mixed term, gain endpoints and no dropped uncertain coefficient")


def mocked_component(A,runtime,dom,se,route,man,proof):
    R,op,L=runtime
    E,F=op.E,op.F
    n=len(proof["global_indices"])
    idx=proof["element"]
    data=bound(man["file"],man["sha256"])
    # No material moments or cache evaluations. The wrappers are observed with exact real input geometry.
    fake_fun=SimpleNamespace(requested_defect_keys=set())
    profile={"normalized_force_total":[arb(0,1)]*n,"two_mass_tail_upper_exact":[Q(0)]*n}
    mock_map=SimpleNamespace()
    with (patch.object(E,"_functional_type",return_value=Mock(return_value=fake_fun)) as ft,
          patch.object(E.U.C,"exact_bulk_maps",return_value=mock_map) as maps,
          patch.object(F,"_force_parts",return_value=profile) as point,
          patch.object(L,"affine_gain_derivative",return_value=[arb(0)]*n) as gain):
        packet=A.bulk_component(runtime,dom,se,route,data,man["sha256"],idx,proof)
        args=point.call_args.args
        current=dom.context.target.nodes()
        require(maps.call_count==64 and args[4].overlaps(L.A(Q(proof["active"])*dom.context.drive))
                and len(args[3])==12
                and all((a-L.A(current[node][j])).contains(0) for a,(node,j) in zip(args[3],[(n,j) for n in proof["local_nodes"] for j in (0,1)]))
                and all(a is b for a,b in zip(gain.call_args.args[-1],[dom.nodes[n][j] for n in proof["local_nodes"] for j in (0,1)]))
                and gain.call_args.args[0]==dom.context.dose and all(m is mock_map for m in args[5])
                and packet["new_Hessian_integrations"] is False and packet["centered_proof_sha256"]==st(proof),
                "explicit same record anchors, full wide coordinates, dose and cache maps for all actual callbacks")
        require(rejected(lambda:A.bulk_component(runtime,dom,se,route,data,man["sha256"],idx,
                                                {**proof,"base_drive_exact":"0"})) and point.call_count==1,
                "changed frozen proof must fail before any numerical callback")


def run():
    stem=BATCH.relative_to(ROOT).as_posix()
    bound(SOURCE,PINS["source_sha256"],False)
    bound(f"{stem}/exact_centered_cohort_bulk_frozen.py",PINS["source_sha256"],False)
    p,pf=[bound(f"{stem}/{n}.json",PINS[n+"_sha256"]) for n in ("protocol","preflight")]
    require(p["status"]=="FROZEN_REMAINING15_CENTERED_HELPER_PREPARATION_REVIEW_PENDING"
            and pf["status"]=="PASS_REMAINING15_CENTERED_HELPER_PREPARATION_ONLY" and yes(pf["checks"]) and len(pf["checks"])==19
            and p["source_sha256"]==pf["source_sha256"]==PINS["source_sha256"] and pf["protocol_sha256"]==PINS["protocol_sha256"]
            and p["actual_new_force_H_or_cache_sum_started"] is pf["new_force_H_or_cache_sum_started"] is False
            and p["completed_record_proof_hashes"]=={f"{stem}/{r}_proofs.json":h for r,h in CD.items()}
            and pf["records_with_bound_144_proofs"]==list(CD),"bound frozen C/D-only prep")
    for name,h in {**p["input_hashes"],**p["failed_attempts"]}.items():
        if "/components/" in name:
            require(any(f"/records/{r}/" in name for r in CD),"non-C/D actual packet in prep manifest")
        bound(name,h,False)
    l_actual=bound(L_ACTUAL,L_ACTUAL_SHA)
    require(l_actual["status"].startswith("PASS") and yes(l_actual["checks"]),"reviewed L primitives evidence")
    sys.path.insert(0,str(ROOT/"scripts"))
    import exact_centered_cohort_bulk as A
    R,op,L,plan,jobs,expect,hashes=A._boot()
    D,E,F=op.D,op.E,op.F
    runtime=R,op,L
    structural(A,L,F)
    ctx.prec=640
    E.core.A.cache_clear()
    synthetic_gain(L)
    ids=["optoGEF2_small_col_"+k for k in "CDEFGHIJKMNOPQR"]
    require([x["record_id"] for x in p["static_profiles"]]==ids,"exact remaining fifteen")
    seen=set()
    for s in p["static_profiles"]:
        record=s["record_id"]
        r=op.C.load_context(record)
        zeroes=[sum(r.target.a["reference"][int(n),0]==0 for n in tri[:3]) for tri in r.target.a["triangles"]]
        require({k:zeroes.count(j) for k,j in (("off_axis",0),("axis_vertex",1),("axis_edge_physical",2))}==s["routes"]
                and r.input_hashes==s["context_input_hashes"] and r.drive==r.dose*r.base["common_gain"]>0
                and r.dose==Q(s["exact_dose"]) and r.old_radius==Q(s["exact_original_radius"])
                and r.base["common_gain"]==Q(s["base_gain_exact"]) and r.row["parameter_cell_sha256"]==s["same_C_sha256"]
                and s["terminal_actual_packets_opened"] is (record in CD),"each isolated actual runtime row")
        seen.add(digest(r.u.tobytes()))
    require(len(seen)==15,"each record has its own anchor")
    blocked=AssertionError("unexpected new numerical work or L-only data loader")
    review={}
    with (patch.object(F,"_geometry",side_effect=blocked),patch.object(F,"_asset",side_effect=blocked),
          patch.object(L,"_source_inputs",side_effect=blocked),patch.object(E.M,"triangle_moment",side_effect=blocked),
          patch.object(E.U,"uniform_cache",side_effect=blocked),patch.object(E.P.Functional,"__init__",side_effect=blocked),
          patch.object(E.core,"bulk_local",side_effect=blocked)):
        for record,h in CD.items():
            b=bound(f"{stem}/{record}_proofs.json",h)
            dom,cases,routes,packets,prior=A.bind_reviewed_record(runtime,plan,jobs,expect,record,A.CD[record],hashes)
            require(b["record_id"]==record and b["actual_direct_review_sha256"]==A.CD[record]
                    and b["closed_cell_sha256"]==dom.context.row["parameter_cell_sha256"]
                    and b["retained_surface_manifest"]==packets["surface"] and b["retained_upper_Hessian_entries"]==6676
                    and b["retained_fullbox_b_exact"]==prior["b_upper_exact"]
                    and b["new_force_or_cache_sums_started"] is False,"retain signed surface, full H and derivative b")
            review[record]=profile(E,L,dom,cases["bulk_selectors"],routes["bulk"],packets["bulk"],b)
            for i in (0,6,72):
                man=packets["bulk"][i]; se=cases["bulk_selectors"][i]; rt=routes["bulk"][i]
                packet=bound(man["file"],man["sha256"])
                require(A.prepare_component(runtime,dom,se,rt,packet,man["sha256"],i)==b["bulk_proofs"][i],
                        "frozen per-record proof replay on each representative chart")
                mocked_component(A,runtime,dom,se,rt,man,b["bulk_proofs"][i])
            print(record,"144 frozen point/projection/tail proofs checked",flush=True)
    checks=dict.fromkeys(("frozen_sources_and_reviewed_transitive_runtime_aliases","all_15_static_contexts_and_distinct_anchors",
                          "completed_CD_only_and_preserved_leaves","exact_same_288_anchor_p2_flows_and_force_tail_norms",
                          "all_six_axial_derivatives_verified_by_independent_barycentric_chain_rule",
                          "affine_gain_dose_parity_fixed_columns_and_signed_state_MVT",
                          "six_mock_callbacks_with_current_full_Y_and_exact_base_drive",
                          "retain_whole_box_H_surface_b_and_exact_original_C",
                          "no_new_moments_force_H_or_cache_jobs"),True)
    report={"status":"PASS",**PINS,"classification":"PASS_CENTERED_COHORT_BULK_MODULE_PREPARATION_ONLY",
            "checks":checks,"review_source_sha256":digest(Path(__file__).read_bytes()),
            "proofs_by_record":CD,"record_chart_coverage":{r:{k:len(v) for k,v in x.items()} for r,x in review.items()},
            "shared_prior_component_engine_review":{"file":L_ACTUAL,"sha256":L_ACTUAL_SHA},
            "new_force_H_or_cache_jobs_started":False,"authorizes_actual_jobs_alone":False,
            "caller_requirements":["Set Arb precision 640 and clear core.A before constructing or binding a Domain.",
                "Use unchanged Domain returned from bind_reviewed_record and exact current packet hash/proof; centered preparation accepts generic Python Domain arguments and relies on this binding.",
                "Bind exact source closure, terminal review and per-record operator proof in a reviewed actual supervisor.",
                "Full corrector derivative uses original full-box H and b. Width or component success alone gives no signed residual or root certificate."],
            "scope":"Pure centered force on authenticated full original nodal box and closed C. Only C and D presently have all reviewed component prerequisites and frozen point proof packets."}
    raw=(json.dumps(report,indent=2,allow_nan=False)+"\n").encode()
    with (BATCH/"independent_review.json").open("xb") as fd:
        fd.write(raw)
    print(json.dumps({"status":report["classification"],"sha256":digest(raw),"checked_proofs":288}))


if __name__=="__main__":
    run()

