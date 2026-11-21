"""Independent cold review of the K/L whole-envelope geometry phase."""
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PILOT=HERE.parent/"exact_loading_tube_pilot_preparation_v1"
FILES={
    "source_sha256":("pilot.py","d39ae6c3bfa3d932eaba03ddd59fc275fd5b3b1c094e7f9a544810d8c10dc8bc"),
    "protocol_sha256":("protocol.json","fe64074d1a169a7672039c73dc07ca35c728cc1c299cf2a28025e0a0e7a97e43"),
    "preflight_sha256":("preflight.json","911f6839bd9cc49dafe626e688cc30d88a283ad87b5aaf92fecc5b629a684881"),
    "plan_sha256":("plan.json","c3b1870e5ad815bc5df5ec94d4731b0228e529ba2e0651f847102b8e0ecaa277"),
    "envelope_sha256":("envelope.json","f50c04b1802176abf4f5494898c8128018d4db8f8d0bb410af7774e5a3fff925"),
}
BINDINGS={k:h for k,(f,h) in FILES.items()}


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def read(p,h=None):
    b=p.read_bytes()
    assert h is None or sha256(b).hexdigest()==h
    return json.loads(b)


def reject(fn, errors=(KeyError,ValueError)):
    try:
        fn()
    except errors:
        return
    raise AssertionError("inert bad input was accepted")


def load_source():
    assert all(sha(PILOT/f)==h for f,h in FILES.values())
    spec=spec_from_file_location("_cold_KL_pilot",PILOT/"pilot.py")
    P=module_from_spec(spec)
    spec.loader.exec_module(P)
    return P


def envelope(P, plan, env):
    k,l=(read(ROOT/p,h) for p,h in zip(plan["anchor_box_files"], plan["anchor_box_sha256"],strict=True))
    assert k["reference_nodes_exact"]==l["reference_nodes_exact"]==env["reference_nodes_exact"]
    assert k["triangles"]==l["triangles"]==env["triangles"]
    assert k["free"]==l["free"]==env["free"] and len(env["free"])==611
    flatten=lambda key,d:[Q(i) for row in d[key] for i in row]
    ref=flatten("reference_nodes_exact",env)
    y0,y1=(flatten("current_centers_exact",d) for d in (k,l))
    mid=flatten("current_centers_exact",env)
    rad=flatten("coordinate_radii_exact",env)
    radius=Q(plan["fiber_radius_exact"])
    do0,do1=map(Q,plan["dose_interval_exact"])
    assert do1-do0==Q(plan["dose_gap_exact"])>0 and radius==Q(1,4096)
    assert len(ref)==len(y0)==len(y1)==len(mid)==len(rad)==650
    free=set(env["free"])
    assert 0<len(free)<650 and all(type(i) is int for i in free)
    for i,(a,b,c,h,x) in enumerate(zip(y0,y1,mid,rad,ref,strict=True)):
        if i in free:
            assert c-h==min(a,b)-radius and c+h==max(a,b)+radius
            # The extreme fibers enclose all affine intermediate fibers.
            for f in (Q(0),Q(1,3),Q(2,3),Q(1)):
                a_f=(1-f)*a+f*b
                assert c-h<=a_f-radius<=a_f+radius<=c+h
        else:
            assert a==b==c==x and h==0
    slopes=[(y1[i]-y0[i])/(do1-do0) for i in env["free"]]
    assert list(map(Q,plan["predictor_slope_exact"]))==slopes
    assert list(map(Q,plan["K_center_displacement_exact"]))==[y0[i]-ref[i] for i in env["free"]]
    assert list(map(Q,plan["L_center_displacement_exact"]))==[y1[i]-ref[i] for i in env["free"]]
    radial_floor=min(mid[i]-rad[i]-ref[i]-Q(float(-(1-1e-7)*float(ref[i]))) for i in env["free"] if i%2==0)
    assert radial_floor>0
    assert env["root_membership_assumed"] is env["any_spatial_derivative_bound_assumed"] is False
    assert plan["cached_or_exact_endpoint_H_used_on_envelope"] is plan["path_certified"] is False
    routes=read(ROOT/plan["K_routes_file"],plan["K_routes_sha256"])
    assert routes["positive_deformed_geometry_proved"] is False
    assert [k["element"] for k in routes["bulk"]]==list(range(144))
    assert [(k["face"],k["edge"]) for k in routes["surface"]]==[(i,j) for i in range(2) for j in range(12)]
    node_map={n:j for j,n in enumerate(env["free"])}
    active=set(plan["active_slots"])
    for row in routes["bulk"]+routes["surface"]:
        nodes=row["nodes"]
        local=[j for j in range(2*len(nodes)) if 2*nodes[j//2]+j%2 in free]
        assert row["free_local"]==local and row["global_indices"]==[node_map[2*nodes[j//2]+j%2] for j in local]
        if "element" in row:
            slot=row["element"]
            axis=sum(ref[2*n]==0 for n in nodes[:3])
            assert nodes==env["triangles"][slot] and row["chart"]==("off_axis","axis_vertex","axis_edge_physical")[axis]
            assert int(row["active"])==(slot in active)
            assert not ({"J_bound","J_bounds_exact","rho_exact"} & set(row))
        else:
            assert Q(row["Rlo_exact"])==ref[2*nodes[0]] and Q(row["Rhi_exact"])==ref[2*nodes[-1]]
            assert not ({"domain","surface_cache_line_geometry"} & set(row))
    return radial_floor,routes


def current_context_spy(P, plan, box):
    verified=[]
    def stopped(A,ctx,marked,actual):
        require=ctx.target.a
        assert isinstance(require,dict) and marked is plan and actual is box
        assert ctx.drive==Q(plan["nominal_drive_exact"])
        assert ctx.target.energy_scale==Q(plan["energy_scale_binary_exact"])
        assert A.a_sha(ctx.R)==plan["K_binding"]["R_binary_sha256"]
        assert A.a_sha(ctx.u)==plan["K_binding"]["center_binary_sha256"]
        assert A.X.arb_context.prec==640
        verified.append(True)
        return "INERT_CONTEXT_ONLY"
    with patch.object(P,"domain_components",side_effect=stopped):
        assert P.geometry_bounds(plan,box)=="INERT_CONTEXT_ONLY"
    assert verified==[True]


def routes_spy(P, routes, plan, box):
    # Only structural dimensions are shared. All mathematical return values are inert tokens.
    calls=[]
    nodes=object()
    failure={}
    indexed={(r["face"],r["edge"]):r for r in routes["surface"]}
    def whole(x,y,tri,eps):
        assert (x,y,tri,eps)==tuple(box[k] for k in ("reference_nodes_exact","current_centers_exact","triangles","coordinate_radii_exact"))
        calls.append("body")
        return {"status":"UNRESOLVED_GEOMETRY_BOX" if failure.get("body") else "CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING",
                "triangles":list(range(144))}
    def common(t,b):
        assert b is box
        calls.append("full_nodes")
        return nodes
    def field(t,Y,slot):
        assert Y is nodes
        r=routes["bulk"][slot]
        calls.append(("bulk",slot))
        return ({**r,"chart":{"kind":r["chart"]},"active":Q(r["active"])},r["nodes"])
    def rho(slot):
        if failure.get("bulk")==slot:
            raise ValueError("inert analytic failure")
        return (Q(1,2),Q(3,2),Q(1,2),Q(1,2))
    def edge(f,k,nodes):
        assert nodes is sentinel
        calls.append(("trace",f,k))
        return dict(indexed[f,k])
    def uniform(fields,max_depth):
        assert max_depth==8
        return {"status":"UNRESOLVED_UNIFORM_DOMAIN" if failure.get("trace")==
                (fields["face"],fields["edge"]) else "CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN"}
    def site(t,b,f,k,i):
        assert b is box
        calls.append(("site",f,k,i))
        if failure.get("site")== (f,k,i):
            raise ValueError("inert same-Y cache line failure")
        return {"status":"INERT_SITE"}
    sentinel=nodes
    fake=SimpleNamespace(X=SimpleNamespace(G=SimpleNamespace(certify_box=whole)),
                          E=SimpleNamespace(_bulk_fields=field,_body_rho=rho,
                          SUR=SimpleNamespace(common_nodes=common,uniform_G_bounds=uniform,local_site_domain=site)))
    ctx=SimpleNamespace(target=SimpleNamespace(surface=edge))
    with patch.object(P,"hashed_json",return_value=routes):
        z=P.domain_components(fake,ctx,plan,box)
        assert z[0]=="PASS_KL_FRESH_WHOLE_ENVELOPE_GEOMETRY_CANDIDATE" and len(z[2])==144 and len(z[3])==24
        assert calls.count("body")==calls.count("full_nodes")==1
        assert sum(type(i) is tuple and i[0]=="site" for i in calls)==288
        assert [x for x in calls if type(x) is tuple and x[0]=="bulk"]==[("bulk",i) for i in range(144)]
        for item, value in (("bulk",0),("trace",(1,4)),("site",(0,7,11))):
            failure.clear();failure[item]=value
            assert P.domain_components(fake,ctx,plan,box)[0]=="UNRESOLVED_KL_WHOLE_ENVELOPE_BULK_OR_TRACE_DOMAIN"
        failure.clear();failure["body"]=True;calls.clear()
        assert P.domain_components(fake,ctx,plan,box)[0]=="UNRESOLVED_KL_WHOLE_ENVELOPE_BODY_OR_RADIAL_BOUND"
        assert calls==["body"]




def lifecycle(P, b, plan, box):
    good={"status":"PASS_INDEPENDENT_KL_LOADING_GEOMETRY_PRELAUNCH","bindings":b,"pilot_id":P.PILOT_ID,
          "checks":dict.fromkeys(P.GEOMETRY_CHECKS,True),"geometry_only_authorized":True,
          "mechanical_kernels_authorized":False}
    P.gate_geometry(good,b)
    for key,bad in (("mechanical_kernels_authorized",True),("geometry_only_authorized",1),
                    ("bindings",{**b,"envelope_sha256":"0"*64}),("pilot_id","some_other_target")):
        reject(lambda:P.gate_geometry({**good,key:bad},b))
    for key in P.GEOMETRY_CHECKS:
        reject(lambda:P.gate_geometry({**good,"checks":{**good["checks"],key:1}},b))
    for outcome in ("pass","unresolved","exception"):
        with TemporaryDirectory(dir=HERE) as where:
            f=Path(where)
            rh=P.save(f/"independent_geometry_prelaunch.json",good)
            candidate={"pass":"PASS_KL_FRESH_WHOLE_ENVELOPE_GEOMETRY_CANDIDATE",
                       "unresolved":"UNRESOLVED_KL_WHOLE_ENVELOPE_BODY_OR_RADIAL_BOUND"}
            def geometry(p,e):
                assert p is plan and e is box
                if outcome=="exception":
                    raise ValueError("inert injected exception")
                return candidate[outcome],{"inert":"body"},[],[],Q(1)
            with patch.object(P,"HERE",f), patch.object(P,"prepared",return_value=(b,plan,box,None)), \
                 patch.object(P,"geometry_bounds",side_effect=geometry):
                reject(lambda:P.run_geometry("0"*64))
                assert not (f/"actual_geometry").exists()
                if outcome=="exception":
                    reject(lambda:P.run_geometry(rh))
                else:
                    P.run_geometry(rh)
                reject(lambda:P.run_geometry(rh),(FileExistsError,))
            out=f/"actual_geometry"
            t=read(out/"terminal.json")
            start=read(out/"start.json")
            assert t["bindings"]==start["bindings"]==b and type(t["pid"]) is int and t["pid"]==start["pid"]>0
            if outcome=="exception":
                assert t["status"]=="KL_GEOMETRY_PROCESS_FAILURE" and sha(out/"failure.json")==t["failure_sha256"]
                assert not (out/"result.json").exists()
            else:
                r=read(out/"result.json",t["result_sha256"])
                assert r["status"]==t["mathematical_status"]==candidate[outcome]
                assert t["status"]=="KL_GEOMETRY_STAGE_FINISHED" and sha(out/"body.json")==r["body_sha256"]
                assert r["root_branch_claimed"] is False and r["new_integrations"]==0


def main():
    P=load_source()
    b,p,e,protocol=P.prepared()
    assert b==BINDINGS and tuple(protocol["geometry_source_required_checks"])==P.GEOMETRY_CHECKS
    rebuilt,box,cross=P.cold_plan()
    assert json.loads(json.dumps(rebuilt,default=str))==p and json.loads(json.dumps(box,default=str))==e
    sources=protocol["geometry_source_tree_sha256"]
    assert len(sources)==34 and all(sha(ROOT/name)==h for name,h in sources.items())
    assert protocol["whole_KL_path_certified"] is False and protocol["geometry_or_mechanics_launch_authorized"] is False
    radial_floor,routes=envelope(P,p,e)
    current_context_spy(P,p,e)
    routes_spy(P,routes,p,e)
    lifecycle(P,b,p,e)
    assert not (PILOT/"actual_geometry").exists() and not (PILOT/"independent_geometry_prelaunch.json").exists()
    r={"status":"PASS_INDEPENDENT_KL_LOADING_GEOMETRY_COLD_FIXTURES","pilot_id":P.PILOT_ID,"bindings":b,
       "checks":dict.fromkeys(P.GEOMETRY_CHECKS,True),
       "source_tree_files":len(sources),"crosswalk_input_pins":cross,"whole_box_coordinate_checks":650,
       "whole_envelope_radial_trial_margin_exact":str(radial_floor),
       "route_checks":{"bulk":144,"trace_edges":24,"inert_same_Y_surface_sites":288},
       "all_actual_nontrivial_geometry_calls_forbidden_in_this_run":True,"geometry_kernels_run":0,
       "mechanical_kernels_run":0,"matched_K_real_context_spy_then_inert_all_168_routing":True,
       "fresh_routing_hold_fixtures":["inert_body_unresolved","inert_bulk_rho_invalid","inert_trace_invalid","inert_surface_line_invalid"],
       "write_once_fixtures":["review_before_output","linked_mathematical_unresolved","linked_runtime_exception","duplicate_run_rejected"],
       "scope":"A geometry-only release may authorize fresh mathematics; no root, branch, or mechanical proof was evaluated."}
    print(P.save(HERE/"audit_result.json",r),flush=True)


if __name__=="__main__":
    main()

