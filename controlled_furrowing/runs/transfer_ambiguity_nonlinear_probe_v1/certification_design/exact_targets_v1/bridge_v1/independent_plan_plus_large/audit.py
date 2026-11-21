"""Recompute the reviewed +large bridge plan with no spatial integrations."""
import argparse
from collections import Counter
from contextlib import ExitStack
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys
from unittest.mock import patch

SOURCE = Path(__file__).resolve()
OUT = SOURCE.parent
ROOT = SOURCE.parents[6]
sys.path.insert(0,str(ROOT/"scripts"))
import transfer_ambiguity_exact_bridge as J

PID = "max_contrast_change_large_03"
SOURCE_PIN = "58ae92f849414f185c4e28bbb2ec3b9d05361391c4d6e3fa71f519162e30e5e7"
GATE = "e9d4f429e467f10cb0d1c3d736fed5b6cd88a195c1873138a8b2d8d7730d416c"
WRAPPER = "269229312793ec101cc301086d634dd4de89c62a5050ccfc60f185cf761dbc7d"
SCOPES = (J.POINT,J.BOX)
E,S,T,C = J.E,J.SUR,J.T,J.C
ORDER = tuple((k,i) for k,n in (("bulk",144),("surface",24)) for i in range(n))


def h(path):
    return sha256(path.read_bytes()).hexdigest()


def digest(obj):
    return sha256(json.dumps(J.EX.doc(obj),sort_keys=True,separators=(",",":")).encode()).hexdigest()


def guard(force,scaled):
    stack = ExitStack()
    stop = AssertionError("force, Hessian, signed moment or AD entered plan audit")
    for mod,name in ((J,"build_bulk"),(J,"build_surface"),(J,"evaluate"),(J,"surface_selector"),
                     (J,"point_maps"),(J,"coefficientwise_inside"),(E,"bulk_component"),
                     (E,"surface_component"),(force,"_force_parts"),(scaled,"bulk_component"),
                     (J.U,"uniform_cache"),(S,"run_edge"),(C,"local_terms"),(C,"assemble_point"),
                     (E.P.Functional,"integral"),(E.P.Functional,"defect"),
                     (S.S.SurfaceFunctional,"integral"),(S.S.SurfaceFunctional,"defect"),
                     (scaled.S.Functional,"integral"),(scaled.S.Functional,"defect")):
        stack.enter_context(patch.object(mod,name,side_effect=stop))
    return stack


def pilot_chain():
    gate_path = J.PREP/"prelaunch_independent_review.json"
    assert h(J.SOURCE)==SOURCE_PIN and h(gate_path)==GATE
    state,b = J.plan_gate(gate_path,GATE)
    assert J.DERIV_ADMISSION_SHA==WRAPPER
    proof = J.derivative_inputs(J.DERIV_ADMISSION,WRAPPER,state[1])[PID]
    links = {**b,"derivative_admission_sha256":WRAPPER,
             "independent_native_derivative_review_sha256":proof["independent_native_sha256"]}
    tail = tuple(Q(proof["planning"][key]) for key in ("per_component_force_entry_tail_cap_exact",
                                                        "per_component_H_entry_tail_cap_exact"))
    assert all(x>0 for x in tail)
    return state,proof,links,tail


def canonical(which, fields,jet):
    def encode(p):
        return [[i,j,list(map(str,S.rational_hull(v)))] for (i,j),v in sorted(p.items())]
    q,B = (fields[k] for k in which)
    return digest({"q":[encode(x) for x in q],"B":[[encode(x) for x in row] for row in B],
                   "J_or_G":encode(jet["value"]),"first":[encode(x) for x in jet["gradient"]],
                   "second":[[i,j,encode(x)] for (i,j),x in sorted(jet["hessian_upper"].items())]})


def inside(which,point,box):
    def flat(f,j):
        q,b = [f[k] for k in which]
        return [*q,*(p for row in b for p in row),j["value"],*j["gradient"],
                *[j["hessian_upper"][k] for k in sorted(j["hessian_upper"])]]
    assert point[1]["hessian_upper"].keys()==box[1]["hessian_upper"].keys()
    count = 0
    for a,b in zip(flat(*point),flat(*box),strict=True):
        for k in a.keys()|b.keys():
            q0,q1 = S.rational_hull(a.get(k,0))
            lo,hi = S.rational_hull(b.get(k,0))
            assert lo <= q0 <= q1 <= hi
            count += 1
    return count


def bulk_orders(identity, body, mass, grad, hess, tail):
    l,u = [Q(body["J_bound"][k]) for k in ("lower_exact","upper_exact")]
    assert body["positive_J_for_entire_box"] is True and Q(1,10**6)<l<=u
    raw = max(abs(l-1),abs(u-1))
    den = 1<<48
    rho = Q(-((-raw.numerator*den)//raw.denominator),den)
    assert raw<=rho<1
    g,d = max(grad),max(hess.values())
    attempts = []
    for n in range(2,33):
        a,b = E.M.volumetric_derivative_tail_bounds(n,rho,1,5000)
        f = mass*a*g
        h = mass*(b*g*g+a*d)
        attempts.append({"degree":n,"phi_prime_tail_upper_exact":str(a),
                         "phi_second_tail_upper_exact":str(b),"force_tail_upper_exact":str(f),
                         "H_entry_tail_upper_exact":str(h)})
        if f<=tail[0] and h<=tail[1]:
            return {**identity,**attempts[-1],"status":"SELECTED_UNIFORM_BULK_ORDER","J_bounds_exact":[str(l),str(u)],
                    "raw_rho_exact":str(raw),"rho_exact":str(rho),"force_tail_budget_exact":str(tail[0]),
                    "hessian_tail_budget_exact":str(tail[1]),"geometry_two_mass_upper_exact":str(mass),
                    "max_DJ_norm_upper_exact":str(g),"max_D2J_norm_upper_exact":str(d),"max_order":32,
                    "order_attempts":attempts}
    raise AssertionError("advertised selected bulk order is unresolved")


def site_map_points(d, index, ids, fields):
    wrapped = J.U.RuleTarget(d["ctx"].target,d["rule"])
    maps = [J.U.C.exact_bulk_maps(wrapped,index,k) for k in range(64)]
    y = [T.rational(c) for i in ids for c in d["point"][i]]
    lower = []
    for pos,m in enumerate(maps):
        assert m["nodes"] == ids and m["global_indices"] == fields["global_indices"]
        assert m["selected_local"] == fields["free_local"]
        assert Q(m["auxiliary_normalized_weight_exact"]) == d["rule"]["volume_weights"][index][pos] > 0
        coeff = list(zip(m["star_full"],m["cache_full"],strict=True))
        states = []
        for s,c in coeff:
            ends = [sum(d*v for d,v in zip(row,y,strict=True)) for row in (c,s)]
            states.append(tuple(sorted(ends)))
        det = min(a*b for a,b in product(states[0],states[3]))-max(a*b for a,b in product(states[1],states[2]))
        q = det*states[4][0]
        assert det > 0 and states[4][0] > 0 and q > Q(1,10**6)
        lower.append(q)
    return {"map_packets_sha256":digest(maps),"positive_line_count":64,
            "minimum_exact_convex_line_J_lower":str(min(lower))}


def bulk(d,index,p,tail,scaled):
    body = d["body"]["triangles"][index]
    route = d["routes"]["bulk"][index]
    assert p["kind"]=="bulk" and type(p["slot"]) is int and p["slot"]==index
    assert p["body_element_sha256"]==digest(body) and p["global_indices"]==route["global_indices"]
    assert p["point_J_guard_scope"]=="Exact point is inside the independently certified same-element box."
    slots = []
    for scope,nodes in zip(SCOPES,(d["point"],d["nodes"]),strict=True):
        C.A.cache_clear()
        # Shared reviewed hydrostatic strain/jet constructors and norm bounds.
        ident,f,_,jet,mass,g,H = E._bulk_setup(d["ctx"],nodes,d["rule"],index)
        assert all(ident[k]==route[k] for k in ("chart","active","global_indices","free_local","auxiliary_element_rule_sha256"))
        assert f["global_indices"]==route["global_indices"] and all(type(j) is int for j in route["global_indices"])
        sel = bulk_orders(ident,body,mass,g,H,tail)
        assert p[scope]["selector"]==sel
        assert p[scope]["polynomial_jets_sha256"]==canonical(("q","B"),f,jet)
        if scope==J.POINT:
            assert p[scope]["same_Y_point_cache_lines"]==site_map_points(d,index,route["nodes"],f)
        else:
            # Existing record-neutral scaled polynomial route proof is a pure helper.
            replay = scaled.prepare_component(d["ctx"],nodes,d["rule"],body,sel,index)
            assert J.EX.doc(replay)==p[scope]["scaled_or_axis_proof"]
            assert replay["record_id"]==d["ctx"].record_id and replay["local_nodes"]==route["nodes"]
            if route["chart"]=="off_axis":
                assert replay["integration_chart"]=="scaled_positive_R_affine_material"
                assert replay["scaled_exact_Qstar_node_count"]==64
                assert Q(replay["scaled_force_tail_max_exact"])<=tail[0]
                assert Q(replay["scaled_H_entry_tail_max_exact"])<=tail[1]
            else:
                assert replay["integration_chart"]==route["chart"] and replay["scaled_polynomial_identity_count"]==0
        slots.append((f,jet))
    count = inside(("q","B"),*slots)
    assert count==p["exact_point_coefficients_enclosed_by_whole_box_count"]
    return count


def trace_orders(d,fields,domain,fi,ei,jet,tail):
    a,b = map(Q,domain["G_bounds_exact"])
    q0 = (a+b)/2
    raw = (b-a)/(b+a)
    assert domain["status"]=="CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN"
    assert 0<a<=b and Q(domain["q0_exact"])==q0 and Q(domain["rho_exact"])==raw
    den=1<<48
    rho=Q(-((-raw.numerator*den)//raw.denominator),den)
    assert 0<=raw<=rho<1
    bounded={**domain,"raw_rho_exact":str(raw),"rho_exact":str(rho)}
    _,_,K,Tn = S.material_boxes(d["ctx"].row,fi)
    coeff=max(abs(t-k) for t,k in product(Tn,K))
    x,y=fields["Rlo"],fields["Rhi"]
    root=Q(str(C.A(q0).sqrt().lower().fmpq()))
    assert root>0 and root**2<=q0
    mass=Q(str(C.upper(d["ctx"].target.normalization()*C.A((y*y-x*x)/2)
                     +C.A(sum(d["rule"]["surface_weights"][fi][ei],Q(0)))).fmpq()))
    g=max(S.poly_norm_upper(v,x,y) for v in jet["gradient"])
    H=max(S.poly_norm_upper(v,x,y) for v in jet["hessian_upper"].values())
    attempts=[]
    for n in range(2,65):
        p=coeff*S.S.binomial_tail(Q(-1,2),n,rho)/(2*root)
        q=coeff*S.S.binomial_tail(Q(-3,2),n,rho)/(4*q0*root)
        f=mass*p*g
        h=mass*(p*H+q*g*g)
        attempts.append({"order":n,"force_tail_upper_exact":str(f),"H_entry_tail_upper_exact":str(h)})
        if f<=tail[0] and h<=tail[1]:
            return n,attempts,bounded
    raise AssertionError("advertised selected trace order is unresolved")


def trace(d,slot,p,tail):
    fi,ei=divmod(slot,12)
    route=d["routes"]["surface"][slot]
    geom=d["geometry"]["trace_geometry_only"][slot]
    assert p["kind"]=="surface" and type(p["slot"]) is int and p["slot"]==slot
    assert (route["face"],route["edge"])==(geom["face"],geom["edge"])==(fi,ei)
    assert p["global_indices"]==route["global_indices"]==geom["global_indices"]
    keys=("q_physical_R","B_physical_R")
    checks=[]
    for scope,nodes,box in zip(SCOPES,(d["point"],d["nodes"]),(d["point_box"],d["box"]),strict=True):
        C.A.cache_clear()
        target=d["ctx"].target
        fields=target.surface(fi,ei,nodes=nodes)
        assert fields["free_local"]==route["free_local"] and fields["global_indices"]==route["global_indices"]
        E._surface_rule(target,d["rule"],fi,ei)
        assert digest([d["rule"]["surface_s"][fi],d["rule"]["surface_weights"][fi][ei]])==route["trace_rule_sha256"]
        domain=S.uniform_G_bounds(fields,max_depth=8)
        sites=[S.local_site_domain(target,box,fi,ei,k) for k in range(12)]
        for k,s in enumerate(sites):
            assert s["quadrature_index"]==k and s["nodes"]==route["nodes"] and s["selected_local"]==route["free_local"]
            assert s["global_indices"]==route["global_indices"] and Q(s["tangent_norm_squared_lower_exact"])>0
            assert Q(s["auxiliary_normalized_weight_exact"])==d["rule"]["surface_weights"][fi][ei][k]>0
        if scope==J.BOX:
            assert domain==geom["domain"] and sites==geom["cache_site_geometry"]
        jet=E.JET.surface_squared_jacobian_jet(fields[keys[0]],fields[keys[1]])
        n,attempts,bounded=trace_orders(d,fields,domain,fi,ei,jet,tail)
        rec={"record_id":d["ctx"].record_id,"face":fi,"edge":ei,"global_indices":route["global_indices"],
             "order":n,"domain":bounded,"cache_site_domains":sites,"order_attempts":attempts,
             "polynomial_jets_sha256":canonical(keys,fields,jet),"node_ball_sha256":E._local_node_digest(nodes,route["nodes"])}
        assert rec==p[scope]
        checks.append((fields,jet))
    count=inside(keys,*checks)
    assert count==p["exact_point_coefficients_enclosed_by_whole_box_count"]
    return count


def prepare():
    _,_,links,tail=pilot_chain()
    assert len(ORDER)==168 and len(set(ORDER))==168
    assert all(x>0 for x in tail)
    # The positive-control coefficient inclusion with a missing exact zero is sound.
    a={"q":[{(0,0):Q(1)}],"B":[[{}]]}
    b={"q":[{(0,0):Q(1)}],"B":[[{(0,0):C.A(0).union(C.A(1))}]]}
    j={"value":{},"gradient":[{}],"hessian_upper":{(0,0):{}}}
    assert inside(("q","B"),(a,j),(b,j))==2
    assert not (OUT/"review.json").exists()
    p={"status":"FROZEN_INDEPENDENT_PLUS_LARGE_PLAN_READBACK","source_sha256":h(SOURCE),
       "producer_bindings":links,"pilot_id":PID,"tail_caps_exact":list(map(str,tail)),
       "requires_externally_observed_real_plan_exit":True,"actual_plan_payloads_opened":False}
    q=J.save(OUT/"protocol.json",p)
    k=J.save(OUT/"preflight.json",{"status":"PASS_PREPARATION_ONLY","source_sha256":h(SOURCE),"protocol_sha256":q,
        "checks":{"all_168_component_keys_and_explicit_pilot":True,"reviewed_frozen_authority_rebound":True,
                  "coefficient_interval_synthetic_containment":True,"no_plan_payload_or_spatial_integrations":True}})
    print(json.dumps({"source_sha256":h(SOURCE),"protocol_sha256":q,"preflight_sha256":k}))


def audit(args):
    pre=J.pinned(OUT/"preflight.json",args.preflight_sha)
    proto=J.pinned(OUT/"protocol.json",args.protocol_sha)
    assert pre["status"]=="PASS_PREPARATION_ONLY" and J.sweep.all_true(pre["checks"])
    assert pre["protocol_sha256"]==args.protocol_sha and pre["source_sha256"]==proto["source_sha256"]==h(SOURCE)
    assert proto["pilot_id"]==PID and proto["actual_plan_payloads_opened"] is False
    assert args.os_exit_code==0 and args.session==17581
    state,deriv,links,tail=pilot_chain()
    assert proto["producer_bindings"]==links and proto["tail_caps_exact"]==list(map(str,tail))
    folder=J.PLAN/PID
    term=J.pinned(folder/"terminal.json",args.terminal_sha)
    started=J.sweep.read(folder/"started.json")
    assert term["status"]=="FINISHED" and type(term["exit_code"]) is int and term["exit_code"]==0
    assert term["pilot_id"]==started["pilot_id"]==PID and term["pid"]==started["pid"]==args.pid
    assert started["bindings"]==links and term["plan_sha256"]==args.plan_sha
    p=J.pinned(folder/"plan.json",args.plan_sha)
    assert p["status"]=="COMPLETE_POINT_FORCE_WHOLE_H_SELECTOR_PLAN_UNREVIEWED" and p["pilot_id"]==PID
    assert p["bindings"]==links and p["tail_caps_exact"]==list(map(str,tail)) and p["derivative"]==deriv
    assert p["point_membership"]==J.POINT and p["H_membership"]==J.BOX
    assert Q(p["whole_box_radius_exact"])==J.RADIUS and p["spatial_component_launch_authorized"] is False
    assert tuple((e["kind"],e["slot"]) for e in p["components"])==ORDER
    for e in p["components"]:
        assert type(e["slot"]) is int and e["file"]==f"{e['kind']}_{e['slot']:03d}.json"
        assert h(folder/e["file"])==e["sha256"]

    force,scaled=J.kernels({})
    lengths=Counter()
    orders={f"{kind}/{what}":Counter() for kind in ("bulk","surface") for what in ("point","box")}
    summary=[]
    with guard(force,scaled):
        d=J.domain(state,PID)
        assert p["context_box_sha256"]==d["entry"]["box_1_sha256"] and p["point_box_sha256"]==digest(d["point_box"])
        assert p["record_id"]==d["ctx"].record_id and C.ctx.prec==640
        for e in p["components"]:
            i=e["slot"]
            proof=J.pinned(folder/e["file"],e["sha256"])
            if e["kind"]=="bulk":
                count=bulk(d,i,proof,tail,scaled)
                info=proof[J.BOX]["scaled_or_axis_proof"]
                lengths[info["integration_chart"]]+=1
                lengths["off_axis_exact_identities"]+=info["scaled_polynomial_identity_count"]
                orders["bulk/point"][proof[J.POINT]["selector"]["degree"]]+=1
                orders["bulk/box"][proof[J.BOX]["selector"]["degree"]]+=1
            else:
                count=trace(d,i,proof,tail)
                lengths["surface"]+=1
                orders["surface/point"][proof[J.POINT]["order"]]+=1
                orders["surface/box"][proof[J.BOX]["order"]]+=1
            lengths["coefficient_inclusion_checks"]+=count
            summary.append({"kind":e["kind"],"slot":i,"plan_leaf_sha256":e["sha256"],"coefficient_checks":count})
            if len(summary) in (36,72,108,144,168):
                print(json.dumps({"plan_components_recomputed":len(summary),"total":168}),flush=True)

    assert lengths["scaled_positive_R_affine_material"]==132
    assert lengths["axis_vertex"]==lengths["axis_edge_physical"]==6 and lengths["surface"]==24
    assert lengths["off_axis_exact_identities"]>0 and lengths["coefficient_inclusion_checks"]>0
    assert len(summary)==168 and not (J.ACTUAL/PID).exists()
    out={"status":"PASS_INDEPENDENT_ACTUAL_SINGLETON_BRIDGE_PLAN","pilot_id":PID,"source_sha256":SOURCE_PIN,
         "plan_sha256":args.plan_sha,"plan_terminal_sha256":args.terminal_sha,
         "plan_observed_os_exit_code":args.os_exit_code,"plan_observed_pid":args.pid,"plan_observed_session":args.session,
         "derivative_admission_sha256":WRAPPER,"audit_source_sha256":h(SOURCE),
         "audit_protocol_sha256":args.protocol_sha,"audit_preflight_sha256":args.preflight_sha,
         "checks":{"root_observed_plan_OS_and_all_leaf_pins":True,"all_168_exact_route_orders_and_tail_budgets_rebuilt":True,
                   "point_and_whole_box_polynomial_and_nodal_fingerprints_recomputed":True,
                   "all_64_point_bulk_and_12_each_trace_sites_matched":True,
                   "whole_box_scaled_and_axis_proofs_recomputed":True,
                   "no_force_H_moment_or_AD_jobs":True},
         "routes_and_checks":dict(lengths),"selected_order_counts":{k:dict(v) for k,v in orders.items()},
         "point_same_Y_bulk_cache_site_count":144*64,"point_and_box_surface_site_count":2*24*12,
         "component_readback":summary,"scope":"Independent selector arithmetic using shared reviewed geometry, sparse polynomial/jet/chart and scalar tail primitives. Whole-box scaled chart identities replay through its reviewed pure preparer. No force, Hessian, moments, AD, root, or equilibrium inference."}
    h0=J.save(OUT/"review.json",out)
    print(json.dumps({"review_sha256":h0,"route_counts":out["routes_and_checks"],"orders":out["selected_order_counts"]}))


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    sub=cli.add_subparsers(dest="stage",required=True)
    sub.add_parser("prepare")
    actual=sub.add_parser("review")
    for n in ("plan-sha","terminal-sha","protocol-sha","preflight-sha"):
        actual.add_argument("--"+n,required=True)
    for n in ("pid","session","os-exit-code"):
        actual.add_argument("--"+n,required=True,type=int)
    args=cli.parse_args()
    if args.stage=="prepare":
        prepare()
    else:
        audit(args)
