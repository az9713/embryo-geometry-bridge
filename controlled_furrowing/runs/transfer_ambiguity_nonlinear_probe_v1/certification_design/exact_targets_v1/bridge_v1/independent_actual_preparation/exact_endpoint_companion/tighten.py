"""Exact tail-radius replay after a genuinely completed component replay."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
import traceback

SOURCE=Path(__file__).resolve()
HERE=SOURCE.parent
PREP=HERE.parent
BRIDGE=PREP.parent
DESIGN=BRIDGE.parents[1]
ROOT=DESIGN.parents[2]
R=Q(1,2**34)
RAW=PREP/"actual_replay"
OUTPUT=HERE/"actual"
PINS={
    PREP/"replay.py":"e504a67f913968888805c0385c6b3fa0501d6b3e7f5345362b453eb1dde50074",
    PREP/"replay_math.py":"e7bb66066d54e131fd009d2bf48e94f3f1fb7bc1955023f66bc96cbe3c6d78e9",
    PREP/"protocol.json":"2ec903e8088596b72d3f59e57ce6700247218b2138fbe9e78edeb1cb046cb9a1",
    PREP/"preflight.json":"effd86b721146884ff4561e6bd24eabe3f32a260d768fd09218ef0432eb5442b",
}
SCOPE=("Tighten already computed independent local enclosures only by removing excess machine "
       "tail inflation and adding the independently justified exact rational tail. Whole-box/point "
       "quantifiers unchanged, actual strict endpoint containment required. No force, Hessian moment, "
       "AD/cache sums, root, stability, continuum or experimental conclusions.")


def name(p):
    return p.relative_to(ROOT).as_posix()


def sha(f):
    return sha256(f.read_bytes()).hexdigest()


def read(f,h):
    raw=f.read_bytes()
    assert sha256(raw).hexdigest()==h,f
    return json.loads(raw)


def save(path,x):
    with path.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(x,f,indent=2,allow_nan=False)
        f.write("\n")
    return sha(path)


def add_back(original, used_radius_exact, sufficient_radius_exact):
    lo,hi=original
    s,t=used_radius_exact,sufficient_radius_exact
    assert 0<=t<=s and lo<=hi and lo+s<=hi-s
    return lo+s-t,hi-s+t


def exact_stored_radius(A, lift):
    t=A.ends(lift)[1]
    assert t>=0
    raw=A.arb(0,A.E.core.upper(lift))
    assert raw.rad().is_exact() and raw.mid().is_zero()
    s=Q(str(raw.rad().fmpq()))
    assert A.ends(raw)==(-s,s) and s>=t
    return s


def environment():
    sys.path.insert(0,str(ROOT/"scripts"))
    sys.path.insert(0,str(PREP))
    import replay_math as A
    import transfer_ambiguity_exact_targets as EX
    import transfer_ambiguity_cached_derivative_extension as DX
    from exact_scaled_cohort_bulk import S
    return A,EX,DX,S


def phase():
    assert all(sha(p)==h for p,h in PINS.items())
    return {"status":"FROZEN_EXACT_ENDPOINT_TAIL_REFINEMENT_REVIEW","source_sha256":sha(SOURCE),
        "input_hashes":{name(f):h for f,h in PINS.items()},"radius_exact":str(R),"precision_bits":640,
        "bound_method":"If J encloses P + [-s,s], then [Jlo+s,Jhi-s] encloses P. "
            "Reinflate by the exact sufficient remainder t <= s. Verify actual new strict containment.",
        "phase_authority":"Sealed raw e504 true terminal, reviewed producer plan and source plus separate source readback before prospective actual.",
        "scope":SCOPE}


def prepare():
    p=phase()
    assert not OUTPUT.exists()
    from flint import arb,ctx,fmpq
    ctx.prec=640
    t=Q(1,10)
    b=arb(0,arb(fmpq(t.numerator,t.denominator)).abs_upper())
    assert b.rad().is_exact()
    s=Q(str(b.rad().fmpq()))
    assert s>t
    f=Q(1,3)
    old=(f-s,f+s)
    new=add_back(old,s,t)
    assert new==(f-t,f+t) and old[0]<new[0]<new[1]<old[1]
    h=save(HERE/"protocol.json",p)
    g=save(HERE/"preflight.json",{"status":"PASS","source_sha256":sha(SOURCE),"protocol_sha256":h,
         "exact_dyadic_stored_radius_and_erosion_fixture":True,
         "actual_raw_component_packets_opened":0,"force_H_cache_or_moment_work":0})
    print(h,g,sha(SOURCE),flush=True)


def bulk_tails(A,S,target,rule,nodes,index,plan,derivative):
    fields,ids=A.E._bulk_fields(target,nodes,index)
    assert fields["global_indices"]==plan["global_indices"]
    degree=plan["degree"]
    a,b=A.E.M.volumetric_derivative_tail_bounds(degree,Q(plan["rho_exact"]),1,5000)
    assert (a,b)==(Q(plan["phi_prime_tail_upper_exact"]),Q(plan["phi_second_tail_upper_exact"]))
    mass=A.E._geometry_mass_upper(target,fields,index,rule)
    assert mass==Q(plan["geometry_two_mass_upper_exact"])
    if derivative and fields["chart"]["kind"]=="off_axis":
        geo=S.chart(fields["corners"])
        q,B=S.numerators(geo,[nodes[j] for j in ids],fields["free_local"])
        jet=A.E.JET.bulk_jacobian_jet(q,B)
        g0,h0=(Q(plan[k]) for k in ("max_DJ_norm_upper_exact","max_D2J_norm_upper_exact"))
        g=[min(g0,A.norm(p)/geo.Rhat_min) for p in jet["gradient"]]
        h={k:min(h0,A.norm(p)/geo.Rhat_min) for k,p in jet["hessian_upper"].items()}
    else:
        jet=A.E.JET.bulk_jacobian_jet(fields["q"],fields["B"])
        g=[A.E._norm(p,fields) for p in jet["gradient"]]
        h={k:A.E._norm(p,fields) for k,p in jet["hessian_upper"].items()}
    if derivative:
        return {k:mass*(b*g[k[0]]*g[k[1]]+a*v) for k,v in h.items()}
    return [mass*a*x for x in g]


def surface_tails(A,target,rule,nodes,index,case,params,derivative):
    K,Tn=map(A.A,params)
    f,e=divmod(index,12)
    fields=target.surface(f,e,nodes=nodes)
    assert fields["global_indices"]==case["global_indices"]
    lo,hi=fields["Rlo"],fields["Rhi"]
    material=target.normalization()*A.A((hi*hi-lo*lo)/2)
    auxiliary=A.A(sum(rule["surface_weights"][f][e],Q(0)))
    mass=A.E.core.upper(material+auxiliary)
    q0=Q(case["domain"]["q0_exact"])
    rho=Q(case["domain"]["rho_exact"])
    order=case["order"]
    root=A.A(q0).sqrt()
    D=max(map(abs,A.SUR.rational_hull(Tn-K)))
    ea=A.A(D*A.SUR.S.binomial_tail(Q(-1,2),order,rho))/(2*root)
    eb=A.A(D*A.SUR.S.binomial_tail(Q(-3,2),order,rho))/(4*A.A(q0)*root)
    J=A.E.JET.surface_squared_jacobian_jet(fields["q_physical_R"],fields["B_physical_R"])
    g=[A.SUR.poly_norm_upper(p,lo,hi) for p in J["gradient"]]
    if derivative:
        return {k:mass*(ea*A.A(A.SUR.poly_norm_upper(p,lo,hi))+eb*A.A(g[k[0]]*g[k[1]]))
                for k,p in J["hessian_upper"].items()}
    return [mass*ea*A.A(v) for v in g]


def tighten(A,v,tail,cap):
    # For an exact bulk tail, avoid rounding the mathematical remainder t.
    # The floating radius s matches the old inflation operation. For the
    # surface tail take an exact upper endpoint of its outward Arb enclosure.
    lift=A.A(tail) if type(tail) is Q else tail
    s=exact_stored_radius(A,lift)
    t=tail if type(tail) is Q else A.ends(tail)[1]
    assert t<=cap and s>=t>=0
    return add_back(v,s,t),s,t


def verify_terminal(args):
    p=RAW/args.pilot_id
    tr=read(p/"terminal.json",args.raw_terminal_sha256)
    o=read(args.raw_os_exit,args.raw_os_exit_sha256)
    rec=read(p/"receipt.json",args.raw_receipt_sha256)
    assert rec["pilot_id"]==args.pilot_id and len(rec["components"])==168
    assert tr["receipt_sha256"]==o["receipt_sha256"]==args.raw_receipt_sha256
    assert o["terminal_sha256"]==args.raw_terminal_sha256
    assert type(o["observed_exit_code"]) is int and o["observed_exit_code"]==tr["exit_code"]
    assert tr["pid"]==rec["pid"]==o["pid"] and o["unified_exec_session"]>0
    assert (tr["status"],rec["status"],o["observed_exit_code"]) in (
        ("FINISHED","PASS_INDEPENDENT_ACTUAL_SINGLETON_BRIDGE_COMPONENT_REPLAY",0),
        ("FINISHED_WITH_HOLD","HOLD_INDEPENDENT_COMPONENT_CONTAINMENT",2))
    assert all(rec[k]==tr[k]==PINS[p] for k,p in (
        ("source_sha256",PREP/"replay.py"),("protocol_sha256",PREP/"protocol.json"),
        ("preflight_sha256",PREP/"preflight.json")))
    for f,h in (rec["checked_packet_sha256"]|rec["plan_packet_sha256"]).items():
        assert sha(ROOT/f)==h
    return rec


def compute(args,raw,out):
    A,EX,DX,S=environment()
    _,cases,manifest,data=DX.authority()
    case,=[c for c in cases if c["pilot_id"]==args.pilot_id]
    ctx,rule,entry=EX.open_context(args.pilot_id,manifest,data)
    stored=read(ROOT/case["box"],case["box_sha256"])
    assert stored==EX.exact_node_box(ctx,R)
    assert raw["box_sha256"]==entry["box_1_sha256"]==case["box_sha256"]
    assert raw["cache_packet_sha256"]==case["cache_packet_sha256"]
    assert raw["record_id"]==ctx.record_id and EX.array_sha(ctx.R)==raw["binary_preconditioner_sha256"]
    plan=read(BRIDGE/"plans"/args.pilot_id/"plan.json",raw["plan_sha256"])
    assert plan["context_box_sha256"]==case["box_sha256"] and Q(raw["radius_exact"])==Q(plan["whole_box_radius_exact"])==R
    caps=list(map(Q,plan["tail_caps_exact"]))
    point=ctx.target.nodes()
    nodes=A.SUR.common_nodes(ctx.target,stored)
    rows=[]
    remaining=[]
    fixed_counts={"point_force":0,"whole_box_H":0}
    for item,selector in zip(raw["components"],plan["components"],strict=True):
        k,i=item["kind"],item["slot"]
        assert type(i) is int and (k,i)==(selector["kind"],selector["slot"])
        d=read(ROOT/item["replay_file"],item["replay_sha256"])
        packet=read(ROOT/item["packet_file"],item["packet_sha256"])
        f=BRIDGE/"plans"/args.pilot_id/selector["file"]
        p=read(f,selector["sha256"])
        assert item["plan_file_sha256"]==selector["sha256"] and item["packet_sha256"]==raw["checked_packet_sha256"][item["packet_file"]]
        assert d["global_indices"]==packet["global_indices"]==p["global_indices"]
        assert all(d[n]==item[n] for n in ("point_force_scope","whole_box_H_scope","plan_file_sha256",
                                         "packet_file","packet_sha256"))
        assert d["point_force_scope"]==A.POINT and d["whole_box_H_scope"]==A.BOX
        assert d["containment_failures"]["direct_point_cache"]==[]
        qold=list(map(A.interval,d["independent_point_F_exact"]))
        hold={(j,n):A.interval((a,b)) for j,n,a,b in d["independent_whole_box_H_exact"]}
        supplied=list(map(A.interval,packet["point_force"]["F_I_minus_F_cache_intervals_exact"]))
        hs={(j,n):A.interval((a,b)) for j,n,a,b in packet["whole_box_hessian"]["H_I_minus_H_cache_upper_exact"]}
        assert len(qold)==len(supplied) and hold.keys()==hs.keys()
        f_bad=[j for j,(a,b) in enumerate(zip(supplied,qold,strict=True)) if not A.contain(a,b)]
        h_bad=[k0 for k0,h in hold.items() if not A.contain(hs[k0],h)]
        assert f_bad==d["containment_failures"]["point_force"]
        assert list(map(list,h_bad))==d["containment_failures"]["whole_box_H"]
        inp=(ctx.target,rule)
        notes={"point_original_excess_rounded_tails":[],"H_original_excess_rounded_tails":[]}
        if f_bad:
            if k=="bulk":
                tf=bulk_tails(A,S,*inp,point,i,p[A.POINT]["selector"],False)
            else:
                face=("basal","apical")[i//12]
                tf=surface_tails(A,*inp,point,i,p[A.POINT],
                     (ctx.base[face+"_area_modulus"],ctx.base[face+"_surface_tension"]),False)
            for j in f_bad:
                qold[j],s,t=tighten(A,qold[j],tf[j],caps[0])
                given=Q(packet["point_force"]["force_tail_upper_exact"][j])
                assert t<=given
                if k=="bulk":
                    assert t==given
                passed=A.contain(supplied[j],qold[j])
                fixed_counts["point_force"]+=int(passed)
                notes["point_original_excess_rounded_tails"].append([j,str(s),str(t),passed])
        if h_bad:
            if k=="bulk":
                th=bulk_tails(A,S,*inp,nodes,i,p[A.BOX]["selector"],True)
            else:
                face=("basal","apical")[i//12]
                th=surface_tails(A,*inp,nodes,i,p[A.BOX],
                     (ctx.base[face+"_area_modulus"],ctx.base[face+"_surface_tension"]),True)
            for j in h_bad:
                hold[j],s,t=tighten(A,hold[j],th[j],caps[1])
                passed=A.contain(hs[j],hold[j])
                fixed_counts["whole_box_H"]+=int(passed)
                notes["H_original_excess_rounded_tails"].append([list(j),str(s),str(t),passed])
        check_f=[j for j,(a,b) in enumerate(zip(supplied,qold,strict=True)) if not A.contain(a,b)]
        check_h=[list(j) for j,b in hold.items() if not A.contain(hs[j],b)]
        if check_f or check_h:
            remaining.append({"kind":k,"slot":i,"point_force":check_f,"whole_box_H":check_h})
        raw_file=item["replay_file"]
        raw_sha=item["replay_sha256"]
        row={**d,"raw_component_replay_file":raw_file,"raw_component_replay_sha256":raw_sha,
            "containment_failures":{"direct_point_cache":[],"point_force":check_f,"whole_box_H":check_h},
            "independent_point_F_exact":[list(map(str,r)) for r in qold],
            "independent_whole_box_H_exact":[[*j,*map(str,r)] for j,r in sorted(hold.items())],**notes}
        file=out/f"{k}_{i:03d}.json"
        h=save(file,row)
        rows.append(item|{"replay_file":name(file),"replay_sha256":h})
    return rows,fixed_counts,remaining


def run(args):
    p=read(HERE/"protocol.json",args.protocol_sha256)
    pf=read(HERE/"preflight.json",args.preflight_sha256)
    assert p==phase() and pf["source_sha256"]==sha(SOURCE) and pf["protocol_sha256"]==args.protocol_sha256
    assert pf["status"]=="PASS" and pf["force_H_cache_or_moment_work"]==pf["actual_raw_component_packets_opened"]==0
    raw=verify_terminal(args)
    out=OUTPUT/args.pilot_id
    out.mkdir(parents=True,exist_ok=False)
    mark={"pid":os.getpid(),"source_sha256":sha(SOURCE),"protocol_sha256":args.protocol_sha256,
          "preflight_sha256":args.preflight_sha256,"raw_component_receipt_sha256":args.raw_receipt_sha256,
          "raw_component_terminal_sha256":args.raw_terminal_sha256,"raw_component_os_exit_sha256":args.raw_os_exit_sha256}
    save(out/"started.json",mark|{"status":"RUNNING"})
    end=mark|{"status":"FAILED","exit_code":1}
    try:
        rows,n,remaining=compute(args,raw,out)
        good=not remaining
        checks=raw["checks"]|{
            "all_point_cache_force_reconstructed_and_contained":True,
            "all_point_force_and_whole_H_contained":good,
            "exact_tail_inflation_erosion_replayed_and_old_result_preserved":True,
        }
        h=save(out/"receipt.json",raw|mark|{
            "status":"PASS_INDEPENDENT_ACTUAL_SINGLETON_BRIDGE_COMPONENT_REPLAY" if good else "HOLD_INDEPENDENT_EXACT_ENDPOINT_CONTAINMENT",
            "checks":checks,"scope":SCOPE,"components":rows,
            "original_raw_component_receipt_file":name(RAW/args.pilot_id/"receipt.json"),
            "exact_endpoint_refinement_counts":n,"remaining_noncontained":remaining,
            "new_force_Hessian_moment_and_AD_jobs":0})
        code=0 if good else 2
        end|={"status":"FINISHED" if good else "FINISHED_WITH_HOLD","exit_code":code,"receipt_sha256":h}
        print(json.dumps({"fixed":n,"remaining":len(remaining)}),flush=True)
        return code
    except BaseException as error:
        save(out/"failure.json",{"error":str(error),"traceback":traceback.format_exc()})
        raise
    finally:
        save(out/"terminal.json",end)
        print(json.dumps(end),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    sub=parser.add_subparsers(dest="action",required=True)
    sub.add_parser("prepare")
    r=sub.add_parser("run")
    r.add_argument("--pilot-id",required=True)
    r.add_argument("--raw-receipt-sha256",required=True)
    r.add_argument("--raw-terminal-sha256",required=True)
    r.add_argument("--raw-os-exit",required=True,type=Path)
    r.add_argument("--raw-os-exit-sha256",required=True)
    r.add_argument("--protocol-sha256",required=True)
    r.add_argument("--preflight-sha256",required=True)
    args=parser.parse_args()
    prepare() if args.action=="prepare" else sys.exit(run(args))
