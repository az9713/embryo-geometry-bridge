"""Independent signed global centered replay and exact rational gain-domain certification."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from math import ceil, floor, lcm
from pathlib import Path
from unittest.mock import patch

from flint import arb, ctx as precision
import review_centered_cohort_packets as P

ROOT, BASE, REDUCTION, CALLER = P.ROOT, P.BASE, P.REDUCTION, P.CALLER
SOURCE = Path(__file__).resolve()
PREP = ROOT / REDUCTION / "independent_global_preparation"
OUT = ROOT / REDUCTION / "independent_global_actual"
PACKET = BASE + "/remaining15_centered_reduction_preparation/independent_packet_preparation"
PINS = {
    **P.PINS,
    "scripts/review_centered_cohort_packets.py": "4f7e6450f57f0c368f3334646df32a9411e7e333402d5c8ec4b05d17397117b3",
    PACKET + "/protocol.json": "0086494f3a69e2587d04ba8824e4a9e3768c46f9671a7f582c72f62fbeefc25b",
    PACKET + "/preflight.json": "3111f11d91f9fe954a4b8ef5b9231cbe79b4abc33d8febf66d95670fbba7113b",
}
ZERO = P.ZERO
RANK = 611
KEYS = ("apical_area_modulus", "basal_area_modulus", "apical_surface_tension", "basal_surface_tension", "common_gain")


def check(ok, msg):
    if not ok:
        raise ValueError(msg)


def power(q):
    need = q.denominator
    check(need.bit_count() == 1, "exact dyadic precision")
    return need.bit_length()-1


def scaled(q, e):
    return q.numerator << (e-power(q))


def binary_matrix(m):
    q = [[P.rational(x) for x in row] for row in m]
    e = max(power(x) for row in q for x in row)
    return [[scaled(x,e) for x in row] for row in q],e


def midrad(intervals):
    e = max((power(x) for a in intervals for x in a),default=0)
    acc = []
    for a,b in intervals:
        check(a<=b,"ordered endpoint matrix")
        x,y=scaled(a,e),scaled(b,e)
        acc.append((x+y,y-x))
    return acc,e+1


def r_force(R,e,values):
    m,q = midrad(values)
    div = 2**(e+q)
    out=[]
    for row in R:
        center = sum(w*x for w,(x,y) in zip(row,m,strict=True))
        extent = sum(abs(w)*y for w,(x,y) in zip(row,m,strict=True))
        out.append((Q(center-extent,div),Q(center+extent,div)))
    return out


def r_absolute(R,e,values):
    ex = max(map(power,values))
    v = [scaled(q,ex) for q in values]
    return [Q(sum(abs(c)*u for c,u in zip(row,v,strict=True)),2**(e+ex)) for row in R]


def r_h_s(R,eR,H,radii,S,eS):
    # From exact sums of midpoints and independent radii for signed RH and RH*S.
    n=len(R)
    check(len(radii)==len(S)==n and all(len(row)==5 for row in S),"full 5-coordinate slope")
    keys=list(H)
    mr,eh=midrad([H[k] for k in keys])
    by_index=[[] for _ in range(n)]
    for (i,j),(m,r) in zip(keys,mr,strict=True):
        by_index[i].append((j,m,r))
        if i!=j:
            by_index[j].append((i,m,r))
    d=lcm(*(q.denominator for q in radii))
    hw=[q.numerator*(d//q.denominator) for q in radii]
    S_cols=list(zip(*S))
    box,coef=[],[]
    for row in R:
        m,u=[0]*n,[0]*n
        for weight,terms in zip(row,by_index,strict=True):
            for j,center,radius in terms:
                m[j]+=weight*center
                u[j]+=abs(weight)*radius
        div=2**(eh+eR)
        box.append(Q(sum((abs(a)+b)*r for a,b,r in zip(m,u,hw,strict=True)),d*div))
        div*=2**eS
        this=[]
        for s in S_cols:
            a=sum(w*x for w,x in zip(s,m,strict=True))
            b=sum(abs(w)*x for w,x in zip(s,u,strict=True))
            this.append((Q(a-b,div),Q(a+b,div)))
        coef.append(this)
    return box,coef


def h_rows(H,n):
    out=[Q(0)]*n
    for (i,j),v in H.items():
        for k in {i,j}:
            out[k]+=P.magnitude(v)
    return out


def force_add(output,ids,values):
    for i,v in zip(ids,values,strict=True):
        output[i]=P.add(output[i],v)


def scatter_H(output,ids,entries):
    expect=set(combinations_with_replacement(range(len(ids)),2))
    seen=set()
    for x in entries:
        ij=(x["i"],x["j"])
        check(ij in expect and ij not in seen,"local upper triangle maps every index once")
        seen.add(ij)
        k=tuple(sorted(ids[c] for c in ij))
        output[k]=P.add(output.get(k,ZERO),P.bounds(x["enclosure"]))
    check(seen==expect,"complete symmetric local Hessian")


def compose(point,RHS,gain,surface,offsets):
    return [P.add(p,s,*(P.mul(P.add(v, g) if j==4 else v, offsets[j]) for j,v in enumerate(row)))
            for p,row,g,s in zip(point,RHS,gain,surface,strict=True)]


def pi_interval():
    with precision.workprec(1024):
        lo,hi=P.packed_ball(arb.pi())
    check(3<lo<=hi<4 and hi-lo<=Q(1,2**1000),"rigorous fine pi hull")
    return lo,hi


def outward_dyadic(a,exponent=800):
    two=2**exponent
    return Q(floor(a[0]*two),two),Q(ceil(a[1]*two),two)


def exact_gain(domain,i,ids,free,pi):
    """Independent exact range on original rational nodal box, with certified pi interval."""
    ctx=domain.context
    a=ctx.target.a
    active=P.rational(a["active"][i,0])
    check(active in (0,1),"original active element mask")
    if not active:
        return [ZERO]*len(free)
    corners=tuple(tuple(map(P.rational,a["reference"][j])) for j in ids[:3])
    d,_,dz=P.physical_derivatives(corners)
    material=P.stiffness_moment(corners,d,dz)
    Qc=dict.fromkeys(material,Q(0))
    weights=domain.rule["volume_weights"][i]
    check(len(weights)==64 and min(weights)>0,"same normalized rule")
    for k,w in enumerate(weights):
        gz=[P.rational(c[1]) for c in a["shape_grad"][i,k]]
        check(w==P.rational(a["volume_weights"][i,k])/ctx.target.energy_scale,"exact binary normalization")
        for j,l in Qc:
            Qc[j,l]+=w*gz[j]*gz[l]
    N=P.mul(pi,(2/ctx.target.energy_scale,)*2)
    coeff={k:P.add(P.mul((v,v),N),(-Qc[k],-Qc[k])) for k,v in material.items()}
    YY=[]
    for j in ids:
        for k in (0,1):
            x,r=Q(domain.box["current_centers_exact"][j][k]),Q(domain.box["coordinate_radii_exact"][j][k])
            YY.append((x-r,x+r))
    return [P.mul((ctx.dose,)*2,P.add(*(P.mul(coeff[tuple(sorted((i//2,j//2)))],YY[j])
                    for j in range(i%2,12,2)))) for i in free]


def minmag(a):
    return min(abs(a[0]),abs(a[1])) if a[0]*a[1]>0 else Q(0)


def strict(eta0,z0,d,b,r):
    eta=eta0+d
    z=z0+b
    lhs=eta+z*r
    ok=z<1 and lhs<r
    return eta,z,lhs,ok,eta/(1-z) if ok else None


def fixtures():
    pi=pi_interval()
    check(P.enclose((Q(314,100),Q(315,100)),pi),"precision-controlled certified pi interval")
    test=(Q(-1,3),Q(1,7))
    check(P.enclose(outward_dyadic(test),test) and not P.enclose((Q(-1,4),Q(1,4)),test),"round outward with negative endpoints")
    R,e=binary_matrix([[.5,-1.],[1.,1.]])
    S,se=binary_matrix([[1.,0.,0.,0.,.5],[.5,.25,0.,0.,-.5]])
    f=[(Q(-2),Q(1)),(Q(-1),Q(3,2))]
    HF={(0,0):(Q(-1),Q(2)),(0,1):(Q(-1,2),Q(1,2)),(1,1):(Q(1),Q(1))}
    r=[Q(1,3),Q(1,9)]
    check(r_force(R,e,f)==[(Q(-5,2),Q(3,2)),(Q(-3),Q(5,2))],"exact signed RF")
    def dot(row,entries):return P.add(*(P.mul((w,w),x) for w,x in zip(row,entries,strict=True)))
    rr=[[Q(v,2**e) for v in a] for a in R]
    rh=[[dot(c,[HF[tuple(sorted((i,j)))] for i in range(2)]) for j in range(2)] for c in rr]
    slopes=[[dot(col,u) for col in zip(*[[Q(v,2**se) for v in vals] for vals in S])] for u in rh]
    expected=[sum(P.magnitude(v)*t for v,t in zip(u,r)) for u in rh]
    check(r_h_s(R,e,HF,r,S,se)==(expected,slopes),"exact signed RH then common predictor columns")
    check(r_absolute(R,e,h_rows(HF,2))==[Q(11,4),Q(4)],"absolute H and preconditioner bounds")
    g=(Q(2),Q(2)); row=[ZERO]*4+[(Q(-2),Q(-2))]
    check(compose([ZERO],[row],[g],[(Q(1),Q(1))],[(Q(-1),Q(1))]*5)==[(Q(1),Q(1))],"same shared gain offset cancels")
    check(strict(Q(1,8),Q(1,4),Q(1,8),Q(1,4),Q(1,2))[-2:]==(False,None),"equality does not pass")
    return {"exact_signed_R_and_RH_S_compared_against_Fraction_oracle":True,"five_shared_offsets_gain_cancellation":True,
            "pi_hull_and_outward_slider_endpoints":True,"strict_equality_unresolved":True}


def prepare():
    for k,v in PINS.items():
        P.checked(k,v,parse=False)
    M,rt,D=P.boot()
    W,H,runtime,A,plan,jobs,expect,hashes,ready=rt
    forbidden=AssertionError("actual packet, gain, root/cached-site, or solve during global preparation")
    with (patch.object(W,"_factory",side_effect=forbidden),patch.object(M,"_calculate",side_effect=forbidden),
          patch.object(M,"_terminal_inputs",side_effect=forbidden),patch.object(H,"bulk_component",side_effect=forbidden),
          patch.object(runtime[1].E.core,"local_terms",side_effect=forbidden)):
        for record in P.RECORDS:
            reference,helper=M._declared_record(D,ready,record,hashes)
            check(len(reference["elements"])==len(helper["bulk_proofs"])==144 and len(helper["retained_surface_manifest"])==24,
                  "same explicit C and D frozen route proofs")
        checks=fixtures()
    checks.update({"same_current_C_D_admission_and_lifetime_pins":True,
                   "no_actual_packets_or_global_sums_in_preparation":True})
    PREP.mkdir(parents=True,exist_ok=True)
    src=SOURCE.read_bytes()
    with (PREP/"review_centered_cohort_global_frozen.py").open("xb") as f:f.write(src)
    p=P.write(PREP/"protocol.json",{"status":"FROZEN_INDEPENDENT_CENTERED_GLOBAL_REPLAY_PREPARATION",
        "source_sha256":P.sha(src),"input_hashes":PINS,"records":list(P.RECORDS),
        "scope":"Rebuild all signed endpoint arrays and unchanged full H from admitted C/D packets. Independent centered preconditioning and exact Fraction gain-range proof on original rational nodal hull with a 1024-bit pi enclosure. Export alternative gain/RG and predictor bound if not contained. Shared reviewed high-order volume/shear, independently replayed point finite cache, geometry and prior H/old Q premises. No mechanics or integration jobs.",
        "review_requires_actual_packet_receipt_and_OS_terminal":True,
        "predictor_only_eta_whole_box_full_H_derivative":True,
        "no_actual_payloads_read":True})
    pf=P.write(PREP/"preflight.json",{"status":"PASS_INDEPENDENT_CENTERED_GLOBAL_PREPARATION_ONLY",
          "source_sha256":P.sha(src),"protocol_sha256":p,"checks":checks,"no_actual_payloads_read":True})
    print(json.dumps({"source_sha256":P.sha(src),"protocol_sha256":p,"preflight_sha256":pf}))


def parse_flag(args,n):
    return getattr(args,n.replace("-","_"))


def terminal_documents(args):
    record=args.record
    packet_path=f"{REDUCTION}/independent_actual/{record}_packet_review.json"
    packet_os=P.checked(f"{REDUCTION}/independent_actual/{record}_packet_review_os_observation.json",args.packet_os_sha256)
    packet=P.checked(packet_path,args.packet_review_sha256)
    result=P.checked(f"{REDUCTION}/{record}/results.json",args.results_sha256)
    term=P.checked(f"{REDUCTION}/{record}/terminal.json",args.terminal_sha256)
    check(packet_os["status"]=="OBSERVED_TRUE_ZERO_EXIT_CENTERED_PACKET_REVIEW"
          and packet_os["observed_exit_code"]==0 and packet_os["review_path"]==packet_path
          and packet_os["review_sha256"]==args.packet_review_sha256
          and packet["status"]=="PASS_CENTERED_PACKET_REPLAY" and packet["components"]==144
          and P.alltrue(packet["checks"]) and packet["admission_sha256"]==args.admission_sha256
          and result["admission_sha256"]==args.admission_sha256 and term["status"]=="FINISHED"
          and type(term["pid"]) is int and term["pid"]>0 and type(term["exit_code"]) is int and term["exit_code"]==0
          and term["results_sha256"]==args.results_sha256
          and result["actual_packet_review_sha256"]==args.packet_review_sha256
          and all(result[k]==packet[k]==v for k,v in P.REDUCER_ID.items())
          and packet["record_id"]==packet_os["record_id"]==result["record_id"]==record
          and result["independent_root_assertion"] is False,"successful matching external producer terminal and independently verified packet review")
    return packet,packet_os,result


def run(args):
    prep=P.checked(P.relative(PREP/"protocol.json"),args.review_protocol_sha256)
    pre=P.checked(P.relative(PREP/"preflight.json"),args.review_preflight_sha256)
    src=SOURCE.read_bytes()
    check(pre["status"]=="PASS_INDEPENDENT_CENTERED_GLOBAL_PREPARATION_ONLY" and P.alltrue(pre["checks"])
          and pre["source_sha256"]==prep["source_sha256"]==P.sha(src)
          and prep["input_hashes"]==PINS and pre["protocol_sha256"]==args.review_protocol_sha256
          and (PREP/"review_centered_cohort_global_frozen.py").read_bytes()==src,"frozen independent replay source")
    for k,v in PINS.items():P.checked(k,v,parse=False)
    record=args.record
    check(record in P.RECORDS and not (OUT/record).exists(),"unique explicit record review")
    packet_review,packet_os,producer=terminal_documents(args)
    # All of the following checks complete before a centered packet is opened.
    M,rt,D=P.boot()
    W,H,runtime,A,plan,jobs,expect,hashes,ready=rt
    admission=P.checked(f"{REDUCTION}/{record}_input_admission.json",args.admission_sha256)
    pins=admission["external_terminal_pins"]
    check(pins["manifest_sha256"]==packet_review["centered_manifest_sha256"]
          and pins["terminal_sha256"]==packet_review["centered_terminal_sha256"]
          and pins["os_terminal_sha256"]==packet_review["os_terminal_sha256"],
          "same admitted worker manifest and OS observation")
    admission,reference,proofs,packets,links,work_pid=P.read_admission(M,D,rt,record,args.admission_sha256,pins)
    domain,actual,routes,prior_files,canonical,nodehash=W._factory(H,runtime,plan,jobs,expect,record,proofs,hashes)
    ctx,p=domain.context,domain.protocol
    check(nodehash==links["factory_nodeball_sha256"] and len(prior_files["surface"])==24
          and len(packet_review["all_packet_details"])==144 and work_pid==packet_review["worker_pid"]
          and producer["same_closed_C_sha256"]==ctx.row["parameter_cell_sha256"]
          and Q(producer["same_exact_source_dose"])==ctx.dose,"same full C, record, dose, nodeballs and prior H")
    f,point,gain,gain_alt,surface=([ZERO]*RANK for _ in range(5))
    Hbulk,Htotal={},{}
    force_count,h_count=[0]*RANK,[0]*RANK
    pi=pi_interval()
    strict_contained,gain_count=0,0
    misses=[]
    for i,entry in enumerate(packets):
        previous=prior_files["bulk"][i]
        proof=proofs["bulk_proofs"][i]
        witness=packet_review["all_packet_details"][i]
        check(witness["element"]==entry["element"]==previous["element"]==i
              and witness["actual_centered_packet_sha256"]==entry["sha256"]
              and witness["same_64_site_point_cache_and_J_verified"] is True,"independent local packet replay binds actual hash")
        current=P.checked(entry["file"],entry["sha256"],packet=True)
        old=P.checked(previous["file"],previous["sha256"])
        ids=proof["global_indices"]
        n=len(ids)
        check(all(current.get(k)==v for k,v in {**proof,**links}.items())
              and current["verified_local_nodeball_sha256"]==old["nodal_hulls_sha256"]==reference["elements"][i]["nodeball_sha256"]
              and current["centered_proof_sha256"]==P.stable(proof)==reference["elements"][i]["centered_proof_sha256"]
              and old["global_indices"]==ids and old["record_id"]==record
              and previous["sha256"]==proof["uniform_packet_sha256"],"same old and new local original box/map")
        oldf=list(map(P.bounds,old["force_integral_minus_cache"]))
        intersections=[P.meet(a,b) for a,b in zip(oldf,map(P.bounds,current["centered_force_enclosure"]),strict=True)]
        check(intersections==[tuple(map(Q,v)) for v in current["force_intersection_with_prior_exact"]],
              "exact signed local intersection")
        raw=list(map(P.bounds,current["affine_gain_derivative_at_Y"]))
        independently=exact_gain(domain,i,proof["local_nodes"],proof["free_local"],pi)
        local_alt=[]
        for j,(x,y) in enumerate(zip(raw,independently,strict=True)):
            check(P.overlap(x,y),"independent rigorous exact-box gain interval does not intersect reviewed formula")
            contained=P.enclose(x,y)
            strict_contained+=int(contained)
            gain_count+=1
            if not contained:
                misses.append({"element":i,"local_coordinate":j,"producer":list(map(str,x)),
                              "independent_outward_dyadic":list(map(str,outward_dyadic(y)))})
            local_alt.append(outward_dyadic(y) if y!=ZERO else ZERO)
        force_add(f,ids,intersections)
        force_add(point,ids,list(map(P.bounds,current["point_force"]["normalized_force_total"])))
        force_add(gain,ids,raw)
        force_add(gain_alt,ids,local_alt)
        scatter_H(Hbulk,ids,old["H_integral_minus_cache_upper"])
        for j in ids:force_count[j]+=1;h_count[j]+=n
        if (i+1)%48==0:print(json.dumps({"record_id":record,"bulk_replayed":i+1,"gain_exact_contained":strict_contained,"gain_coordinates":gain_count}),flush=True)
    Htotal.update(Hbulk)
    for i,entry in enumerate(prior_files["surface"]):
        q=P.checked(entry["file"],entry["sha256"])
        face,edge=divmod(i,12)
        local=ctx.target.a[f"surface_{face}_edges"][edge]
        ids=[ctx.target.index[k] for j in local for k in (2*int(j),2*int(j)+1) if k in ctx.target.index]
        check(q["face"]==face and q["edge"]==edge and q["global_indices"]==ids and q["record_id"]==record
              and q["point_cache_sha256"]==ctx.row["point_cache_sha256"] and q["parameter_cell_sha256"]==ctx.row["parameter_cell_sha256"],
              "same verified retained C/box surface data")
        values=list(map(P.bounds,q["force_integral_minus_cache"]))
        force_add(f,ids,values)
        force_add(surface,ids,values)
        scatter_H(Htotal,ids,q["H_integral_minus_cache_upper"])
        for j in ids:force_count[j]+=1;h_count[j]+=len(ids)
    check(force_count==p["force_incidence"] and h_count==p["H_row_incidence"] and len(Htotal)==6676
          and len(f)==611 and len(prior_files["surface"])==24,"complete signed original FE topology")
    prior_audit=P.checked(f"{H.DIRECT}/independent_actual/{record}.json",reference["operator_review_sha256"])
    prior=P.checked(f"{H.DIRECT}/{record}/results.json",prior_audit["producer_pins"]["results_sha256"])
    entries=[[i,j,*Htotal[i,j]] for i,j in sorted(Htotal)]
    check(P.stable(entries)==P.stable(prior["global_H_upper_intervals_exact"]),"original full H entries unchanged")
    R,eR=binary_matrix(ctx.R)
    S,eS=binary_matrix(ctx.S)
    hr=h_rows(Htotal,RANK)
    b=max(r_absolute(R,eR,hr))
    norm=Q(max(sum(map(abs,row)) for row in R),2**eR)
    check(b==Q(prior["b_upper_exact"])==Q(reference["unchanged_fullbox_b_exact"])
          and norm==Q(prior["R_infinity_norm_exact"]),"unchanged whole-box Hessian defect and exact R norm")
    rf,rp,rg,rs,alt_rg=[r_force(R,eR,rows) for rows in (f,point,gain,surface,gain_alt)]
    radii=[Q(domain.box["coordinate_radii_exact"][int(j)//2][int(j)%2]) for j in ctx.free]
    st,RHS=r_h_s(R,eR,Hbulk,radii,S,eS)
    order=tuple(runtime[1].C.ORDER)
    offsets=[tuple(x-ctx.base[k] for x in ctx.cell[k]) for k in order]
    check(order==KEYS,"exact same 5 ordered closed offset intervals")
    sep=[P.add(a,(-b,b),P.mul(c,offsets[-1]),d) for a,b,c,d in zip(rp,st,rg,rs,strict=True)]
    joint=[P.meet(a,b) for a,b in zip(rf,sep,strict=True)]
    pred=compose(rp,RHS,rg,rs,offsets)
    final=[P.meet(a,b) for a,b in zip(joint,pred,strict=True)]
    alt_pred=compose(rp,RHS,alt_rg,rs,offsets)
    # Prior direct whole-box RF is independently reviewed and available even when a gain comparison is only overlapping.
    old_rf=[tuple(map(Q,a)) for a in prior["preconditioned_RF_intervals_exact"]]
    alt_final=[P.meet(a,b) for a,b in zip(alt_pred,old_rf,strict=True)]
    arrays={"global_F_intervals_exact":f,"retained_global_H_upper_intervals_exact":entries,
            "direct_whole_box_RF_intervals_exact":rf,"separated_whole_box_RF_intervals_exact":sep,
            "whole_box_RF_intersections_exact":joint,"coupled_predictor_RF_only_intervals_exact":pred,
            "final_predictor_RF_only_intersections_exact":final,"point_bulk_RF_intervals_exact":rp,
            "RH_bulk_S_signed_intervals_exact":RHS,"affine_gain_RF_intervals_exact":rg,
            "retained_surface_RF_intervals_exact":rs,"node_path_bulk_only_RH_radii_exact":st,
            "parameter_order":order,"closed_five_offsets_exact":offsets}
    digests={k:P.stable(v) for k,v in arrays.items()}
    check(all(P.stable(producer[k])==h for k,h in digests.items()),"exact agreement of all signed/whole H/predictor arrays")
    eta0,z0,r,provider=D._old_premises(ctx,p,hashes)
    d=max(map(P.magnitude,final))
    d_lo=max(map(minmag,final))
    alt_d=max(map(P.magnitude,alt_final))
    eta,z,lhs,ok,error=strict(eta0,z0,d,b,r)
    eta_alt,z_alt,lhs_alt,ok_alt,error_alt=strict(eta0,z0,alt_d,b,r)
    vals={"exact_R_infinity_norm":norm,"same_fullbox_b_exact":b,"old_cached_eta_exact":eta0,"old_cached_Z_exact":z0,
          "same_radius_exact":r,"RF_predictor_defect_upper_exact":d,"RF_predictor_defect_lower_exact":d_lo,
          "total_eta_upper_exact":eta,"total_Z_upper_exact":z,"strict_selfmap_lhs_exact":lhs}
    check(all(Q(producer[k])==v for k,v in vals.items()) and producer["sufficient_conditions_passed"] is ok
          and producer["conditional_root_error_exact"]==(None if error is None else str(error))
          and producer["status"]==("CONDITIONAL_UNREVIEWED_CENTERED_SUFFICIENT_TEST_PASSES" if ok
                                   else "UNREVIEWED_CENTERED_SUFFICIENT_TEST_UNRESOLVED"),"strict root arithmetic readback")
    all_contained=strict_contained==gain_count
    final_ok=ok and (all_contained or ok_alt)
    curvature=None
    if final_ok:
        proof=D._minimum(ctx,hashes,z0,z,r,norm,provider)
        check(proof["strict_local_minimum"] is True and Q(proof["normalized_curvature_lower_exact"])==(1-z)/norm,
              "same-record SPD homotopy")
        curvature=Q(proof["normalized_curvature_lower_exact"])
    future={"status":"INDEPENDENT_EXACT_CENTERED_INPUT_ARRAYS","record_id":record,
        "same_closed_C_sha256":ctx.row["parameter_cell_sha256"],"point_bulk_RF_intervals_exact":rp,
        "RH_bulk_S_signed_intervals_exact":RHS,"retained_surface_RF_intervals_exact":rs,
        "independent_affine_gain_RF_intervals_exact":alt_rg,"verified_producer_gain_RF_intervals_exact":rg if all_contained else None,
        "raw_gain_operand_proven_to_cover_exact_box":all_contained,"closed_five_offsets_exact":offsets,"parameter_order":order,
        "independent_coupled_predictor_RF_intervals_exact":alt_pred,"independent_final_predictor_RF_intervals_exact":alt_final,
        "same_full_box_b_exact":b,"same_total_Z_upper_exact":z,"same_radius_exact":r,
        "normalization_norm_exact":norm,"gain_rational_pi_enclosure":pi,
        "declared_outward_dyadic_bits":800,"noncontaining_raw_gain_diagnostics":misses,
        "reducer_exact_array_digest":digests}
    folder=OUT/record
    folder.mkdir(parents=True,exist_ok=False)
    a=P.write(folder/"signed_inputs.json",future)
    report={"status":"PASS_INDEPENDENT_CENTERED_ROOT_AND_MINIMUM" if final_ok else "PASS_INDEPENDENT_CENTERED_REDUCTION_UNRESOLVED",
        "record_id":record,"root_certified":final_ok,"strict_local_minimum_certified":final_ok,
        "actual_admission_sha256":args.admission_sha256,"actual_packet_review_sha256":args.packet_review_sha256,
        "actual_packet_review_os_sha256":args.packet_os_sha256,"producer_results_sha256":args.results_sha256,
        "producer_terminal_sha256":args.terminal_sha256,"review_source_sha256":P.sha(src),
        "review_protocol_sha256":args.review_protocol_sha256,"review_preflight_sha256":args.review_preflight_sha256,
        "exact_array_sha256":digests,"global_F_entries":611,"H_upper_entries":6676,"first_144_bulk_and_retained_24_surfaces":True,
        "gain_slots":gain_count,"producer_gain_contains_exact_rational_nodal_gain_slots":strict_contained,
        "producer_gain_overlapping_only_slots":gain_count-strict_contained,
        "gain_check_scope":"The independent dyadic range follows from rational basis and full original exact nodal box times certified pi; an overlap alone never proves producer containment.",
        "future_signed_inputs_sha256":a,"scope":"Original closed C and original tube. Full H derivative b unchanged; signed predictor-only force bound. Local high-order volume/shear and prior H/geometry premises shared with reviewed components. No new mechanical evaluation or spatial high-order integration.",
        "scalar_exact_readback":vals,"raw_candidate_sufficient_test_passed":ok,
        "independent_gain_predictor_defect":alt_d,"independent_gain_predictor_lhs":lhs_alt,
        "independent_gain_predictor_passes":ok_alt,
        "certified_root_error_exact":(error if all_contained else error_alt) if final_ok else None,
        "normalized_curvature_lower_exact":curvature,"no_new_actual_mechanical_jobs":True}
    h=P.write(folder/"receipt.json",report)
    print(json.dumps({"status":report["status"],"record_id":record,"review_sha256":h,"signed_inputs_sha256":a,
                      "raw_gain_contains":strict_contained,"gain_slots":gain_count,
                      "lhs":float(lhs),"alternate_lhs":float(lhs_alt),"r":float(r),"z":float(z),
                      "defect_lower":float(d_lo),"defect_upper":float(d)}))


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    sub=a.add_subparsers(dest="mode",required=True)
    sub.add_parser("prepare")
    job=sub.add_parser("actual")
    job.add_argument("record",choices=P.RECORDS)
    for flag in ("review-protocol-sha256","review-preflight-sha256","admission-sha256","packet-review-sha256",
                 "packet-os-sha256","results-sha256","terminal-sha256"):
        job.add_argument("--"+flag,required=True)
    args=a.parse_args()
    prepare() if args.mode=="prepare" else run(args)

