"""Independent local centered packet replay on the explicitly admitted C or D domain."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import importlib
from itertools import combinations_with_replacement, product
import json
from math import factorial
from pathlib import Path, PurePosixPath
from unittest.mock import patch

from flint import arb, ctx as precision, fmpq

ROOT = Path(__file__).resolve().parents[1]
BASE = "runs/exact_integration_cohort_v1"
REDUCTION = BASE + "/remaining15_centered_reduction_preparation"
CALLER = BASE + "/remaining15_centered_actual_preparation"
PREP = ROOT / REDUCTION / "independent_packet_preparation"
OUTPUT = ROOT / REDUCTION / "independent_actual"
REDUCER = "scripts/exact_centered_cohort_reduction.py"
SOURCE = "scripts/review_centered_cohort_packets.py"
RECORDS = tuple("optoGEF2_small_col_" + x for x in "CD")
REDUCER_ID = {
    "source_sha256": "af64f05e165a119d14be761ea79f258d10ac371f6f97de653b85c1c5f79f18a2",
    "protocol_sha256": "ec5419359e9ecfe21ba73a01ae1a1b72fed8561942920eaea43e27921236999d",
    "preflight_sha256": "984954626ade96dbb6ad802f7e465d674aed9087555ee3b419970adc85c22b82",
}
REDUCER_GATE = "d7911fc87fa703fea389823f4386bfbb3728cd6d4a8af6fd6463bf96e3b3b64a"
CALLER_GATE = "31ed5e5a15e73242414a44ed3b9412679d01d0e590eedbf4d94973dcbb556ca2"
PINS = {
    REDUCER: REDUCER_ID["source_sha256"],
    REDUCTION + "/protocol.json": REDUCER_ID["protocol_sha256"],
    REDUCTION + "/preflight.json": REDUCER_ID["preflight_sha256"],
    REDUCTION + "/prelaunch_independent_review.json": REDUCER_GATE,
    CALLER + "/prelaunch_independent_review.json": CALLER_GATE,
    CALLER + "/manifest.json": "8cf34d5a9fded6cc6c5a761d57e157da292c315d871bf34c2ec352d707b9ce65",
    BASE + "/remaining15_centered_preparation/optoGEF2_small_col_C_proofs.json":
        "9c3ff54d940735bae7db651056dfbf44c4dd893596520c4dfda2e276cb773d87",
    BASE + "/remaining15_centered_preparation/optoGEF2_small_col_D_proofs.json":
        "13d51f9626f125eeb34f429a55cea05298ff2494643da679fc048d66ee772269",
}
ZERO = (Q(0), Q(0))
EDGES = ((0, 1), (1, 2), (2, 0))


def require(ok, why):
    if not ok:
        raise ValueError(why)


def sha(raw):
    return sha256(raw).hexdigest()


def stable(x):
    return sha(json.dumps(x, default=str, separators=(",", ":"), sort_keys=True).encode())


def relative(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def checked(path, expected, *, packet=False, parse=True):
    local = (ROOT / path).resolve()
    require(local.is_relative_to(ROOT) and relative(local) == path, "canonical explicit artifact required")
    require(packet or f"{CALLER}/actual/" not in path or "/components/" not in path,
            "no centered packet read before complete independent admission")
    raw = local.read_bytes()
    require(sha(raw) == expected, f"changed artifact: {path}")
    return json.loads(raw) if parse else raw


def write(path, item):
    raw = (json.dumps(item, default=str, indent=2, allow_nan=False) + "\n").encode()
    with path.open("xb") as out:
        out.write(raw)
    return sha(raw)


def alltrue(checks):
    return isinstance(checks, dict) and bool(checks) and all(x is True for x in checks.values())


def bounds(packet):
    v = tuple(Q(packet[k]) for k in ("lower_exact", "upper_exact"))
    require(v[0] <= v[1] and all(x.denominator.bit_count() == 1 for x in v),
            "finite ordered dyadic stored endpoints required")
    return v


def enclose(outer, inner):
    return outer[0] <= inner[0] <= inner[1] <= outer[1]


def overlap(left, right):
    return max(left[0], right[0]) <= min(left[1], right[1])


def meet(left, right):
    require(overlap(left, right), "disjoint separately valid intervals")
    return max(left[0], right[0]), min(left[1], right[1])


def add(*rows):
    return tuple(sum((v[i] for v in rows), Q(0)) for i in (0, 1))


def mul(x, y):
    values = [a*b for a, b in product(x, y)]
    return min(values), max(values)


def magnitude(v):
    return max(map(abs, v))


def rational(x):
    return Q.from_float(float(x))


def ball(x):
    return arb(fmpq(x.numerator, x.denominator))


def packed_ball(x):
    require(x.is_finite(), "finite evaluated independent ball required")
    return tuple(Q(str(t.fmpq())) for t in (x.lower(), x.upper()))


def physical_derivatives(corners):
    r1, z1 = [corners[1][i] - corners[0][i] for i in (0, 1)]
    r2, z2 = [corners[2][i] - corners[0][i] for i in (0, 1)]
    det = r1*z2-r2*z1
    require(det > 0, "positive original reference triangle area")
    return det, ((z1-z2)/det, z2/det, -z1/det), ((r2-r1)/det, -r2/det, r1/det)


def p2(lam, deriv):
    return ([(4*x-1)*a for x,a in zip(lam,deriv)]
            + [4*(lam[i]*deriv[j] + lam[j]*deriv[i]) for i,j in EDGES])


def dot(l, r):
    return sum((a*b for a,b in zip(l,r,strict=True)), Q(0))


def exact_strain(rad, ax, dr, dz, shape, radius):
    return [dot(dr,rad), dot(dz,rad), dot(dr,ax), dot(dz,ax), dot(shape,rad)/radius]


def line_J(qs, qc):
    line = [tuple(sorted((s,c))) for s,c in zip(qs,qc,strict=True)]
    ad, bc = mul(line[0],line[3]), mul(line[1],line[2])
    det = ad[0]-bc[1]
    require(det > 0 and line[4][0] > 0 and det*line[4][0] > Q(1,10**6),
            "exact same-Y point strain line must lie in the positive-J analytic domain")
    return det*line[4][0]


def bulk_gradient(q, drive):
    """Independent analytic point strain gradient, evaluated on exact real point inputs."""
    a,b,c,d,h = map(ball,q)
    det = a*d - b*c
    J = det*h
    require(det > 0 and h > 0 and J > ball(Q(1,10**6)), "analytic cache derivative domain")
    penalty = 5000*J.log() - 1
    gain = ball(drive)
    return [a+penalty*d/det, (1+gain)*b-penalty*c/det,
            c-penalty*b/det, (1+gain)*d+penalty*a/det, h+penalty/h]


def stiffness_moment(corners, det, dz):
    """Exact integral of R times the product of physical axial P2 derivatives."""
    vertices = [[Q(int(a==b)) for a in range(3)] for b in range(3)]
    D = list(zip(*(p2(v,dz) for v in vertices)))
    radii = [c[0] for c in corners]
    def simplex(a,b,c):
        return Q(factorial(3) if a==b==c else 2 if len({a,b,c})==2 else 1, factorial(5))
    W = [[det*sum((r*simplex(i,j,k) for i,r in enumerate(radii)),Q(0))
          for k in range(3)] for j in range(3)]
    return {(i,j):sum((W[k][l]*D[i][k]*D[j][l] for k,l in product(range(3),repeat=2)),Q(0))
            for i,j in combinations_with_replacement(range(6),2)}


def independent_point_and_gain(domain, i, ids, selected):
    """Rebuild all 64 exact site maps, point cache correction and affine gain coefficients."""
    target = domain.context.target
    a = target.a
    corners = [tuple(map(rational, a["reference"][j])) for j in ids[:3]]
    det, dR, dZ = physical_derivatives(corners)
    center = target.nodes()
    rad, ax = zip(*(center[j] for j in ids))
    qrule = domain.rule
    exact_weights = [rational(x)/target.energy_scale for x in a["volume_weights"][i]]
    require(exact_weights == qrule["volume_weights"][i] and min(exact_weights)>0 and len(exact_weights)==64,
            "exact positive prior normalized rule weights")
    active = rational(a["active"][i,0])
    require(active in (0,1), "element constant active mask")
    drive = active*domain.context.drive
    min_J = None
    cache_force = [arb(0) for _ in selected]
    cached_stiffness = dict.fromkeys(combinations_with_replacement(range(6),2),Q(0))
    for k,(lam,w) in enumerate(zip(qrule["volume_barycentric"],exact_weights,strict=True)):
        Ns = [x*(2*x-1) for x in lam] + [4*lam[u]*lam[v] for u,v in EDGES]
        bsR,bsZ = p2(lam,dR),p2(lam,dZ)
        R = dot(lam,[v[0] for v in corners])
        Rc = rational(a["qr"][i,k])
        Nc = list(map(rational,a["shape"][k]))
        bcR,bcZ = zip(*(tuple(map(rational,z)) for z in a["shape_grad"][i,k]))
        xi = Nc[1]+(Nc[3]+Nc[4])/2
        eta = Nc[2]+(Nc[4]+Nc[5])/2
        require(R>0 and Rc>0 and tuple(lam)==(1-xi-eta,xi,eta) and min(lam)>0,
                "exact interior Qstar and binary cached shape sites")
        qs = exact_strain(rad,ax,bsR,bsZ,Ns,R)
        qc = exact_strain(rad,ax,bcR,bcZ,Nc,Rc)
        lo = line_J(qs,qc)
        min_J = lo if min_J is None else min(min_J,lo)
        gs,gc = bulk_gradient(qs,drive),bulk_gradient(qc,drive)
        ww = ball(w)
        for p, column in enumerate(selected):
            j, component = divmod(column,2)
            s = ball(bsR[j])*gs[2*component]+ball(bsZ[j])*gs[2*component+1]
            c = ball(bcR[j])*gc[2*component]+ball(bcZ[j])*gc[2*component+1]
            if component==0:
                s += ball(Ns[j]/R)*gs[4]
                c += ball(Nc[j]/Rc)*gc[4]
            cache_force[p] += ww*(s-c)
        if active:
            for u,v in cached_stiffness:
                cached_stiffness[u,v] += w*bcZ[u]*bcZ[v]
    if not active:
        return min_J, max(exact_weights), list(map(packed_ball,cache_force)), [ZERO]*len(selected)
    M = stiffness_moment(corners,det,dZ)
    N = 2*arb.pi()/ball(target.energy_scale)
    K = {key:N*ball(m)-ball(cached_stiffness[key]) for key,m in M.items()}
    y = [domain.nodes[j][a] for j in ids for a in (0,1)]
    out = []
    for i in selected:
        terms = (K[tuple(sorted((i//2,j//2)))]*y[j] for j in range(i%2,12,2))
        out.append(ball(domain.context.dose)*sum(terms,arb(0)))
    return min_J,max(exact_weights),list(map(packed_ball,cache_force)),list(map(packed_ball,out))


def read_admission(M,D,rt,record,admission_sha,external):
    W,H,runtime,A,plan,jobs,expect,hashes,ready = rt
    path = f"{REDUCTION}/{record}_input_admission.json"
    a = checked(path,admission_sha)
    require(a["status"]=="BOUND_144_CENTERED_PACKETS_BEFORE_LOCAL_READ" and a["record_id"]==record
            and a["admission_schema_version"]==2 and a["local_centered_payloads_opened"] is False
            and a["review_sha256"]==REDUCER_GATE and a["external_terminal_pins"]==external,
            "exact root-supplied admission/OS manifest pins and version required")
    require(M._own_gate(D,REDUCER_GATE,hashes)==REDUCER_ID,"same current centered reduction gate")
    reference,helper = M._declared_record(D,ready,record,hashes)
    files,links = M._terminal_inputs(D,record,reference,external,hashes)
    require(a["expected_packets"]==files and a[M.ADMITTED_WORKER_LINKS]==links,
            "same full externally observed worker manifest and links")
    terminal = checked(f"{CALLER}/actual/{record}/terminal.json",external["terminal_sha256"])
    sup = PurePosixPath(external["supervisor_dir"])
    require(sup.parent.as_posix()==CALLER+"/launches" and len(sup.name)==64
            and set(sup.name)<=set("0123456789abcdef"),"canonical supervisor directory")
    wait = checked(f"{sup}/{record}.os_terminal.json",external["os_terminal_sha256"])
    require(type(wait["exit_code"]) is type(terminal["exit_code"]) is int
            and wait["exit_code"]==terminal["exit_code"]==0
            and type(terminal["pid"]) is int and terminal["pid"]==wait["pid"]>0
            and terminal["results_sha256"]==external["manifest_sha256"],"true witnessed successful OS and worker exits")
    prefix = f"{CALLER}/actual/{record}/components"
    require(len(files)==144 and all(type(x["element"]) is int and x["element"]==i
            and x["file"]==f"{prefix}/element_{i:03d}.json" for i,x in enumerate(files)),
            "all and only canonical 144 centered packets")
    for f,h in a["input_hashes"].items():
        checked(f,h,parse=False)
    return a,reference,helper,files,links,terminal["pid"]


def packet_endpoints(p, proof, old, independently_cached, independently_gain):
    """Exact endpoint reconstruction using reviewed volume/shear kernel premises."""
    point = p["point_force"]
    n = len(proof["global_indices"])
    def field(v,k):
        x = list(map(bounds,v[k]))
        require(len(x)==n, f"complete {k} signed coordinates")
        return x
    vol, she, ca, F = (field(point,k) for k in ("volumetric_truncation","quadratic_strain_defect",
                       "finite_same_Y_cache_correction","normalized_force_total"))
    G,dg,h,wide = (field(p,k) for k in ("affine_gain_derivative_at_Y","affine_gain_increment",
                                       "uniform_H_times_nodal_displacement","centered_force_enclosure"))
    prior = field(old,"force_integral_minus_cache")
    inds = set(combinations_with_replacement(range(n),2))
    H = {}
    for row in old["H_integral_minus_cache_upper"]:
        key=(row["i"],row["j"])
        require(key in inds and key not in H and all(type(x) is int for x in key),"same signed local H coordinates")
        H[key]=bounds(row["enclosure"])
    require(H.keys()==inds,"all retained local H entries")
    tails = list(map(Q,point["two_mass_tail_upper_exact"]))
    r = list(map(Q,proof["original_local_radii_exact"]))
    gdelta = tuple(map(Q,proof["common_gain_shift_exact"]))
    out = [tuple(map(Q,x)) for x in p["force_intersection_with_prior_exact"]]
    require(len(r)==len(tails)==len(out)==len(independently_cached)==len(independently_gain)==n
            and len(gdelta)==2 and gdelta[0]<=0<=gdelta[1] and min(tails)>=0
            and max(tails)==Q(proof["point_force_tail_max_exact"])<=Q(proof["force_tail_budget_exact"]),
            "same exact point tails and verified original nodal/gain segment")
    counts=Counter()
    for i in range(n):
        require(enclose(ca[i],independently_cached[i]), "producer cache enclosure does not contain independent direct same-Y gradient")
        # Both gain evaluations use different outward grouping of the same exact affine expression.
        require(overlap(G[i],independently_gain[i]),"analytic independent physical/cached gain expression disagrees")
        if proof["active"]=="0":
            require(G[i]==dg[i]==independently_gain[i]==ZERO,"inactive gain must be exact zero")
        point_components=add(vol[i],she[i],ca[i],(-tails[i],tails[i]))
        require(overlap(F[i],point_components),"signed point force components inconsistent")
        u=sum((magnitude(H[min(i,j),max(i,j)])*v for j,v in enumerate(r)),Q(0))
        require(enclose(h[i],(-u,u)) and enclose(dg[i],mul(G[i],gdelta)),"same-box H and affine gain increments")
        centered_components=add(F[i],h[i],dg[i])
        require(overlap(wide[i],centered_components),"signed centered state+gain composition inconsistent")
        require(out[i]==meet(prior[i],wide[i]),"exact original-vs-centered scalar interval intersection")
        counts["coordinates"]+=1
        counts["point_contains_serialized_sum"]+=int(enclose(F[i],point_components))
        counts["centered_contains_serialized_sum"]+=int(enclose(wide[i],centered_components))
        counts["gain_contains_independent_ball"]+=int(enclose(G[i],independently_gain[i]))
    return counts


def boot():
    for p,h in PINS.items():
        checked(p,h,parse=False)
    M = importlib.import_module("exact_centered_cohort_reduction")
    require(relative(M.__file__)==REDUCER,"same frozen centered reducer path")
    rt=M._boot()
    D=rt[2][1].D
    require(M.ADMISSION_VERSION==2 and M.ADMITTED_WORKER_LINKS=="expected_worker_links"
            and tuple(rt[0].RECORDS)==RECORDS,"only explicit authenticated C/D v2 jobs")
    M._own_gate(D,REDUCER_GATE,{})
    return M,rt,D


def fixture_tests():
    precision.prec=640
    corners=((Q(1),Q(0)),(Q(3),Q(0)),(Q(1),Q(2)))
    det,dR,dZ=physical_derivatives(corners)
    M=stiffness_moment(corners,det,dZ)
    require(sum(M[i,j]*(1 if i==j else 2) for i,j in M)==0
            and M[2,2]==Q(7,10),"exact weighted P2 axial stiffness and constant kernel")
    q=[Q(1),Q(0),Q(0),Q(1),Q(1)]
    require(list(map(packed_ball,bulk_gradient(q,Q(1,4))))==[ZERO,ZERO,ZERO,(Q(1,4),Q(1,4)),ZERO]
            and line_J(q,q)==1,"independent analytic derivative at identity and positive cache line")
    bad=q.copy();bad[4]=Q(0)
    try:
        line_J(q,bad)
    except ValueError:
        pass
    else:
        raise AssertionError("nonpositive synthetic line accepted")
    require(mul((Q(-2),Q(3)),(Q(-1),Q(4)))==(Q(-8),Q(12))
            and meet((Q(-3),Q(1,2)),(Q(-1),Q(5)))==(Q(-1),Q(1,2)),
            "signed exact product and intersection")
    def box(a, b):
        return {"lower_exact": str(a), "upper_exact": str(b)}

    local = {"global_indices":[7], "original_local_radii_exact":["1/8"], "common_gain_shift_exact":["-1/16","1/8"],
             "point_force_tail_max_exact":"1/16", "force_tail_budget_exact":"1/4", "active":"1"}
    old = {"force_integral_minus_cache":[box(Q(-1,2),Q(1,2))],
           "H_integral_minus_cache_upper":[{"i":0,"j":0,"enclosure":box(2,3)}]}
    point = {"volumetric_truncation":[box(Q(1,8),Q(1,4))], "quadratic_strain_defect":[box(Q(-1,4),Q(1,4))],
             "finite_same_Y_cache_correction":[box(0,Q(1,8))], "normalized_force_total":[box(Q(-3,16),Q(11,16))],
             "two_mass_tail_upper_exact":["1/16"]}
    packet = {"point_force":point, "affine_gain_derivative_at_Y":[box(Q(-1,2),Q(1,4))],
              "affine_gain_increment":[box(Q(-1,16),Q(1,32))], "uniform_H_times_nodal_displacement":[box(Q(-3,8),Q(3,8))],
              "centered_force_enclosure":[box(Q(-5,8),Q(35,32))], "force_intersection_with_prior_exact":[["-1/2","1/2"]]}
    cache, gain = [(Q(1,16),Q(1,16))], [(Q(-1,4),Q(1,8))]
    require(packet_endpoints(packet,local,old,cache,gain)
            == {"coordinates":1,"point_contains_serialized_sum":1,"centered_contains_serialized_sum":1,
                "gain_contains_independent_ball":1},"signed decomposition and independent analytic point/gain containment")
    try:
        packet_endpoints({**packet, "uniform_H_times_nodal_displacement":[box(Q(-1,4),Q(1,4))]},local,old,cache,gain)
    except ValueError:
        pass
    else:
        raise AssertionError("narrowed wrong Hessian segment accepted")
    try:
        checked(CALLER+"/actual/"+RECORDS[0]+"/components/element_000.json","0"*64)
    except ValueError as error:
        require("no centered packet read" in str(error),"unauthorized read reaches artifact access")
    else:
        raise AssertionError("pre-admission component read accepted")
    return {"synthetic_independent_exact_axial_stiffness":True,"synthetic_analytic_point_cache_gradient":True,
            "synthetic_nonanalytic_strain_rejected":True,"signed_exact_interval_and_intersection":True,
            "synthetic_mixed_force_composition_and_narrowed_segment_rejection":True,
            "pre_admission_centered_packet_read_rejected":True}


def prepare():
    M,rt,D=boot()
    W,H,runtime,A,plan,jobs,expect,hashes,ready=rt
    blocked = AssertionError("actual force, moments, cache sums or centered I/O during reviewer preparation")
    with (patch.object(M,"_terminal_inputs",side_effect=blocked),
          patch.object(W,"_factory",side_effect=blocked),
          patch.object(H,"bulk_component",side_effect=blocked),
          patch.object(runtime[1].F,"_force_parts",side_effect=blocked),
          patch.object(runtime[1].E.core,"local_terms",side_effect=blocked)):
        for r in RECORDS:
            ref,proof=M._declared_record(D,ready,r,hashes)
            require(len(proof["bulk_proofs"])==144==len(ref["elements"]) and len(ref["retained_surface_manifest"])==24,
                    "reviewer targets same complete authenticated 144 proofs per record")
        checks=fixture_tests()
    checks.update({"only_C_D_declared_proofs_read":True,"no_actual_packets_sites_cache_sums_or_force_integrals_evaluated":True})
    PREP.mkdir(parents=True,exist_ok=True)
    src=Path(__file__).read_bytes()
    with (PREP/"review_centered_cohort_packets_frozen.py").open("xb") as f:
        f.write(src)
    ps=write(PREP/"protocol.json",{"status":"FROZEN_INDEPENDENT_CENTERED_PACKET_REVIEW_PREPARATION",
        "source_sha256":sha(src),"frozen_input_hashes":PINS,"records":list(RECORDS),
        "scope":"Local endpoint and cache/gain replay only. Previously reviewed high-order volume/shear moment evaluations and the prior whole-box H/geometry kernels are shared premises. Reconstruct signed sums, exact scalar tails, same-Y analytic point finite cache correction, 64 positive exact strain-line domains per element, and affine physical dZ stiffness at explicit 640-bit nodes.",
        "review_uses_no_new_spatial_force_moments_or_equilibrium":True,
        "packet_read_gate":"Explicit root-provided SHA of v2 immutable admission and matching complete worker results, terminal, true wait receipt, and canonical supervisor. Twelve axis and 132 off-axis cases per record; each worker must have true terminal before any new packet read.",
        "endpoint_policy":"Exact containment for independently reassembled state and gain increments and direct cache bounds; exact intersections with prior. Point/component and affine gain comparisons also record containment counts. Separately serialized Arb balls can widen different terms, so consistency requires overlap in those comparisons and shared validated primitive enclosures."})
    pf=write(PREP/"preflight.json",{"status":"PASS_INDEPENDENT_CENTERED_PACKET_PREPARATION_ONLY",
        "source_sha256":sha(src),"protocol_sha256":ps,"checks":checks,
        "actual_component_packets_read":False,"actual_cache_sites_replayed":False,"new_spatial_integrals":False})
    print(json.dumps({"source_sha256":sha(src),"protocol_sha256":ps,"preflight_sha256":pf,"checks":len(checks)}))


def actual(args):
    src=Path(__file__).read_bytes()
    pre=checked(relative(PREP/"preflight.json"),args.review_preflight_sha256)
    p=checked(relative(PREP/"protocol.json"),args.review_protocol_sha256)
    require(pre["status"]=="PASS_INDEPENDENT_CENTERED_PACKET_PREPARATION_ONLY" and alltrue(pre["checks"])
            and p["source_sha256"]==pre["source_sha256"]==sha(src)
            and pre["protocol_sha256"]==args.review_protocol_sha256 and p["records"]==list(RECORDS)
            and p["frozen_input_hashes"]==PINS
            and (PREP/"review_centered_cohort_packets_frozen.py").read_bytes()==src,"same frozen independent review program")
    record=args.record
    require(record in RECORDS and not (OUTPUT/f"{record}_packet_review.json").exists(),"unique explicit current C/D review")
    supplied={"worker_review_sha256":CALLER_GATE,"terminal_sha256":args.centered_terminal_sha256,
              "manifest_sha256":args.centered_manifest_sha256,"os_terminal_sha256":args.os_terminal_sha256,
              "supervisor_dir":args.supervisor_dir}
    M,rt,D=boot()
    W,H,runtime,A,plan,jobs,expect,hashes,ready=rt
    admission,reference,helper,manifest,links,pid=read_admission(M,D,rt,record,args.admission_sha256,supplied)
    domain,domains,routes,prior,canonical,key=W._factory(H,runtime,plan,jobs,expect,record,helper,hashes)
    require(key==reference["factory_nodeball_sha256"]==links["factory_nodeball_sha256"] and precision.prec==640
            and prior["surface"]==reference["retained_surface_manifest"],"same unchanged exact full factory and all retained surfaces")
    counters=Counter()
    charts=Counter()
    details=[]
    for i,entry in enumerate(manifest):
        receipt,old,selector,route,proof=W._element_guard(H,runtime,domain,canonical,domains,routes,prior,helper,i,hashes)
        require(receipt==reference["elements"][i],"same frozen center/route metadata and prior H")
        packet=checked(entry["file"],entry["sha256"],packet=True)
        require(all(packet.get(k)==v for k,v in {**proof,**links}.items())
                and packet["status"]=="CENTERED_RECORD_BULK_FORCE_UNREVIEWED" and packet["precision_bits"]==640
                and packet["centered_proof_sha256"]==stable(proof)==receipt["centered_proof_sha256"]
                and packet["verified_local_nodeball_sha256"]==receipt["nodeball_sha256"]
                and packet["new_Hessian_integrations"] is False,"strict packet identity and no new Hessian")
        E=runtime[1].E
        center=domain.context.target.nodes()
        fields,ids=E._bulk_fields(domain.context.target,center,i)
        jet=E.JET.bulk_jacobian_jet(fields["q"],fields["B"])
        t,_=E.M.volumetric_derivative_tail_bounds(selector["degree"],Q(selector["rho_exact"]),1,5000)
        mass=Q(selector["geometry_two_mass_upper_exact"])
        tails=[mass*t*E._norm(x,fields) for x in jet["gradient"]]
        pt=packet["point_force"]
        require(fields["global_indices"]==proof["global_indices"] and fields["free_local"]==proof["free_local"]
                and fields["chart"]["kind"]==proof["chart"] and ids==proof["local_nodes"]
                and t==Q(selector["phi_prime_tail_upper_exact"])
                and tails==list(map(Q,pt["two_mass_tail_upper_exact"])),"exact per-coordinate first-derivative tails and original projection")
        min_J,w,cache,gain=independent_point_and_gain(domain,i,ids,fields["free_local"])
        J=bounds(pt["cache_min_J_bound"])
        require(Q(1,10**6)<J[0]<=J[1]<=min_J and Q(pt["largest_cached_weight_exact"])==w,
                "all jointly reconstructed exact positive same-Y cache lines and same normalized rule")
        scalar=packet_endpoints(packet,proof,old,cache,gain)
        counters.update(scalar)
        charts[proof["chart"]]+=1
        details.append({"element":i,"actual_centered_packet_sha256":entry["sha256"],
            "same_64_site_point_cache_and_J_verified":True,
            "independent_min_cache_line_J_exact":str(min_J),"signed_endpoint_counts":dict(scalar)})
        if (i+1)%24==0:
            print(json.dumps({"record_id":record,"independent_packets":i+1,"of":144}),flush=True)
    require(charts=={"off_axis":132,"axis_vertex":6,"axis_edge_physical":6} and len(details)==144,
            "full explicit centered coverage")
    checks=dict.fromkeys(("externally_pinned_full_worker_manifest_and_true_OS_exit",
        "admission_before_actual_payload_reads","same_exact_base_C_full_Y_nodeballs_and_prior_H",
        "all_144_packet_metadata_and_signed_scalar_endpoint_replays",
        "all_9216_exact_analytic_same_Y_cache_domains_and_direct_point_correction",
        "exact_degree3_radial_weight_axial_gain_and_frozen_binary_cache_weights",
        "held_original_24_surface_manifests_and_full_box_derivative",
        "no_new_high_order_spatial_integrals_or_equilibrium_solves"),True)
    report={"status":"PASS_CENTERED_PACKET_REPLAY",**REDUCER_ID,"record_id":record,"components":144,
        "admission_sha256":args.admission_sha256,"centered_manifest_sha256":args.centered_manifest_sha256,
        "centered_terminal_sha256":args.centered_terminal_sha256,"os_terminal_sha256":args.os_terminal_sha256,
        "worker_pid":pid,"checks":checks,"review_source_sha256":sha(src),
        "review_protocol_sha256":args.review_protocol_sha256,"review_preflight_sha256":args.review_preflight_sha256,
        "signed_endpoint_counts":dict(counters),"route_counts":dict(charts),"all_packet_details":details,
        "scope":"Local packet verification, including independently recomputed analytic point cache correction and exact low-degree axial gain. Rounded serialized scalar accumulations checked by exact containment where required and overlap where independently rounded. Shared reviewed high-order spatial volume/shear primitives plus prior whole-box H and geometric certificate. No assembled root or physics claim."}
    OUTPUT.mkdir(parents=True,exist_ok=True)
    pin=write(OUTPUT/f"{record}_packet_review.json",report)
    M._packet_review(D,record,pin,args.admission_sha256,supplied,REDUCER_ID,{})
    print(json.dumps({"status":report["status"],"record_id":record,"packet_review_sha256":pin,
                      "signed_endpoint_counts":dict(counters)},allow_nan=False))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_subparsers(dest="mode",required=True)
    modes.add_parser("prepare")
    run=modes.add_parser("actual")
    run.add_argument("record",choices=RECORDS)
    for k in ("review-protocol-sha256","review-preflight-sha256","admission-sha256",
              "centered-manifest-sha256","centered-terminal-sha256","os-terminal-sha256","supervisor-dir"):
        run.add_argument("--"+k,required=True)
    ns=parser.parse_args()
    prepare() if ns.mode=="prepare" else actual(ns)

