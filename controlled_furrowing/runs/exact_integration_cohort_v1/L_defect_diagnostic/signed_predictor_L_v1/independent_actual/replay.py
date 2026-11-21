"""Independent 256-bit replay of L's signed full integral predictor residual."""
import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
import os
from pathlib import Path
import sys
from time import monotonic
from unittest.mock import patch

from flint import arb, arb_mat, ctx

OUT = Path(__file__).resolve().parent
STAGE = OUT.parent
ROOT = STAGE.parents[3]
PREP = OUT / "preparation"
GATE = STAGE / "actual_independent_review.json"
CENTERED = STAGE.parent / "centered_bulk_L_v1"
PREC = 256
PID = 25152
N = 611
ORDER = ("apical_area_modulus","basal_area_modulus","apical_surface_tension",
         "basal_surface_tension","common_gain")
PINS = {
    "scripts/algebraic_root_feasibility.py":"7321d38bf17cf354cc03be2871f3fbaa45bcebbcbfca7be90d0f867b80078acf",
    "scripts/review_cohort_exact_record.py":"67db748a5d40f31471ddd3d6c84fae280ff7ead650e6eb4e950a4766e54e0285",
    "scripts/exact_cohort_context.py":"c22612a1de1ccfd34e8b0f0859f6413c4da09c765ba4945d060340e93e928f9f",
    "scripts/exact_signed_predictor_L.py":"00cebb37ed7cd23152ee46a659ff9878a32f62733faa2964d3b051fa7678dd1b",
    "preparation/protocol.json":"4e7152f01104da4c5f09a8134f54ac3228ae6467b7ae34d862f50265b5a003af",
    "preparation/preflight.json":"e91f11f8f7f82aae6670fe4215750df2506e59f7501205c8faeec76577ce99e3",
    "preparation/prelaunch_independent_review.json":"d67a83d4d4e5b11c8ca6c6d81a5ad01c07bfb3a10b704d4cf7a2b85c7313500f",
    "actual/started.json":"9f0f1d5ffde4c2896a9e3b151cde4ee2b41362ddfb5d1bbf716c06ba76880284",
    "actual/terminal.json":"f3d3e4b2bc6c9d71e0815683dae52f7c46ca02cc1f51ffbc3b186c98d721d42c",
    "actual/results.json":"4af745c36cc3150b94451b1df39c1442c5d3efaf17a2781cfbfc4bbce301a382",
}
CENTER_PINS = {
    "reduction_preparation/results.json":"5f07f0aea44bc176d615de9e593e16695cc478082b7cf9fc9847aada2a0604e8",
    "reduction_preparation/actual_independent_review.json":"cce6f01147b6ed436e3f3a8620e94bbd9df84a8fc1aa45e9bbed9cba8c1f8c4a",
    "independent_actual/component_and_composition_replay.json":"cf918c5c12a06437aee8c0de9ca76de0da48db67ceaf22ea1ea36c542a9665c8",
}
USE = ("separated_point_bulk_RF_intervals_exact","separated_retained_surface_RF_intervals_exact",
       "coupled_RH_bulk_S_intervals_exact","separated_gain_derivative_RF_intervals_exact",
       "exact_closed_parameter_offsets")


def read_set():
    hashes = {}

    def read(name,expected=None,parse=True):
        f = (ROOT/name).resolve()
        assert f.is_relative_to(ROOT)
        raw = f.read_bytes()
        h = sha256(raw).hexdigest()
        assert expected is None or h == expected, name
        assert hashes.setdefault(f.relative_to(ROOT).as_posix(),h) == h
        return json.loads(raw) if parse else raw
    return hashes,read


def fingerprint(value):
    return sha256(json.dumps(value,default=str,sort_keys=True,separators=(",",":")).encode()).hexdigest()


def pair(x):
    lo,hi = map(Q,x)
    assert lo <= hi
    return lo,hi


def ball_pair(x):
    assert x.is_finite()
    return Q(str(x.lower().fmpq())),Q(str(x.upper().fmpq()))


def aloft(x):
    assert x.is_finite()
    return x.abs_upper()


def add(*xs):
    return tuple(sum((x[i] for x in xs),Q(0)) for i in (0,1))


def mul(a,b):
    candidates = [x*y for x,y in product(a,b)]
    return min(candidates),max(candidates)


def inside(x,ball):
    a,b = ball_pair(ball)
    assert a <= x <= b


def overlap(a,b):
    return max(a[0],b[0]) <= min(a[1],b[1])


def write(path,d):
    with path.open("x",encoding="utf8",newline="\n") as f:
        f.write(json.dumps(d,default=str,indent=2)+"\n")
    return sha256(path.read_bytes()).hexdigest()


def stamp():
    return datetime.now(timezone.utc).isoformat()


def setup(hashes,read):
    for path,h in PINS.items():
        if path.startswith("scripts/"):
            read(path,h,False)
    landing = {path:read(STAGE/path,h) for path,h in PINS.items() if not path.startswith("scripts/")}
    frozen,pre,gate = [landing[f"preparation/{k}.json"] for k in
                       ("protocol","preflight","prelaunch_independent_review")]
    assert pre["status"] == "PASS_SIGNED_PREDICTOR_PREPARATION_ONLY"
    assert gate["status"] == "PASS" and all(v is True for v in gate["checks"].values())
    assert len(pre["checks"]) == 17 and all(v is True for v in pre["checks"].values())
    assert pre["actual_cached_9504_site_sums"] is frozen["actual_cached_point_or_variation_run"] is False
    for path,h in {**frozen["input_hashes"],**frozen["failed_attempts"]}.items():
        read(path,h,False)
    start,term,prod = [landing[f"actual/{k}.json"] for k in ("started","terminal","results")]
    assert start["pid"] == term["pid"] == PID and term["status"] == "FINISHED"
    assert type(term["exit_code"]) is int and term["exit_code"] == 0
    assert term["results_sha256"] == PINS["actual/results.json"]
    assert prod["status"] == "CONDITIONAL_UNREVIEWED_SIGNED_PREDICTOR_SUFFICIENT_TEST_PASSES"
    assert prod["independently_reviewed_root_assertion"] is False
    for name,key in (("scripts/exact_signed_predictor_L.py","source_sha256"),
                     ("preparation/protocol.json","protocol_sha256"),("preparation/preflight.json","preflight_sha256"),
                     ("preparation/prelaunch_independent_review.json","prelaunch_sha256")):
        assert prod[key] == start[key] == PINS[name]
    assert 0 < prod["seconds"] <= (datetime.fromisoformat(term["utc"])-datetime.fromisoformat(start["utc"])).total_seconds()+1
    assert all(prod["input_hashes"][path] == h for path,h in frozen["input_hashes"].items())
    for path,h in prod["input_hashes"].items():
        read(path,h,False)

    sys.path.insert(0,str(ROOT/"scripts"))
    import review_cohort_exact_record as D
    import exact_cohort_context as C
    import algebraic_root_feasibility as core
    assert Path(core.__file__).resolve() == (ROOT/"scripts/algebraic_root_feasibility.py").resolve()
    record = C.load_context("optoGEF2_small_col_L")
    assert record.record_id == prod["record_id"] == frozen["record_id"]
    assert record.row["parameter_cell_sha256"] == prod["same_closed_C_sha256"]
    assert frozen["parameter_order"] == list(ORDER) == list(C.ORDER)
    assert {k:Q(v) for k,v in prod["base_parameters_exact"].items()} == record.base
    assert record.drive == record.base["common_gain"] * record.dose > 0
    for path,h in record.input_hashes.items():
        read(path,h,False)
    source = f"{D._record_dir(record.record_id)}/protocol.json"
    original = read(source,frozen["input_hashes"][source])
    assert original["row"] == record.row
    prefs = frozen["point_kernel_arguments"]
    assert prefs["point_certificate_file"] == record.row["point_cache_file"]
    assert Path(prefs["predictor_file"]).resolve() == (ROOT/record.row["predictor_file"]).resolve()
    assert {k:Q(v) for k,v in prefs["base_parameters_exact"].items()} == record.base
    assert Q(prefs["dose_exact"]) == record.dose
    oldeta,oldZ,r,ref = D._old_premises(record,original,hashes)
    Jcenter = read(record.row["root_certificate_file"],record.row["root_certificate_sha256"])
    root_pre = Jcenter["uniform_cell_result"]["radius_attempts"]
    assert any(Q(x["radius_exact"]) == r and x["strict_selfmap_contraction"] is True for x in root_pre)
    assert r == record.old_radius == Q(prod["r_exact"])
    assert Q(prod["cached_Z_upper_exact"]) == oldZ
    centered = {p:read(CENTERED/p,h) for p,h in CENTER_PINS.items()}
    center = centered["reduction_preparation/results.json"]
    centerreview = centered["reduction_preparation/actual_independent_review.json"]
    central_replay = centered["independent_actual/component_and_composition_replay.json"]
    assert centerreview["status"] == "PASS_ACTUAL_CENTERED_ORIGINAL_CUBE_REDUCTION_REVIEW"
    assert centerreview["reducer_results_sha256"] == CENTER_PINS["reduction_preparation/results.json"]
    assert centerreview["component_replay_sha256"] == CENTER_PINS["independent_actual/component_and_composition_replay.json"]
    assert all(center[k] == value for k,value in centerreview["decision"].items())
    for key in USE:
        assert fingerprint(center[key]) == central_replay["expected_reducer_array_digests"][key]
    assert center["component_count"] == 168 and center["retained_surface_components"] == 24
    b = Q(center["b_upper_exact"])
    z = oldZ+b
    assert z == Q(center["z_upper_exact"]) == Q(prod["z_upper_exact"]) < 1
    assert Q(prod["strict_eta_budget_exact"]) == r*(1-z)
    for face,index in product(range(2),range(12)):
        path = f"{D._record_dir(record.record_id)}/surface/face{face}_edge{index:02d}.json"
        read(path,center["source_inputs_and_packets"][path],False)
    delta = [(record.cell[k][0]-record.base[k],record.cell[k][1]-record.base[k]) for k in ORDER]
    assert list(map(pair,center[USE[4]])) == list(map(pair,prod["same_closed_parameter_offsets"])) == delta
    assert all(a <= 0 <= x for a,x in delta)
    norm = max(sum((abs(Q.from_float(float(v))) for v in row),Q(0)) for row in record.R)
    assert norm == Q(center["R_infinity_norm_exact"]) > 0
    curv = D._minimum(record,hashes,oldZ,z,r,norm,ref)
    assert curv["strict_local_minimum"] is True
    assert Q(curv["normalized_curvature_lower_exact"]) == (1-z)/norm > 0
    return D,core,record,prod,center,delta,b,z,r,oldZ,norm,curv


def site_gradient_columns(term,dose,active,basal,A):
    q = term["q"]
    F = arb_mat(len(q),5)
    if term["kind"] == "bulk":
        F[1,4],F[3,4] = A(dose*active)*q[1],A(dose*active)*q[3]
    else:
        rr,x,y = q
        scale = term["parameters"][2]
        length = (x*x+y*y).sqrt()
        Js = scale*rr*length
        columns = (1,3) if basal else (0,2)
        dJs = (scale*length,scale*rr*x/length,scale*rr*y/length)
        for i,v in enumerate(dJs):
            F[i,columns[0]],F[i,columns[1]] = (Js-1)*v,v
    return F


def build_independently(core,record,delta):
    """Form H·S locally; never call the producer point/variation or build dense H."""
    A = core.A
    model,u,R = core.restore_cache(record.target.packet)
    assert u.tobytes() == record.u.tobytes() and R.tobytes() == record.R.tobytes()
    assert model.free.tobytes() == record.free.tobytes()
    assert Q(model.parameters.mu) == 1 and Q(model.parameters.bulk) == 5000
    assert Q(model.energy_scale) == record.target.energy_scale > 0
    for name in ORDER[:4]:
        assert Q(getattr(model.parameters,name)) == record.base[name]
    box_theta = [A(record.cell[k][0]).union(A(record.cell[k][1])) for k in ORDER]
    # All five exact offsets contain zero. This is the same interval hull of
    # t*S*delta on the theta-first, state-second path t in [0,1].
    interval_delta = arb_mat(5,1,[A(a).union(A(b)) for a,b in delta])
    g = [arb(0) for _ in range(N)]
    columns = [[arb(0) for _ in range(5)] for _ in range(N)]
    rem = [arb(0) for _ in range(N)]
    counts = Counter()
    counts_active = Counter()
    basal_count = len(model.surfaces[0][0])*len(model.surfaces[0][1])
    nb = len(model.shape)
    minJ0 = minJpath = minR = minTangent = None
    for term in core.local_terms(model,u,record.drive):
        kind = term["kind"]
        ids = term["indices"]
        q0 = term["q"]
        n = len(q0)
        p0 = term["parameters"]
        w = term["weight"]
        assert ids and len(set(ids)) == len(ids) and w > 0
        B = arb_mat(n,len(ids),[v for row in term["B"] for v in row])
        BT = B.transpose()
        S = arb_mat(len(ids),5,[A(float(x)) for pid in ids for x in record.S[pid]])
        BS = B*S
        delta_q = BS*interval_delta
        qpath = [q+delta_q[i,0] for i,q in enumerate(q0)]
        if kind == "bulk":
            element,site = divmod(counts[kind],nb)
            mask = int(record.target.a["active"][element,0])
            assert mask in (0,1)
            counts_active[mask] += 1
            f0,H0,J0 = core.bulk_local(q0,*p0)
            _,Ht,Jt = core.bulk_local(qpath,p0[0],p0[1],A(mask*record.dose)*box_theta[4])
            minJ0 = J0.lower() if minJ0 is None else min(minJ0,J0.lower())
            minJpath = Jt.lower() if minJpath is None else min(minJpath,Jt.lower())
            C = site_gradient_columns(term,record.dose,mask,False,A)
        else:
            is_basal = counts[kind] < basal_count
            i,j = (1,3) if is_basal else (0,2)
            f0,H0,_ = core.surface_local(q0,*p0)
            _,Ht,_ = core.surface_local(qpath,box_theta[i],box_theta[j],p0[2])
            assert p0[2] > 0
            rho = qpath[0].lower()
            tangent = (qpath[1]*qpath[1]+qpath[2]*qpath[2]).lower()
            assert rho > 0 and tangent > 0
            minR = rho if minR is None else min(minR,rho)
            minTangent = tangent if minTangent is None else min(minTangent,tangent)
            C = site_gradient_columns(term,record.dose,0,is_basal,A)
        counts[kind] += 1
        f_sub = BT*arb_mat(n,1,f0)
        cols = BT*(arb_mat(H0)*BS+C)
        uncertainty = [aloft(delta_q[i,0]) for i in range(n)]
        rest = [sum((aloft(Ht[i][j]-H0[i][j])*v for j,v in enumerate(uncertainty)),arb(0))
                for i in range(n)]
        for j,full in enumerate(ids):
            g[full] += w*f_sub[j,0]
            for k in range(5):
                columns[full][k] += w*cols[j,k]
            rem[full] += w*sum((aloft(B[i,j])*rest[i] for i in range(n)),arb(0))
        done = counts["bulk"]+counts["surface"]
        if done in (3072,6144,9216,9504):
            print(json.dumps({"cached_sites_checked":done,"total":9504}),flush=True)
    assert counts == {"bulk":9216,"surface":288} and counts_active[1] == 24*64
    R_arb = arb_mat(N,N,[A(float(x)) for x in R.ravel()])
    R_mag = arb_mat(N,N,[A(abs(float(x))) for x in R.ravel()])
    gr = R_arb*arb_mat(N,1,g)
    M = R_arb*arb_mat(columns)
    remainder = R_mag*arb_mat(N,1,[v.upper() for v in rem])
    assert minJ0 > A(Q(1,10**6)) and minJpath > A(Q(1,10**6))
    rg = [ball_pair(gr[i,0]) for i in range(N)]
    affine = [[ball_pair(M[i,k]) for k in range(5)] for i in range(N)]
    radius = [ball_pair(remainder[i,0])[1] for i in range(N)]
    assert min(radius) >= 0
    return rg,affine,radius,{"site_counts":dict(counts),"active_site_counts":dict(counts_active),
                              "cached_point_min_J":ball_pair(minJ0)[0],"cached_path_min_J":ball_pair(minJpath)[0],
                              "cached_path_min_surface_radius":ball_pair(minR)[0],
                              "cached_path_min_surface_tangent_squared":ball_pair(minTangent)[0]}


def predictor_enclosures(base,columns,rem,center,delta):
    anchors,coeffs,rows = [],[],[]
    point,surf,hs,gain = (center[key] for key in USE[:4])
    for i,(p,M,r) in enumerate(zip(base,columns,rem)):
        an = add(p,pair(point[i]),pair(surf[i]))
        coeff = [add(a,pair(h),pair(gain[i])) if k==4 else add(a,pair(h)) for k,(a,h) in enumerate(zip(M,hs[i]))]
        bounds = add(an,*[mul(a,s) for a,s in zip(coeff,delta)],(-r,r))
        anchors.append(an)
        coeffs.append(coeff)
        rows.append(bounds)
    return anchors,coeffs,rows


def synthetic(core):
    A = core.A
    q = list(map(A,(1,Q(1,8),0,1,1)))
    dose = Q(5,4)
    mask = 1
    term = {"kind":"bulk","q":q,"parameters":(A(1),A(5000),A(0))}
    C = site_gradient_columns(term,dose,mask,False,A)
    inside(Q(5,32),C[1,4])
    inside(dose,C[3,4])
    assert C[0,4].is_zero() and C[4,4].is_zero()
    assert all(v.is_zero() for v in site_gradient_columns(term,dose,0,False,A).entries())
    q = list(map(A,(2,3,4)))
    term = {"kind":"surface","q":q,"parameters":(A(1),A(0),A(Q(1,7)))}
    for basal in (False,True):
        C = site_gradient_columns(term,dose,0,basal,A)
        k,t = (1,3) if basal else (0,2)
        Js = Q(10,7)
        for i,v in enumerate((Q(5,7),Q(6,35),Q(8,35))):
            inside((Js-1)*v,C[i,k])
            inside(v,C[i,t])
    # Exact differentiable theta-first example: gradient of
    # E(u,g)=u^3/3+(1+g)u^2/2, u0=1, g0=2, S=-1/4.
    # The independent final-theta average H−H0 times Sδ equals the omitted term.
    u0,g0,S = Q(1),Q(2),Q(-1,4)
    Qforce = lambda u,g:u*u+(1+g)*u
    for d in (Q(-1,6),Q(0),Q(1,4)):
        move = S*d
        H0 = 2*u0+1+g0
        H_delta_avg = move+d
        cached = Qforce(u0,g0)+(H0*S+u0)*d+H_delta_avg*move
        assert cached == Qforce(u0+move,g0+d)
    zero = (Q(0),Q(0))
    center = {USE[0]:[(-1,-1)],USE[1]:[zero],USE[2]:[[zero,zero,zero,zero,(-2,-2)]],USE[3]:[(1,1)]}
    an,C,ret = predictor_enclosures([(Q(1),Q(1))],[[zero]*4+[(Q(1),Q(1))]],[Q(1,10)],
                                     center,[(Q(-1,8),Q(3,7))]*5)
    assert an == [zero] and C == [[zero]*5] and ret == [(Q(-1,10),Q(1,10))]
    assert not (Q(1,4)+Q(3,4) < 1)


def prepare():
    hashes,read = read_set()
    # No original cached objective, assembly or variation on this path.
    data = setup(hashes,read)
    _,core,record,_,_,_,b,z,r,_,norm,curv = data
    ctx.prec = PREC
    core.A.cache_clear()
    forbidden = AssertionError("preparation attempted an actual cached site")
    with patch.object(core,"restore_cache",side_effect=forbidden), \
            patch.object(core,"local_terms",side_effect=forbidden), \
            patch.object(core,"assemble_point",side_effect=forbidden):
        synthetic(core)
    f = Path(__file__).relative_to(ROOT)
    raw = read(f,parse=False)
    h = sha256(raw).hexdigest()
    PREP.mkdir(parents=True,exist_ok=True)
    snapshot = PREP / "replay_frozen.py"
    with snapshot.open("xb") as out:
        out.write(raw)
    protocol = {"status":"FROZEN_INDEPENDENT_SIGNED_Q_REPLAY_PREPARATION","record_id":record.record_id,
                "source_file":f.as_posix(),"source_sha256":h,"input_hashes":dict(sorted(hashes.items())),
                "precision_bits":PREC,"radius_exact":r,"b_exact":b,"z_exact":z,"R_infinity_norm_exact":norm,
                "same_record_SPD_curvature_lower_exact":curv["normalized_curvature_lower_exact"],
                "math":"Rebuild all cached g_Q0 and five H_Q0*S+B_Q0 columns by local chart derivatives, applying H to S before global assembly. Rebuild a rowwise theta-first and final-theta state-path remainder using H_path-H_0 at radius zero. Combine independently with existing reviewed signed defect arrays.",
                "shared_premises":"Frozen numeric cache restoration, site maps, local energy derivatives, original whole-box cached contraction, original/scaled integral H bound and reviewed centered force/gain arrays.",
                "actual_cached_sites_run":False}
    psha = write(PREP/"protocol.json",protocol)
    pre = {"status":"PASS_PREPARATION_ONLY","source_sha256":h,"protocol_sha256":psha,
           "checks":{"frozen_signed_job_terminal_and_exact_inputs_current":True,
                     "synthetic_signed_bulk_surface_parameter_columns":True,
                     "exact_theta_first_final_theta_path_and_five_shared_offsets":True,
                     "independent_local_H_S_signed_row_plan":True,
                     "cached_site_assembly_blocked_during_preparation":True,
                     "old_cached_SPD_and_same_R_H_radius_cell_on_original_box":True},
           "actual_cached_sites_run":False}
    pre_sha = write(PREP/"preflight.json",pre)
    print(json.dumps({"status":pre["status"],"source_sha256":h,"protocol_sha256":psha,
                      "preflight_sha256":pre_sha,"actual_cached_sites_run":False}),flush=True)


def replay(args):
    prep_hash,read_prep = read_set()
    protocol = read_prep(PREP/"protocol.json",args.protocol_sha)
    pre = read_prep(PREP/"preflight.json",args.preflight_sha)
    p = Path(__file__).relative_to(ROOT)
    s = read_prep(p,protocol["source_sha256"],False)
    assert s == read_prep(PREP/"replay_frozen.py",protocol["source_sha256"],False)
    assert pre["status"] == "PASS_PREPARATION_ONLY" and all(v is True for v in pre["checks"].values())
    assert pre["source_sha256"] == protocol["source_sha256"] and pre["protocol_sha256"] == args.protocol_sha
    hashes,read = read_set()
    D,core,record,prod,center,delta,b,z,r,oldZ,norm,curv = setup(hashes,read)
    read(p,protocol["source_sha256"],False)
    assert hashes == protocol["input_hashes"]
    for name,h in prep_hash.items():
        read(name,h,False)
    ctx.prec = PREC
    core.A.cache_clear()
    import torch
    torch.set_num_threads(1)
    t0 = monotonic()
    write(OUT/"started.json",{"utc":stamp(),"pid":os.getpid(),"protocol_sha256":args.protocol_sha,
                              "preflight_sha256":args.preflight_sha})
    try:
        rg,M,rem,physical = build_independently(core,record,delta)
        anchors,coefficients,rows = predictor_enclosures(rg,M,rem,center,delta)
        assert len(rg) == len(M) == len(rows) == N
        assert all(len(prod[key]) == N for key in ("signed_cached_anchor_gR","signed_cached_affine_M",
                                                  "signed_full_integral_RF_predictor_intervals"))
        assert all(len(row) == 5 for row in prod["signed_cached_affine_M"])
        assert all(overlap(a,pair(b)) for a,b in zip(rg,prod["signed_cached_anchor_gR"]))
        assert all(overlap(a,pair(b)) for x,y in zip(M,prod["signed_cached_affine_M"]) for a,b in zip(x,y))
        assert all(overlap(a,pair(b)) for a,b in zip(rows,prod["signed_full_integral_RF_predictor_intervals"]))
        eta = max(max(map(abs,x)) for x in rows)
        lhs = eta+z*r
        passed = bool(z < 1 and lhs < r)
        err = eta/(1-z) if passed else None
        assert not passed or err < r
        for name,h in tuple(hashes.items()):
            read(name,h,False)
        result = {
            "status":"INDEPENDENT_SIGNED_PREDICTOR_STRICT_TEST_PASSES" if passed else "INDEPENDENT_SIGNED_PREDICTOR_STRICT_TEST_UNRESOLVED",
            "scope":"Frozen exact spatial finite-element model, L only, selected original cube for every theta in the same closed C. Independent cached reconstruction, shared reviewed local derivatives and previously reviewed integral-defect and box-contraction premises.",
            "record_id":record.record_id,"same_closed_C_sha256":record.row["parameter_cell_sha256"],
            "precision_bits":PREC,"independent_cached_gR":rg,"independent_cached_affine_columns":M,
            "independent_cached_remainders":rem,"shared_parameter_offsets":delta,
            "independent_full_signed_anchors":anchors,"independent_full_signed_coefficients":coefficients,
            "independent_full_RF_predictor_intervals":rows,"strict_signed_force_budget_exact":r*(1-z),
            "independent_eta_I_exact":eta,"lhs_exact":lhs,"r_exact":r,"z_exact":z,"old_Z_Q_exact":oldZ,
            "b_exact":b,"R_infinity_norm_exact":norm,"strict_conditions_pass":passed,"root_error_upper_exact":err,
            "strict_local_minimum_curvature":curv if passed else None,
            "all_producer_rows_overlap_independent_enclosures":True,"no_producer_build_point_variation_or_combination":True,
            "domains":physical,"input_hashes":dict(sorted(hashes.items())),"seconds":monotonic()-t0,
        }
        digest = write(OUT/"results.json",result)
        write(OUT/"terminal.json",{"status":"FINISHED","utc":stamp(),"pid":os.getpid(),"exit_code":0,"results_sha256":digest})
        print(json.dumps({"status":result["status"],"result_sha256":digest,
                          "eta_I":float(eta),"lhs":float(lhs),"error":float(err) if passed else None,"r":float(r),
                          "z":float(z),"curvature":float(Q(curv["normalized_curvature_lower_exact"])) if passed else None}),flush=True)
    except BaseException as error:
        write(OUT/"failure.json",{"error":str(error),"type":type(error).__name__})
        write(OUT/"terminal.json",{"status":"FAILED","exit_code":1,"pid":os.getpid(),"utc":stamp()})
        raise


def seal(args):
    hashes,read = read_set()
    terminal = read(OUT/"terminal.json",args.terminal_sha)
    result = read(OUT/"results.json",args.results_sha)
    started = read(OUT/"started.json")
    assert terminal["status"] == "FINISHED" and terminal["exit_code"] == 0
    assert terminal["pid"] == started["pid"] and terminal["results_sha256"] == args.results_sha
    for name,h in result["input_hashes"].items():
        read(name,h,False)
    eta,lhs,r,z = (Q(result[k]) for k in ("independent_eta_I_exact","lhs_exact","r_exact","z_exact"))
    assert lhs == eta+z*r
    passed = z < 1 and lhs < r
    assert passed is result["strict_conditions_pass"]
    if passed:
        assert Q(result["root_error_upper_exact"]) == eta/(1-z) < r
        assert Q(result["strict_local_minimum_curvature"]["normalized_curvature_lower_exact"]) == (1-z)/Q(result["R_infinity_norm_exact"]) > 0
    review = {"status":"PASS_INDEPENDENT_SIGNED_ORIGINAL_CUBE_ROOT_CERTIFICATE" if passed else "PASS_INDEPENDENT_SIGNED_REPLAY_UNRESOLVED",
              "record_id":result["record_id"],"same_closed_C_sha256":result["same_closed_C_sha256"],
              "scope":result["scope"],"independent_actual_results_sha256":args.results_sha,
              "independent_actual_terminal_sha256":args.terminal_sha,
              "producer_results_sha256":PINS["actual/results.json"],"producer_terminal_sha256":PINS["actual/terminal.json"],
              "checks":{"matching_external_process_exit_and_pinned_results":True,"all_inputs_frozen":True,
                        "authoritative_9504_cache_sites_exact_signed_rows_and_independent_remainder_replayed":True,
                        "all_cached_anchor_affine_and_integral_rows_overlap_producer":True,
                        "independent_strict_conditions_and_same_record_curvature_recomputed":True},
              "sufficient_test_pass":passed,"root_certified":passed,"strict_local_minimum_certified":passed,
              "uniform_eta_I_upper_exact":str(eta),"uniform_selfmap_lhs_exact":str(lhs),
              "original_radius_exact":str(r),"same_box_z_exact":str(z),
              "root_error_upper_exact":result["root_error_upper_exact"],
              "uniform_curvature_lower_exact":result["strict_local_minimum_curvature"]["normalized_curvature_lower_exact"] if passed else None,
              "site_domain_and_normalization":result["domains"],"input_hashes":dict(sorted(hashes.items()))}
    h = write(GATE,review)
    print(json.dumps({"status":review["status"],"review_sha256":h,"root_certified":passed}),flush=True)


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Replay requires enabled assertions")
    root = argparse.ArgumentParser()
    sub = root.add_subparsers(dest="mode",required=True)
    sub.add_parser("prepare")
    rep = sub.add_parser("replay")
    rep.add_argument("--protocol-sha",required=True)
    rep.add_argument("--preflight-sha",required=True)
    post = sub.add_parser("seal")
    post.add_argument("--terminal-sha",required=True)
    post.add_argument("--results-sha",required=True)
    args = root.parse_args()
    if args.mode == "prepare":
        prepare()
    elif args.mode == "replay":
        replay(args)
    else:
        seal(args)
