"""Audit the complete centered L actual job and independently compose exact endpoints."""
import argparse
from datetime import datetime
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement, product
import json
from math import lcm
from pathlib import Path
import sys

OUT = Path(__file__).resolve().parent
STAGE = OUT.parent
ROOT = STAGE.parents[3]
RED = STAGE / "reduction_preparation"
RECEIPT = OUT / "component_and_composition_replay.json"
ACTUAL_GATE = STAGE / "actual_independent_review.json"
RED_GATE = RED / "actual_independent_review.json"
TERMINAL = "dc27f3e080bb0ecc49aaf023ca301b3dc7a5ee25e852eef13b6002f1cda33eba"
MANIFEST = "c1628822b0b20e278e2438802c5c721a3f7db1e168dd1e7a60914fde9372cb5b"
PRELAUNCH = "2e70e234fd719985326a009ee41d857da96bad45e029a4db6bba91ff0fba1807"
PID = 29844
N = 611
ZERO = (Q(0),Q(0))


def reader():
    hashes = {}

    def read(name,sha=None,parse=True):
        f = (ROOT/name).resolve()
        assert f.is_relative_to(ROOT)
        path = f.relative_to(ROOT).as_posix()
        raw = f.read_bytes()
        digest = sha256(raw).hexdigest()
        assert sha is None or sha == digest, path
        assert hashes.setdefault(path,digest) == digest, path
        return json.loads(raw) if parse else raw
    return hashes,read


def digest(value):
    return sha256(json.dumps(value,default=str,separators=(",",":"),sort_keys=True).encode()).hexdigest()


def pack(p):
    a,b = (Q(p[k]) for k in ("lower_exact","upper_exact"))
    assert a <= b
    return a,b


def add(*args):
    return tuple(sum((a[k] for a in args),Q(0)) for k in (0,1))


def mul(a,b):
    terms = [x*y for x,y in product(a,b)]
    return min(terms),max(terms)


def enclosing(outer,inner):
    return outer[0] <= inner[0] <= inner[1] <= outer[1]


def meet(a,b):
    u,v = max(a[0],b[0]),min(a[1],b[1])
    assert u <= v
    return u,v


def mag(a):
    return max(map(abs,a))


def accum(dst,idx,values):
    assert len(idx) == len(values) == len(set(idx))
    for j,val in zip(idx,values):
        dst[j] = add(dst[j],val)


def accum_H(H,idx,local):
    assert set(local) == set(combinations_with_replacement(range(len(idx)),2))
    for (k,l),pair in local.items():
        ij = tuple(sorted((idx[k],idx[l])))
        H[ij] = add(H.get(ij,ZERO),pair)


def cache_geometry(data,index,fields,nodes):
    """Exact material/cached point strain-line domain replay without force or Hessian sums."""
    arr = data.context.target.a
    vertices = fields["corners"]
    a,b = (vertices[j][0]-vertices[0][0] for j in (1,2))
    c,d = (vertices[j][1]-vertices[0][1] for j in (1,2))
    det = a*d-b*c
    assert det > 0
    lam_R = [(c-d)/det,d/det,-c/det]
    lam_Z = [(b-a)/det,-b/det,a/det]
    full_r,full_z = zip(*(data.center[node] for node in nodes))
    values = []
    edges = ((0,1),(1,2),(2,0))
    rational = lambda v:Q.from_float(float(v))
    dot = lambda u,v:sum((x*y for x,y in zip(u,v)),Q(0))
    for iq,L in enumerate(data.rule["volume_barycentric"]):
        sv = [l*(2*l-1) for l in L]+[4*L[j]*L[k] for j,k in edges]
        basis = lambda dl:[(4*l-1)*t for l,t in zip(L,dl)]+[4*(L[j]*dl[k]+L[k]*dl[j]) for j,k in edges]
        br,bz = basis(lam_R),basis(lam_Z)
        rad = sum(l*p[0] for l,p in zip(L,vertices))
        assert rad > 0
        qS = [dot(br,full_r),dot(bz,full_r),dot(br,full_z),dot(bz,full_z),dot(sv,full_r)/rad]
        mesh_grads = arr["shape_grad"][index,iq]
        cr,cz = zip(*[(rational(x),rational(y)) for x,y in mesh_grads])
        old_radius = rational(arr["qr"][index,iq])
        cached = list(map(rational,arr["shape"][iq]))
        assert old_radius > 0
        qO = [dot(cr,full_r),dot(cz,full_r),dot(cr,full_z),dot(cz,full_z),dot(cached,full_r)/old_radius]
        line = [(min(x,y),max(x,y)) for x,y in zip(qO,qS)]
        ad,bc = mul(line[0],line[3]),mul(line[1],line[2])
        determinant = (ad[0]-bc[1],ad[1]-bc[0])
        assert determinant[0] > 0 and line[4][0] > 0
        values.append(mul(determinant,line[4])[0])
    assert len(values) == 64
    return min(values)


def verify_packet(A,data,get,pin,proof,expected):
    i = proof["element"]
    rebuilt,fields,jet,nodes,case,H,prior,radii = A.G._element(data,i)
    assert digest(rebuilt) == digest(proof)
    packet = get(pin["file"],pin["sha256"])
    assert all(packet[k] == v for k,v in expected.items())
    assert all(packet[k] == v for k,v in proof.items())
    assert packet["local_proof_sha256"] == digest(proof)
    assert packet["record_id"] == data.context.record_id and Q(packet["dose_exact"]) == data.context.dose
    assert packet["parameter_cell_sha256"] == data.context.row["parameter_cell_sha256"]
    assert packet["status"] == "CENTERED_BULK_FORCE_COMPONENT_UNREVIEWED"
    assert packet["new_Hessian_integrations"] is False and packet["precision_bits"] == 640
    assert packet["chart"] == fields["chart"]["kind"] and packet["order"] == case["degree"]
    norm = [A.G.E._norm(p,fields) for p in jet["gradient"]]
    t,_ = A.G.M.volumetric_derivative_tail_bounds(case["degree"],Q(case["rho_exact"]),1,5000)
    assert t == Q(case["phi_prime_tail_upper_exact"])
    mass = Q(case["geometry_two_mass_upper_exact"])
    point = packet["point_force"]
    tails = list(map(Q,point["two_mass_tail_upper_exact"]))
    assert tails == [mass*t*n for n in norm]
    assert 0 <= min(tails) <= max(tails) == Q(proof["point_force_tail_max_exact"])
    assert max(tails) <= Q(proof["original_force_tail_budget_exact"])
    min_J = cache_geometry(data,i,fields,nodes)
    J = pack(point["cache_min_J_bound"])
    assert Q(1,10**6) < J[0] <= J[1] <= min_J
    pt = [list(map(pack,point[key])) for key in ("volumetric_truncation","quadratic_strain_defect",
                 "finite_same_Y_cache_correction","normalized_force_total")]
    gd,gi,hs,both,cpr = [list(map(pack,packet[key])) for key in
                        ("affine_gain_derivative_at_Y","affine_gain_increment",
                        "same_original_box_Hessian_times_displacement","centered_force_enclosure",
                        "force_intersection_with_prior")]
    n = len(proof["global_indices"])
    assert all(len(v)==n for v in [tails,*pt,gd,gi,hs,both,cpr,prior])
    if proof["active"] == "0":
        assert gd == gi == [ZERO]*n
    assert len(H) == n*(n+1)//2
    for j in range(n):
        vol,shear,cache,point_force = [v[j] for v in pt]
        assert enclosing(point_force,add(vol,shear,cache,(-tails[j],tails[j])))
        state_bound = sum(mag(H[tuple(sorted((j,k)))])*r for k,r in enumerate(radii))
        assert enclosing(hs[j],(-state_bound,state_bound))
        assert enclosing(gi[j],mul(gd[j],data.gain_bounds))
        assert enclosing(both[j],add(point_force,hs[j],gi[j]))
        assert cpr[j] == meet(both[j],pack(prior[j]))
    return cpr,pt[-1],gd,H,{"element":i,"point_force_entries":n,"cache_site_domains_checked":64,
                             "point_reported_cache_J_lower_exact":str(J[0]),"independent_cache_J_lower_exact":str(min_J)}


def exponent(v):
    d = v.denominator
    assert d & (d-1) == 0
    return d.bit_length()-1


def binary_matrix(array):
    rows = [[Q.from_float(float(v)) for v in row] for row in array]
    e = max(exponent(v) for row in rows for v in row)
    return [[v.numerator << (e-exponent(v)) for v in row] for row in rows],e


def midpoint_radius_scaled(pairs):
    e = max((exponent(v) for pair in pairs for v in pair),default=0)
    mid,rad = [],[]
    for a,b in pairs:
        x = a.numerator << (e-exponent(a))
        y = b.numerator << (e-exponent(b))
        mid.append(x+y)
        rad.append(y-x)
    return mid,rad,e+1


def matvec(R,eR,intervals):
    mid,rad,e = midpoint_radius_scaled(intervals)
    den = 2**(e+eR)
    result = []
    for row in R:
        m = sum(w*v for w,v in zip(row,mid))
        r = sum(abs(w)*v for w,v in zip(row,rad))
        result.append((Q(m-r,den),Q(m+r,den)))
    return result


def signed_RH(R,eR,H,radii,S,eS):
    """Midpoint/radius exact oracle for staged RH and RH·S interval products."""
    keys = list(H)
    mids,rads,e = midpoint_radius_scaled([H[k] for k in keys])
    entries = [[] for _ in range(N)]
    for (i,j),m,r in zip(keys,mids,rads):
        entries[i].append((j,m,r))
        if i != j:
            entries[j].append((i,m,r))
    d = lcm(*(x.denominator for x in radii))
    rs = [x.numerator*(d//x.denominator) for x in radii]
    rfinal,columns = [],[]
    S_cols = list(zip(*S))
    for Rrow in R:
        m,r = [0]*N,[0]*N
        for w,local in zip(Rrow,entries):
            u = abs(w)
            for j,a,b in local:
                m[j] += w*a
                r[j] += u*b
        rfinal.append(Q(sum((abs(x)+y)*d0 for x,y,d0 in zip(m,r,rs)),d*2**(e+eR)))
        coeffs = []
        for s in S_cols:
            mid = sum(x*y for x,y in zip(s,m))
            rad = sum(abs(x)*y for x,y in zip(s,r))
            coeffs.append((Q(mid-rad,2**(e+eR+eS)),Q(mid+rad,2**(e+eR+eS))))
        columns.append(coeffs)
    return rfinal,columns


def write(path,report):
    with path.open("x",encoding="utf8",newline="\n") as file:
        file.write(json.dumps(report,indent=2,default=str)+"\n")
    return sha256(path.read_bytes()).hexdigest()


def compose():
    hashes,read = reader()
    terminal = read(STAGE/"actual/terminal.json",TERMINAL)
    started = read(STAGE/"actual/started.json")
    manifest = read(STAGE/"actual/results.json",MANIFEST)
    assert type(terminal["exit_code"]) is int and terminal["exit_code"] == 0 and terminal["status"] == "FINISHED"
    assert started["pid"] == terminal["pid"] == PID and terminal["results_sha256"] == MANIFEST
    duration = (datetime.fromisoformat(terminal["utc"])-datetime.fromisoformat(started["utc"])).total_seconds()
    assert 0 < manifest["elapsed_seconds"] <= duration+1
    gate = read(RED/"prelaunch_independent_review.json",PRELAUNCH)
    assert gate["status"] == "PASS" and all(v is True for v in gate["checks"].values())
    sys.path[:0] = [str(RED),str(ROOT/"scripts")]
    import reduce_centered_L as A
    pins,protocol,proofs,expected = A._prepared()
    assert manifest["status"] == "ALL_144_CENTERED_BULK_FORCE_PACKETS_UNREVIEWED"
    assert manifest["record_id"] == A.RECORD and manifest["new_Hessian_integrations"] is manifest["root_asserted"] is False
    assert all(started[k] == manifest[k] == v for k,v in expected.items())
    for path,h in pins.items():
        read(path,h,False)
    assert [s["element"] for s in manifest["files"]] == list(range(144))
    assert [s["file"] for s in manifest["files"]] == protocol["actual_expected_files"]
    data = A.G._source_inputs()
    assert data.hashes == protocol["input_hashes"]
    assert digest(data.center) == protocol["point_center_sha256"]
    assert data.context.base["common_gain"] == Q(protocol["exact_base_gain"])
    assert data.context.drive == data.context.dose*Q(protocol["exact_base_gain"])
    assert data.gain_bounds == tuple(map(Q,protocol["closed_gain_shift_exact"]))
    f,point,gain,surface = ([ZERO]*N for _ in range(4))
    Hbulk,Htotal = {},{}
    count,hcount = [0]*N,[0]*N
    locals_ = []

    def counts(idx):
        for j in idx:
            count[j] += 1
            hcount[j] += len(idx)

    for proof,pin in zip(proofs,manifest["files"]):
        i = proof["element"]
        original = f"{A.ORIGINAL}/bulk/element_{i:03d}.json"
        read(original,data.old_files[original],False)
        if i in data.scaled_files:
            prior = data.scaled_files[i]
            read(prior["file"],prior["sha256"],False)
        result,pt,gd,H,local = verify_packet(A,data,read,pin,proof,expected)
        idx = proof["global_indices"]
        accum(f,idx,result)
        accum(point,idx,pt)
        accum(gain,idx,gd)
        accum_H(Hbulk,idx,H)
        counts(idx)
        locals_.append(local)
    Htotal.update(Hbulk)
    old = read(A.B.OLD+"_input_admission.json",A.B.PINS[A.B.OLD+"_input_admission.json"])
    surf_packets = old["expected_packets"]["surface"]
    assert len(surf_packets) == 24
    for i,binding in enumerate(surf_packets):
        p = read(binding["file"],binding["sha256"])
        idx = A.D._packet_meta(data.context,"surface",i,data.domains["surface_cases"][i],p,data.body,data.protocol)
        g = list(map(pack,p["force_integral_minus_cache"]))
        accum(surface,idx,g)
        accum(f,idx,g)
        h = {(v["i"],v["j"]):pack(v["enclosure"]) for v in p["H_integral_minus_cache_upper"]}
        accum_H(Htotal,idx,h)
        counts(idx)
    assert count == data.protocol["force_incidence"] and hcount == data.protocol["H_row_incidence"]
    assert len(Htotal) == 6676 and len(f) == 611
    h_hash = digest([(i,j,Htotal[i,j]) for i,j in sorted(Htotal)])
    previous_path = f"{A.G.SCALED}/reduction_preparation/results.json"
    oldred = read(previous_path,A.G.PIN[previous_path])
    assert h_hash == oldred["H_endpoint_sha256"]

    R,eR = binary_matrix(data.context.R)
    S,eS = binary_matrix(data.context.S)
    norm = Q(max(sum(map(abs,row)) for row in R),2**eR)
    rows = [Q(0)]*N
    for (i,j),v in Htotal.items():
        rows[i] += mag(v)
        if i != j:
            rows[j] += mag(v)
    b = max(x[1] for x in matvec([list(map(abs,row)) for row in R],eR,[(v,v) for v in rows]))
    assert b == Q(oldred["b_upper_exact"]) and norm == Q(oldred["R_infinity_norm_exact"])
    direct,point_R,gain_R,surface_R = [matvec(R,eR,v) for v in (f,point,gain,surface)]
    radii = [Q(data.box["coordinate_radii_exact"][int(do)//2][int(do)%2]) for do in data.context.free]
    state,coupling = signed_RH(R,eR,Hbulk,radii,S,eS)
    gain_term = [mul(v,data.gain_bounds) for v in gain_R]
    separated = [add(pt,(-st,st),g,sf) for pt,st,g,sf in zip(point_R,state,gain_term,surface_R)]
    joint = [meet(a,b) for a,b in zip(direct,separated)]
    keys = list(A.G.T.ORDER)
    ds = [(v[0]-data.context.base[k],v[1]-data.context.base[k]) for k,v in ((k,data.context.cell[k]) for k in keys)]
    k_gain = keys.index("common_gain")
    predictor = [add(p,us,*[mul(add(x,g) if i==k_gain else x,d) for i,(x,d) in enumerate(zip(row,ds))])
                 for p,row,g,us in zip(point_R,coupling,gain_R,surface_R)]
    final = [meet(a,b) for a,b in zip(joint,predictor)]
    force = max(map(mag,final))
    eta0,z0,r,provider = A.D._old_premises(data.context,data.protocol,hashes)
    eta,z = eta0+force,z0+b
    lhs = eta+z*r
    passed = z<1 and lhs<r
    error = eta/(1-z) if passed else None
    curv = A.D._minimum(data.context,hashes,z0,z,r,norm,provider) if passed else None
    omitted = {key:digest(array) for key,array in {
        "direct_whole_box_RF_intervals_exact":direct,"separated_whole_box_RF_intervals_exact":separated,
        "whole_box_RF_intersections_exact":joint,"coupled_predictor_only_RF_intervals_exact":predictor,
        "final_predictor_only_RF_intersections_exact":final,"separated_point_bulk_RF_intervals_exact":point_R,
        "separated_state_bulk_only_radius_exact":state,"separated_gain_derivative_RF_intervals_exact":gain_R,
        "separated_gain_term_exact":gain_term,"separated_retained_surface_RF_intervals_exact":surface_R,
        "coupled_RH_bulk_S_intervals_exact":coupling,"exact_closed_parameter_offsets":ds}.items()}
    includes_zero = sum(a <= 0 <= b for a,b in final)
    lower = max(Q(0) if a <= 0 <= b else min(abs(a),abs(b)) for a,b in final)
    scalars = {
        "record_id":data.context.record_id,"component_count":168,"centered_bulk_components":144,
        "retained_surface_components":24,"force_entries":611,"upper_Hessian_entries":6676,
        "F_endpoint_sha256":digest(f),"uncentered_H_endpoint_sha256":h_hash,
        "bulk_only_H_endpoint_sha256":digest([(i,j,Hbulk[i,j]) for i,j in sorted(Hbulk)]),
        "predictor_parameter_order":keys,
        "coupled_predictor_only_RF_uniform_upper_exact":max(map(mag,predictor)),
        "direct_RF_uniform_upper_exact":max(map(mag,direct)),
        "separated_RF_uniform_upper_exact":max(map(mag,separated)),
        "whole_box_RF_uniform_upper_exact":max(map(mag,joint)),
        "RF_uniform_upper_exact":force,"predictor_only_RF_uniform_lower_exact":lower,
        "R_infinity_norm_exact":norm,"b_upper_exact":b,"cached_eta_upper_exact":eta0,"cached_Z_upper_exact":z0,
        "eta_upper_exact":eta,"z_upper_exact":z,"selfmap_lhs_exact":lhs,"normal_radius_exact":r,
        "sufficient_conditions_pass":passed,"conditional_root_error_upper_exact":error,
        "conditional_curvature":curv,"same_point_cell_and_cube":True,"root_assertion_independently_reviewed":False,
        "status":"CONDITIONAL_UNREVIEWED_CENTERED_ORIGINAL_CUBE_SUFFICIENT_TEST_PASSES" if passed
            else "UNREVIEWED_CENTERED_ORIGINAL_CUBE_SUFFICIENT_TEST_UNRESOLVED",
    }
    own = Path(__file__).relative_to(ROOT)
    read(own,parse=False)
    for name,h in tuple(hashes.items()):
        read(name,h,False)
    report = {"status":"PASS_ACTUAL_CENTERED_PACKET_AND_COMPOSITION_REPLAY",
              "scope":"Independent endpoint and signed midpoint/radius arithmetic, original box. Point force/gain evaluation and original/scaled Hessian moment/AD premises remain shared. Analytic site checks evaluated strain/J bounds only. No fresh nonlinear solve, force integral, or cache force sum.",
              "checks":{
                  "actual_144_manifest_and_terminal_match_supplied_hashes_and_PID":True,
                  "all_144_points_proof_center_base_drive_original_charts_order_rho_tail":True,
                  "exact_signed_local_F_parts_and_state_gain_sums_contained":True,
                  "independent_exact_all9216_point_site_strain_line_J_bounds":True,
                  "all_local_intersections_and_611_6676_exact_signed_global_assembly":True,
                  "three_compositions_with_independent_exact_midpoint_radius_matrix_arithmetic":True,
                  "predictor_only_eta_and_same_full_box_H_Z_and_strict_conditional_gate":True,
              },"source_and_inputs":dict(sorted(hashes.items())),"local_metadata":locals_,
              "expected_reducer_array_digests":omitted,"expected_reducer_scalar_fields":scalars,
              "final_predictor_RF_zero_containing_coordinates":includes_zero,
              "preconditioned_point_bulk_bound_exact":max(map(mag,point_R)),"actual_new_integrations":False}
    component_sha = write(RECEIPT,report)
    gate = {key:v for key,v in report.items() if key != "local_metadata"}
    gate["support_file"] = RECEIPT.relative_to(ROOT).as_posix()
    gate["support_sha256"] = component_sha
    actual_sha = write(ACTUAL_GATE,gate)
    print(json.dumps({"status":report["status"],"component_sha":component_sha,"review_sha":actual_sha,
                      "point_bulk":float(report["preconditioned_point_bulk_bound_exact"]),"RF":float(force),
                      "ratio":float(lhs/r),"pass":passed,"RF_includes_0":includes_zero}),flush=True)


def finalize(args):
    hashes,read = reader()
    receipt = read(RECEIPT,args.component_sha)
    for p,h in receipt["source_and_inputs"].items():
        read(p,h,False)
    adm = read(RED/"input_admission.json",args.admission_sha)
    term = read(RED/"terminal.json",args.terminal_sha)
    start = read(RED/"started.json")
    actual = read(RED/"results.json",args.result_sha)
    assert term["status"] == "FINISHED" and type(term["exit_code"]) is int and term["exit_code"] == 0
    assert term["pid"] == start["pid"] and term["results_sha256"] == args.result_sha
    assert start["admission_sha256"] == actual["admission_sha256"] == args.admission_sha
    assert adm["local_centered_payloads_opened"] is False
    assert adm["centered_terminal_sha256"] == TERMINAL and adm["centered_manifest_sha256"] == MANIFEST
    assert adm["prelaunch_review_sha256"] == PRELAUNCH
    assert adm["expected_centered_packets"] == read(STAGE/"actual/results.json",MANIFEST)["files"]
    for p,h in actual["source_inputs_and_packets"].items():
        read(p,h,False)
    for k,v in receipt["expected_reducer_scalar_fields"].items():
        assert actual[k] == v, k
    for k,h in receipt["expected_reducer_array_digests"].items():
        assert digest(actual[k]) == h, k
    report = {"status":"PASS_ACTUAL_CENTERED_ORIGINAL_CUBE_REDUCTION_REVIEW",
              "checks":{"terminal_and_manifest_admission_current":True,
                        "all_reducer_force_matrix_array_digests_and_exact_scalars_match":True,
                        "all_bound_inputs_current":True},
              "component_replay_sha256":args.component_sha,"reducer_results_sha256":args.result_sha,
              "admission_sha256":args.admission_sha,"terminal_sha256":args.terminal_sha,
              "decision":receipt["expected_reducer_scalar_fields"],
              "final_predictor_RF_zero_containing_coordinates":receipt["final_predictor_RF_zero_containing_coordinates"],
              "scope":"Independent exact endpoint replay for the frozen original L cube only. Uses reviewed shared local moment, AD and cache kernels as premises.",
              "input_hashes":dict(sorted(hashes.items()))}
    value = write(RED_GATE,report)
    print(json.dumps({"review_sha256":value,"status":report["status"],
                      "sufficient_conditions_pass":actual["sufficient_conditions_pass"]}),flush=True)


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Replay requires enabled assertions")
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="action",required=True)
    sub.add_parser("compose")
    final = sub.add_parser("finalize")
    for name in ("component","admission","terminal","result"):
        final.add_argument(f"--{name}-sha",required=True)
    args = parser.parse_args()
    compose() if args.action == "compose" else finalize(args)
