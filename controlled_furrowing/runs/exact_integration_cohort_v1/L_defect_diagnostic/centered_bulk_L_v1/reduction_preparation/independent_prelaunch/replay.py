"""Independent prelaunch audit of the three centered-L force bounds."""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement, product
from io import StringIO
import json
from pathlib import Path
from random import Random
import sys
from contextlib import redirect_stdout
from unittest.mock import patch

OUT = Path(__file__).resolve().parent
JOB = OUT.parent
ROOT = JOB.parents[4]
GATE = JOB / "prelaunch_independent_review.json"
EXPECTED = {
    "source_sha256":"6e35f18b3d6f51360667074b312ca094f1aa544a17c06be8548cdd9f6c7bd1cc",
    "protocol_sha256":"d57e6d58cc3520560f108768a5063b01603c9ebac66f042bd1622815e78f2ddc",
    "preflight_sha256":"ac46b00d8f8d1c9b566685950d6bafe993efaae5c22f08c323eae7a0018d760f",
}
DESIGN = {"README.md":"5ede3ed1b759442f6b6fbfcd25bcb1819f3de8b69acc4ebafa0ea17c62506c6b",
           "receipt.json":"1ba6f141868f887f4e3dcdaea5c313de4e3591e5f9b81a7375d03fd3d048a3dd"}


def reject(action,text):
    try:
        action()
    except ValueError as error:
        assert text in str(error), str(error)
    else:
        raise AssertionError("invalid synthetic fixture admitted")


def mul(a,b):
    candidates = [x*y for x,y in product(a,b)]
    return min(candidates),max(candidates)


def isum(values):
    entries = tuple(values)
    return tuple(sum((v[k] for v in entries),Q(0)) for k in (0,1))


def direct_RH(R,H):
    n = len(R)
    return [[isum(mul((Q(w),Q(w)),H.get(tuple(sorted((i,j))),(Q(0),Q(0))))
                  for i,w in enumerate(row)) for j in range(n)] for row in R]


def linear_oracles(A):
    rng = Random(623190)
    count = 0
    for n in (1,2,3,4):
        for _ in range(12):
            R = [[rng.randint(-24,24)/2**rng.randrange(6) for _ in range(n)] for _ in range(n)]
            slope = [[rng.randint(-12,12)/2**rng.randrange(8) for _ in range(5)] for _ in range(n)]
            rad = [Q(rng.randint(1,11),rng.randrange(2,13)) for _ in range(n)]
            H = {}
            for i,j in combinations_with_replacement(range(n),2):
                if rng.random() < 0.25:
                    continue
                a = Q(rng.randrange(-24,15),16)
                H[i,j] = a,a+Q(rng.randrange(12),32)
            RH = direct_RH(R,H)
            want_perm = [sum(max(map(abs,v))*d for v,d in zip(row,rad)) for row in RH]
            want_slope = [[isum(mul(pair,(Q(s[k]),Q(s[k]))) for pair,s in zip(row,slope)) for k in range(5)] for row in RH]
            perm,staged = A._bulk_RH_effects(R,H,rad,slope)
            assert perm == want_perm and staged == want_slope
            count += 1
    # Exact nonlinear scalar corner replay for the relaxed common-parameter
    # representation. Here H is uncertain symmetric and varies independently.
    R = [[1.5,-0.5],[-1.25,0.25]]
    S = [[1.0,-1.5,.25,0.0,-1.0],[.5,1.0,-.25,1.0,2.0]]
    H = {(0,0):(Q(-2),Q(-1)),(0,1):(Q(-1),Q(1,2)),(1,1):(Q(1),Q(2))}
    gain = [(Q(-3),Q(2)),(Q(1),Q(3))]
    point = [(Q(-1),Q(1,2)),(Q(1),Q(3))]
    surf = [(Q(-1,2),Q(1,2)),(Q(-3,2),Q(-1))]
    dp = [(Q(-1,3),Q(1,9)),(Q(0),Q(1,5)),(Q(-1,7),Q(3,7)),(Q(0),Q(1,2)),(Q(-3,4),Q(1,4))]
    _,Scols = A._bulk_RH_effects(R,H,[Q(1),Q(2)],S)
    enclosure = A._coupled_predictor_force(point,Scols,gain,surf,dp,4)
    for dtheta in product(*dp):
        dY = [sum(Q(v)*d for v,d in zip(row,dtheta)) for row in S]
        for Ha,Hb,Hc in product(*H.values()):
            Hx = [Ha*dY[0]+Hb*dY[1],Hb*dY[0]+Hc*dY[1]]
            for i,row in enumerate(R):
                state = sum(Q(w)*v for w,v in zip(row,Hx))
                for k,p,s in product(gain[i],point[i],surf[i]):
                    value = p+s+state+k*dtheta[4]
                    assert enclosure[i][0] <= value <= enclosure[i][1]
    # Same gain in both increments can cancel exactly. Independent
    # relaxation would create a nonzero interval unwarranted here.
    z = (Q(0),Q(0))
    actual = A._coupled_predictor_force([z],[[z,z,z,z,(Q(2),Q(2))]],
                                       [(Q(-2),Q(-2))],[z],dp,4)
    assert actual == [z]
    return count


def local_packet_oracle(A,expected):
    pack = lambda a,b:{"lower_exact":str(a),"upper_exact":str(b)}
    proof = {"element":11,"global_indices":[4,1],"active":"1","point_force_tail_max_exact":"1/16",
             "original_force_tail_budget_exact":"1/8","free_local":[2,3]}
    row = {"parameter_cell_sha256":"d"*64}
    p = {**expected,**proof,"status":"CENTERED_BULK_FORCE_COMPONENT_UNREVIEWED","record_id":A.RECORD,
         "dose_exact":"3/2","parameter_cell_sha256":row["parameter_cell_sha256"],"precision_bits":640,
         "local_proof_sha256":A.D._digest(proof),"new_Hessian_integrations":False}
    p["point_force"] = {"volumetric_truncation":[pack(-1,-1),pack(1,1)],
                        "quadratic_strain_defect":[pack(0,0),pack(-1,-1)],
                        "finite_same_Y_cache_correction":[pack(0,0),pack(-1,-1)],
                        "normalized_force_total":[pack("-17/16","-15/16"),pack(-1,-1)],
                        "two_mass_tail_upper_exact":["1/16","0"],"cache_min_J_bound":pack("1/2",1)}
    p["affine_gain_derivative_at_Y"] = [pack(-1,2),pack(2,3)]
    p["affine_gain_increment"] = [pack(-1,"1/2"),pack("-3/2","3/4")]
    p["same_original_box_Hessian_times_displacement"] = [pack(-2,2),pack(-2,2)]
    p["centered_force_enclosure"] = [pack(-5,3),pack(-5,3)]
    old = [pack(-4,4),pack(-3,4)]
    p["force_intersection_with_prior"] = [pack(-4,3),pack(-3,3)]
    H = {(0,0):(Q(1),Q(1)),(0,1):(Q(-2),Q(-1)),(1,1):(Q(-1),Q(1))}
    f,center,gain = A._centered_packet(p,proof,expected,row,Q(3,2),old,H,[Q(1,8),Q(1,4)],(Q(-1,2),Q(1,4)))
    assert f == [(Q(-4),Q(3)),(Q(-3),Q(3))]
    assert center[0] == (Q(-17,16),Q(-15,16)) and gain[1] == (Q(2),Q(3))
    verify = lambda packet:A._centered_packet(packet,proof,expected,row,Q(3,2),old,H,[Q(1,8),Q(1,4)],(Q(-1,2),Q(1,4)))
    reject(lambda:verify({**p,"same_original_box_Hessian_times_displacement":[pack(0,0),pack(-2,2)]}),"same-box")
    reject(lambda:verify({**p,"affine_gain_increment":[pack(0,0),p["affine_gain_increment"][1]]}),"same-box")
    reject(lambda:verify({**p,"force_intersection_with_prior":[pack(-5,3),pack(-3,3)]}),"intersection changed")


def admission_oracle(A,plan,expected):
    files = [{"element":i,"file":path,"sha256":sha256(path.encode()).hexdigest()}
             for i,path in enumerate(plan["actual_expected_files"])]
    manifest = {**expected,"record_id":A.RECORD,"status":"ALL_144_CENTERED_BULK_FORCE_PACKETS_UNREVIEWED",
                "new_Hessian_integrations":False,"root_asserted":False,"files":files}
    terminal_sha,manifest_sha,review_sha = "a"*64,"b"*64,"c"*64
    terminal = {"status":"FINISHED","exit_code":0,"pid":4242,"results_sha256":manifest_sha}
    touched,outputs = [],[]
    def real_path_only(path,pin,hashes):
        assert path in (f"{A.STAGE}/actual/terminal.json",f"{A.STAGE}/actual/results.json")
        is_term = path.endswith("terminal.json")
        assert pin == (terminal_sha if is_term else manifest_sha)
        touched.append(path)
        hashes[path] = pin
        return terminal if is_term else manifest
    def fake_write(path,value):
        assert path.name == "input_admission.json"
        outputs.append(value)
        return "d"*64
    with patch.object(A,"_prepared",return_value=({},plan,[],expected)), \
            patch.object(A,"_check_gate",side_effect=lambda h,b:b.update(gate=h)), \
            patch.object(A.D,"_bound",side_effect=real_path_only), \
            patch.object(A.D,"_write",side_effect=fake_write),redirect_stdout(StringIO()):
        A.admit(review_sha,terminal_sha,manifest_sha)
        assert len(touched) == 2 and len(outputs) == 1
        assert outputs[0]["expected_centered_packets"] == files
        assert outputs[0]["local_centered_payloads_opened"] is False
        terminal["exit_code"] = 2
        reject(lambda:A.admit(review_sha,terminal_sha,manifest_sha),"terminal")
        terminal["exit_code"] = 0
        terminal["results_sha256"] = "0"*64
        reject(lambda:A.admit(review_sha,terminal_sha,manifest_sha),"binding")
        terminal["results_sha256"] = manifest_sha
        manifest["files"] = files[::-1]
        reject(lambda:A.admit(review_sha,terminal_sha,manifest_sha),"canonical")
        manifest["files"] = files
        manifest["root_asserted"] = True
        reject(lambda:A.admit(review_sha,terminal_sha,manifest_sha),"reviewed preparation")
        assert len(outputs) == 1


def run():
    if not __debug__:
        raise RuntimeError("Replay requires enabled assertions")
    hashes = {}

    def read(name,digest=None,parse=True):
        path = (ROOT/name).resolve()
        assert path.is_relative_to(ROOT)
        assert not path.is_relative_to(JOB.parent/"actual")
        raw = path.read_bytes()
        value = sha256(raw).hexdigest()
        assert digest is None or digest == value, name
        assert hashes.setdefault(path.relative_to(ROOT).as_posix(),value) == value
        return json.loads(raw) if parse else raw

    rel = JOB.relative_to(ROOT).as_posix()
    src = read(f"{rel}/reduce_centered_L.py",EXPECTED["source_sha256"],False)
    assert src == read(f"{rel}/reduce_centered_L_frozen.py",EXPECTED["source_sha256"],False)
    proto,pre = [read(f"{rel}/{key}.json",EXPECTED[f"{key}_sha256"]) for key in ("protocol","preflight")]
    assert pre["status"] == "PASS_CENTERED_L_COMPOSITION_PREPARATION_ONLY"
    assert len(pre["checks"]) == pre["check_count"] == 21 and all(v is True for v in pre["checks"].values())
    assert pre["source_sha256"] == proto["source_sha256"] == EXPECTED["source_sha256"]
    assert pre["protocol_sha256"] == EXPECTED["protocol_sha256"]
    assert pre["actual_reduction_run"] is proto["actual_reduction_run"] is proto["actual_centered_packet_payloads_opened"] is False
    assert proto["component_count"] == 168 and proto["force_entries"] == 611 and proto["upper_Hessian_entries"] == 6676
    for path,digest in {**proto["input_hashes"],**proto["failed_attempts"]}.items():
        read(path,digest,False)
    for file,digest in DESIGN.items():
        read(JOB.parent/"composition_design"/file,digest,False)
    original_open = Path.open
    def guarded_open(path,*args,**kwargs):
        assert not path.resolve().is_relative_to(JOB.parent/"actual"), path
        return original_open(path,*args,**kwargs)

    sys.path.insert(0,str(JOB))
    import reduce_centered_L as A
    with patch.object(Path,"open",guarded_open):
        prep_hash,plan,proofs,expected = A._prepared()
        assert prep_hash == proto["input_hashes"] and len(proofs) == 144
        assert len(A.tests()) == 20
        data = A.G._source_inputs()
        assert data.hashes == plan["input_hashes"]
        c = data.context
        offsets = [(lo-c.base[k],hi-c.base[k]) for k in A.G.T.ORDER for lo,hi in (c.cell[k],)]
        assert len(offsets) == 5 and all(lo <= 0 <= hi for lo,hi in offsets)
        assert offsets[A.G.T.ORDER.index("common_gain")] == data.gain_bounds
        assert len(c.S) == len(c.free) == 611
        for i,free in enumerate(c.free):
            node,axis = divmod(int(free),2)
            lo,hi = isum(mul((Q.from_float(float(v)),)*2,d) for v,d in zip(c.S[i],offsets))
            a,b = map(Q,data.box["predictor_hulls_exact"][i])
            assert (lo,hi) == (a,b)
            rr = Q(data.box["coordinate_radii_exact"][node][axis])
            assert -rr <= lo-c.old_radius <= hi+c.old_radius <= rr
        cached = {}
        eta,z,r,origin = A.D._old_premises(c,data.protocol,cached)
        source = f"{A.G.SCALED}/reduction_preparation/results.json"
        previous = read(source,A.G.PIN[source])
        b = Q(previous["b_upper_exact"])
        v = 1-z-b-eta/r
        assert v > 0 and Q(previous["z_upper_exact"]) == z+b
        assert not A.D._strict_test(eta,z,v*r,b,r)[3]
        assert A.D._strict_test(eta,z,v*r/2,b,r)[3]
        existing = A.D._minimum(c,cached,z,z+b,r,Q(previous["R_infinity_norm_exact"]),origin)
        assert existing["strict_local_minimum"] is True
        assert Q(existing["normalized_curvature_lower_exact"]) == (1-z-b)/Q(previous["R_infinity_norm_exact"])
        for path,digest in cached.items():
            read(path,digest,False)
        random_count = linear_oracles(A)
        local_packet_oracle(A,expected)
        admission_oracle(A,plan,expected)

    own = Path(__file__).relative_to(ROOT)
    read(own,parse=False)
    for path,digest in tuple(hashes.items()):
        read(path,digest,False)
    report = {
        **EXPECTED,"status":"PASS",
        "scope":"Independent prelaunch audit of three signed force enclosures on frozen L. Uses prior Hessians and prepared states only. No live centered terminal, manifest, payload, reduction or root claim.",
        "checks":{
            "frozen_source_protocol_preflight_parent_checks_and_design_current":True,
            "true_exit_and_ordered144_manifest_admission_before_packet_read_verified_synthetically":True,
            "same_signed_RH_first_symmetric_H_bulk_radii_and_5_parameter_matrix_dense_Fraction_oracle":True,
            "coupled_same_gain_column_and_all_parameter_corners_synthetic_oracle":True,
            "all_611_predictor_hulls_plus_corrector_cube_fit_original_box_and_same_gain_C":True,
            "same_cached_residual_gate_uses_predictor_only_upper_and_full_total_box_derivative":True,
            "signed_local_F_parts_and_mixed_H_state_gain_interval_checks":True,
            "strict_test_and_conditional_same_R_cached_SPD_homotopy":True,
            "live_centered_input_open_blocked_throughout_replay":True,
        },
        "small_dense_RH_fixtures":random_count,
        "whole_box_separated_state_policy":"Source uses the original symmetric radii. The design's optional asymmetric predictor-only enclosure is not the separated whole-box strategy.",
        "joint_predictor_scope":"Centered RH_bulk S and R-gain derivative share the same delta gain, along the parameter-matched predictor only. Full-box total H still controls z.",
        "local_packet_scope":"Point/centered force component overlap is a consistency diagnostic under the shared prelaunched numerical-kernel premise. Local state/gain arithmetic is required to enclose exact reconstructed signed endpoints.",
        "force_room_at_retained_derivative_exact":str(v*r),
        "input_hashes":dict(sorted(hashes.items())),
        "review_source_file":own.as_posix(),"review_source_sha256":hashes[own.as_posix()],
        "actual_centered_payloads_opened":False,"actual_reduction_run":False,
    }
    if "--freeze" in sys.argv:
        receipt = OUT/"replay.json"
        assert not receipt.exists() and not GATE.exists()
        with receipt.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(report,indent=2)+"\n")
        report["support_file"] = receipt.relative_to(ROOT).as_posix()
        report["support_sha256"] = sha256(receipt.read_bytes()).hexdigest()
        with GATE.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":report["status"],"checks":len(report["checks"]),
                      "small_dense_RH_fixtures":random_count,"force_room":float(v*r)}),flush=True)


if __name__ == "__main__":
    run()
