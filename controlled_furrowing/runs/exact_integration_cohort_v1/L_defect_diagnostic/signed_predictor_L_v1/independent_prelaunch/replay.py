"""Prelaunch check of the signed cached-plus-integral predictor composition."""
from base64 import b64decode
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
from itertools import product
import json
from pathlib import Path
from random import Random
import sys
from unittest.mock import patch

from flint import arb, arb_mat, ctx
import numpy as np

OUT = Path(__file__).resolve().parent
STAGE = OUT.parent
ROOT = STAGE.parents[3]
PREP = STAGE / "preparation"
SOURCE = "scripts/exact_signed_predictor_L.py"
EXPECTED = {
    "source_sha256":"00cebb37ed7cd23152ee46a659ff9878a32f62733faa2964d3b051fa7678dd1b",
    "protocol_sha256":"4e7152f01104da4c5f09a8134f54ac3228ae6467b7ae34d862f50265b5a003af",
    "preflight_sha256":"e91f11f8f7f82aae6670fe4215750df2506e59f7501205c8faeec76577ce99e3",
}


def reject(call,word):
    try:
        call()
    except ValueError as error:
        assert word in str(error),str(error)
    else:
        raise AssertionError("malformed signed predictor fixture accepted")


def corner_oracle(A):
    rng = Random(192276)
    interval = lambda: sorted((Q(rng.randint(-9,5),16),Q(rng.randint(-5,12),8)))
    count = 0
    for _ in range(48):
        offsets = [(Q(-rng.randint(0,10),9),Q(rng.randint(0,10),17)) for _ in range(5)]
        g = [tuple(interval()) for _ in range(3)]
        d = [tuple(interval()) for _ in range(3)]
        surface = [tuple(interval()) for _ in range(3)]
        cached = [[tuple(interval()) for _ in range(5)] for j in range(3)]
        defect = [[tuple(interval()) for _ in range(5)] for j in range(3)]
        gain = [tuple(interval()) for _ in range(3)]
        remainder = [Q(rng.randint(0,12),23) for _ in range(3)]
        gain_col = rng.randrange(5)
        origins = [A._add(a,b,c) for a,b,c in zip(g,d,surface)]
        coeffs = A.combine_slopes(cached,defect,gain,gain_col)
        result = A.affine_rows(origins,coeffs,offsets,remainder)
        for i,(a,b) in enumerate(result):
            expect = []
            for theta in product(*offsets):
                def end(direction):
                    part = g[i][direction]+d[i][direction]+surface[i][direction]
                    for j,t in enumerate(theta):
                        terms = (cached[i][j],defect[i][j],gain[i]) if j==gain_col else (cached[i][j],defect[i][j])
                        part += sum(t*term[direction if t>=0 else 1-direction] for term in terms)
                    return part
                expect.append((end(0),end(1)))
            assert (a,b) == (min(x for x,y in expect)-remainder[i],max(y for x,y in expect)+remainder[i])
        count += 1
    return count


def q_cache_affine_oracle(A):
    """Independent physical coefficients for a positive synthetic bulk/surface site."""
    to_arb = A.P.core.A

    def encl(ball,value):
        a,b = A._arb_bounds(ball)
        assert a <= value <= b

    qbulk = [to_arb(v) for v in (1,Q(1,8),0,1,1)]
    dose = Q(4,3)
    gs = (Q(5,7),Q(9,8))
    drives = [A.P.core.bulk_local(qbulk,to_arb(1),to_arb(5000),to_arb(dose*g))[0] for g in gs]
    for j,expected in enumerate((Q(0),dose*(gs[1]-gs[0])/8,Q(0),dose*(gs[1]-gs[0]),Q(0))):
        encl(drives[1][j]-drives[0][j],expected)

    # Cached [R,x,z] are raw linear strains. Scale fixes the old binary
    # normalization of Js; nonzero tangential components make all rows matter.
    scale = Q(1,7)
    raw = (Q(2),Q(3),Q(4))
    material = ((Q(11,5),Q(1,3)),(Q(9,7),Q(8,9)))
    J = scale*raw[0]*5
    dJ = [5*scale,scale*raw[0]*raw[1]/5,scale*raw[0]*raw[2]/5]
    grades = []
    for K,T in material:
        grad,_,_ = A.P.core.surface_local(list(map(to_arb,raw)),to_arb(K),to_arb(T),to_arb(scale))
        for v,d in zip(grad,dJ):
            encl(v,(K*(J-1)+T)*d)
        grades.append(grad)
    dK,dT = tuple(y-x for x,y in zip(*material))
    for a,b,w in zip(*grades,dJ):
        encl(b-a,(dK*(J-1)+dT)*w)


def extract_mock(A,context,args,offsets):
    """Check the 611-row boundary with stubbed point/cache paths and remainders."""
    aq = A.P.core.A
    n = 611
    fake = {"n":n,"S":context.S,"cache":[{"kind":"bulk"}]*9216 + [{"kind":"surface"}]*288,
            "base":[context.base[k] for k in A.P.NAMES],"dose":aq(context.dose),
            "Z0":aq(Q(1,4)),"gR":arb_mat(n,1,[aq(Q(j,256)) for j in range(n)]),
            "affine_M":arb_mat(n,5,[aq(Q(j-300,256)) for j in range(n*5)])}
    calls = []

    def variation(point,cell,radius,rem):
        assert point is fake and radius.is_zero() and rem is True and cell["method"] == "affine"
        assert all(A.D._contains(A._arb_bounds(v),p) for v,p in zip(cell["delta"],offsets))
        calls.append(True)
        return [aq(Q(j,1024)) for j in range(n)],aq(Q(1,2))

    with patch.object(A.P,"build_point",return_value=fake) as builder, \
            patch.object(A.P,"variation",side_effect=variation):
        shape,M,bounds,J,Z = A._point_and_remainder(context,args,offsets)
        assert calls == [True] and builder.call_args.args == (args,)
        assert shape[35] == (Q(35,256),)*2 and M[143][4] == (Q(419,256),)*2
        assert bounds[-1] == Q(610,1024) and J == (Q(1,2),)*2 and Z == (Q(1,4),)*2
        wrong = context.S.copy()
        wrong[0,0] += 1
        fake["S"] = wrong
        reject(lambda:A._point_and_remainder(context,args,offsets),"same frozen decoded")
        fake["S"] = context.S
        base = fake["base"]
        fake["base"] = list(reversed(base))
        reject(lambda:A._point_and_remainder(context,args,offsets),"same frozen decoded")
        fake["base"] = base
        foreign = list(offsets)
        foreign[-1] = (Q(-1),Q(1))
        reject(lambda:A._point_and_remainder(context,args,foreign),"enclose the exact offsets")
        assert calls == [True]
    with patch.object(A.P,"build_point",return_value=fake), \
            patch.object(A.P,"variation",return_value=([aq(0)]*611,aq(Q(1,10**6)))):
        reject(lambda:A._point_and_remainder(context,args,offsets),"predictor path and remainder")


def run():
    if not __debug__:
        raise RuntimeError("Review requires enabled assertions")
    inputs = {}

    def read(name,h=None,parse=True):
        path = (ROOT/name).resolve()
        assert path.is_relative_to(ROOT) and not path.is_relative_to(STAGE/"actual")
        raw = path.read_bytes()
        digest = sha256(raw).hexdigest()
        assert h is None or h == digest, name
        assert inputs.setdefault(path.relative_to(ROOT).as_posix(),digest) == digest
        return json.loads(raw) if parse else raw

    rel = PREP.relative_to(ROOT).as_posix()
    source = read(SOURCE,EXPECTED["source_sha256"],False)
    assert source == read(f"{rel}/exact_signed_predictor_L_frozen.py",EXPECTED["source_sha256"],False)
    p,pre = [read(f"{rel}/{key}.json",EXPECTED[f"{key}_sha256"]) for key in ("protocol","preflight")]
    assert pre["status"] == "PASS_SIGNED_PREDICTOR_PREPARATION_ONLY"
    assert len(pre["checks"]) == 17 and all(v is True for v in pre["checks"].values())
    assert pre["source_sha256"] == p["source_sha256"] == EXPECTED["source_sha256"]
    assert pre["protocol_sha256"] == EXPECTED["protocol_sha256"]
    assert p["actual_cached_point_or_variation_run"] is pre["actual_cached_9504_site_sums"] is False
    for path,h in {**p["input_hashes"],**p["failed_attempts"]}.items():
        read(path,h,False)
    attempt = read(f"{rel}/development_revision_0/attempt.json")
    assert attempt["status"] == "REJECTED_BEFORE_PREFLIGHT_INPUT_LOAD" and attempt["cached_site_sums_started"] is False
    assert p["failed_attempts"][f"{rel}/development_revision_0/exact_signed_predictor_L.py"] == attempt["source_sha256"]
    sys.path.insert(0,str(ROOT/"scripts"))
    import exact_signed_predictor_L as A
    ctx.prec = 640
    A.P.core.A.cache_clear()
    blocked = AssertionError("actual 9504 cached sites or variation attempted during prelaunch")
    with patch.object(A.P,"build_point",side_effect=blocked), \
            patch.object(A.P,"variation",side_effect=blocked), \
            patch.object(A.P.core,"assemble_point",side_effect=blocked), \
            patch.object(A.P.core,"local_terms",side_effect=blocked):
        bindings,c,original,old,bulk,surf,gain,bulkS,delta,b,z,r,old_Z,provider,args = A._inputs()
        assert bindings == p["input_hashes"]
        assert c.row == original["row"]
        assert r == Q(p["radius_exact"]) == c.old_radius and z == Q(p["unchanged_z_upper_exact"]) == old_Z+b
        assert Q(p["strict_signed_force_budget_exact"]) == r*(1-z)
        assert tuple(A.P.NAMES) == tuple(A.C.ORDER) == tuple(p["parameter_order"])
        assert args == p["point_kernel_arguments"]
        assert (A.P.ROOT/args["point_certificate_file"]).resolve() == (ROOT/c.row["point_cache_file"]).resolve()
        assert Path(args["predictor_file"]).is_absolute()
        assert (A.P.OUT/args["predictor_file"]).resolve() == (ROOT/c.row["predictor_file"]).resolve()
        assert args["dose_exact"] == str(c.dose) and c.drive == c.dose * c.base["common_gain"]
        assert c.base == {k:Q(v) for k,v in args["base_parameters_exact"].items()}
        assert delta == [tuple(map(Q,v)) for v in p["same_exact_closed_offsets"]]
        assert all(o <= 0 <= u for o,u in delta)
        packet = read(c.row["point_cache_file"],c.row["point_cache_sha256"])
        pars = packet["parameters"]
        assert Q(pars["mu"]) == 1 and Q(pars["bulk"]) == 5000
        for name in A.P.NAMES[:4]:
            assert Q.from_float(float(pars[name])) == c.base[name]
        binary = b64decode(packet["frozen_binary_cache"]["data"],validate=True)
        assert sha256(binary).hexdigest() == c.row["frozen_binary_cache_sha256"]
        with np.load(BytesIO(binary),allow_pickle=False) as saved:
            assert saved["u"].tobytes() == c.u.tobytes()
            assert saved["preconditioner"].tobytes() == c.R.tobytes()
            assert saved["free"].tobytes() == c.free.tobytes()
            assert all(np.isfinite(saved[name]).all() and np.all(saved[name]>0)
                       for name in ("energy_scale","volume_weights","qr","surface_0_weights","surface_1_weights",
                                    "surface_0_r0","surface_1_r0","surface_0_l0","surface_1_l0"))
        with np.load(A.P.OUT/args["predictor_file"],allow_pickle=False) as predicted:
            assert predicted["S"].tobytes() == c.S.tobytes()
        assert sum(f"{A.ORIGINAL}/surface/" in key for key in bindings) == 24
        assert len(bulk) == len(surf) == len(gain) == len(bulkS) == 611
        assert len(A.tests()) == 12

    q_cache_affine_oracle(A)
    extract_mock(A,c,args,delta)
    random_count = corner_oracle(A)
    if A.strict_bound(r*(1-z),z,r)[1]:
        raise AssertionError("self-map equality promoted")
    assert A.strict_bound(r*(1-z)/2,z,r)[1]
    own = Path(__file__).relative_to(ROOT).as_posix()
    read(own,parse=False)
    for path,digest in tuple(inputs.items()):
        read(path,digest,False)
    report = {
        **EXPECTED,"status":"PASS","record_id":c.record_id,
        "scope":"Signed cached-plus-defect predictor residual preparation for frozen L and C. Shared cached kernels reviewed, no actual cached assembly/variation, exact FE integration, nonlinear solve or root claim.",
        "checks":{
            "frozen_input_source_preflight_history_and_retained_surface_arrays":True,
            "verified_exact_binary_u_R_S_free_order_positive_measures_and_exact_five_base_parameters":True,
            "fixed_cached_mu_kappa_same_normalization_and_exact_active_gain_times_dose":True,
            "authoritative_cache_join_relative_to_same_root_and_absolute_predictor_join":True,
            "synthetic_affine_at_fixed_Y0_bulk_and_surface_physical_derivative_oracles":True,
            "theta_first_nodal_path_uses_affine_predictor_only_zero_radius_remainder_then_common_signed_columns":True,
            "independent_5_parameter_exact_signed_coefficient_corner_oracle":True,
            "mixed_decoded_predictor_base_cell_and_bad_cache_J_mock_rejected":True,
            "strict_equality_unresolved_and_prior_whole_box_z_r_unchanged":True,
            "actual_cached_9504_site_sums_blocked_for_all_prelaunch_inputs":True,
        },
        "static_path_argument":"At fixed Y0, F_Q is affine in the five material/drive inputs. The subsequent Y0-to-p(theta) segment at final theta has exact affine strain with C-hull and radius zero. Its H_Q−H_Q0 term is covered by the rowwise predictor remainder bound. Combine signs before applying common parameter increments.",
        "independent_signed_corner_fixtures":random_count,
        "strict_signed_force_budget_exact":str(r*(1-z)),
        "actual_cached_9504_site_sums_run":False,"root_asserted":False,
        "input_hashes":dict(sorted(inputs.items())),"review_source_file":own,"review_source_sha256":inputs[own],
    }
    if "--freeze" in sys.argv:
        receipt,gate = OUT/"replay.json",PREP/"prelaunch_independent_review.json"
        assert not receipt.exists() and not gate.exists()
        with receipt.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(report,indent=2)+"\n")
        report["support_file"] = receipt.relative_to(ROOT).as_posix()
        report["support_sha256"] = sha256(receipt.read_bytes()).hexdigest()
        with gate.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":report["status"],"checks":len(report["checks"]),
                      "strict_signed_force_budget":float(r*(1-z))}),flush=True)


if __name__ == "__main__":
    run()
