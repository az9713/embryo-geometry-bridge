"""Independent mathematical and input admission review of the frozen pilot launch."""
from contextlib import ExitStack
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
PRE = HERE.parent
PRIMARY = PRE.parents[2]
ROOT = PRIMARY.parents[1]
NATIVE = PRIMARY / "independent_actual_preparation/actual_replay/receipt.json"
MANIFEST = PRIMARY / "field_manifest.json"
BUDGET = PRE.parents[1] / "budgets_v1"
PRIMARY_PROTOCOL = PRIMARY / "preparation/protocol.json"
SRC = ROOT / "scripts/transfer_ambiguity_cached_pilots.py"
PINS = {
    SRC: "1b22ebc307c0cee51bf4193ff3999242b1cfe278d07c10cc63619f2795489c58",
    PRE / "protocol.json": "83c6e54dfcf775dbf1ca5b5eff64aacc37f5d3db074bba43bf6f029b59fe9fef",
    PRE / "preflight.json": "83b34a2f06c0579ad8c70da26accbbde7ef3248047b0446a118ff36635dfd8d6",
    BUDGET / "actual/results.json": "689948766db950ee4ddf399e6731aea08890cc7a1f84c5615d19e7ded68b0872",
    BUDGET / "independent_actual/receipt.json": "978a2a4095bb8735103a5ee2c11ad266a6b0e052e22af2b7a36ef29f1f38265c",
    MANIFEST: "308ca4c05e476ba6a14cbc7b8d34018983b8f11263e2c224488672b35eb38ac0",
    NATIVE: "c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd",
    PRIMARY_PROTOCOL: "9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94",
}


def digest(f):
    return sha256(f.read_bytes()).hexdigest()


def get(f):
    return json.loads(f.read_bytes())


def save(f, value):
    with f.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2)
        handle.write("\n")
    return digest(f)


def rejects(f):
    try:
        f()
    except ValueError:
        return True
    return False


def exact_fixture(m, core):
    import numpy as np
    # Distinct rows expose a column/row or maximum-entry versus row-sum error.
    raw = [[-2., 0.75], [0.125, -3.25]]
    R = [list(map(Q, row)) for row in raw]
    H = [[Q(-1, 2), Q(1, 8)], [Q(1, 8), Q(-1, 4)]]
    g = [Q(-1, 16), Q(1, 32)]
    expected_rf = [sum(a*b for a,b in zip(row,g,strict=True)) for row in R]
    expected_z0 = max(sum(abs(int(i==j)-sum(R[i][k]*H[k][j] for k in range(2))) for
                          j in range(2)) for i in range(2))
    expected_rnorm = max(sum(map(abs,row)) for row in R)
    with patch.object(m,"NFREE",2):
        packed = lambda rows: core.arb_mat(len(rows),len(rows[0]),[core.A(x) for a in rows for x in a])
        eta, z0, abs_r, info = m.preconditioned_point(core,packed([[x] for x in g]),packed(H),np.asarray(raw))
        assert [Q(x["upper_exact"]) for x in info["preconditioned_center_gradient_by_row"]] == expected_rf
        assert eta == max(map(abs,expected_rf)) and z0 == expected_z0
        assert Q(info["preconditioner_infinity_norm_upper_exact"]) == expected_rnorm
        r, egoal, reta, rz0 = Q(1,64), Q(1,16), Q(1,256), Q(1,8)
        rows = [Q(1,32), Q(1,16)]
        names = ("bulk_det_min","bulk_hoop_min","surface_radius_min",
                 "surface_tangent_squared_min","normalized_weight_min")
        domain={"term_counts":m.TERMS,"strict_lower_margins":{x:core.pack(core.A(1)) for x in names}}
        def toy_variation(cache,n,radius):
            assert cache == "toy" and n==2 and core.pack(radius)["upper_exact"] == str(r)
            return packed([[d] for d in rows]),core.A(1),domain
        with patch.object(core,"variation_rows",side_effect=toy_variation):
            output = m.test_radius(core,"toy",abs_r,reta,rz0,r,egoal)
            assert not m.test_radius(core,"toy",abs_r,reta,rz0,egoal,egoal)[
                "strict_selfmap_contraction_and_observer_radius"]
        v = max(sum(abs(a)*b for a,b in zip(row,rows,strict=True)) for row in R)
        z = rz0 + v
        assert z < 1 and reta+z*r < r
        assert output["status"] == "STRICT_CACHED_QUADRATURE_ROOT_CANDIDATE"
        assert Q(output["preconditioned_variation_upper_exact"]) == v
        assert Q(output["z_upper_exact"]) == z
        assert Q(output["eta_plus_zr_upper_exact"]) == reta + z*r
        assert Q(output["posterior_error_upper_exact"]) == reta/(1-z) < r < egoal
    return {"synthetic_binary_R_norm_exact":str(expected_rnorm),"synthetic_Z0_exact":str(expected_z0),
            "synthetic_variation_exact":str(v),"synthetic_root_displacement_upper_exact":str(reta/(1-z))}


def protocol_fixture(m, prototype, bindings):
    assert m.prepare_protocol()==prototype
    bad_gate = {"status":"PASS","checks":{"nonliteral":1},**bindings}
    for label, example in (("nonliteral",bad_gate),
                           ("mixed",{**bad_gate,"checks":{"x":True},"budget_actual_sha256":"a"*64}),
                           ("empty",{**bad_gate,"checks":{}}),
                           ("fail",{**bad_gate,"status":"FAIL","checks":{"x":True}})):
        gate=HERE/f"synthetic_{label}_rejection.json"
        dg=save(gate,example)
        assert rejects(lambda:m.prelaunch(gate,dg))
        assert rejects(lambda:m.prelaunch(gate,"d"*64))
    identity={"pid":32255,"supervisor_pid":54321,"pilot_id":"test","bindings":bindings}
    terminal={**identity,"status":"FINISHED"}
    m.admit_child(0,terminal,identity)
    for code in (False,1,0.0,None):
        assert rejects(lambda:m.admit_child(code,terminal,identity))
    for bad in ({**terminal,"pid":9},{**terminal,"supervisor_pid":2},
                {**terminal,"bindings":{}},{**terminal,"pilot_id":"elsewhere"},
                {**terminal,"status":"FAILED"}):
        assert rejects(lambda:m.admit_child(0,bad,identity))
    assert not m.contraction(Q(1,8),Q(1,2),Q(1,4),Q(1,2))[
        "strict_selfmap_contraction_and_observer_radius"]
    assert not m.contraction(Q(1,16),Q(1,2),Q(1,4),Q(1,4))[
        "strict_selfmap_contraction_and_observer_radius"]
    assert m.contraction(Q(1,16),Q(1),Q(1,4),Q(1,2))["posterior_error_upper_exact"] is None


def main():
    assert all(digest(f)==h for f,h in PINS.items())
    prep, pre, native, manifest, primary, reviewed, budget = map(get,(
        PRE/"protocol.json",PRE/"preflight.json",NATIVE,MANIFEST,PRIMARY_PROTOCOL,
        BUDGET/"independent_actual/receipt.json",BUDGET/"actual/results.json"))
    assert digest(PRE/"transfer_ambiguity_cached_pilots_frozen.py")==PINS[SRC]
    assert pre["source_sha256"]==prep["source_sha256"]==PINS[SRC]
    assert pre["protocol_sha256"]==PINS[PRE/"protocol.json"]
    assert pre["status"]=="PASS" and pre["checks"] and all(x is True for x in pre["checks"].values())
    for k in ("actual_preconditioners_recovered","actual_primary_field_payloads_opened",
              "actual_cached_sites_assembled","new_nonlinear_work"):
        assert pre[k]==0
    assert not (PRE.parent/"actual").exists()
    hashes=prep["input_hashes"] | prep["source_and_prerequisite_hashes"]
    assert all(digest(ROOT/f)==h for f,h in hashes.items())
    assert len(manifest["fields"])==native["checkpoints"]==44
    assert manifest["primary_results_sha256"]==native["actual_results_sha256"]
    assert manifest["independent_primary_review_sha256"]==PINS[NATIVE]
    assert manifest["independent_replay_observed_os_exit_code"]==0
    assert all(row[f"{label}_sha256"]==native["all_hashes"][row[label]]
               for row in manifest["fields"] for label in ("field","verification"))
    assert native["source_compatible_witness"] is native["exact_equilibrium_position_bounded"] is False
    assert native["heldout_response_tables_read"] is False and native["new_equilibrium_solves"] is False
    assert prep["precision_bits"]==256 and prep["torch_cpu_threads"]==1
    assert prep["expected_site_counts"]=={"bulk":144*8**2,"surface":24*12}
    radii=tuple(map(Q,prep["descending_dyadic_radii_exact"]))
    assert radii==tuple(Q(1,2**k) for k in (20,24,28,32,36))
    cases={row["endpoint_id"]:row for row in primary["endpoints"]}
    full={(row["endpoint_id"],row["checkpoint"]):row for row in manifest["fields"]}
    identities=[]
    for p,b,a in zip(prep["ordered_pilots"],budget["four_pilot_selection"],
                     reviewed["four_pilot_fields"],strict=True):
        keys=("endpoint_id","checkpoint","record_id","field","field_sha256","e_goal_dimensionless_exact")
        assert all(p[k]==a[k]==b[k] for k in keys)
        j=p["endpoint_id"],p["checkpoint"]
        identities.append(j)
        assert all(p[k]==full[j][k] for k in ("field","field_sha256","step","verification","verification_sha256"))
        assert p["pilot_id"]==f"{j[0]}_{j[1]}"
        assert len(p["free_map"])==len(set(p["free_map"]))==611
        assert len(set(p["active_elements"]))==24
        design=cases[j[0]]
        material=design["materials"]["small" if j[1]=="small_17" else "large"]
        assert material==p["material"]
        optical=Q(p["step"]["dose_exact"])
        gain=Q(material["gain_exact"])
        drive=gain*optical
        rounded=Q(float(drive))
        assert rounded!=drive>0
        assert Q(p["step"]["requested_strength_exact"])==drive
        assert Q(p["step"]["solver_strength_binary64_exact"])==rounded
        assert {k:Q(v) for k,v in p["exact_drive"].items()}=={
            "optical_dose_exact":optical,"common_gain_exact":gain,"activation_strength_exact":drive,
            "activation_strength_binary_exact":rounded,"binary_input_minus_exact_strength":rounded-drive}
        assert 0<radii[-1]<=radii[0]<Q(p["e_goal_dimensionless_exact"])/2
    assert identities==[(k,c) for k in cases for c in ("small_17","large_03")]
    assert prep["rounded_solver_input_used_as_exact_proof_drive"] is False
    assert prep["new_nonlinear_equilibrium_solves"]==0
    assert prep["archived_fixture_cache_or_preconditioner_reused"] is False

    binds=dict(zip(("source_sha256","protocol_sha256","preflight_sha256","budget_actual_sha256",
                    "budget_independent_actual_review_sha256","whole_44_field_manifest_sha256"),
                   (PINS[x] for x in (SRC,PRE/"protocol.json",PRE/"preflight.json",
                                     BUDGET/"actual/results.json",BUDGET/"independent_actual/receipt.json",MANIFEST)),strict=True))
    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_cached_pilots as m
    import algebraic_root_feasibility as core
    def forbid(*_a,**_k):
        raise AssertionError("actual work forbidden in prelaunch readback")
    original=Path.open
    def no_payload(path,*args,**kw):
        if path.resolve().is_relative_to(PRIMARY/"actual") and path.name in ("field.json","verification.json"):
            forbid()
        return original(path,*args,**kw)

    with ExitStack() as plain:
        for module, funcs in ((m,("model_input","sparse_proposal","calculate","worker","launch")),
                              (core,("local_terms","freeze_cache","replay","attempt"))):
            for name in funcs:
                plain.enter_context(patch.object(module,name,side_effect=forbid))
        plain.enter_context(patch.object(Path,"open",no_payload))
        plain.enter_context(patch.object(subprocess,"Popen",side_effect=forbid))
        protocol_fixture(m,prep,binds)
        ctx0=core.ctx.prec
        core.ctx.prec=256
        core.A.cache_clear()
        try:
            values=exact_fixture(m,core)
            for p in prep["ordered_pilots"]:
                altered=dict(p["step"],requested_strength_exact=p["exact_drive"]["activation_strength_binary_exact"])
                assert rejects(lambda:m.exact_drive(altered,p["material"]))
                assert m.exact_drive(p["step"],p["material"])==p["exact_drive"]
            calls=[]
            def local(model,u,strength):
                calls.append((model,u,strength,type(strength)))
                yield {"kind":"bulk"}
            toy=SimpleNamespace(local_terms=local,assemble_point=lambda terms,n:(0,0,list(terms),0))
            pending=Q(1,10)
            m.point_assembly(toy,"toy",[],pending,2,{"bulk":1})
            assert calls==[("toy",[],pending,Q)]
            assert rejects(lambda:m.point_assembly(toy,"toy",[],float(pending),2,{"bulk":1}))
            assert rejects(lambda:m.point_assembly(toy,"toy",[],pending,2,{"bulk":2}))
        finally:
            core.ctx.prec=ctx0
            core.A.cache_clear()

    assert not (PRE.parent/"actual").exists()
    checks={
      "all_frozen_bindings_and_native_44_metadata_verified":True,
      "four_authoritative_observer_budget_fields_and_exact_gain_times_optical_dose":True,
      "high_precision_bootstrap_and_fresh_binary_normalized_model_R_source_review":True,
      "exact_9504_site_force_Hessian_and_whole_cube_variation_source_review":True,
      "weighted_absolute_rows_and_matrix_norms_checked_with_nondiagonal_synthetic_R":True,
      "strict_radii_selfmap_posterior_and_closed_cube_math_checked":True,
      "failed_mixed_nonliteral_and_wrong_child_exit_rejected":True,
      "serial_os_exit_receipts_exclusive_output_and_failure_paths_source_reviewed":True,
      "no_actual_field_payload_force_hessian_cache_or_worker_opened":True,
    }
    support=HERE/"support.json"
    sh=save(support,{"status":"PASS","checks":checks,"source_and_input_hashes":hashes,
        "audit_source_sha256":digest(Path(__file__)),"small_matrix_fixture":values,
        "all_radii_below_half_each_budget":True,
        "premise_reuse":"The previously independently reviewed analytic local and whole-cube kernel is reused. This review exercises new wrapper calculations only on synthetic data, and verifies metadata and source.",
        "proof":"For fixed exact drive, local strain maps are affine in u on the closed cube. The positive cached determinant, hoop, tangent, radius and weights keep the asserted law analytic with penalty inactive at all 9504 sites. If t_i bounds sum_j |H(u)-H(c)|_ij, then the i-th preconditioned difference row sum is bounded by (|R|t)_i, by the triangle inequality. The independent max row defect bound Z0 gives z. T(u)=u-R F(u) has Lipschitz constant at most z, while eta+z*r<r maps the closed cube into its interior. At the center RH is invertible via the Neumann series for I-(I-RH), giving invertible square R. The unique local fixed point is a zero, with error at most eta/(1-z)<r<e_goal. No positive-definiteness implication is made.",
        "scope":"Exact real potential from rounded binary quadrature caches with exact rational optical drive. A pilot is separate from continuum, exact spatial integration, source compatibility and minimum/whole-body certificates.",
        "actual_point_R_H_or_cache_work":0})
    gate=PRE/"prelaunch_independent_review.json"
    gate_hash=save(gate,{"status":"PASS","checks":checks,**binds,
        "support_file":support.relative_to(ROOT).as_posix(),"support_sha256":sh})
    print(gate_hash,sh)


if __name__=="__main__":
    main()
