"""Independent replay of four exact-drive binary-cache root candidates."""
import argparse
import base64
from contextlib import ExitStack
from fractions import Fraction as Q
from hashlib import sha256
import importlib.metadata
import io
import json
import os
from pathlib import Path
import re
import sys
import time
import traceback
from unittest.mock import patch

SOURCE = Path(__file__).resolve()
HERE = SOURCE.parent
PILOTS = HERE.parent
PRIMARY = PILOTS.parents[1]
ROOT = PRIMARY.parents[1]
OLD = PILOTS / "preparation"
ACTUAL = PILOTS / "actual"
REPLAY = HERE / "actual_replay"
N = 611
PREC = 256
PARTS = ("child_started", "worker_started", "observed_os_exit", "terminal", "results")
PRODUCER_RESULTS_SHA = "3138dd8c6c72dd086a9864d6ac3e662d8e1a70ecbd17d401c688687e0beb0757"
AD_SOURCE = ROOT / "scripts/review_algebraic_root.py"
AD_PRIOR_REVIEW = ROOT / "evidence/algebraic_root_independent_review.json"
PINS = {
    ROOT / "scripts/transfer_ambiguity_cached_pilots.py": "1b22ebc307c0cee51bf4193ff3999242b1cfe278d07c10cc63619f2795489c58",
    OLD / "protocol.json": "83c6e54dfcf775dbf1ca5b5eff64aacc37f5d3db074bba43bf6f029b59fe9fef",
    OLD / "preflight.json": "83b34a2f06c0579ad8c70da26accbbde7ef3248047b0446a118ff36635dfd8d6",
    OLD / "prelaunch_independent_review.json": "84c13d4a3289ac9b6d4a03cf0760549f9439c732b3eb0880e973ec34cf4b693c",
    AD_SOURCE: "6af2311d0667a28af263b6fb6237063481aa8e65209ae9722c2b0246a8cff27a",
    AD_PRIOR_REVIEW: "f327f2fab76f9fdcab3a83882adbb1a4ae3f0cf1e5c269e7335fb909691085d5",
}
COUNTS = {"bulk": 9216, "surface": 288}
RADIUS = Q(1, 2**36)
SCOPE = ("Unique local zeros only of the four exact-real finite binary-cache quadrature potentials, "
         "at their exact requested fixed rational strengths. No exact spatial integration, energy "
         "minimum, full-body geometry or PDE or source-data inference.")
MARGINS = {"bulk_det_min", "bulk_hoop_min", "bulk_J_min", "surface_radius_min",
           "surface_tangent_squared_min", "surface_J_min", "normalized_weight_min"}


def name(f):
    return f.relative_to(ROOT).as_posix()


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def pinned(path, expected):
    assert isinstance(expected,str) and re.fullmatch("[0-9a-f]{64}", expected), path
    data = path.read_bytes()
    assert sha256(data).hexdigest() == expected, path
    return json.loads(data)


def save(path, doc):
    with path.open("x", encoding="utf-8", newline="\n") as out:
        json.dump(doc, out, indent=2, allow_nan=False)
        out.write("\n")
    return sha(path)


def verified(f, h):
    assert sha(f) == h, f


def source_protocol():
    for f, h in PINS.items():
        verified(f, h)
    old = pinned(OLD / "protocol.json", PINS[OLD / "protocol.json"])
    prior = pinned(AD_PRIOR_REVIEW, PINS[AD_PRIOR_REVIEW])
    assert prior["status"] == "PASS" and prior["reviewer_source_sha256"] == PINS[AD_SOURCE]
    assert prior["checks"] and all(x is True for x in prior["checks"].values())
    assert old["expected_site_counts"] == COUNTS
    refs = old["input_hashes"] | old["source_and_prerequisite_hashes"]
    for f, h in refs.items():
        verified(ROOT/f, h)
    ids = [p["pilot_id"] for p in old["ordered_pilots"]]
    assert len(ids) == len(set(ids)) == 4 and all(re.fullmatch("[a-z0-9_]+", x) for x in ids)
    names = [name(ACTUAL/f"{kind}.json") for kind in ("started", "terminal", "results")]
    names += [name(ACTUAL / ident / f"{kind}.json") for ident in ids for kind in PARTS]
    assert len(names) == len(set(names)) == 23
    return {"source_sha256":sha(SOURCE), "immutable_pins":{name(f):h for f,h in PINS.items()} | refs,
            "producer_top_results_sha256":PRODUCER_RESULTS_SHA,"required_supervisor_pid":27828,
            "required_supervisor_session":46941,"required_supervisor_exit_code":0,
            "external_parent_pins":names,"ordered_pilot_ids":ids,"precision_bits":PREC,
            "accepted_radius_to_replay_exact":str(RADIUS), "local_site_counts":COUNTS,
            "runtime_packages":old["runtime_packages"],
            "proof_target":"Exact binary X+injection(u), frozen binary quadrature and weights; exact requested Fraction gain*dose.",
            "independence":"Pinned prior independent forward interval AD and assembler, regenerated current binary model arrays, independent exact norms and Banach inequalities. No producer cache assembler/variation calls. Existing cache R is the shared map by design.",
            "scope":SCOPE}


def load_ad():
    sys.path.insert(0,str(ROOT/"scripts"))
    import review_algebraic_root as ad
    assert Path(ad.__file__).resolve()==AD_SOURCE and sha(AD_SOURCE)==PINS[AD_SOURCE]
    ad.ctx.prec=PREC
    ad.A.cache_clear()
    return ad


def admission_map(admission, plan):
    assert admission["supervisor_pid"] == plan["required_supervisor_pid"]
    assert admission["unified_exec_session"] == plan["required_supervisor_session"]
    code = admission["observed_exit_code"]
    assert type(code) is int and code == plan["required_supervisor_exit_code"]
    h = admission["sha256"]
    assert set(h)==set(plan["external_parent_pins"])
    assert h[name(ACTUAL/"results.json")]==plan["producer_top_results_sha256"]
    assert all(isinstance(value,str) and re.fullmatch("[0-9a-f]{64}",value) for value in h.values())
    return h


def admit(path, digest, plan, old):
    admission = pinned(path, digest)
    h = admission_map(admission, plan)
    docs = {key:pinned(ROOT/key,value) for key,value in h.items()}
    def pick(directory, kind):
        return docs[name(directory/f"{kind}.json")]
    start, terminal, whole = (pick(ACTUAL,kind) for kind in ("started","terminal","results"))
    assert start["pid"] == terminal["pid"] == plan["required_supervisor_pid"]
    assert terminal["status"] == "FINISHED" and start["utc"] < terminal["utc"]
    assert terminal["results_sha256"] == plan["producer_top_results_sha256"]
    assert whole["status"] == "COMPLETE_FOUR_CACHED_QUADRATURE_PILOT_ATTEMPTS"
    assert whole["independent_actual_review_required"] is True
    ids=plan["ordered_pilot_ids"]
    assert start["ordered_pilot_ids"] == ids
    assert [x["pilot_id"] for x in terminal["workers"]] == [x["pilot_id"] for x in whole["pilots"]] == ids
    assert len(terminal["workers"])==len({x["pid"] for x in terminal["workers"]})==4
    gate = pinned(OLD/"prelaunch_independent_review.json",PINS[OLD/"prelaunch_independent_review.json"])
    bindings = {k:gate[k] for k in ("source_sha256","protocol_sha256","preflight_sha256","budget_actual_sha256",
        "budget_independent_actual_review_sha256","whole_44_field_manifest_sha256")}
    bindings["prelaunch_review_sha256"] = PINS[OLD/"prelaunch_independent_review.json"]
    assert whole["bindings"] == terminal["bindings"] == start["bindings"] == bindings
    output = []
    prior_utc = start["utc"]
    for ident, roster, result_meta, prototype in zip(ids,terminal["workers"],whole["pilots"],
                                                     old["ordered_pilots"],strict=True):
        f = ACTUAL/ident
        parent_start, worker_start, worker_exit, worker_terminal, result = (pick(f,kind) for kind in PARTS)
        assert roster == worker_exit | {"os_exit_receipt_sha256":h[name(f/"observed_os_exit.json")]}
        pid = worker_exit["pid"]
        assert type(pid) is int and pid > 0 and pid != start["pid"]
        common = {"pid":pid,"supervisor_pid":start["pid"],"pilot_id":ident,"bindings":bindings}
        for observation in (parent_start,worker_start,worker_exit,worker_terminal):
            assert all(observation[k]==v for k,v in common.items())
        assert type(worker_exit["observed_exit_code"]) is int and worker_exit["observed_exit_code"] == 0
        assert worker_terminal["status"]=="FINISHED"
        assert prior_utc <= parent_start["utc"] <= worker_exit["utc"] <= terminal["utc"]
        assert prior_utc <= worker_start["utc"] <= worker_terminal["utc"] <= worker_exit["utc"]
        prior_utc = worker_exit["utc"]
        assert worker_exit["terminal_sha256"]==h[name(f/"terminal.json")]
        assert worker_terminal["results_sha256"]==result_meta["results_sha256"]==h[name(f/"results.json")]
        assert result["pilot"] == prototype and result["bindings"] == bindings
        assert result["status"] == result_meta["status"] == "CANDIDATE_CERTIFIED_CACHED_QUADRATURE_LOCAL_ROOT"
        assert result["cache_packet"] == result_meta["cache_packet"] == name(f/"cache_packet.json")
        assert result["cache_packet_sha256"] == result_meta["cache_packet_sha256"]
        assert result["new_nonlinear_equilibrium_solves"] == 0
        assert result["independent_actual_result_review_required"] is True
        assert result["exact_spatial_integration_established"] is result["whole_body_geometry_established"] is result["local_minimum_established"] is False
        output.append((f,result))
    return h, output


def cache_extract(meta, packet, pilot, f, hashes):
    import numpy as np
    from transfer_pair_pilot import AxisymmetricModel, IntegrationPolicy, Parameters, SolvePolicy, SparsePolicy
    from dataclasses import asdict

    integration = meta["numerical_policy"]
    expected_material = pilot["material"]["parameters"]
    assert asdict(Parameters(**expected_material)) == packet["parameters"] == expected_material
    storage = pilot["material"]["parameter_storage_exact"]
    assert {k:Q(v) for k,v in expected_material.items()} == {k:Q(v) for k,v in storage.items()}
    assert packet["material_parameter_storage_exact"] == storage
    assert packet["precision_bits"] == 256 and packet["free_dof"] == N
    assert packet["integration_policy"] == integration["solve_quadrature"]
    assert packet["pilot_id"] == pilot["pilot_id"] and packet["record_id"] == pilot["record_id"]
    assert packet["input_state"] == pilot["field"]
    assert packet["input_state_sha256"] == pilot["field_sha256"]
    assert packet["field_verification_sha256"] == pilot["verification_sha256"]
    hashes[pilot["field"]] = pilot["field_sha256"]
    hashes[pilot["verification"]] = pilot["verification_sha256"]
    state = pinned(ROOT/pilot["field"],pilot["field_sha256"])
    verified(ROOT/pilot["verification"],pilot["verification_sha256"])
    assert state["reached_target"] is state["converged"] is True
    assert state["integration_policy"] == integration["solve_quadrature"]
    assert state["parameters"] == expected_material

    dose = Q(pilot["step"]["dose_exact"])
    gain = Q(pilot["material"]["gain_exact"])
    strength = dose * gain
    binary_strength = Q(state["activation_strength"])
    assert strength == Q(pilot["step"]["requested_strength_exact"]) > 0 and binary_strength != strength
    assert binary_strength == Q(float(strength)) == Q(pilot["step"]["solver_strength_binary64_exact"])
    drive = {"optical_dose_exact":dose,"common_gain_exact":gain,"activation_strength_exact":strength,
             "activation_strength_binary_exact":binary_strength,"binary_input_minus_exact_strength":binary_strength-strength}
    assert drive == {k:Q(v) for k,v in pilot["exact_drive"].items()} == {k:Q(packet[k]) for k in drive}
    assert Q(pilot["step"]["binary64_rounding_difference_exact"])==binary_strength-strength

    blob = packet["frozen_binary_cache"]
    raw = base64.b64decode(blob["data"],validate=True)
    assert len(raw)==blob["bytes"] and sha256(raw).hexdigest()==blob["sha256"]
    with np.load(io.BytesIO(raw),allow_pickle=False) as archive:
        arrays = {k:archive[k] for k in archive.files}
        assert len(arrays)==len(archive.files)
    assert blob["array_shapes"] == {k:list(v.shape) for k,v in arrays.items()}
    assert all(a.dtype.kind in ("i","u","f") and np.isfinite(a).all() for a in arrays.values())
    model = AxisymmetricModel(Parameters(**expected_material),SolvePolicy(**integration["solve"]),
                             IntegrationPolicy(**integration["solve_quadrature"]),SparsePolicy(**integration["sparse"]))
    live = {"reference":model.reference,"free":model.free,"triangles":model.triangles,
            "shape":model.shape.numpy(),"shape_grad":model.shape_grad.numpy(),"qr":model.qr.numpy(),
            "volume_weights":model.integration.numpy(),"active":model.active.numpy(),
            "energy_scale":np.array([model.energy_scale]),"u":np.asarray(state["free_displacement"])}
    for side, (edges,S,dS,r0,l0,w,K,T) in enumerate(model.surfaces):
        k = "basal" if side==0 else "apical"
        assert K==expected_material[f"{k}_area_modulus"] and T==expected_material[f"{k}_surface_tension"]
        for attr, value in zip(("edges","N","dN","r0","l0","weights"),(edges,S,dS,r0,l0,w),strict=True):
            live[f"surface_{side}_{attr}"]=value.numpy()
    assert set(arrays)==set(live)|{"preconditioner"}
    assert all(a.dtype==arrays[k].dtype and a.shape==arrays[k].shape
               and a.tobytes()==arrays[k].tobytes() for k,a in live.items())
    assert live["reference"].shape==(325,2) and live["triangles"].shape==(144,6)
    assert arrays["u"].shape==arrays["free"].shape==(N,)
    assert arrays["free"].tolist()==pilot["free_map"]
    assert arrays["reference"].tolist()==state["reference_nodes"] and arrays["triangles"].tolist()==state["triangles"]
    assert np.flatnonzero(model.activation_elements).tolist()==pilot["active_elements"]
    assert model.activation_elements.tolist()==state["activation_elements"]
    nf=set(pilot["free_map"])
    ref=[Q(x) for x in arrays["reference"].ravel().tolist()]
    stored=[Q(x) for row in state["current_nodes"] for x in row]
    u=dict(zip(pilot["free_map"],map(Q,arrays["u"].tolist()),strict=True))
    assert len(nf)==len(u)==N and len(ref)==len(stored)==650
    gaps=[abs(y-x-u.get(i,0)) for i,(x,y) in enumerate(zip(ref,stored,strict=True))]
    assert max(gaps)==Q(pilot["xi_dimensionless_exact"]) and all(g==0 for i,g in enumerate(gaps) if i not in nf)
    assert Q(float(model.energy_scale))==Q(packet["normalization_binary_exact"])>0
    assert arrays["preconditioner"].dtype==np.dtype("<f8") and arrays["preconditioner"].shape==(N,N)

    hfile=f/"fresh_sparse_hessian.npz"
    assert packet["hessian_archive"]==name(hfile) and packet["numeric_proposal_only"] is True
    verified(hfile,packet["hessian_archive_sha256"])
    hashes[name(hfile)]=packet["hessian_archive_sha256"]
    with np.load(hfile,allow_pickle=False) as h:
        from scipy.sparse import csr_matrix
        assert set(h.files)=={"H_data","H_indices","H_indptr","H_shape","free_indices","free_vector","gradient"}
        assert np.array_equal(h["free_indices"],arrays["free"]) and np.array_equal(h["free_vector"],arrays["u"])
        assert tuple(h["H_shape"])==(N,N)
        assert h["gradient"].shape==(N,) and np.isfinite(h["gradient"]).all() and np.isfinite(h["H_data"]).all()
        proposal=csr_matrix((h["H_data"],h["H_indices"],h["H_indptr"]),shape=(N,N))
        assert proposal.nnz==packet["matrix_nonzero_count"]
    return arrays, strength, proposal


def qhi(ball):
    return Q(str(ball.upper().fmpq()))


def norm_rows(ad, a, n):
    return max(qhi(sum((a[i,j].abs_upper() for j in range(a.ncols())),ad.arb(0))) for i in range(n))


def replay_bound(ad, g, H, cache, raw_R, radius, counts):
    n = len(raw_R)
    rr = [[Q(x) for x in row] for row in raw_R.tolist()]
    norm_r = max(sum(map(abs,row)) for row in rr)
    mat = ad.arb_mat(n,n,[ad.A(x) for row in rr for x in row])
    abs_r=ad.arb_mat(n,n,[ad.A(abs(x)) for row in rr for x in row])
    Id = ad.arb_mat(n,n)
    for i in range(n):
        Id[i,i]=1
    rg=mat*g
    eta = norm_rows(ad,rg,n)
    z0 = norm_rows(ad,Id-mat*H,n)
    rows, domain = ad.variation(cache,n,ad.A(radius))
    assert domain["counts"] == counts
    dz = norm_rows(ad,abs_r*rows,n)
    return rg, eta, z0, dz, norm_r, domain


def producer_assertions(folder, result, p, hashes):
    attempts = result["radius_attempts"]
    radii = tuple(Q(1,2**k) for k in (20,24,28,32,36))
    assert len(attempts)==len(radii)
    for k,(row,r) in enumerate(zip(attempts,radii,strict=True)):
        f = folder/f"radius_{k:02d}.json"
        assert row["path"]==f.name
        hashes[name(f)] = row["sha256"]
        raw=pinned(f,row["sha256"])
        assert raw=={a:b for a,b in row.items() if a not in ("path","sha256")}
        assert Q(row["radius_requested_exact"])==r
        assert row["strict_selfmap_contraction_and_observer_radius"] is (r==RADIUS)
        assert Q(p["e_goal_dimensionless_exact"])>r
        if row["status"]=="UNRESOLVED_QUADRATURE_ANALYTIC_DOMAIN":
            assert r!=RADIUS
            continue
        z,eta=Q(row["z_upper_exact"]),Q(row["eta_upper_exact"])
        assert z==Q(result["center_Z0_upper_exact"])+Q(row["preconditioned_variation_upper_exact"])
        assert eta==Q(result["eta_upper_exact"])
        assert Q(row["eta_plus_zr_upper_exact"])==eta+z*r
        assert Q(row["strict_selfmap_gap_lower_exact"])==r-eta-z*r
        if z<1:
            assert Q(row["posterior_error_upper_exact"])==eta/(1-z)
        else:
            assert row["posterior_error_upper_exact"] is None
        assert (z<1 and eta+z*r<r) is (r==RADIUS)
    assert result["accepted_radius"] == {a:b for a,b in attempts[-1].items() if a not in ("path","sha256")}
    return attempts[-1]


def replay_one(ad, meta, folder, result, inputs):
    import numpy as np
    start=time.perf_counter()
    p=result["pilot"]
    c=folder/"cache_packet.json"
    centerfile=folder/"point_enclosure.json"
    inputs[name(c)]=result["cache_packet_sha256"]
    inputs[name(centerfile)]=result["point_enclosure_sha256"]
    packet=pinned(c,inputs[name(c)])
    center=pinned(centerfile,inputs[name(centerfile)])
    assert packet["bindings"]==result["bindings"]
    assert packet["frozen_binary_cache"]["sha256"] == result["embedded_binary_cache_sha256"]
    assert center["expected_and_rebuilt_site_counts"]==COUNTS
    assert Q(center["eta_upper_exact"])==Q(result["eta_upper_exact"])
    assert Q(center["center_Z0_upper_exact"])==Q(result["center_Z0_upper_exact"])
    baseline=producer_assertions(folder,result,p,inputs)
    arrays, strength, proposal=cache_extract(meta,packet,p,folder,inputs)
    ad.A.cache_clear()
    g,H,cache,counts = ad.assemble(ad.terms_from_arrays(arrays,packet["parameters"],strength),N)
    assert counts==COUNTS and len(cache)==sum(COUNTS.values())
    rg,eta,z0,dz,rnorm,domain=replay_bound(ad,g,H,cache,arrays["preconditioner"],RADIUS,COUNTS)
    assert set(domain["strict_lower_margins"])==MARGINS
    assert all(Q(v["lower_exact"])>0 for v in domain["strict_lower_margins"].values())
    assert Q(domain["strict_lower_margins"]["bulk_J_min"]["lower_exact"])>Q(1,1_000_000)
    z=z0+dz
    goal=Q(p["e_goal_dimensionless_exact"])
    assert 0<=eta and 0<=z<1
    lhs=eta+z*RADIUS
    posterior=eta/(1-z)
    assert lhs<RADIUS and 0<=posterior<RADIUS<goal
    assert rnorm<=Q(center["preconditioner_infinity_norm_upper_exact"])
    own_force=[ad.pack(rg[i,0]) for i in range(N)]
    compared=center["preconditioned_center_gradient_by_row"]
    overlaps=0
    for left,right in zip(own_force,compared,strict=True):
        overlaps+=int(Q(left["upper_exact"])>=Q(right["lower_exact"])
                      and Q(right["upper_exact"])>=Q(left["lower_exact"]))
    assert overlaps==N
    sample=proposal.tocoo()
    diff=max(abs(float(H[a,b].mid())-float(x)) for a,b,x in zip(sample.row,sample.col,sample.data,strict=True))
    scale=max(1.,float(np.max(np.abs(sample.data))))
    return {
        "status":"PASS_INDEPENDENT_CACHED_LOCAL_ROOT","pilot_id":p["pilot_id"],"record_id":p["record_id"],
        "field":p["field"],"field_sha256":p["field_sha256"],"cache_packet_sha256":inputs[name(c)],
        "cache_array_contents_sha256":packet["frozen_binary_cache"]["sha256"],"site_counts":counts,
        "exact_strength":str(strength),"rounded_strength_used_in_proof":False,
        "npz_current_model_equal_all_except_R":True,"preconditioner_infinity_norm_exact":str(rnorm),
        "radius_exact":str(RADIUS),"science_radius_budget_exact":str(goal),
        "independent_eta_upper_exact":str(eta),"independent_Z0_upper_exact":str(z0),
        "independent_variation_upper_exact":str(dz),"independent_z_upper_exact":str(z),
        "independent_selfmap_lhs_upper_exact":str(lhs),"strict_selfmap_gap_lower_exact":str(RADIUS-lhs),
        "independent_posterior_error_upper_exact":str(posterior),
        "display":{"eta":float(eta),"Z0":float(z0),"z":float(z),"lhs":float(lhs),
                   "error_dimensionless":float(posterior),"error_micrometers":float(Q(67,2)*posterior),
                   "observer_budget_fraction":float(posterior/goal)},
        "signed_center_force_intervals":own_force,"analytic_domain":domain,
        "producer_signed_center_force_overlap_rows":overlaps,
        "producer_variation_bound_upper_exact":baseline["preconditioned_variation_upper_exact"],
        "producer_four_larger_unresolved_attempts_retained":True,
        "floating_proposal_on_stored_nonzero_entries_relative_difference_diagnostic":diff/scale,
        "seconds":time.perf_counter()-start,"scope":SCOPE,
    }


def fixture(ad, plan):
    import numpy as np
    exact=Q(1,10)
    q=list(map(ad.A,(1,Q(1,8),0,1,1)))
    g,H,_,_=ad.bulk_AD(q,ad.A(1),ad.A(5),ad.A(exact))
    bg,bH,_,_=ad.bulk_AD(q,ad.A(1),ad.A(5),ad.A(float(exact)))
    assert ad.A(Q(1,8)*(exact-Q(float(exact)))) in g[1]-bg[1]
    assert ad.A(exact-Q(float(exact))) in H[1][1]-bH[1][1]
    term={"kind":"bulk","q":q,"parameters":(ad.A(1),ad.A(5),ad.A(exact)),
          "B":[[ad.A(1 if j==1 else 0)] for j in range(5)],"indices":[0],"weight":ad.A(1)}
    g1,H1,cache,d=ad.assemble([term],1)
    _,eta,z0,dz,rnorm,domain=replay_bound(ad,g1,H1,cache,np.asarray([[0.5]]),Q(1,256),d)
    assert z0==Q(9,20) or abs(z0-Q(9,20))<Q(1,2**245)
    assert eta>=Q(11,160) and eta-Q(11,160)<Q(1,2**245)
    assert dz>=0 and rnorm==Q(1,2) and domain["counts"]=={"bulk":1,"surface":0}
    map_test=dict.fromkeys(plan["external_parent_pins"],"a"*64)
    map_test[name(ACTUAL/"results.json")]=PRODUCER_RESULTS_SHA
    a={"supervisor_pid":27828,"unified_exec_session":46941,"observed_exit_code":0,"sha256":map_test}
    assert admission_map(a,plan)==map_test
    for faulty in (a|{"observed_exit_code":False},a|{"supervisor_pid":12},
                   a|{"sha256":{}},a|{"sha256": map_test|{name(ACTUAL/"results.json"):"b"*64}}):
        try:
            admission_map(faulty,plan)
        except AssertionError:
            continue
        raise AssertionError("admitted a mixed or failed producer")
    return {"nondyadic_drive_changes_active_AD_gradient_and_Hessian":True,
            "synthetic_assembly_and_weighted_rows_at_256_bits":True,
            "requires_exact_23_hash_pins_and_true_supervisor_exit":True,
            "actual_producer_or_primary_field_payloads_opened":0,
            "actual_site_count":0}


def prepare():
    plan=source_protocol()
    assert not REPLAY.exists()
    old=Path.open
    def exclude(path,*a,**k):
        resolved=path.resolve()
        if resolved.is_relative_to(ACTUAL) or (resolved.is_relative_to(PRIMARY/"actual")
                and resolved.name in ("field.json","verification.json")):
            raise AssertionError("preparation payload access excluded")
        return old(path,*a,**k)
    with patch.object(Path,"open",exclude):
        checks=fixture(load_ad(),plan)
    assert not REPLAY.exists()
    protocol=HERE/"protocol.json"
    ph=save(protocol,plan)
    f=HERE/"preflight.json"
    h=save(f,{"status":"PASS","source_sha256":sha(SOURCE),"protocol_sha256":ph,
              "synthetic":checks,"scope":"Terminal-admission and independent interval-AD fixture only.",
              "force_or_root_replays_launched":0})
    print(ph,h,sha(SOURCE),flush=True)


def run(args):
    ph,fh=args.protocol_sha256,args.preflight_sha256
    plan=pinned(HERE/"protocol.json",ph)
    pre=pinned(HERE/"preflight.json",fh)
    assert pre["status"]=="PASS" and pre["source_sha256"]==sha(SOURCE) and pre["protocol_sha256"]==ph
    assert pre["force_or_root_replays_launched"]==0
    assert plan==source_protocol()
    old=pinned(OLD/"protocol.json",PINS[OLD/"protocol.json"])
    inputs,cases=admit(args.parent_admission.resolve(),args.parent_admission_sha256,plan,old)
    REPLAY.mkdir(exist_ok=False)
    stamp={"pid":os.getpid(),"protocol_sha256":ph,"preflight_sha256":fh,"source_sha256":sha(SOURCE),
           "parent_admission":str(args.parent_admission.resolve()),"parent_admission_sha256":args.parent_admission_sha256}
    save(REPLAY/"started.json",stamp|{"status":"RUNNING"})
    done={"status":"FAILED","exit_code":1,**stamp}
    try:
        ad=load_ad()
        import torch
        torch.set_num_threads(1)
        assert torch.get_num_threads()==1
        for package, version in plan["runtime_packages"].items():
            assert importlib.metadata.version(package)==version
        rows=[]
        for f,result in cases:
            x=replay_one(ad,old,f,result,inputs)
            ch=save(REPLAY/f"{x['pilot_id']}.json",x)
            rows.append({k:x[k] for k in ("pilot_id","status","record_id","display","field","field_sha256")}
                        | {"receipt_sha256":ch})
            print(json.dumps(rows[-1]),flush=True)
        assert len(rows)==4
        for f,h in inputs.items():
            verified(ROOT/f,h)
        assert plan==source_protocol()
        h=save(REPLAY/"receipt.json",{"status":"PASS_INDEPENDENT_FOUR_CACHED_LOCAL_ROOTS",
               "scope":SCOPE,"checks":{"all_four_independently_strict":True,"all_four_within_observer_error_budgets":True,
                        "all_four_exact_nondyadic_drives":True,"true_parent_and_child_exit_admission":True},
               "parents":stamp,"producer_and_payload_sha256":inputs,"field_results":rows,
               "precision_bits":PREC,"separate_forward_interval_AD_source_sha256":PINS[AD_SOURCE],
               "root_certificates_cached_quadrature_only":4,"root_certificates_exact_spatial_integral":0,
               "new_equilibrium_solves":0,"stability_or_whole_body_embedding_certified":False,
               "source_compatibility_or_PDE_certified":False})
        done.update(status="FINISHED",exit_code=0,receipt_sha256=h)
    except BaseException as error:
        save(REPLAY/"failure.json",{"message":str(error),"traceback":traceback.format_exc()})
        raise
    finally:
        save(REPLAY/"terminal.json",done)
        print(json.dumps(done),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    commands=parser.add_subparsers(dest="action",required=True)
    commands.add_parser("prepare")
    r=commands.add_parser("replay")
    r.add_argument("--parent-admission",required=True,type=Path)
    r.add_argument("--parent-admission-sha256",required=True)
    r.add_argument("--protocol-sha256",required=True)
    r.add_argument("--preflight-sha256",required=True)
    args=parser.parse_args()
    if args.action=="prepare":
        prepare()
    else:
        run(args)
