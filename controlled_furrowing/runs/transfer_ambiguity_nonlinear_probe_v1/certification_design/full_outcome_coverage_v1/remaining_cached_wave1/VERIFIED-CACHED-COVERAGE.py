"""Verify and consolidate the two record M cached-quadrature results."""
import base64
from fractions import Fraction as Q
from hashlib import sha256
import io
import json
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
DESIGN=HERE.parents[1]
ROOT=DESIGN.parents[2]
MAIN=HERE/"independent_actual_preparation/actual_replay"
NARROW=HERE/"negative_M_narrow_v1"
FOUR=DESIGN/"cached_pilots_v1/independent_actual_preparation/actual_replay/receipt.json"
LEDGER=HERE.parent/"receipt.json"
PLUS="max_contrast_change_small_09"
MINUS="max_negative_contrast_change_small_09"
PINS={
    HERE/"preparation/protocol.json":"c1543cbbd2520664313b56114b66a885faaf4d75869f56dbb1aa83a43a78fe3a",
    HERE/"independent_actual_preparation/replay.py":"1895a08aa52b3d007e871e8e4d0f15d032fcb9abf209228fd4e37137b5283a15",
    HERE/"independent_actual_preparation/protocol.json":"b04ccb595239fc9715f8e67b819725ede554cf7d602188bc025f2c95ed63315b",
    HERE/"independent_actual_preparation/preflight.json":"49e4a5ec1ec7202db731191285d1c08beb37c507627c250603520ae1ff192331",
    MAIN/"receipt.json":"9854b5033a055ea33d7af818cf8da27b8eee35b3f4bb57606fcbf9bca96e3b56",
    MAIN/"terminal.json":"745374287982959495edb20217fbf8d2a6828acaceae9b87b51d4f7f7b917fde",
    MAIN/"os_exit_observation.json":"42e11a4b7065b07eca318e5127951e884e8bd7703773b37df184bc9b747706ef",
    MAIN/f"{PLUS}.json":"c5c78d779e4537859c6a0aaf4b206a1339ba6cc7d1392a99fdcb4d7255d007b4",
    MAIN/f"{MINUS}.json":"17f68e45c56df1a19bfa8233520e23a5ea6453128d3095780d017deb7b61214a",
    NARROW/"replay.py":"f8e310e95bc68beb691960f0d77de60a9ac31326b9ee15bfdad5f07dc8fd1e24",
    NARROW/"protocol.json":"4a967869eb203f6b6d0f7f5f231567e21c05cb5a9338d7c528ee81f534e9ae47",
    NARROW/"preflight.json":"0b29b819a6ace56b8cbfd4f45d5dcf01482b6302e146889a3bc77ded2dbf9071",
    NARROW/"root_source_readback.json":"460d4320277c26517e5d46d499ac59d46e7b5387d7d739403d19996baf7119b6",
    NARROW/"actual/receipt.json":"7aa4c765a6df459e802e2d1e99c7fb007ab533e3b29ef07cee2d92df3f52d4f1",
    NARROW/"actual/terminal.json":"5face2fc79564c12b842de61ffad2493686fb026bb93d05ed1ad36408ba9a93d",
    NARROW/"actual/os_exit_observation.json":"3fa3737870a585f58ae7ca846b9d696c2ee98df87820c44f503529214273a6b8",
    FOUR:"4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17",
    LEDGER:"235e8336c735b1bf1a75889018ae6a18586fd0905a65d38cb683c5ed236df19a",
}
COUNTS={"bulk":9216,"surface":288}
SCOPE=("Two additional unique local zeros of the prescribed exact-real binary cached-quadrature potentials. "
       "No exact spatial-integration, stability or whole-body geometry certificate, stable loading path, "
       "continuum or experimental validation.")


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def name(p):
    return p.relative_to(ROOT).as_posix()


def os_exit(directory,pid,exec_session):
    t=read(directory/"terminal.json")
    o=read(directory/"os_exit_observation.json")
    assert t["status"]=="FINISHED" and type(t["exit_code"]) is int and t["exit_code"]==0
    assert t["pid"]==o["pid"]==pid and o["unified_exec_session"]==exec_session
    assert type(o["observed_exit_code"]) is int and o["observed_exit_code"]==0
    assert t["receipt_sha256"]==o["receipt_sha256"]==PINS[directory/"receipt.json"]
    assert o["terminal_sha256"]==PINS[directory/"terminal.json"]
    return o


def proof_sum(row,base,c,cache,expected_radius):
    assert row["case_id"]==base["case_id"]==c["case_id"]==cache["case_id"]
    assert row["record_id"]==base["record_id"]==cache["record_id"]==c["record_id"]=="optoGEF2_small_col_M"
    assert row["field"]==base["field"]==c["field"]==cache["input_state"]
    assert row["field_sha256"]==base["field_sha256"]==c["field_sha256"]==cache["input_state_sha256"]
    assert c["verification_sha256"]==cache["field_verification_sha256"]
    assert cache["parameters"]==c["material"]["parameters"]
    assert cache["material_parameter_storage_exact"]==c["material"]["parameter_storage_exact"]
    assert cache["closed_parameter_cell_exact"]==c["closed_parameter_cell_exact"]
    assert all(a==b for a,b in c["closed_parameter_cell_exact"].values())
    assert Q(c["step"]["dose_exact"])==Q(78208913,100000000)
    exact=Q(c["material"]["gain_exact"])*Q(c["step"]["dose_exact"])
    assert Q(row["exact_strength"])==Q(base["exact_strength"])==Q(cache["activation_strength_exact"])==exact
    assert Q(float(exact))==Q(cache["activation_strength_binary_exact"])!=exact
    assert base["rounded_strength_used_in_proof"] is False
    blob=cache["frozen_binary_cache"]
    raw=base64.b64decode(blob["data"],validate=True)
    assert len(raw)==blob["bytes"] and sha256(raw).hexdigest()==blob["sha256"]==base["cache_array_contents_sha256"]
    with np.load(io.BytesIO(raw),allow_pickle=False) as archive:
        arrays={k:archive[k] for k in archive.files}
        assert {k:list(v.shape) for k,v in arrays.items()}==blob["array_shapes"]
        assert {k:sha256(v.tobytes(order="C")).hexdigest() for k,v in arrays.items()}==base["cache_tensor_content_sha256"]
        assert arrays["free"].tolist()==c["free_map"] and arrays["u"].shape==(611,)
        assert arrays["preconditioner"].shape==(611,611) and arrays["preconditioner"].dtype==np.dtype("<f8")
        R=max(sum(abs(Q(x)) for x in line) for line in arrays["preconditioner"].tolist())
    assert R==Q(base["preconditioner_infinity_norm_exact"])==Q(row["preconditioner_infinity_norm_exact"])
    assert base["site_counts"]==COUNTS
    eta=Q(base["independent_eta_upper_exact"])
    z0=Q(base["independent_Z0_upper_exact"])
    if row is not base:
        assert row["preconditioner_binary_sha256"]==base["cache_tensor_content_sha256"]["preconditioner"]
        assert eta==Q(row["eta_upper_exact"]) and z0==Q(row["center_Z0_upper_exact"])
        assert row["original_unresolved_test_preserved"] is True
        assert row["original_unresolved_radius_exact"]==base["radius_exact"]
        assert row["reused_signed_cached_anchor_sha256"]==PINS[MAIN/f"{MINUS}.json"]
        assert row["reused_signed_cached_anchor_file"]==name(MAIN/f"{MINUS}.json")
    signed=base["signed_center_force_intervals"]
    assert len(signed)==611 and all(-eta<=Q(p["lower_exact"])<=Q(p["upper_exact"])<=eta for p in signed)
    r=Q(row["radius_exact"])
    z=Q(row["independent_ZQ_upper_exact"])
    lhs=Q(row["strict_selfmap_lhs_upper_exact"])
    post=Q(row["independent_posterior_error_upper_exact" if row is base else "posterior_error_upper_exact"])
    goal=Q(row["observer_goal_exact"])
    assert r==expected_radius and goal==Q(c["e_goal_dimensionless_exact"])
    assert z==z0+Q(row["independent_variation_upper_exact"])
    assert 0<=eta and 0<=z<1 and lhs==eta+z*r<r<goal
    assert post==eta/(1-z)<r
    geo=row["analytic_domain" if row is base else "independent_whole_cube_domains"]
    assert geo["counts"]==COUNTS and all(Q(v["lower_exact"])>0 for v in geo["strict_lower_margins"].values())
    return {"case_id":c["case_id"],"record_id":c["record_id"],
        "field":c["field"],"field_sha256":c["field_sha256"],"cache_packet_sha256":row["cache_packet_sha256"],
        "preconditioner_binary_sha256":base["cache_tensor_content_sha256"]["preconditioner"],
        "exact_requested_drive":str(exact),"radius_exact":str(r),"cached_ZQ_upper_exact":str(z),
        "posterior_error_upper_exact":str(post),"observer_goal_exact":str(goal),
        "display":{"radius":float(r),"cached_ZQ":float(z),"posterior_dimensionless_displacement":float(post),
                   "observer_goal":float(goal),"fraction_of_observer_goal":float(post/goal)}}


def main():
    assert all(sha(p)==v for p,v in PINS.items())
    old=read(MAIN/"receipt.json")
    narrow=read(NARROW/"actual/receipt.json")
    main_exit=os_exit(MAIN,34020,29833)
    narrow_exit=os_exit(NARROW/"actual",8136,3263)
    assert narrow_exit["root_source_readback_sha256"]==PINS[NARROW/"root_source_readback.json"]
    assert old["status"]=="INDEPENDENT_TWO_M_CACHED_REPLAY_WITH_UNRESOLVED_CASE"
    assert old["cached_local_roots"]==1 and old["spatial_root_count"]==0
    assert narrow["status"]=="PASS_INDEPENDENT_NEGATIVE_M_SMALLER_CACHED_ROOT"
    assert old["producer_supervisor_pid"]==1996 and old["producer_supervisor_exec_session"]==22696
    assert type(old["producer_observed_exit_code"]) is int and old["producer_observed_exit_code"]==0
    plus,minus=[read(MAIN/f"{p}.json") for p in (PLUS,MINUS)]
    assert plus["status"]=="PASS_INDEPENDENT_M_CACHED_LOCAL_ROOT"
    assert minus["status"]=="UNRESOLVED_INDEPENDENT_M_CACHED_ROOT"
    assert Q(minus["radius_exact"])==Q(1,2**32) and Q(minus["independent_ZQ_upper_exact"])>1
    assert minus["independent_posterior_error_upper_exact"] is None
    assert old["results"][0]["native_sha256"]==PINS[MAIN/f"{PLUS}.json"]
    assert old["results"][1]["native_sha256"]==PINS[MAIN/f"{MINUS}.json"]
    fixed=read(HERE/"preparation/protocol.json")
    cold=read(HERE/"independent_actual_preparation/protocol.json")
    cold2=read(NARROW/"protocol.json")
    hashes=old["producer_result_bytes_and_payloads"]|cold["input_hashes"]|cold2["input_sha256"]
    for file,h in hashes.items():
        assert sha(ROOT/file)==h,file
    cplus,cminus=fixed["cases"]
    roots=[]
    for b,c,cert,path,r in (
        (plus,cplus,plus,MAIN/f"{PLUS}.json",Q(1,2**36)),
        (minus,cminus,narrow,NARROW/"actual/receipt.json",Q(1,2**34)),
    ):
        cachefile=HERE/"actual"/c["case_id"]/"cache_packet.json"
        assert hashes[name(cachefile)]==b["cache_packet_sha256"]==cert["cache_packet_sha256"]
        roots.append(proof_sum(cert,b,c,read(cachefile),r)|{"native_proof_file":name(path),"native_proof_sha256":PINS[path]})
    four=read(FOUR)
    assert four["status"]=="PASS_INDEPENDENT_FOUR_CACHED_LOCAL_ROOTS" and four["root_certificates_cached_quadrature_only"]==4
    ledger=read(LEDGER)
    assert ledger["status"]=="PASS_READ_ONLY_36_FIELD_COVERAGE_INVENTORY"
    assert all(sha(p)==v for p,v in PINS.items())
    output={"status":"PASS_VERIFIED_TWO_ADDITIONAL_M_CACHED_LOCAL_ROOTS",
        "scope":SCOPE,"new_cached_local_root_count":2,"prior_other_cached_local_root_count":4,
        "exact_spatial_root_count":0,"ad_or_equilibrium_jobs_in_readback":0,
        "case_results":roots,"original_negative_M_unresolved":{"native_file":name(MAIN/f"{MINUS}.json"),
            "native_sha256":PINS[MAIN/f"{MINUS}.json"],"radius_exact":minus["radius_exact"],
            "derivative_bound_upper_exact":minus["independent_ZQ_upper_exact"]},
        "true_exit_observations":{"original_independent_pair":main_exit,"separate_negative_M_smaller_cube":narrow_exit},
        "checks":{"all_authoritative_packet_producer_and_native_sources_rehashed":True,
                  "same_exact_drive_cached_binary_arrays_preconditioner_fields_and_observer_budgets":True,
                  "strict_rational_selfmap_and_posterior_for_each_certified_radius":True,
                  "unresolved_larger_minus_M_and_prior_four_cached_rows_preserved":True},
        "source_and_input_sha256":hashes|{name(p):h for p,h in PINS.items()},
        "review_source_sha256":sha(Path(__file__))}
    path=HERE/"VERIFIED-CACHED-COVERAGE.json"
    with path.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(output,f,indent=2,allow_nan=False)
        f.write("\n")
    print(sha(path),len(output["source_and_input_sha256"]),[r["display"] for r in roots])


if __name__=="__main__":
    main()
