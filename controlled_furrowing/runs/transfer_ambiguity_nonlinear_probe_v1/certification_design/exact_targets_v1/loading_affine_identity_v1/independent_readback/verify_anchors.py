"""Metadata and exact scalar checks for the existing affine loading memo and endpoint anchors."""
import base64
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import re
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent
EX=HERE.parents[1]
ROOT=next(p for p in EX.parents if (p/"PLAN.md").is_file())
MEMO=HERE.parent/"READBACK.md"
MEMO_SHA="61e1a3"
links={}


def digest(b):
    return sha256(b).hexdigest()


def read(path,h=None,parse=True):
    p=(ROOT/path).resolve()
    name=p.relative_to(ROOT).as_posix()
    data=p.read_bytes()
    d=digest(data)
    if h is not None:
        assert d.startswith(h) if len(h)<64 else d==h, name
    assert links.setdefault(name,d)==d,name
    return json.loads(data) if parse else data


def members(cache):
    z=cache["frozen_binary_cache"]
    data=base64.b64decode(z["data"],validate=True)
    assert digest(data)==z["sha256"] and len(data)==z["bytes"]
    with ZipFile(BytesIO(data)) as f:
        # These streams are byte digests only. No matrix values are decoded or multiplied.
        return {n:digest(f.read(n)) for n in f.namelist() if n not in ("u.npy","preconditioner.npy")}


def main():
    source=read(MEMO,MEMO_SHA,False).decode()
    locators=re.findall(r"\]\(([^)]+)\)",source)
    for link in locators:
        path,_,span=link.partition("#")
        data=read(MEMO.parent/path,parse=False)
        if span.startswith("L"):
            assert 1<=int(span[1:])<=len(data.splitlines())
    prev=read("evidence/energy_gain_cell_independent_review.json")
    geom=read("evidence/energy_gain_geometry_independent_review.json")
    gain=read("evidence/energy_gain_cell.json")
    body=read("evidence/energy_gain_geometry.json")
    assert prev["status"]==geom["status"]=="PASS" and all(v is True for v in geom["checks"].values())
    assert Q(gain["gain_relative_halfwidth_exact"])==Q(1,1000)
    read("runs/exact_loading_tube_design_v1/METHOD.md",parse=False)
    old=read(EX/"preparation/target_manifest.json","cb25fbeaef753cea700c187cd8b377b7dcf4aa1fe2b0b246ed17134a07a89bc3")
    new=read(EX/"remaining32_targets_v1/preparation/manifest.json","dfea134640f99b45e9ca75d1ce718043d67f2917375506997d2465abda8c1242")
    cases=("max_negative_contrast_change_small_16","max_negative_contrast_change_small_17")
    k,= (r for r in new["targets"] if r["case_id"]==cases[0])
    l,= (r for r in old["targets"] if r["pilot_id"]==cases[1])
    ck,cl=(read(r["model_contract"],r["model_contract_sha256"]) for r in (k,l))
    pk,pl=(read(c["point_cache"],c["point_cache_sha256"]) for c in (ck,cl))
    for r in (k,l):
        read(r["rule"],r["rule_sha256"],False)
        read(r["routes"],r["routes_sha256"],False)
    assert k["rule_sha256"]==l["rule_sha256"] and ck["parameter_cell"]==cl["parameter_cell"]
    assert ck["material_parameters_binary_exact"]==cl["material_parameters_binary_exact"]==pk["material_parameter_storage_exact"]==pl["material_parameter_storage_exact"]
    mk,ml=members(pk),members(pl)
    assert mk==ml and len(mk)==21
    assert Q(pk["normalization_binary_exact"])==Q(pl["normalization_binary_exact"])==Q(ck["energy_scale_binary_exact"])==Q(cl["energy_scale_binary_exact"])>0
    doses,gains,drive=zip(*[tuple(Q(c["exact_requested_drive"][s]) for s in ("optical_dose_exact","common_gain_exact","activation_strength_exact")) for c in (ck,cl)],strict=True)
    assert doses==(Q(128898317,10**8),Q(2360271,1562500)) and doses[1]-doses[0]==Q(22159027,10**8)
    assert gains[0]==gains[1]>0 and drive==tuple(i*gains[0] for i in doses)
    assert all(Q(p["activation_strength_exact"])==a!=Q(p["activation_strength_binary_exact"]) for p,a in zip((pk,pl),drive,strict=True))
    base=EX/"remaining32_bridge_v1/first_wave_reduction_v1/independent_actual_preparation/actual_replay"
    kr=read(base/"receipt.json","e73ea05a4703256235a46ce7b914298d3f8258e9fef029f4dbe952b04bb9798a")
    ko=read(base/"os_exit_observation.json","aa3c90101c64ccefe1403af2c296bcb22cafad3b2e35d73539aa7e3c5227f281")
    kt=read(base/"terminal.json",ko["terminal_sha256"])
    assert type(ko["observed_exit_code"]) is int and ko["observed_exit_code"]==kt["exit_code"]==0
    assert kr["pid"]==ko["pid"]==kt["pid"] and kr["exact_spatial_roots"]==kr["strict_local_minima"]==4
    ka,= (r for r in kr["case_results"] if r["case_id"]==cases[0])
    kp=read(ka["native_file"],ka["native_sha256"])
    assert kp["status"]=="PASS_INDEPENDENT_FIRST_WAVE_EXACT_SPATIAL_ROOT" and ka["strict_local_minimum"] is True
    assert kp["target"]["model_contract_sha256"]==k["model_contract_sha256"]
    assert kp["target"]["point_cache_sha256"]==ck["point_cache_sha256"] and Q(kp["target"]["exact_drive"])==drive[0]
    lr=read(EX/"bridge_reduction_v1/independent_actual_preparation/all_four_readback.json",
            "859f507ad3762cfef44d53672bbe15a43592ce54d6e08ba8828bf50a925ec340")
    la,= (r for r in lr["cases"] if r["pilot_id"]==cases[1])
    lp=read(la["independent_proof_file"],la["independent_proof_sha256"])
    lo=read(la["independent_os_observation_file"],la["independent_os_observation_sha256"])
    assert type(lo["observed_exit_code"]) is int and lo["observed_exit_code"]==0 and lo["pid"]==lp["pid"]
    assert Q(la["whole_box_normalized_H_min_exact"])>0 and la["selfmap"]["strict_selfmap"] is True
    assert lp["status"]=="PASS_INDEPENDENT_EXACT_SPATIAL_ROOT_AND_STRICT_LOCAL_MINIMUM"
    assert lp["target_identity"]==la["target_identity"] and lp["target_identity"]["model_contract_sha256"]==l["model_contract_sha256"]
    assert lp["target_identity"]["cache_packet_sha256"]==cl["point_cache_sha256"] and Q(lp["target_identity"]["drive_exact"])==drive[1]
    out={"status":"PASS_INDEPENDENT_AFFINE_LOADING_SOURCE_AND_K_L_ANCHOR_READBACK","checked_links":len(locators),
         "unique_memo_link_paths":len({x.split("#")[0] for x in locators}),"source_sha256":digest(Path(__file__).read_bytes()),
         "bound_inputs_sha256":links,"verified_binary_operator_member_hashes_common_K_L":mk,
         "anchors":list(cases),"exact_optical_doses":list(map(str,doses)),"exact_dose_gap":str(doses[1]-doses[0]),
         "common_scale_exact":ck["energy_scale_binary_exact"],"same_normalization":True,
         "new_mechanical_or_matrix_jobs":0,"exact_spatial_loading_interval_certified":False}
    p=HERE/"receipt.json"
    with p.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(out,f,indent=2)
        f.write("\n")
    print(digest(p.read_bytes()))


if __name__=="__main__":
    main()
