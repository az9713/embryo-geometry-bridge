"""Verify the public context-open path on each native proof class, without mechanics."""
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
sys.path.insert(0, str(ROOT/"scripts"))
import transfer_ambiguity_remaining32_exact_targets as a

pin = "dfea134640f99b45e9ca75d1ce718043d67f2917375506997d2465abda8c1242"
plan = a.read(a.PREP/"manifest.json",pin)
assert plan["geometry_or_variation_actual_launch_authorized"] is False
assert len(plan["targets"]) == 32
assert a.reject(lambda:a.open_context(a.IDS[0],"0"*64))
assert a.reject(lambda:a.open_context("max_contrast_change_small_17",pin))
cases = ["max_contrast_change_small_01", "max_contrast_change_small_08",
         "max_contrast_change_small_09", "max_negative_contrast_change_small_09",
         "max_contrast_change_small_13", "max_negative_contrast_change_small_16"]
sample = []
for case in cases:
    ctx,rule,d = a.open_context(case,pin)
    bound = a.read(a.canonical(d["model_contract"]),d["model_contract_sha256"])
    assert ctx.record_id == "transfer_exact_"+case and len(rule["volume_weights"])==144 and len(ctx.u)==611
    assert a.a_sha(ctx.u)==bound["center_binary_sha256"] and a.a_sha(ctx.R)==bound["binary_R_sha256"]
    assert a.node_box(ctx,a.Q(d["planned_geometry_radius_exact"]))["current_centers_exact"]==a.doc(ctx.target.nodes())
    assert (a.Q(bound["archived_shift_exact"]) > 0) == (case in a.RECENTERED)
    assert bound["root_claimed"] is False
    sample.append({"case_id":case,"model_contract_sha256":d["model_contract_sha256"]})
doc={"status":"PASS_PUBLIC_CONTEXT_READBACK","manifest_sha256":pin,"cases":sample,
     "wrong_external_digest_and_foreign_case_rejected":True,"new_geometry_AD_integrations_or_solves":0,
     "source_sha256":sha256(Path(__file__).read_bytes()).hexdigest()}
print(a.X.save(HERE/"context_readback.json",doc))
