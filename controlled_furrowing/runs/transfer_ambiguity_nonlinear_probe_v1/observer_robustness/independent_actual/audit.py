"""Independent postrun readback of the four-observer stored-field sweep."""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import importlib.util
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
SWEEP=HERE.parent
PRIMARY=SWEEP.parent
ROOT=HERE.parents[3]
A=SWEEP/"actual"
P=SWEEP/"preparation"
PK=PRIMARY/"independent_actual_preparation"
SOURCE=ROOT/"scripts/transfer_ambiguity_observer_sweep.py"
PINS={
 A/"results.json":"53a17296ddae08b4d1fe6b7e9b20ad0a377bf8dc962142bbfb9eb9877a381647",
 A/"terminal.json":"65e6a2158094aba9c565232814e3a5ade303fde68476c37b63bca0765624177e",
 SOURCE:"9eacb1de680cefb84d90d603cc1aca4dd6ef4892372a91cffa3480651c5874a8",
 P/"protocol.json":"dc48cee2cdf5f0da3b49f6972dc37191cf482a238f37b3185f8ceaf1f93bbeea",
 P/"preflight.json":"1e9e318c7f6ae06d99ddb1a12eb64c78f101812dc16631b56cc582e474f5247b",
 P/"prelaunch_independent_review.json":"144ca9889393133d1cb7174a7c11b88e89a36a55ce497c7d0763336e2af9407c",
 PRIMARY/"field_manifest.json":"308ca4c05e476ba6a14cbc7b8d34018983b8f11263e2c224488672b35eb38ac0",
 PK/"actual_replay/receipt.json":"c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd",
 PK/"replay_frozen.py":"e1ff63fddd6981f9ef8f05183c884537546cb66d59b9ec063911b22a1c50d27d",
 PRIMARY/"actual/results.json":"2ddc2fb3f0da426551bf74b6e2aac7bb489dfe848f8d9e73cfc42fddecbe99d6",
 PRIMARY/"actual/terminal.json":"84aafea79a60878e139ed034c199e4875e91bac21fa650354f12f944e87fd0b9",
 PRIMARY/"preparation/protocol.json":"9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94",
}
VARIANTS=("current_disk__spatial_max","current_disk__regional_mean","current_radial_line__spatial_max","current_radial_line__regional_mean")
OBS=("maximum_shortening_um","region_averaged_shortening_um","apical_max_shortening_um","basal_max_shortening_um")
ENDPOINTS=("max_contrast_change","max_negative_contrast_change")


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_bytes())


def binding(p,expected,hashes):
    rel=p.relative_to(ROOT).as_posix()
    assert sha(p)==expected==hashes.setdefault(rel,expected),rel


def yes(x):
    return bool(x) and all(v is True for v in x.values())


def as_q(v):
    return Q(v["numerator"],v["denominator"])


def pair(x):
    a,b=map(Q,x)
    assert a<=b
    return a,b


def scalar_error(y,v,d):
    a,b=v
    lo=max(a-y,y-b,Q(0))
    hi=max(abs(y-a),abs(y-b))
    return lo,hi,"IN" if hi<=d else ("OUT" if lo>d else "UNRESOLVED")


def maxima(rows):
    return tuple(max(x[i] for x in rows) for i in (0,1))


def verdict(bounds,delta):
    return ("ALLOWANCE_SATISFIED" if bounds[1]<=delta else
            ("ALLOWANCE_EXCEEDED" if bounds[0]>delta else "UNRESOLVED"))


def keypair(x):
    return x["endpoint_id"],x["checkpoint"]


def serialized_bounds(data):
    assert data["lower_closed"] is data["upper_closed"] is True
    return as_q(data["lower"]),as_q(data["upper"])


def line_mean(edges,radius,graph):
    """Independently integrate Z(t) R'(t) on the exact radius-bracketed segments."""
    total=(Q(0),Q(0))
    for e in edges:
        if e[0][0]>=radius:
            break
        rad,z=graph.power(e,0),graph.power(e,1)
        assert min(rad[1],rad[1]+2*rad[2])>0
        cut=graph.bracket(rad,min(radius,e[-1][0]))
        coefficients=graph.conv(z,[rad[1],2*rad[2]])
        primitive=[Q(0)]+[value/(i+1) for i,value in enumerate(coefficients)]
        total=graph.add(total,graph.interval_horner(primitive,cut))
    return graph.scale(total,1/radius)


def replay_graph(state,packet,graph,geometry):
    L=Q(67,2)
    controls=geometry.geometric_controls(state,L)
    physical=graph.calc_small(state,controls,geometry)
    H=L*Q(state["parameters"]["height"])
    radius=L*Q(state["parameters"]["patch_radius"])
    contract=packet["observer_contract"]
    assert tuple(map(as_q,contract["current_radius_ROI_um"]))==(0,radius)
    assert as_q(contract["reference_apical_um"])==H and as_q(contract["reference_basal_um"])==0
    assert contract["weighting"]=="current_disk" and contract["component_aggregation"]=="spatial_max"
    assert contract["root_bracketing_steps"]==80
    a,b=controls["apical"],controls["basal"]
    disk_a=graph.disk_and_extrema(a,radius)[0]
    disk_b=graph.disk_and_extrema(b,radius)[0]
    line_a,line_b=line_mean(a,radius,graph),line_mean(b,radius,graph)
    assert graph.add((H,H),graph.minus(disk_b,disk_a))==physical[OBS[1]]
    cert=packet["maximum_shortening_certificate"]
    mx=serialized_bounds(cert["interval"])
    assert graph.overlap(mx,physical[OBS[0]])
    assert cert["requested_width_achieved"] is True and mx[1]-mx[0]<=Q(1,1000)
    assert Q(packet["positive_thickness_lower_um_exact"])>0
    variants={}
    for label,ma,mb in (("current_disk",disk_a,disk_b),("current_radial_line",line_a,line_b)):
        common={OBS[0]:physical[OBS[0]],OBS[1]:graph.add((H,H),graph.minus(mb,ma))}
        variants[label+"__spatial_max"]={**common,OBS[2]:physical[OBS[2]],OBS[3]:physical[OBS[3]]}
        variants[label+"__regional_mean"]={**common,OBS[2]:graph.minus((H,H),ma),OBS[3]:mb}
    return variants,mx


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    hashes={}
    for f,h in PINS.items():
        binding(f,h,hashes)
    terminal,start,result=[read(A/f) for f in ("terminal.json","started.json","results.json")]
    binding(A/"started.json",terminal["start_sha256"],hashes)
    assert terminal["status"]=="FINISHED" and terminal["pid"]==start["pid"]==33196
    assert terminal["results_sha256"]==PINS[A/"results.json"] and terminal["bindings"]==start["bindings"]==result["bindings"]
    assert start["utc"]<terminal["utc"] and result["status"]=="COMPLETE_SECONDARY_STATIC_OBSERVER_SWEEP"
    expected={"source_sha256":SOURCE,"protocol_sha256":P/"protocol.json","preflight_sha256":P/"preflight.json",
              "primary_protocol_sha256":PRIMARY/"preparation/protocol.json",
              "secondary_prelaunch_review_sha256":P/"prelaunch_independent_review.json",
              "primary_actual_review_sha256":PK/"actual_replay/receipt.json",
              "whole_44_field_manifest_sha256":PRIMARY/"field_manifest.json",
              "primary_results_sha256":PRIMARY/"actual/results.json",
              "primary_terminal_sha256":PRIMARY/"actual/terminal.json"}
    assert result["bindings"]=={k:PINS[path] for k,path in expected.items()}
    gate,pre,prot=[read(P/f) for f in ("prelaunch_independent_review.json","preflight.json","protocol.json")]
    assert gate["status"]==pre["status"]=="PASS" and yes(gate["checks"]) and yes(pre["checks"])
    assert prot["source_sha256"]==PINS[SOURCE] and prot["variant_order"]==list(VARIANTS)
    for name,h in (prot["source_hashes"]|prot["input_hashes"]).items():
        binding(ROOT/name,h,hashes)
    binding(P/"transfer_ambiguity_observer_sweep_frozen.py",PINS[SOURCE],hashes)
    binding(P/"independent_prelaunch/audit.py",gate["review_source_sha256"],hashes)
    primary,native,manifest,t1,primary_protocol=[read(x) for x in (
        PRIMARY/"actual/results.json",PK/"actual_replay/receipt.json",PRIMARY/"field_manifest.json",
        PRIMARY/"actual/terminal.json",PRIMARY/"preparation/protocol.json")]
    assert t1["status"]=="FINISHED" and t1["results_sha256"]==PINS[PRIMARY/"actual/results.json"]
    assert all(x["observed_exit_code"]==0 for x in t1["workers"]) and len(t1["workers"])==2
    nt=read(PK/"actual_replay/terminal.json")
    assert nt["status"]=="FINISHED" and nt["exit_code"]==0 and nt["receipt_sha256"]==PINS[PK/"actual_replay/receipt.json"]
    binding(PK/"actual_replay/terminal.json",sha(PK/"actual_replay/terminal.json"),hashes)
    assert native["checkpoints"]==44 and native["independently_recomputed_small_residuals"]==136
    assert native["actual_results_sha256"]==manifest["primary_results_sha256"]==PINS[PRIMARY/"actual/results.json"]
    assert manifest["primary_terminal_sha256"]==PINS[PRIMARY/"actual/terminal.json"]
    assert manifest["independent_primary_review_sha256"]==PINS[PK/"actual_replay/receipt.json"]
    assert manifest["independent_replay_observed_os_exit_code"]==0
    binding(PRIMARY/"freeze_field_manifest.py",manifest["source_sha256"],hashes)
    for k,h in native["all_hashes"].items():
        binding(ROOT/k,h,hashes)
    rows=manifest["fields"]
    actual=[dict(endpoint_id=case["endpoint_id"],**row) for case in primary["endpoints"] for row in case["checkpoints"]]
    assert len(rows)==len(actual)==44 and len({keypair(x) for x in rows})==44
    for row,old,plan in zip(rows,actual,prot["frozen_checkpoint_plan"],strict=True):
        assert keypair(row)==keypair(old)==keypair(plan)
        for k in ("role","step","field","verification"):
            assert row[k]==old[k]==plan[k]
        for file in ("field","verification"):
            assert row[file+"_sha256"]==old[file+"_sha256"]
            binding(ROOT/row[file],row[file+"_sha256"],hashes)
    assert result["primary_replay_exactly_matches"] is True
    for case in primary["endpoints"]:
        assert result["by_endpoint_and_variant"][case["endpoint_id"]][VARIANTS[0]]==case["small_residual"]
    assert result["retuning"] is False and result["additional_equilibria"]==result["additional_force_evaluations"]==0
    assert prot["heldout_large_response_read"] is False
    delta=Q(prot["conditional_source_allowance_um_exact"])
    assert delta==Q(result["frozen_primary_source_allowance_um_exact"])==Q(primary_protocol["selection"]["conditional_source_allowance_um_exact"])
    small={keypair(x):x for x in rows if x["role"]=="small_record"}
    checkpoints={}
    matches=0
    counts=Counter()
    graph=module(PK/"replay_frozen.py","primary_frozen")
    rr=ROOT/"runs/transfer_pair_pilot_v1/independent_actual/replay.py"
    assert native["all_hashes"][rr.relative_to(ROOT).as_posix()]==sha(rr)
    sys.path.insert(0,str(ROOT/"scripts"))
    geometry=module(rr,"transfer_graph_reviewer")
    for doc in result["field_observer_receipts"]:
        fp=ROOT/doc["path"]
        binding(fp,doc["sha256"],hashes)
        t=read(fp)
        k=keypair(t)
        row=small[k]
        assert fp==A/k[0]/(k[1]+".json") and k not in checkpoints
        assert all(row[name]==t[name] for name in ("role","step","field","verification","field_sha256","verification_sha256"))
        state=read(ROOT/row["field"])
        case=next(x for x in primary_protocol["endpoints"] if x["endpoint_id"]==k[0])
        assert state["parameters"]==case["materials"]["small"]["parameters"] and state["converged"] is state["reached_target"] is True
        assert Q(state["activation_strength"])==Q(t["step"]["solver_strength_binary64_exact"])
        assert t["record_id"]==row["step"]["record_id"]
        independent,shared_max=replay_graph(state,t,graph,geometry)
        archive={v:{name:pair(x) for name,x in vals.items()} for v,vals in t["variant_observables_um_exact"].items()}
        assert set(archive)==set(independent)==set(VARIANTS)
        for v,name in product(VARIANTS,OBS):
            assert graph.overlap(archive[v][name],independent[v][name])
            matches+=(archive[v][name]==independent[v][name])
            counts[name]+=(archive[v][name]==independent[v][name])
        assert all(archive[v][OBS[0]]==shared_max for v in VARIANTS)
        for v in VARIANTS:
            # Fixed-ROI means use independent antiderivative and must agree.
            assert archive[v][OBS[1]]==independent[v][OBS[1]]
            if v.endswith("regional_mean"):
                assert all(archive[v][name]==independent[v][name] for name in OBS[2:])
        v1=read(ROOT/row["verification"])
        assert {k:pair(x) for k,x in v1["observer_intervals_um_exact"].items()}==archive[VARIANTS[0]]
        checkpoints[k]=(row["step"],archive,independent)
    assert len(checkpoints)==len(small)==result["measured_small_fields"]==34
    totals={}
    recomposed={}
    for v in VARIANTS:
        branch=[]
        for ep in ENDPOINTS:
            rec=[checkpoints[k] for k in checkpoints if k[0]==ep]
            measurement=result["by_endpoint_and_variant"][ep][v]
            table=measurement["all68_residuals"]
            assert Q(measurement["declared_allowance_um_exact"])==delta
            assert len(rec)==table["records"]==17 and len(table["constraint_rows"])==table["constraints"]==68
            actual_rows,independent_rows=[],[]
            for (r,k),t in zip(product(rec,OBS),table["constraint_rows"],strict=True):
                step,ar,ind=r
                y=Q(step["observed_exact"][k])
                assert t["record_id"]==step["record_id"] and t["observable"]==k and Q(t["observed_exact"])==y
                assert Q(t["discrepancy_allowance_exact"])==delta
                a,b,tag=scalar_error(y,ar[v][k],delta)
                assert pair(t["prediction_interval_exact"])==ar[v][k]
                assert a==Q(t["minimum_interval_residual_exact"]) and b==Q(t["maximum_interval_residual_exact"])
                assert tag==t["status"]
                actual_rows.append((a,b,k,step["record_id"]))
                da,db,_=scalar_error(y,ind[v][k],delta)
                independent_rows.append((da,db,k,step["record_id"]))
            for k in OBS:
                ar=maxima([a for a in actual_rows if a[2]==k])
                assert ar==(Q(table["discrepancy_lower_requirements_exact"][k]),Q(table["discrepancy_sufficient_requirements_exact"][k]))
            ar,ir=maxima(actual_rows),maxima(independent_rows)
            assert graph.overlap(ar,ir)
            assert ar==pair(measurement["required_allowance_um_exact"])
            assert verdict(ar,delta)==verdict(ir,delta)==measurement["status"]
            assert {"ALLOWANCE_SATISFIED":"IN","ALLOWANCE_EXCEEDED":"OUT","UNRESOLVED":"UNRESOLVED"}[measurement["status"]]==table["status"]
            worst=max(actual_rows,key=lambda r:r[1])
            branch.append((ar,ir,(ep,worst[2],worst[3])))
        ar,ir=maxima([x[0] for x in branch]),maxima([x[1] for x in branch])
        closed=result["common_pair_requirement_by_variant"][v]
        st=verdict(ar,delta)
        assert ar==pair(closed["required_allowance_um_exact"]) and st==verdict(ir,delta)==closed["status_at_frozen_allowance"]
        totals[v]=ar
        recomposed[v]={"required_allowance_um_exact":list(map(str,ar)),"independent_required_allowance_um_exact":list(map(str,ir)),
                       "required_allowance_um_display":list(map(float,ar)),"status_at_original_delta":st,
                       "limiting_endpoint_observable_record":max(branch,key=lambda x:x[0][1])[2]}
    all4=result["common_pair_requirement_across_all_four_conventions"]
    worst=maxima(list(totals.values()))
    assert worst==pair(all4["required_allowance_um_exact"]) and verdict(worst,delta)==all4["status_at_frozen_allowance"]=="ALLOWANCE_EXCEEDED"
    assert recomposed[VARIANTS[0]]["status_at_original_delta"]=="ALLOWANCE_SATISFIED"
    assert all(recomposed[v]["status_at_original_delta"]=="ALLOWANCE_EXCEEDED" for v in VARIANTS[1:])
    assert result["small_constraint_evaluations"]==544
    assert all(sha(ROOT/k)==h for k,h in hashes.items())
    receipt={
      "status":"PASS_INDEPENDENT_SECONDARY_STORED_OBSERVER_SWEEP",
      "scope":"Archived static stored fields, fixed current-radius ROI and flat reference. Changing observer weight or component aggregation is a contract sensitivity check, without force/solve or biological inference.",
      "checks":{"terminal_and_all_current_sources_bound":True,"native_44_field_manifest_and_primary_receipts_unchanged":True,
                "independent_P2_controls_disk_and_line_means_and_component_extrema":True,
                "all_544_exact_interval_residuals_reconstructed":True,"all_four_common_thresholds_independently_recomputed":True},
      "primary_exact_result_kept":True,"source_configured_delta_um_exact":str(delta),
      "original_delta_um_display":float(delta),"scope_of_common_allowance":"One identical scalar discrepancy bound within each stated observer variant.",
      "by_observer_common_requirements":recomposed,"all_four_required_common_allowance_um_exact":list(map(str,worst)),
      "stored_field_constraints":544,"independent_interval_exact_matches":matches,"exact_matches_by_observable":dict(counts),
      "geometry_premises":"Independent boundary topology and exact physical controls from reviewed prior replay. Disk and line antiderivatives and component extrema independently rebuilt and paired with the reused exact same-radius total-force-free graph maximum and thickness kernels.",
      "OS_scope":"Parent directly observed secondary process exit zero. Static reviewer matched PID 33196 and pinned terminal/start/results; no independent operating-system observation is claimed.",
      "independent_force_or_equilibrium_jobs_run":0,"source_compatibility_established":False,"robust_across_conventions_at_frozen_delta":False,
      "heldout_large_responses_opened":False,"input_and_output_sha256":hashes,"audit_source_sha256":sha(Path(__file__))
    }
    path=HERE/"receipt.json"
    with path.open("x",encoding="utf-8") as stream:
        json.dump(receipt,stream,indent=2)
        stream.write("\n")
    print(sha(path),len(hashes),matches)
    for v,row in recomposed.items():
        print(v,row["status_at_original_delta"],row["required_allowance_um_display"],row["limiting_endpoint_observable_record"])


if __name__=="__main__":
    main()

