"""Deterministic cold wave profiles. Importing this module cannot launch work."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
BRIDGE=HERE.parent
ROOT=BRIDGE.parents[4]
LEDGER=BRIDGE.parent.parent/"full_outcome_coverage_v1/receipt.json"
RAW=BRIDGE/"actual_first_wave"
FIRST=("max_contrast_change_small_01","max_contrast_change_small_13",
       "max_negative_contrast_change_small_01","max_negative_contrast_change_small_16")
PINS={
 BRIDGE/"preparation/manifest.json":"d42b7260283ba02ccc34b9e9ba69ce76d58507115c22615c90d28fa38020cb41",
 LEDGER:"235e8336c735b1bf1a75889018ae6a18586fd0905a65d38cb683c5ed236df19a",
 RAW/"root_terminal_admission.json":"bb9f059699e6d98fe626491f252edddb1664fbf83d79cb74258084ffd399c880",
 ROOT/"scripts/transfer_ambiguity_remaining32_exact_bridge.py":"4d52e8fb8dbc01ba58c38052bbff537443693234499a9a3b9b3e138c82d9b51b",
 ROOT/"scripts/transfer_ambiguity_remaining32_spatial_worker.py":"61b0082a68a297f3ae287232704a9b72ed4fbf53e96090db990ba0d8a8458f2d",
 ROOT/"scripts/transfer_ambiguity_remaining32_geometry_variation.py":"68bca869ee486d7dc25b34cae3ceee1ddc2adc53500a20301c0038c259a2b26a",
 BRIDGE/"independent_component_replay_preparation/replay.py":"283044f9d213bdfa82d1333ae9be1e86b8d2dc94f5ceb7d8186e9c6e1a1c0720",
 BRIDGE/"independent_component_replay_preparation/hull_math.py":"f8b4e1e8b4fb9069d235303f3fd0a45ad629ccd9d0bf37bf30ceca509c96b38b",
 ROOT/"scripts/transfer_ambiguity_remaining32_first_wave_reduce.py":"58de9f5e96b43005b0ba83dd0f4d8cd1b35f0690eef58290a75ef5abce692d92",
 BRIDGE/"first_wave_reduction_v1/independent_actual_preparation/replay_global.py":"4682165bf59bd633092c042a7df619c9bc1ad272dc80d1755404c2999a65287f",
}
PARTS={"bulk":144,"surface":24}
R34=Fraction(1,2**34)


def rel(path):
    return path.resolve().relative_to(ROOT).as_posix()


def pinned(path,digest):
    b=path.read_bytes()
    if sha256(b).hexdigest()!=digest:
        raise ValueError(f"immutable input changed: {rel(path)}")
    return b


def get(path,d=None):
    return json.loads(pinned(path,d if d is not None else PINS[path]))


def require(value,msg):
    if not value:
        raise ValueError(msg)


def bound_targets(cases):
    out={}
    for e in cases:
        case=e["case_id"]
        require(type(case) is str and case not in out and e["record_id"]=="transfer_exact_"+case
                and Fraction(e["trial"]["trial_radius_exact"])==R34,"same ordered immutable target")
        cat=e["component_catalogue"]
        scopes=[(k,i) for k,v in PARTS.items() for i in range(v)]
        require(len(cat)==168 and all(type(row["slot"]) is int for row in cat)
                and [(row["kind"],row["slot"]) for row in cat]==scopes,"full exact 168 local identities")
        t=e["trial"]
        out[case]={"case_id":case,"record_id":e["record_id"],
            "source_record_id":e["source_record_id"],"cached_input_kind":e["cached_input_kind"],
            "original_observer_field_sha256":e["original_observer_field_sha256"],
            "original_addition_gap_exact":e["original_addition_gap_exact"],
            "observer_goal_exact":e["observer_goal_exact"],
            **{k:t[k] for k in ("trial_box_sha256","model_contract_sha256","point_cache_sha256",
                "R_binary_sha256","center_binary_sha256","archived_center_shift_exact","exact_drive",
                "energy_scale_binary_exact","trial_radius_exact","routes_sha256","rule_sha256")}}
    require(len(out)==32 and set(FIRST).issubset(out),"full all32 bridge catalogue")
    return out


def waves(targets,priority):
    """Freeze paired checkpoint grouping for complete coverage, regardless of outcomes."""
    require(len(priority)==len(set(priority))==len(targets)==32 and set(priority)==set(targets),"full original priority set")
    left=set(targets)-set(FIRST)
    suffixes=list(dict.fromkeys(case.rsplit("_",1)[1] for case in priority))
    paired=[]
    singles=[]
    for number in suffixes:
        members=[case for case in priority if case in left and case.rsplit("_",1)[1]==number]
        if len(members)==2:
            paired.append(members)
        elif len(members)==1:
            singles.extend(members)
    require(len(paired)==13 and len(singles)==2,"13 complete pairs and two remaining counterparts")
    units=paired+[singles]
    out=[]
    for i in range(0,len(units),2):
        here=set(units[i]+units[i+1])
        order=[k for k in targets if k in here]
        out.append({"wave_id":f"wave_{i//2+2:02}","case_order":order,
                    "targets":{case:targets[case] for case in order},
                    "part_jobs":[[case,part] for case in order for part in PARTS],
                    "producer_components":672,"independent_components":672,
                    "properly_observed_raw_artifact_count":715,"max_owned_native_workers":4,
                    "selector_and_raw_actual_authorized":False})
    validate(out,targets)
    return out


def validate(profiles,targets):
    require(len(profiles)==7 and [p["wave_id"] for p in profiles]==[f"wave_{k:02}" for k in range(2,9)],
            "seven distinct later phase IDs")
    listed=[]
    for p in profiles:
        case=p["case_order"]
        require(type(case) is list and len(case)==len(set(case))==4 and all(k in targets and k not in FIRST for k in case)
                and case==[k for k in targets if k in case],"exact 4-case canonical later scope")
        require(p["targets"]=={k:targets[k] for k in case} and p["part_jobs"]==[[k,part] for k in case for part in PARTS]
                and type(p["max_owned_native_workers"]) is int and p["max_owned_native_workers"]==4
                and p["selector_and_raw_actual_authorized"] is False,"same case metadata and eight parts")
        listed.extend(case)
    require(len(listed)==len(set(listed))==28 and set(listed)==set(targets)-set(FIRST),"all 28 exactly once")


def load_profile(path,digest,wave):
    """Data-only input for a separately reviewed future caller."""
    m=get(path,digest)
    require(m["status"]=="FROZEN_COLD_REMAINING28_EXECUTION_PROFILES"
            and m["new_actual_jobs_authorized"] is False and m["source_sha256"]==sha256(Path(__file__).read_bytes()).hexdigest(),
            "explicit cold profile manifest digest")
    targets=bound_targets(get(BRIDGE/"preparation/manifest.json")["cases"])
    validate(m["waves"],targets)
    chosen=[p for p in m["waves"] if p["wave_id"]==wave]
    require(len(chosen)==1,"known explicit later-wave profile ID")
    return chosen[0]


def timing_metadata():
    obs=get(RAW/"root_terminal_admission.json")
    require(obs["status"]=="PASS_ROOT_OBSERVED_REMAINING32_FIRST_WAVE_SPATIAL_TERMINALS"
            and type(obs["observed_exit_code"]) is int and obs["observed_exit_code"]==0
            and type(obs["pid"]) is int and obs["pid"]>0 and len(obs["sha256"])==715,"same completed first raw wave")
    pins=obs["sha256"]
    def duration(folder):
        s,t=folder/"started.json",folder/"terminal.json"
        for p in (s,t):
            pinned(p,pins[rel(p)])
        # File creation through last write is a wall estimate, not CPU time.
        return Fraction(t.stat().st_mtime_ns-s.stat().st_mtime_ns,10**9)
    elapsed=duration(RAW)
    parts=[{"case_id":c,"part":p,"wall_seconds_filesystem_exact":str(duration(RAW/"cases"/c/p))} for c in FIRST for p in PARTS]
    raw_sizes={}
    for case in FIRST:
        total=0
        for kind,n in PARTS.items():
            for i in range(n):
                p=RAW/"cases"/case/kind/f"{kind}_{i:03d}.json"
                require(rel(p) in pins,"same expected raw canonical filename")
                total+=p.stat().st_size
        raw_sizes[case]=total
    times=[]
    for c in FIRST:
        folder=BRIDGE/"plans"/c
        root_ob=folder/"root_os_exit_observation.json"
        obs_bytes=root_ob.read_bytes()
        osr=json.loads(obs_bytes)
        require(type(osr["observed_exit_code"]) is int and osr["observed_exit_code"]==0 and osr["case_id"]==c,
                "first-wave completed selector")
        terminal=folder/"terminal.json"
        pinned(terminal,osr["terminal_sha256"])
        plan=Fraction(terminal.stat().st_mtime_ns-(folder/"started.json").stat().st_mtime_ns,10**9)
        times.append({"case_id":c,"selector_observation_sha256":sha256(obs_bytes).hexdigest(),"elapsed_seconds_exact":str(plan)})
    return {"source":"Completed first-wave file timestamps/sizes, no raw numerical packet bytes read",
            "raw_wave_elapsed_seconds_exact":str(elapsed),"raw_part_elapsed_seconds":parts,
            "raw_packet_bytes_first_wave":raw_sizes,
            "raw_packet_bytes_remaining28_if_similar":7*sum(raw_sizes.values()),
            "seven_similar_raw_waves_elapsed_seconds_estimate_exact":str(7*elapsed),
            "selector_times_first_four":times,"seven_similar_four_wide_selector_groups_seconds_estimate_exact":
                 str(7*max(Fraction(x["elapsed_seconds_exact"]) for x in times)),
            "independent_component_hull_replay_timing_assumed":False}


def fixtures(profiles,targets):
    import copy
    for altered in ("repeat","foreign","moved","too_many","payload"):
        p=copy.deepcopy(profiles)
        if altered=="repeat":
            p[0]["case_order"][1]=p[0]["case_order"][0]
        elif altered=="foreign":
            p[0]["case_order"][0]=FIRST[0]
        elif altered=="moved":
            p[0]["case_order"].reverse()
        elif altered=="too_many":
            p[0]["max_owned_native_workers"]=5
        else:
            sample=p[0]["case_order"][0]
            p[0]["targets"][sample]["R_binary_sha256"]="0"*64
        try:
            validate(p,targets)
        except ValueError:
            continue
        raise ValueError(f"unsafe wave profile fixture accepted {altered}")


def prepare():
    for path,digest in PINS.items():
        pinned(path,digest)
    m=get(BRIDGE/"preparation/manifest.json")
    ledger=get(LEDGER)
    original_order=[case for wave in ledger["recommended_batch_waves"] for case in wave["case_ids"]]
    targets=bound_targets(m["cases"])
    profiles=waves(targets,original_order)
    fixtures(profiles,targets)
    counts=Counter(targets[k]["cached_input_kind"] for profile in profiles for k in profile["case_order"])
    result={"status":"FROZEN_COLD_REMAINING28_EXECUTION_PROFILES","source_sha256":sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256":{rel(k):v for k,v in PINS.items()},"first_wave_outside_scope":list(FIRST),
        "first_wave_final_independent_global_outcome_review_required":True,"new_actual_jobs_authorized":False,
        "one_global_spatial_process_limit":4,"later_cases":28,"later_cached_authority_counts":dict(counts),
        "new_spatial_producer_components_planned":4704,"new_independent_component_replays_planned":4704,
        "new_selector_components_planned":4704,"new_integrations_run":0,
        "numerical_payloads_opened":0,"failed_scope_fixtures_rejected":True,"waves":profiles,
        "resource_estimate":timing_metadata()}
    with (HERE/"profiles.json").open("x",encoding="utf8",newline="\n") as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write("\n")
    digest=sha256((HERE/"profiles.json").read_bytes()).hexdigest()
    for row in profiles:
        require(load_profile(HERE/"profiles.json",digest,row["wave_id"])==row,"public profile readback")
    print(digest)


if __name__=="__main__":
    prepare()
