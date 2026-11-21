"""Pin true per-part exits, then run the frozen same-record direct and independent replays."""
import argparse
from datetime import datetime,timezone
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
BASE=OUT.parent.relative_to(ROOT).as_posix()
LAUNCH=ROOT/BASE/"remaining15_integrals_preparation"
REDUCE=ROOT/BASE/"remaining15_reduction_preparation"
SOURCE=ROOT/"scripts/exact_cohort_record_reduction.py"
REPLAY=REDUCE/"independent_actual/replay.py"
GATE="66c9a028bc45c204d5c889c4c698fef28511c63ab9df6bddf0cd02e8452d425a"
PROTOCOL="a4bb9eccf03ce91a32eae4e4dfe086f933e867d4ccd7b88a9679c331ec4d8c09"
PREFLIGHT="a1efe382446c96d222a1dd37deb5d78046179c4ca2bd4651088255105f666a58"
RECORDS=tuple("optoGEF2_small_col_"+x for x in "EFGHIJKMNOPQR")


def require(ok, message):
    if not ok:
        raise ValueError(message)


def utc():
    return datetime.now(timezone.utc).isoformat()


def get(p):
    raw=p.read_bytes()
    return json.loads(raw),sha256(raw).hexdigest()


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def rel(p):
    return p.relative_to(ROOT).as_posix()


def save(p, d):
    raw=(json.dumps(d,indent=2,allow_nan=False)+"\n").encode()
    with p.open("xb") as f:
        f.write(raw)
    return sha256(raw).hexdigest()


def pinned():
    pins={
        SOURCE:"8e3177434299fbffeee2fa5161dd4025fc0e9cd0b7f46be85fceaeccae204ed2",
        REDUCE/"prelaunch_independent_review.json":GATE,
        REPLAY:"d60b0a952ea55938f70f9c67c85207a3208eb888a904e5b063f05ad3c4d975bd",
        REPLAY.parent/"protocol.json":PROTOCOL,
        REPLAY.parent/"preflight.json":PREFLIGHT,
        LAUNCH/"prelaunch_independent_review.json":"e89922fe616b089f1723757b7b4f79f2d01da8b7675b32d62c81e17bbb1f0542",
        LAUNCH/"manifest.json":"e4cbc8826b763d235e24ff827d684a00526b79e18bef63b51bc74f727dcdd363",
    }
    require(all(digest(k)==v for k,v in pins.items()),"current frozen reducer/replay and launch gate required")
    return {rel(p):h for p,h in pins.items()}


def part(record, name):
    label=f"{record}_{name}"
    candidates=list((LAUNCH/"launches").glob(f"*/{label}.os_terminal.json"))
    if not candidates:
        return None
    ok=[]
    for p in candidates:
        terminal, terminal_sha=get(p)
        if type(terminal["exit_code"]) is int and terminal["exit_code"]==0:
            ok.append((p,terminal,terminal_sha))
    if not ok:
        return None
    require(len(ok)==1,"ambiguous successful OS part")
    receipt,os_doc,os_sha=ok[0]
    folder=receipt.parent
    starter,starter_sha=get(folder/"started.json")
    lead,lead_sha=get(folder/f"{label}.started.json")
    terminal_root=ROOT/BASE/"records"/record/"scaled_integrals_v1"/"actual"/name
    term,term_sha=get(terminal_root/"terminal.json")
    actual,actual_sha=get(terminal_root/"results.json")
    key=sha256(json.dumps([starter["manifest_sha256"],starter["review_sha256"],starter["explicit_jobs"]]).encode()).hexdigest()
    require(folder.name==key and starter["manifest_sha256"]==digest(LAUNCH/"manifest.json")
        and starter["review_sha256"]==digest(LAUNCH/"prelaunch_independent_review.json")
        and f"{record}:{name}" in starter["explicit_jobs"] and lead["command"][-1]==f"{record}:{name}"
        and type(os_doc["pid"]) is int and os_doc["pid"]>0 and os_doc["pid"]==lead["pid"]==term["pid"]
        and type(term["exit_code"]) is int and term["exit_code"]==0 and term["status"]=="FINISHED"
        and term["results_sha256"]==actual_sha and actual["record_id"]==record and actual["part"]==name
        and actual["status"]=="COMPONENT_PART_COMPLETE_UNREVIEWED" and actual["prelaunch_sha256"]==starter["review_sha256"],
        f"{record}:{name} wrong OS receipt or manifest")
    files=actual["files"]
    expected = ([f"element_{i:03d}.json" for i in range(144)] if name=="bulk"
                else [f"face{f}_edge{e:02d}.json" for f in (0,1) for e in range(12)])
    require([x["file"] for x in files]==[rel(terminal_root/"components"/x) for x in expected]
           and all(len(x["sha256"])==64 for x in files),"canonical component manifest")
    return {"terminal_sha256":term_sha,"manifest_sha256":actual_sha,"os_terminal_sha256":os_sha,
            "supervisor_dir":rel(folder),"pid":os_doc["pid"],"components":len(files),
            "os_terminal_file":rel(receipt),"start_digests":{rel(folder/"started.json"):starter_sha,
              rel(folder/f"{label}.started.json"):lead_sha}}


def run(cmd, receipts, key, receipt_dir):
    args=[sys.executable,"-B","-u",*cmd]
    with subprocess.Popen(args,cwd=ROOT) as proc:
        print(json.dumps({"starting":key,"pid":proc.pid}),flush=True)
        start=utc()
        status=proc.wait()
    receipt={"command":args,"pid":proc.pid,"exit_code":status,"started_utc":start,"terminal_utc":utc()}
    receipts[key]={"file":rel(receipt_dir/f"{key}.os_terminal.json"),
                   "sha256":save(receipt_dir/f"{key}.os_terminal.json",receipt)}
    require(status==0,f"true OS exit {status} for {key}; evidence preserved")


def main(args):
    source_pins=pinned()
    selected=tuple("optoGEF2_small_col_"+x for x in args.record_letters)
    require(len(set(selected))==len(selected) and selected and all(x in RECORDS for x in selected),
            "explicit remaining record letters; C and D preserved")
    tag="".join(x[-1] for x in selected)
    records=[]
    incomplete=[]
    for r in selected:
        p=REDUCE/r
        require(not p.exists() and not (REDUCE/f"{r}_input_admission.json").exists() and not (REPLAY.parent/f"{r}.json").exists(),
                f"{r} has prior production or independent output")
        parts={name:part(r,name) for name in ("bulk","surface")}
        if all(parts.values()):
            records.append({"record_id":r,"parts":parts})
        else:
            incomplete.append(r)
    require(bool(records),"no newly complete selected records")
    obs=OUT/f"{tag}_terminal_observation.json"
    osh=save(obs,{"status":"COMPLETE_RECORD_PARTS_OBSERVED_BEFORE_ANY_REDUCTION_PAYLOAD_READ",
                  "utc":utc(),"records":records,"incomplete_selected_records":incomplete,
                  "current_source_pins":source_pins,"observer_source_sha256":digest(Path(__file__)),
                  "component_payloads_opened":False})
    print(json.dumps({"observation_file":rel(obs),"observation_sha256":osh,
                      "ready":[x["record_id"] for x in records],"incomplete":incomplete}),flush=True)
    dest=OUT/f"{tag}_reducer_os_receipts"
    dest.mkdir(exist_ok=False)
    children, results={},[]
    for target in records:
        r=target["record_id"]
        code=r[-1]
        admit=[rel(SOURCE),"admit",r,"--review-sha256",GATE]
        for p,v in target["parts"].items():
            for k in ("terminal_sha256","manifest_sha256","os_terminal_sha256","supervisor_dir"):
                admit.extend(["--"+p+"-"+k.replace("_","-"),v[k]])
        run(admit,children,f"{code}_admit",dest)
        admission=digest(REDUCE/f"{r}_input_admission.json")
        run([rel(SOURCE),"reduce",r,"--admission-sha256",admission],children,f"{code}_reduce",dest)
        res,rsha=get(REDUCE/r/"results.json")
        terminal,tsha=get(REDUCE/r/"terminal.json")
        reducer_exit,_=get(dest/f"{code}_reduce.os_terminal.json")
        require(terminal["status"]=="FINISHED" and terminal["exit_code"]==0
                and reducer_exit["pid"]==terminal["pid"] and terminal["results_sha256"]==rsha
                and res["admission_sha256"]==admission and res["record_id"]==r,"terminal and result pin")
        run([rel(REPLAY),"replay",r,"--protocol-sha256",PROTOCOL,"--preflight-sha256",PREFLIGHT,
             "--observation",rel(obs),"--observation-sha256",osh,"--admission-sha256",admission,
             "--terminal-sha256",tsha,"--results-sha256",rsha],children,f"{code}_independent",dest)
        v,h=get(REPLAY.parent/f"{r}.json")
        entry={"record_id":r,"producer":{"admission_sha256":admission,"terminal_sha256":tsha,"results_sha256":rsha},
               "independent":{"file":rel(REPLAY.parent/f"{r}.json"),"sha256":h,
                    "status":v["status"],"root_certified":v["root_certified"],
                    "strict_local_minimum_certified":v["strict_local_minimum_certified"],
                    "RF_entries_containing_zero":v["RF_entries_containing_zero"],"exact_scalars":v["exact_scalars"]}}
        results.append(entry)
        save(dest/f"{code}_result_binding.json",entry)
    h=save(dest/"batch_readback.json",{"status":"FINISHED_PER_RECORD_DIRECT_REDUCTION_AND_INDEPENDENT_READBACK",
        "observation":{"file":rel(obs),"sha256":osh},"children":children,"records":results,"utc":utc()})
    print(json.dumps({"batch_readback_sha256":h,"file":rel(dest/"batch_readback.json"),"records":len(results)}),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("record_letters")
    main(parser.parse_args())

