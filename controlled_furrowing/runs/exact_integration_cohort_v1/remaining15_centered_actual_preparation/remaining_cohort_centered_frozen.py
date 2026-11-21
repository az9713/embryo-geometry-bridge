"""Prepare and supervise explicit centered force jobs on authenticated record boxes."""
import argparse
from datetime import datetime, timezone
from hashlib import sha256
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
from time import monotonic
from unittest.mock import patch

from flint import arb, ctx as precision


ROOT = Path(__file__).resolve().parents[1]
BASE = "runs/exact_integration_cohort_v1"
PREPARED = BASE + "/remaining15_centered_preparation"
OUT = ROOT / BASE / "remaining15_centered_actual_preparation"
HELPER = "scripts/exact_centered_cohort_bulk.py"
RECORDS = ("optoGEF2_small_col_C", "optoGEF2_small_col_D")
PIN = {
    HELPER: "524c22266a7e7e8138c06c8e4b444f44bac887d8404296d67d58a3541726d02c",
    PREPARED + "/protocol.json": "dbe461e20104cf60cee9de78010d2c37d086663c420e3fb47869ce5ad2ea339d",
    PREPARED + "/preflight.json": "9c4ea9ed77c22fda2a2b36ada0bc7ac3c1101f100f8d07d68a9a3bbb30af790c",
    PREPARED + "/independent_review.json": "ffcebc86827a4e12cdbec8e0a7fd2b5be84e884eb83f88078d9e61356e7b3769",
    PREPARED + f"/{RECORDS[0]}_proofs.json": "9c3ff54d940735bae7db651056dfbf44c4dd893596520c4dfda2e276cb773d87",
    PREPARED + f"/{RECORDS[1]}_proofs.json": "13d51f9626f125eeb34f429a55cea05298ff2494643da679fc048d66ee772269",
}


def require(ok, why):
    if not ok:
        raise ValueError(why)


def _relative(path):
    return Path(path).relative_to(ROOT).as_posix()


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _checked(path, digest, parse=True):
    raw = (ROOT / path).read_bytes()
    require(sha256(raw).hexdigest() == digest, f"frozen source/input changed: {path}")
    return json.loads(raw) if parse else None


def _utc():
    return datetime.now(timezone.utc).isoformat()


def _boot():
    # Verify the declared project-source closure before loading the helper.
    for path, digest in PIN.items():
        _checked(path, digest, parse=False)
    protocol = _checked(PREPARED + "/protocol.json", PIN[PREPARED + "/protocol.json"])
    sources = {p: h for p, h in protocol["input_hashes"].items() if p.endswith(".py")}
    for path, digest in sources.items():
        require((ROOT / path).resolve().is_relative_to(ROOT), "source escaped the study")
        _checked(path, digest, parse=False)
    H = importlib.import_module("exact_centered_cohort_bulk")
    require(Path(H.__file__).resolve() == (ROOT / HELPER).resolve() and _sha(H.__file__) == PIN[HELPER],
            "wrong generic centered module loaded")
    R, op, L, plan, jobs, expect, hashes = H._boot()
    D = op.D
    for path, digest in {**PIN, **sources}.items():
        D._bound(path, digest, hashes, parse=False)
    pre = _checked(PREPARED + "/preflight.json", PIN[PREPARED + "/preflight.json"])
    review = _checked(PREPARED + "/independent_review.json", PIN[PREPARED + "/independent_review.json"])
    require(review["status"] == "PASS" and D._yes(review["checks"])
            and review["source_sha256"] == PIN[HELPER] and protocol["source_sha256"] == PIN[HELPER]
            and review["protocol_sha256"] == pre["protocol_sha256"] == PIN[PREPARED + "/protocol.json"]
            and review["preflight_sha256"] == PIN[PREPARED + "/preflight.json"]
            and pre["status"] == "PASS_REMAINING15_CENTERED_HELPER_PREPARATION_ONLY"
            and D._yes(pre["checks"]) and set(review["proofs_by_record"]) == set(RECORDS)
            and all(review["proofs_by_record"][r] == PIN[PREPARED + f"/{r}_proofs.json"] for r in RECORDS),
            "same independent reviewed centered helper and actual C/D proof inputs required")
    return H, (R, op, L), plan, jobs, expect, hashes


def _precision_guard():
    require(precision.prec == 640, "centered factory, nodeball proof and actual calls require 640-bit precision")


def _factory(H, runtime, plan, jobs, expect, record, proof, hashes):
    require(record in RECORDS and proof["record_id"] == record, "explicit currently reviewed C or D record required")
    R, op, _ = runtime
    precision.prec = 640
    op.E.core.A.cache_clear()
    _precision_guard()
    require(op.E.core.A.cache_info().currsize == 0, "clear cached conversions before authenticating domain")
    domain, actual, routes, packets, retained = H.bind_reviewed_record(
        runtime, plan, jobs, expect, record, proof["actual_direct_review_sha256"], hashes
    )
    ctx = domain.context
    original = op.C.exact_node_box(ctx, ctx.old_radius)
    require(all(op.D._digest(domain.box[key]) == op.D._digest(v) for key, v in original.items())
            and ctx.record_id == record and ctx.row["parameter_cell_sha256"] == proof["closed_cell_sha256"]
            and packets["surface"] == proof["retained_surface_manifest"]
            and retained["upper_Hessian_entries"] == proof["retained_upper_Hessian_entries"] == 6676
            and retained["b_upper_exact"] == proof["retained_fullbox_b_exact"]
            and len(packets["bulk"]) == 144 and len(packets["surface"]) == 24,
            "same authenticated full domain, retained surface and Hessian prerequisites required")
    canonical_nodes = op.E.SUR.common_nodes(ctx.target, domain.box)
    require(len(domain.nodes) == len(canonical_nodes) == 325 and all(len(p) == 2 for p in domain.nodes),
            "complete authenticated domain nodes required")
    all_nodes = list(range(325))
    canonical_hash = op.E._local_node_digest(canonical_nodes, all_nodes)
    require(op.E._local_node_digest(domain.nodes, all_nodes) == canonical_hash,
            "factory returned a changed or enlarged domain")
    return domain, actual, routes, packets, canonical_nodes, canonical_hash


def _node_guard(E, supplied, canonical, ids, expected):
    _precision_guard()
    actual = E._local_node_digest(supplied, ids)
    known = E._local_node_digest(canonical, ids)
    require(actual == known == expected, "changed or enlarged realized local nodeball domain")
    return actual


def _element_guard(H, runtime, domain, canonical, actual, routes, packets, frozen, index, hashes):
    _precision_guard()
    _, op, _ = runtime
    D = op.D
    ctx = domain.context
    raw = packets["bulk"][index]
    route = routes["bulk"][index]
    selector = actual["bulk_selectors"][index]
    proof = frozen["bulk_proofs"][index]
    ids = list(map(int, ctx.target.a["triangles"][index]))
    require(type(index) is int and index == raw["element"] == route["element"] == selector["element"] == proof["element"]
            and proof["record_id"] == route["record_id"] == ctx.record_id
            and ids == route["local_nodes"] == proof["local_nodes"]
            and proof["original_route_proof_sha256"] == D._digest(route)
            and proof["uniform_packet_sha256"] == raw["sha256"]
            and selector["nodal_hulls_sha256"] == route["nodal_hulls_sha256"],
            "same exact record/element/uniform route proof required")
    realized = _node_guard(op.E, domain.nodes, canonical, ids, route["nodal_hulls_sha256"])
    uniform = D._bound(raw["file"], raw["sha256"], hashes)
    require(uniform["nodal_hulls_sha256"] == realized
            and uniform["record_component_proof_sha256"] == D._digest(route), "actual uniform source domain differs")
    reconstructed = H.prepare_component(runtime, domain, selector, route, uniform, raw["sha256"], index)
    require(D._digest(reconstructed) == D._digest(proof), "rederived point, H segment and gain proof differs from frozen proof")
    receipt = {
        "record_id": ctx.record_id, "element": index, "nodeball_sha256": realized,
        "prior_uniform_file": raw["file"], "prior_uniform_sha256": raw["sha256"],
        "centered_proof_sha256": D._digest(proof), "original_route_sha256": D._digest(route),
        "chart": proof["chart"], "point_force_tail_upper_exact": proof["point_force_tail_max_exact"],
    }
    return receipt, uniform, selector, route, proof


def _jobs(requested, eligible):
    require(bool(requested) and len(requested) == len(set(requested))
            and all(r in RECORDS and r in eligible for r in requested),
            "explicit unique currently reviewed record authorization required")
    return requested


def _tests(runtime):
    _, op, _ = runtime
    E = op.E
    c = [[arb(1), arb(2)]]
    expected = E._local_node_digest(c, [0])
    result = {"rigorous_same_local_box_fingerprint_accepted": _node_guard(E, c, c, [0], expected) == expected}

    def reject(name, act):
        try:
            act()
        except ValueError:
            result[name] = True
        else:
            raise AssertionError(name)

    reject("enlarged_ball_with_same_center_rejected",
           lambda: _node_guard(E, [[arb(1, 1), arb(2)]], c, [0], expected))
    reject("different_center_ball_rejected",
           lambda: _node_guard(E, [[arb(2), arb(2)]], c, [0], expected))
    reject("stale_uniform_nodeball_hash_rejected", lambda: _node_guard(E, c, c, [0], "0" * 64))
    reject("duplicate_record_job", lambda: _jobs([RECORDS[0]] * 2, RECORDS))
    reject("unreviewed_record_job", lambda: _jobs(["optoGEF2_small_col_E"], RECORDS))
    reject("wrong_record_scope", lambda: _jobs([RECORDS[0]], [RECORDS[1]]))
    precision.prec = 512
    try:
        reject("wrong_runtime_precision_rejected_before_use", _precision_guard)
    finally:
        precision.prec = 640
        E.core.A.cache_clear()
    return result


def prepare():
    H, runtime, plan, jobs, expect, hashes = _boot()
    _, op, _ = runtime
    D, E = op.D, op.E
    no = AssertionError("centered actual preparation attempted force, moment, cache or child launch")
    receipts = {}
    with (
        patch.object(H, "bulk_component", side_effect=no),
        patch.object(E.U.C, "exact_bulk_maps", side_effect=no),
        patch.object(E.P.Functional, "__init__", side_effect=no),
        patch.object(E.U, "uniform_cache", side_effect=no),
        patch.object(subprocess, "Popen", side_effect=no),
    ):
        for record in RECORDS:
            path = PREPARED + f"/{record}_proofs.json"
            p = D._bound(path, PIN[path], hashes)
            domain, actual, routes, packets, canonical, key = _factory(H, runtime, plan, jobs, expect, record, p, hashes)
            entries = [_element_guard(H, runtime, domain, canonical, actual, routes, packets, p, i, hashes)[0]
                       for i in range(144)]
            receipts[record] = {
                "record_id": record, "factory_nodeball_sha256": key, "helper_proof_file": path,
                "helper_proof_sha256": PIN[path], "operator_review_sha256": p["actual_direct_review_sha256"],
                "elements": entries, "retained_surface_manifest": p["retained_surface_manifest"],
                "unchanged_full_H_upper_entries": p["retained_upper_Hessian_entries"],
                "unchanged_fullbox_b_exact": p["retained_fullbox_b_exact"],
            }
    checks = _tests(runtime)
    checks.update({
        "both_exact_640bit_factories_checked": set(receipts) == set(RECORDS),
        "all_288_actual_local_ball_fingerprints_and_replayed_exact_proofs": all(
            [r["element"] for r in p["elements"]] == list(range(144)) for p in receipts.values()),
        "all_48_retained_surface_bindings_and_full_H": all(len(p["retained_surface_manifest"]) == 24
            and p["unchanged_full_H_upper_entries"] == 6676 for p in receipts.values()),
        "no_actual_centered_force_sum_or_child_launch": True,
    })
    require(D._yes(checks), "all centered actual caller guards required")
    OUT.mkdir(parents=True, exist_ok=True)
    names = ("remaining_cohort_centered_frozen.py", "protocol.json", "preflight.json", "manifest.json")
    require(not any((OUT / n).exists() for n in names), "preserve write-once centered supervisor preparation")
    with (OUT / names[0]).open("xb") as f:
        f.write(Path(__file__).read_bytes())
    psha = D._write(OUT / names[1], {
        "status": "FROZEN_CENTERED_SUPERVISOR_PREPARATION_REVIEW_PENDING", "records": list(RECORDS),
        "source_file": _relative(__file__), "source_sha256": _sha(__file__), "input_hashes": hashes,
        "requirements": "Set precision 640 and clear core.A before factory bind, authenticate the unchanged exact factory "
                        "Domain, independently derive local fingerprints from supplied boxes versus canonical 640-bit nodes "
                        "and frozen proof/route before each component. Only the two current fully reviewed records are eligible.",
        "scope": "All 144 bulk point forces, frozen H path and affine gain on same original C; keep original 24 surfaces "
                 "and full 6676 upper Hessian for future independent reduction. No new math proof or root assertion here.",
        "gate": "prelaunch_independent_review.json binds source_sha256, protocol_sha256, preflight_sha256, manifest_sha256, "
                "helper_review_sha256 and explicit records. Operator resources are scheduled by the supervising caller.",
        "actual_centered_force_or_child_launch": False,
    })
    fsha = D._write(OUT / names[2], {
        "status": "PASS_CENTERED_ACTUAL_CALLER_PREPARATION_ONLY", "source_sha256": _sha(__file__),
        "protocol_sha256": psha, "checks": checks, "actual_centered_force_or_child_launch": False,
    })
    msha = D._write(OUT / names[3], {
        "status": "READY_FOR_CENTERED_LAUNCHER_INDEPENDENT_PRELAUNCH", "source_sha256": _sha(__file__),
        "protocol_sha256": psha, "preflight_sha256": fsha, "records": receipts,
        "actual_centered_force_or_child_launch": False,
    })
    print(json.dumps({"source_sha256": _sha(__file__), "protocol_sha256": psha, "preflight_sha256": fsha,
                      "manifest_sha256": msha, "checks": len(checks)}), flush=True)


def _gate(D, manifest_sha, review_sha, records):
    p, f, m = (OUT / k for k in ("protocol.json", "preflight.json", "manifest.json"))
    expected = {"source_sha256": _sha(__file__), "protocol_sha256": _sha(p),
                "preflight_sha256": _sha(f), "manifest_sha256": manifest_sha}
    gate = D._bound(_relative(OUT / "prelaunch_independent_review.json"), review_sha, {})
    manifest = D._bound(_relative(m), manifest_sha, {})
    require(gate["status"] == "PASS" and D._yes(gate["checks"])
            and all(gate[k] == v for k, v in expected.items())
            and all(manifest[k] == v for k, v in expected.items() if k != "manifest_sha256")
            and gate["helper_review_sha256"] == PIN[PREPARED + "/independent_review.json"]
            and manifest["status"] == "READY_FOR_CENTERED_LAUNCHER_INDEPENDENT_PRELAUNCH"
            and len(gate["records"]) == len(set(gate["records"])) and set(gate["records"]) <= set(RECORDS),
            "current bound independent centered caller prelaunch required")
    _jobs(records, gate["records"])
    for record in records:
        require(not (OUT / "actual" / record).exists(), "write-once actual centered record already exists")
    return manifest, {**expected, "prelaunch_sha256": review_sha}


def worker(manifest_sha, review_sha, record):
    H, runtime, plan, jobs, expect, hashes = _boot()
    R, op, _ = runtime
    D = op.D
    manifest, links = _gate(D, manifest_sha, review_sha, [record])
    expected = manifest["records"][record]
    path = PREPARED + f"/{record}_proofs.json"
    require(expected["helper_proof_file"] == path and expected["helper_proof_sha256"] == PIN[path],
            "same current record helper proof file required")
    frozen = D._bound(path, PIN[path], hashes)
    domain, actual, routes, packets, canonical, key = _factory(H, runtime, plan, jobs, expect, record, frozen, hashes)
    require(key == expected["factory_nodeball_sha256"] and len(expected["elements"]) == 144
            and expected["operator_review_sha256"] == frozen["actual_direct_review_sha256"]
            and expected["retained_surface_manifest"] == frozen["retained_surface_manifest"],
            "same unchanged independently reviewed factory box and full proof set required")
    out = OUT / "actual" / record
    (out / "components").mkdir(parents=True, exist_ok=False)
    pid, start = os.getpid(), monotonic()
    record_link = {**links, "record_id": record, "factory_nodeball_sha256": key, "operator_review_sha256": expected["operator_review_sha256"]}
    D._write(out / "started.json", {**record_link, "pid": pid, "utc": _utc()})
    completed = []
    try:
        for i in range(144):
            receipt, uniform, s, route, proof = _element_guard(H, runtime, domain, canonical, actual, routes, packets, frozen, i, hashes)
            require(receipt == expected["elements"][i], "changed launcher element proof")
            value = H.bulk_component(runtime, domain, s, route, uniform, receipt["prior_uniform_sha256"], i, proof)
            require(value["record_id"] == record and value["element"] == i and value["new_Hessian_integrations"] is False
                    and value["centered_proof_sha256"] == receipt["centered_proof_sha256"],
                    "returned centered component disagrees with current explicit record")
            f = out / "components" / f"element_{i:03d}.json"
            h = D._write(f, {**value, **record_link, "verified_local_nodeball_sha256": receipt["nodeball_sha256"]})
            completed.append({"element": i, "file": _relative(f), "sha256": h})
            print(json.dumps({"record": record, "complete": i + 1, "total": 144}), flush=True)
        rsha = D._write(out / "results.json", {**record_link, "status": "COMPLETE_144_CENTERED_FORCE_PACKETS_UNREVIEWED",
                                               "files": completed, "retained_surface_manifest": packets["surface"],
                                               "retained_full_H_upper_entries": expected["unchanged_full_H_upper_entries"],
                                               "retained_fullbox_b_exact": expected["unchanged_fullbox_b_exact"],
                                               "new_Hessian_integrations": False, "root_asserted": False,
                                               "seconds": monotonic() - start})
        D._write(out / "terminal.json", {"status": "FINISHED", "exit_code": 0, "pid": pid, "utc": _utc(),
                                         "results_sha256": rsha})
    except BaseException as error:
        D._write(out / "failure.json", {**record_link, "error": str(error), "completed": completed})
        D._write(out / "terminal.json", {"status": "FAILED", "exit_code": 1, "pid": pid, "utc": _utc()})
        raise


def launch(manifest_sha, review_sha, records):
    _, runtime, _, _, _, _ = _boot()
    D = runtime[1].D
    _gate(D, manifest_sha, review_sha, records)
    key = sha256(json.dumps([manifest_sha, review_sha, records]).encode()).hexdigest()
    out = OUT / "launches" / key
    out.mkdir(parents=True, exist_ok=False)
    D._write(out / "started.json", {"pid": os.getpid(), "utc": _utc(), "explicit_records": records,
                                    "manifest_sha256": manifest_sha, "review_sha256": review_sha})
    complete = []
    try:
        for record in records:
            args = [sys.executable, "-B", "-u", str(Path(__file__).resolve()), "worker",
                    "--manifest-sha256", manifest_sha, "--review-sha256", review_sha, "--record", record]
            with ((out / f"{record}.stdout.log").open("xb") as stdout,
                  (out / f"{record}.stderr.log").open("xb") as stderr):
                child = subprocess.Popen(args, cwd=ROOT, stdout=stdout, stderr=stderr)
                D._write(out / f"{record}.started.json", {"pid": child.pid, "utc": _utc(), "command": args})
                print(json.dumps({"record": record, "started_pid": child.pid}), flush=True)
                code = child.wait()
            D._write(out / f"{record}.os_terminal.json", {"pid": child.pid, "utc": _utc(), "exit_code": code})
            require(code == 0, "centered child failed; inspect retained OS receipt and logs")
            location = OUT / "actual" / record
            terminal = json.loads((location / "terminal.json").read_bytes())
            rsha = _sha(location / "results.json")
            D._terminal(terminal, rsha)
            require(terminal["pid"] == child.pid, "actual centered OS/worker PID mismatch")
            complete.append({"record_id": record, "pid": child.pid, "exit_code": 0, "manifest_sha256": rsha,
                             "terminal_sha256": _sha(location / "terminal.json")})
        sha = D._write(out / "results.json", {"status": "EXPLICIT_CENTERED_RECORDS_COMPLETE_UNREVIEWED", "records": complete})
        D._write(out / "terminal.json", {"status": "FINISHED", "exit_code": 0, "pid": os.getpid(),
                                         "results_sha256": sha, "utc": _utc()})
    except BaseException as error:
        D._write(out / "failure.json", {"error": str(error), "completed": complete})
        D._write(out / "terminal.json", {"status": "FAILED", "exit_code": 1, "pid": os.getpid(), "utc": _utc()})
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("prepare")
    for mode in ("launch", "worker"):
        call = commands.add_parser(mode)
        call.add_argument("--manifest-sha256", required=True)
        call.add_argument("--review-sha256", required=True)
        call.add_argument("--record", action="append", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare()
    elif args.command == "launch":
        launch(args.manifest_sha256, args.review_sha256, args.record)
    else:
        require(len(args.record) == 1, "one explicit worker record required")
        worker(args.manifest_sha256, args.review_sha256, args.record[0])
