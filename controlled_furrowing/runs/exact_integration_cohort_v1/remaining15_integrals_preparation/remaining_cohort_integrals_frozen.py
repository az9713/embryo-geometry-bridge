"""Prepare and explicitly authorize per-record bulk/surface workers for the remaining cohort."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
import importlib
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
from time import monotonic
from types import SimpleNamespace
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
BASE = "runs/exact_integration_cohort_v1"
OUT = ROOT / BASE / "remaining15_integrals_preparation"
DOMAIN_JOB = BASE + "/remaining15_preparation/domains_job"
ADAPTER = BASE + "/remaining15_scaled_preparation"
SCALED = BASE + "/L_defect_diagnostic/scaled_chart_v1/scaled_element21.py"
AUX_FORCE = BASE + "/L_defect_diagnostic/element21_force.py"
JET = "runs/exact_integration_bridge_v1/safe_kernel_v2/polynomial_jet_safe.py"
RECORDS = ["optoGEF2_small_col_" + c for c in "CDEFGHIJKMNOPQR"]
AXIS = [*range(6), *range(72, 78)]
PARTS = ("bulk", "surface")
PIN = {
    "scripts/exact_scaled_cohort_bulk.py": "6639bb2fae3c3892ec44115ff90e27242586744a2900daa16177cc5a091f1526",
    "scripts/cohort_exact_pilot.py": "a21b680117f236cc38efc46e086be1d9aff80cfb09ae613691b1642db8bed25e",
    AUX_FORCE: "9a5f93bb3d4a39020e28bee165ba606d1d17f63e44a17ee758959aa05109d329",
    ADAPTER + "/protocol.json": "5103212ad5870d370a9fa51000c9bfcdd1496f99652b065622c19c5a262fb656",
    ADAPTER + "/preflight.json": "8cd2dcd90a06b029a4b454b036889943e51649c160a8a0bf762fb8fb2ae37641",
    BASE + "/remaining15_preparation/manifest.json": "93a2f5858a450ca9a8db9f9bc944b918e2a7446abba614f97de5e47fe434495d",
    DOMAIN_JOB + "/results.json": "ecf1fa4ff2dbff45e70116e2ce1f156cf868d333f97963d028a6795ecd8f9c2e",
    DOMAIN_JOB + "/terminal.json": "a1bd414f6a413e6cea710fa86fb1d9b52f014b8040c68a339d3496759e9171ab",
    DOMAIN_JOB + "/started.json": "dd25ebff04f1d68aece549effb48199b99084a9609b9a84c795ea3cadf51e2b8",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return Path(path).relative_to(ROOT).as_posix()


def checked(path, expected, parse=True):
    canonical = PurePosixPath(path)
    require(canonical.as_posix() == path and not canonical.is_absolute()
            and ".." not in canonical.parts and "\\" not in path and ":" not in path
            and (ROOT / canonical).resolve().is_relative_to(ROOT), "canonical local artifact path required")
    raw = (ROOT / canonical).read_bytes()
    require(sha256(raw).hexdigest() == expected, f"frozen bytes changed: {path}")
    return json.loads(raw) if parse else None


def writes(path, value):
    with path.open("x", encoding="utf8", newline="\n") as out:
        json.dump(value, out, default=str, indent=2, allow_nan=False)
        out.write("\n")
    return digest(path)


def utc():
    return datetime.now(timezone.utc).isoformat()


def _yes(checks):
    return isinstance(checks, dict) and bool(checks) and all(x is True for x in checks.values())


def _boot():
    # Verify the complete reviewed project-source closure before importing any of it.
    adapter = checked(ADAPTER + "/protocol.json", PIN[ADAPTER + "/protocol.json"])
    pinned = {**{n: h for n, h in adapter["input_hashes"].items() if n.endswith(".py")},
              **{n: h for n, h in PIN.items() if n.endswith(".py")}}
    for n, h in pinned.items():
        checked(n, h, parse=False)
    I = importlib.import_module("exact_scaled_cohort_bulk")
    E, D, S = I.E, I.D, I.S
    C, F = D.C, S.F
    modules = {
        "scripts/exact_scaled_cohort_bulk.py": I, "scripts/exact_cohort_integrator.py": E,
        "scripts/review_cohort_exact_record.py": D, "scripts/exact_cohort_context.py": C,
        SCALED: S, AUX_FORCE: F, "scripts/exact_uniform_bulk_defect_B_safe.py": E.U,
        "scripts/exact_bulk_bridge_safe.py": E.P, "scripts/exact_fe_target_safe.py": E.T,
        "scripts/algebraic_root_feasibility.py": E.core, JET: E.JET,
        "scripts/exact_fe_moments.py": E.M, "scripts/exact_axis_edge_physical.py": E.AX,
        "scripts/exact_uniform_axis_edge_bridge_B_safe.py": E.AXbridge,
        "scripts/exact_cache_consistency.py": E.U.C, "scripts/exact_fe_target.py": E.U.C.target_core,
        "scripts/exact_uniform_surface_bridge_B_safe.py": E.SUR, "scripts/exact_surface_bridge_safe.py": E.SUR.S,
        "scripts/exact_surface_cache_consistency.py": E.SUR.L,
    }
    for n, module in modules.items():
        require(Path(module.__file__).resolve() == (ROOT / n).resolve()
                and digest(module.__file__) == pinned[n], f"import alias origin changed: {n}")
    require(S.E is F.E is E and S.D is F.D is D and F.C is C and C.T is E.T
            and S.U is E.U and S.M is E.M and S.JET is E.JET and S.core is E.core
            and E.SUR.core is E.core and E.SUR.T is E.T and E.SUR.JET is E.JET and E.SUR.M is E.M
            and E.AXbridge.U is E.U and E.AXbridge.AX is E.AX
            and E.SUR.L.target is E.U.C.target_core and E.T.core is E.core,
            "shared runtime functions resolve to different aliases")
    return SimpleNamespace(I=I, E=E, D=D, C=C, S=S, F=F, pinned=pinned, adapter=adapter)


def _record_dir(record):
    require(record in RECORDS, "explicit remaining record required, B and L excluded")
    return ROOT / BASE / "records" / record


def _terminal(document, expected):
    require(document["status"] == "FINISHED" and type(document["exit_code"]) is int
            and document["exit_code"] == 0 and type(document["pid"]) is int and document["pid"] > 0
            and document["results_sha256"] == expected, "successful matching terminal required")


def _domains(runtime):
    prepped = checked(BASE + "/remaining15_preparation/manifest.json", PIN[BASE + "/remaining15_preparation/manifest.json"])
    completed = checked(DOMAIN_JOB + "/results.json", PIN[DOMAIN_JOB + "/results.json"])
    terminal = checked(DOMAIN_JOB + "/terminal.json", PIN[DOMAIN_JOB + "/terminal.json"])
    started = checked(DOMAIN_JOB + "/started.json", PIN[DOMAIN_JOB + "/started.json"])
    _terminal(terminal, PIN[DOMAIN_JOB + "/results.json"])
    require(terminal["pid"] == started["pid"] and started["manifest_sha256"] == PIN[BASE + "/remaining15_preparation/manifest.json"]
            and started["driver_sha256"] == PIN["scripts/cohort_exact_pilot.py"]
            and completed["status"] == "COMPLETE_15_DOMAINS_ONLY_UNREVIEWED"
            and completed["force_H_integrals_run"] is False
            and [p["record_id"] for p in completed["records"]] == RECORDS
            and [p["record_id"] for p in prepped["records"]] == RECORDS
            and [p["record_id"] for p in runtime.adapter["record_profiles"]] == RECORDS,
            "full terminal remaining15-only domain set required")
    return [(p, result, gate) for p, result, gate in zip(prepped["records"], completed["records"],
                                             started["separately_verified_domains_only_gates"])]


def _record_data(runtime, row, terminal_row, gate):
    E, C, D = runtime.E, runtime.C, runtime.D
    record = row["record_id"]
    original = _record_dir(record) / "original_cube_v1"
    hashes = {value["file"]: value["sha256"] for value in row["files"].values()}
    os_path = f"{DOMAIN_JOB}/{record}.os_terminal.json"
    require(gate["record_id"] == record and gate["file"] == relative(original / "domains_prelaunch_independent_review.json"),
            "cross-record domain prelaunch gate")
    hashes[gate["file"]] = gate["sha256"]
    hashes[relative(original / "domains.json")] = terminal_row["domains_sha256"]
    hashes[relative(original / "domains_terminal.json")] = terminal_row["terminal_sha256"]
    hashes[os_path] = digest(ROOT / os_path)
    docs = {p: checked(p, h) for p, h in hashes.items()}
    protocol = docs[relative(original / "protocol.json")]
    box = docs[relative(original / "coordinate_box.json")]
    domains = docs[relative(original / "domains.json")]
    true_terminal = docs[relative(original / "domains_terminal.json")]
    _terminal(true_terminal, hashes[relative(original / "domains.json")])
    os_terminal = docs[os_path]
    require(domains["record_id"] == os_terminal["record_id"] == record
            and domains["status"] == "PASS_FRESH_DOMAINS_AND_SELECTORS_UNREVIEWED"
            and domains["actual_force_H_integrals_run"] is False
            and domains["protocol_sha256"] == hashes[relative(original / "protocol.json")]
            and domains["coordinate_box_sha256"] == hashes[relative(original / "coordinate_box.json")]
            and true_terminal["pid"] == terminal_row["pid"] == os_terminal["pid"]
            and type(os_terminal["exit_code"]) is int and os_terminal["exit_code"] == 0
            and docs[gate["file"]]["status"] == "PASS" and _yes(docs[gate["file"]]["checks"]),
            "completed original-domain identity or OS process receipt changed")
    body_path = relative(original / "body_geometry.json")
    body = checked(body_path, domains["body_sha256"])
    hashes[body_path] = domains["body_sha256"]
    require(body["status"] == "CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING"
            and len(body["triangles"]) == len(domains["bulk_selectors"]) == 144
            and len(domains["surface_cases"]) == 24, "complete original domains required")
    ctx = C.load_context(record)
    require(ctx.row == protocol["row"] == row["selected_context"]
            and Q(protocol["normal_radius_exact"]) == ctx.old_radius
            and all(D._digest(box[k]) == D._digest(v) for k, v in C.exact_node_box(ctx, ctx.old_radius).items()),
            "unchanged selected record context, C and exact declared cube required")
    for p, h in ctx.input_hashes.items():
        prior = hashes.setdefault(p, h)
        require(prior == h, "conflicting selected context binding")
    from flint import ctx as precision
    precision.prec = 640
    E.core.A.cache_clear()
    nodes = E.SUR.common_nodes(ctx.target, box)
    return SimpleNamespace(ctx=ctx, nodes=nodes, rule=ctx.target.auxiliary_rule(), protocol=protocol,
                           domains=domains, body=body, hashes=hashes, original=original)


def _surface_proof(runtime, data, face, edge, case):
    E, D = runtime.E, runtime.D
    require((case["face"], case["edge"]) == (face, edge) and case["record_id"] == data.ctx.record_id,
            "explicit record surface case required")
    E._surface_rule(data.ctx.target, data.rule, face, edge)
    fields = data.ctx.target.surface(face, edge, nodes=data.nodes)
    require(fields["global_indices"] == case["global_indices"] and len(case["cache_site_domains"]) == 12
            and [p["quadrature_index"] for p in case["cache_site_domains"]] == list(range(12)),
            "same explicit original force/H surface map and cache domains required")
    E.SUR.verify_frozen_domain(fields, case["domain"])
    D._order_history(case["order"], case["order_attempts"], data.protocol["max_surface_order"],
                     Q(data.protocol["force_tail_budget_exact"]), Q(data.protocol["H_entry_tail_budget_exact"]))
    return {"face": face, "edge": edge, "case_sha256": D._digest(case),
            "global_indices": fields["global_indices"], "unchanged_order": case["order"],
            "unchanged_rho_exact": case["domain"]["rho_exact"]}


def _validate_jobs(values, allowed):
    parts = [tuple(x.split(":")) for x in values]
    require(bool(parts) and len(parts) == len(set(parts)) and all(len(job) == 2
            and job[0] in RECORDS and job[1] in PARTS and job in allowed for job in parts),
            "explicit unique authorized remaining record:part jobs required")
    return parts


def _unused(record, part):
    root = _record_dir(record)
    original = root / "original_cube_v1"
    require(not (root / "scaled_integrals_v1" / "actual" / part).exists()
            and not (original / f"{part}_started.json").exists() and not (original / part).exists(),
            "existing record-part output must be preserved")


def _guards():
    eligible = {(RECORDS[0], "bulk"), (RECORDS[0], "surface")}
    values = [RECORDS[0] + ":bulk", RECORDS[0] + ":surface"]
    good = _validate_jobs(values, eligible)
    results = {"complete_explicit_distinct_record_parts": good == [(RECORDS[0], p) for p in PARTS]}
    bad = [
        ("duplicate_part", lambda: _validate_jobs(values[:1] * 2, eligible)),
        ("wrong_record_part", lambda: _validate_jobs([RECORDS[1] + ":bulk"], eligible)),
        ("completed_L", lambda: _validate_jobs(["optoGEF2_small_col_L:bulk"], eligible)),
        ("explicit_only", lambda: _validate_jobs(["all:bulk"], eligible)),
        ("malformed_route", lambda: _validate_jobs([RECORDS[0] + ":bulk:surface"], eligible)),
        ("failed_terminal_rejected", lambda: _terminal({"status": "FAILED", "exit_code": 1, "pid": 1,
                                                         "results_sha256": "a" * 64}, "a" * 64)),
        ("wrong_terminal_digest", lambda: _terminal({"status": "FINISHED", "exit_code": 0, "pid": 1,
                                                     "results_sha256": "a" * 64}, "b" * 64)),
        ("supplemental_alias_digest_rejected", lambda: checked(AUX_FORCE, "0" * 64, False)),
    ]
    for name, check in bad:
        try:
            check()
        except ValueError:
            results[name] = True
        else:
            results[name] = False
    require(_yes(results), "all launch routing and terminal guards required")
    return results


def prepare():
    runtime = _boot()
    packs = _domains(runtime)
    OUT.mkdir(parents=True, exist_ok=True)
    for name in ("protocol.json", "preflight.json", "manifest.json", "started.json", "terminal.json"):
        require(not (OUT / name).exists(), "existing launcher preparation must be preserved")
    require(all(not (_record_dir(r) / "scaled_integrals_v1" / "preparation").exists() for r in RECORDS),
            "earlier per-record launcher preparation must be preserved")
    source = Path(__file__).resolve()
    with (OUT / "remaining_cohort_integrals_frozen.py").open("xb") as file:
        file.write(source.read_bytes())
    psha = writes(OUT / "protocol.json", {
        "status": "FROZEN_REMAINING15_INTEGRAL_LAUNCHER_PREPARATION", "source_file": relative(source),
        "source_sha256": digest(source), "pins": PIN, "transitive_sources": runtime.pinned,
        "records": RECORDS, "precision_bits": 640,
        "scope": "Pure route proof for each successful terminal original domain. Retain original twelve axis kernels and "
                 "twenty-four surface kernels per record. Launch only explicitly authorized record/part pairs after "
                 "bound adapter, per-record actual domain and launcher independent reviews.",
        "authorization": "A launcher review with source_sha256, protocol_sha256, preflight_sha256 and manifest_sha256, "
                         "adapter_review {file,sha256}, and explicit jobs [{record_id,part,domain_review_sha256}].",
        "force_hessian_or_cache_sums_run": False,
    })
    checks = _guards()
    writes(OUT / "started.json", {"pid": os.getpid(), "utc": utc(), "protocol_sha256": psha, "phase": "PURE_ROUTE_PREPARATION"})
    receipts, failures = [], {}
    forbidden = AssertionError("preparation attempted actual integration or cache sums")
    I, E, S = runtime.I, runtime.E, runtime.S
    with (
        patch.object(I, "bulk_component", side_effect=forbidden),
        patch.object(E, "surface_component", side_effect=forbidden),
        patch.object(S.Functional, "__init__", side_effect=forbidden),
        patch.object(E.P.Functional, "__init__", side_effect=forbidden),
        patch.object(E.U, "uniform_cache", side_effect=forbidden),
        patch.object(E.U.C, "exact_bulk_maps", side_effect=forbidden),
    ):
        for row, terminal, gate in packs:
            record = row["record_id"]
            try:
                data = _record_data(runtime, row, terminal, gate)
                d = data.domains
                bulk = [I.prepare_component(data.ctx, data.nodes, data.rule, data.body["triangles"][i],
                                            d["bulk_selectors"][i], i) for i in range(144)]
                axis = [i for i, item in enumerate(bulk) if item["integration_chart"] != "scaled_positive_R_affine_material"]
                require(axis == AXIS and [r["element"] for r in bulk] == list(range(144))
                        and all(item["record_id"] == record for item in bulk), "all 144 routes and twelve original axes required")
                surface = [_surface_proof(runtime, data, *divmod(i, 12), case) for i, case in enumerate(d["surface_cases"])]
                directory = _record_dir(record) / "scaled_integrals_v1" / "preparation"
                directory.mkdir(parents=True, exist_ok=False)
                bsha = writes(directory / "proofs.json", {"bulk": bulk, "surface": surface})
                rsha = writes(directory / "protocol.json", {
                    "status": "FROZEN_RECORD_ROUTE_PROOFS_REVIEW_PENDING", "record_id": record,
                    "batch_protocol_sha256": psha, "source_sha256": digest(source),
                    "domain_and_record_input_hashes": data.hashes, "proofs_sha256": bsha,
                    "independently_reviewed_actual_domains_path": f"{BASE}/remaining15_domains_review/{record}.json",
                    "required_parts": list(PARTS), "original_axis_indices": axis,
                })
                fsha = writes(directory / "preflight.json", {
                    "status": "PASS_ROUTE_PREPARATION_ONLY", "record_id": record, "protocol_sha256": rsha,
                    "source_sha256": digest(source), "proofs_sha256": bsha, "bulk_count": 144, "surface_count": 24,
                    "checks": {"same_original_box_and_rule": True, "all_scaled_and_axis_routes_replayed": True,
                               "all_original_orders_and_tail_budgets_reviewed": True},
                    "force_hessian_or_cache_sums_run": False,
                })
                receipts.append({"record_id": record, "files": {relative(directory / k): v for k, v in (
                    ("protocol.json", rsha), ("preflight.json", fsha), ("proofs.json", bsha))}, "bulk": 144, "surface": 24})
                print(json.dumps({"record": record, "pure_routes_prepared": len(receipts), "total": 15}), flush=True)
            except Exception as error:
                failures[record] = {"type": type(error).__name__, "error": str(error)}
                print(json.dumps({"record": record, **failures[record]}), flush=True)
    ok = not failures and [r["record_id"] for r in receipts] == RECORDS
    checks.update({"all2160_bulk_and360_surface_proofs": ok, "transitive_aliases_verified": True,
                   "no_actual_integral_or_cache_calls": True})
    fsha = writes(OUT / "preflight.json", {"status": "PASS_PREPARATION_ONLY" if ok else "PARTIAL_PRESERVED",
                                             "source_sha256": digest(source), "protocol_sha256": psha,
                                             "checks": checks, "force_hessian_or_cache_sums_run": False})
    msha = writes(OUT / "manifest.json", {
        "status": "READY_FOR_INDEPENDENT_INTEGRAL_PRELAUNCH" if ok else "PARTIAL_PRESERVED",
        "source_sha256": digest(source), "protocol_sha256": psha, "preflight_sha256": fsha,
        "records": receipts, "rejected_attempts": failures, "force_hessian_or_cache_sums_run": False,
    })
    writes(OUT / "terminal.json", {"status": "FINISHED" if ok else "FAILED", "exit_code": 0 if ok else 1,
                                  "pid": os.getpid(), "results_sha256": msha, "utc": utc()})
    print(json.dumps({"prepared": len(receipts), "source_sha256": digest(source), "protocol_sha256": psha,
                      "preflight_sha256": fsha, "manifest_sha256": msha}), flush=True)
    return 0 if ok else 1


def _gate(manifest_sha, review_sha, jobs):
    manifest = checked(relative(OUT / "manifest.json"), manifest_sha)
    psha, prsha = (digest(OUT / n) for n in ("protocol.json", "preflight.json"))
    expect = {"source_sha256": digest(__file__), "protocol_sha256": psha,
              "preflight_sha256": prsha, "manifest_sha256": manifest_sha}
    review = checked(relative(OUT / "prelaunch_independent_review.json"), review_sha)
    require(manifest["status"] == "READY_FOR_INDEPENDENT_INTEGRAL_PRELAUNCH"
            and [p["record_id"] for p in manifest["records"]] == RECORDS
            and review["status"] == "PASS" and _yes(review["checks"])
            and all(review[k] == v for k, v in expect.items())
            and all(manifest[k] == v for k, v in expect.items() if k != "manifest_sha256"),
            "current actual launcher/source/proof manifest review required")
    a = review["adapter_review"]
    approved = checked(a["file"], a["sha256"])
    require(approved["status"].startswith("PASS") and _yes(approved["checks"])
            and approved["source_sha256"] == PIN["scripts/exact_scaled_cohort_bulk.py"]
            and approved["protocol_sha256"] == PIN[ADAPTER + "/protocol.json"]
            and approved["preflight_sha256"] == PIN[ADAPTER + "/preflight.json"],
            "independent reviewed generic adapter required")
    eligible = {(p["record_id"], p["part"]): p["domain_review_sha256"] for p in review["jobs"]}
    require(len(eligible) == len(review["jobs"]), "duplicate independent job authorization")
    select = _validate_jobs(jobs, eligible)
    selected = {}
    for record, part in select:
        _unused(record, part)
        row, = [p for p in manifest["records"] if p["record_id"] == record]
        folder = _record_dir(record) / "scaled_integrals_v1" / "preparation"
        all_docs = {str(Path(p).name): checked(p, h) for p, h in row["files"].items()}
        require(set(all_docs) == {"protocol.json", "preflight.json", "proofs.json"}
                and all(relative(folder / name) in row["files"] for name in all_docs),
                "record-specific preparation paths required")
        p = all_docs["protocol.json"]
        pre = all_docs["preflight.json"]
        require(p["record_id"] == pre["record_id"] == record and p["batch_protocol_sha256"] == psha
                and pre["protocol_sha256"] == row["files"][relative(folder / "protocol.json")]
                and p["proofs_sha256"] == pre["proofs_sha256"] == row["files"][relative(folder / "proofs.json")]
                and pre["status"] == "PASS_ROUTE_PREPARATION_ONLY" and _yes(pre["checks"]), "bound reviewed record proof required")
        review_path = p["independently_reviewed_actual_domains_path"]
        domain_review = checked(review_path, eligible[record, part])
        require(domain_review["status"] == "PASS_RECORD_DOMAIN_REPLAY" and _yes(domain_review["checks"])
                and domain_review["record_id"] == record
                and all(domain_review["input_sha256"].get(k) == v for k, v in p["domain_and_record_input_hashes"].items()
                        if k.startswith(relative(_record_dir(record) / "original_cube_v1")) or k.endswith(".os_terminal.json")),
                "same terminal geometry and independent actual domain replay required")
        selected[record, part] = (p, all_docs["proofs.json"], eligible[record, part])
    return expect, selected


def _worker(manifest_sha, review_sha, job):
    expect, selected = _gate(manifest_sha, review_sha, [job])
    runtime = _boot()
    (record, part), = selected
    protocol, proofs, domain_sha = selected[record, part]
    row, = [row for row in _domains(runtime) if row[0]["record_id"] == record]
    data = _record_data(runtime, *row)
    require(data.hashes == protocol["domain_and_record_input_hashes"], "worker received changed record/domain hashes")
    out = _record_dir(record) / "scaled_integrals_v1" / "actual" / part
    (out / "components").mkdir(parents=True, exist_ok=False)
    pid, start = os.getpid(), monotonic()
    links = {**expect, "prelaunch_sha256": review_sha, "independent_domain_review_sha256": domain_sha,
             "record_id": record, "part": part, "record_proofs_sha256": protocol["proofs_sha256"]}
    writes(out / "started.json", {**links, "pid": pid, "utc": utc()})
    files = []
    try:
        if part == "bulk":
            for i, proof in enumerate(proofs["bulk"]):
                case = data.domains["bulk_selectors"][i]
                packet = runtime.I.bulk_component(data.ctx, data.nodes, data.rule, data.body["triangles"][i], case, i, proof)
                path = out / "components" / f"element_{i:03d}.json"
                files.append({"element": i, "file": relative(path), "sha256": writes(path, packet)})
                print(json.dumps({"record": record, "part": part, "complete": i + 1, "total": 144}), flush=True)
        else:
            for i, proof in enumerate(proofs["surface"]):
                case = data.domains["surface_cases"][i]
                face, edge = divmod(i, 12)
                require(_surface_proof(runtime, data, face, edge, case) == proof, "surface proof changed")
                packet = runtime.E.surface_component(data.ctx, data.nodes, data.rule, case)
                path = out / "components" / f"face{face}_edge{edge:02d}.json"
                files.append({"face": face, "edge": edge, "file": relative(path), "sha256": writes(path, packet)})
                print(json.dumps({"record": record, "part": part, "complete": i + 1, "total": 24}), flush=True)
        require(len(files) == (144 if part == "bulk" else 24), "incomplete record part")
        rsha = writes(out / "results.json", {**links, "status": "COMPONENT_PART_COMPLETE_UNREVIEWED",
                                             "files": files, "seconds": monotonic() - start})
        writes(out / "terminal.json", {"status": "FINISHED", "exit_code": 0, "pid": pid, "utc": utc(),
                                      "results_sha256": rsha})
    except BaseException as error:
        writes(out / "failure.json", {**links, "error": str(error), "completed": files})
        writes(out / "terminal.json", {"status": "FAILED", "exit_code": 1, "pid": pid, "utc": utc()})
        raise


def launch(manifest_sha, review_sha, jobs):
    _, selected = _gate(manifest_sha, review_sha, jobs)
    key = sha256(json.dumps([manifest_sha, review_sha, jobs]).encode()).hexdigest()
    out = OUT / "launches" / key
    out.mkdir(parents=True, exist_ok=False)
    writes(out / "started.json", {"pid": os.getpid(), "utc": utc(), "manifest_sha256": manifest_sha,
                                  "review_sha256": review_sha, "explicit_jobs": jobs})
    completed = []
    try:
        for record, part in selected:
            _unused(record, part)
            label = record + "_" + part
            args = [sys.executable, "-B", "-u", str(Path(__file__).resolve()), "worker", "--manifest-sha256",
                    manifest_sha, "--review-sha256", review_sha, "--job", record + ":" + part]
            with (
                (out / f"{label}.stdout.log").open("xb") as stdout,
                (out / f"{label}.stderr.log").open("xb") as stderr,
            ):
                process = subprocess.Popen(args, cwd=ROOT, stdout=stdout, stderr=stderr)
                writes(out / f"{label}.started.json", {"pid": process.pid, "utc": utc(), "command": args})
                print(json.dumps({"record": record, "part": part, "started_pid": process.pid}), flush=True)
                code = process.wait()
            writes(out / f"{label}.os_terminal.json", {"pid": process.pid, "utc": utc(), "exit_code": code})
            require(code == 0, f"{record}:{part} child failed; retained OS receipt and logs")
            folder = _record_dir(record) / "scaled_integrals_v1" / "actual" / part
            terminal = json.loads((folder / "terminal.json").read_bytes())
            result_hash = digest(folder / "results.json")
            _terminal(terminal, result_hash)
            require(terminal["pid"] == process.pid, "child terminal PID differs from OS receipt")
            completed.append({"record_id": record, "part": part, "pid": process.pid, "exit_code": 0,
                              "terminal_file": relative(folder / "terminal.json"),
                              "terminal_sha256": digest(folder / "terminal.json"), "results_sha256": result_hash})
        msha = writes(out / "results.json", {"status": "EXPLICIT_RECORD_PARTS_COMPLETE_UNREVIEWED", "jobs": completed})
        writes(out / "terminal.json", {"status": "FINISHED", "exit_code": 0, "pid": os.getpid(),
                                      "results_sha256": msha, "utc": utc()})
    except BaseException as error:
        writes(out / "failure.json", {"error": str(error), "completed": completed})
        writes(out / "terminal.json", {"status": "FAILED", "exit_code": 1, "pid": os.getpid(), "utc": utc()})
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest="command", required=True)
    command.add_parser("prepare")
    for name in ("launch", "worker"):
        run = command.add_parser(name)
        run.add_argument("--manifest-sha256", required=True)
        run.add_argument("--review-sha256", required=True)
        run.add_argument("--job", required=True, action="append")
    args = parser.parse_args()
    if args.command == "prepare":
        raise SystemExit(prepare())
    if args.command == "launch":
        launch(args.manifest_sha256, args.review_sha256, args.job)
    else:
        require(len(args.job) == 1, "one exact worker target required")
        _worker(args.manifest_sha256, args.review_sha256, args.job[0])
