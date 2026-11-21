"""Read back original-cube domain geometry/selectors. Run prepare, replay, then seal."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
from time import perf_counter
import traceback

HERE = Path(__file__).resolve().parent
COHORT = HERE.parent
ROOT = COHORT.parents[1]
BATCH = COHORT / "remaining15_preparation"
JOB = BATCH / "domains_job"
LREVIEW = COHORT / "driver_independent_review"
DRIVER = ROOT / "scripts/cohort_exact_pilot.py"
DRIVER_SHA = "a21b680117f236cc38efc46e086be1d9aff80cfb09ae613691b1642db8bed25e"
PINS = {
    DRIVER: DRIVER_SHA,
    JOB / "results.json": "ecf1fa4ff2dbff45e70116e2ce1f156cf868d333f97963d028a6795ecd8f9c2e",
    JOB / "terminal.json": "a1bd414f6a413e6cea710fa86fb1d9b52f014b8040c68a339d3496759e9171ab",
    BATCH / "manifest.json": "93a2f5858a450ca9a8db9f9bc944b918e2a7446abba614f97de5e47fe434495d",
    BATCH / "independent_prelaunch/batch_independent_review.json": "34bc25ad801e579ff49d12b992b1dff53a1f2d19d7796e060d08d23ef4a59e26",
    LREVIEW / "review_integral_prelaunch.py": "88dd02b67e8005b06d930aa09872de6e2d61ab97c06a19146b0c50b117cdcdad",
}
LETTERS = "CDEFGHIJKMNOPQR"
RECORDS = ["optoGEF2_small_col_" + suffix for suffix in LETTERS]
SCOPE = "Original-cube geometry and scalar tail selectors only. No force or Hessian integration, cache sums, global assembly or root assertion."


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def rel(path):
    return path.relative_to(ROOT).as_posix()


def read(path):
    return json.loads(path.read_bytes())


def write(path, payload):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2)
        stream.write("\n")


def utc():
    return datetime.now(timezone.utc).isoformat()


def yes(checks):
    return bool(checks) and all(value is True for value in checks.values())


def outdir(record):
    return COHORT / "records" / record / "original_cube_v1"


def read_bindings():
    """Bind the upstream producer terminal, domain inputs and record artifacts."""
    bound = {}

    def bind(path, expected=None):
        name = rel(path)
        h = sha(path)
        assert (expected is None or h == expected) and bound.get(name, h) == h, name
        bound[name] = h
        return h

    for path, h in PINS.items():
        bind(path, h)
    bind(JOB / "started.json")
    supervisor = read(JOB / "started.json")
    bind(BATCH / "run_domains.py", supervisor["supervisor_sha256"])
    assert supervisor["driver_sha256"] == DRIVER_SHA
    assert supervisor["manifest_sha256"] == PINS[BATCH / "manifest.json"]
    terminal = read(JOB / "terminal.json")
    run = read(JOB / "results.json")
    assert terminal["status"] == "FINISHED" and terminal["pid"] == supervisor["pid"] == 31848
    assert type(terminal["exit_code"]) is int and terminal["exit_code"] == 0
    assert terminal["results_sha256"] == PINS[JOB / "results.json"]
    assert run["status"] == "COMPLETE_15_DOMAINS_ONLY_UNREVIEWED" and run["force_H_integrals_run"] is False
    manifest = read(BATCH / "manifest.json")
    batch_audit = read(BATCH / "independent_prelaunch/batch_independent_review.json")
    bind(BATCH / "independent_prelaunch/audit.py", batch_audit["audit_source_sha256"])
    assert batch_audit["status"] == "PASS_REMAINING15_DOMAIN_PRELAUNCH_ONLY" and yes(batch_audit["checks"])
    assert batch_audit["integration_authorized"] is False
    assert manifest["expected_record_ids"] == RECORDS
    assert manifest["excluded_completed_records"] == ["optoGEF2_small_col_B", "optoGEF2_small_col_L"]
    collections = (manifest["records"], run["records"], batch_audit["approved_gates"],
                   supervisor["separately_verified_domains_only_gates"])
    assert all([item["record_id"] for item in rows] == RECORDS for rows in collections)
    for row, actual, audit, active in zip(*collections, strict=True):
        record = row["record_id"]
        out = outdir(record)
        for name in ("protocol.json", "preflight.json", "coordinate_box.json"):
            item = row["files"][name]
            assert item["file"] == rel(out / name)
            bind(out / name, item["sha256"])
        protocol = read(out / "protocol.json")
        pre = read(out / "preflight.json")
        assert pre["status"] == "PASS_PREPARATION_ONLY" and yes(pre["checks"])
        assert pre["source_sha256"] == DRIVER_SHA and pre["protocol_sha256"] == bound[rel(out / "protocol.json")]
        assert protocol["record_id"] == record and protocol["row"] == row["selected_context"]
        for file, h in protocol["source_and_input_hashes"].items():
            bind(ROOT / file, h)
        old_gate = out / "domains_prelaunch_independent_review.json"
        assert audit["gate_file"] == active["file"] == rel(old_gate)
        assert bind(old_gate, audit["sha256"]) == active["sha256"]
        gate = read(old_gate)
        assert gate["status"] == "PASS" and yes(gate["checks"])
        assert gate["source_sha256"] == DRIVER_SHA and gate["protocol_sha256"] == pre["protocol_sha256"]
        assert gate["preflight_sha256"] == bound[rel(out / "preflight.json")]
        for path in (JOB / f"{record}.started.json", JOB / f"{record}.os_terminal.json",
                     out / "domains_started.json"):
            bind(path)
        child, os_exit, domain_start = [read(path) for path in (
            JOB / f"{record}.started.json", JOB / f"{record}.os_terminal.json", out / "domains_started.json")]
        assert os_exit == {k: actual[k] for k in os_exit} and child["command"][1:3] == ["-B", "-u"]
        assert Path(child["command"][3]).resolve() == DRIVER and child["command"][4:] == ["domains", record]
        bind(out / "domains_terminal.json", actual["terminal_sha256"])
        bind(out / "domains.json", actual["domains_sha256"])
        domain_terminal = read(out / "domains_terminal.json")
        assert domain_terminal["status"] == "FINISHED"
        assert domain_terminal["pid"] == domain_start["pid"] == actual["pid"] == child["pid"]
        assert all(type(obj["exit_code"]) is int and obj["exit_code"] == 0 for obj in (os_exit, domain_terminal))
        domain = read(out / "domains.json")
        assert domain_terminal["results_sha256"] == actual["domains_sha256"]
        assert domain["status"] == "PASS_FRESH_DOMAINS_AND_SELECTORS_UNREVIEWED"
        assert domain["record_id"] == record and domain["actual_force_H_integrals_run"] is False
        assert domain["protocol_sha256"] == bound[rel(out / "protocol.json")]
        assert protocol["coordinate_box_file"] == rel(out / "coordinate_box.json")
        assert protocol["coordinate_box_sha256"] == domain["coordinate_box_sha256"] == bound[rel(out / "coordinate_box.json")]
        bind(out / "body_geometry.json", domain["body_sha256"])
        assert not any((out / f).exists() for f in ("domains_failure.json", "bulk_started.json", "surface_started.json",
                                                    "integration_prelaunch_independent_review.json"))
    # Pin shared executable source before the domain review imports it.
    for filename in ("rebind_radius_fix.py", "review_prelaunch.py", "review_repaired.py",
                     "check_selector_smoke.py", "check_driver_gate.py"):
        bind(LREVIEW / filename)
    bind(Path(__file__).resolve())
    return dict(sorted(bound.items())), manifest


def prepare():
    bindings, _ = read_bindings()
    protocol = {"status": "FROZEN_REMAINING15_DOMAIN_READBACK", "records": RECORDS,
                "source_sha256": sha(Path(__file__)), "input_sha256": bindings, "scope": SCOPE,
                "prior_geometry_and_selector_premises": "Same reviewed geometry/Arb/polynomial engines, replayed against every input. Exact tail and cache-site arithmetic from the L independent review."}
    write(HERE / "protocol.json", protocol)
    write(HERE / "preflight.json", {
        "status": "PASS_READBACK_PREPARATION_ONLY", "source_sha256": protocol["source_sha256"],
        "protocol_sha256": sha(HERE / "protocol.json"),
        "checks": {"all_15_supervisor_and_child_exit_receipts_bound": True,
                   "all_record_box_body_domains_preflight_protocol_and_prior_gate_hashes_bound": True,
                   "all_frozen_upstream_inputs_current": True, "no_original_integral_job_started": True},
        "geometry_or_selectors_replayed_during_preparation": False,
    })
    print(json.dumps({"status": "FROZEN", "inputs": len(bindings), "protocol": sha(HERE/"protocol.json"),
                      "preflight": sha(HERE/"preflight.json")}), flush=True)


def frozen():
    protocol, pre = [read(HERE / name) for name in ("protocol.json", "preflight.json")]
    assert pre["status"] == "PASS_READBACK_PREPARATION_ONLY" and yes(pre["checks"])
    assert protocol["status"] == "FROZEN_REMAINING15_DOMAIN_READBACK" and protocol["records"] == RECORDS
    assert protocol["source_sha256"] == pre["source_sha256"] == sha(Path(__file__))
    assert pre["protocol_sha256"] == sha(HERE / "protocol.json")
    bindings, manifest = read_bindings()
    assert bindings == protocol["input_sha256"]
    return bindings, manifest


def check_one(row):
    started = perf_counter()
    record = row["record_id"]
    sys.path.insert(0, str(LREVIEW))
    import review_integral_prelaunch as check
    from review_cohort_exact_record import _old_premises

    P = check.P
    assert P.ROOT == ROOT
    check.RECORD = record
    out = outdir(record)
    protocol, box, body, domains = [read(out / f) for f in
                                   ("protocol.json", "coordinate_box.json", "body_geometry.json", "domains.json")]
    P.precision.prec = 640
    P.core.A.cache_clear()
    assert protocol["precision_bits"] == P.E.PRECISION_BITS == 640
    assert protocol["max_bulk_order"] == 32 and protocol["max_surface_order"] == 64
    assert protocol["surface_domain_depth"] == 8
    context = P.C.load_context(record)
    assert context.row == protocol["row"] == row["selected_context"]
    assert protocol["parameter_cell_sha256"] == context.row["parameter_cell_sha256"]
    assert all(protocol["source_and_input_hashes"][path] == h for path, h in context.input_hashes.items())
    assert context.base == {k: Q(v) for k, v in row["exact_base"].items()}
    assert context.dose == Q(row["exact_dose"]) and context.drive == Q(row["exact_drive"])
    for name in ("R", "S", "u"):
        assert sha256(getattr(context, name).tobytes(order="C")).hexdigest() == row[f"{name}_binary_sha256"]
    eta, z, radius, provider = _old_premises(context, protocol, {})
    assert provider == row["audited_radius_and_cached_eta_Z_provider"]
    assert eta == Q(row["exact_cached_eta_upper"]) and z == Q(row["exact_cached_Z_upper"])
    assert radius == context.old_radius == Q(row["selected_corrector_radius_exact"])
    same_box = P.serial(P.C.exact_node_box(context, radius))
    same_box["domain_id"] = "exactFE:" + record + ":original_cube_v1"
    assert same_box == box and box["root_membership_assumed"] is False
    assert all(Q(x) == 0 for x in box["shift_exact"])
    norm = max(sum((abs(Q.from_float(float(v))) for v in Rrow), Q(0)) for Rrow in context.R)
    assert norm == Q(protocol["R_infinity_norm_exact"]) == Q(row["exact_R_infinity_norm"])
    force_incidence, h_incidence = [0]*611, [0]*611
    a = context.target.a
    for topology in (a["triangles"], a["surface_0_edges"], a["surface_1_edges"]):
        for ids in topology:
            _, proj = check.projection(context, ids)
            for j in proj:
                force_incidence[j] += 1
                h_incidence[j] += len(proj)
    assert force_incidence == protocol["force_incidence"] and h_incidence == protocol["H_row_incidence"]
    cf, ch = max(force_incidence), max(h_incidence)
    assert (cf, ch) == (protocol["max_force_incidence"], protocol["max_H_row_incidence"])
    assert min(cf,ch,norm)>0 and protocol["margin_divisor"]==10000
    ft, ht = Q(protocol["force_tail_budget_exact"]), Q(protocol["H_entry_tail_budget_exact"])
    assert ft == (1-z)*radius/(10000*norm*cf) == Q(row["force_tail_budget_exact"])
    assert ht == (1-z)/(10000*norm*ch) == Q(row["H_entry_tail_budget_exact"])
    assert norm*cf*ft == Q(protocol["global_force_tail_budget_exact"])
    assert norm*ch*ht == Q(protocol["global_derivative_tail_budget_exact"])

    nodes = P.S.common_nodes(context.target, box)
    rule = context.target.auxiliary_rule()
    recomputed = P.G.certify_box(box["reference_nodes_exact"], box["current_centers_exact"],
                                 box["triangles"], box["coordinate_radii_exact"])
    assert recomputed == body
    assert body["status"] == "CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING"
    assert body["axis_preserved_for_entire_box"] is True and body["global_injectivity_for_entire_box"] is True
    assert len(body["triangles"]) == len(domains["bulk_selectors"]) == 144
    charts = Counter()
    bulk_attempts = 0
    for i, case in enumerate(domains["bulk_selectors"]):
        steps, chart = check.bulk_check(context, nodes, rule, case, body["triangles"][i], i, ft, ht)
        charts[chart] += 1
        bulk_attempts += steps
    assert dict(charts) == {"off_axis": 132, "axis_vertex": 6, "axis_edge_physical": 6}
    assert len(domains["surface_cases"]) == 24
    surface_attempts = leaves = cache_sites = 0
    for i, case in enumerate(domains["surface_cases"]):
        steps, width = check.surface_check(context, nodes, box, rule, case, i//12, i%12, ft, ht)
        surface_attempts += steps
        leaves += width
        cache_sites += len(case["cache_site_domains"])
    assert cache_sites == 288
    paths = [out / f for f in ("protocol.json", "preflight.json", "coordinate_box.json", "body_geometry.json",
                              "domains.json", "domains_terminal.json", "domains_prelaunch_independent_review.json")]
    paths += [JOB / f"{record}.os_terminal.json"]
    return {"status":"PASS_RECORD_DOMAIN_REPLAY", "record_id":record,
            "checks":{"original_selected_record_context_and_authoritative_box":True,
                      "old_proof_center_cached_contraction_and_exact_incidence_budgets":True,
                      "whole_box_embedding_and_144_J_domains_replayed_with_reviewed_engine":True,
                      "all_144_exact_bulk_selector_budget_attempts_recomputed":True,
                      "all_24_closed_surface_G_domains_replayed_and_hulls_checked":True,
                      "all_24_exact_surface_selector_budget_attempts_recomputed":True,
                      "all_288_exact_unweighted_surface_cache_line_domains_reconstructed":True},
            "same_C_sha256": context.row["parameter_cell_sha256"],
            "proof_center_state":context.row["proof_center_state_file"],
            "selected_cached_bound_provider":provider, "dose_exact":str(context.dose),
            "input_sha256":{rel(path):sha(path) for path in paths},
            "body_J":body["J_bound"], "bulk_charts":dict(charts),
            "bulk_elements":144, "bulk_attempts":bulk_attempts,
            "bulk_selected_orders":sorted({x["degree"] for x in domains["bulk_selectors"]}),
            "surface_edges":24, "surface_attempts":surface_attempts, "surface_leaves":leaves,
            "surface_selected_orders":sorted({x["order"] for x in domains["surface_cases"]}),
            "surface_cache_sites":cache_sites,"seconds":perf_counter()-started,
            "scope":SCOPE}


def replay(workers):
    inputs, manifest = frozen()
    started = perf_counter()
    write(HERE / "started.json", {"pid": os.getpid(), "utc": utc(), "workers": workers,
                                  "protocol_sha256": sha(HERE / "protocol.json"),
                                  "preflight_sha256": sha(HERE / "preflight.json")})
    results = {}
    try:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            jobs = {pool.submit(check_one, row): row["record_id"] for row in manifest["records"]}
            for future in as_completed(jobs):
                record = jobs[future]
                try:
                    item = future.result()
                except Exception:
                    print(f"FAIL {record}", flush=True)
                    raise
                assert all(inputs[p] == h for p,h in item["input_sha256"].items())
                name = HERE / (record + ".json")
                write(name, item)
                results[record] = {"record_id":record,"file":rel(name),"sha256":sha(name),
                                   "bulk_elements":item["bulk_elements"],"surface_edges":item["surface_edges"],
                                   "surface_cache_sites":item["surface_cache_sites"]}
                print(f"PASS {record} {item['seconds']:.2f}s", flush=True)
        frozen()
        rows = [results[k] for k in RECORDS]
        totals = {k:sum(item[k] for item in rows) for k in ("bulk_elements","surface_edges","surface_cache_sites")}
        assert totals == {"bulk_elements":2160,"surface_edges":360,"surface_cache_sites":4320}
        write(HERE / "results.json", {"status":"PASS_REMAINING15_DOMAINS_REPLAY", "scope":SCOPE,
                                     "records":rows,"totals":totals,"seconds":perf_counter()-started})
        write(HERE / "terminal.json", {"status":"FINISHED","pid":os.getpid(),"exit_code":0,"utc":utc(),
                                       "results_sha256":sha(HERE / "results.json")})
        print(json.dumps(totals), flush=True)
    except BaseException:
        write(HERE / "failure.json", {"completed":list(results), "traceback":traceback.format_exc()})
        write(HERE / "terminal.json", {"status":"FAILED","pid":os.getpid(),"exit_code":1,"utc":utc()})
        raise


def seal():
    inputs, _ = frozen()
    terminal, results, started = [read(HERE / f"{k}.json") for k in ("terminal","results","started")]
    assert terminal["status"] == "FINISHED" and terminal["exit_code"] == 0 and terminal["pid"] == started["pid"]
    assert terminal["results_sha256"] == sha(HERE / "results.json")
    assert started["protocol_sha256"] == sha(HERE / "protocol.json")
    assert started["preflight_sha256"] == sha(HERE / "preflight.json")
    assert results["status"] == "PASS_REMAINING15_DOMAINS_REPLAY"
    assert [x["record_id"] for x in results["records"]] == RECORDS
    pending = []
    review_hash = sha(HERE / "results.json")
    for row in results["records"]:
        record = row["record_id"]
        path = ROOT / row["file"]
        assert sha(path) == row["sha256"]
        item = read(path)
        assert item["status"] == "PASS_RECORD_DOMAIN_REPLAY" and item["record_id"] == record and yes(item["checks"])
        assert all(inputs[k] == h for k,h in item["input_sha256"].items())
        out = outdir(record)
        gate_name = out / "integration_prelaunch_independent_review.json"
        assert not gate_name.exists()
        pending.append((gate_name, {
            "status":"PASS","record_id":record,"checks":item["checks"], "source_sha256":DRIVER_SHA,
            "protocol_sha256":inputs[rel(out / "protocol.json")],
            "preflight_sha256":inputs[rel(out / "preflight.json")],
            "domains_result_sha256":inputs[rel(out / "domains.json")],
            "body_sha256": inputs[rel(out / "body_geometry.json")],
            "coordinate_box_sha256": inputs[rel(out / "coordinate_box.json")],
            "support_file":row["file"], "support_sha256":row["sha256"],
            "review_source_sha256":sha(Path(__file__)), "batch_results_sha256":review_hash,
            "scope":SCOPE + " This receipt applies to the original-driver domain gate only.",
        }))
    for path, gate in pending:
        write(path, gate)
    write(HERE / "issued_gates.json", {"status":"PASS_15_ORIGINAL_DRIVER_DOMAIN_GATES",
                                      "results_sha256":review_hash, "totals":results["totals"],
                                      "gates":[{"file":rel(path),"sha256":sha(path)} for path,_ in pending],
                                      "adapter_specific_review_performed":False,
                                      "force_H_integral_jobs_run":False})
    print(sha(HERE / "issued_gates.json"), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("prepare","replay","seal"))
    parser.add_argument("--workers", type=int, choices=(1,2,3), default=3)
    args = parser.parse_args()
    if args.phase == "prepare":
        prepare()
    elif args.phase == "replay":
        replay(args.workers)
    else:
        seal()

