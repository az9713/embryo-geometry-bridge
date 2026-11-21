"""Exact signed direct reduction for one completed original-cube cohort record."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
import importlib
from itertools import combinations_with_replacement
import json
import os
from pathlib import Path, PurePosixPath
from time import monotonic
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
BASE = "runs/exact_integration_cohort_v1"
OUT = ROOT / BASE / "remaining15_reduction_preparation"
LAUNCH = BASE + "/remaining15_integrals_preparation"
LAUNCHER = "scripts/remaining_cohort_integrals.py"
RECORDS = ["optoGEF2_small_col_" + c for c in "CDEFGHIJKMNOPQR"]
PARTS = ("bulk", "surface")
PINS = {
    LAUNCHER: "c9654ab19be69aa70a0100ede0919bc8a87040728ca38487758f55ad3047dcfc",
    LAUNCH + "/protocol.json": "ce18ae97d8e9a31db1bb12460371700eb0b17a5cb90e05daf0248584e7d8585e",
    LAUNCH + "/preflight.json": "a6b5939d7e7aefecc150f45067764e30d87e9b8551355a2b5e223e4d0e992d50",
    LAUNCH + "/manifest.json": "e4cbc8826b763d235e24ff827d684a00526b79e18bef63b51bc74f727dcdd363",
    LAUNCH + "/prelaunch_independent_review.json": "e89922fe616b089f1723757b7b4f79f2d01da8b7675b32d62c81e17bbb1f0542",
    BASE + "/remaining15_scaled_preparation/independent_review.json":
        "821b3730305b1d55da9a1cb119780a1717adf99b953539e2357535d81e597fe2",
}


def require(ok, text):
    if not ok:
        raise ValueError(text)


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _relative(path):
    return Path(path).relative_to(ROOT).as_posix()


def _yes(v):
    return isinstance(v, dict) and bool(v) and all(x is True for x in v.values())


def _record_dir(record):
    require(record in RECORDS, "explicit remaining record required; completed B/L excluded")
    return ROOT / BASE / "records" / record


def _boot():
    require(_sha(ROOT / LAUNCHER) == PINS[LAUNCHER], "frozen launcher changed before import")
    L = importlib.import_module("remaining_cohort_integrals")
    require(Path(L.__file__).resolve() == (ROOT / LAUNCHER).resolve()
            and _sha(L.__file__) == PINS[LAUNCHER], "wrong launcher import path")
    runtime = L._boot()
    copies = {}
    docs = {k: runtime.D._bound(k, h, copies, parse=k.endswith(".json")) for k, h in PINS.items()}
    plan, gate = docs[f"{LAUNCH}/manifest.json"], docs[f"{LAUNCH}/prelaunch_independent_review.json"]
    expect = {"source_sha256": PINS[LAUNCHER], "protocol_sha256": PINS[f"{LAUNCH}/protocol.json"],
              "preflight_sha256": PINS[f"{LAUNCH}/preflight.json"], "manifest_sha256": PINS[f"{LAUNCH}/manifest.json"],
              "prelaunch_sha256": PINS[f"{LAUNCH}/prelaunch_independent_review.json"]}
    require(plan["status"] == "READY_FOR_INDEPENDENT_INTEGRAL_PRELAUNCH"
            and gate["status"] == "PASS" and _yes(gate["checks"])
            and all(gate[k] == v for k, v in expect.items() if k != "prelaunch_sha256")
            and [r["record_id"] for r in plan["records"]] == RECORDS
            and docs[gate["adapter_review"]["file"]]["status"] == "PASS"
            and gate["adapter_review"]["sha256"] == PINS[gate["adapter_review"]["file"]],
            "same independently approved remaining15 launcher and adapter required")
    jobs = {(j["record_id"], j["part"]): j["domain_review_sha256"] for j in gate["jobs"]}
    require(len(jobs) == len(gate["jobs"]) == 30 and set(jobs) == {(r, p) for r in RECORDS for p in PARTS},
            "all 30 distinct independently eligible record-part jobs required")
    copies.update(runtime.pinned)
    return runtime, plan, jobs, expect, copies


def _load_record(runtime, plan, jobs, expect, record, hashes):
    D = runtime.D
    directory = _record_dir(record)
    original = directory / "original_cube_v1"
    stage = directory / "scaled_integrals_v1"
    row, = [r for r in plan["records"] if r["record_id"] == record]
    expected_names = {_relative(stage / "preparation" / n) for n in ("preflight.json", "protocol.json", "proofs.json")}
    require(set(row["files"]) == expected_names, "same canonical route-proof files required")
    prepared = {Path(path).name: D._bound(path, h, hashes) for path, h in row["files"].items()}
    protocol, pre, routes = (prepared[k] for k in ("protocol.json", "preflight.json", "proofs.json"))
    require(pre["status"] == "PASS_ROUTE_PREPARATION_ONLY" and _yes(pre["checks"])
            and protocol["record_id"] == pre["record_id"] == record
            and protocol["batch_protocol_sha256"] == expect["protocol_sha256"]
            and protocol["proofs_sha256"] == pre["proofs_sha256"] == row["files"][_relative(stage / "preparation/proofs.json")]
            and [p["element"] for p in routes["bulk"]] == list(range(144)) and len(routes["surface"]) == 24,
            "complete same record original route proofs required")
    review_path = protocol["independently_reviewed_actual_domains_path"]
    review_sha = jobs[record, "bulk"]
    require(review_sha == jobs[record, "surface"], "bulk and surface must share the one reviewed original domain")
    review = D._bound(review_path, review_sha, hashes)
    dom_hashes = protocol["domain_and_record_input_hashes"]
    for path, digest in dom_hashes.items():
        D._bound(path, digest, hashes, parse=False)
    require(review["status"] == "PASS_RECORD_DOMAIN_REPLAY" and _yes(review["checks"])
            and review["record_id"] == record
            and all(review["input_sha256"].get(n) == h for n, h in dom_hashes.items()
                    if n.startswith(_relative(original)) or n.endswith(".os_terminal.json")),
            "actual independent domain descriptor or source bindings changed")
    original_protocol, box, domains, body = [
        D._bound(_relative(original / n), dom_hashes[_relative(original / n)], hashes)
        for n in ("protocol.json", "coordinate_box.json", "domains.json", "body_geometry.json")
    ]
    ctx = runtime.C.load_context(record)
    require(ctx.row == original_protocol["row"]
            and all(dom_hashes.get(n) == h for n, h in ctx.input_hashes.items())
            and domains["record_id"] == record and domains["body_sha256"] == dom_hashes[_relative(original / "body_geometry.json")]
            and domains["coordinate_box_sha256"] == dom_hashes[_relative(original / "coordinate_box.json")]
            and original_protocol["coordinate_box_file"] == _relative(original / "coordinate_box.json")
            and Q(original_protocol["normal_radius_exact"]) == ctx.old_radius
            and all(D._digest(box[k]) == D._digest(v) for k, v in runtime.C.exact_node_box(ctx, ctx.old_radius).items()),
            "same exact old proof center, predictor, R, full closed C and nodal box required")
    links = {**expect, "record_id": record, "independent_domain_review_sha256": review_sha,
             "record_proofs_sha256": protocol["proofs_sha256"]}
    return ctx, original_protocol, domains, body, routes, links, stage


def _manifest(document, links, part, components):
    require(document["status"] == "COMPONENT_PART_COMPLETE_UNREVIEWED" and document["part"] == part
            and all(document.get(k) == v for k, v in links.items()), "mixed record, proof or prelaunch manifest")
    files = document["files"]
    if part == "bulk":
        expected = [("element", i, f"{components}/element_{i:03d}.json") for i in range(144)]
        require(len(files) == 144 and all(type(row.get(k)) is int and row[k] == i and row["file"] == path
                    for row, (k, i, path) in zip(files, expected)), "complete unique canonical 144-bulk manifest required")
    else:
        expected = [(face, edge, f"{components}/face{face}_edge{edge:02d}.json") for face in (0, 1) for edge in range(12)]
        require(len(files) == 24 and all(type(row.get("face")) is type(row.get("edge")) is int
                    and (row["face"], row["edge"], row["file"]) == entry for row, entry in zip(files, expected)),
                "complete unique canonical 24-surface manifest required")
    require(len({p["file"] for p in files}) == len(files) and all(
        type(p["sha256"]) is str and len(p["sha256"]) == 64 and set(p["sha256"]) <= set("0123456789abcdef")
        for p in files), "one exact digest for every distinct owned packet required")
    return files


def _terminal_part(D, stage, links, record, part, expected, bindings):
    folder = stage / "actual" / part
    terminal = D._bound(_relative(folder / "terminal.json"), expected["terminal_sha256"], bindings)
    D._terminal(terminal, expected["manifest_sha256"])
    supervisor = PurePosixPath(expected["supervisor_dir"])
    require(supervisor.parent.as_posix() == LAUNCH + "/launches" and D._sha(supervisor.name) == supervisor.name,
            "canonical launched supervisor directory required")
    label = record + "_" + part
    os_exit = D._bound(f"{supervisor}/{label}.os_terminal.json", expected["os_terminal_sha256"], bindings)
    started_path = f"{supervisor}/{label}.started.json"
    started = D._bound(started_path, _sha(ROOT / started_path), bindings)
    main_path = f"{supervisor}/started.json"
    main = D._bound(main_path, _sha(ROOT / main_path), bindings)
    key = sha256(json.dumps([links["manifest_sha256"], links["prelaunch_sha256"], main["explicit_jobs"]]).encode()).hexdigest()
    require(key == supervisor.name and main["manifest_sha256"] == links["manifest_sha256"]
            and main["review_sha256"] == links["prelaunch_sha256"] and record + ":" + part in main["explicit_jobs"]
            and type(os_exit["exit_code"]) is int and os_exit["exit_code"] == 0
            and os_exit["pid"] == terminal["pid"] == started["pid"]
            and started["command"][-1] == record + ":" + part,
            "same separately observed OS and worker terminal exits required")
    manifest = D._bound(_relative(folder / "results.json"), expected["manifest_sha256"], bindings)
    return _manifest(manifest, links, part, _relative(folder / "components"))


def _prelaunch(D, review_sha, hashes):
    path = _relative(OUT / "prelaunch_independent_review.json")
    gate = D._bound(path, review_sha, hashes)
    p, pre = (OUT / n for n in ("protocol.json", "preflight.json"))
    cases = {"source_sha256": _sha(__file__), "protocol_sha256": _sha(p), "preflight_sha256": _sha(pre)}
    saved = D._bound(_relative(pre), cases["preflight_sha256"], hashes)
    require(gate["status"] == "PASS" and _yes(gate["checks"]) and all(gate[k] == v for k, v in cases.items())
            and saved["status"] == "PASS_GENERIC_DIRECT_REDUCTION_PREPARATION_ONLY"
            and saved["source_sha256"] == cases["source_sha256"] and saved["protocol_sha256"] == cases["protocol_sha256"]
            and _yes(saved["checks"]), "independent reviewed exact reducer prelaunch required")
    D._bound(_relative(p), cases["protocol_sha256"], hashes, parse=False)
    D._bound(_relative(__file__), cases["source_sha256"], hashes, parse=False)


def admit(record, review_sha, pins):
    runtime, plan, jobs, expect, hashes = _boot()
    D = runtime.D
    _prelaunch(D, review_sha, hashes)
    _, _, _, _, _, links, stage = _load_record(runtime, plan, jobs, expect, record, hashes)
    paths = {part: _terminal_part(D, stage, links, record, part, pins[part], hashes) for part in PARTS}
    target = OUT / f"{record}_input_admission.json"
    sha = D._write(target, {"status": "ALL_168_EXPECTED_PACKETS_BOUND_BEFORE_LOCAL_READ",
                           "record_id": record, "review_sha256": review_sha, "external_part_pins": pins,
                           "expected_packets": paths, "source_and_manifest_hashes": hashes,
                           "local_component_payloads_opened": False})
    print(json.dumps({"record_id": record, "admission_sha256": sha, "components": 168}), flush=True)


def _exact_tail(D, packet_bound, limit, enclosure):
    value = D._fraction(packet_bound)
    require(0 <= value <= limit and D._contains(D._packed(enclosure), (value, value)),
            "exact tail exceeds unchanged local budget or its rounded enclosure")
    return value


def _identity(ctx, case, packet, ids):
    row = ctx.row
    expected = {
        "record_id": ctx.record_id, "dose_exact": str(ctx.dose), "precision_bits": 640,
        "parameter_cell_sha256": row["parameter_cell_sha256"], "point_cache_sha256": row["point_cache_sha256"],
        "model_observation_contract_sha256": row["model_observation_contract_sha256"],
        "predictor_sha256": row["predictor_sha256"], "global_indices": ids,
    }
    require(all(packet.get(k) == v for k, v in expected.items())
            and case["record_id"] == ctx.record_id and case["global_indices"] == ids,
            "mixed record, exact dose, C, predictor, cache or same-Y free map")


def _scaled_packet(D, ctx, i, case, body_piece, proof, packet, protocol):
    nodes = ctx.target.a["triangles"][i]
    ids, free = D._mesh_map(ctx, nodes)
    _identity(ctx, case, packet, ids)
    require(all(ctx.target.a["reference"][int(i), 0] > 0 for i in nodes[:3])
            and packet["status"] == "UNIFORM_COHORT_BULK_COMPONENT_UNREVIEWED"
            and packet["chart"] == case["chart"] == "off_axis"
            and packet["integration_chart"] == proof["integration_chart"] == "scaled_positive_R_affine_material"
            and packet["element"] == case["element"] == proof["element"] == i
            and packet["free_local"] == proof["free_local"] == case["free_local"] == free
            and packet["active"] == case["active"] == str(int(ctx.target.a["active"][i, 0]))
            and packet["record_component_proof_sha256"] == D._digest(proof)
            and proof["original_selector_sha256"] == D._digest(case)
            and proof["body_element_sha256"] == D._digest(body_piece)
            and packet["nodal_hulls_sha256"] == case["nodal_hulls_sha256"] == proof["nodal_hulls_sha256"]
            and packet["auxiliary_element_rule_sha256"] == case["auxiliary_element_rule_sha256"] == proof["auxiliary_element_rule_sha256"]
            and packet["degree"] == case["degree"] == proof["unchanged_degree"]
            and packet["rho_exact"] == case["rho_exact"] == proof["unchanged_rho_exact"],
            "same reviewed scaled route, nodal/auxiliary geometry and local ordered projection required")
    J = D._packed(body_piece["J_bound"], dyadic=False)
    rho = max(abs(J[0] - 1), abs(J[1] - 1))
    require(J[0] > Q(1, 10**6) and body_piece["positive_J_for_entire_box"] is True
            and list(map(Q, case["J_bounds_exact"])) == list(J) and Q(case["raw_rho_exact"]) == rho,
            "same independently reviewed analytic material J bounds required")
    D._rho(Q(case["rho_exact"]), rho)
    fmax, hmax = Q(protocol["force_tail_budget_exact"]), Q(protocol["H_entry_tail_budget_exact"])
    D._order_history(case["degree"], case["order_attempts"], protocol["max_bulk_order"], fmax, hmax)
    ft = packet["force_two_mass_tail_upper_exact"]
    fp = packet["force_two_mass_tail_upper"]
    hp = packet["H_integral_minus_cache_upper"]
    require(len(ft) == len(fp) == len(ids) and len(hp) == len(ids) * (len(ids) + 1) // 2,
            "full local scaled tail and Hessian coverage required")
    f = [_exact_tail(D, a, fmax, b) for a, b in zip(ft, fp)]
    h = [_exact_tail(D, p["two_mass_tail_upper_exact"], hmax, p["two_mass_tail_upper"]) for p in hp]
    require(max(f) == Q(proof["scaled_force_tail_max_exact"]) and max(h) == Q(proof["scaled_H_entry_tail_max_exact"]),
            "scaled exact tail maxima differ from the bound route proof")
    cache = packet["cache_consistency"]
    require(cache["all64_convex_strain_lines_positive"] is True
            and D._packed(cache["minimum_cache_convex_line_J_lower"])[0] > Q(1, 10**6)
            and D._packed(cache["max_coefficientwise_delta_q_upper"])[0] >= 0,
            "same-node cache strain segments must stay in reviewed domain")
    return ids


def _reduce(D, ctx, protocol, domains, body, routes, parts, hashes):
    f, H, counts, hcounts = [D.ZERO] * 611, {}, [0] * 611, [0] * 611
    for part, cases in (("bulk", domains["bulk_selectors"]), ("surface", domains["surface_cases"])):
        for i, (case, item) in enumerate(zip(cases, parts[part])):
            packet = D._bound(item["file"], item["sha256"], hashes)
            proof = routes[part][i]
            if part == "bulk" and proof["integration_chart"] == "scaled_positive_R_affine_material":
                ids = _scaled_packet(D, ctx, i, case, body["triangles"][i], proof, packet, protocol)
            else:
                ids = D._packet_meta(ctx, part, i, case, packet, body, protocol)
                if part == "bulk":
                    require(packet["record_component_proof_sha256"] == D._digest(proof)
                            and packet["integration_chart"] == proof["integration_chart"] == case["chart"]
                            and case["chart"] in ("axis_vertex", "axis_edge_physical")
                            and proof["body_element_sha256"] == D._digest(body["triangles"][i]),
                            "same reviewed axis route required")
                else:
                    require(proof["case_sha256"] == D._digest(case) and proof["global_indices"] == ids,
                            "same independent surface proof required")
            D._scatter(f, H, ids, packet["force_integral_minus_cache"], packet["H_integral_minus_cache_upper"])
            for j in ids:
                counts[j] += 1
                hcounts[j] += len(ids)
    require(len(parts["bulk"]) == 144 and len(parts["surface"]) == 24 and len(H) == 6676
            and min(counts) > 0 and counts == protocol["force_incidence"] and hcounts == protocol["H_row_incidence"],
            "strict complete 168-component force/Hessian topology required")
    Hrows = D._h_rows(611, H)
    RF, b_rows, norm = D._dyadic_reduce(ctx.R, f, Hrows)
    require(norm == Q(protocol["R_infinity_norm_exact"]), "same exact binary preconditioner norm required")
    d, b = max(map(D._abs_upper, RF)), max(b_rows)
    eta0, z0, r, provider = D._old_premises(ctx, protocol, hashes)
    eta, z, lhs, passes, error = D._strict_test(eta0, z0, d, b, r)
    curvature = D._minimum(ctx, hashes, z0, z, r, norm, provider) if passes else None
    return {
        "status": "CONDITIONAL_UNREVIEWED_DIRECT_SUFFICIENT_TEST_PASSES" if passes
                  else "UNREVIEWED_DIRECT_SUFFICIENT_TEST_UNRESOLVED",
        "record_id": ctx.record_id, "component_count": 168, "force_entries": 611, "upper_Hessian_entries": len(H),
        "exact_source_dose": ctx.dose, "same_closed_C_sha256": ctx.row["parameter_cell_sha256"],
        "global_F_intervals_exact": f, "global_H_upper_intervals_exact": [[i, j, *H[i, j]] for i, j in sorted(H)],
        "preconditioned_RF_intervals_exact": RF, "H_absolute_row_bounds_exact": Hrows,
        "abs_R_H_row_bounds_exact": b_rows, "R_infinity_norm_exact": norm,
        "RF_defect_uniform_upper_exact": d, "b_upper_exact": b, "old_cached_eta_upper_exact": eta0,
        "old_cached_Z_upper_exact": z0, "old_radius_exact": r, "old_cached_bound_provider": provider,
        "total_eta_upper_exact": eta, "total_Z_upper_exact": z, "strict_selfmap_lhs_exact": lhs,
        "sufficiency_passed": passes, "conditional_root_error_exact": error, "conditional_curvature": curvature,
        "root_assertion_independently_reviewed": False,
        "scope": "Direct sufficient test on the unchanged original cube and closed C. Independent actual packet and exact "
                 "reduction review required before claiming a root. A failed bound leaves sufficiency unresolved.",
        "method": "Original reviewed signed interval scatter and integer-dyadic R force/absolute Hessian row kernels.",
    }


def reduce(record, admission_sha):
    runtime, plan, jobs, expect, hashes = _boot()
    D = runtime.D
    admission = D._bound(_relative(OUT / f"{record}_input_admission.json"), admission_sha, hashes)
    require(admission["status"] == "ALL_168_EXPECTED_PACKETS_BOUND_BEFORE_LOCAL_READ"
            and admission["record_id"] == record and admission["local_component_payloads_opened"] is False,
            "successful immutable admission required before any packet reads")
    _prelaunch(D, admission["review_sha256"], hashes)
    ctx, protocol, domains, body, routes, links, stage = _load_record(runtime, plan, jobs, expect, record, hashes)
    for n, h in admission["source_and_manifest_hashes"].items():
        D._bound(n, h, hashes, parse=False)
    require(all(_terminal_part(D, stage, links, record, part, admission["external_part_pins"][part], hashes)
                == admission["expected_packets"][part] for part in PARTS), "completed part expectation changed")
    day = lambda: datetime.now(timezone.utc).isoformat()
    start = monotonic()
    prefix = OUT / record
    prefix.mkdir(exist_ok=False)
    D._write(prefix / "started.json", {"record_id": record, "pid": os.getpid(),
                                       "admission_sha256": admission_sha, "utc": day()})
    try:
        result = _reduce(D, ctx, protocol, domains, body, routes, admission["expected_packets"], hashes)
        result.update({"admission_sha256": admission_sha, "input_hashes": hashes, "seconds": monotonic() - start})
        h = D._write(prefix / "results.json", result)
        D._write(prefix / "terminal.json", {"status": "FINISHED", "exit_code": 0, "pid": os.getpid(),
                                           "results_sha256": h, "utc": day()})
        print(json.dumps({"record_id": record, "status": result["status"], "results_sha256": h}), flush=True)
    except BaseException as error:
        D._write(prefix / "failure.json", {"record_id": record, "message": str(error)})
        D._write(prefix / "terminal.json", {"status": "FAILED", "exit_code": 1, "pid": os.getpid(), "utc": day()})
        raise


def tests(D):
    checks = {}

    def reject(name, action):
        try:
            action()
        except ValueError:
            checks[name] = True
        else:
            raise AssertionError(name)

    links = {"record_id": RECORDS[0], "record_proofs_sha256": "b" * 64}
    files = [{"face": f, "edge": e, "file": f"local/face{f}_edge{e:02d}.json", "sha256": "a" * 64}
             for f in (0, 1) for e in range(12)]
    doc = {**links, "part": "surface", "status": "COMPONENT_PART_COMPLETE_UNREVIEWED", "files": files}
    read = lambda d: _manifest(d, links, "surface", "local")
    checks["complete_exact_part_manifest_only"] = read(doc) == files
    reject("missing_manifest_item", lambda: read({**doc, "files": files[:-1]}))
    reject("duplicate_packet_in_manifest", lambda: read({**doc, "files": files[:-1] + files[:1]}))
    reject("wrong_record_manifest", lambda: read({**doc, "record_id": RECORDS[1]}))
    reject("mixed_proof_manifest", lambda: read({**doc, "record_proofs_sha256": "c" * 64}))
    reject("completed_record_rejected", lambda: _record_dir("optoGEF2_small_col_L"))
    reject("failed_terminal", lambda: D._terminal({"status": "FAILED", "pid": 1, "exit_code": 1,
                                                 "results_sha256": "0" * 64}))
    pack = lambda a, b: {"lower_exact": str(a), "upper_exact": str(b)}
    checks["exact_tail_obeys_budget_even_if_packed_rounding_outward_upper_is_larger"] = (
        _exact_tail(D, "1/3", Q(1, 3), pack(Q(5, 16), Q(3, 8))) == Q(1, 3))
    reject("negative_tail_rejected", lambda: _exact_tail(D, "-1/3", Q(1), pack(-1, 0)))
    reject("overbudget_exact_tail_rejected", lambda: _exact_tail(D, "1/3", Q(1, 4), pack(0, 1)))
    reject("noncontaining_tail_rejected", lambda: _exact_tail(D, "1/3", Q(1), pack(0, Q(1, 4))))
    force, H = [D.ZERO] * 2, {}
    hess = [{"i": i, "j": j, "enclosure": pack(-3, -1)} for i, j in combinations_with_replacement(range(2), 2)]
    D._scatter(force, H, [1, 0], [pack(-3, -1), pack(-2, 4)], hess)
    rf, rows, _ = D._dyadic_reduce([[0.5, -1.0], [0., -0.5]], force, D._h_rows(2, H))
    checks["negative_force_and_negative_Hessian_preserved_until_preconditioning"] = (
        rf == [(Q(0), Q(5)), (Q(1, 2), Q(3, 2))] and rows == [Q(9), Q(3)])
    reject("missing_Hessian_coordinate", lambda: D._scatter([D.ZERO] * 2, {}, [1, 0],
                                                            [pack(0, 0)] * 2, hess[:-1]))
    checks["strict_root_gate_equality_unresolved"] = D._strict_test(Q(1, 4), Q(1, 4), Q(0), Q(1, 4), Q(1, 2))[-2:] == (False, None)
    checks["strict_pass_root_radius_bound_exact"] = D._strict_test(Q(1, 8), Q(1, 4), Q(1, 8), Q(1, 8), Q(1))[-2:] == (True, Q(2, 5))
    require(_yes(checks), "unsigned/substituted or incomplete test fixture unexpectedly accepted")
    return checks


def prepare():
    runtime, plan, jobs, expect, hashes = _boot()
    D = runtime.D
    no = AssertionError("preparation attempted actual packet reduction or cached numerical work")
    with (
        patch.object(D, "_scatter", side_effect=no),
        patch.object(D, "_dyadic_reduce", side_effect=no),
        patch.object(runtime.E.core, "assemble_point", side_effect=no),
        patch.object(runtime.E.core, "local_terms", side_effect=no),
        patch.object(runtime.E.U, "uniform_cache", side_effect=no),
    ):
        for record in RECORDS:
            ctx, p, _, _, _, _, _ = _load_record(runtime, plan, jobs, expect, record, hashes)
            D._old_premises(ctx, p, hashes)
    checks = tests(D)
    checks.update({"all15_same_record_domain_and_strict_old_radius_provider_bound": True,
                   "no_actual_component_files_opened_or_global_sums": True})
    OUT.mkdir(parents=True, exist_ok=True)
    names = ("exact_cohort_record_reduction_frozen.py", "protocol.json", "preflight.json")
    require(not any((OUT / n).exists() for n in names), "preserve write-once direct reduction preparation")
    source = Path(__file__).resolve()
    with (OUT / names[0]).open("xb") as f:
        f.write(source.read_bytes())
    psha = D._write(OUT / names[1], {
        "status": "FROZEN_REMAINING15_DIRECT_REDUCTION_PREPARATION_REVIEW_PENDING",
        "source_file": _relative(source), "source_sha256": _sha(source), "input_hashes": hashes, "records": RECORDS,
        "binding": "For a single explicit record, pin both expected worker terminals and manifests, plus their observed "
                   "OS exit digests and matching supervisor folders. Neither the other records' progress nor a partial "
                   "packet set authorizes reading any component payload.",
        "math": "Signed exact interval scatter for all 144+24 components, recorded 611F and 6676 upper H entries, "
                "integer-dyadic R*F and |R|*absolute Hessian rows. Direct etaQ+d and zQ+b on same original r/cell. "
                "Only a strict pass can conditionally use the identical-record cached SPD review.",
        "independence": "The new reducer reuses the previously independently reviewed endpoint primitives. "
                        "Actual component and actual reduction review remain separate requirements.",
        "actual_component_files_opened": False, "actual_global_sums_run": False,
        "failed_attempts": {_relative(p): _sha(p) for p in sorted(OUT.glob("development_revision_*/*")) if p.is_file()},
    })
    fsha = D._write(OUT / names[2], {
        "status": "PASS_GENERIC_DIRECT_REDUCTION_PREPARATION_ONLY", "source_sha256": _sha(source),
        "protocol_sha256": psha, "checks": checks, "records": len(RECORDS),
        "actual_component_files_opened": False, "actual_global_sums_run": False,
    })
    print(json.dumps({"source_sha256": _sha(source), "protocol_sha256": psha, "preflight_sha256": fsha,
                      "checks": len(checks), "records": len(RECORDS)}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest="command", required=True)
    command.add_parser("prepare")
    admission = command.add_parser("admit")
    admission.add_argument("record")
    admission.add_argument("--review-sha256", required=True)
    for part in PARTS:
        for flag in ("terminal-sha256", "manifest-sha256", "os-terminal-sha256", "supervisor-dir"):
            admission.add_argument(f"--{part}-{flag}", required=True)
    computation = command.add_parser("reduce")
    computation.add_argument("record")
    computation.add_argument("--admission-sha256", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare()
    elif args.command == "admit":
        keys = ("terminal_sha256", "manifest_sha256", "os_terminal_sha256", "supervisor_dir")
        admit(args.record, args.review_sha256, {part: {k: getattr(args, f"{part}_{k}") for k in keys} for part in PARTS})
    else:
        reduce(args.record, args.admission_sha256)
