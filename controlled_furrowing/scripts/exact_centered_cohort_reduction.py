"""Record-neutral exact reduction of completed independently replayed centered force."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
import importlib
import importlib.util
import json
import os
from pathlib import Path, PurePosixPath
from time import monotonic
from types import SimpleNamespace
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
BASE = "runs/exact_integration_cohort_v1"
OUT = ROOT / BASE / "remaining15_centered_reduction_preparation"
CALLER = BASE + "/remaining15_centered_actual_preparation"
HELPER = BASE + "/remaining15_centered_preparation"
ARITHMETIC = BASE + "/L_defect_diagnostic/centered_bulk_L_v1/reduction_preparation/reduce_centered_L.py"
ADMITTED_WORKER_LINKS = "expected_worker_links"
ADMISSION_VERSION = 2
PIN = {
    "scripts/remaining_cohort_centered.py": "bd384270f40fd08962d1c82c33ce5203af7c3f40c79654b8275801ffa22e397c",
    CALLER + "/protocol.json": "1a4ffdd07db4b1073efc2715eb6a548f2011da9d44e286d115713fe7f1d2b845",
    CALLER + "/preflight.json": "3b10c934c20f6c802355257b7cc710663b873c34f39201816bce8ca4ac88d77f",
    CALLER + "/manifest.json": "8cf34d5a9fded6cc6c5a761d57e157da292c315d871bf34c2ec352d707b9ce65",
    HELPER + "/independent_review.json": "ffcebc86827a4e12cdbec8e0a7fd2b5be84e884eb83f88078d9e61356e7b3769",
    ARITHMETIC: "6e35f18b3d6f51360667074b312ca094f1aa544a17c06be8548cdd9f6c7bd1cc",
}


def require(ok, why):
    if not ok:
        raise ValueError(why)


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _relative(path):
    return Path(path).relative_to(ROOT).as_posix()


def _utc():
    return datetime.now(timezone.utc).isoformat()


def _boot():
    for n, h in PIN.items():
        require(_sha(ROOT / n) == h, f"changed prerequisite: {n}")
    W = importlib.import_module("remaining_cohort_centered")
    require(Path(W.__file__).resolve() == ROOT / "scripts/remaining_cohort_centered.py", "centered caller alias changed")
    H, runtime, plan, jobs, expect, hashes = W._boot()
    _, op, L = runtime
    D = op.D
    spec = importlib.util.spec_from_file_location("reviewed_signed_centered_arithmetic", ROOT / ARITHMETIC)
    A = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(A)
    require(Path(A.__file__).resolve() == (ROOT / ARITHMETIC).resolve() and A.D is D and A.G is L,
            "reviewed signed arithmetic aliases changed")
    for n, h in PIN.items():
        D._bound(n, h, hashes, parse=False)
    ready = D._bound(CALLER + "/manifest.json", PIN[CALLER + "/manifest.json"], hashes)
    require(ready["status"] == "READY_FOR_CENTERED_LAUNCHER_INDEPENDENT_PRELAUNCH"
            and set(ready["records"]) == set(W.RECORDS)
            and all(len(r["elements"]) == 144 and len(r["retained_surface_manifest"]) == 24
                    and r["unchanged_full_H_upper_entries"] == 6676 for r in ready["records"].values()),
            "complete current centered caller preparation required")
    return W, H, runtime, A, plan, jobs, expect, hashes, ready


def _declared_record(D, ready, record, hashes):
    require(record in ready["records"], "explicit currently reviewed centered record required")
    reference = ready["records"][record]
    helper = D._bound(reference["helper_proof_file"], reference["helper_proof_sha256"], hashes)
    require(reference["record_id"] == helper["record_id"] == record
            and [x["element"] for x in reference["elements"]] == [x["element"] for x in helper["bulk_proofs"]] == list(range(144))
            and reference["retained_surface_manifest"] == helper["retained_surface_manifest"]
            and all(a["record_id"] == b["record_id"] == record and a["centered_proof_sha256"] == D._digest(b)
                    for a, b in zip(reference["elements"], helper["bulk_proofs"])),
            "one complete record-scoped centered proof manifest required")
    return reference, helper


def _own_gate(D, review_sha, hashes):
    gate_path = _relative(OUT / "prelaunch_independent_review.json")
    gate = D._bound(gate_path, review_sha, hashes)
    identity = {"source_sha256": _sha(__file__), "protocol_sha256": _sha(OUT / "protocol.json"),
                "preflight_sha256": _sha(OUT / "preflight.json")}
    require(gate["status"] == "PASS" and D._yes(gate["checks"]) and all(gate[k] == v for k, v in identity.items()),
            "current independent centered composition prelaunch required")
    for n, h in ((_relative(__file__), identity["source_sha256"]),
                 (_relative(OUT / "protocol.json"), identity["protocol_sha256"]),
                 (_relative(OUT / "preflight.json"), identity["preflight_sha256"])):
        D._bound(n, h, hashes, parse=False)
    return identity


def _worker_gate(D, expected, record, hashes):
    gate = D._bound(CALLER + "/prelaunch_independent_review.json", expected, hashes)
    identity = {"source_sha256": PIN["scripts/remaining_cohort_centered.py"], "protocol_sha256": PIN[CALLER + "/protocol.json"],
                "preflight_sha256": PIN[CALLER + "/preflight.json"], "manifest_sha256": PIN[CALLER + "/manifest.json"]}
    require(gate["status"] == "PASS" and D._yes(gate["checks"])
            and all(gate[k] == v for k, v in identity.items())
            and gate["helper_review_sha256"] == PIN[HELPER + "/independent_review.json"] and record in gate["records"],
            "matching independent centered caller launch permission required")
    return {**identity, "prelaunch_sha256": expected}


def _complete_files(D, value, folder, links):
    require(value["status"] == "COMPLETE_144_CENTERED_FORCE_PACKETS_UNREVIEWED"
            and value["new_Hessian_integrations"] is False and value["root_asserted"] is False
            and all(value.get(k) == v for k, v in links.items()), "mixed centered terminal manifest")
    nodes = value["files"]
    require(len(nodes) == 144 and all(type(x["element"]) is int and x["element"] == i
                and x["file"] == f"{folder}/components/element_{i:03d}.json" for i, x in enumerate(nodes)),
            "strict canonical complete 144-element centered coverage required")
    for x in nodes:
        D._sha(x["sha256"])
    return nodes


def _terminal_inputs(D, record, reference, pins, hashes):
    launch = _worker_gate(D, pins["worker_review_sha256"], record, hashes)
    folder = f"{CALLER}/actual/{record}"
    terminal = D._bound(folder + "/terminal.json", pins["terminal_sha256"], hashes)
    D._terminal(terminal, pins["manifest_sha256"])
    supervisor = PurePosixPath(pins["supervisor_dir"])
    require(supervisor.parent.as_posix() == CALLER + "/launches" and D._sha(supervisor.name) == supervisor.name,
            "exact authorized centered OS supervisor folder required")
    os_exit = D._bound(f"{supervisor}/{record}.os_terminal.json", pins["os_terminal_sha256"], hashes)
    child_path, super_path = f"{supervisor}/{record}.started.json", f"{supervisor}/started.json"
    child = D._bound(child_path, _sha(ROOT / child_path), hashes)
    history = D._bound(super_path, _sha(ROOT / super_path), hashes)
    supervisor_digest = sha256(json.dumps([launch["manifest_sha256"], launch["prelaunch_sha256"],
                                           history["explicit_records"]]).encode()).hexdigest()
    require(supervisor.name == supervisor_digest and history["manifest_sha256"] == launch["manifest_sha256"]
            and history["review_sha256"] == launch["prelaunch_sha256"]
            and record in history["explicit_records"] and child["command"][-1] == record
            and os_exit["pid"] == child["pid"] == terminal["pid"] and type(os_exit["exit_code"]) is int
            and os_exit["exit_code"] == 0, "observed successful OS/worker terminal pair required")
    current = D._bound(folder + "/results.json", pins["manifest_sha256"], hashes)
    record_links = {**launch, "record_id": record, "factory_nodeball_sha256": reference["factory_nodeball_sha256"],
                    "operator_review_sha256": reference["operator_review_sha256"]}
    files = _complete_files(D, current, folder, record_links)
    require(current["retained_surface_manifest"] == reference["retained_surface_manifest"]
            and current["retained_full_H_upper_entries"] == reference["unchanged_full_H_upper_entries"] == 6676
            and current["retained_fullbox_b_exact"] == reference["unchanged_fullbox_b_exact"],
            "unchanged full-box H and retained surface binding required")
    return files, record_links


def admit(record, review_sha, pins):
    _, _, runtime, _, _, _, _, hashes, ready = _boot()
    D = runtime[1].D
    _own_gate(D, review_sha, hashes)
    reference, _ = _declared_record(D, ready, record, hashes)
    files, links = _terminal_inputs(D, record, reference, pins, hashes)
    h = D._write(OUT / f"{record}_input_admission.json", {
        "status": "BOUND_144_CENTERED_PACKETS_BEFORE_LOCAL_READ", "record_id": record,
        "admission_schema_version": ADMISSION_VERSION,
        "review_sha256": review_sha, "external_terminal_pins": pins, ADMITTED_WORKER_LINKS: links,
        "expected_packets": files, "input_hashes": hashes, "local_centered_payloads_opened": False,
    })
    print(json.dumps({"record_id": record, "admission_sha256": h, "expected_centered_packets": len(files)}), flush=True)


def _packet_review(D, record, expected_sha, admission_sha, pins, identity, hashes):
    result = D._bound(_relative(OUT / "independent_actual" / f"{record}_packet_review.json"), expected_sha, hashes)
    require(result["status"] == "PASS_CENTERED_PACKET_REPLAY" and D._yes(result["checks"])
            and result["record_id"] == record and result["components"] == 144
            and result["admission_sha256"] == admission_sha
            and result["centered_manifest_sha256"] == pins["manifest_sha256"]
            and result["centered_terminal_sha256"] == pins["terminal_sha256"]
            and all(result[k] == v for k, v in identity.items()),
            "independent actual replay of this admitted complete centered packet set required")


def _packet_metadata(D, packet, proof, receipt, links, record):
    require(packet["status"] == "CENTERED_RECORD_BULK_FORCE_UNREVIEWED" and packet["record_id"] == proof["record_id"] == record
            and packet["new_Hessian_integrations"] is False and packet["precision_bits"] == 640
            and packet["verified_local_nodeball_sha256"] == receipt["nodeball_sha256"]
            and packet["centered_proof_sha256"] == receipt["centered_proof_sha256"] == D._digest(proof)
            and all(packet.get(k) == v for k, v in {**proof, **links}.items()),
            "mixed actual centered packet, exact record proof or local uniform source")


def _bounds(D, packed, n):
    require(len(packed) == n, "complete signed local centered coordinates required")
    return list(map(D._packed, packed))


def _packet_force(D, A, packet, proof, receipt, links, record, uniform_H, uniform_force):
    _packet_metadata(D, packet, proof, receipt, links, record)
    n = len(proof["global_indices"])
    point = packet["point_force"]
    f0 = _bounds(D, point["normalized_force_total"], n)
    volume, shear, cache = [_bounds(D, point[k], n) for k in (
        "volumetric_truncation", "quadratic_strain_defect", "finite_same_Y_cache_correction")]
    tail = list(map(Q, point["two_mass_tail_upper_exact"]))
    require(len(tail) == n and min(tail) >= 0 and max(tail) == Q(proof["point_force_tail_max_exact"])
            and max(tail) <= Q(proof["force_tail_budget_exact"])
            and D._packed(point["cache_min_J_bound"])[0] > Q(1, 10**6), "same exact point tail and analytic cache required")
    gain = _bounds(D, packet["affine_gain_derivative_at_Y"], n)
    offset = tuple(map(Q, proof["common_gain_shift_exact"]))
    require(len(offset) == 2 and offset[0] <= 0 <= offset[1]
            and (proof["active"] != "0" or all(x == D.ZERO for x in gain)), "same gain path and inactive exact zero required")
    state, parameter, centered = [_bounds(D, packet[k], n) for k in (
        "uniform_H_times_nodal_displacement", "affine_gain_increment", "centered_force_enclosure")]
    intersections = [tuple(map(Q, p)) for p in packet["force_intersection_with_prior_exact"]]
    radii = list(map(Q, proof["original_local_radii_exact"]))
    require(len(intersections) == len(uniform_force) == len(radii) == n and min(radii) >= 0, "same local complete retained bundle required")
    for i in range(n):
        A._intersect(A._sum(volume[i], shear[i], cache[i], (-tail[i], tail[i])), f0[i])
        u = sum((D._abs_upper(uniform_H[min(i, j), max(i, j)]) * r for j, r in enumerate(radii)), Q(0))
        require(D._contains(state[i], (-u, u)) and D._contains(parameter[i], A._multiply(gain[i], offset)),
                "centered nodal and gain increments differ from exact replay")
        A._intersect(centered[i], A._sum(f0[i], state[i], parameter[i]))
        require(intersections[i] == A._intersect(centered[i], uniform_force[i]), "signed local force intersection changed")
    return intersections, f0, gain


def _calculate(W, H, runtime, A, plan, jobs, expect, record, reference, helper, admission, hashes):
    R, op, _ = runtime
    D = op.D
    domain, actual, routes, packets, canonical, sig = W._factory(H, runtime, plan, jobs, expect, record, helper, hashes)
    require(sig == reference["factory_nodeball_sha256"], "same authenticated original domain mandatory")
    ctx, p = domain.context, domain.protocol
    f = [D.ZERO] * 611
    Hbulk, Htotal = {}, {}
    point, gain, surface = ([D.ZERO] * 611 for _ in range(3))
    counts, hcounts = [0] * 611, [0] * 611
    for i, bound in enumerate(admission["expected_packets"]):
        receipt, uniform, _, _, proof = W._element_guard(
            H, runtime, domain, canonical, actual, routes, packets, helper, i, hashes
        )
        require(receipt == reference["elements"][i] and bound["element"] == i, "changed or crossed launcher element proof")
        payload = D._bound(bound["file"], bound["sha256"], hashes)
        ids = proof["global_indices"]
        rows = uniform["H_integral_minus_cache_upper"]
        local_H = {k: D._packed(v["enclosure"]) for k, v in D._entry_map(rows, len(ids)).items()}
        force, anchor, gain_rows = _packet_force(D, A, payload, proof, receipt, admission[ADMITTED_WORKER_LINKS],
                    record, local_H, list(map(D._packed, uniform["force_integral_minus_cache"])))
        A._force_scatter(point, ids, anchor)
        A._force_scatter(gain, ids, gain_rows)
        D._scatter(f, Hbulk, ids, [{"lower_exact": str(a), "upper_exact": str(b)} for a, b in force], rows)
        for j in ids:
            counts[j] += 1
            hcounts[j] += len(ids)
    Htotal.update(Hbulk)
    for i, bound in enumerate(packets["surface"]):
        old = D._bound(bound["file"], bound["sha256"], hashes)
        case = actual["surface_cases"][i]
        ids = D._packet_meta(ctx, "surface", i, case, old, domain.body, p)
        force = old["force_integral_minus_cache"]
        A._force_scatter(surface, ids, list(map(D._packed, force)))
        D._scatter(f, Htotal, ids, force, old["H_integral_minus_cache_upper"])
        for j in ids:
            counts[j] += 1
            hcounts[j] += len(ids)
    require(counts == p["force_incidence"] and hcounts == p["H_row_incidence"] and len(Htotal) == 6676
            and len(packets["surface"]) == 24, "same full exact 168 finite-element incidences required")
    audit = D._bound(f"{H.DIRECT}/independent_actual/{record}.json", reference["operator_review_sha256"], hashes)
    prior = D._bound(f"{H.DIRECT}/{record}/results.json", audit["producer_pins"]["results_sha256"], hashes)
    entries = [[i, j, *Htotal[i, j]] for i, j in sorted(Htotal)]
    require(D._digest(entries) == D._digest(prior["global_H_upper_intervals_exact"]),
            "retained full-box H changed")
    Rmat = ctx.R
    RF, br, norm = D._dyadic_reduce(Rmat, f, D._h_rows(611, Htotal))
    b = max(br)
    require(b == Q(reference["unchanged_fullbox_b_exact"]) == Q(prior["b_upper_exact"])
            and norm == Q(prior["R_infinity_norm_exact"]), "same full-box Hessian defect b and preconditioner required")
    zero = [Q(0)] * 611
    rp, rg, rs = [D._dyadic_reduce(Rmat, g, zero)[0] for g in (point, gain, surface)]
    radii = [Q(domain.box["coordinate_radii_exact"][int(j) // 2][int(j) % 2]) for j in ctx.free]
    bulk_state, RHS = A._bulk_RH_effects(Rmat, Hbulk, radii, ctx.S)
    order = op.C.ORDER
    offsets = [(ctx.cell[k][0] - ctx.base[k], ctx.cell[k][1] - ctx.base[k]) for k in order]
    gain_col = order.index("common_gain")
    gain_term = [A._multiply(row, offsets[gain_col]) for row in rg]
    separated = [A._sum(x, (-y, y), z, w) for x, y, z, w in zip(rp, bulk_state, gain_term, rs)]
    joint = [A._intersect(x, y) for x, y in zip(RF, separated)]
    predictor = A._coupled_predictor_force(rp, RHS, rg, rs, offsets, gain_col)
    final = [A._intersect(a, z) for a, z in zip(joint, predictor)]
    eta0, z0, r, provider = D._old_premises(ctx, p, hashes)
    defect = max(map(D._abs_upper, final))
    eta, z, lhs, ok, error = D._strict_test(eta0, z0, defect, b, r)
    curvature = D._minimum(ctx, hashes, z0, z, r, norm, provider) if ok else None
    return {
        "status": "CONDITIONAL_UNREVIEWED_CENTERED_SUFFICIENT_TEST_PASSES" if ok
                  else "UNREVIEWED_CENTERED_SUFFICIENT_TEST_UNRESOLVED",
        "record_id": record, "same_closed_C_sha256": ctx.row["parameter_cell_sha256"], "same_exact_source_dose": ctx.dose,
        "bulk_components": 144, "retained_surface_components": 24, "F_entries": 611, "upper_H_entries": len(Htotal),
        "global_F_intervals_exact": f, "retained_global_H_upper_intervals_exact": entries,
        "direct_whole_box_RF_intervals_exact": RF, "separated_whole_box_RF_intervals_exact": separated,
        "whole_box_RF_intersections_exact": joint, "coupled_predictor_RF_only_intervals_exact": predictor,
        "final_predictor_RF_only_intersections_exact": final,
        "point_bulk_RF_intervals_exact": rp, "RH_bulk_S_signed_intervals_exact": RHS,
        "affine_gain_RF_intervals_exact": rg, "retained_surface_RF_intervals_exact": rs,
        "node_path_bulk_only_RH_radii_exact": bulk_state, "parameter_order": order, "closed_five_offsets_exact": offsets,
        "exact_R_infinity_norm": norm, "same_fullbox_b_exact": b, "old_cached_eta_exact": eta0,
        "old_cached_Z_exact": z0, "same_radius_exact": r, "RF_predictor_defect_upper_exact": defect,
        "RF_predictor_defect_lower_exact": max(map(A._active_bound, final)),
        "total_eta_upper_exact": eta, "total_Z_upper_exact": z, "strict_selfmap_lhs_exact": lhs,
        "sufficient_conditions_passed": ok, "conditional_root_error_exact": error, "conditional_curvature": curvature,
        "independent_root_assertion": False,
        "scope": "All 168 components, original closed C and tube. The final force bound is predictor-only. The derivative "
                 "bound remains the full-box Hessian including original surfaces. Actual global replay still needed.",
    }


def reduce(record, admission_sha, packet_review_sha):
    W, H, runtime, A, plan, jobs, expect, hashes, ready = _boot()
    D = runtime[1].D
    prep = D._bound(_relative(OUT / f"{record}_input_admission.json"), admission_sha, hashes)
    require(prep["status"] == "BOUND_144_CENTERED_PACKETS_BEFORE_LOCAL_READ" and prep["record_id"] == record
            and prep["admission_schema_version"] == ADMISSION_VERSION
            and prep["local_centered_payloads_opened"] is False, "same immutable completed centered admission required")
    identity = _own_gate(D, prep["review_sha256"], hashes)
    reference, helper = _declared_record(D, ready, record, hashes)
    for path, digest in prep["input_hashes"].items():
        D._bound(path, digest, hashes, parse=False)
    files, links = _terminal_inputs(D, record, reference, prep["external_terminal_pins"], hashes)
    require(files == prep["expected_packets"] and links == prep[ADMITTED_WORKER_LINKS], "changed admitted centered manifest")
    _packet_review(D, record, packet_review_sha, admission_sha, prep["external_terminal_pins"], identity, hashes)
    folder = OUT / record
    folder.mkdir(exist_ok=False)
    start = monotonic()
    D._write(folder / "started.json", {"pid": os.getpid(), "utc": _utc(), "admission_sha256": admission_sha,
                                       "actual_packet_review_sha256": packet_review_sha})
    try:
        result = _calculate(W, H, runtime, A, plan, jobs, expect, record, reference, helper, prep, hashes)
        h = D._write(folder / "results.json", {**result, **identity, "admission_sha256": admission_sha,
                                             "actual_packet_review_sha256": packet_review_sha, "input_hashes": hashes,
                                             "seconds": monotonic() - start})
        D._write(folder / "terminal.json", {"pid": os.getpid(), "utc": _utc(), "status": "FINISHED", "exit_code": 0,
                                            "results_sha256": h})
        print(json.dumps({"record_id": record, "status": result["status"], "results_sha256": h}), flush=True)
    except BaseException as error:
        D._write(folder / "failure.json", {"message": str(error)})
        D._write(folder / "terminal.json", {"pid": os.getpid(), "utc": _utc(), "status": "FAILED", "exit_code": 1})
        raise


def _mocked_admit_reduce_roundtrip(D):
    """Run the real admission and reduction entry points over in-memory artifact I/O."""
    record = "mocked_roundtrip_record"
    stage = OUT / "synthetic_mocked_roundtrip"
    admission_file = _relative(stage / f"{record}_input_admission.json")
    review_file = _relative(stage / "independent_actual" / f"{record}_packet_review.json")
    result_file = _relative(stage / record / "results.json")
    terminal_file = _relative(stage / record / "terminal.json")
    identity = {"source_sha256": "1" * 64, "protocol_sha256": "2" * 64, "preflight_sha256": "3" * 64}
    pins = {"manifest_sha256": "4" * 64, "terminal_sha256": "5" * 64}
    worker_links = {"record_id": record, "operator_review_sha256": "6" * 64}
    files = [{"element": i, "file": f"mocked_components/element_{i:03d}.json", "sha256": "7" * 64}
             for i in range(144)]
    documents = {}

    def write(path, value):
        name = _relative(path)
        require(name not in documents, "mocked artifacts must be write-once")
        documents[name] = value
        return D._digest(value)

    def read(path, expected, bindings, parse=True):
        value = documents[path]
        require(D._digest(value) == expected, "mocked artifact digest changed")
        bindings[path] = expected
        return value if parse else None

    def calculate(_, __, ___, ____, _____, ______, _______, actual_record, reference, helper, admission, hashes):
        require(actual_record == record and admission["expected_packets"] == files
                and admission["expected_worker_links"] == worker_links,
                "real admission document does not reach the actual reducer schema")
        require(review_file in hashes, "independent actual packet gate must run before any calculation")
        return {"status": "SYNTHETIC_CALCULATION_ONLY", "record_id": actual_record}

    def mkdir(self, **_):
        require(self == stage / record, "unexpected filesystem path in mocked roundtrip")

    runtime = (None, SimpleNamespace(D=D), None)
    replacements = {
        "OUT": stage,
        "_boot": lambda: (None, None, runtime, None, None, None, None, {}, {}),
        "_own_gate": lambda *_: identity,
        "_declared_record": lambda *_: ({}, {}),
        "_terminal_inputs": lambda *_: (files, worker_links),
        "_calculate": calculate,
    }
    with (
        patch.dict(globals(), replacements),
        patch.object(D, "_bound", side_effect=read),
        patch.object(D, "_write", side_effect=write),
        patch.object(Path, "mkdir", mkdir),
    ):
        admit(record, "8" * 64, pins)
        admitted = documents[admission_file]
        require(admitted["expected_worker_links"] == worker_links and admitted["admission_schema_version"] == 2
                and admitted["local_centered_payloads_opened"] is False, "new admission schema was not serialized")
        h = D._digest(admitted)
        documents[review_file] = {
            **identity, "status": "PASS_CENTERED_PACKET_REPLAY", "checks": {"mocked": True}, "record_id": record,
            "components": 144, "admission_sha256": h, "centered_manifest_sha256": pins["manifest_sha256"],
            "centered_terminal_sha256": pins["terminal_sha256"],
        }
        review_sha = D._digest(documents[review_file])
        reduce(record, h, review_sha)
        require(documents[result_file]["admission_sha256"] == h
                and documents[result_file]["actual_packet_review_sha256"] == review_sha
                and documents[terminal_file]["results_sha256"] == D._digest(documents[result_file])
                and documents[terminal_file]["exit_code"] == 0, "mocked end-to-end terminal binding failed")
        documents[admission_file] = {**admitted, "expected_worker_links": {**worker_links, "record_id": "other"}}
        try:
            reduce(record, D._digest(documents[admission_file]), review_sha)
        except ValueError as error:
            require(str(error) == "changed admitted centered manifest", "wrong mixed-link rejection")
        else:
            raise AssertionError("replayed admission accepted mixed worker links")
    return True


def tests(D, A):
    record = "synthetic"
    paths = [{"element": i, "file": f"toy/components/element_{i:03d}.json", "sha256": "a" * 64} for i in range(144)]
    doc = {"status": "COMPLETE_144_CENTERED_FORCE_PACKETS_UNREVIEWED", "record_id": record, "files": paths,
           "new_Hessian_integrations": False, "root_asserted": False}
    links = {"record_id": record}
    results = {"canonical_complete_expected_set": _complete_files(D, doc, "toy", links) == paths}

    def reject(name, call):
        try:
            call()
        except ValueError:
            results[name] = True
        else:
            raise AssertionError(name)

    reject("partial_manifest", lambda: _complete_files(D, {**doc, "files": paths[:-1]}, "toy", links))
    reject("duplicate_manifest", lambda: _complete_files(D, {**doc, "files": paths[:-1] + paths[:1]}, "toy", links))
    reject("mixed_record", lambda: _complete_files(D, {**doc, "record_id": "other"}, "toy", links))
    reject("root_claim_in_worker", lambda: _complete_files(D, {**doc, "root_asserted": True}, "toy", links))
    reject("nonterminal_inputs", lambda: D._terminal({"status": "FAILED", "pid": 5, "exit_code": 1, "results_sha256": "a" * 64}))
    results["strict_gate_equality_is_unresolved"] = (
        D._strict_test(Q(1, 4), Q(1, 4), Q(0), Q(1, 4), Q(1, 2))[-2:] == (False, None))
    proof = {"record_id": record, "uniform_packet_sha256": "a" * 64}
    receipt = {"nodeball_sha256": "b" * 64, "centered_proof_sha256": D._digest(proof)}
    packet = {**proof, "status": "CENTERED_RECORD_BULK_FORCE_UNREVIEWED", "new_Hessian_integrations": False,
              "precision_bits": 640, "verified_local_nodeball_sha256": "b" * 64, "centered_proof_sha256": D._digest(proof)}
    _packet_metadata(D, packet, proof, receipt, links, record)
    results["same_exact_centered_proof_accepted"] = True
    reject("changed_uniform_packet", lambda: _packet_metadata(D, {**packet, "uniform_packet_sha256": "c" * 64},
                                                            proof, receipt, links, record))
    reject("changed_realized_nodeball", lambda: _packet_metadata(D, {**packet, "verified_local_nodeball_sha256": "c" * 64},
                                                              proof, receipt, links, record))
    results["signed_RH_bulk_only_unchanged"] = A._bulk_RH_effects(
        [[1.0, 1.0], [-1.0, 1.0]], {(0, 0): (Q(1), Q(1)), (0, 1): (Q(-1), Q(-1)), (1, 1): (Q(1), Q(1))},
        [Q(1, 3), Q(1, 7)], [[1., 0., 0., 0., 1.], [0., 0., 0., 0., -1.]])[0] == [Q(0), Q(20, 21)]
    results["real_admit_then_reduce_roundtrip_with_mocked_artifacts_and_mixed_link_rejection"] = (
        _mocked_admit_reduce_roundtrip(D))
    require(D._yes(results), "all centered composition guard tests required")
    return results


def prepare():
    _, H, runtime, A, _, _, _, hashes, ready = _boot()
    D = runtime[1].D
    no = AssertionError("unexpected actual centered data, global reduction or integration during preparation")
    with (patch.object(D, "_scatter", side_effect=no), patch.object(D, "_dyadic_reduce", side_effect=no),
          patch.object(H, "bulk_component", side_effect=no), patch.object(H, "bind_reviewed_record", side_effect=no)):
        for record in ready["records"]:
            _declared_record(D, ready, record, hashes)
    checks = tests(D, A)
    checks.update({"same_current_explicit_C_D_proof_inputs": True, "no_actual_packets_or_global_real_sums": True})
    OUT.mkdir(parents=True, exist_ok=True)
    names = ("exact_centered_cohort_reduction_frozen.py", "protocol.json", "preflight.json")
    require(not any((OUT / n).exists() for n in names), "preserve write-once centered reduction preparation")
    with (OUT / names[0]).open("xb") as output:
        output.write(Path(__file__).read_bytes())
    psha = D._write(OUT / names[1], {
        "status": "FROZEN_CENTERED_COHORT_REDUCTION_PREPARATION_REVIEW_PENDING",
        "source_file": _relative(__file__), "source_sha256": _sha(__file__), "input_hashes": hashes,
        "admission_schema_version": ADMISSION_VERSION, "admitted_worker_links_field": ADMITTED_WORKER_LINKS,
        "currently_eligible_fully_reviewed_records": list(ready["records"]),
        "gates": "Independent reduction prelaunch, externally pinned complete worker and observed OS terminal, then "
                 "independent actual centered packet replay binding admission/terminal/manifest before any actual sum.",
        "arithmetic": "All 144 centered signed intersections with original reviewed force, original 24 surface F/H, exact "
                      "signed R*F, unchanged full absolute H b, separated point-bulk RF, signed RH_bulk*S, shared "
                      "preconditioned gain column and retained surface RF. Five closed increments enter once at each predictor.",
        "scope": "Only actual independently reviewed centered record packets of the same original model/cell/radius. "
                 "Global replay remains independently reviewable. No new actual integrations or point/packet work in preparation.",
        "actual_component_payloads_opened_or_global_sums": False,
        "failed_attempts": {_relative(p): _sha(p) for p in sorted(OUT.glob("development_revision_*/*")) if p.is_file()},
    })
    fsha = D._write(OUT / names[2], {"status": "PASS_GENERIC_CENTERED_REDUCTION_PREPARATION_ONLY",
                                    "source_sha256": _sha(__file__), "protocol_sha256": psha, "checks": checks,
                                    "actual_component_payloads_opened_or_global_sums": False})
    print(json.dumps({"source_sha256": _sha(__file__), "protocol_sha256": psha, "preflight_sha256": fsha,
                      "checks": len(checks)}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("prepare")
    ar = commands.add_parser("admit")
    ar.add_argument("record")
    ar.add_argument("--review-sha256", required=True)
    PIN_KEYS = ("worker-review-sha256", "terminal-sha256", "manifest-sha256", "os-terminal-sha256", "supervisor-dir")
    for key in PIN_KEYS:
        ar.add_argument("--" + key, required=True)
    run = commands.add_parser("reduce")
    run.add_argument("record")
    run.add_argument("--admission-sha256", required=True)
    run.add_argument("--packet-review-sha256", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare()
    elif args.command == "admit":
        admit(args.record, args.review_sha256, {k.replace("-", "_"): getattr(args, k.replace("-", "_")) for k in PIN_KEYS})
    else:
        reduce(args.record, args.admission_sha256, args.packet_review_sha256)
