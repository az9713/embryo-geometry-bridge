"""Independent exact endpoint reduction of a bound original-cube cohort run."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from math import ceil
from pathlib import Path, PurePosixPath

import exact_cohort_context as C


ROOT = Path(__file__).resolve().parents[1]
BASE = "runs/exact_integration_cohort_v1"
OUT = ROOT / BASE / "reduction_preparation"
ENERGY_RESULT = "runs/cohort_energy_stability_v1/results.json"
ENERGY_REVIEW = "evidence/cohort_energy_stability_independent_review.json"
PREPARATION_PINS = {
    ENERGY_RESULT: "5f590ed1bea078e2e88d17bb39f8030d144a95fb73bf96ab10b5823951883164",
    ENERGY_REVIEW: "d4f4c4d30d16bbbf2141cabf1cee18b47fcc35ab958a2a8a203000f0a7acaf95",
    BASE + "/context_preparation/protocol.json": "c1834232a5f92ef1e3d997aeba3bc6a543ef2fe4a2e1b0d20995b4a4fd9f2409",
    BASE + "/integrator_preparation/protocol.json": "cced8eb0a3b7b15d96fe0c4f43da042c807ecd109bcaf1dd37e0a0de1e979096",
}
PRIMARY_PINS = {
    "protocol": "protocol.json", "domains": "domains.json", "domains_terminal": "domains_terminal.json",
    "integration_review": "integration_prelaunch_independent_review.json",
    "bulk_terminal": "bulk_terminal.json", "surface_terminal": "surface_terminal.json",
    "assembled": "assembled.json",
}
ZERO = (Q(0), Q(0))
DYADIC_RHO_DENOM = 2**48


def require(condition, why):
    if not condition:
        raise ValueError(why)


def _yes(checks):
    return isinstance(checks, dict) and bool(checks) and all(value is True for value in checks.values())


def _sha(value):
    require(type(value) is str and len(value) == 64 and all(c in "0123456789abcdef" for c in value),
            "explicit lowercase expected SHA-256 required")
    return value


def _digest(value):
    return sha256(json.dumps(value, default=str, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _write(path, value):
    with path.open("x", encoding="utf8", newline="\n") as stream:
        json.dump(value, stream, default=str, indent=2, allow_nan=False)
        stream.write("\n")
    return sha256(path.read_bytes()).hexdigest()


def _bound(path, digest, bindings, parse=True):
    relative = PurePosixPath(path)
    require(relative.as_posix() == path and not relative.is_absolute() and ".." not in relative.parts
            and "\\" not in path and ":" not in path, "canonical study-relative artifact path required")
    source = ROOT / relative
    require(source.resolve().is_relative_to(ROOT) and bindings.get(path, digest) == _sha(digest),
            f"conflicting or out-of-study artifact: {path}")
    raw = source.read_bytes()
    require(sha256(raw).hexdigest() == digest, f"artifact changed: {path}")
    bindings[path] = digest
    return json.loads(raw) if parse else None


def _fraction(value):
    require(type(value) in (str, int, Q), "exact rational endpoint required")
    return Q(value)


def _exp(value):
    d = value.denominator
    require(d & (d - 1) == 0, "binary or dyadic endpoint required")
    return d.bit_length() - 1


def _packed(value, dyadic=True):
    lo, hi = (_fraction(value[k]) for k in ("lower_exact", "upper_exact"))
    require(lo <= hi, "ordered exact endpoints required")
    if dyadic:
        _exp(lo)
        _exp(hi)
    return lo, hi


def _contains(outer, inner):
    return outer[0] <= inner[0] <= inner[1] <= outer[1]


def _subset(outer, inner, label):
    require(_contains(outer, inner), f"producer fails endpoint containment: {label}")


def _entry_map(rows, n):
    entries = {}
    for row in rows:
        i, j = row["i"], row["j"]
        require(type(i) is type(j) is int and 0 <= i <= j < n and (i, j) not in entries,
                "unique indexed upper-triangular entries required")
        entries[i, j] = row
    return entries


def _plus(a, b):
    return a[0] + b[0], a[1] + b[1]


def _abs_upper(interval):
    return max(abs(interval[0]), abs(interval[1]))


def _scatter(g, H, ids, force, hessian):
    n = len(ids)
    require(n > 0 and len(set(ids)) == n and all(type(i) is int and 0 <= i < len(g) for i in ids),
            "unique valid free projection required")
    entries = _entry_map(hessian, n)
    require(len(force) == n and entries.keys() == set(combinations_with_replacement(range(n), 2)),
            "complete local force and upper-triangle coverage required")
    for index, value in zip(ids, force):
        g[index] = _plus(g[index], _packed(value))
    for (i, j), value in entries.items():
        a, b = ids[i], ids[j]
        key = (min(a, b), max(a, b))
        H[key] = _plus(H.get(key, ZERO), _packed(value["enclosure"]))


def _h_rows(n, H):
    sums = [Q(0)] * n
    for (i, j), value in H.items():
        absolute = _abs_upper(value)
        sums[i] += absolute
        if i != j:
            sums[j] += absolute
    return sums


def _dyadic_reduce(R, g, H_rows):
    """Signed interval Rg and positive |R| rows via scaled integers."""
    n = len(g)
    require(n > 0 and len(R) == len(H_rows) == n and all(len(row) == n for row in R),
            "matching square preconditioner, force and Hessian row dimensions required")
    r_pairs = [[float(value).as_integer_ratio() for value in row] for row in R]
    er = max(den.bit_length() - 1 for row in r_pairs for _, den in row)
    eg = max(_exp(value) for row in g for value in row)
    eh = max(map(_exp, H_rows))
    require(all(lo <= hi for lo, hi in g) and min(H_rows) >= 0, "valid endpoint/row direction required")
    dlo, dhi = zip(*[(a.numerator << (eg - _exp(a)), b.numerator << (eg - _exp(b))) for a, b in g])
    h = [value.numerator << (eh - _exp(value)) for value in H_rows]
    Rg, Rh, row_norms = [], [], []
    for row in r_pairs:
        lo = hi = abs_row = norm = 0
        for i, (num, den) in enumerate(row):
            k = num << (er - (den.bit_length() - 1))
            left, right = (dlo[i], dhi[i]) if k >= 0 else (dhi[i], dlo[i])
            lo += k * left
            hi += k * right
            abs_row += abs(k) * h[i]
            norm += abs(k)
        Rg.append((Q(lo, 2**(er + eg)), Q(hi, 2**(er + eg))))
        Rh.append(Q(abs_row, 2**(er + eh)))
        row_norms.append(Q(norm, 2**er))
    return Rg, Rh, max(row_norms)


def _strict_test(old_eta, old_z, rf, b, radius):
    require(min(old_eta, old_z, rf, b) >= 0 and 0 <= old_z < 1 and radius > 0,
            "positive radius and valid retained sufficient-bound inputs required")
    eta, z = old_eta + rf, old_z + b
    lhs = eta + z * radius
    ok = z < 1 and lhs < radius
    return eta, z, lhs, ok, eta / (1 - z) if ok else None


def _order_history(order, attempts, maximum, force_budget, H_budget):
    require(type(order) is type(maximum) is int and 2 <= order <= maximum and len(attempts) == order - 1,
            "complete bounded order history required")
    for i, a in enumerate(attempts, 2):
        force, H = (_fraction(a[key]) for key in ("force_tail_upper_exact", "H_entry_tail_upper_exact"))
        require(type(a["degree"]) is int and a["degree"] == i and min(force, H) >= 0
                and ((force <= force_budget and H <= H_budget) == (i == order)),
                "order history must choose the first order passing both budgets")


def _manifest(packet_dir, files, names):
    names = {f"{packet_dir}/{name}" for name in names}
    paths = [item["file"] for item in files]
    require(len(paths) == len(names) and len(set(paths)) == len(paths) and set(paths) == names,
            "complete, distinct canonical component manifest required")
    for item in files:
        _sha(item["sha256"])
    return files


def _terminal(terminal, results_hash=None):
    require(terminal["status"] == "FINISHED" and type(terminal["exit_code"]) is int
            and terminal["exit_code"] == 0 and type(terminal["pid"]) is int and terminal["pid"] > 0,
            "successful terminal process required")
    actual = _sha(terminal["results_sha256"])
    require(results_hash is None or actual == results_hash, "terminal result binding mismatch")
    return actual


def _prepared():
    p = json.loads((OUT / "protocol.json").read_bytes())
    pre = json.loads((OUT / "preflight.json").read_bytes())
    require(pre["status"] == "PASS_REDUCTION_PREPARATION_ONLY" and _yes(pre["checks"])
            and pre["protocol_sha256"] == sha256((OUT / "protocol.json").read_bytes()).hexdigest(),
            "matching exact reduction preparation required")
    bindings = {}
    for name, value in p["input_hashes"].items():
        _bound(name, value, bindings, parse=False)
    require(bindings[Path(__file__).relative_to(ROOT).as_posix()] == pre["source_sha256"],
            "reducer source/preflight digest mismatch")
    return p, bindings


def _record_dir(record):
    require(record.startswith("optoGEF2_small_col_") and len(record) == len("optoGEF2_small_col_L")
            and record[-1] in "BCDEFGHIJKLMNOPQR", "explicit cohort record ID required")
    return f"{BASE}/records/{record}/original_cube_v1"


def admit(record, pins):
    """Write expected packet hashes before opening a local component."""
    _, bindings = _prepared()
    folder = _record_dir(record)
    require(set(pins) == set(PRIMARY_PINS), "all externally supplied terminal and assembly pins required")
    primary = {key: _bound(f"{folder}/{name}", pins[key], bindings) for key, name in PRIMARY_PINS.items()}
    _terminal(primary["domains_terminal"], pins["domains"])
    domains, protocol, reviewed = (primary[k] for k in ("domains", "protocol", "integration_review"))
    require(domains["status"] == "PASS_FRESH_DOMAINS_AND_SELECTORS_UNREVIEWED"
            and domains["actual_force_H_integrals_run"] is False and domains["protocol_sha256"] == pins["protocol"]
            and all(document["record_id"] == record for document in (domains, protocol, reviewed))
            and reviewed["status"] == "PASS" and _yes(reviewed["checks"])
            and reviewed["domains_result_sha256"] == pins["domains"]
            and reviewed["protocol_sha256"] == pins["protocol"] and reviewed["body_sha256"] == domains["body_sha256"],
            "current independently reviewed record domain required")
    for name, digest in protocol["source_and_input_hashes"].items():
        _bound(name, digest, bindings, parse=False)
    _bound(f"{folder}/body_geometry.json", domains["body_sha256"], bindings, parse=False)
    _bound(f"{folder}/preflight.json", reviewed["preflight_sha256"], bindings, parse=False)
    _bound(reviewed["support_file"], reviewed["support_sha256"], bindings, parse=False)
    require(bindings["scripts/cohort_exact_pilot.py"] == reviewed["source_sha256"],
            "integration gate and bound producer source differ")
    packets = {}
    for part, names in (
        ("bulk", [f"element_{i:03d}.json" for i in range(144)]),
        ("surface", [f"face{face}_edge{edge:02d}.json" for face in (0, 1) for edge in range(12)]),
    ):
        result_hash = _terminal(primary[f"{part}_terminal"])
        manifest = _bound(f"{folder}/{part}_results.json", result_hash, bindings)
        require(manifest["record_id"] == record and manifest["part"] == part
                and manifest["domains_sha256"] == pins["domains"]
                and manifest["status"] == "COMPONENT_PART_COMPLETE_UNREVIEWED", "valid completed part required")
        packets[part] = _manifest(f"{folder}/{part}", manifest["files"], names)
    assembled = primary["assembled"]
    require(assembled["record_id"] == record and assembled["protocol_sha256"] == pins["protocol"]
            and assembled["domains_sha256"] == pins["domains"]
            and assembled["component_packets"] == packets["bulk"] + packets["surface"],
            "assembled packet list must exactly match the terminal manifests")
    output = OUT / f"{record}_input_admission.json"
    digest = _write(output, {"status": "BOUND_COMPONENT_INPUTS_BEFORE_LOCAL_READ", "record_id": record,
                             "run_dir": folder, "external_pins": pins, "source_and_manifest_hashes": bindings,
                             "expected_packets": packets, "local_component_payloads_opened": False})
    print(json.dumps({"admission_file": str(output), "admission_sha256": digest, "packets": 168}), flush=True)


def _old_premises(context, protocol, bindings):
    row = context.row
    root = _bound(row["root_certificate_file"], row["root_certificate_sha256"], bindings)
    review = _bound(row["independent_root_capsule_file"], row["independent_root_capsule_sha256"], bindings)
    (wrapper,) = [a for a in review["record_audits"] if a["record_id"] == context.record_id]
    audit = wrapper["audit"]
    proof = _bound(wrapper["review_file"], wrapper["review_sha256"], bindings)
    require(audit == proof and _yes(audit["checks"]), "same original cached root independent capsule required")
    if "bounds" in audit:
        source = wrapper["review_file"]
        b = audit["bounds"]
    else:
        source = audit["reuse_review_file"]
        nested = _bound(source, audit["reuse_review_sha256"], bindings)
        n = nested["independent_bounds"]
        require(nested["status"].startswith("PASS") and _yes(nested["checks"]) and _yes(nested["root_audit"]["checks"])
                and nested["root_audit"]["branch_domain_id"] == row["branch_domain_id"]
                and nested["root_audit"]["parameter_cell_sha256"] == row["parameter_cell_sha256"]
                and _fraction(nested["base_strength_exact"]) == context.drive, "matching cached M reuse required")
        b = {"eta": n["eta_uniform"], "Z": n["Z_uniform"], "radius": n["radius"]}
    eta = max(_packed(b["eta"])[1], _packed(root["uniform_cell_result"]["eta_uniform"])[1])
    z = _packed(b["Z"])[1]
    radius = context.old_radius
    require(_contains(_packed(b["radius"]), (radius, radius)) and eta >= 0 and 0 <= z < 1
            and eta + z * radius < radius
            and protocol["cached_bound_provider"] == source and _fraction(protocol["cached_eta_upper_exact"]) == eta
            and _fraction(protocol["cached_Z_upper_exact"]) == z and _fraction(protocol["normal_radius_exact"]) == radius,
            "same audited cached eta, Z and exact original radius required")
    return eta, z, radius, source


def _mesh_map(context, ids):
    index = {int(dof): pos for pos, dof in enumerate(context.free)}
    dofs = [2 * int(node) + axis for node in ids for axis in (0, 1)]
    selected = [i for i, dof in enumerate(dofs) if dof in index]
    return [index[dofs[i]] for i in selected], selected


def _tail(value, budget):
    lo, hi = _packed(value)
    require(0 <= lo <= hi <= budget, "local nonnegative tail must fit the reviewed budget")


def _rho(actual, raw):
    require(actual == Q(ceil(raw * DYADIC_RHO_DENOM), DYADIC_RHO_DENOM)
            and 0 <= raw <= actual < 1, "exact outward convergent rho required")


def _packet_meta(context, part, index, case, packet, body, protocol):
    target = context.target
    row = context.row
    labels = {
        "record_id": context.record_id, "dose_exact": str(context.dose), "precision_bits": 640,
        "parameter_cell_sha256": row["parameter_cell_sha256"], "point_cache_sha256": row["point_cache_sha256"],
        "model_observation_contract_sha256": row["model_observation_contract_sha256"],
        "predictor_sha256": row["predictor_sha256"],
    }
    require(all(packet[key] == value for key, value in labels.items()) and case["record_id"] == context.record_id,
            "same record/cache/model/dose and 640-bit packet required")
    if part == "bulk":
        i = index
        nodes = target.a["triangles"][i]
        ids, selected = _mesh_map(context, nodes)
        axis_nodes = sum(target.a["reference"][int(node), 0] == 0 for node in nodes[:3])
        chart = {0: "off_axis", 1: "axis_vertex", 2: "axis_edge_physical"}[axis_nodes]
        J = body["triangles"][i]
        low, high = _packed(J["J_bound"], dyadic=False)
        rho, raw = _fraction(case["rho_exact"]), max(abs(low - 1), abs(high - 1))
        require(J["positive_J_for_entire_box"] is True and low > Q(1, 10**6)
                and list(map(_fraction, case["J_bounds_exact"])) == [low, high]
                and _fraction(case["raw_rho_exact"]) == raw, "element rho must bind to its positive body J")
        _rho(rho, raw)
        require(packet["status"] == "UNIFORM_COHORT_BULK_COMPONENT_UNREVIEWED"
                and case["status"] == "SELECTED_UNIFORM_BULK_ORDER" and packet["element"] == case["element"] == i
                and packet["chart"] == case["chart"] == chart and packet["free_local"] == case["free_local"] == selected
                and packet["active"] == case["active"] == str(int(target.a["active"][i, 0]))
                and packet["nodal_hulls_sha256"] == case["nodal_hulls_sha256"]
                and packet["auxiliary_element_rule_sha256"] == case["auxiliary_element_rule_sha256"]
                and packet["rho_exact"] == case["rho_exact"] and packet["degree"] == case["degree"],
                "bulk packet chart, positivity, selected case or mapping mismatch")
        cache = packet["cache_consistency"]
        require(cache["all64_convex_strain_lines_positive"] is True
                and _packed(cache["minimum_cache_convex_line_J_lower"])[0] > Q(1, 10**6)
                and _packed(cache["max_coefficientwise_delta_q_upper"])[0] >= 0,
                "positive same-node bulk cache line certificate required")
        order, maximum = case["degree"], protocol["max_bulk_order"]
    else:
        face, edge = divmod(index, 12)
        nodes = target.a[f"surface_{face}_edges"][edge]
        ids, selected = _mesh_map(context, nodes)
        domain = case["domain"]
        low, high = map(_fraction, domain["G_bounds_exact"])
        require(domain["status"] == "CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN"
                and 0 < low <= high and _fraction(domain["hoop_lower_exact"]) > 0, "positive whole-edge trace domain required")
        rho, raw = _fraction(domain["rho_exact"]), (high - low) / (high + low)
        _rho(rho, raw)
        require(_fraction(domain["raw_rho_exact"]) == raw and _fraction(domain["q0_exact"]) == (low + high) / 2
                and packet["status"] == "UNIFORM_COHORT_SURFACE_COMPONENT_UNREVIEWED"
                and packet["face"] == case["face"] == face and packet["edge"] == case["edge"] == edge
                and packet["order"] == case["order"] and packet["rho_exact"] == domain["rho_exact"]
                and packet["q0_exact"] == domain["q0_exact"], "surface case or outward expansion mismatch")
        require(len(packet["actual_cache_sites"]) == len(case["cache_site_domains"]) == 12,
                "all twelve surface cache sites required")
        for iq, (v, s) in enumerate(zip(packet["actual_cache_sites"], case["cache_site_domains"])):
            require(v["quadrature_index"] == s["quadrature_index"] == iq and s["global_indices"] == ids
                    and s["selected_local"] == selected and s["nodes"] == list(map(int, nodes))
                    and s["status"] == "CERTIFIED_LOCAL_FULLY_COEFFICIENTWISE_LINE_DOMAIN_NO_CACHE_SUM"
                    and _fraction(s["convex_line_exact"][0][0]) > 0
                    and _fraction(s["tangent_norm_squared_lower_exact"]) > 0
                    and _fraction(v["weight_applied_exactly_once"]) == _fraction(s["auxiliary_normalized_weight_exact"]) > 0,
                    "surface cache site map, positive domain or measure mismatch")
            require(all(_packed(v[key])[0] >= 0 for key in ("unweighted_F_max", "unweighted_H_max")),
                    "surface cache norms must be nonnegative")
        require(all(_packed(packet[key])[0] > 0 for key in ("material_mass", "auxiliary_mass")),
                "positive surface material and auxiliary measures required")
        order, maximum = case["order"], protocol["max_surface_order"]
    require(packet["global_indices"] == case["global_indices"] == ids and bool(ids), "same free-column map required")
    fbudget, hbudget = (_fraction(protocol[k]) for k in ("force_tail_budget_exact", "H_entry_tail_budget_exact"))
    _order_history(order, case["order_attempts"], maximum, fbudget, hbudget)
    require(len(packet["force_two_mass_tail_upper"]) == len(ids), "full force tails required")
    for value in packet["force_two_mass_tail_upper"]:
        _tail(value, fbudget)
    for value in packet["H_integral_minus_cache_upper"]:
        _tail(value["two_mass_tail_upper"], hbudget)
    return ids


def _verify_producer(producer, g, H, RF, H_rows, rf, b, eta0, z0, radius):
    n = len(g)
    for label, values, exact in (
        ("F", producer["force_integral_minus_cache"], g),
        ("RF", producer["preconditioned_force_defect"], RF),
    ):
        require(len(values) == len(exact) == n, f"producer complete {label} required")
        for i, (outer, inner) in enumerate(zip(values, exact)):
            _subset(_packed(outer), inner, f"{label} {i}")
    rows = producer["H_absolute_row_bounds"]
    require(len(rows) == n, "producer complete H row upper bounds required")
    for i, (value, exact) in enumerate(zip(rows, H_rows)):
        lo, hi = _packed(value)
        # The producer sums absolute upper bounds of widened balls. Even the
        # lower endpoint of that sum may exceed the exact endpoint reduction.
        require(0 <= lo and exact <= hi, f"producer H row {i} upper bound fails endpoint domination")
    published = _entry_map(producer["H_integral_minus_cache_upper"], n)
    require(published.keys() == H.keys(), "producer upper Hessian sparsity changed")
    for key, exact in H.items():
        _subset(_packed(published[key]["enclosure"]), exact, f"H {key}")
    producer_rf, producer_b = (_fraction(producer[key]) for key in ("RF_defect_uniform_upper_exact", "b_upper_exact"))
    require(producer_rf >= rf and producer_b >= b, "producer norm scalars fail endpoint domination")
    pe, pz, plhs, ok, perror = _strict_test(eta0, z0, producer_rf, producer_b, radius)
    require(tuple(_fraction(producer[k]) for k in ("eta_upper_exact", "z_upper_exact", "selfmap_lhs_exact", "r_exact"))
            == (pe, pz, plhs, radius) and producer["sufficient_conditions_pass"] is ok
            and producer["status"] == ("CANDIDATE_UNIFORM_EXACT_FE_ROOT_UNREVIEWED" if ok else "UNRESOLVED_SUFFICIENT_ROOT_TEST")
            and producer["root_assertion_independently_reviewed"] is False
            and (producer["root_error_upper_exact"] is None if perror is None
                 else _fraction(producer["root_error_upper_exact"]) == perror), "producer sufficient-test arithmetic mismatch")
    return ok


def _minimum(context, bindings, old_z, z, radius, norm, old_source):
    row = context.row
    summaries = _bound(ENERGY_RESULT, PREPARATION_PINS[ENERGY_RESULT], bindings)
    reviews = _bound(ENERGY_REVIEW, PREPARATION_PINS[ENERGY_REVIEW], bindings)
    require(reviews["status"] == "PASS" and _yes(reviews["checks"]), "cached SPD independent review required")
    (entry,) = [r for r in summaries["records"] if r["record_id"] == context.record_id]
    (wrapped,) = [r for r in reviews["record_audits"] if r["record_id"] == context.record_id]
    cached = _bound(entry["file"], entry["sha256"], bindings)
    audit = _bound(wrapped["review_file"], wrapped["review_sha256"], bindings)
    require(audit == wrapped["audit"] and audit["status"] == "PASS" and _yes(audit["checks"])
            and audit["producer_result_sha256"] == entry["sha256"]
            and _fraction(audit["uniform_Hessian_lower_bound_exact"]) > 0
            and cached["status"] == "CERTIFIED_UNIFORM_STRICT_LOCAL_ENERGY_MINIMUM"
            and cached["whole_tube_Hessian_positive"] is True, "same-record cached SPD premise required")
    root = cached["root_binding"]
    for key in ("record_id", "root_certificate_sha256", "parameter_cell_sha256",
                "model_observation_contract_sha256", "branch_domain_id"):
        expected = row[key]
        require(audit[key] == root[key] == expected, f"cached SPD {key} binding mismatch")
    for key in ("point_cache_file", "point_cache_sha256", "frozen_binary_cache_sha256"):
        require(root[key] == row[key], f"cached SPD {key} mismatch")
    require(root["preconditioner_binary_sha256"] == sha256(context.R.tobytes(order="C")).hexdigest()
            and _fraction(audit["exact_R_infinity_norm"]) == _fraction(cached["exact_R_infinity_norm"]) == norm
            and _fraction(audit["independent_fulltube_Z_upper_exact"]) == old_z
            and _fraction(audit["normal_root_radius_exact"]) == radius
            and root["independent_Z_source"]["file"] == old_source
            and _fraction(cached["uniform_free_coordinate_Hessian_lower_bound_exact"])
                == _fraction(audit["uniform_Hessian_lower_bound_exact"]) > 0
            and 0 <= z < 1 and norm > 0, "same R, Z, radius and cached curvature premise required")
    return {"strict_local_minimum": True, "normalized_curvature_lower_exact": str((1 - z) / norm),
            "cached_independent_review_sha256": wrapped["review_sha256"],
            "argument": "At each parameter-matched integral root, the symmetric homotopy H_Q+t(H_I-H_Q) "
                        "satisfies ||I-RH_t||∞<=z<1. The cached SPD endpoint fixes inertia. "
                        "The inverse infinity-norm bound ||R||∞/(1-z) gives the curvature bound."}


def review(record, admission_sha256):
    _, bindings = _prepared()
    folder = _record_dir(record)
    admission_path = (OUT / f"{record}_input_admission.json").relative_to(ROOT).as_posix()
    admission = _bound(admission_path, admission_sha256, bindings)
    require(admission["status"] == "BOUND_COMPONENT_INPUTS_BEFORE_LOCAL_READ" and admission["record_id"] == record
            and admission["run_dir"] == folder and admission["local_component_payloads_opened"] is False,
            "prewritten terminal input admission required")
    accepted = admission["source_and_manifest_hashes"]
    for name, digest in accepted.items():
        _bound(name, digest, bindings, parse=False)
    protocol, domains, body, assembled = [
        _bound(f"{folder}/{name}", accepted[f"{folder}/{name}"], bindings)
        for name in ("protocol.json", "domains.json", "body_geometry.json", "assembled.json")
    ]
    ctx = C.load_context(record)
    require(protocol["row"] == ctx.row and protocol["precision_bits"] == 640
            and protocol["parameter_cell_sha256"] == ctx.row["parameter_cell_sha256"]
            and body["status"] == "CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING"
            and all(body[key] is True for key in ("all_element_J_positive", "axis_preserved_for_entire_box",
                                                  "uniform_strong_monotonicity", "global_injectivity_for_entire_box"))
            and len(body["triangles"]) == len(domains["bulk_selectors"]) == 144 and len(domains["surface_cases"]) == 24,
            "same full explicit record context and complete fresh domains required")
    box = _bound(protocol["coordinate_box_file"], protocol["coordinate_box_sha256"], bindings)
    expected_box = C.exact_node_box(ctx, ctx.old_radius)
    require(domains["coordinate_box_sha256"] == protocol["coordinate_box_sha256"]
            and all(_digest(box[key]) == _digest(value) for key, value in expected_box.items()),
            "coordinate box must equal the original unshifted exact predictor box")
    eta0, z0, radius, old_source = _old_premises(ctx, protocol, bindings)
    g, H, counts, rowcounts = [ZERO] * 611, {}, [0] * 611, [0] * 611
    for part, cases in (("bulk", domains["bulk_selectors"]), ("surface", domains["surface_cases"])):
        files = admission["expected_packets"][part]
        by_name = {item["file"]: item["sha256"] for item in files}
        for index, case in enumerate(cases):
            face, edge = divmod(index, 12)
            short = f"element_{index:03d}.json" if part == "bulk" else f"face{face}_edge{edge:02d}.json"
            name = f"{folder}/{part}/{short}"
            packet = _bound(name, by_name[name], bindings)
            ids = _packet_meta(ctx, part, index, case, packet, body, protocol)
            _scatter(g, H, ids, packet["force_integral_minus_cache"], packet["H_integral_minus_cache_upper"])
            for i in ids:
                counts[i] += 1
                rowcounts[i] += len(ids)
    require(min(counts) > 0 and counts == protocol["force_incidence"] and rowcounts == protocol["H_row_incidence"]
            and max(counts) == protocol["max_force_incidence"] and max(rowcounts) == protocol["max_H_row_incidence"],
            "full record component incidence required")
    H_rows = _h_rows(611, H)
    RF, b_rows, norm = _dyadic_reduce(ctx.R, g, H_rows)
    require(norm == _fraction(protocol["R_infinity_norm_exact"]), "exact binary preconditioner norm mismatch")
    divisor = _fraction(protocol["margin_divisor"])
    fbudget, hbudget = (_fraction(protocol[k]) for k in ("force_tail_budget_exact", "H_entry_tail_budget_exact"))
    require(divisor > 0 and norm * max(counts) * fbudget == (1 - z0) * radius / divisor
            and norm * max(rowcounts) * hbudget == (1 - z0) / divisor, "exact global tail budget accounting mismatch")
    rf, b = max(map(_abs_upper, RF)), max(b_rows)
    producer_pass = _verify_producer(assembled, g, H, RF, H_rows, rf, b, eta0, z0, radius)
    eta, z, lhs, passed, error = _strict_test(eta0, z0, rf, b, radius)
    minimum = _minimum(ctx, bindings, z0, z, radius, norm, old_source) if passed else None
    output = OUT / f"{record}_reduction_review.json"
    report = {
        "status": "PASS_CERTIFIED_UNIFORM_EXACT_FE_ROOT" if passed else "PASS_REDUCTION_UNRESOLVED_SUFFICIENT_ROOT_TEST",
        "record_id": record, "input_admission_sha256": admission_sha256, "producer_sufficient_test_pass": producer_pass,
        "assembly_enclosure_containment_verified": True, "whole_cell_original_predictor_cube": True,
        "component_count": 168, "force_entries": len(g), "upper_Hessian_entries": len(H),
        "F_endpoint_sha256": _digest(g), "H_endpoint_sha256": _digest([(i, j, H[i, j]) for i, j in sorted(H)]),
        "RF_endpoint_sha256": _digest(RF), "abs_R_H_row_sha256": _digest(b_rows),
        "R_infinity_norm_exact": norm, "RF_uniform_upper_exact": rf, "b_upper_exact": b,
        "cached_Z_upper_exact": z0, "cached_eta_upper_exact": eta0,
        "eta_upper_exact": eta, "z_upper_exact": z, "selfmap_lhs_exact": lhs, "normal_radius_exact": radius,
        "sufficient_conditions_pass": passed, "root_error_upper_exact": error, "energy": minimum,
        "source_input_and_component_hashes": bindings,
        "independence": "Exact endpoint assembly and integer dyadic reduction are independent of producer Arb assembly. "
                        "The checked record/context loader, reviewed geometry and same local integration kernels are shared premises.",
    }
    digest = _write(output, report)
    print(json.dumps({"status": report["status"], "file": str(output), "sha256": digest}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest="action", required=True)
    admission = actions.add_parser("admit")
    admission.add_argument("record")
    for key in PRIMARY_PINS:
        admission.add_argument("--" + key.replace("_", "-") + "-sha256", required=True)
    replay = actions.add_parser("review")
    replay.add_argument("record")
    replay.add_argument("--admission-sha256", required=True)
    args = parser.parse_args()
    if args.action == "admit":
        admit(args.record, {key: getattr(args, key + "_sha256") for key in PRIMARY_PINS})
    else:
        review(args.record, args.admission_sha256)
