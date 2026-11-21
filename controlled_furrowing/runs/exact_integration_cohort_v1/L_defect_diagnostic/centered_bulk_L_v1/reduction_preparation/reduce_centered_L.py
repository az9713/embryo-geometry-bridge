"""Exact composition of the complete centered bulk force and retained L components."""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import lcm
import os
from pathlib import Path
import sys
from time import monotonic


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[4]
sys.path.insert(0, str(ROOT / "scripts"))
import exact_centered_bulk_L as G


D, B = G.D, G.B
RECORD, ORIGINAL = G.RECORD, G.ORIGINAL
STAGE = OUT.parent.relative_to(ROOT).as_posix()
SOURCE_PIN = "a83c89ffc302be341df220fcb97c5be82f191a70e317fea9af839d2917a87492"
PINS = {
    "scripts/exact_centered_bulk_L.py": SOURCE_PIN,
    f"{STAGE}/preparation/protocol.json": "466048ec456499aa904dd50daa744245aff2eecb6fc2a30fb15bec3d97d23f9b",
    f"{STAGE}/preparation/preflight.json": "2ed81adb1f8dd7e52d7c49685be49e508af4ee395911c0e292a74c67148632e0",
    f"{STAGE}/preparation/selectors.json": "8fbb0e8844d7f47f119b795e453032ae46d24a53ea489d2ce9fe2119576a9bf6",
    f"{STAGE}/preparation/prelaunch_independent_review.json": "c915d5b2f54fe16b7e72d3f1e8bb4d10b94fca9a1c42209f21e7c2b032c16297",
}


def _relative(path):
    return Path(path).relative_to(ROOT).as_posix()


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _multiply(a, b):
    product = [x * y for x in a for y in b]
    return min(product), max(product)


def _intersect(a, b):
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    D.require(lo <= hi, "independent signed interval bounds are disjoint")
    return lo, hi


def _sum(*values):
    return sum((a for a, _ in values), Q(0)), sum((b for _, b in values), Q(0))


def _active_bound(box):
    return min(map(abs, box)) if box[0] * box[1] > 0 else Q(0)


def _prepared():
    hashes = {}
    docs = {p: D._bound(p, h, hashes, parse=p.endswith(".json")) for p, h in PINS.items()}
    protocol = docs[f"{STAGE}/preparation/protocol.json"]
    pre = docs[f"{STAGE}/preparation/preflight.json"]
    selectors = docs[f"{STAGE}/preparation/selectors.json"]
    gate = docs[f"{STAGE}/preparation/prelaunch_independent_review.json"]
    expected = {
        "source_sha256": SOURCE_PIN,
        "protocol_sha256": PINS[f"{STAGE}/preparation/protocol.json"],
        "preflight_sha256": PINS[f"{STAGE}/preparation/preflight.json"],
        "selectors_sha256": PINS[f"{STAGE}/preparation/selectors.json"],
        "prelaunch_sha256": PINS[f"{STAGE}/preparation/prelaunch_independent_review.json"],
    }
    D.require(_sha(G.__file__) == SOURCE_PIN and gate["status"] == "PASS" and D._yes(gate["checks"])
              and gate["record_id"] == protocol["record_id"] == RECORD
              and all(gate[key] == value for key, value in expected.items() if key != "prelaunch_sha256")
              and pre["status"] == "PASS_CENTERED_BULK_L_PREPARATION_ONLY" and D._yes(pre["checks"])
              and protocol["elements"] == 144 and [p["element"] for p in selectors] == list(range(144)),
              "same reviewed 144-element centered preparation required")
    for path, digest in protocol["input_hashes"].items():
        D._bound(path, digest, hashes, parse=False)
    D._bound(gate["support_file"], gate["support_sha256"], hashes, parse=False)
    return hashes, protocol, selectors, expected


def _check_gate(expected_sha, hashes):
    gate_path = _relative(OUT / "prelaunch_independent_review.json")
    source, protocol, pre = map(_sha, (Path(__file__), OUT / "protocol.json", OUT / "preflight.json"))
    review = D._bound(gate_path, expected_sha, hashes)
    D.require(review["status"] == "PASS" and D._yes(review["checks"])
              and (review["source_sha256"], review["protocol_sha256"], review["preflight_sha256"])
                  == (source, protocol, pre), "matching independent prelaunch of this exact composition required")
    for path, digest in ((_relative(__file__), source), (_relative(OUT / "protocol.json"), protocol),
                         (_relative(OUT / "preflight.json"), pre)):
        D._bound(path, digest, hashes, parse=False)


def _manifest(terminal_sha, manifest_sha, protocol, expected, hashes):
    terminal = D._bound(f"{STAGE}/actual/terminal.json", terminal_sha, hashes)
    D._terminal(terminal, manifest_sha)
    manifest = D._bound(f"{STAGE}/actual/results.json", manifest_sha, hashes)
    fields = {**expected, "record_id": RECORD, "status": "ALL_144_CENTERED_BULK_FORCE_PACKETS_UNREVIEWED",
              "new_Hessian_integrations": False, "root_asserted": False}
    return B._validate_manifest(manifest, range(144), protocol["actual_expected_files"], fields)


def admit(review_sha, terminal_sha, manifest_sha):
    hashes, protocol, _, expected = _prepared()
    _check_gate(review_sha, hashes)
    files = _manifest(terminal_sha, manifest_sha, protocol, expected, hashes)
    out = {"status": "BOUND_COMPLETE_CENTERED_MANIFEST_BEFORE_LOCAL_READ", "record_id": RECORD,
           "prelaunch_review_sha256": review_sha, "centered_terminal_sha256": terminal_sha,
           "centered_manifest_sha256": manifest_sha, "input_and_manifest_hashes": hashes,
           "expected_centered_packets": files, "local_centered_payloads_opened": False}
    digest = D._write(OUT / "input_admission.json", out)
    print(json.dumps({"admission_sha256": digest, "centered_bulk_components": len(files)}), flush=True)


def _packet_metadata(packet, proof, expected, row, dose):
    values = {**expected, **proof, "status": "CENTERED_BULK_FORCE_COMPONENT_UNREVIEWED",
              "record_id": RECORD, "dose_exact": str(dose), "parameter_cell_sha256": row["parameter_cell_sha256"],
              "precision_bits": 640, "local_proof_sha256": D._digest(proof), "new_Hessian_integrations": False}
    D.require(all(packet.get(key) == value for key, value in values.items()),
              "mixed record, center, original/scaled evidence, chart or local geometry packet")


def _gvalues(array, count):
    D.require(len(array) == count, "complete local force coverage required")
    return list(map(D._packed, array))


def _centered_packet(packet, proof, expected, row, dose, prior, H, radii, gain_bounds):
    _packet_metadata(packet, proof, expected, row, dose)
    count = len(proof["global_indices"])
    point = packet["point_force"]
    vol, shear, cache, posted_point = [
        _gvalues(point[key], count) for key in ("volumetric_truncation", "quadratic_strain_defect",
                  "finite_same_Y_cache_correction", "normalized_force_total")
    ]
    tail = list(map(Q, point["two_mass_tail_upper_exact"]))
    D.require(len(tail) == count and min(tail) >= 0 and max(tail) == Q(proof["point_force_tail_max_exact"])
              and max(tail) <= Q(proof["original_force_tail_budget_exact"])
              and D._packed(point["cache_min_J_bound"])[0] > Q(1, 10**6), "same point tails and analytic cache required")
    for v, a, b, c, t in zip(posted_point, vol, shear, cache, tail):
        _intersect(v, _sum(a, b, c, (-t, t)))
    point_force = posted_point
    gain = _gvalues(packet["affine_gain_derivative_at_Y"], count)
    if proof["active"] == "0":
        D.require(all(b == D.ZERO for b in gain), "inactive element gain effect must be zero")
    state = _gvalues(packet["same_original_box_Hessian_times_displacement"], count)
    gain_piece = _gvalues(packet["affine_gain_increment"], count)
    centered = _gvalues(packet["centered_force_enclosure"], count)
    intersection = _gvalues(packet["force_intersection_with_prior"], count)
    old = _gvalues(prior, count)
    for i in range(count):
        u = sum((D._abs_upper(H[min(i, j), max(i, j)]) * x for j, x in enumerate(radii)), Q(0))
        D.require(D._contains(state[i], (-u, u)) and D._contains(gain_piece[i], _multiply(gain[i], gain_bounds)),
                  "local centered increments fail exact same-box signed reconstruction")
        _intersect(centered[i], _sum(point_force[i], state[i], gain_piece[i]))
        D.require(intersection[i] == _intersect(centered[i], old[i]), "local signed force intersection changed")
    return intersection, point_force, gain


def _force_scatter(array, ids, values):
    D.require(len(ids) == len(values), "exact mapped force coverage required")
    for i, value in zip(ids, values):
        array[i] = D._plus(array[i], value)


def _bulk_RH_effects(R, H, radii, predictor):
    """Exact signed RH_bulk, then independent nodal radii or shared predictor columns."""
    n = len(radii)
    D.require(n and len(R) == n and all(len(row) == n for row in R) and min(radii) >= 0,
              "matching nonnegative nodal radii and square R required")
    half = [[] for _ in range(n)]
    shift_h = max((D._exp(v) for value in H.values() for v in value), default=0)
    for (i, j), (lo, hi) in H.items():
        a, b = lo.numerator << (shift_h - D._exp(lo)), hi.numerator << (shift_h - D._exp(hi))
        half[j].append((i, a, b))
        if i != j:
            half[i].append((j, a, b))
    pairs = [[float(v).as_integer_ratio() for v in row] for row in R]
    shift_r = max(den.bit_length() - 1 for row in pairs for _, den in row)
    common = lcm(*(r.denominator for r in radii))
    weights = [r.numerator * (common // r.denominator) for r in radii]
    p = len(predictor[0])
    D.require(len(predictor) == n and p and all(len(row) == p for row in predictor),
              "same explicit predictor projection required")
    s_pairs = [[float(v).as_integer_ratio() for v in row] for row in predictor]
    shift_s = max(den.bit_length() - 1 for row in s_pairs for _, den in row)
    slope = [[num << (shift_s - (den.bit_length() - 1)) for num, den in row] for row in s_pairs]
    results, coupled = [], []
    denominator = common * 2**(shift_h + shift_r)
    predictor_denominator = 2**(shift_h + shift_r + shift_s)
    for row in pairs:
        coef = [num << (shift_r - (den.bit_length() - 1)) for num, den in row]
        value = 0
        lower, upper = [0] * p, [0] * p
        for j, column in enumerate(half):
            lo = hi = 0
            for i, a, b in column:
                k = coef[i]
                lo += k * (a if k >= 0 else b)
                hi += k * (b if k >= 0 else a)
            value += max(abs(lo), abs(hi)) * weights[j]
            for a, k in enumerate(slope[j]):
                lower[a] += k * (lo if k >= 0 else hi)
                upper[a] += k * (hi if k >= 0 else lo)
        results.append(Q(value, denominator))
        coupled.append([(Q(a, predictor_denominator), Q(b, predictor_denominator)) for a, b in zip(lower, upper)])
    return results, coupled


def _coupled_predictor_force(point, A, gain, surface, delta, gain_column):
    """Reuse one delta-g coordinate for the shape predictor and affine gain."""
    D.require(len(point) == len(A) == len(gain) == len(surface)
              and 0 <= gain_column < len(delta) and all(len(row) == len(delta) for row in A),
              "same preconditioned predictor and gain dimensions required")
    rows = []
    for p, v, g, s in zip(point, A, gain, surface):
        terms = [(_sum(value, g) if j == gain_column else value) for j, value in enumerate(v)]
        rows.append(_sum(p, s, *(_multiply(a, b) for a, b in zip(terms, delta))))
    return rows


def _actual(admission, protocol, proofs, expected, hashes):
    data = G._source_inputs()
    D.require(data.hashes == protocol["input_hashes"], "same prior reviewed combined H and frozen box binding required")
    f = [D.ZERO] * 611
    Hbulk, H = {}, {}
    point, gain, surface = ([D.ZERO] * 611 for _ in range(3))
    count, hcount = [0] * 611, [0] * 611

    def scatter(ids, values, local_H, all_H):
        D._scatter(f, all_H, ids, values, local_H)
        for i in ids:
            count[i] += 1
            hcount[i] += len(ids)

    for proof, binding in zip(proofs, admission["expected_centered_packets"]):
        i = proof["element"]
        rebuilt, _, _, _, _, local_H, prior, radii = G._element(data, i)
        D.require(D._digest(rebuilt) == D._digest(proof) and binding["element"] == i, "changed or crossed point proof")
        packet = D._bound(binding["file"], binding["sha256"], hashes)
        interval, pf, ag = _centered_packet(packet, proof, expected, data.context.row, data.context.dose,
                                           prior, local_H, radii, data.gain_bounds)
        ids = proof["global_indices"]
        _force_scatter(point, ids, pf)
        _force_scatter(gain, ids, ag)
        Hrows = [{"i": a, "j": b, "enclosure": B._packed(*value)} for (a, b), value in sorted(local_H.items())]
        scatter(ids, [B._packed(*p) for p in interval], Hrows, Hbulk)
    H.update(Hbulk)
    old = D._bound(B.OLD + "_input_admission.json", B.PINS[B.OLD + "_input_admission.json"], hashes)
    for i, case in enumerate(data.domains["surface_cases"]):
        face, edge = divmod(i, 12)
        name = f"{ORIGINAL}/surface/face{face}_edge{edge:02d}.json"
        binding = next(p for p in old["expected_packets"]["surface"] if p["file"] == name)
        packet = D._bound(name, binding["sha256"], hashes)
        ids = D._packet_meta(data.context, "surface", i, case, packet, data.body, data.protocol)
        g = packet["force_integral_minus_cache"]
        scatter(ids, g, packet["H_integral_minus_cache_upper"], H)
        _force_scatter(surface, ids, list(map(D._packed, g)))
    D.require(len(proofs) == 144 and len(data.domains["surface_cases"]) == 24 and len(f) == 611 and len(H) == 6676
              and count == data.protocol["force_incidence"] and hcount == data.protocol["H_row_incidence"],
              "incomplete centered 168-component assembly")
    old_name = f"{G.SCALED}/reduction_preparation/results.json"
    scaled = D._bound(old_name, G.PIN[old_name], hashes)
    Hhash = D._digest([(i, j, H[i, j]) for i, j in sorted(H)])
    D.require(Hhash == scaled["H_endpoint_sha256"], "retained full Hessian endpoints changed")
    R = data.context.R
    rf_direct, b_rows, norm = D._dyadic_reduce(R, f, D._h_rows(611, H))
    b = max(b_rows)
    D.require(b == Q(scaled["b_upper_exact"]) and norm == Q(scaled["R_infinity_norm_exact"]),
              "unchanged reviewed normalization and Hessian bound required")
    zeros = [Q(0)] * 611
    rp, rg, rs = [D._dyadic_reduce(R, v, zeros)[0] for v in (point, gain, surface)]
    nodal_radii = [Q(data.box["coordinate_radii_exact"][int(dof) // 2][int(dof) % 2]) for dof in data.context.free]
    state, coupled = _bulk_RH_effects(R, Hbulk, nodal_radii, data.context.S)
    gain = [_multiply(v, data.gain_bounds) for v in rg]
    separated = [_sum(p, (-u, u), g, s) for p, u, g, s in zip(rp, state, gain, rs)]
    joint = [_intersect(a, b) for a, b in zip(rf_direct, separated)]
    theta_order = G.T.ORDER
    delta = [(a - data.context.base[key], b - data.context.base[key])
             for key in theta_order for a, b in [data.context.cell[key]]]
    predictor_only = _coupled_predictor_force(rp, coupled, rg, rs, delta, theta_order.index("common_gain"))
    predictor_joint = [_intersect(a, b) for a, b in zip(joint, predictor_only)]
    force = max(map(D._abs_upper, predictor_joint))
    eta0, z0, r, provider = D._old_premises(data.context, data.protocol, hashes)
    eta, z, lhs, passes, error = D._strict_test(eta0, z0, force, b, r)
    curvature = D._minimum(data.context, hashes, z0, z, r, norm, provider) if passes else None
    return {
        "status": "CONDITIONAL_UNREVIEWED_CENTERED_ORIGINAL_CUBE_SUFFICIENT_TEST_PASSES" if passes
                  else "UNREVIEWED_CENTERED_ORIGINAL_CUBE_SUFFICIENT_TEST_UNRESOLVED",
        "record_id": RECORD, "component_count": 168, "centered_bulk_components": 144,
        "retained_surface_components": 24, "force_entries": 611, "upper_Hessian_entries": 6676,
        "F_endpoint_sha256": D._digest(f), "uncentered_H_endpoint_sha256": Hhash,
        "bulk_only_H_endpoint_sha256": D._digest([(i, j, Hbulk[i, j]) for i, j in sorted(Hbulk)]),
        "direct_whole_box_RF_intervals_exact": rf_direct, "separated_whole_box_RF_intervals_exact": separated,
        "whole_box_RF_intersections_exact": joint, "coupled_predictor_only_RF_intervals_exact": predictor_only,
        "final_predictor_only_RF_intersections_exact": predictor_joint,
        "separated_point_bulk_RF_intervals_exact": rp, "separated_state_bulk_only_radius_exact": state,
        "separated_gain_derivative_RF_intervals_exact": rg, "separated_gain_term_exact": gain,
        "separated_retained_surface_RF_intervals_exact": rs,
        "coupled_RH_bulk_S_intervals_exact": coupled, "predictor_parameter_order": list(theta_order),
        "exact_closed_parameter_offsets": delta, "coupled_predictor_only_RF_uniform_upper_exact": max(map(D._abs_upper, predictor_only)),
        "direct_RF_uniform_upper_exact": max(map(D._abs_upper, rf_direct)),
        "separated_RF_uniform_upper_exact": max(map(D._abs_upper, separated)),
        "whole_box_RF_uniform_upper_exact": max(map(D._abs_upper, joint)),
        "RF_uniform_upper_exact": force,
        "predictor_only_RF_uniform_lower_exact": max(map(_active_bound, predictor_joint)),
        "RF_uniform_upper_scope": "All parameter-matched exact predictors p(theta), theta in unchanged C.",
        "R_infinity_norm_exact": norm, "b_upper_exact": b, "cached_eta_upper_exact": eta0, "cached_Z_upper_exact": z0,
        "eta_upper_exact": eta, "z_upper_exact": z, "selfmap_lhs_exact": lhs, "normal_radius_exact": r,
        "sufficient_conditions_pass": passes, "conditional_root_error_upper_exact": error, "conditional_curvature": curvature,
        "same_point_cell_and_cube": True, "root_assertion_independently_reviewed": False,
        "primitives": "Reuses original exact endpoint scatter and integer-dyadic R-force and |R|-absolute-H-row kernels. "
                      "Separated reduction uses signed RH_bulk first, componentwise fixed nodal radii, shared gain applied "
                      "after signed RG and retained surface force. Predictor-only bound instead uses RH_bulk S "
                      "and adding RG in the same gain column before applying all five shared parameter intervals. "
                      "Surface H contributes to b only.",
    }


def reduce(admission_sha):
    hashes, protocol, proofs, expected = _prepared()
    admission = D._bound(_relative(OUT / "input_admission.json"), admission_sha, hashes)
    D.require(admission["status"] == "BOUND_COMPLETE_CENTERED_MANIFEST_BEFORE_LOCAL_READ"
              and admission["record_id"] == RECORD and admission["local_centered_payloads_opened"] is False,
              "complete frozen centered admission required before opening packets")
    _check_gate(admission["prelaunch_review_sha256"], hashes)
    for path, digest in admission["input_and_manifest_hashes"].items():
        D._bound(path, digest, hashes, parse=False)
    D.require(_manifest(admission["centered_terminal_sha256"], admission["centered_manifest_sha256"],
                        protocol, expected, hashes) == admission["expected_centered_packets"], "same admitted packet set required")
    stamp = lambda: datetime.now(timezone.utc).isoformat()
    start = monotonic()
    D._write(OUT / "started.json", {"pid": os.getpid(), "utc": stamp(), "admission_sha256": admission_sha})
    try:
        result = _actual(admission, protocol, proofs, expected, hashes)
        result.update({"admission_sha256": admission_sha, "source_inputs_and_packets": hashes, "seconds": monotonic() - start})
        digest = D._write(OUT / "results.json", result)
        D._write(OUT / "terminal.json", {"status": "FINISHED", "exit_code": 0, "pid": os.getpid(),
                                        "results_sha256": digest, "utc": stamp()})
        print(json.dumps({"status": result["status"], "results_sha256": digest}), flush=True)
    except BaseException as error:
        D._write(OUT / "failure.json", {"message": str(error), "admission_sha256": admission_sha})
        D._write(OUT / "terminal.json", {"status": "FAILED", "exit_code": 1, "pid": os.getpid(), "utc": stamp()})
        raise


def tests():
    checks = {}

    def check(key, expected):
        D.require(expected, key)
        checks[key] = True

    def reject(key, action):
        try:
            action()
        except ValueError:
            checks[key] = True
        else:
            raise AssertionError(key)

    x, y, z = (Q(-3), Q(-1)), (Q(-2), Q(2)), (Q(0), Q(1))
    check("strict_negative_and_zero_crossing_intersection", _intersect(x, y) == (Q(-2), Q(-1)))
    reject("disjoint_signed_intersection_rejected", lambda: _intersect(x, z))
    check("signed_scalar_gain_product", _multiply(x, y) == (Q(-6), Q(6)))
    H = {(0, 0): (Q(1), Q(1)), (0, 1): (Q(-1), Q(-1)), (1, 1): (Q(1), Q(1))}
    R = [[1.0, 1.0], [-1.0, 1.0]]
    S = [[1.0, 0.0, -0.5, 0.5, 1.0], [0.0, 1.0, 0.5, 0.0, -1.0]]
    radii, coupling = _bulk_RH_effects(R, H, [Q(1, 3), Q(1, 7)], S)
    check("signed_RH_bulk_cancels_before_radius_multiplication", radii == [Q(0), Q(20, 21)])
    check("five_exact_signed_predictor_columns_retained",
          coupling[0] == [D.ZERO] * 5 and coupling[1] == [(Q(v), Q(v)) for v in (-2, 2, 2, -1, -4)])
    delta = [(Q(0), Q(0)), (Q(-1), Q(1)), (Q(0), Q(1)), (Q(1), Q(1)), (Q(-2), Q(2))]
    coupled = _coupled_predictor_force([D.ZERO], [coupling[1]], [(Q(4), Q(4))], [D.ZERO], delta, 4)
    check("same_gain_column_cancels_before_mixed_five_parameter_corners",
          coupled == [(Q(-3), Q(3))] and coupled != [(Q(-19), Q(19))])
    H[0, 1] = (Q(-1), Q(1))
    check("RH_bulk_interval_crossing_with_nondyadic_radii",
          _bulk_RH_effects(R, H, [Q(1, 3), Q(1, 7)], S)[0] == [Q(20, 21)] * 2)
    # An exactly canceling surface H is inadmissible in this bulk-only segment.
    check("bulk_state_nonzero_even_if_full_H_zero",
          _bulk_RH_effects([[1.0]], {(0, 0): (Q(2), Q(2))}, [Q(1, 4)], [[1.0]])[0] == [Q(1, 2)])
    reject("negative_radius_rejected", lambda: _bulk_RH_effects([[1.0]], {}, [Q(-1)], [[1.0]]))
    check("strict_original_cube_sufficient_test", D._strict_test(Q(1, 8), Q(1, 4), Q(1, 8), Q(1, 8), Q(1))[-2:] == (True, Q(2, 5)))
    check("selfmap_at_equality_unresolved", D._strict_test(Q(1, 4), Q(1, 4), Q(0), Q(1, 4), Q(1, 2))[-2:] == (False, None))
    files = [{"element": 1, "file": "a", "sha256": "a" * 64}, {"element": 2, "file": "b", "sha256": "b" * 64}]
    v = lambda f: B._validate_manifest({"files": f}, [1, 2], ["a", "b"], {})
    check("complete_packet_manifest_accepted", v(files) == files)
    reject("partial_manifest_rejected", lambda: v(files[:1]))
    reject("duplicate_manifest_rejected", lambda: v(files[:1] * 2))
    reject("mixed_packet_element_rejected", lambda: v([{**files[0], "element": 0}, files[1]]))
    reject("substituted_path_rejected", lambda: v([{**files[0], "file": "x"}, files[1]]))
    reject("failed_terminal_rejected",
           lambda: D._terminal({"status": "FAILED", "pid": 1, "exit_code": 1, "results_sha256": "a" * 64}))
    proof = {"element": 1, "original_packet_sha256": "a" * 64}
    row = {"parameter_cell_sha256": "b" * 64}
    packet = {**proof, "status": "CENTERED_BULK_FORCE_COMPONENT_UNREVIEWED", "record_id": RECORD,
              "dose_exact": "2", "parameter_cell_sha256": "b" * 64, "precision_bits": 640,
              "local_proof_sha256": D._digest(proof), "new_Hessian_integrations": False}
    _packet_metadata(packet, proof, {}, row, Q(2))
    check("complete_exact_packet_metadata", True)
    reject("mixed_record_packet_rejected", lambda: _packet_metadata({**packet, "record_id": "other"}, proof, {}, row, Q(2)))
    reject("mixed_local_proof_rejected", lambda: _packet_metadata({**packet, "original_packet_sha256": "c" * 64}, proof, {}, row, Q(2)))
    return checks


def prepare():
    hashes, _, _, _ = _prepared()
    checks = tests()
    checks["frozen_144_and_prior_reviewed_corrected_domain_current"] = True
    source = Path(__file__).resolve()
    names = ("reduce_centered_L_frozen.py", "protocol.json", "preflight.json")
    D.require(not any((OUT / n).exists() for n in names), "preserve write-once composition preparation")
    with (OUT / names[0]).open("xb") as fp:
        fp.write(source.read_bytes())
    p = {
        "status": "FROZEN_CENTERED_L_REDUCTION_PREPARATION_REVIEW_PENDING", "record_id": RECORD,
        "source_file": _relative(source), "source_sha256": _sha(source), "input_hashes": hashes,
        "force": "Each of the 144 centered packets binds exact local proof metadata and signed force fields. "
                 "Use exact local centered/prior intersections, then assemble all signed entries. "
                 "Independently compose a separate R-force bound with point-bulk, signed RH_bulk acting on original "
                 "node radii, shared exact gain after signed R times bulk gain derivative, and original surface force.",
        "predictor_only": "Also enclose only parameter-matched predictor F. Form the signed 611 by 5 interval RH_bulk*S. "
                          "Add the signed preconditioned gain derivative to its gain column before multiplying by "
                          "the same five closed parameter increments. Add point bulk and retained surface. "
                          "Only the resulting predictor-only force norm enters eta. Contraction still uses full box H.",
        "Hessian": "Reuse full same 132 scaled/original intersected + 12 original axis + 24 original surface Hessian. "
                   "Verify unchanged endpoints and b. Use only the bulk portion for the separated bulk state segment.",
        "gate": "Terminal manifest and independent reduction prelaunch expected SHA supplied externally before payload reads. "
                "Even a passing strict arithmetic test is conditional pending independent actual review.",
        "component_count": 168, "force_entries": 611, "upper_Hessian_entries": 6676,
        "actual_centered_packet_payloads_opened": False, "actual_reduction_run": False,
        "failed_attempts": {_relative(p): _sha(p) for p in sorted(OUT.glob("development_revision_*/*")) if p.is_file()},
    }
    psha = D._write(OUT / names[1], p)
    pre = {"status": "PASS_CENTERED_L_COMPOSITION_PREPARATION_ONLY", "checks": checks, "check_count": len(checks),
           "source_sha256": p["source_sha256"], "protocol_sha256": psha, "actual_reduction_run": False}
    presha = D._write(OUT / names[2], pre)
    print(json.dumps({"status": pre["status"], "source_sha256": p["source_sha256"], "protocol_sha256": psha,
                      "preflight_sha256": presha, "checks": len(checks)}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest="command", required=True)
    command.add_parser("prepare")
    freeze = command.add_parser("admit")
    for flag in ("review", "terminal", "manifest"):
        freeze.add_argument(f"--{flag}-sha256", required=True)
    run = command.add_parser("reduce")
    run.add_argument("--admission-sha256", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare()
    elif args.command == "admit":
        admit(args.review_sha256, args.terminal_sha256, args.manifest_sha256)
    else:
        reduce(args.admission_sha256)
