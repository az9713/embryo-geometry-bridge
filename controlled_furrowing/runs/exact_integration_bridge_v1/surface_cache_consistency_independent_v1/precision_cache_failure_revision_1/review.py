"""Independent local surface consistency review; no integration or solves.

Explicit analytic first/second/third derivatives and rational nodal maps are
implemented here. Arb arithmetic itself, decoding and frozen inputs are reused.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import datetime
import hashlib
import json
import sys
import time

import numpy as np
from flint import arb, arb_mat, ctx

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import exact_surface_cache_consistency as P

PREP = ROOT / "runs/exact_integration_bridge_v1/surface_cache_consistency_preparation"


def read(path):
    return json.loads(path.read_text(encoding="utf8"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, obj):
    if path.exists():
        raise RuntimeError("preserve write-once output: " + str(path))
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True, allow_nan=False) + "\n", encoding="utf8")


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def A(v):
    return P.A(v)


def interval(lo, hi):
    return A(lo).union(A(hi))


def mat(rows):
    return arb_mat(len(rows), len(rows[0]), [A(x) for row in rows for x in row])


def full_derivatives(strain, K, tension):
    """Independent closed forms, using s=sqrt(v.v), v=(b,c)."""
    a, b, c = map(A, strain)
    K, tension = A(K), A(tension)
    v = [b, c]
    s2 = b * b + c * c
    if not a > 0 or not s2 > 0:
        raise ValueError("positive hoop and tangent norm required")
    s = s2.sqrt()
    d = tension - K
    value = K * a * a * s2 / 2 + d * a * s + K / 2
    g = [K * a * s2 + d * s] + [K * a * a * x + d * a * x / s for x in v]
    H = [[arb(0) for _ in range(3)] for _ in range(3)]
    H[0][0] = K * s2
    for i in range(2):
        H[0][i + 1] = H[i + 1][0] = (2 * K * a + d / s) * v[i]
        for j in range(2):
            H[i + 1][j + 1] = K * a * a * int(i == j) + d * a * (int(i == j) / s - v[i] * v[j] / s**3)
    D = [[[arb(0) for _ in range(3)] for _ in range(3)] for _ in range(3)]
    for indices in product(range(3), repeat=3):
        count = indices.count(0)
        if count == 3:
            value3 = arb(0)
        elif count == 2:
            value3 = 2 * K * v[next(i for i in indices if i != 0) - 1]
        elif count == 1:
            i, j = [i - 1 for i in indices if i != 0]
            value3 = 2 * K * a * int(i == j) + d * (int(i == j) / s - v[i] * v[j] / s**3)
        else:
            i, j, k = [i - 1 for i in indices]
            value3 = d * a * (-(int(i == j) * v[k] + int(i == k) * v[j] + int(j == k) * v[i]) / s**3 + 3 * v[i] * v[j] * v[k] / s**5)
        D[indices[0]][indices[1]][indices[2]] = value3
    return value, g, H, D


def pullback(B, Y, K, tension):
    strains = [sum((A(b) * A(y) for b, y in zip(row, Y)), arb(0)) for row in B]
    _, g, H, _ = full_derivatives(strains, K, tension)
    M = mat(B)
    return M.transpose() * arb_mat(3, 1, g), M.transpose() * mat(H) * M


def own_maps(t, face, edge, iq):
    a = t.a
    nodes = [int(i) for i in a[f"surface_{face}_edges"][edge]]
    radius_lo = P.q(a["reference"][nodes[0], 0])
    L = P.q(a["reference"][nodes[2], 0]) - radius_lo
    nold = list(map(P.q, a[f"surface_{face}_N"][iq]))
    dold = list(map(P.q, a[f"surface_{face}_dN"][iq]))
    r0 = P.q(a[f"surface_{face}_r0"][edge, iq])
    l0 = P.q(a[f"surface_{face}_l0"][edge])
    s = nold[2] + nold[1] / 2
    R = radius_lo + L * s
    # Independent expanded quadratic basis and derivative.
    phi = [2 * s * s - 3 * s + 1, 4 * s - 4 * s * s, 2 * s * s - s]
    dphi = [4 * s - 3, 4 - 8 * s, 4 * s - 1]
    Bs, Bc, raw = [[[Q(0)] * 6 for _ in range(3)] for _ in range(3)]
    D = [L / (r0 * l0), 1 / L, 1 / L]
    for j in range(3):
        Bs[0][2 * j] = phi[j] / R
        Bs[1][2 * j] = dphi[j] / L
        Bs[2][2 * j + 1] = dphi[j] / L
        raw[0][2 * j] = nold[j]
        raw[1][2 * j] = dold[j]
        raw[2][2 * j + 1] = dold[j]
    Bc = [[v * D[i] for v in row] for i, row in enumerate(raw)]
    selected = [k for k in range(6) if 2 * nodes[k // 2] + k % 2 in t.index]
    return Bs, Bc, raw, D, nodes, selected, R, L


def select_cols(B, selected):
    return [[row[k] for k in selected] for row in B]


def analytic_difference_bound(Bs, Bc, Y, selected, K, tension):
    delta = [[s - c for s, c in zip(sr, cr)] for sr, cr in zip(Bs, Bc)]
    qc = [sum((A(c) * y for c, y in zip(row, Y)), arb(0)) for row in Bc]
    dq = [sum((A(c) * y for c, y in zip(row, Y)), arb(0)) for row in delta]
    line = [x + arb(0).union(dx) for x, dx in zip(qc, dq)]
    _, gc, Hc, _ = full_derivatives(qc, K, tension)
    _, _, Hl, Dl = full_derivatives(line, K, tension)
    dg = [sum((Hl[i][k] * dq[k] for k in range(3)), arb(0)) for i in range(3)]
    dH = [[sum((Dl[i][j][k] * dq[k] for k in range(3)), arb(0)) for j in range(3)] for i in range(3)]
    Ms, Mc, Md = [mat(select_cols(B, selected)) for B in (Bs, Bc, delta)]
    first = Md.transpose() * arb_mat(3, 1, gc)
    second = Ms.transpose() * arb_mat(3, 1, dg)
    H1 = Md.transpose() * mat(Hc) * Mc
    H2 = Ms.transpose() * mat(Hc) * Md
    H3 = Ms.transpose() * mat(dH) * Ms
    return first + second, H1 + H2 + H3, dq, (H1, H2, H3)


def prepare():
    p = read(PREP / "protocol.json")
    delivery = read(PREP / "preparation_delivery.json")
    bindings = dict(p["source_and_input_hashes"])
    for name, digest in delivery["file_hashes"].items():
        file = PREP / name
        if sha(file) != digest:
            raise ValueError("producer delivery hash mismatch: " + name)
        bindings[file.relative_to(ROOT).as_posix()] = digest
    for file in (Path(__file__), PREP / "preparation_delivery.json"):
        bindings[file.relative_to(ROOT).as_posix()] = sha(file)
    for name, digest in bindings.items():
        if sha(ROOT / name) != digest:
            raise ValueError("changed source/input: " + name)
    protocol = {
        "status": "PREPARED_INDEPENDENT_LOCAL_REVIEW_ONLY",
        "utc": now(),
        "source_sha256": sha(Path(__file__)),
        "source_and_input_hashes": bindings,
        "precision_bits": 512,
        "actual_local_sites": [[0, 0, 0], [1, 11, 11]],
        "policy": "Independent analytic derivatives, symbolic rational and Torch AD checks, zero-crossing/near-cancellation/invalid-domain adversaries; 64 local nodal vertices times 4 shared material endpoints at each of two actual common-box sites plus center. No optimization, integration or whole24 assembly.",
        "scope": "Unweighted point formula and two supplied common nodal-box sites; no root membership or five-gate completion.",
    }
    write(OUT / "protocol.json", protocol)
    print(json.dumps({"status": protocol["status"], "bindings": len(bindings), "source_sha256": protocol["source_sha256"], "protocol_sha256": sha(OUT / "protocol.json")}))


def run():
    p = read(OUT / "protocol.json")
    if p["source_sha256"] != sha(Path(__file__)):
        raise ValueError("review source changed after freeze")
    for name, digest in p["source_and_input_hashes"].items():
        if sha(ROOT / name) != digest:
            raise ValueError("changed frozen source/input: " + name)
    ctx.prec = p["precision_bits"]
    start = time.monotonic()
    checks = {}
    diagnostics = {}

    def ck(name, value):
        if value is not True:
            raise AssertionError(name)
        checks[name] = True

    def contains(ball, val):
        return P.contains_exact(ball, val)

    def all_contains(outer, inner):
        return all(x.contains(y) for x, y in zip(outer.entries(), inner.entries()))

    # This direct helper has safe semantics despite its use of == rather than
    # is_zero: Arb ==0 is true only for exact zero, so not(==0) retains overlap.
    for exponent in (1, 80, 300):
        x = arb(0, P.core.upper(A(Q(1, 10**exponent))))
        ck(f"zero_crossing_{exponent}_retained", P.ZERO in P.jclean({P.ZERO: x}))
        ck(f"zero_crossing_{exponent}_eq_false", (x == 0) is False)
        near = P.jadd(P.jconstant(arb(1, P.core.upper(A(Q(1, 10**exponent))))), P.jconstant(-1))
        ck(f"near_cancelled_{exponent}_retained", P.ZERO in near and not near[P.ZERO].is_zero())
        multiplied = P.jmul(near, P.jvariable(1, 1))
        ck(f"near_cancelled_{exponent}_jet_product_retained", (0, 1, 0) in multiplied and not multiplied[(0, 1, 0)].is_zero())
    ck("exact_zero_jet_coefficient_removed", P.jconstant(0) == {})
    for label, ball in (("positive_infinity", arb("+inf")), ("indeterminate", arb("nan"))):
        try:
            P.jclean({P.ZERO: ball})
            ck("nonfinite_rejected_" + label, False)
        except ValueError:
            ck("nonfinite_rejected_" + label, True)

    import sympy as sy
    a, b, c, K, T = sy.symbols("a b c K T", real=True)
    xs = (a, b, c)
    W = K * a**2 * (b**2 + c**2) / 2 + (T - K) * a * sy.sqrt(b**2 + c**2) + K / 2
    points = [([Q(6, 5), Q(4, 5), Q(3, 5)], Q(3), Q(7)),
              ([Q(9, 10), Q(-3, 5), Q(4, 5)], Q(2), Q(1)),
              ([Q(5, 4), Q(0), Q(-7, 5)], Q(7, 3), Q(2, 3)),
              ([Q(1), Q(1, 100), Q(0)], Q(2), Q(3))]
    for ip, (vals, kk, tt) in enumerate(points):
        truth = full_derivatives(vals, kk, tt)
        actual = P.surface_derivatives(vals, kk, tt)
        substitutions = dict(zip(xs, map(sy.Rational, vals))) | {K: sy.Rational(kk), T: sy.Rational(tt)}
        for order in range(4):
            indices_list = [()] if order == 0 else product(range(3), repeat=order)
            for indices in indices_list:
                exact = Q(str(sy.diff(W, *[xs[i] for i in indices]).subs(substitutions))) if indices else Q(str(W.subs(substitutions)))
                got = truth[order]
                produced = actual[("value", "gradient", "Hessian", "third")[order]]
                for i in indices:
                    got, produced = got[i], produced[i]
                label = f"symbolic_rational_{ip}_{order}_{''.join(map(str, indices))}"
                ck(label, contains(got, exact) and contains(produced, exact))

    # Numeric AD is an additional diagnostic; the exact rational tests above
    # and the analytic FTC derivation supply the mathematical checks.
    import torch
    torch.set_default_dtype(torch.float64)
    rng = np.random.default_rng(4531)
    max_ad_relative = 0.0
    for ip in range(8):
        x = torch.tensor([rng.uniform(.5, 1.7), rng.uniform(.3, 1.4), rng.uniform(-.8, .8)], requires_grad=True)
        kk, tt = rng.uniform(.5, 4), rng.uniform(0, 5)
        def potential(z):
            aa, bb, cc = z
            return kk * aa**2 * (bb**2 + cc**2) / 2 + (tt - kk) * aa * torch.sqrt(bb**2 + cc**2) + kk / 2
        g = torch.autograd.functional.jacobian(potential, x)
        H = torch.autograd.functional.hessian(potential, x)
        D = torch.autograd.functional.jacobian(lambda z: torch.autograd.functional.hessian(potential, z, create_graph=True), x)
        own = full_derivatives([Q.from_float(v) for v in x.detach().numpy()], Q.from_float(float(kk)), Q.from_float(float(tt)))
        produced = P.surface_derivatives(x.detach().numpy(), float(kk), float(tt))
        for label, numeric, analytic, formal in (("g", g.detach().numpy(), own[1], produced["gradient"]), ("H", H.detach().numpy(), own[2], produced["Hessian"]), ("third", D.detach().numpy(), own[3], produced["third"])):
            for index in np.ndindex(numeric.shape):
                va, vf = analytic, formal
                for i in index:
                    va, vf = va[i], vf[i]
                error = max(abs(float(va.mid()) - numeric[index]), abs(float(vf.mid()) - numeric[index])) / max(1.0, abs(numeric[index]))
                max_ad_relative = max(max_ad_relative, error)
                ck(f"numeric_AD_{ip}_{label}_{'_'.join(map(str, index))}", bool(error < 1e-12))
    diagnostics["numeric_AD_max_relative_error"] = max_ad_relative

    # Fixed-coordinate adversary: selected columns omit the only map change.
    # An implementation that forms deltaq from selected/free Y alone fails.
    Bc = [[Q(1), 0, 0, 0, Q(1), 0], [0, 0, Q(1), 0, 0, 0], [0, 0, 0, Q(1), 0, 0]]
    Bs = [r.copy() for r in Bc]
    Bs[0][0] += Q(1, 2**42)
    Y = [Q(1, 2), 0, 1, 0, Q(1, 2), 0]
    selected = [2, 3, 4, 5]
    bound = P.consistency_bound(Bs, Bc, Y, selected, 2, 3)
    fullF, fullH = pullback(Bs, Y, 2, 3)
    oldF, oldH = pullback(Bc, Y, 2, 3)
    dF, dH = fullF - oldF, fullH - oldH
    ck("fixed_node_deltaq_nonzero", not bound["delta_q"][0].contains(0))
    ck("fixed_node_change_selected_deltaB_zero", all(v == 0 for row in bound["delta_B_full"] for j, v in enumerate(row) if j in selected))
    ck("fixed_node_effect_force_nonzero", any(not dF[j, 0].contains(0) for j in selected))
    ck("fixed_node_effect_H_nonzero", any(not dH[i, j].contains(0) for i in selected for j in selected))
    ck("fixed_node_fullY_force_enclosed", all(bound["force_difference"][i, 0].contains(dF[j, 0]) for i, j in enumerate(selected)))
    ck("fixed_node_fullY_H_enclosed", all(bound["H_difference"][i, j].contains(dH[u, v]) for i, u in enumerate(selected) for j, v in enumerate(selected)))

    # Exact normalized raw Js identity, independent of whether L_I equals L0c.
    for ir, (aa, bb, cc, L, alpha) in enumerate([(Q(6, 5), Q(4, 5), Q(3, 5), Q(7, 3), Q(11, 5)), (Q(2), Q(-5, 13), Q(12, 13), Q(1, 10), Q(3, 8))]):
        raw = [aa, bb, cc]
        normalized = [alpha * L * aa, bb / L, cc / L]
        ck(f"normalized_Js_squared_exact_{ir}", normalized[0]**2 * (normalized[1]**2 + normalized[2]**2) == alpha**2 * aa**2 * (bb**2 + cc**2))
        oldg, oldH, _ = P.core.surface_local(list(map(A, raw)), A(3), A(7), A(alpha))
        _, ng, nh, _ = full_derivatives(normalized, 3, 7)
        Ddiag = [alpha * L, 1 / L, 1 / L]
        for i in range(3):
            ck(f"normalized_raw_gradient_identity_{ir}_{i}", oldg[i].overlaps(A(Ddiag[i]) * ng[i]))
            for j in range(3):
                ck(f"normalized_raw_H_identity_{ir}_{i}_{j}", oldH[i][j].overlaps(A(Ddiag[i] * Ddiag[j]) * nh[i][j]))

    producer_protocol = read(PREP / "protocol.json")
    row = producer_protocol["selected_record"]
    boxfile = ROOT / "runs/exact_integration_bridge_v1/retube_B_geometry/coordinate_boxes.json"
    box = read(boxfile)
    ck("same_cell_commonbox_binding", box["parameter_cell_sha256"] == row["parameter_cell_sha256"] and sha(boxfile) == producer_protocol["new_common_geometry_box_sha256"])
    ck("commonbox_does_not_assume_root_membership", box["root_membership_assumed"] is False)
    target = P.target.Target(read(ROOT / row["point_cache_file"]))
    # Importing legacy target is safe ONLY along this inspected decode/rational
    # container path. A raising sentinel proves these local fixtures do not call
    # its invalid polynomial cleanup. No production source is modified.
    saved_clean = P.target.clean
    def reject_unsafe_cleanup(*args, **kwargs):
        raise AssertionError("legacy polynomial cleanup reached")
    P.target.clean = reject_unsafe_cleanup
    site_receipts = []
    try:
        for face, edge, iq in p["actual_local_sites"]:
            label = f"{face}_{edge}_{iq}"
            Bs, Bc, raw, D, nodes, selected, R, L = own_maps(target, face, edge, iq)
            pm = P.actual_maps(target, face, edge, iq)
            ck("own_rational_full_maps_" + label, Bs == pm["star_full"] and Bc == pm["cache_normalized_full"])
            ck("own_normalization_and_selected_columns_" + label, D == pm["normalization_D_exact"] and selected == pm["selected_local"])
            centers = [Q(box["current_centers_exact"][n][c]) for n in nodes for c in (0, 1)]
            radii = [Q(box["coordinate_radii_exact"][n][c]) for n in nodes for c in (0, 1)]
            Ybox = [interval(x - r, x + r) for x, r in zip(centers, radii)]
            prefix = "basal" if face == 0 else "apical"
            kl, kh = map(Q, row["parameter_cell"][prefix + "_area_modulus"])
            tl, th = map(Q, row["parameter_cell"][prefix + "_surface_tension"])
            Kbox, Tbox = interval(kl, kh), interval(tl, th)
            fresh = P.consistency_bound(Bs, Bc, Ybox, selected, Kbox, Tbox)
            independentF, independentH, independentdq, hterms = analytic_difference_bound(Bs, Bc, Ybox, selected, Kbox, Tbox)
            ck("explicit_analytic_F_enclosure_overlap_" + label, all(x.overlaps(y) for x, y in zip(independentF.entries(), fresh["force_difference"].entries())))
            ck("explicit_analytic_H_enclosure_overlap_" + label, all(x.overlaps(y) for x, y in zip(independentH.entries(), fresh["H_difference"].entries())))
            ck("coefficientwise_dq_matches_" + label, all(x.overlaps(y) for x, y in zip(independentdq, fresh["delta_q"])))
            for it, (ours, name) in enumerate(zip(hterms, ("H_deltaBT_Hc_Bc", "H_BstarT_Hc_deltaB", "H_BstarT_dH_Bstar"))):
                ck(f"each_Hessian_term_{label}_{it}", all(x.overlaps(y) for x, y in zip(ours.entries(), fresh[name].entries())))
            vertex_count = 0
            for signs in product((-1, 1), repeat=6):
                vertex = [x + sign * r for x, sign, r in zip(centers, signs, radii)]
                for kval, tval in product((kl, kh), (tl, th)):
                    Fs, Hs = pullback(Bs, vertex, kval, tval)
                    Fc, Hc = pullback(Bc, vertex, kval, tval)
                    FF, HH = Fs - Fc, Hs - Hc
                    ck(f"actual_vertex_force_{label}_{vertex_count}", all(fresh["force_difference"][i, 0].contains(FF[j, 0]) for i, j in enumerate(selected)))
                    ck(f"actual_vertex_H_{label}_{vertex_count}", all(fresh["H_difference"][i, j].contains(HH[u, v]) for i, u in enumerate(selected) for j, v in enumerate(selected)))
                    vertex_count += 1
            # A bit-for-bit producer replay at its own declared precision is
            # separate from independent analytic tests at 512 bits.
            ctx.prec = producer_protocol["precision_bits"]
            replay = P.actual_fixture(target, row, box, face, edge, iq)
            file = PREP / f"local_face{face}_edge{edge}_q{iq}.json"
            ck("fresh_producer_fixture_reproduction_" + label, P.serial(replay) == read(file))
            ctx.prec = p["precision_bits"]
            site = {
                "face": face, "edge": edge, "quadrature_index": iq,
                "scope": "Uniform LOCAL unweighted supplied common-box site, not integration/root.",
                "producer_fixture_sha256": sha(file),
                "nodes": nodes, "selected_local": selected,
                "canonical_R_exact": str(R), "canonical_L_exact": str(L),
                "normalization_D_exact": list(map(str, D)),
                "same_fullY_vertex_cases": vertex_count,
                "independent_force_difference": P.pack_matrix(independentF),
                "independent_H_difference": P.pack_matrix(independentH),
                "fresh_jet_force_difference": P.pack_matrix(fresh["force_difference"]),
                "fresh_jet_H_difference": P.pack_matrix(fresh["H_difference"]),
                "legacy_sparse_cleanup_calls": 0,
            }
            name = f"site_face{face}_edge{edge}_q{iq}.json"
            write(OUT / name, site)
            site_receipts.append({"file": name, "sha256": sha(OUT / name)})
            print(json.dumps({"site_complete": [face, edge, iq], "vertex_cases": vertex_count}), flush=True)
    finally:
        P.target.clean = saved_clean
    ck("two_actual_sites_completed", len(site_receipts) == 2)
    ck("legacy_sparse_cleanup_not_used_by_actual_local_path", True)
    for name, digest in p["source_and_input_hashes"].items():
        ck("binding_current_" + name, sha(ROOT / name) == digest)
    receipt = {
        "status": "PASS_INDEPENDENT_LOCAL_SURFACE_CACHE_CONSISTENCY_REVIEW",
        "utc": now(), "checks": checks,
        "source_sha256": sha(ROOT / "scripts/exact_surface_cache_consistency.py"),
        "review_source_sha256": sha(Path(__file__)),
        "review_protocol_sha256": sha(OUT / "protocol.json"),
        "producer_protocol_sha256": sha(PREP / "protocol.json"),
        "source_and_input_hashes": p["source_and_input_hashes"],
        "site_receipts": site_receipts,
        "diagnostics": diagnostics,
        "elapsed_seconds": time.monotonic() - start,
        "third_derivative_proof": "Finite degree3 nilpotent jet of analytic sqrt on positive norm box; explicit independent closed forms and exact rational symbolic derivatives checked. No spatial Taylor truncation remainder.",
        "formula_proof": "dq=(Bs-Bc)Y full six coordinates. Fundamental theorem of calculus bounds dg by H on qc+[0,1]dq and dH by Wthird there. Fdiff=deltaB.T gc+Bs.T dg. Hdiff=deltaB.T Hc Bc+Bs.T Hc deltaB+Bs.T dH Bs.",
        "safe_cleanup_scope": "Direct jclean uses not(v==0), retains every tested nonexact zero-crossing finite ball, rejects nonfinite via A. Legacy target import is limited here to rational/decode/container/bind; actual two-site path also passed with legacy clean replaced by a raising sentinel. No approval of legacy integration arithmetic.",
        "actual_uniform_local_sites": 2,
        "site_vertex_cases_each": 256,
        "all24_surfaces_assembled": False,
        "spatial_integration_run": False,
        "uniform_integral_or_root_bound_proved": False,
        "old_unsafe_point_certificate_restored": False,
        "original_five_gates_complete": False,
    }
    write(OUT / "results.json", receipt)
    write(ROOT / "evidence/exact_surface_cache_consistency_independent_review.json", receipt)
    print(json.dumps({"status": receipt["status"], "checks": len(checks), "receipt_sha256": sha(ROOT / "evidence/exact_surface_cache_consistency_independent_review.json"), "elapsed_seconds": receipt["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    if sys.argv[1:] == ["prepare"]:
        prepare()
    elif sys.argv[1:] == ["run"]:
        run()
    else:
        raise SystemExit("prepare | run")
