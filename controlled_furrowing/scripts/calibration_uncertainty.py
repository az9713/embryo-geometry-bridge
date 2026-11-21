"""Exact small-only shape-constrained calibration compatibility audit.

No spatial solver or larger-outcome file is read. Exact statements concern the
published decimal numbers, not precision of the biological measurements.
"""
from pathlib import Path
from decimal import Decimal, getcontext
from fractions import Fraction as F
import csv
import hashlib
import json
import platform
import xml.etree.ElementTree as ET
import zipfile

import fitz
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy.optimize import linprog

getcontext().prec = 60
ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs/calibration_uncertainty"
OBS = ["maximum_shortening_um", "region_averaged_shortening_um",
       "apical_max_shortening_um", "basal_max_shortening_um"]
LOCATORS = ["maximum_source", "average_source", "apical_source", "basal_source"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rat(q):
    q = F(q)
    return {"numerator": q.numerator, "denominator": q.denominator,
            "decimal_approximation": str(Decimal(q.numerator) / Decimal(q.denominator))}


def decimal_pairs(x, y):
    """Direct all-pairs Decimal calculation; independent of the scan below."""
    drops = [(y[i] - y[j], i, j) for i in range(len(x)) for j in range(len(x)) if x[i] <= x[j]]
    largest = max(Decimal(0), max(t[0] for t in drops))
    return largest / 2, [(i, j) for d, i, j in drops if d == largest]


def fraction_group_scan(x, y):
    """Sorted group-prefix scan, handles ties without imposing record order."""
    groups = {}
    for a, b in zip(x, y):
        groups.setdefault(a, []).append(b)
    high = None
    drop = F(0)
    for a in sorted(groups):
        high = max(groups[a]) if high is None else max(high, max(groups[a]))
        drop = max(drop, high - min(groups[a]))
    return drop / 2


def envelopes(x, y, delta):
    lower = [max(y[j] - delta for j in range(len(x)) if x[j] <= a) for a in x]
    upper = [min(y[j] + delta for j in range(len(x)) if x[j] >= a) for a in x]
    return lower, upper


def feasible_greedy(x, y, delta):
    """Independent interval-feasibility algorithm; no pair-drop formula."""
    last = None
    for a in sorted(set(x)):
        group = [b for xx, b in zip(x, y) if xx == a]
        lo, hi = max(group) - delta, min(group) + delta
        last = lo if last is None else max(last, lo)
        if last > hi:
            return False
    return True


def verify_curve(x, y, f, delta):
    assert all(abs(a - b) <= delta for a, b in zip(y, f))
    assert all(f[i] <= f[j] for i in range(len(x)) for j in range(len(x)) if x[i] <= x[j])


def source_decimal_check(rows):
    """Read only Fig2c/d/e XML worksheets; do not inspect larger outcomes."""
    path = ROOT / "data/raw/41467_2025_62483_MOESM16_ESM.xlsx"
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    relns = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    with zipfile.ZipFile(path) as archive:
        relations = {a.get("Id"): a.get("Target") for a in ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))}
        maps = {}
        for sheet in ET.fromstring(archive.read("xl/workbook.xml")).findall("m:sheets/m:sheet", ns):
            if sheet.get("name") in ["Fig 2c", "Fig 2d", "Fig 2e"]:
                target = relations[sheet.get("{" + relns + "}id")]
                target = target.lstrip("/") if target.startswith("/") else "xl/" + target
                cells = {}
                for cell in ET.fromstring(archive.read(target)).findall(".//m:c", ns):
                    value = cell.find("m:v", ns)
                    if value is not None and cell.get("t") != "s":
                        cells[cell.get("r")] = value.text
                maps[sheet.get("name")] = cells
    checks = []
    raw_rows = []
    for row in rows:
        raw_row = {"record_id": row["record_id"]}
        for key, locator in [("activation_magnitude", "dose_source")] + list(zip(OBS, LOCATORS)):
            sheet, cell = row[locator].split("!")
            original = maps[sheet][cell]
            # CSV uses canonical shortest float representation, while OOXML
            # includes decimal expansions of binary numbers. Retain both.
            assert float(row[key]) == float(original)
            assert abs(Decimal(row[key]) - Decimal(original)) <= Decimal("2e-15")
            raw_row[key] = original
            checks.append({"record": row["record_id"], "observable": key,
                           "locator": row[locator], "canonical_csv_decimal": row[key],
                           "exact_OOXML_numeric_token": original,
                           "OOXML_minus_CSV_decimal": str(Decimal(original) - Decimal(row[key]))})
        raw_rows.append(raw_row)
    assert len(checks) == 85
    return checks, raw_rows


def synthetic_checks():
    cases = [([0, 1, 2], [1, 2, 3], "0"),
             ([0, 1, 2], [3, 2, 1], "1"),
             ([0, 0, 1], [0, 4, 5], "2"),
             ([0, 1, 2, 3], [-1, -5, -2, -8], "3.5"),
             ([0], [7], "0")]
    for xx, yy, expected in cases:
        xd, yd = list(map(lambda a: Decimal(str(a)), xx)), list(map(lambda a: Decimal(str(a)), yy))
        d, _ = decimal_pairs(xd, yd)
        assert d == Decimal(expected)
        x, y = list(map(F, xx)), list(map(F, yy))
        assert fraction_group_scan(x, y) == F(d)
        assert feasible_greedy(x, y, F(d))
        if d > 0:
            assert not feasible_greedy(x, y, F(d) - F(1, 10**30))
        low, high = envelopes(x, y, F(d))
        verify_curve(x, y, low, F(d)); verify_curve(x, y, high, F(d))
    return {"known_cases": len(cases), "coverage": ["monotone", "descending", "dose ties", "negative observations", "singleton"], "pass": True}


def concave_lp_matrix(x, y):
    n = len(x)
    A, b, labels = [], [], []
    for i, yy in enumerate(y):
        for sign in [1, -1]:
            row = [F(0)] * (n + 1); row[i] = sign; row[-1] = -1
            A.append(row); b.append(sign * yy); labels.append(["observation", i, sign])
    for i in range(n):
        row = [F(0)] * (n + 1); row[i] = -1
        if i:
            row[i - 1] = 1
        A.append(row); b.append(F(0)); labels.append(["monotonicity", i])
    grid = [F(0)] + x
    for j in range(1, n):
        d1, d2 = grid[j] - grid[j - 1], grid[j + 1] - grid[j]
        row = [F(0)] * (n + 1)
        if j > 1:
            row[j - 2] = 1 / d1
        row[j - 1], row[j] = -1 / d1 - 1 / d2, 1 / d2
        A.append(row); b.append(F(0)); labels.append(["concavity", j])
    row = [F(0)] * (n + 1); row[-1] = -1
    A.append(row); b.append(F(0)); labels.append(["nonnegative_delta"])
    return A, b, labels


def main():
    RUN.mkdir(parents=True, exist_ok=True)
    data = ROOT / "data/derived/calibration_small_18um.csv"
    assert sha(data) == "1aa53b6e6ce0973faf0b038ae8d4059059b05c02bb87ed9698e24d1dea7b1518"
    source_paths = [data, ROOT / "data/raw/primary.xml",
                    ROOT / "data/raw/41467_2025_62483_MOESM16_ESM.xlsx",
                    ROOT / "data/raw/41467_2025_62483_MOESM15_ESM.pdf",
                    ROOT / "data/raw/41467_2025_62483_MOESM1_ESM.pdf"]
    before = {p.relative_to(ROOT).as_posix(): sha(p) for p in source_paths}
    with data.open(encoding="utf8", newline="") as handle:
        rows = sorted(csv.DictReader(handle), key=lambda r: Decimal(r["activation_magnitude"]))
    assert len(rows) == 17 and len(set(r["record_id"] for r in rows)) == 17
    source_checks, raw_rows = source_decimal_check(rows)
    xd = [Decimal(r["activation_magnitude"]) for r in rows]
    x = list(map(F, xd))
    results = {}
    for key in OBS:
        yd = [Decimal(r[key]) for r in rows]; y = list(map(F, yd))
        delta, pairs = decimal_pairs(xd, yd)
        raw_delta, _ = decimal_pairs([Decimal(r["activation_magnitude"]) for r in raw_rows],
                                     [Decimal(r[key]) for r in raw_rows])
        assert F(delta) == fraction_group_scan(x, y)
        assert feasible_greedy(x, y, F(delta))
        assert not feasible_greedy(x, y, F(delta) - F(1, 10**30))
        contracts = []
        for extra in [F(0), F(1, 2), F(1), F(delta)]:
            dd = F(delta) + extra; lower, upper = envelopes(x, y, dd)
            midpoint = [(a + b) / 2 for a, b in zip(lower, upper)]
            for curve in [lower, upper, midpoint]:
                verify_curve(x, y, curve, dd)
            table = [{"record_id": r["record_id"], "dose": r["activation_magnitude"],
                      "observed_um": r[key], "lower_um_exact": str(lo), "upper_um_exact": str(hi),
                      "midpoint_um_exact": str(mid)}
                     for r, lo, hi, mid in zip(rows, lower, upper, midpoint)]
            contracts.append({"delta": rat(dd), "extra_over_floor": rat(extra), "at_observed_doses": table})
        results[key] = {"monotonicity": "declared conditional common-curve assumption; not a theorem for the full spatial FE maximum",
                        "delta_star_decimal_exact": str(delta), "delta_star": rat(F(delta)),
                        "OOXML_token_delta_star_decimal_exact": str(raw_delta),
                        "witness_pairs": [{"lower_dose_record": rows[i]["record_id"],
                                           "higher_dose_record": rows[j]["record_id"],
                                           "doses_exact": [str(xd[i]), str(xd[j])],
                                           "responses_um_exact": [str(yd[i]), str(yd[j])],
                                           "response_drop_um_exact": str(yd[i] - yd[j])}
                                          for i, j in pairs], "contracts": contracts,
                        "below_floor_exact_rejection_test": True,
                        "pair_scan_interval_algorithms_agree_exactly": True}

    y = [F(r[OBS[0]]) for r in rows]
    # Sharper shape class: nonnegative increasing concave f with f(0)=0.
    j = next(i for i, r in enumerate(rows) if r["record_id"] == "optoGEF2_small_col_I")
    k = next(i for i, r in enumerate(rows) if r["record_id"] == "optoGEF2_small_col_M")
    beta = x[j] / x[k]
    dcon = (beta * y[k] - y[j]) / (1 + beta)
    gain = (y[k] - dcon) / x[k]
    cap = max(yy - dcon for yy in y)
    curve = [min(gain * xx, cap) for xx in x]
    A, b, labels = concave_lp_matrix(x, y)
    z = curve + [dcon]
    slacks = [rhs - sum(a * zz for a, zz in zip(row, z)) for row, rhs in zip(A, b)]
    assert min(slacks) >= 0 and max(abs(yy - ff) for yy, ff in zip(y, curve)) == dcon
    # Exact dual certificate: observation upper at j, lower at k, and the
    # origin-concavity consequence -f_j + beta*f_k <= 0.
    weights = [1 / (1 + beta), beta / (1 + beta), 1 / (1 + beta)]
    dual_rows = [[F(0)] * 18 for _ in range(3)]
    dual_rows[0][j] = 1; dual_rows[0][-1] = -1
    dual_rows[1][k] = -1; dual_rows[1][-1] = -1
    dual_rows[2][j] = -1; dual_rows[2][k] = beta
    combined = [sum(w * row[i] for w, row in zip(weights, dual_rows)) for i in range(18)]
    assert combined == [F(0)] * 17 + [F(-1)]
    assert weights[0] * y[j] - weights[1] * y[k] == -dcon
    objective = np.zeros(18); objective[-1] = 1
    lp = linprog(objective, A_ub=np.asarray(A, float), b_ub=np.asarray(b, float),
                 bounds=[(None, None)] * 18, method="highs")
    assert lp.success and abs(lp.fun - float(dcon)) < 1e-9
    concave = {"class": "nonnegative nondecreasing concave f on x>=0 with f(0)=0",
               "delta_star": rat(dcon), "witness_records": [rows[j]["record_id"], rows[k]["record_id"]],
               "beta": rat(beta), "dual_weights": list(map(rat, weights)),
               "exact_dual_combined_lhs": "-delta", "exact_dual_combined_rhs": rat(-dcon),
               "matching_curve": {"formula": "f(x)=min(g*x,C)", "g": rat(gain), "C_um": rat(cap),
                                  "saturation_dose": rat(cap / gain)},
               "exact_primal_constraint_count": len(slacks), "all_exact_slacks_nonnegative": True,
               "max_exact_residual": rat(max(abs(yy - ff) for yy, ff in zip(y, curve))),
               "numeric_independent_lp_delta": lp.fun, "numeric_lp_is_not_exact_proof": True,
               "restriction_warning": "Sharp for relaxed shape class; only necessary for the restricted column family, not full spatial FE"}
    raw_x = [F(r["activation_magnitude"]) for r in raw_rows]
    raw_y = [F(r[OBS[0]]) for r in raw_rows]
    raw_beta = raw_x[j] / raw_x[k]
    raw_dcon = (raw_beta * raw_y[k] - raw_y[j]) / (1 + raw_beta)
    raw_gain = (raw_y[k] - raw_dcon) / raw_x[k]
    raw_cap = max(yy - raw_dcon for yy in raw_y)
    raw_curve = [min(raw_gain * xx, raw_cap) for xx in raw_x]
    raw_A, raw_b, _ = concave_lp_matrix(raw_x, raw_y)
    assert all(rhs - sum(a * zz for a, zz in zip(row, raw_curve + [raw_dcon])) >= 0
               for row, rhs in zip(raw_A, raw_b))
    assert max(abs(yy - ff) for yy, ff in zip(raw_y, raw_curve)) == raw_dcon
    concave["OOXML_token_variant"] = {"delta_star": rat(raw_dcon), "matching_g": rat(raw_gain),
                                      "matching_C_um": rat(raw_cap), "all_68_exact_constraints_verified": True}

    # Optional homogeneous linear null: a concrete stricter calibration contract.
    dlin = max(F(0), max((y[i] * x[jj] - y[jj] * x[i]) / (x[i] + x[jj])
                         for i in range(17) for jj in range(17)))
    gl = max(F(0), max((yy - dlin) / xx for xx, yy in zip(x, y)))
    gu = min((yy + dlin) / xx for xx, yy in zip(x, y))
    assert gl == gu and max(abs(gl * xx - yy) for xx, yy in zip(x, y)) == dlin
    linear = {"class": "f(I)=g*I, g>=0; a null comparator, not the FE active-stress gain",
              "delta_star": rat(dlin), "attaining_gain_um_per_normalized_dose": rat(gl),
              "gain_contract": "g in [max(0,max_i((y_i-delta)/I_i)), min_i((y_i+delta)/I_i)]"}

    # Independent symbolic differentiation using d/dS=(x^2/h)d/dx.
    xx, hh, aa = sp.symbols("x h a", positive=True)
    nu = sp.symbols("nu", real=True)
    inverse = ((1 - 2 * nu) * xx ** sp.Rational(3, 2)
               + nu * xx ** sp.Rational(5, 2) - (1 - nu)) / aa
    first = sp.diff(inverse, xx) * xx**2 / hh
    second = sp.diff(first, xx) * xx**2 / hh
    expected_first = xx ** sp.Rational(5, 2) * (3 * (1 - 2 * nu) + 5 * nu * xx) / (2 * aa * hh)
    expected_second = 5 * xx ** sp.Rational(7, 2) * (3 * (1 - 2 * nu) + 7 * nu * xx) / (4 * aa * hh**2)
    assert sp.simplify(first - expected_first) == 0
    assert sp.simplify(second - expected_second) == 0

    xmlroot = ET.parse(ROOT / "data/raw/primary.xml").getroot()
    texts = {s.get("id"): " ".join(s.itertext()) for s in xmlroot.iter("sec")
             if s.get("id") in ["Sec14", "Sec15", "Sec16", "Sec17", "Sec22"]}
    assert "All measurements were taken from distinct samples." in texts["Sec22"]
    with fitz.open(ROOT / "data/raw/41467_2025_62483_MOESM15_ESM.pdf") as pdf:
        report_text = pdf[1].get_text()
    assert "Individual embryos were randomly selected." in report_text
    assert before == {p.relative_to(ROOT).as_posix(): sha(p) for p in source_paths}

    plt.style.use("dark_background")
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for axis, key, title in zip(axes.flat, OBS, ["Maximum shortening", "Region-averaged shortening",
                                               "Separate apical maximum (conditional)", "Separate basal maximum (conditional)"]):
        yy = np.asarray([float(r[key]) for r in rows]); dd = F(results[key]["delta_star_decimal_exact"])
        lower, upper = envelopes(x, list(map(F, [r[key] for r in rows])), dd)
        axis.fill_between(list(map(float, x)), list(map(float, lower)), list(map(float, upper)),
                          color="#4ea6b7", alpha=.35, label="Sharp bounds at observed doses")
        axis.plot(list(map(float, x)), [(float(a) + float(b)) / 2 for a, b in zip(lower, upper)], color="#76d4e9", lw=1)
        axis.scatter(list(map(float, x)), yy, s=24, color="#e8e8e8", zorder=3, label="Published small records")
        if key == OBS[0]:
            xx = np.linspace(0, float(max(x)), 300)
            axis.plot(xx, np.minimum(float(gain) * xx, float(cap)), color="#f5bc66", ls="--",
                      label=f"Concave fit (larger residual {float(dcon):.3f} µm)")
        axis.set_title(f"{title}\nmonotone floor = {float(dd):.9g} µm", fontsize=10)
        axis.set_ylabel("Shortening (µm)"); axis.grid(alpha=.15)
    for axis in axes[1]:
        axis.set_xlabel("Normalized activation magnitude")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Small-only deterministic compatibility: not confidence intervals\n"
                 "Bounds evaluated at observed doses; connectors interpolate the display grid", fontsize=11)
    fig.tight_layout(); fig.savefig(RUN / "compatibility_envelopes.png", dpi=170); plt.close(fig)

    result = {"status": "PASS exact deterministic bounds and independent checks; not a statistical noise estimate",
              "scope": "17 small-only records; no larger outcomes or spatial solver read",
              "script_sha256": sha(Path(__file__)), "software": {"python": platform.python_version(),
                    "numpy": np.__version__}, "source_hashes_before_after_equal": before,
              "published_decimal_workbook_crosschecks": source_checks,
              "synthetic_checks": synthetic_checks(), "monotone_observables": results,
              "concave_anchored_maximum_shortening": concave, "homogeneous_linear_maximum_shortening": linear,
              "column_inverse_derivatives": {"symbolic_first_derivative_verified": True,
                   "symbolic_second_derivative_verified": True,
                   "first_bracket_lower_bound_for_x_ge1_nu_0to1": "3-nu >=2",
                   "second_bracket_lower_bound_for_x_ge1_nu_0to1": "3+nu >=3",
                   "property": "Restricted column law increasing and strictly concave; not generalized to spatial FE"},
              "nested_set_peak_theorem": {"formula": "M(I;theta)=sup_{0<=A<=gI} O(q(A;theta))",
                   "proof": "For I1<=I2 and g>=0 the first amplitude set is contained in the second",
                   "common_branch_scaled_cap_full_path_required": True,
                   "instantaneous_FE_observable_need_not_be_monotone": True,
                   "source_spatial_temporal_maximum": "Methods Sec16",
                   "actual_experimental_path_equivalence_proven": False,
                   "restrictions": ["fixed measured doses", "common material/geometry/reference law", "common branch", "full nested path coverage", "same observation operator"]},
              "dose_order_caveat": {"maximum_witness_dose_gap": "0.00731338",
                   "uniform_dose_error_allowance_preserving_witness_order_inclusive": "0.00365669",
                   "larger_error": "Can reverse this witness order; does not prove all other compatibility bounds disappear"},
              "primary_unit_audit": {"Sec22_distinct_samples_explicit": True,
                   "Sec22_OptoGEF2_17_embryos_explicit": True,
                   "Fig2_literal_each_marker_is_individual_embryo_sentence": "not found in inspected caption",
                   "reporting_summary_page2_random_selection_and_blinding": True,
                   "missing_specimen_ID_is_not_evidence_of_dependence": True,
                   "IID_noise_or_independent_measurement_errors_proven": False,
                   "shared_dose_normalization": "Sec17 reference curve averaged across experiments; common calibration may couple errors",
                   "per_record_quantitative_error_bounds": "not found in inspected methods, worksheets, supplement captions and reporting summary",
                   "sampling_resolution_is_not_segmentation_uncertainty": True,
                   "locator_text_sha256": {k: hashlib.sha256(t.encode()).hexdigest() for k, t in texts.items()},
                   "reporting_page2_extracted_text_sha256": hashlib.sha256(report_text.encode()).hexdigest()},
              "interpretation": {"not_novel": "Standard order/shape-constrained deterministic inference",
                    "not_confidence_interval": True, "not_measurement_noise_estimate": True,
                    "not_FE_monotonicity_or_concavity_theorem": True,
                    "full_physical_identified_set_not_computed": True}}
    (ROOT / "evidence/calibration_uncertainty_bounds.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf8")
    print(json.dumps({"status": result["status"], "monotone_floor_um": results[OBS[0]]["delta_star_decimal_exact"],
                      "concave_floor_um": float(dcon), "linear_floor_um": float(dlin),
                      "workbook_checks": len(source_checks), "concave_exact_constraints": len(slacks)}))


if __name__ == "__main__":
    main()
