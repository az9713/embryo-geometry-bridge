from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[4]
RUN = Path("runs/transfer_ambiguity_nonlinear_probe_v1")
REP = RUN / "report"
OUTPUT = Path(__file__).with_name("receipt.json")
LINEAR = Path("runs/transfer_ambiguity_linear_v1/actual/results.json")
PINS = {
    REP / "report_build_receipt.json": "79cae6792ea8fbfef89bd9d69a3c87f14f419fb1a6c03384f781caf92eeff7c9",
    REP / "REPORT.html": "d33f3d594406ed498d364c3462450b62cb27054f0146f4bbca3a1c2ae22b41f6",
    REP / "browser_readback_v2/receipt.json": "c0236dd4a464cd29b2b4925a3232fb4017ca0ddc65952cefb839118005a52d8d",
    Path("scripts/report_transfer_ambiguity_nonlinear.py"): "4240680fce5745ead9e7d2539b6d04498183b56559da29dbe5dbf325d5568e22",
    RUN / "actual/results.json": "2ddc2fb3f0da426551bf74b6e2aac7bb489dfe848f8d9e73cfc42fddecbe99d6",
    RUN / "independent_actual_preparation/actual_replay/receipt.json": "c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd",
    RUN / "observer_robustness/actual/results.json": "53a17296ddae08b4d1fe6b7e9b20ad0a377bf8dc962142bbfb9eb9877a381647",
    RUN / "observer_robustness/independent_actual/receipt.json": "703b574bc51a083b10c7328d9622862e939ddabe791b18ce9e1817f07718f69f",
    RUN / "originality_followup/receipt.json": "6797ad5d05f71535f13edc519276e57bc5176d7d5b0e6273b1aa91eba48653d0",
    RUN / "field_manifest.json": "308ca4c05e476ba6a14cbc7b8d34018983b8f11263e2c224488672b35eb38ac0",
    RUN / "preparation/protocol.json": "9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94",
    RUN / "actual/terminal.json": "84aafea79a60878e139ed034c199e4875e91bac21fa650354f12f944e87fd0b9",
    LINEAR: "68ddc210523cb40001c99babb251916ac22c62eecd836d976fac1ab78dc0fe3f",
}
VARIANTS = (
    "current_disk__spatial_max", "current_disk__regional_mean",
    "current_radial_line__spatial_max", "current_radial_line__regional_mean",
)
LABELS = (
    "Disk · component maxima", "Disk · component means",
    "Radial line · component maxima", "Radial line · component means",
)
HASHES = {}


def bytes_of(path, digest=None):
    path = Path(path)
    data = (ROOT / path).read_bytes()
    observed = sha256(data).hexdigest()
    if digest is not None:
        assert observed == digest, (path, observed, digest)
    HASHES[path.as_posix()] = observed
    return data


def read(path, digest=None):
    path = Path(path)
    return json.loads(bytes_of(path, digest or PINS.get(path)))


def pair(values):
    lo, hi = map(F, values)
    assert lo <= hi
    return lo, hi


def minus(a, b):
    return a[0] - b[1], a[1] - b[0]


def midpoint(a):
    return sum(a) / 2


def max_pair(a):
    return max(v[0] for v in a), max(v[1] for v in a)


def bound(rows, observations, delta):
    lows = {}
    highs = {}
    assert len(rows) == 68
    seen = set()
    for row in rows:
        key = (row["record_id"], row["observable"])
        assert key not in seen
        seen.add(key)
        y = F(row["observed_exact"])
        assert y == observations[key]
        lo, hi = (x - y for x in pair(row["prediction_interval_exact"]))
        low = F(0) if lo <= 0 <= hi else min(abs(lo), abs(hi))
        high = max(abs(lo), abs(hi))
        assert low == F(row["minimum_interval_residual_exact"])
        assert high == F(row["maximum_interval_residual_exact"])
        assert delta == F(row["discrepancy_allowance_exact"])
        label = key[1]
        lows[label] = max(low, lows.get(label, 0))
        highs[label] = max(high, highs.get(label, 0))
    assert seen == observations.keys()
    return (max(lows.values()), max(highs.values())), lows, highs


def compare_case(case, observations, delta):
    a = case["all68_residuals"]
    envelope, low, high = bound(a["constraint_rows"], observations, delta)
    assert low == {k: F(v) for k, v in a["discrepancy_lower_requirements_exact"].items()}
    assert high == {k: F(v) for k, v in a["discrepancy_sufficient_requirements_exact"].items()}
    assert envelope == pair(case["required_allowance_um_exact"])
    return envelope


class Tags(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.links = []
        self.images = []
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "a" and "href" in attrs:
            self.links.append(attrs["href"])
        if tag == "img":
            self.images.append(attrs["src"])


def audit():
    for path, digest in PINS.items():
        bytes_of(path, digest)
    build = read(REP / "report_build_receipt.json")
    assert build["source_sha256"] == PINS[Path(build["source"])]
    assert len(build["input_hashes"]) == 19 and len(build["report_and_figure_hashes"]) == 6
    for bindings in (build["input_hashes"], build["report_and_figure_hashes"]):
        for path, digest in bindings.items():
            bytes_of(path, digest)

    env = read(RUN / "preparation/protocol.json")
    a = read(RUN / "actual/results.json")
    i = read(RUN / "independent_actual_preparation/actual_replay/receipt.json")
    s = read(RUN / "observer_robustness/actual/results.json")
    si = read(RUN / "observer_robustness/independent_actual/receipt.json")
    linear = read(LINEAR)
    manifest = read(RUN / "field_manifest.json")
    terminal = read(RUN / "actual/terminal.json")
    assert i["actual_results_sha256"] == PINS[RUN / "actual/results.json"]
    assert si["input_and_output_sha256"][(RUN / "observer_robustness/actual/results.json").as_posix()] == PINS[RUN / "observer_robustness/actual/results.json"]
    assert terminal["status"] == "FINISHED" and terminal["results_sha256"] == i["actual_results_sha256"]
    assert len(terminal["workers"]) == 2 and all(r["observed_exit_code"] == 0 for r in terminal["workers"])
    assert manifest["independent_replay_observed_os_exit_code"] == 0
    assert manifest["independent_primary_review_sha256"] == PINS[RUN / "independent_actual_preparation/actual_replay/receipt.json"]
    assert manifest["primary_terminal_sha256"] == PINS[RUN / "actual/terminal.json"]
    assert manifest["primary_results_sha256"] == i["actual_results_sha256"]

    selection = env["selection"]
    baseline = F(linear["source_baseline_residual_um_exact"])
    delta = baseline + F(1, 100)
    row, = [r for r in linear["rows"] if r["mode"] == selection["mode"] and F(r["budget_exact"]) == F(1, 100)]
    assert F(selection["affine_baseline_source_residual_um_exact"]) == baseline
    assert F(selection["conditional_source_allowance_um_exact"]) == F(row["allowance_um_exact"]) == delta
    assert F(selection["increment_um_exact"]) == F(1, 100)
    assert F(selection["scale_each_frozen_step_exact"]) == F(1, 2)
    assert F(s["frozen_primary_source_allowance_um_exact"]) == F(si["source_configured_delta_um_exact"]) == delta
    assert env["status"] == "FROZEN_OUTCOME_FREE_NONLINEAR_DIAGNOSTIC"
    gradient = list(map(F, linear["rounded_affine_contrast_gradient_exact"]))
    params = env["matched_input"]
    dose = F(params["matched_dose_exact"])
    diameter = F(params["large_reported_diameter_um_exact"])
    shoulder = F(env["large_observable"]["shoulder_radius_um_exact"])
    assert diameter == 2 * shoulder and dose == F("0.89767374")

    md = bytes_of(REP / "REPORT.md").decode()
    shifts = []
    obs = {}
    endpoint_ids = []
    md_primary_rows = []
    for symbol, case in zip(("+", "−"), env["endpoints"], strict=True):
        name = case["endpoint_id"]
        endpoint_ids.append(name)
        full = list(map(F, row[name]["step_exact"]))
        half = [t / 2 for t in full]
        assert list(map(F, case["full_linear_step_exact"])) == full
        assert list(map(F, case["half_step_exact"])) == half
        assert row[name]["exact_primal_feasible"] is True and row[name]["exact_weak_duality"] is True
        shift = sum(g * h for g, h in zip(gradient, half, strict=True))
        assert shift == F(case["affine_large_Q_change_um_exact"])
        shifts.append(shift)
        assert f"| {symbol} | " + " | ".join(f"{float(x):+.8f}" for x in half) + " |" in md

        gains = [F(case["materials"][mode]["gain_exact"]) for mode in ("small", "large")]
        assert gains[0] == gains[1]
        assert [F(r["dose_exact"]) for r in case["large_checkpoints"]] == [dose * F(k, 3) for k in range(4)]
        small = case["small_records"]
        doses = [F(r["dose_exact"]) for r in small]
        assert len(small) == 17 and doses == sorted(doses)
        for checkpoint in [case["small_zero"], *small, *case["large_checkpoints"]]:
            strength = F(checkpoint["dose_exact"]) * gains[0]
            rounded = F.from_float(float(strength))
            assert F(checkpoint["requested_strength_exact"]) == strength
            assert F(checkpoint["solver_strength_binary64_exact"]) == rounded
            assert abs(F(checkpoint["binary64_rounding_difference_exact"])) == abs(rounded - strength)
        obs[name] = {(r["record_id"], k): F(v) for r in small for k, v in r["observed_exact"].items()}
        assert len(obs[name]) == 68
    assert obs[endpoint_ids[0]] == obs[endpoint_ids[1]]
    affine_gap = shifts[0] - shifts[1]
    assert 2 * affine_gap == pair(row["contrast_range_width_um_exact"])[0]
    assert affine_gap == F(selection["half_pair_affine_large_Q_range_um_exact"]) == F(i["half_step_affine_Q_difference_um_exact"])

    q = {}
    primary = {}
    for symbol, e, ie in zip(("+", "−"), a["endpoints"], i["endpoints"], strict=True):
        name = e["endpoint_id"]
        assert name == ie["endpoint_id"]
        q[name] = pair(e["large_Q_um_exact"])
        assert q[name] == pair(ie["large_Q_um_exact"])
        primary[name] = compare_case(e["small_residual"], obs[name], delta)
        assert primary[name] == pair(ie["required_allowance_um_exact"])
        assert primary[name][1] < delta
        table_row = (f"| {symbol} | {float(midpoint(q[name])):+.11f} | "
                     f"{float(primary[name][1]):.11f} | {float(delta - primary[name][1]):.11f} |")
        assert table_row in md
        md_primary_rows.append(table_row)
    g = minus(*(q[name] for name in endpoint_ids))
    departure = minus(g, (affine_gap, affine_gap))
    assert 0 < g[0] and g == pair(a["large_Q_positive_direction_minus_negative_direction_um_exact"])
    assert g == pair(i["stored_large_Q_signed_difference_um_exact"])
    assert departure == pair(i["nonlinear_minus_affine_signed_difference_um_exact"])
    assert max_pair(list(primary.values())) == pair(a["both_endpoints_required_allowance_um_exact"])

    variants = {}
    md_secondary_rows = []
    for variant, label in zip(VARIANTS, LABELS, strict=True):
        envelopes = []
        for e in a["endpoints"]:
            name = e["endpoint_id"]
            case = s["by_endpoint_and_variant"][name][variant]
            if variant == VARIANTS[0]:
                assert case == e["small_residual"]
            envelopes.append(compare_case(case, obs[name], delta))
        common = max_pair(envelopes)
        sc = s["common_pair_requirement_by_variant"][variant]
        sic = si["by_observer_common_requirements"][variant]
        assert common == pair(sc["required_allowance_um_exact"])
        assert common == pair(sic["independent_required_allowance_um_exact"])
        assert common == pair(sic["required_allowance_um_exact"])
        if variant == VARIANTS[0]:
            assert common[1] < delta
            decision = "Both satisfy"
        else:
            assert common[0] > delta
            decision = "At least one exceeds"
        table_row = f"| {label} | {float(common[1]):.10f} | {float(common[1] - delta):+.10f} | {decision} |"
        assert table_row in md
        md_secondary_rows.append(table_row)
        variants[variant] = common
    all_four = max_pair(list(variants.values()))
    assert all_four == pair(s["common_pair_requirement_across_all_four_conventions"]["required_allowance_um_exact"])
    assert all_four == pair(si["all_four_required_common_allowance_um_exact"])
    assert s["primary_replay_exactly_matches"] is True and si["primary_exact_result_kept"] is True
    assert si["robust_across_conventions_at_frozen_delta"] is False
    assert s["additional_equilibria"] == s["additional_force_evaluations"] == 0 and s["retuning"] is False

    roles = Counter(f["role"] for f in manifest["fields"])
    assert roles == {"small_record": 34, "large_final": 2, "large_checkpoint": 4, "small_zero": 2, "large_zero": 2}
    assert manifest["status"] == "COMPLETE_44_PRIMARY_CHECKPOINTS"
    assert a["verified_outcome_fields"] == 36 and a["small_observer_intervals"] == i["independently_recomputed_small_residuals"] == 136
    assert s["measured_small_fields"] == roles["small_record"]
    assert s["small_constraint_evaluations"] == si["independent_interval_exact_matches"] == 544
    assert s["measured_small_fields"] * 4 * len(VARIANTS) == 544

    profiles = {}
    for f in manifest["fields"]:
        for kind in ("field", "verification"):
            path, digest = f[kind], f[f"{kind}_sha256"]
            assert i["all_hashes"][path] == digest
            bytes_of(path, digest)
        if f["role"] != "large_final":
            continue
        name = f["endpoint_id"]
        proof = read(f["verification"], f["verification_sha256"])["surface"]
        ps = {}
        for sample in proof["sampled_full_current_profiles"]:
            r = F(sample["radius_um_exact"])
            value = sample["apical_height_um"]
            assert r not in ps and value["lower_closed"] is True and value["upper_closed"] is True
            ps[r] = pair([F(value[k]["numerator"], value[k]["denominator"]) for k in ("lower", "upper")])
        assert F(proof["shoulder_radius_um_exact"]) == shoulder
        assert proof["usable_stored_graph"] is True
        assert minus(ps[0], ps[shoulder]) == q[name] == pair(proof["Q_um_exact"])
        profiles[name] = ps
    pp, pm = (profiles[name] for name in endpoint_ids)
    radii = sorted(r for r in pp.keys() & pm.keys() if r <= 54)
    assert radii[0] == 0 and shoulder in radii
    dz0, dza = (midpoint(pp[r]) - midpoint(pm[r]) for r in (0, shoulder))
    assert dz0 - dza == midpoint(g)
    plot_rows = {
        "radius_um_exact": list(map(str, radii)),
        "apical_height_plus_midpoints_um_exact": [str(midpoint(pp[r])) for r in radii],
        "apical_height_minus_midpoints_um_exact": [str(midpoint(pm[r])) for r in radii],
        "height_difference_midpoints_um_exact": [str(midpoint(pp[r]) - midpoint(pm[r])) for r in radii],
    }

    checkpoints = [cp for e in i["endpoints"] for cp in e["checkpoints"]]
    assert len(checkpoints) == i["checkpoints"] == 44
    assert i["numerical_force_evaluations"] == 2 * len(checkpoints) == 88
    assert i["stored_exact_element_checks"] == 144 * len(checkpoints) == 6336
    q8 = max(cp["Q8"] for cp in checkpoints)
    q12 = max(cp["Q12"] for cp in checkpoints)
    j = min(F(cp["J_lower_exact"]) for cp in checkpoints)
    assert 0 < j and max(q8, q12) < 1e-7
    assert all(a[k] is False for k in ("source_witness_certified", "exact_integral_FE_root_certified", "continuum_accuracy_certified", "heldout_response_tables_read"))
    assert i["exact_equilibrium_position_bounded"] is False
    margin = delta - max(v[1] for v in primary.values())
    displayed = [
        f"{float(delta):.11f}", f"{float(baseline):.11f}", f"{float(margin):.11f}",
        f"{float(shoulder):.8f}", f"{float(diameter):.7f}", f"{float(dose):.8f}",
        f"{float(midpoint(g)):.12f}", f"{float(affine_gap):.12f}", f"{float(midpoint(departure)):.12f}",
        f"{float(all_four[1]):.10f}", f"{float(j):.6f}", f"{q8:.3e}", f"{q12:.3e}", f"{float(g[1]-g[0]):.2e}",
    ]
    assert all(label in md for label in displayed)
    assert f"{float(1000 * margin):.8f}" == "4.81253906"

    novelty = read(RUN / "originality_followup/receipt.json")
    assert novelty["novelty_confirmed"] is False
    nrefs = {r["citation"]: r for r in novelty["sources"]}
    cited = build["primary_source_links"]
    assert nrefs["Liu, Maini and Baker"]["doi"] in cited["Liu, Maini and Baker (2026)"]
    assert nrefs["Kumar et al."]["url"] == cited["Kumar et al. (2024)"]
    assert nrefs["Countryman et al."]["doi"].split("/")[-1] in cited["Countryman et al. (2025)"]
    assert nrefs["Foucart and Liao"]["url"] == cited["Foucart and Liao"]
    assert nrefs["DeVore, Foucart, Petrova and Wojtaszczyk"]["url"] == cited["DeVore et al."]
    assert "A bounded literature search cannot establish priority" in md
    assert "conditional model discrepancy threshold" in md and "Experimental noise" in md

    html = bytes_of(REP / "REPORT.html", PINS[REP / "REPORT.html"])
    tags = Tags(html.decode())
    linked = [href for href in tags.links if not urlsplit(href).scheme]
    assert len(linked) == 18
    assert len(tags.images) == 3
    for link in linked + tags.images:
        path = (ROOT / REP / unquote(urlsplit(link).path)).resolve()
        path.relative_to(ROOT)
        assert path.is_file(), path
        bytes_of(path.relative_to(ROOT))
    assert set(tags.links) - set(linked) == set(cited.values()) - {cited["Liu et al. author preprint"]}

    browser = read(REP / "browser_readback_v2/receipt.json")
    assert browser["status"] == "PASS_BROWSER_READBACK" and browser["protected_documents_unchanged"] is True
    for path, digest in browser["document_hashes"].items():
        bytes_of(ROOT.parent / path, digest)
    assert [c["width"] for c in browser["checks"]] == [1440, 390]
    for check in browser["checks"]:
        assert check["status"] == "PASS" and check["width"] == check["document_width"]
        assert check["math"] == check["rendered"] == 34
        assert not any(check[k] for k in ("page_errors", "failed_requests", "math_errors", "raw_tex"))
        assert all(i["complete"] and i["width"] > 0 for i in check["images"])
        assert [i["src"] for i in check["images"]] == tags.images
        assert [x["href"] for x in check["local_links"]] == linked
        assert all(x["status"] == 200 for x in check["local_links"])
        for screenshot in check["screenshots"]:
            bytes_of(REP / "browser_readback_v2" / screenshot["file"], screenshot["sha256"])
    rev = build["revision"]
    archive = read(rev["preservation_manifest"], rev["preservation_manifest_sha256"])
    for path, data in archive["files"].items():
        old = bytes_of(data["archive"], data["sha256"])
        if path == (REP / "REPORT.html").as_posix() or path == build["source"]:
            assert old.count(b"table{border-collapse:") == 1
            assert old.replace(b"table{border-collapse:", b"table{position:relative;border-collapse:", 1) == bytes_of(path)

    return {
        "status": "PASS_INDEPENDENT_SCIENTIFIC_REPORT_READBACK",
        "audit_source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": [
            "Exact endpoint, budget, selection, contrast and threshold arithmetic from retained inputs and reviewed stored-field intervals.",
            "Existing independent primary and secondary model evaluations are premises. This readback performs no force evaluation, solve, integration, or observer recomputation from nodes.",
            "Sampled profile midpoint reconstruction and audited plot source; the connected lines carry no continuous error or root enclosure.",
            "Previously observed Edge screenshots and receipts were rehashed at both widths. Local HTML targets also checked directly.",
            "Literature credit checked against retained source audit. External links were not fetched again and originality is unestablished.",
        ],
        "counts": {"final_fields": 36, "all_checkpoints": 44, "primary_comparisons": 136,
                   "secondary_comparisons_including_primary": 544, "existing_force_replays": 88,
                   "existing_geometry_checks": 6336, "browser_widths": [1440, 390],
                   "browser_math_nodes_each": 34, "browser_images_each": 3, "browser_local_links_each": 18},
        "primary_allowance_um_exact": str(delta),
        "minimum_primary_margin_um_exact": str(margin),
        "minimum_primary_margin_nm_exact": str(1000 * margin),
        "signed_large_gap_um_exact": list(map(str, g)),
        "affine_half_step_gap_um_exact": str(affine_gap),
        "nonlinear_minus_affine_um_exact": list(map(str, departure)),
        "primary_table_rows": md_primary_rows,
        "secondary_table_rows": md_secondary_rows,
        "observer_pair_allowances_um_exact": {k: list(map(str, v)) for k, v in variants.items()},
        "profile_chart_samples": plot_rows,
        "rehashed_files": HASHES,
    }


if __name__ == "__main__":
    receipt = audit()
    with OUTPUT.open("x", encoding="utf-8") as f:
        json.dump(receipt, f, indent=2, ensure_ascii=False)
        f.write("\n")
