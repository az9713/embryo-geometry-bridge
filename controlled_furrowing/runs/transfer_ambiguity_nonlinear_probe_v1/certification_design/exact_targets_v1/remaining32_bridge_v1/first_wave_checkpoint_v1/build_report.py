"""Build the local first-wave checkpoint from exact, hash-pinned receipts."""
import csv
from collections import Counter
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from hashlib import sha256
from html import escape
from importlib.util import module_from_spec, spec_from_file_location
import json
import os
import re
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, PercentFormatter

HERE = Path(__file__).resolve().parent
BR = HERE.parent
EX = BR.parent
DESIGN = EX.parent
RUN = DESIGN.parent
ROOT = RUN.parents[1]
SITE = ROOT.parent
RED = BR / "first_wave_reduction_v1/independent_actual_preparation/actual_replay"
OBS = BR / "first_wave_observer_v1/independent_readback/actual"
GEOM = EX / "remaining32_geometry_v1/independent_actual_preparation/actual_replay"
HULL = BR / "independent_component_replay_preparation/actual_replay"
PRIOR = EX / "bridge_reduction_v1/independent_actual_preparation/all_four_readback.json"
PRIOR_PAGE = EX / "four_spatial_checkpoint_v1/SPATIAL.html"
CACHED = DESIGN / "full_outcome_coverage_v1/observer_application_v1/independent_readback/receipt.json"
CACHED_PAGE = DESIGN / "full_outcome_coverage_v1/observer_checkpoint_v1/OBSERVERS.html"
POLICY = BR / "admission/INDEPENDENT-HULL-IMPLEMENTATION-POLICY.md"
HELPER = EX / "certification_checkpoint_v1/build_checkpoint.py"
CONVERT = ROOT / "scripts/build_report.py"
ASSETS = SITE / "docs/assets/katex"
FAVICON = RUN / "report/favicon.png"
TITLE = "Exact spatial mechanics in the first new wave"
FIGURE = HERE / "first_wave_bounds.png"
LABELS = {
    "max_contrast_change_small_01": "+B",
    "max_contrast_change_small_13": "+E",
    "max_negative_contrast_change_small_01": "−B",
    "max_negative_contrast_change_small_16": "−K",
}
VARIANTS = {
    "current_disk__spatial_max": "Disk, spatial maximum (primary)",
    "current_disk__regional_mean": "Disk, regional mean",
    "current_radial_line__spatial_max": "Radial line, spatial maximum",
    "current_radial_line__regional_mean": "Radial line, regional mean",
}
PINS = {}
STATIC = {
    PRIOR: "859f507ad3762cfef44d53672bbe15a43592ce54d6e08ba8828bf50a925ec340",
    PRIOR_PAGE: "6dd9100643191f6617f5032d2a9d4e35ef3e179c9826115c841b041ce92c9e4d",
    CACHED: "a06549f048abba4eb2b4d1c1fb398f06e8f524196c92a673bd0493c527c8f40d",
    CACHED_PAGE: "dc7e2aa40c2511944f4ba0e226da1da844d866c85cf436c165ffa0e6a78ae557",
    HELPER: "50524d01db5003b0839e873efb37d41bb0552ad0f0c162f73bd668ad0e58d31b",
    CONVERT: "90148d53d74d8ea5cc8a665e6a3660abf2db72487a1fa59586caba98f9cc38e2",
    ASSETS / "katex.min.css": "5a0ff74f11887f7260dc16b216906f51352b84ca29f21abb9a14b62a2ca6eb15",
    ASSETS / "katex.min.js": "e6bfe5deebd4c7ccd272055bab63bd3ab2c73b907b6e6a22d352740a81381fd4",
    FAVICON: "961e6e85ec16258652bf3c3f322d336fb598b90bccc49821047babfef3925530",
}


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def pinned(path, expected):
    path = path.resolve()
    assert len(expected) == 64 and sha(path) == expected, path
    PINS[path] = expected
    return path


def load(path, h):
    return json.loads(pinned(path, h).read_bytes())


def authority(path, outer_hash, outer_status, native_status):
    outer = load(path / "os_exit_observation.json", outer_hash)
    arts = load(path / "receipt.json", outer["receipt_sha256"])
    term = load(path / "terminal.json", outer["terminal_sha256"])
    assert outer["status"] == outer_status and arts["status"] == native_status
    assert type(outer["pid"]) is int and outer["pid"] > 0 and type(term["pid"]) is int and term["pid"] == outer["pid"]
    assert type(term["exit_code"]) is type(outer["observed_exit_code"]) is int
    assert term["exit_code"] == outer["observed_exit_code"] == 0 and term["status"] == "FINISHED"
    assert term["receipt_sha256"] == outer["receipt_sha256"]
    return arts


def readback():
    for path, value in STATIC.items():
        pinned(path, value)
    pinned(POLICY, "b77eef177dbee7a61823e9f975d1b838a532f9768e4553c1d6a0e7978a05da1b")
    red = authority(RED, "aa3c90101c64ccefe1403af2c296bcb22cafad3b2e35d73539aa7e3c5227f281",
                    "PASS_ROOT_OBSERVED_FIRST_WAVE_INDEPENDENT_GLOBAL_EXIT",
                    "COMPLETE_INDEPENDENT_FIRST_WAVE_SIGNED_ENDPOINT_REPLAY")
    observer = authority(OBS, "6832d622e2f3d45084aa74242e955adbb83f58e026711daa6dca184a144c553c",
                         "PASS_OBSERVED_INDEPENDENT_FIRST_WAVE_SPATIAL_OBSERVER_READBACK_EXIT",
                         "PASS_INDEPENDENT_ACTUAL_FIRST_WAVE_SPATIAL_OBSERVER_READBACK")
    geometry = authority(GEOM, "d39b2d5e97ad6d15bfb68f3f8360e0214c45d11b37dbbd3eaae48313a2e8fc0c",
                         "PASS_ROOT_OBSERVED_REMAINING32_GEOMETRY_REPLAY_EXIT",
                         "PASS_INDEPENDENT_ACTUAL_REMAINING32_SINGLETON_GEOMETRY")
    assert red["exact_spatial_roots"] == red["strict_local_minima"] == 4
    assert observer["all_64_exact_comparisons_recreated"] is True
    assert geometry["successfully_replayed_geometry_cases"] == 32
    old = json.loads(PRIOR.read_bytes())
    cached = json.loads(CACHED.read_bytes())
    assert old["status"] == "PASS_FOUR_INDEPENDENT_FIXED_INPUT_EXACT_SPATIAL_ROOTS_AND_MINIMA"
    assert cached["status"] == "PASS_INDEPENDENT_ACTUAL_36_CACHED_ROOT_OBSERVERS"
    assert list(LABELS) == [d["case_id"] for d in red["case_results"]] == [d["case_id"] for d in observer["cases"]]
    gs = {d["pilot_id"]: d for d in geometry["geometry_records"]}
    counts = {v: Counter() for v in VARIANTS}
    rows, comparisons, hull_noncontained, hull_disjoint = [], [], Counter(), Counter()
    delta = Q(observer["allowance_um_exact"])
    for rr, ob in zip(red["case_results"], observer["cases"], strict=True):
        case = rr["case_id"]
        g_row = gs[case]
        nat = load(ROOT / rr["native_file"], rr["native_sha256"])
        g = load(GEOM / (case + ".json"), g_row["sha256"])
        h = load(HULL / case / "receipt.json", nat["independent_trusted_hull_native_receipt"])
        assert nat["status"] == rr["status"] == "PASS_INDEPENDENT_FIRST_WAVE_EXACT_SPATIAL_ROOT"
        assert g["status"] == "PASS_INDEPENDENT_REMAINING32_WHOLE_BOX_GEOMETRY_ONLY"
        assert h["status"] == "PASS_INDEPENDENT_ACTUAL_REMAINING32_COMPONENT_HULLS"
        assert h["full_independent_hulls_computed"] is True
        assert g["all_144_body_and_24_trace_domains"] is True
        assert h["overlap_or_producer_containment_used_as_proof"] is False
        assert len(h["case_components"]) == 168
        assert Counter(d["kind"] for d in h["case_components"]) == {"bulk": 144, "surface": 24}
        assert all(value == 0 for value in h["disjoint_counts"].values())
        hull_noncontained.update(h["producer_containment_miss_counts"])
        hull_disjoint.update(h["disjoint_counts"])
        a, t, m = (nat[k] for k in ("arithmetic", "target", "observer_and_optional_minimum"))
        assert nat["case_id"] == t["case_id"] == h["case_id"] == g["pilot_id"] == case
        assert t["record_id"] == h["record_id"] == g["record_id"]
        assert ob["record_id"] == "optoGEF2_small_col_" + LABELS[case][1:]
        assert t["trial_box_sha256"] == h["trial_box_sha256"] == g["box_sha256"] == g_row["box_sha256"]
        assert t["point_cache_sha256"] == h["cache_packet_sha256"] == g["point_cache_sha256"]
        assert t["R_binary_sha256"] == h["R_binary_sha256"] == g["R_binary_sha256"]
        assert t["center_binary_sha256"] == h["center_binary_sha256"] == g["center_binary_sha256"]
        assert t["model_contract_sha256"] == h["model_contract_sha256"] == g["model_contract_sha256"]
        assert t["exact_drive"] == h["exact_drive"] == g["exact_drive"]
        assert t["original_observer_field_sha256"] == g["archived_field_sha256"] == ob["original_field_sha256"]
        assert ob["root_proof_sha256"] == rr["native_sha256"]
        r = Q(nat["radius_exact"])
        assert r == Q(1, 2**34) == Q(h["radius_exact"]) == Q(g["radius_exact"])
        zq, bh, eta, z = (Q(a[k]) for k in ("ZQ_upper_exact", "b_H_upper_exact", "eta_upper_exact", "z_upper_exact"))
        q = z + eta / r
        slack, posterior = r * (1 - q), eta / (1 - z)
        sigma, xi, goal = (Q(t[k]) for k in ("archived_center_shift_exact", "original_addition_gap_exact", "observer_goal_exact"))
        curvature = (1 - z) / Q(nat["same_R_norm_exact"])
        eps = Q(67, 2) * (posterior + sigma + xi)
        assert z == zq + bh and 0 <= z <= q < 1 and a["strict_selfmap"] is True
        assert r*q == Q(a["selfmap_lhs_upper_exact"]) and slack == Q(a["selfmap_margin_lower_exact"])
        assert posterior == Q(a["posterior_error_upper_exact"]) == Q(ob["posterior_exact"])
        assert sigma == Q(ob["archived_sigma_exact"]) == Q(g["archived_shift_exact"]) == Q(h["original_archived_shift_exact"])
        assert xi == Q(ob["original_xi_exact"]) == Q(g["archived_addition_gap_exact"])
        assert sigma + posterior < goal == Q(g["archived_observer_goal_exact"])
        assert m["strict_minimum_proved"] is True and m["archived_state_budget_proved"] is True
        assert curvature == Q(m["normalized_H_euclidean_lower_exact"]) > 0
        assert eps == Q(m["physical_nodal_error_um_upper_exact"]) == Q(ob["physical_nodal_error_um_exact"])
        j, G, hoop, monotone = (Q(g["body_J_bounds_exact"][0]), Q(g["min_surface_G_exact"]),
                                Q(g["min_trace_hoop_exact"]), Q(g["strong_monotonicity_lower_exact"]))
        assert min(j, G, hoop, monotone) > 0
        rows.append(dict(case_id=case, label=LABELS[case], record_id=t["record_id"], archived_record_id=ob["record_id"], radius=r, ZQ_upper=zq,
                         b_H_upper=bh, eta_upper=eta, z_upper=z, selfmap_fraction_upper=q,
                         normalized_selfmap_slack_lower=1-q, selfmap_slack_lower=slack, posterior_upper=posterior,
                         sigma=sigma, xi=xi, archived_observer_goal=goal, physical_nodal_error_um_upper=eps,
                         normalized_H_euclidean_lower=curvature, body_J_lower=j, trace_G_lower=G,
                         trace_hoop_lower=hoop, body_monotonicity_lower=monotone,
                         root_proof_sha256=rr["native_sha256"], geometry_proof_sha256=g_row["sha256"],
                         hull_proof_sha256=nat["independent_trusted_hull_native_receipt"]))
        seen = set()
        for item in ob["comparisons"]:
            variant, observable = item["variant"], item["observable"]
            assert variant in VARIANTS and (variant, observable) not in seen
            seen.add((variant, observable))
            stored_lo, stored_hi = map(Q, item["stored_interval_um_exact"])
            infl = Q(item["extra_observer_error_um_exact"])
            lo, hi = map(Q, item["exact_spatial_root_observer_interval_um_exact"])
            obs = Q(item["observed_um_exact"])
            res_lo = max(lo-obs, obs-hi, 0)
            res_hi = max(abs(lo-obs), abs(hi-obs))
            standing = "IN" if res_hi <= delta else "OUT" if res_lo > delta else "UNRESOLVED"
            assert lo == stored_lo - infl and hi == stored_hi + infl and infl >= 0
            assert (res_lo, res_hi) == tuple(map(Q, item["absolute_residual_um_exact"]))
            assert standing == item["standing"]
            counts[variant][standing] += 1
            comparisons.append(dict(case_id=case, label=LABELS[case], variant=variant, observable_key=observable,
                                    stored_lower_um=stored_lo, stored_upper_um=stored_hi, observer_inflation_um=infl,
                                    root_observer_lower_um=lo, root_observer_upper_um=hi, observed_um=obs,
                                    absolute_residual_lower_um=res_lo, absolute_residual_upper_um=res_hi,
                                    allowance_um=delta, standing=standing,
                                    strict_OUT_margin_um=res_lo-delta if standing == "OUT" else ""))
        assert len(seen) == 16
    assert counts == observer["counts"]
    assert len(comparisons) == 64
    out, = (c for c in comparisons if c["standing"] == "OUT")
    saved, = observer["strict_out"]
    assert (out["case_id"], out["variant"], out["observable_key"]) == (saved["case_id"], saved["variant"], saved["observable"])
    assert out["variant"] == "current_disk__regional_mean" and out["label"] == "−K"
    assert out["observable_key"] == "apical_max_shortening_um"
    assert out["strict_OUT_margin_um"] == Q(saved["strict_OUT_margin_um_exact"])
    return rows, comparisons, counts, out, dict(hull_noncontained), dict(hull_disjoint)


def rel(path):
    return Path(os.path.relpath(path, HERE)).as_posix()


def link(label, path):
    return f"[{label}]({rel(path)})"


def dec(q, digits=5):
    with localcontext() as c:
        c.prec = max(digits + 4, 22)
        val = Decimal(q.numerator) / Decimal(q.denominator)
        return format(val, f".{digits}g").replace("e-", "e−")


def directed(q, places, upper):
    a = q * 10**places
    integer = -(-a.numerator // a.denominator) if upper else a.numerator // a.denominator
    return f"{integer // 10**places}.{integer % 10**places:0{places}d}"


def save(path, data):
    with path.open("x", encoding="utf-8", newline="\n") as f:
        f.write(data)


def dump(path, data):
    save(path, json.dumps(data, ensure_ascii=False, indent=2) + "\n")


def csv_exact(path, rows):
    with path.open("x", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with path.open(encoding="utf-8", newline="") as f:
        assert list(csv.DictReader(f)) == [{k: str(v) for k, v in d.items()} for d in rows]


def exact_json(rows):
    return [{k: str(v) if isinstance(v, Q) else v for k, v in d.items()} for d in rows]


def plot(rows):
    css = ["#315a83", "#b6502b", "#dde4ea", "#5a778c"]
    fields = ("z_upper", "eta_over_radius", "normalized_selfmap_slack_lower",
              "posterior_1e9_um", "shift_1e9_um", "addition_gap_1e9_um")
    plot_rows = []
    for r in rows:
        data = (r["z_upper"], r["eta_upper"]/r["radius"], r["normalized_selfmap_slack_lower"],
                *(r[k] * Q(67, 2) * 10**9 for k in ("posterior_upper", "sigma", "xi")))
        assert sum(data[:3]) == 1 and sum(data[3:]) == r["physical_nodal_error_um_upper"] * 10**9
        plot_rows.append(dict(case_id=r["case_id"], label=r["label"], **dict(zip(fields, data, strict=True))))
    dump(HERE / "figure_data.json", {
        "full_case_order": list(LABELS),
        "axis_a": "dimensionless fraction of self-map radius",
        "axis_b": "physical nodal error contribution in 10^-9 micrometres",
        "rounding": "Exact rational channels projected to binary float only for drawing.",
        "data_exact": exact_json(plot_rows),
    })
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.spines.left": False, "axes.labelcolor": "#192b3b",
                         "axes.edgecolor": "#a4b1bc", "xtick.color": "#314759", "ytick.color": "#192b3b"})
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.2), layout="constrained")
    labels = [r["label"] for r in plot_rows]
    panels = [
        (fields[:3], [r"$z$", r"$\eta/r$", "Remaining headroom"],
         css[:3], "A   Certified self-map budget", r"Fraction of radius $r$"),
        (fields[3:], ["Posterior", "Archived-center shift", "Binary addition gap"],
         [css[0], css[1], css[3]], "B   Error transferred to the archived coordinates",
         r"Physical nodal error ($10^{-9}\,\mathrm{\mu m}$)"),
    ]
    for i, (ax, (keys, names, colors, title, xlabel)) in enumerate(zip(axes, panels, strict=True)):
        left = [0.0] * len(plot_rows)
        for key, name, color in zip(keys, names, colors, strict=True):
            vals = [float(r[key]) for r in plot_rows]
            ax.barh(labels, vals, left=left, height=.48, color=color, label=name)
            left = [x + y for x, y in zip(left, vals, strict=True)]
        ax.invert_yaxis()
        ax.set_title(title, loc="left", fontsize=13, fontweight="bold", color="#192b3b", pad=14)
        ax.set_xlabel(xlabel, labelpad=4)
        ax.tick_params(axis="y", length=0, pad=6)
        ax.set_axisbelow(True)
        ax.grid(axis="x", color="#e5e9ed", linewidth=.75)
        ax.legend(loc="upper left", bbox_to_anchor=(-.01, -.21), frameon=False, ncol=3,
                  fontsize=9, handlelength=1, columnspacing=1.2)
        if i == 0:
            ax.set_xlim(0, 1.025)
            ax.xaxis.set_major_locator(MultipleLocator(.25))
            ax.xaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
            ax.axvline(1, color="#192b3b", linestyle="--", linewidth=.9)
        else:
            ax.set_xlim(0, 3.45)
            ax.xaxis.set_major_locator(MultipleLocator(.5))
    with FIGURE.open("xb") as f:
        fig.savefig(f, format="png", dpi=210, facecolor="#ffffff")
    plt.close(fig)


def md_table(header, body):
    return "\n".join(["| " + " | ".join(header) + " |", "| " + " | ".join("---" for _ in header) + " |"]
                     + ["| " + " | ".join(map(str, r)) + " |" for r in body])


def markdown(rows, counts, out, noncontained):
    title = "# " + TITLE
    certificate = md_table(("Case", r"$z$", "Headroom", "Curvature lower"),
                           [(d["label"], dec(d["z_upper"]), dec(d["normalized_selfmap_slack_lower"]*100, 4)+"%",
                             dec(d["normalized_H_euclidean_lower"], 6)) for d in rows])
    transfer = md_table(("Case", "Posterior", "Center shift", "Nodal error (µm)"),
                        [(d["label"], dec(d["posterior_upper"]), dec(d["sigma"]), dec(d["physical_nodal_error_um_upper"])) for d in rows])
    geom = md_table(("Case", r"Body $J$ lower", r"Trace $G$ lower", "Hoop lower"),
                    [(d["label"], directed(d["body_J_lower"], 5, False), directed(d["trace_G_lower"], 5, False),
                      directed(d["trace_hoop_lower"], 5, False)) for d in rows])
    obs = md_table(("Convention", "IN", "OUT"), [(name, counts[key]["IN"], counts[key]["OUT"]) for key, name in VARIANTS.items()])
    d = directed(out["strict_OUT_margin_um"], 9, False)
    res_lo = directed(out["absolute_residual_lower_um"], 9, False)
    res_hi = directed(out["absolute_residual_upper_um"], 9, True)
    allowance = directed(out["allowance_um"], 9, True)
    return rf"""{title}

Four additional fixed-input finite-element equilibria have a unique root in their respective certified cubes and a strict local energy minimum. The cumulative count is **8 of 36** roots with exact spatial integration. All **36 of 36** also have independently certified cached-quadrature roots. The four earlier exact-spatial certificates are documented separately. {link('Earlier exact-spatial checkpoint', PRIOR_PAGE)} · {link('Cached-quadrature checkpoint', CACHED_PAGE)}

## Mechanics and contraction

The fixed cubes have radius $r=2^{{-34}}$ in normalized free coordinates. The signed residual gives $\eta\geq\lVert RF_I(u_c)\rVert_\infty$ at the center. The whole-cube bound $z=Z_Q+b_H$ controls $\lVert I-RH_I\rVert_\infty$ using the cached derivative bound and the spatial Hessian correction. All four satisfy the strict self-map inequality.

$$
\eta/r+z<1.
$$

Headroom is $1-z-\eta/r$, shown as a percentage of the cube radius. The positive symmetric part of the same $R$ and the whole-cube derivative bound give the following lower bound on the normalized energy Hessian.

$$
\lambda_{{\min}}(H_I) \ge (1-z)/\lVert R\rVert_\infty > 0.
$$

{certificate}

Here **+** and **−** identify the two fixed parameter endpoints. B, E and K are archived small-field record labels. Numeric tables display rounded values, with the geometry lower bounds rounded down. {link('Exact rational bounds and case identities', HERE/'exact_bounds.csv')}

![Four horizontal case bars in each panel. The first panel shows positive self-map headroom in every case. The second shows that the +E and −K archived-center shifts contribute most of their physical error.](first_wave_bounds.png)

## Geometry and archived observer coordinates

Independent geometry checks covered all 144 body elements and 24 surface traces in each complete cube. They give positive body Jacobian, positive trace metric, positive trace hoop, and a positive body monotonicity bound. The table shows the first three lower bounds.

{geom}

The posterior is $e=\eta/(1-z)$. The center shift $\sigma$ is zero for the two B evaluations. For +E and −K, the proof center shifts while observer comparisons continue to use the original archived fields. The transfer bound includes the original binary coordinate-addition gap $\xi$.

$$
\epsilon = 33.5(e+\sigma+\xi)\,\mu\mathrm{{m}}.
$$
 The posterior and shift columns are in normalized nodal coordinates.

{transfer}

{link('Independent geometry closure', GEOM/'receipt.json')} · {link('Independent global signed reduction', RED/'receipt.json')} · {link('Figure data with exact rational endpoints', HERE/'figure_data.json')}

## Four observers per convention

The same frozen model-discrepancy allowance applies to all comparisons. The very small nodal bounds above quantify numerical error in reaching a root. The allowance governs mismatch between the model observer and the archived measurement. These bounds leave physical accuracy and experimental detectability open. Each of the four fields contributes four observables to each convention. All **16 primary comparisons are IN**. Of **48 alternative comparisons**, 47 are IN and one is OUT.

{obs}

The OUT result is the **apical regional mean at −K** under disk averaging. Its absolute residual lies in **[{res_lo}, {res_hi}] µm** after outward decimal rounding. The allowance is at most **{allowance} µm** at the displayed precision. The exact difference between the residual lower endpoint and the exact allowance is positive, at least **{d} µm** after rounding down. The archived data key names the apical channel; the disk regional-mean convention determines the statistic in this row. {link('All 64 comparisons and the exact OUT deficit', HERE/'exact_observer_comparisons.csv')} · {link('Independent observer readback', OBS/'receipt.json')}

## Independent intervals and scope

The independent component replay enclosed the signed point force correction, the whole-cube Hessian correction, and the same-coordinate cache correction separately. It retained the full interval hull of each producer and independently valid enclosure for the signed reduction. Among the 672 components, {noncontained['point_force']} force entries and {noncontained['whole_box_H']} Hessian entries had independent enclosures that extended outside the corresponding producer interval. No channel had a disjoint pair. The independent validity of the intervals supports the reduction. Their overlap is a consistency check. {link('Signed interval hull policy', POLICY)}

The certificates cover eight fixed inputs and their stated cubes. Twenty-eight exact-spatial targets remain outside this checkpoint. The observer outcome is sensitive to the stated averaging convention. The fixed-input cubes and observer inequalities leave source compatibility, loading stability, and continuum adequacy open.
"""


CSS = """
:root{color-scheme:dark}
*{box-sizing:border-box}
body{background:#111923;color:#e8edf3;margin:0;font:1rem/1.58 system-ui,sans-serif}
main{width:min(100% - 32px,900px);margin:0 auto;padding:35px 0 50px}
nav{display:flex;flex-wrap:wrap;gap:.5rem 1.1rem;margin:0 0 1.4rem;font-size:.9rem}
h1,h2{line-height:1.16;letter-spacing:-.025em;color:#fff}
h1{font-size:clamp(1.95rem,5vw,3.45rem);margin:.3rem 0 1.25rem}
h2{font-size:clamp(1.5rem,4vw,2.05rem);margin:2.3rem 0 .85rem}
a{color:#9bd1f7;text-underline-offset:.19em;overflow-wrap:anywhere}
strong{color:#fff}
img{display:block;width:100%;height:auto;border-radius:.35rem;margin:1.1rem 0}
table{border-collapse:collapse;width:100%;margin:1em 0 1.2em;font-variant-numeric:tabular-nums}
th,td{text-align:left;padding:.6rem .55rem;border-bottom:1px solid #3c4a5c}
th{color:#fff;font-weight:650;background:#202d3b}
tbody tr:nth-child(even){background:#16212d}
@media(max-width:520px){main{padding-top:22px}th,td{font-size:.83rem;padding:.45rem .24rem}nav{gap:.35rem .7rem}}
"""


def render(md):
    spec = spec_from_file_location("checkpoint_katex", HELPER)
    helper = module_from_spec(spec)
    spec.loader.exec_module(helper)
    sys.path.insert(0, str(ROOT / "scripts"))
    import build_report
    assert Path(build_report.__file__).resolve() == CONVERT
    md = re.sub(r'\$\$(.*?)\$\$', lambda m: r'\[' + m[1] + r'\]', md, flags=re.S)
    md = re.sub(r'(?<!\$)\$([^$\n]+)\$(?!\$)', lambda m: r'\(' + m[1] + r'\)', md)
    body, count = helper.compile_math(build_report.convert(md))
    page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width,initial-scale=1"><title>{escape(TITLE)}</title>'
            f'<link rel="icon" type="image/png" href="{rel(FAVICON)}">'
            f'<link rel="stylesheet" href="{rel(ASSETS/"katex.min.css")}"><style>{CSS}</style></head>'
            f'<body><main><nav aria-label="Report formats"><a href="REPORT.md">Markdown</a>'
            f'<a href="exact_bounds.csv">Exact bounds CSV</a><a href="exact_observer_comparisons.csv">Observer CSV</a>'
            f'</nav>{body}</main></body></html>\n')
    parser = helper.LocalLinks()
    parser.feed(page)
    links = {}
    for v in parser.paths:
        target = (HERE / v).resolve()
        assert target.is_relative_to(SITE) and target.is_file(), v
        links[v] = target.relative_to(SITE).as_posix()
    return page, count, links


def main():
    rows, comparisons, counts, out, noncontained, disjoint = readback()
    csv_exact(HERE / "exact_bounds.csv", rows)
    csv_exact(HERE / "exact_observer_comparisons.csv", comparisons)
    plot(rows)
    md = markdown(rows, counts, out, noncontained)
    save(HERE / "REPORT.md", md)
    html, count, links = render(md)
    save(HERE / "REPORT.html", html)
    assert all(sha(d) == h for d, h in PINS.items())
    readback_path = HERE / "scientific_readback.json"
    dump(readback_path, {
        "status": "PASS_EXACT_CHECKPOINT_REPORT_SCALARS_AND_PROOF_BINDINGS",
        "cases": exact_json(rows), "observer_comparison_counts": counts, "strict_OUT": exact_json([out])[0],
        "number_of_observer_comparisons_recomputed_from_exact_endpoints": 64,
        "producer_noncontainment_counts": noncontained, "disjoint_counts": disjoint,
        "cumulative_exact_spatial_roots": 8, "total_cached_roots": 36,
        "figure_projection_channels_each_case": 6, "shared_premises": "Native prior mathematical proofs are reused.",
        "full_five_gate_complete": False,
    })
    binding = HERE / "report_build_receipt.json"
    dump(binding, {
        "status": "PASS_BUILD_EXACT_SPATIAL_FIRST_NEW_WAVE_CHECKPOINT",
        "source_sha256": sha(Path(__file__)), "rendered_math_count": count, "local_links": links,
        "input_sha256": {p.relative_to(SITE).as_posix(): digest for p, digest in PINS.items()},
        "file_sha256": {name: sha(HERE/name) for name in ("REPORT.md", "REPORT.html", "exact_bounds.csv",
                        "exact_observer_comparisons.csv", FIGURE.name, "figure_data.json", readback_path.name)},
        "expected_table_rows": 16,
    })
    print(json.dumps({"build_receipt_sha256": sha(binding), "readback_sha256": sha(readback_path),
                      "compiled_math": count, "local_links": len(links), "checks": "4 cases and 64 observer comparisons"}))


if __name__ == "__main__":
    main()

