"""Build a standalone, source-pinned checkpoint for the cached observer result."""
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction as Q
from hashlib import sha256
import html
import json
from pathlib import Path

import markdown
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DESIGN = HERE.parents[1]
RUN = DESIGN.parent
RECEIPT = HERE.parent/"observer_application_v1/independent_readback/receipt.json"
ALT = RUN/"observer_robustness/independent_actual/receipt.json"
SPATIAL = DESIGN/"exact_targets_v1/four_spatial_checkpoint_v1/SPATIAL.html"
PINS = {
    RECEIPT:"a06549f048abba4eb2b4d1c1fb398f06e8f524196c92a673bd0493c527c8f40d",
    ALT:"703b574bc51a083b10c7328d9622862e939ddabe791b18ce9e1817f07718f69f",
    SPATIAL:"6dd9100643191f6617f5032d2a9d4e35ef3e179c9826115c841b041ce92c9e4d",
}
ROOT = next(p for p in HERE.parents if (p/"PLAN.md").is_file())
ENDPOINTS = ("max_contrast_change","max_negative_contrast_change")


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def text(p,value):
    with p.open("x",encoding="utf-8",newline="\n") as sink:
        sink.write(value)


def directed(q, digits, direction):
    with localcontext() as ctx:
        ctx.prec=90
        n=Decimal(q.numerator)/Decimal(q.denominator)
        rounded=n.quantize(Decimal(10)**(-digits),rounding=direction)
        assert Q(rounded)<=q if direction==ROUND_FLOOR else Q(rounded)>=q
        return str(rounded)


def display(q, digits=12):
    with localcontext() as ctx:
        ctx.prec=90
        return f"{Decimal(q.numerator)/Decimal(q.denominator):.{digits}f}"


def show(interval,digits=11):
    lo,hi=interval
    return f"[{directed(lo,digits,ROUND_FLOOR)}, {directed(hi,digits,ROUND_CEILING)}]"


def pair(values):
    a,b=map(Q,values)
    assert a<=b
    return a,b


def plot(q, gap, path):
    bg="#121c2b"
    foreground="#eff5fd"
    f, ax=plt.subplots(figsize=(7.4,2.65),dpi=180,facecolor=bg)
    ax.set_facecolor(bg)
    values=[float(sum(q[k])/2) for k in ENDPOINTS]
    colors=["#62dfc6","#c4b2ff"]
    for value,y,color in zip(values,(1,0),colors,strict=True):
        ax.scatter(value,y,s=94,color=color,edgecolor=bg,linewidth=1.4,zorder=4)
        ax.text(value+0.009,y,f"{value:.6f}",color=color,fontsize=11,va="center")
    ax.plot(values,[1,0],color="#8192ac",linewidth=1.1,linestyle=":",zorder=2)
    ax.set_yticks([1,0],["Positive direction","Negative direction"])
    ax.set_xlim(-4.315,-4.075)
    ax.set_ylim(-.42,1.42)
    ax.set_xticks([-4.3,-4.25,-4.2,-4.15,-4.1])
    ax.set_xlabel("Apical center minus shoulder, Q (µm)",color=foreground,labelpad=10)
    ax.set_title(f"Signed gap ≈ {float(sum(gap)/2):.6f} µm",loc="left",color=foreground,fontsize=14,pad=10)
    ax.tick_params(axis="both",colors=foreground,length=0,pad=8,labelsize=10)
    ax.grid(axis="x",color="#3b4b64",alpha=.55,linewidth=.75)
    for spine in ax.spines.values():
        spine.set_visible(False)
    f.subplots_adjust(left=.25,right=.98,bottom=.27,top=.79)
    assert not path.exists()
    f.savefig(path,facecolor=bg,metadata={"Title":"Cached quadrature large-region apical contrast"})
    plt.close(f)


def main():
    for p,want in PINS.items():
        assert digest(p)==want,p
    reviewed=json.loads(RECEIPT.read_bytes())
    alt=json.loads(ALT.read_bytes())
    assert reviewed["status"]=="PASS_INDEPENDENT_ACTUAL_36_CACHED_ROOT_OBSERVERS"
    assert reviewed["cached_fixed_parameter_roots"]==36 and reviewed["primary_signed_residuals"]==136
    assert reviewed["scope"]=="current_disk__spatial_max" and all(
        reviewed["counts"][k]=={"IN":68,"OUT":0,"UNRESOLVED":0} for k in ENDPOINTS)
    assert reviewed["new_novelty_claim"] is reviewed["source_compatibility_or_PDE"] is reviewed["full_five_gate_goal_achieved"] is False
    assert reviewed["exact_spatial_observer_applications_from_this_readback"]==0
    delta=Q(reviewed["frozen_model_discrepancy_um_exact"])
    assert delta==Q(alt["source_configured_delta_um_exact"])
    q={row["endpoint_id"]:pair(row["Q_cached_um_exact"]) for row in reviewed["large_Q_cached"]}
    gap=pair(reviewed["cache_root_signed_large_gap_um_exact"])
    p,n=(q[k] for k in ENDPOINTS)
    assert gap==(p[0]-n[1],p[1]-n[0]) and gap[0]>0
    slack=Q(reviewed["minimum_slack"]["strict_slack_um_exact"])
    assert slack==min(Q(r["strict_slack_um_exact"]) for r in reviewed["all_rows"])>0
    assert reviewed["minimum_slack"]["case_id"]=="max_negative_contrast_change_small_17"
    assert reviewed["minimum_slack"]["observable"]=="region_averaged_shortening_um"
    labels={
        "current_disk__regional_mean":"Disk weighting, component regional means",
        "current_radial_line__spatial_max":"Radial-line weighting, component spatial maxima",
        "current_radial_line__regional_mean":"Radial-line weighting, component regional means",
    }
    alternate_rows=[]
    for k,label in labels.items():
        row=alt["by_observer_common_requirements"][k]
        lo,hi=pair(row["independent_required_allowance_um_exact"])
        assert delta<lo<=hi and row["status_at_original_delta"]=="ALLOWANCE_EXCEEDED"
        alternate_rows.append(f"| {label} | {display(hi,6)} µm |")
    plot(q,gap,HERE/"large_contrast.png")
    markdown_text=f"""# Observer comparisons at 36 cached-quadrature equilibria

Both frozen parameter choices satisfy all 68 primary small-region comparisons at the same conditional allowance of **{display(delta)} µm**. They still differ on the model’s large-region apical contrast. The certified positive-direction minus negative-direction gap is **{show(gap)} µm**. [Independent observer readback](../observer_application_v1/independent_readback/receipt.json)

## Calibration and prediction

Each candidate contributes 17 small-region equilibria and one large-region equilibrium. At each small state the primary observer compares maximum shortening, disk-averaged shortening and the separate apical and basal spatial maxima on a fixed region in current radius.

| Frozen candidate | Primary comparisons within allowance | Large-region Q (µm) |
| --- | --- | --- |
| Positive direction | 68 / 68 | {show(p)} |
| Negative direction | 68 / 68 | {show(n)} |

The smallest slack across all 136 comparisons is about **{display(slack)} µm**. It occurs at the negative-direction L state for region-averaged shortening. Slack is the frozen allowance minus the upper endpoint of the absolute residual interval.

![Large-region Q centers at minus 4.136003 and minus 4.290758 micrometers. The positive-direction value is greater.](large_contrast.png)

The plot marks the centers of the enclosures. The certified intervals are narrower than the plot resolution. The large observable **Q = Zₐ(0) − Zₐ(Rₛ)** compares the apical height on the axis with the height at the same fixed current shoulder radius Rₛ in each field. Both Q values are negative. The signed gap is positive because the positive-direction candidate has the higher Q. These values are model predictions at the prescribed large input.

## From a root bound to an observer bound

Let **e** be the certified root error measured from the verification center in dimensionless free coordinates. Let **ξ** be the original binary storage gap between the exact nodal sum and the saved coordinates. For verification at the original center, each nodal coordinate differs from the archived field by at most **33.5(e + ξ) µm**. For the two recentered fields, positive E and negative K, the actual binary center shift **σ** gives the bound **33.5(σ + e + ξ) µm**. The comparison uses the original saved observer and the original ξ in both cases.

A strictly positive radial derivative bound and fixed axis and rim radii keep the quadratic traces invertible over the same region. The nodal error bounds the apical and basal heights at fixed current radius. Their error bounds add for total shortening. Each component maximum uses its own surface bound. The bound on Q adds the axis-node vertical error and the same-radius apical shoulder error. Every comparison uses the worst absolute residual over the full inflated interval. [Observer transfer method](../observer_application_v1/METHOD.md)

## Dependence on observer convention

The primary convention uses current-radius disk weighting and spatial maxima for the separate components. In the earlier stored-field comparison, the three alternate conventions each exceeded the same frozen allowance for the pair.

| Alternate stored-field convention | Required common allowance, approximately |
| --- | --- |
{chr(10).join(alternate_rows)}

The allowance is a chosen model discrepancy bound. Its calibration as experimental noise and the experimental observation mapping remain unresolved. This comparison covers fixed states of the finite, cached-quadrature model. Exact-spatial integration has separate certificates for four fixed-input equilibria and minima. The full 36-state exact-spatial observer comparison, loading-path validity, continuum adequacy and novelty remain open. [Alternate-observer review](../../observer_robustness/independent_actual/receipt.json) · [Four exact-spatial certificates](../ex-href)
"""
    markdown_text=markdown_text.replace("../ex-href","../../exact_targets_v1/four_spatial_checkpoint_v1/SPATIAL.html").replace("../../observer_robustness/","../../../observer_robustness/")
    marked=markdown.markdown(markdown_text,extensions=["tables"])
    marked=marked.replace("<table>",'<div class="table-wrap"><table>').replace("</table>","</table></div>")
    marked=marked.replace('<p><img ', '<figure><img ').replace('src="large_contrast.png" /></p>', 'src="large_contrast.png" /></figure>')
    css="""
:root {color-scheme:dark;font:17px/1.62 system-ui,sans-serif;background:#0d1420;color:#eaf1fb}
* {box-sizing:border-box}
body {margin:0}
main {max-width:900px;margin:0 auto;padding:clamp(20px,5vw,60px)}
h1,h2 {line-height:1.12;text-wrap:balance}
h1 {font-size:clamp(2rem,5.4vw,3.25rem);margin:0 0 1rem;letter-spacing:-.04em}
h2 {font-size:clamp(1.35rem,3.1vw,1.9rem);margin:2.25rem 0 .8rem;border-top:1px solid #364963;padding-top:1rem}
p {margin:.75rem 0 1rem}
strong {color:#fff}
a {color:#84dec9;text-underline-offset:.2em}
a:focus-visible {outline:2px solid #fff;outline-offset:3px}
.table-wrap {margin:1.3rem 0; overflow-wrap:anywhere}
table {width:100%;border-collapse:collapse;table-layout:fixed;font-variant-numeric:tabular-nums;background:#121c2b}
td,th {padding:.72rem;text-align:left;vertical-align:top;border-bottom:1px solid #30425b}
th {color:#b9c9df;font-size:.9rem}
figure {margin:1.3rem 0;border:1px solid #364963;border-radius:14px;background:#121c2b;overflow:hidden}
img {display:block;width:100%;height:auto}
@media(max-width:540px){:root{font-size:16px}td,th{padding:.42rem}table{font-size:.82rem}}
"""
    page=f"""<!doctype html>
<html lang="en"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape("Observer comparisons at 36 cached-quadrature equilibria")}</title><style>{css}</style></head>
<body><main>{marked}</main></body></html>
"""
    text(HERE/"OBSERVERS.md",markdown_text)
    text(HERE/"OBSERVERS.html",page)
    receipt={
        "status":"BUILT_SCOPED_CACHED_OBSERVER_CHECKPOINT",
        "source_sha256":digest(Path(__file__)),"source_inputs":{p.relative_to(ROOT).as_posix():h for p,h in PINS.items()},
        "artifacts":{n:digest(HERE/n) for n in ("OBSERVERS.md","OBSERVERS.html","large_contrast.png")},
        "counts":reviewed["counts"],"frozen_delta_exact":str(delta),"min_slack_exact":str(slack),
        "gap_um_exact":list(map(str,gap)),"Q_um_exact":{k:list(map(str,v)) for k,v in q.items()},
        "build_scope":"Static report; 36 cached-quadrature equilibria and the primary observer. No new computation of equilibrium or predictions.",
        "full_goal_achieved":False,
    }
    text(HERE/"build_receipt.json",json.dumps(receipt,indent=2)+"\n")


if __name__=="__main__":
    main()

