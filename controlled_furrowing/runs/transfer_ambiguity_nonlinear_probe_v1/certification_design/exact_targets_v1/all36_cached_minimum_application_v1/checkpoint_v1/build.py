"""Build a local report and scientific plot from the independent exact receipt."""
import csv
from fractions import Fraction
from hashlib import sha256
from html import escape
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
READBACK=HERE.parent/"independent_readback"
RECEIPT="88cd552e7240bffd7f1c6596eb6671fa2a9faf5836a9ae086ed7e96f538a7aaf"
OS="68c81d84539b334a7e826f12cc57f473b49c7d9068d6a2c242891a808ef0d957"
TITLE="36 strict minima of the cached-quadrature energies"
R34=Fraction(1,2**34)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def verified(path,h):
    assert sha(path)==h, path
    return json.loads(path.read_bytes())


def write(path,s):
    with path.open("x",encoding="utf-8",newline="\n") as out:
        out.write(s)


def label(case):
    return ("−" if "negative" in case else "+")+case.rsplit("_change_",1)[-1].replace("_"," ")


def value(s):
    return float(Fraction(s))


def series(x):
    return f"{value(x):.8e}"


def plot(rows,order):
    fig,ax=plt.subplots(2,1,figsize=(8.7,6.4),sharex=True)
    ink="#273c4a"
    fig.patch.set_facecolor("#faf7ed")
    coords={name:i for i,name in enumerate(order)}
    coords["large_03"]+=.85
    groups=(("Positive endpoint","#006767","o",0.10),("Negative endpoint","#a4342a","s",-0.10))
    cols=(("cached_H_eigenvalue_lower_exact",10**4,"Curvature lower bound (×10⁻⁴)"),
          ("strict_cube_interior_margin_exact",10**11,"Cube interior margin (×10⁻¹¹)"))
    for axes,(key,scale,title) in zip(ax,cols,strict=True):
        axes.set_facecolor("#ffffff")
        for lab,col,symbol,offset in groups:
            data=[r for r in rows if ("negative" in r["pilot_id"])==(offset<0)]
            xy=[(coords[r["pilot_id"].rsplit("_change_",1)[-1]]+offset,float(Fraction(r[key])*scale)) for r in data]
            handle=axes.scatter(*zip(*xy,strict=True),s=38,label=lab,color=col,marker=symbol,edgecolors="white",linewidths=.6,zorder=3)
            assert [list(v) for v in handle.get_offsets().tolist()] == [list(v) for v in xy]
        axes.set_title(title,loc="left",fontsize=14,fontweight="semibold",color=ink,pad=9)
        axes.set_xlim(-.55,len(order)+.5)
        axes.set_ylim(0,max(float(Fraction(r[key])*scale) for r in rows)+.4)
        axes.tick_params(colors=ink,labelsize=10)
        axes.grid(axis="y",color="#dce3df",linewidth=.85)
        for spine in ("top","right"):
            axes.spines[spine].set_visible(False)
        for spine in ("left","bottom"):
            axes.spines[spine].set_color("#879b9f")
        axes.axvline(16.9,color="#879b9f",linestyle=":",linewidth=.9)
    ticks=[0,2,4,6,8,10,12,14,16,coords["large_03"]]
    ax[1].set_xticks(ticks,["01","03","05","07","09","11","13","15","17","Large 03"])
    ax[1].set_xlabel("Source record",fontsize=11,color=ink)
    handles,text=ax[0].get_legend_handles_labels()
    fig.legend(handles,text,loc="upper center",ncol=2,frameon=False,labelcolor=ink,fontsize=11)
    fig.subplots_adjust(left=.12,right=.975,top=.89,bottom=.115,hspace=.35)
    img=HERE/"cached_minima.png"
    assert not img.exists()
    fig.savefig(img,dpi=210,facecolor=fig.get_facecolor())
    plt.close(fig)


INTRO="Each of the 36 fixed parameter and drive cases has a unique minimum of its normalized cached-quadrature energy within a matched closed cube. Each minimum lies strictly inside its cube."
THEOREM="With the same preconditioner R and whole-cube defect Z<sub>Q</sub>, the lower Hessian bound is <i>m</i> = (1 − Z<sub>Q</sub>)/‖R‖<sub>∞</sub> &gt; 0. The proof also requires an independently verified positive-definite symmetric part of R and a smooth cached energy. For every u in the cube with minimizer u<sub>*</sub>, the energy gap is at least (m/2)‖u − u<sub>*</sub>‖<sub>2</sub><sup>2</sup>."
AUX="For +small 02, the native posterior of 9.05979 × 10⁻¹¹ came from the larger radius 2⁻³². The unchanged point force and a fresh whole-cube derivative bound give a strict self-map on the smaller cube, with posterior 6.59641 × 10⁻¹². The two cubes share the same center, so uniqueness on the larger cube identifies the same root. The frozen observer error uses the original native posterior. The two recentered proofs retain their actual shifts σ and archived binary addition gaps ξ."
SCOPE="These minima concern the frozen binary cached-quadrature energies. Exact-spatial finite-element integration, continuum equations, loading stability, and experimental compatibility require separate evidence."
CAPTION="For each of the 18 source records, the two markers show the positive and negative endpoint. Curvature is for the normalized cached energy. The interior margin is r − e in dimensionless free coordinates, for r = 2⁻³⁴."


def output(rows,doc):
    emin=doc["minimum_curvature_case"]
    imin=doc["minimum_interior_case"]
    table=[(r["pilot_id"],label(r["pilot_id"]),series(r["cached_H_eigenvalue_lower_exact"]),series(r["strict_cube_interior_margin_exact"])) for r in rows]
    assert len(table)==len(set(t[0] for t in table))==36
    with (HERE/"bounds.csv").open("x",encoding="utf-8",newline="") as f:
        w=csv.writer(f)
        fields=("pilot_id","cached_H_eigenvalue_lower_exact","strict_cube_interior_margin_exact")
        w.writerow(fields)
        w.writerows([r[k] for k in fields] for r in rows)
    links=(("Exact scalar readback","../independent_readback/receipt.json"),
           ("Proof and interpretation","../independent_readback/READBACK.md"),
           ("Observed exit","../independent_readback/os_exit_observation.json"),
           ("Exact plot data (CSV)","bounds.csv"))
    metric1=f"{value(doc['minimum_curvature_exact'])*1e4:.5f} × 10⁻⁴"
    metric2=f"{value(doc['minimum_interior_margin_exact'])*1e11:.5f} × 10⁻¹¹"
    md=[f"# {TITLE}","",INTRO,"",
        f"The smallest normalized Hessian eigenvalue lower bound is **{metric1}** ({label(emin)}). The smallest certified interior margin is **{metric2}** ({label(imin)}).",
        "", "For a whole-cube defect $Z_Q$ and the fixed preconditioner $R$, the lower Hessian bound is $m=(1-Z_Q)/\\|R\\|_\\infty>0$. The same $R$ has an independently verified positive-definite symmetric part. Smoothness and the lower Hessian bound give an energy gap of at least $(m/2)\\|u-u_*\\|_2^2$ for every $u$ in the cube.",
        "", "![Normalized curvature lower bounds and positive cube interior margins for all 36 cases](cached_minima.png)","",CAPTION,"",
        "## Root localization and observer errors","",AUX,"",SCOPE,"",
        "## All 36 bounds","","| Case | Normalized Hessian lower bound | Cube interior margin |","|---|---:|---:|",
        *[f"| {name} | {m} | {margin} |" for _,name,m,margin in table],"",
        " · ".join(f"[{a}]({b})" for a,b in links),""]
    write(HERE/"MINIMA.md","\n".join(md))
    td="".join(f'<tr data-case="{escape(key)}"><th scope="row">{escape(name)}</th><td>{m}</td><td>{margin}</td></tr>' for key,name,m,margin in table)
    linkhtml="".join(f'<li><a href="{href}">{escape(name)}</a></li>' for name,href in links)
    page=f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{TITLE}</title><style>
:root{{color-scheme:dark;background:#101c23;color:#e9efed;font:16px/1.62 system-ui,sans-serif}}
*{{box-sizing:border-box}}body{{margin:0}}main{{max-width:950px;margin:auto;padding:40px 30px 50px}}
h1{{font-weight:670;line-height:1.08;font-size:clamp(2.2rem,5vw,3.7rem);letter-spacing:-.045em;margin:.3em 0}}
h2{{font-weight:630;line-height:1.24;margin:34px 0 10px}}p{{margin:0 0 16px}}
.metrics{{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin:22px 0}}
.metrics div{{background:#172932;border:1px solid #32515d;padding:13px 16px;border-radius:8px}}
.metrics strong{{display:block;color:#90dddd;font-size:1.55em;letter-spacing:-.02em}}
.metrics span,figcaption{{color:#c8d7d7}}figure{{margin:26px 0}}figure img{{display:block;width:100%;height:auto;border-radius:8px}}
figcaption{{font-size:.92em;margin-top:8px}}.scope{{border-left:2px solid #6bb9b6;padding-left:13px}}
table{{width:100%;border-collapse:collapse;font-size:.94rem;font-variant-numeric:tabular-nums}}
thead{{background:#203640}}th,td{{padding:7px 9px;border-bottom:1px solid #3a535f;text-align:right}}tr th:first-child{{text-align:left}}
tbody tr:nth-child(2n){{background:#15262e}}a{{color:#9ae4e1}}ul{{padding-left:21px}}
@media(max-width:540px){{main{{padding:24px 15px 30px}}.metrics{{grid-template-columns:1fr}}th,td{{padding:6px 5px;font-size:.83rem}}}}
</style></head><body><main>
<h1>{TITLE}</h1><p>{INTRO}</p>
<div class="metrics"><div><strong>{metric1}</strong><span>Smallest curvature lower bound · {escape(label(emin))}</span></div>
<div><strong>{metric2}</strong><span>Smallest interior margin · {escape(label(imin))}</span></div></div>
<p>{THEOREM}</p>
<figure><img src="cached_minima.png" alt="Normalized Hessian lower bounds and strict cube interior margins for all 36 cases. The exact numerical values follow in a table."><figcaption>{CAPTION}</figcaption></figure>
<section><h2>Root localization and observer errors</h2><p>{AUX}</p><p class="scope">{SCOPE}</p></section>
<section><h2>All 36 bounds</h2>
<table aria-label="Normalized cached Hessian lower bounds and strict cube interior margins"><thead><tr><th scope="col">Case</th><th scope="col">Curvature bound</th><th scope="col">Interior margin</th></tr></thead><tbody>{td}</tbody></table></section>
<section><h2>Evidence files</h2><ul>{linkhtml}</ul></section>
</main></body></html>
"""
    write(HERE/"MINIMA.html",page)


def main():
    doc=verified(READBACK/"receipt.json",RECEIPT)
    os=verified(READBACK/"os_exit_observation.json",OS)
    assert doc["status"]=="PASS_INDEPENDENT_ACTUAL_ALL36_CACHED_ENERGY_LOCAL_MINIMA_SCALAR_READBACK" and os["observed_exit_code"]==0
    assert doc["verified_cached_local_minima"]==36 and os["receipt_sha256"]==RECEIPT
    order=[f"small_{i:02d}" for i in range(1,18)]+["large_03"]
    rows=sorted(doc["certificate_cases"],key=lambda r:(order.index(r["pilot_id"].rsplit("_change_",1)[-1]),"negative" in r["pilot_id"]))
    assert all(Fraction(r["cached_H_eigenvalue_lower_exact"])*Fraction(r["norm_R_exact"])==1-Fraction(r["z_r34_exact"])
               and Fraction(r["strict_cube_interior_margin_exact"])+Fraction(r["localized_posterior_exact"])==R34 for r in rows)
    plot(rows,order)
    output(rows,doc)
    links={fn:sha(HERE/fn) for fn in ("MINIMA.md","MINIMA.html","cached_minima.png","bounds.csv")}
    receipt={"status":"PASS_EXACT_INPUTS_AND_RENDERED_PLOT_SERIES","independent_minima_sha256":RECEIPT,
             "observed_readback_exit_sha256":OS,"build_source_sha256":sha(Path(__file__)),
             "artifacts":links,"cases":36,"plotted_scalar_values":72,"new_mechanical_evaluations":0}
    write(HERE/"build_receipt.json",json.dumps(receipt,indent=2)+"\n")
    print(sha(HERE/"build_receipt.json"))


if __name__=="__main__":
    main()
