"""Render the independently verified four exact-spatial observer comparisons."""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from html import escape
from importlib.util import module_from_spec, spec_from_file_location
import json
import os
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


HERE=Path(__file__).resolve().parent
EXACT=HERE.parent
DESIGN=EXACT.parent
RUN=DESIGN.parent
ROOT=RUN.parents[1]
SITE=ROOT.parent
ASSETS=SITE/'docs/assets/katex'
AUTH=HERE/'independent_actual_readback/receipt.json'
OS=HERE/'independent_actual_readback/observed_os_exit.json'
SPATIAL=EXACT/'four_spatial_checkpoint_v1/SPATIAL.html'
FOUR=EXACT/'bridge_reduction_v1/independent_actual_preparation/all_four_readback.json'
CACHED=DESIGN/'full_outcome_coverage_v1/observer_checkpoint_v1/OBSERVERS.html'
STUDY=RUN/'report/REPORT.html'
SWEEP=RUN/'observer_robustness/independent_actual/receipt.json'
HELPER=EXACT/'certification_checkpoint_v1/build_checkpoint.py'
CONVERT=ROOT/'scripts/build_report.py'
FAVICON=RUN/'report/favicon.png'
FIGURE=HERE/'exact_spatial_contrast.png'
TITLE='Observer bounds at four exact-spatial roots'
PLUS,MINUS='max_contrast_change','max_negative_contrast_change'
VARIANTS=('current_disk__spatial_max','current_disk__regional_mean',
          'current_radial_line__spatial_max','current_radial_line__regional_mean')
LABELS=('Disk · component maxima','Disk · component means','Radial line · component maxima','Radial line · component means')
PINS={
    AUTH:'e57c6e18e4e585b002ebab8aa0782d063e281c201ab1449cfac484ef325287e7',
    OS:'285d93d875619ed1f4e82457320a80a4dfd399be29531a843f0f632898053c4c',
    AUTH.with_name('readback.py'):'a5b8484c26303022e216675580f5bc7bf491023136a885c03b9abd373f9503c6',
    HERE/'results.json':'df5c68d82b27748c5c82a883e963df8e55b47fc43ad64e54c9957afb0e216c9c',
    HERE/'observed_os_exit.json':'e70cfb92af9f47dc44d6b1682ec520a663e032850d68c59ada456b64245526c1',
    SPATIAL:'6dd9100643191f6617f5032d2a9d4e35ef3e179c9826115c841b041ce92c9e4d',
    FOUR:'859f507ad3762cfef44d53672bbe15a43592ce54d6e08ba8828bf50a925ec340',
    CACHED:'dc7e2aa40c2511944f4ba0e226da1da844d866c85cf436c165ffa0e6a78ae557',
    STUDY:'d33f3d594406ed498d364c3462450b62cb27054f0146f4bbca3a1c2ae22b41f6',
    SWEEP:'703b574bc51a083b10c7328d9622862e939ddabe791b18ce9e1817f07718f69f',
    HELPER:'50524d01db5003b0839e873efb37d41bb0552ad0f0c162f73bd668ad0e58d31b',
    CONVERT:'90148d53d74d8ea5cc8a665e6a3660abf2db72487a1fa59586caba98f9cc38e2',
    FAVICON:'961e6e85ec16258652bf3c3f322d336fb598b90bccc49821047babfef3925530',
    ASSETS/'katex.min.css':'5a0ff74f11887f7260dc16b216906f51352b84ca29f21abb9a14b62a2ca6eb15',
    ASSETS/'katex.min.js':'e6bfe5deebd4c7ccd272055bab63bd3ab2c73b907b6e6a22d352740a81381fd4',
}


def check(cond,msg):
    if not cond:
        raise ValueError(msg)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def link(path):
    return Path(os.path.relpath(path,HERE)).as_posix()


def parsed(path):
    return json.loads(path.read_bytes())


def directed(x,places,up):
    factor=10**places
    v=x*factor
    n=v.__ceil__() if up else v.__floor__()
    s='−' if n<0 else ''
    a,b=divmod(abs(n),factor)
    return f'{s}{a}.{b:0{places}d}'


def exact_interval(w):
    lo,hi=map(Q,w)
    check(lo<=hi,'ordered exact endpoints required')
    return lo,hi


def numbers():
    for path,h in PINS.items():
        check(sha(path)==h,f'changed evidence, report or renderer: {path}')
    a,o,f,s=map(parsed,(AUTH,OS,FOUR,SWEEP))
    check(a['status']=='PASS_INDEPENDENT_EXACT_FOUR_SPATIAL_OBSERVER_READBACK'
          and a['independently_matched_spatial_cases']==4 and a['independently_matched_L_alternative_comparisons']==32
          and a['producer_functions_executed']==0 and a['all_136_primary_comparisons_claimed'] is False
          and a['producer_results_sha256']==o['producer_results_sha256']==PINS[HERE/'results.json']
          and a['producer_observed_os_exit_sha256']==o['producer_observed_os_exit_sha256']==PINS[HERE/'observed_os_exit.json']
          and o['status']=='PASS_AGENT_OBSERVED_INDEPENDENT_EXACT_FOUR_OBSERVER_EXIT'
          and type(o['observed_exit_code']) is int and o['observed_exit_code']==0
          and o['source_sha256']==a['independent_source_sha256']==PINS[AUTH.with_name('readback.py')]
          and o['native_readback_sha256']==PINS[AUTH] and a['spatial_proof_closure_sha256']==PINS[FOUR]
          and [r['spatial_root_proof_sha256'] for r in a['cases']]==[r['independent_proof_sha256'] for r in f['cases']]
          and len(f['cases'])==4 and f['full_five_gate_goal_achieved'] is False,
          'incorrect proof, completed readback, or scope')
    large={r['endpoint_id']:exact_interval(r['exact_spatial_root_Q_um_exact']) for r in a['cases'] if r['role']=='large_final'}
    plus,minus=(large[k] for k in (PLUS,MINUS))
    gap=plus[0]-minus[1],plus[1]-minus[0]
    check(plus[1]<0 and minus[1]<0 and gap[0]>0 and len(large)==2
          and gap==exact_interval(a['large_signed_gap_um_exact']) and a['strict_positive_gap'] is True,
          'signed gap or individual sign mismatch')
    small=[(p['endpoint_id'],r) for p in a['cases'] if p['role']=='small_record' for r in p['comparisons']]
    counts={k:Counter(r['status'] for _,r in small if r['variant']==k) for k in VARIANTS}
    delta=Q(a['conditional_allowance_um_exact'])
    minimums={}
    for endpoint,row in small:
        lo,hi=exact_interval(row['absolute_residual_um_exact'])
        check((row['status']=='IN' and hi<=delta) or (row['status']=='OUT' and lo>delta),
              'wrong bounded conditional comparison')
        if row['status']=='OUT':
            check(row['variant'] in VARIANTS[2:] and row['observable']=='region_averaged_shortening_um',
                  'radial-line mean interpretation changed')
            minimums.setdefault(endpoint,set()).add(lo)
    check(len(small)==32 and a['comparison_counts']=={'IN':28,'OUT':4}
          and list(counts.values())==[{'IN':8},{'IN':8},{'IN':6,'OUT':2},{'IN':6,'OUT':2}]
          and all(len(x)==1 for x in minimums.values()) and set(minimums)=={PLUS,MINUS}
          and delta==Q(s['source_configured_delta_um_exact']) and s['robust_across_conventions_at_frozen_delta'] is False,
          'wrong primary/alternative grouping or conditional allowance')
    return plus,minus,gap,counts,delta,{k:next(iter(v)) for k,v in minimums.items()}


def make_figure(plus,minus,gap):
    bg,fg,muted,grid='#0d1927','#eaf3fc','#b6cada','#32475d'
    cold,warm='#80d6ea','#efa781'
    p,n=(float((min(pair)+max(pair))/2) for pair in (plus,minus))
    fig,ax=plt.subplots(figsize=(8.3,3.45),dpi=185,facecolor=bg)
    ax.set_facecolor(bg)
    ax.grid(axis='x',color=grid,linewidth=.7)
    ax.set_axisbelow(True)
    ax.axhline(0,color='#526b82',lw=.8)
    ax.plot([n,p],[0,0],lw=2.5,color=fg,solid_capstyle='round')
    ax.scatter([n,p],[0,0],s=80,c=[warm,cold],zorder=3)
    ax.annotate('+ direction',xy=(p,0),xytext=(p,.21),ha='center',color=cold,fontsize=11,weight='medium')
    ax.annotate('− direction',xy=(n,0),xytext=(n,.21),ha='center',color=warm,fontsize=11,weight='medium')
    ax.text(n,-.20,f'{n:.6f}'.replace('-','−'),ha='center',color=warm,fontsize=10)
    ax.text(p,-.20,f'{p:.6f}'.replace('-','−'),ha='center',color=cold,fontsize=10)
    ax.annotate(f'+{float(sum(gap)/2):.6f} µm',xy=((p+n)/2,.025),ha='center',va='bottom',color=fg,
                fontsize=12,weight='medium',bbox=dict(facecolor=bg,edgecolor='none',pad=3))
    ax.set_ylim(-.32,.40)
    ax.set_xlim(-4.36,-4.06)
    ax.set_yticks([])
    ax.xaxis.set_major_locator(MultipleLocator(.06))
    ax.tick_params(axis='x',colors=muted)
    ax.set_xlabel('Large-region apical contrast Q (µm)',color=muted,labelpad=8)
    for side in ('top','left','right'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color(grid)
    fig.text(.065,.90,'Both contrasts are negative',color=fg,fontsize=19,weight='semibold')
    fig.text(.065,.799,'The plus candidate has the higher apical contrast',color=muted,fontsize=11)
    fig.subplots_adjust(left=.095,right=.94,bottom=.225,top=.68)
    with FIGURE.open('xb') as handle:
        fig.savefig(handle,format='png',dpi=185,facecolor=bg)
    plt.close(fig)


BODY=r"""# Observer bounds at four exact-spatial roots

At the matched large input, the two frozen finite-element models have a positive signed difference in the apical contrast. The independently verified gap is **[[GAP_LO]] to [[GAP_HI]] µm**. Both individual contrasts are negative. The two verified small-L roots account for eight primary comparisons and 24 alternative comparisons. They provide **28 IN and 4 OUT** at the frozen allowance of approximately **[[DELTA]] µm**.

## Same sign, different magnitude

The large observable is \(Q=Z_a(0)-Z_a(R_s)\), with the axis at zero and the shoulder at current radius \(R_s=24.94698765\,\mu\mathrm m\) in both fields. A negative value means the axis lies below the shoulder in this axial coordinate. The labels + and − identify the two original affine search directions. The contrast result gives a difference in magnitude and no evidence of opposite folding directions or biological significance.

![On a magnified range of the negative Q axis, the minus candidate is at minus 4.290758 micrometers and the plus candidate is at minus 4.136003 micrometers. The difference is about plus 0.154755 micrometers.](exact_spatial_contrast.png)

| Large-state value | Lower bound (µm) | Upper bound (µm) |
|:--|--:|--:|
[[LARGE_TABLE]]

The interval widths are smaller than the plot resolution. The numeric bounds round outward. The signed difference adds the uncertainty from both root observables through interval subtraction:

\[
[Q_{+,l}-Q_{-,u},\ Q_{+,u}-Q_{-,l}] .
\]

## Comparison at the two small L states

Each line includes both candidates and four small-L observables. IN requires the entire residual interval to stay at or below the allowance. OUT requires its lower bound to exceed the allowance.

| Static convention | IN | OUT |
|:--|--:|--:|
[[SMALL_TABLE]]

All four OUT rows are the regional mean of total shortening under radial-line weighting. This mean is shared by the two radial-line conventions. The exact residual lower bounds at L are **[[LOW_PLUS]] µm** for + and **[[LOW_MINUS]] µm** for −. The allowance is a conditional model discrepancy threshold. The experiment’s measurement uncertainty, time aggregation and reference registration remain open.

## How the root certificate reaches an observation

Let \(e\) bound the free-coordinate distance of the spatial root from the exact binary center and let \(\xi\) bound the coordinate addition error in the saved field. The independent readback reconstructed both quadratic traces from the stored nodes, confirmed fixed axis and rim radii, and checked the radial derivative against \(8\epsilon\) on every edge for the new \(\epsilon\). A face with stored slope bound \(M_s\) then has Eulerian height error \(E_s\) given by

\[
\epsilon=33.5(e+\xi)\,\mu\mathrm m,
\qquad E_s=\tfrac54(1+M_s)\epsilon .
\]

For total shortening, the two face errors add. A component maximum or mean uses its own face bound. The axis-node error \(\epsilon\) adds to the apical shoulder error \(E_a\) for \(Q\). This retains the fixed observer, including evaluation at the same current radius as the curve changes.

## What the gap would imply under full compatibility

If both models belong to the genuinely calibration-compatible set, every single-point prediction \(q\) for their large observable incurs worst-case absolute error at least half their gap. The triangle inequality gives the classical two-value bound

\[
\max\{|q-Q_+|,|q-Q_-|\}\ge\tfrac12|Q_+-Q_-| .
\]

With these two Q bounds that conditional lower bound would be **at least [[HALF_LO]] µm**. The four spatial certificates alone cover only eight of 136 primary comparisons. The other 128 and their associated 32 exact-spatial roots remain unverified at this scope. The separate 36-state result uses cached quadrature. Source compatibility, a stable loading path, continuum adequacy and scientific novelty remain unresolved.

## Evidence

- [Independent four-root observer receipt]([[AUTH]]) · [Observed independent replay exit]([[OS]])
- [Four exact-spatial root and minimum certificates]([[SPATIAL]]) · [Spatial proof readback]([[FOUR]])
- [Full 36-field cached-quadrature observer checkpoint]([[CACHED]]) · [Nonlinear stored-field report]([[STUDY]]) · [Independent stored-field convention sweep]([[SWEEP]])
"""


def markdown(plus,minus,gap,counts,delta,fail):
    def row(label,values):
        lo,hi=values
        return f'| {label} | {directed(lo,12,False)} | {directed(hi,12,True)} |'
    paths={'AUTH':AUTH,'OS':OS,'SPATIAL':SPATIAL,'FOUR':FOUR,'CACHED':CACHED,'STUDY':STUDY,'SWEEP':SWEEP}
    terms={**{k:link(v) for k,v in paths.items()},'GAP_LO':directed(gap[0],14,False),
           'GAP_HI':directed(gap[1],14,True),'HALF_LO':directed(gap[0]/2,14,False),'DELTA':f'{float(delta):.9f}',
           'LARGE_TABLE':'\n'.join((row('Candidate +',plus),row('Candidate −',minus),row('Signed gap',gap))),
           'SMALL_TABLE':'\n'.join(f"| {label} | {counts[k]['IN']} | {counts[k]['OUT']} |"
                                   for label,k in zip(LABELS,VARIANTS,strict=True)),
           'LOW_PLUS':directed(fail[PLUS],8,False),'LOW_MINUS':directed(fail[MINUS],8,False)}
    text=BODY
    for k,v in terms.items():
        text=text.replace(f'[[{k}]]',v)
    check('[[' not in text,'unfilled report parameter')
    return text


CSS='''
*{box-sizing:border-box}html{color-scheme:dark}body{margin:0;background:#091422;color:#e5f0fb;font:1rem/1.65 system-ui,sans-serif}
main{width:min(100%,52.5rem);min-width:0;margin:auto;padding:2rem clamp(.95rem,3vw,2.4rem) 4rem}
nav{display:flex;flex-wrap:wrap;gap:.22rem 1.05rem;margin:0 0 1.55rem;font-size:.91rem}
h1,h2{line-height:1.18;letter-spacing:-.025em;text-wrap:balance}h1{font-size:clamp(2rem,5.5vw,3.1rem);margin:.6rem 0 1rem}
h2{font-size:clamp(1.37rem,4vw,1.85rem);margin:2rem 0 .8rem}h1+p{font-size:1.12rem;color:#cbdeed}
a{color:#8ae5fa;text-underline-offset:.18em}a:focus-visible{outline:2px solid #ffd0a8;outline-offset:3px}
img{display:block;width:100%;height:auto;border-radius:8px}li{margin:.42rem 0}
table{width:100%;border-collapse:collapse;font-size:.91rem;text-align:right;font-variant-numeric:tabular-nums}
td,th{padding:.56rem .58rem;border-bottom:1px solid #36536d}th{color:#c2dcec}td:first-child,th:first-child{text-align:left}
.math.display{margin:1.1rem 0;max-width:100%;overflow-x:auto}.katex-display{margin:0}.math{color:#f1f7ff}
@media(max-width:500px){main{padding-top:1.08rem}td,th{padding:.44rem .23rem;font-size:.76rem;line-height:1.3}
.math.display .katex{font-size:.93em}ul{padding-left:1.2rem}}
'''


def render(text):
    spec=spec_from_file_location('pinned_static_katex',HELPER)
    helper=module_from_spec(spec)
    spec.loader.exec_module(helper)
    sys.path.insert(0,str(ROOT/'scripts'))
    import build_report
    check(Path(build_report.__file__).resolve()==CONVERT,'wrong Markdown converter')
    body,count=helper.compile_math(build_report.convert(text))
    page=('<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
          '<meta name="viewport" content="width=device-width,initial-scale=1">'
          f'<title>{escape(TITLE)}</title><link rel="icon" type="image/png" href="{link(FAVICON)}">'
          f'<link rel="stylesheet" href="{link(ASSETS/"katex.min.css")}"><style>{CSS}</style></head>'
          '<body><main><nav aria-label="Related reports">'
          f'<a href="{link(SPATIAL)}">Spatial certificates</a><a href="{link(CACHED)}">Cached observer results</a>'
          '<a href="OBSERVERS.md">Markdown</a></nav>'+body+'</main></body></html>\n')
    parser=helper.LocalLinks(); parser.feed(page)
    paths={}
    for v in parser.paths:
        target=(HERE/v).resolve()
        check(target.is_relative_to(SITE) and (target.is_file() or target==HERE/'OBSERVERS.md'),f'broken local link: {v}')
        paths[v]=target.relative_to(SITE).as_posix()
    return page,count,paths


def save(path,data):
    with path.open('x',encoding='utf-8',newline='\n') as f:
        f.write(data)


def main():
    plus,minus,gap,counts,delta,fail=numbers()
    make_figure(plus,minus,gap)
    md=markdown(plus,minus,gap,counts,delta,fail)
    html,math,paths=render(md)
    save(HERE/'OBSERVERS.md',md); save(HERE/'OBSERVERS.html',html)
    receipt={'status':'PASS_BUILD_FROM_INDEPENDENT_EXACT_SPATIAL_OBSERVER_READBACK',
             'source_sha256':sha(Path(__file__)),'rendered_math_count':math,'local_links':paths,
             'input_sha256':{p.relative_to(SITE).as_posix():h for p,h in PINS.items()},
             'file_sha256':{n:sha(HERE/n) for n in ('OBSERVERS.md','OBSERVERS.html',FIGURE.name)},
             'signed_comparison_um_exact':{k:list(map(str,b)) for k,b in zip((PLUS,MINUS,'gap'),(plus,minus,gap),strict=True)},
             'conditional_half_gap_lower_exact':str(gap[0]/2),'conditional_model_allowance_um_exact':str(delta),
             'observer_comparison_counts':counts,'primary_comparisons_covered':8,'alternative_comparisons_covered':24,
             'source_compatibility_or_136_exact_spatial_comparisons_claimed':False}
    path=HERE/'observer_report_build_receipt.json'
    save(path,json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'build_receipt_sha256':sha(path),'math':math,'files':receipt['file_sha256']}))


if __name__=='__main__':
    main()
