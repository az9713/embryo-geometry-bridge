"""Build the independently joined whole-box geometry and cached derivative snapshot."""
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
from matplotlib.lines import Line2D


HERE=Path(__file__).resolve().parent
EXACT=HERE.parent
DESIGN=EXACT.parent
RUN=DESIGN.parent
ROOT=RUN.parents[1]
SITE=ROOT.parent
ASSETS=SITE/'docs/assets/katex'
JOIN=EXACT/'remaining32_bridge_v1/admission/bridge_input_admission.json'
GEO=EXACT/'remaining32_geometry_v1/independent_actual_preparation/actual_replay/receipt.json'
VAR=EXACT/'remaining32_derivative_extension_v1/independent_actual_preparation/actual_replay/receipt.json'
TARGET=EXACT/'remaining32_targets_v1/independent_prelaunch/review.json'
SPATIAL=EXACT/'four_spatial_checkpoint_v1/SPATIAL.html'
OBSERVER=EXACT/'observer_application_v2/OBSERVERS.html'
HELPER=EXACT/'certification_checkpoint_v1/build_checkpoint.py'
CONVERT=ROOT/'scripts/build_report.py'
FAVICON=RUN/'report/favicon.png'
FIGURE=HERE/'geometry_and_cached_defect.png'
TITLE='Verified inputs for 32 spatial integrations'
ENDPOINTS=('max_contrast_change','max_negative_contrast_change')
ORDER=tuple(f'{p}_small_{i:02d}' for p in ENDPOINTS for i in range(1,17))
PINS={
    JOIN:'2e2e388e798ed509bb52bc72a6169b8b6ee5db55173b3529ac73c3812b23a6cb',
    JOIN.with_name('os_exit_observation.json'):'636d312911742e22722377ee540e13ed725c116de6d1853f137349d03b9dee96',
    GEO:'74179dde74571f792b6f3949a06d6e26076b8c746e6700c355264d0de3debf04',
    GEO.with_name('os_exit_observation.json'):'d39b2d5e97ad6d15bfb68f3f8360e0214c45d11b37dbbd3eaae48313a2e8fc0c',
    VAR:'2a1be56b3050215e05da8e11b4785a65cdc8a229210fa82ebee37f1aa8236d05',
    VAR.with_name('os_exit_observation.json'):'69ededa4a104e07ac85e251e6a019b93fc11016cdd83954aaee9f7a987f83325',
    TARGET:'67d33bbc3171839b1c578b02b5cfbf9ea6ba5d6cb351d0655dd14aaf99df75f9',
    SPATIAL:'6dd9100643191f6617f5032d2a9d4e35ef3e179c9826115c841b041ce92c9e4d',
    OBSERVER:'1a0e24c8369db76af550e7182ed87310e6125f78f253f197d15f3c08d4efa6eb',
    HELPER:'50524d01db5003b0839e873efb37d41bb0552ad0f0c162f73bd668ad0e58d31b',
    CONVERT:'90148d53d74d8ea5cc8a665e6a3660abf2db72487a1fa59586caba98f9cc38e2',
    FAVICON:'961e6e85ec16258652bf3c3f322d336fb598b90bccc49821047babfef3925530',
    ASSETS/'katex.min.css':'5a0ff74f11887f7260dc16b216906f51352b84ca29f21abb9a14b62a2ca6eb15',
    ASSETS/'katex.min.js':'e6bfe5deebd4c7ccd272055bab63bd3ab2c73b907b6e6a22d352740a81381fd4',
}


def require(condition,message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def link(path):
    return Path(os.path.relpath(path,HERE)).as_posix()


def canon(name):
    p=(ROOT/name).resolve()
    require(p.is_relative_to(ROOT) and p.relative_to(ROOT).as_posix()==name,'invalid study reference')
    return p


def parsed(path):
    return json.loads(path.read_bytes())


def dataset():
    for path,h in PINS.items():
        require(sha(path)==h,f'changed fixed source or proof: {path.name}')
    join,geo,var,t=map(parsed,(JOIN,GEO,VAR,TARGET))
    require(join['status']=='PASS_INDEPENDENT_REMAINING32_SAME_TARGET_BRIDGE_INPUTS'
            and join['spatial_workers_authorized'] is False
            and geo['status']=='PASS_INDEPENDENT_ACTUAL_REMAINING32_SINGLETON_GEOMETRY'
            and geo['successfully_replayed_geometry_cases']==32 and geo['root_asserted'] is False
            and var['status']=='PASS_INDEPENDENT_ACTUAL_TEN_FRESH_CACHED_DERIVATIVES'
            and t['status']=='PASS_INDEPENDENT_REMAINING32_TARGET_CONTEXTS'
            and join['geometry_receipt_sha256']==PINS[GEO] and join['variation_receipt_sha256']==PINS[VAR]
            and join['geometry_observed_exit_sha256']==PINS[GEO.with_name('os_exit_observation.json')]
            and join['variation_observed_exit_sha256']==PINS[VAR.with_name('os_exit_observation.json')]
            and tuple(c['case_id'] for c in join['cases'])==tuple(r['pilot_id'] for r in geo['geometry_records'])
            ==tuple(t['ordered_case_ids'])==ORDER,'complete proof join and verified process scope required')
    for receipt in (JOIN,GEO,VAR):
        status=parsed(receipt.with_name('os_exit_observation.json'))
        require(type(status['observed_exit_code']) is int and status['observed_exit_code']==0,
                'observed true serial completion required')
    links=join['source_and_input_sha256']
    fresh={r['case_id']:r for r in var['case_results']}
    rows=[]
    for shared,gr,ctx in zip(join['cases'],geo['geometry_records'],t['per_case'],strict=True):
        case=shared['case_id']; kind=shared['cached_input_kind']
        file=canon(gr['receipt'])
        require(gr['sha256']==links[gr['receipt']]==shared['geometry_replay_case_sha256']==sha(file),
                'one geometry case proof changed')
        g=parsed(file)
        j,z,n=(Q(shared[k]) for k in ('geometry_J_lower_exact','reviewed_ZQ_upper_exact','R_inf_exact'))
        require(gr['status']==g['status']=='PASS_INDEPENDENT_REMAINING32_WHOLE_BOX_GEOMETRY_ONLY'
                and case==gr['pilot_id']==g['pilot_id']==ctx['case_id']
                and Q(gr['radius_exact'])==Q(g['radius_exact'])==Q(1,2**34)
                and 0<j==Q(g['body_J_bounds_exact'][0]) and 0<=z<1 and n>0
                and shared['trial_box_sha256']==g['box_sha256']
                and Q(g['strong_monotonicity_lower_exact'])>0 and Q(g['min_surface_G_exact'])>0
                and Q(g['min_trace_hoop_exact'])>0 and Q(g['min_site_tangent_squared_exact'])>0
                and g['all_144_body_and_24_trace_domains'] is True and g['root_asserted'] is False
                and g['raw_surface_cache_geometry_lines']==288 and g['closed_full_nodal_coordinates_rebuilt']==650
                and Q(g['archived_shift_exact'])==Q(ctx['observer_shift_sigma_exact'])
                and Q(g['archived_addition_gap_exact'])==Q(ctx['observer_archived_xi_exact']),
                'geometry or archive authority not independently validated')
        if kind=='FRESH_R34':
            d=fresh[case]
            require(d['status']=='PASS_REVIEWED_FRESH_R34_CACHED_DERIVATIVE'
                    and d['native_sha256']==shared['cached_derivative_case_sha256']==links[d['native_file']]==sha(canon(d['native_file']))
                    and Q(d['reviewed_ZQ_upper_exact'])==z, 'fresh cached derivative mismatch')
        else:
            require(kind=='NATIVE_R34' and case not in fresh and Q(ctx['native_cached_radius_exact'])==Q(1,2**34)
                    and Q(ctx['native_ZQ_upper_exact'])==z, 'native cached derivative scope mismatch')
        rows.append({'id':case,'input_type':kind,'J_lower':j,'Z_Q_upper':z,'R_inf':n,
                     'strong_monotonicity_lower':Q(g['strong_monotonicity_lower_exact']),
                     'bulk_rho_upper':Q(g['max_bulk_rho_outward_exact']),
                     'trace_rho_upper':Q(g['max_surface_rho_raw_exact']),
                     'archived_sigma':Q(g['archived_shift_exact'])})
    require(Counter(r['input_type'] for r in rows)=={'FRESH_R34':10,'NATIVE_R34':22}
            and {r['id'] for r in rows if r['archived_sigma']>0}
            =={'max_contrast_change_small_13','max_negative_contrast_change_small_16'},'wrong scope counts or recentered cases')
    return rows


def numbers(rows):
    return {'min_J':min(r['J_lower'] for r in rows),'min_monotonicity':min(r['strong_monotonicity_lower'] for r in rows),
            'max_bulk_rho':max(r['bulk_rho_upper'] for r in rows),'max_trace_rho':max(r['trace_rho_upper'] for r in rows),
            'min_ZQ':min(r['Z_Q_upper'] for r in rows),'max_ZQ':max(r['Z_Q_upper'] for r in rows),
            'min_headroom':min(1-r['Z_Q_upper'] for r in rows)}


def decimal(x,places,up=False):
    factor=10**places
    f=x*factor; n=f.__ceil__() if up else f.__floor__()
    a,b=divmod(n,factor)
    return f'{a}.{b:0{places}d}'


def plot(rows):
    bg,fg,muted,grid='#0d1927','#e9f3fb','#b5cadd','#30465d'
    colors=('#83d7e9','#efa885')
    fig,axs=plt.subplots(2,1,sharex=True,figsize=(8.6,5.9),dpi=175,facecolor=bg)
    for i,(ax,y,title) in enumerate(zip(axs,('J_lower','Z_Q_upper'),
                                       (r'$J$ lower bound',r'$\bar Z_Q$ upper bound'),strict=True)):
        ax.set_facecolor(bg)
        ax.set_axisbelow(True)
        ax.yaxis.grid(color=grid,linewidth=.7)
        ax.tick_params(colors=muted,labelsize=9)
        ax.set_ylabel(title,color=fg,labelpad=10,fontsize=9)
        for side in ('left','right','top'):
            ax.spines[side].set_visible(False)
        ax.spines['bottom'].set_color(grid)
        for e,(kind,col,marker) in enumerate(zip(ENDPOINTS,colors,('o','s'),strict=True)):
            subset=[(int(r['id'][-2:]),r) for r in rows if r['id'].startswith(kind+'_')]
            for new in (True,False) if i else (False,):
                shown=[(k,r) for k,r in subset if i==0 or (r['input_type']=='FRESH_R34')==new]
                ax.scatter([k+(-.12 if e==0 else .12) for k,_ in shown],[float(r[y]) for _,r in shown],
                           s=34,marker=marker,edgecolors=col,facecolors='none' if new else col,linewidths=1.25)
        ax.set_ylim(.957,1.001) if i==0 else ax.set_ylim(.203,.317)
    axs[1].set_xticks(range(1,17),[f'{i:02d}' for i in range(1,17)])
    axs[1].set_xlabel('Small input index within each candidate',color=muted,labelpad=8)
    handles=[Line2D([],[],marker=m,linestyle='none',markerfacecolor=c,markeredgecolor=c,label=t)
             for m,c,t in zip(('o','s'),colors,('Candidate +','Candidate −'),strict=True)]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.067,.87),ncol=2,frameon=False,labelcolor=fg,fontsize=9)
    fig.text(.068,.925,'Geometry and cached derivative bounds',color=fg,fontsize=19,weight='semibold')
    fig.text(.068,.87,'Open markers in the lower panel indicate fresh derivative bounds.',color=muted,fontsize=9)
    fig.subplots_adjust(left=.16,right=.945,top=.79,bottom=.145,hspace=.24)
    with FIGURE.open('xb') as f:
        fig.savefig(f,format='png',facecolor=bg)
    plt.close(fig)


BODY=r"""# Verified inputs for 32 spatial integrations

The remaining 32 fixed-input states now have independently replayed whole-box geometry and cached derivative bounds on the same closed free-displacement radius \(r=2^{-34}\). The geometry covers all 32 trial boxes. Ten cached derivative bounds were recomputed at the matching radius, and 22 already-reviewed bounds match the same trial cubes.

![Two panels show 32 input bounds, indexed 01–16 for each of two candidates. Every reported whole-box J lower bound exceeds 0.961. Every cached derivative defect upper bound is below 0.309. Open markers identify 10 fresh cached derivative bounds.](geometry_and_cached_defect.png)

The trial boxes use each case’s fixed exact material values, gain times dose and frozen binary preconditioner. The plotted input index labels records. The points do not track a loading history.

## Whole-box geometry

Each case used an admissible fixed-reference quadratic finite-element mesh. The reviewer independently reconstructed all 650 exact nodal coordinate intervals, checked 144 bulk charts and 24 traces, and replayed the Bernstein domain bounds. Across the 32 cases, the reported lower bounds on the deformation Jacobian \(J\) are all at least **[[J]]**. The certified lower bound on strong monotonicity is at least **[[M]]** in the frozen coordinate scaling. Bounds on the bulk and trace series-domain ratios are **[[BULK]]** and **[[TRACE]]** respectively. These are whole domain checks for the stated finite-element model. They place no bound on continuum error.

## Cached derivative bounds

With cached quadrature Hessian \(H_Q\) and the recorded preconditioner \(R\), the stated quantity is a whole-cube upper bound

\[
\bar Z_Q\ge\sup_{u\text{ in cube}}\|I-RH_Q(u)\|_\infty .
\]

| Cached derivative route | Cases | Matching evidence |
|:--|--:|:--|
| Existing exact-radius bound | 22 | Native bound already covers this trial cube. |
| Fresh exact-radius bound | 10 | Seven smaller trial cubes and three enlarged from radius \(2^{-36}\). |

The largest reported \(\bar Z_Q\) is at most **[[ZMAX]]**, leaving at least **[[ROOM]]** below one in this cached derivative bound. The independently reviewed ten-case extension and the original 22-case bounds enter one same-target manifest. Both recentered states (+ small 13 and − small 16) retain their exact shifts \(\sigma\) from the archived observer fields and the stored-coordinate addition gaps \(\xi\). If \(e\) is a verified spatial root error relative to the proof center, the coordinate difference against the originally saved observer field is bounded by

\[
\epsilon_{\text{nodes}}\le 33.5(\sigma+e+\xi)\,\mu\mathrm m .
\]

## Remaining proof obligations

The spatial root test still needs a signed enclosure of the center residual \(RF_I(u_0)\) to obtain its norm bound \(\eta_I\). It also needs a whole-cube norm bound \(b\) on \(R(H_I-H_Q)\). With \(z_I=\bar Z_Q+b\), the sufficient test requires

\[
z_I<1,\qquad \eta_I+z_I r<r .
\]

The current geometry and cached derivative evidence does not supply either spatial defect. Minimum proofs additionally require the applicable positive symmetric-part premise for \(R\). The earlier four exact-spatial certificates remain distinct. Full exact-spatial outcome coverage, source compatibility, loading-path validity and continuum adequacy remain open.

## Evidence

- [32-case geometry replay]([[GEO]]) · [Ten fresh derivative replays]([[VAR]]) · [Joined input authority]([[JOIN]])
- [Exact target-context review]([[TARGET]]) · [Earlier four spatial proofs]([[SPATIAL]]) · [Earlier four spatial observer bounds]([[OBSERVER]])
"""


def write(path,text):
    with path.open('x',encoding='utf-8',newline='\n') as f:
        f.write(text)


def markdown(n):
    vals={**{k:link(v) for k,v in (('GEO',GEO),('VAR',VAR),('JOIN',JOIN),('TARGET',TARGET),('SPATIAL',SPATIAL),('OBSERVER',OBSERVER))},
          'J':decimal(n['min_J'],6),'M':decimal(n['min_monotonicity'],6),'BULK':decimal(n['max_bulk_rho'],6,True),
          'TRACE':decimal(n['max_trace_rho'],6,True),'ZMAX':decimal(n['max_ZQ'],6,True),'ROOM':decimal(n['min_headroom'],6)}
    text=BODY
    for k,v in vals.items():
        text=text.replace(f'[[{k}]]',v)
    require('[[' not in text,'empty document slot')
    return text


CSS='''
*{box-sizing:border-box}html{color-scheme:dark}body{margin:0;background:#091422;color:#e5f0fb;font:1rem/1.65 system-ui,sans-serif}
main{width:min(100%,55rem);min-width:0;margin:auto;padding:2rem clamp(.95rem,3vw,2.4rem) 4rem}
nav{display:flex;flex-wrap:wrap;gap:.25rem 1.15rem;margin:0 0 1.6rem;font-size:.91rem}
h1,h2{line-height:1.2;letter-spacing:-.025em;text-wrap:balance}h1{font-size:clamp(2rem,5.5vw,3.15rem);margin:.6rem 0 1rem}
h2{margin:2rem 0 .8rem;font-size:clamp(1.38rem,4vw,1.84rem)}h1+p{font-size:1.12rem;color:#c7dbea}
a{color:#8ae5fa;text-underline-offset:.18em}a:focus-visible{outline:2px solid #ffd0a8;outline-offset:3px}
img{display:block;width:100%;height:auto;border-radius:8px}table{width:100%;border-collapse:collapse;font-size:.91rem}
td,th{padding:.54rem .56rem;border-bottom:1px solid #36536d}th{color:#c2dcec}li{margin:.42rem 0}
.math.display{margin:1.1rem 0;max-width:100%;overflow-x:auto}.katex-display{margin:0}.math{color:#f1f7ff}
@media(max-width:500px){main{padding-top:1.05rem}td,th{padding:.44rem .25rem;font-size:.77rem;line-height:1.4}ul{padding-left:1.2rem}}
'''


def render(text):
    spec=spec_from_file_location('fixed_katex_renderer',HELPER)
    helper=module_from_spec(spec); spec.loader.exec_module(helper)
    sys.path.insert(0,str(ROOT/'scripts'))
    import build_report
    require(Path(build_report.__file__).resolve()==CONVERT,'wrong Markdown converter')
    body,count=helper.compile_math(build_report.convert(text))
    page=('<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
          '<meta name="viewport" content="width=device-width,initial-scale=1">'
          f'<title>{escape(TITLE)}</title><link rel="icon" type="image/png" href="{link(FAVICON)}">'
          f'<link rel="stylesheet" href="{link(ASSETS/"katex.min.css")}"><style>{CSS}</style></head>'
          '<body><main><nav aria-label="Related reports">'
          f'<a href="{link(SPATIAL)}">Four spatial proofs</a><a href="{link(OBSERVER)}">Spatial observer bounds</a>'
          '<a href="INPUTS.md">Markdown</a></nav>'+body+'</main></body></html>\n')
    links=helper.LocalLinks(); links.feed(page)
    linked={}
    for url in links.paths:
        path=(HERE/url).resolve()
        require(path.is_relative_to(SITE) and (path.is_file() or path==HERE/'INPUTS.md'),f'broken local evidence link {url}')
        linked[url]=path.relative_to(SITE).as_posix()
    return page,count,linked


def main():
    rows=dataset(); n=numbers(rows)
    plot(rows)
    md=markdown(n); page,count,links=render(md)
    write(HERE/'INPUTS.md',md); write(HERE/'INPUTS.html',page)
    receipt={'status':'PASS_BUILD_FROM_INDEPENDENT_REMAINING32_JOINED_INPUTS','source_sha256':sha(Path(__file__)),
             'input_sha256':{p.relative_to(SITE).as_posix():h for p,h in PINS.items()},'rendered_math_count':count,
             'file_sha256':{f:sha(HERE/f) for f in ('INPUTS.md','INPUTS.html',FIGURE.name)},'local_links':links,
             'recomputed_exact_extrema':{k:str(v) for k,v in n.items()},
             'all32_exact_plot_data':[{k:str(v) if isinstance(v,Q) else v for k,v in r.items()} for r in rows],
             'new_mechanics_jobs':0,'root_or_minimum_claimed':False}
    path=HERE/'build_receipt.json'
    write(path,json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'receipt_sha256':sha(path),'files':receipt['file_sha256'],'math':count}))


if __name__=='__main__':
    main()
