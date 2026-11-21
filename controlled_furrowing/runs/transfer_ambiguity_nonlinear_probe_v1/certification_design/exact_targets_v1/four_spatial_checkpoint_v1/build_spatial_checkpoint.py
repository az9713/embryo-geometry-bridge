"""Build the four fixed-input spatial certificates report from independent receipts."""

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


HERE = Path(__file__).resolve().parent
EXACT = HERE.parent
DESIGN = EXACT.parent
RUN = DESIGN.parent
ROOT = RUN.parents[1]
SITE = ROOT.parent
ASSETS = SITE/'docs/assets/katex'
REDUCTION = EXACT/'bridge_reduction_v1/independent_actual_preparation'
COVERAGE = DESIGN/'full_outcome_coverage_v1'
AUTH = REDUCTION/'all_four_readback.json'
CONVERT = ROOT/'scripts/build_report.py'
HELPER = EXACT/'certification_checkpoint_v1/build_checkpoint.py'
FAVICON = RUN/'report/favicon.png'
PAST = EXACT/'cached_minimum_application_v1/MINIMA.html'
REPORT = RUN/'report/REPORT.html'
SWEEP = RUN/'observer_robustness/independent_actual/receipt.json'
ORIGINALITY = RUN/'originality_followup/receipt.json'
PLAN = ROOT/'PLAN.md'
FIGURE = HERE/'selfmap_budget.png'
TITLE = 'Four exact-spatial finite-element minima'
PILOTS = tuple(f'{endpoint}_{field}' for endpoint in ('max_contrast_change','max_negative_contrast_change')
               for field in ('small_17','large_03'))
LABELS = ('+ · small L','+ · large','− · small L','− · large')
RADIUS = Q(1, 2**34)
PINS = {
    AUTH:'859f507ad3762cfef44d53672bbe15a43592ce54d6e08ba8828bf50a925ec340',
    REDUCTION/'METHOD.md':'258048dcff70265a95f83eeb949d9cc2e4571fb2927017717efe5637a18d494e',
    EXACT/'preconditioner_spd_v1/independent_actual_preparation/actual_replay/receipt.json':
        '769a65a128aaa79018c050c7bd61d38999166621742946cd0f6a13876854ba83',
    EXACT/'independent_actual/replay/receipt.json':
        '884e0bc76a27465c86b62a202cab491d764aedf030623592238e09adee86d77d',
    COVERAGE/'COVERAGE.md':'e580b217050841f32abcc162cb989ab5975e7b8a8da0c6b5724b34b7405fa3b4',
    COVERAGE/'receipt.json':'235e8336c735b1bf1a75889018ae6a18586fd0905a65d38cb683c5ed236df19a',
    PAST:'6415becb909a5889eaa1f9175581dbb88a0e0b523a140cfb64fb000d8283b8cf',
    REPORT:'d33f3d594406ed498d364c3462450b62cb27054f0146f4bbca3a1c2ae22b41f6',
    SWEEP:'703b574bc51a083b10c7328d9622862e939ddabe791b18ce9e1817f07718f69f',
    ORIGINALITY:'6797ad5d05f71535f13edc519276e57bc5176d7d5b0e6273b1aa91eba48653d0',
    PLAN:'993240af5d30c8e5b72ab31a1278a524314d38ede48486b677b021d84bcde500',
    CONVERT:'90148d53d74d8ea5cc8a665e6a3660abf2db72487a1fa59586caba98f9cc38e2',
    HELPER:'50524d01db5003b0839e873efb37d41bb0552ad0f0c162f73bd668ad0e58d31b',
    ASSETS/'katex.min.css':'5a0ff74f11887f7260dc16b216906f51352b84ca29f21abb9a14b62a2ca6eb15',
    ASSETS/'katex.min.js':'e6bfe5deebd4c7ccd272055bab63bd3ab2c73b907b6e6a22d352740a81381fd4',
    FAVICON:'961e6e85ec16258652bf3c3f322d336fb598b90bccc49821047babfef3925530',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    with path.open('rb') as handle:
        hasher = sha256()
        for block in iter(lambda: handle.read(2**20), b''):
            hasher.update(block)
        return hasher.hexdigest()


def link(path):
    return Path(os.path.relpath(path,HERE)).as_posix()


def source_path(name):
    path=(ROOT/name).resolve()
    require(path.is_relative_to(ROOT) and path.relative_to(ROOT).as_posix()==name,
            'noncanonical study-relative evidence')
    return path


def load(path):
    return json.loads(path.read_bytes())


def reviewed_numbers():
    for path,h in PINS.items():
        require(sha(path)==h,f'changed report input: {path.name}')
    receipt=load(AUTH)
    hashes=receipt['current_file_sha256']
    require(receipt['status']=='PASS_FOUR_INDEPENDENT_FIXED_INPUT_EXACT_SPATIAL_ROOTS_AND_MINIMA'
            and receipt['independent_exact_spatial_root_count']==4
            and receipt['stronger_loaded_or_continuum_claim'] is False
            and receipt['full_five_gate_goal_achieved'] is False
            and tuple(c['pilot_id'] for c in receipt['cases'])==PILOTS
            and len(hashes)==receipt['current_file_binding_count']==3006,
            'four-case exact-spatial independent authority required')
    for name,h in hashes.items():
        require(sha(source_path(name))==h,f'changed independent evidence: {name}')
    rows=[]
    for case in receipt['cases']:
        n,o=(case[f'independent_{role}_file'] for role in ('proof','os_observation'))
        require(hashes[n]==case['independent_proof_sha256'] and hashes[o]==case['independent_os_observation_sha256'],
                'native proof or observation pin mismatch')
        native,observed=map(load,map(source_path,(n,o)))
        terminal_path=source_path(n).with_name('terminal.json')
        terminal=load(terminal_path)
        norm=Q(native['exact_R_infinity_norm'])
        r=Q(case['radius_exact'])
        comparison=case['selfmap']
        eta,z,lhs,margin,error=(Q(comparison[k]) for k in ('eta_upper_exact','z_upper_exact','selfmap_lhs_upper_exact',
                  'selfmap_margin_lower_exact','posterior_error_upper_exact'))
        m=Q(case['whole_box_normalized_H_min_exact'])
        goal=Q(case['state_observer_design_goal_exact'])
        target=case['target_identity']
        require(native['status']=='PASS_INDEPENDENT_EXACT_SPATIAL_ROOT_AND_STRICT_LOCAL_MINIMUM'
                and observed['status']=='PASS_OBSERVED_INDEPENDENT_GLOBAL_REDUCTION_TRUE_OS_EXIT'
                and all(type(w) is int and w==0 for w in (observed['observed_exit_code'],terminal['exit_code']))
                and type(native['pid']) is int and native['pid']>0 and native['pid']==observed['pid']==terminal['pid']
                and observed['receipt_sha256']==terminal['receipt_sha256']==hashes[n]
                and observed['terminal_sha256']==hashes[terminal_path.relative_to(ROOT).as_posix()]
                and native['pilot_id']==observed['pilot_id']==case['pilot_id']
                and native['target_identity']==target and Q(target['dose_exact'])*Q(target['gain_exact'])==Q(target['drive_exact'])
                and Q(target['energy_scale_binary_exact'])>0
                and native['counts']=={'bulk':144,'surface':24,'free':611}
                and case['bits']==native['precision_bits']==1024 and Q(native['radius_exact'])==r==RADIUS
                and native['conservative_compared_bounds']==comparison and comparison['strict_selfmap'] is True
                and eta>0 and 0<z<1 and lhs==eta+z*r<r and margin==r-lhs>0
                and norm>0 and m==Q(native['whole_box_normalized_H_euclidean_lower_exact'])==(1-z)/norm>0
                and error==eta/(1-z)<r<goal==Q(native['observer_state_design_goal_exact'])
                and native['posterior_below_observer_state_design_goal'] is True,
                f"failed exact linked theorem premises: {case['pilot_id']}")
        rows.append({'pilot_id':case['pilot_id'],'eta':eta,'z':z,'lhs':lhs,'margin':margin,'root_error':error,'m':m,
                     'state_goal':goal,'norm_R':norm,'dose':Q(target['dose_exact']), 'proof':source_path(n),'OS':source_path(o)})
    coverage=load(COVERAGE/'receipt.json')
    require((coverage['outcome_state_count'],coverage['total_primary_small_constraints'],
             coverage['primary_small_constraints_if_four_exact_pilots_succeed'],
             coverage['still_uncovered_after_four_successes'])==(36,136,8,128),
            'cohort denominator and selected-pilot accounting changed')
    return rows,Q(coverage['conditional_allowance_um_exact']),len(hashes)


def decimal(value,places,up):
    den=10**places
    scaled=value*den
    scaled_int=scaled.__ceil__() if up else scaled.__floor__()
    a,b=divmod(scaled_int,den)
    return f'{a}.{b:0{places}d}'


def figure(rows):
    bg,fg,muted,grid='#0c1826','#e9f3fb','#b2cadc','#31465c'
    point,derivative='#4bb2cb','#dba17d'
    fig,ax=plt.subplots(figsize=(8.2,4.25),dpi=180,facecolor=bg)
    ax.set_facecolor(bg)
    ax.set_axisbelow(True)
    ax.grid(axis='x',color=grid,linewidth=.7)
    ys=list(range(4))
    a=[float(row['eta']/RADIUS) for row in rows]
    b=[float(row['z']) for row in rows]
    ax.barh(ys,a,color=point,height=.55,label='Center residual  η/r')
    ax.barh(ys,b,left=a,color=derivative,height=.55,label='Whole-cube contraction  z')
    ax.axvline(1,color=fg,linewidth=1.1,linestyle='--')
    for i,row in enumerate(rows):
        ax.text(a[i]+b[i]+.018,i,decimal(row['lhs']/RADIUS,3,True),color=fg,va='center',fontsize=10)
    ax.set_xlim(0,1.04)
    ax.set_yticks(ys,LABELS,color=fg,fontsize=11)
    ax.invert_yaxis()
    ax.tick_params(axis='x',colors=muted)
    ax.tick_params(axis='y',length=0)
    ax.set_xlabel('Fraction of the proof radius',color=muted,labelpad=9)
    for name in ('left','right','top'):
        ax.spines[name].set_visible(False)
    ax.spines['bottom'].set_color(grid)
    ax.legend(ncols=2,loc='upper center',bbox_to_anchor=(.48,1.21),frameon=False,labelcolor=fg,fontsize=9)
    fig.text(.06,.91,'Every self-map bound stays below 1',color=fg,fontsize=18,weight='semibold')
    fig.text(.06,.844,'Four fixed inputs · exact-spatial finite-element equations',color=muted,fontsize=10)
    fig.subplots_adjust(left=.19,right=.94,bottom=.19,top=.68)
    with FIGURE.open('xb') as handle:
        fig.savefig(handle,format='png',dpi=180,facecolor=bg)
    plt.close(fig)


TEXT = r"""# Four exact-spatial finite-element minima

For each of two frozen material candidates, an independent calculation certifies one stationary root and a strict local minimum of the exactly integrated finite-element energy at small record L and at the matched large input. Each result applies to one closed cube in 611 free displacement coordinates at one fixed input. All four roots lie strictly inside their cubes.

## Numerical certificates

The common radius is \(r=2^{-34}\approx5.821\times10^{-11}\) in normalized displacement coordinates. Candidate + maximizes the frozen affine contrast direction and candidate − maximizes its negative. The final material points for this test were both frozen at half of the chosen affine steps. The table rounds upper bounds upward and the curvature lower bound downward.

![Four stacked horizontal bars give bounds on the center residual divided by the radius plus the contraction bound. Totals from top to bottom are at most 0.732, 0.381, 0.677 and 0.341. All totals fall below the dashed threshold at one.](selfmap_budget.png)

| Case | \(z\) | Self-map ratio | Error ×10⁻¹¹ | Curvature ×10⁻⁴ |
|:--|--:|--:|--:|--:|
[[TABLE]]

The self-map ratio bounds \((\eta+zr)/r\). Error bounds the free-coordinate infinity distance from the frozen center to the exact-spatial root. Curvature bounds the smallest eigenvalue of the **normalized** Hessian on the entire cube in frozen free coordinates. It has no direct interpretation as a measured modulus. Exact values and each observed native proof are linked below.

## Why the certificates hold

Let \(E_I\) be the normalized scalar energy evaluated by exact spatial integration on the fixed quadratic axisymmetric mesh. Its gradient and Hessian are \(F_I\) and \(H_I\). The matching cached operator has derivatives \(F_Q,H_Q\) and serves as an auxiliary calculation. The same case-specific geometry, linear free map, singleton material parameters, exact gain times dose and energy normalization bind both calculations. Independently reviewed body and trace domains provide admissibility throughout each cube. Strict positive domain margins and smooth integrands on the compact reference domains give a \(C^2\) energy near the cube.

The proofs separately enclose the signed point defect \(F_I(u_0)-F_Q(u_0)\) and the symmetric whole-cube derivative defect \(H_I-H_Q\). The independent reduction checks the 144 bulk and 24 surface component maps and uses separate exact endpoint scatter and directed interval products. It adds the already-preconditioned cached center force once, **before** taking norms. With the frozen binary matrix \(R\), the resulting bounds are

\[
\begin{gathered}
\eta\ge\|RF_I(u_0)\|_\infty,\\
z\ge\sup_{u\text{ in cube}}\|I-RH_I(u)\|_\infty.
\end{gathered}
\]

The derivative bound combines the cached defect bound with the preconditioned Hessian defect bound. Since \(z<1\), \(RH_I\) and the square matrix \(R\) are invertible. The strict inequality \(\eta+zr<r\) makes \(u\mapsto u-RF_I(u)\) a contraction into the cube interior. Banach's fixed-point theorem gives a unique stationary point on the cube and

\[
\|u_*-u_0\|_\infty\le e=\frac{\eta}{1-z}<r .
\]

For the minimum, the separately certified symmetric part \(P=(R+R^{\mathsf T})/2\) is positive definite. All eigenvalues of \(A=RH_I\) have positive real part. With \(S=H_I^{-1}\), symmetry gives \(AS+SA^{\mathsf T}=2P\). Integrating the derivative of \(e^{-tA}Se^{-tA^{\mathsf T}}\) from zero to infinity gives a positive definite integral representation of \(S\). Thus \(H_I\) is positive definite. The Neumann bound on the inverse then gives

\[
\lambda_{\min}(H_I)\ge m=\frac{1-z}{\|R\|_\infty}>0 .
\]

Taylor integration along each segment from the interior stationary point gives an energy gap of at least \(m\|u-u_*\|_2^2/2\). The point is a strict local minimum and the unique minimum on its declared cube.

## Coverage and scientific meaning

These four proofs cover **4 of 36** final outcome fields. They cover the two large contrast inputs and the L field for each candidate. The two L fields span **8 of 136** primary small comparisons. The other 32 fields and their **128** small comparisons still require exact-spatial root certificates. The root errors pass the prior state-error design budgets for these four fields. A reported observation bound also requires transfer through the chosen observer, including stored-coordinate rounding and moving evaluation points.

The separate stored-field report finds a large contrast gap of about 0.154755 µm under the predeclared scenario. The frozen [[DELTA]] µm allowance is a conditional model residual budget. Both saved candidates fit it under the primary disk/spatial-maximum convention. The same candidates exceed it under the other three audited static conventions. Experimental compatibility, dynamics along a stable loading path, the physical identified set, continuum and nonaxisymmetric accuracy, and literature novelty remain unresolved. The original five-gate objective remains active.

## Evidence

- [Consolidated independent four-case certificate]([[AUTH]]) · [Exact spatial reduction method]([[METHOD]])
- [[PROOFS]]
- [[OS]]
- [Whole-cube geometry replay]([[GEOMETRY]]) · [Positive symmetric parts of the four preconditioners]([[SPD]])
- [Fixed-field coverage ledger]([[COVERAGE]]) · [Conditional nonlinear report]([[REPORT]]) · [Observer sensitivity replay]([[SWEEP]])
- [Earlier cached-energy minima]([[PAST]]) · [Originality audit]([[ORIGINALITY]]) · [Five-gate research plan]([[PLAN]])
"""


def research_text(rows,delta):
    table=[]
    for label,row in zip(LABELS,rows,strict=True):
        values=(row['z'],row['lhs']/RADIUS,row['root_error']*10**11,row['m']*10**4)
        rounded=[decimal(v,4,i!=3) for i,v in enumerate(values)]
        table.append('| '+' | '.join([label,*rounded])+' |')
    paths={'AUTH':AUTH,'METHOD':REDUCTION/'METHOD.md','GEOMETRY':EXACT/'independent_actual/replay/receipt.json',
           'SPD':EXACT/'preconditioner_spd_v1/independent_actual_preparation/actual_replay/receipt.json',
           'COVERAGE':COVERAGE/'COVERAGE.md','REPORT':REPORT,'SWEEP':SWEEP,'PAST':PAST,
           'ORIGINALITY':ORIGINALITY,'PLAN':PLAN}
    values={**{key:link(path) for key,path in paths.items()},'TABLE':'\n'.join(table),
            'DELTA':f'{float(delta):.8f}'}
    for token,key in (('PROOFS','proof'),('OS','OS')):
        values[token]=' · '.join(f'[{label} {key}]({link(row[key])})' for label,row in zip(LABELS,rows,strict=True))
    text=TEXT
    for key,value in values.items():
        text=text.replace(f'[[{key}]]',value)
    require('[[' not in text,'unfilled document field')
    return text


CSS='''
*{box-sizing:border-box}html{color-scheme:dark}
body{margin:0;background:#091422;color:#e5f0fb;font:1rem/1.65 system-ui,sans-serif}
main{width:min(100%,55rem);min-width:0;margin:auto;padding:2rem clamp(.95rem,3vw,2.4rem) 4rem}
nav{display:flex;flex-wrap:wrap;gap:.25rem 1.15rem;margin:0 0 1.6rem;font-size:.91rem}
h1,h2{line-height:1.2;letter-spacing:-.025em;text-wrap:balance}
h1{font-size:clamp(2rem,5.5vw,3.15rem);margin:.6rem 0 1rem}
h2{margin:2rem 0 .8rem;font-size:clamp(1.38rem,4vw,1.84rem)}
h1+p{font-size:1.12rem;color:#c7dbea}
a{color:#8ae5fa;text-underline-offset:.18em}a:focus-visible{outline:2px solid #ffd0a8;outline-offset:3px}
img{display:block;width:100%;height:auto;border-radius:8px}
table{width:100%;border-collapse:collapse;font-size:.91rem;text-align:right;font-variant-numeric:tabular-nums}
td,th{padding:.54rem .56rem;border-bottom:1px solid #36536d}
th{color:#c2dcec}td:first-child,th:first-child{text-align:left}
.math.display{margin:1.1rem 0;max-width:100%;overflow-x:auto}.katex-display{margin:0}.math{color:#f1f7ff}
li{margin:.42rem 0}
@media(max-width:500px){main{padding-top:1.05rem}td,th{padding:.4rem .24rem;font-size:.75rem;line-height:1.35}
.math.display .katex{font-size:1em}ul{padding-left:1.2rem}}
'''


def render(text):
    spec=spec_from_file_location('reviewed_static_math_renderer',HELPER)
    helper=module_from_spec(spec)
    spec.loader.exec_module(helper)
    sys.path.insert(0,str(ROOT/'scripts'))
    import build_report
    require(Path(build_report.__file__).resolve()==CONVERT,'wrong frozen Markdown converter')
    body,count=helper.compile_math(build_report.convert(text))
    page=('<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
          '<meta name="viewport" content="width=device-width,initial-scale=1">'
          f'<title>{escape(TITLE)}</title><link rel="icon" type="image/png" href="{link(FAVICON)}">'
          f'<link rel="stylesheet" href="{link(ASSETS/"katex.min.css")}"><style>{CSS}</style></head>'
          '<body><main><nav aria-label="Related reports">'
          f'<a href="{link(REPORT)}">Nonlinear study</a><a href="{link(PAST)}">Cached minima</a>'
          '<a href="SPATIAL.md">Markdown</a></nav>'+body+'</main></body></html>\n')
    html_links=helper.LocalLinks()
    html_links.feed(page)
    targets={}
    for url in html_links.paths:
        path=(HERE/url).resolve()
        require(path.is_relative_to(SITE) and (path.is_file() or path==HERE/'SPATIAL.md'),f'missing local resource: {url}')
        targets[url]=path.relative_to(SITE).as_posix()
    return page,count,targets


def write(path,contents):
    with path.open('x',encoding='utf-8',newline='\n') as handle:
        handle.write(contents)


def main():
    rows,delta,checked=reviewed_numbers()
    figure(rows)
    md=research_text(rows,delta)
    page,count,targets=render(md)
    write(HERE/'SPATIAL.md',md)
    write(HERE/'SPATIAL.html',page)
    receipt={'status':'PASS_BUILD_FROM_FOUR_INDEPENDENT_EXACT_SPATIAL_CERTIFICATES',
             'source_sha256':sha(Path(__file__)),'rendered_math_count':count,'local_links':targets,
             'file_sha256':{name:sha(HERE/name) for name in ('SPATIAL.md','SPATIAL.html',FIGURE.name)},
             'input_sha256':{p.relative_to(SITE).as_posix():h for p,h in PINS.items()},
             'independent_current_file_bindings_rehashed':checked,
             'four_exact_rows':[{k:str(v) if isinstance(v,Q) else link(v) if isinstance(v,Path) else v
                                 for k,v in row.items()} for row in rows],
             'radius_exact':str(RADIUS),'conditional_allowance_um_exact':str(delta),
             'exact_spatial_roots':4,'exact_spatial_strict_local_minima':4,'outcome_fields':36,
             'remaining_exact_spatial_root_certificates':32,'remaining_small_comparisons':128,
             'full_five_gate_goal_achieved':False,'new_mechanics_jobs':0}
    path=HERE/'build_receipt.json'
    write(path,json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'receipt_sha256':sha(path),'rendered_math':count,'files':receipt['file_sha256']}),flush=True)


if __name__=='__main__':
    main()
