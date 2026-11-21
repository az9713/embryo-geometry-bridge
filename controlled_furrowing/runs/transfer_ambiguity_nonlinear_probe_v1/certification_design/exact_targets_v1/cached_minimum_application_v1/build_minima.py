"""Build the reviewed cached-minima note, exact-data chart and static KaTeX HTML."""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
import os
from pathlib import Path
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


HERE = Path(__file__).resolve().parent
EXACT = HERE.parent
DESIGN = EXACT.parent
RUN = DESIGN.parent
ROOT = EXACT.parents[3]
SITE = ROOT.parent
ASSET = SITE/'docs/assets/katex'
FIGURE = HERE/'curvature_bounds.png'
CHECKPOINT = EXACT/'certification_checkpoint_v1'
S = EXACT/'preconditioner_spd_v1/independent_actual_preparation/actual_replay/receipt.json'
D = EXACT/'cached_derivative_extension_v1/independent_actual_preparation/actual_replay_v2/receipt.json'
C = DESIGN/'cached_pilots_v1/independent_actual_preparation/actual_replay/receipt.json'
G = EXACT/'independent_actual/replay/receipt.json'
FAVICON = RUN/'report/favicon.png'
CONVERT = ROOT/'scripts/build_report.py'
HELPER = CHECKPOINT/'build_checkpoint.py'
PINS = {
    HERE/'independent_readback.json':'a7bfee53506f54837c6a0c04f2c67b4ae0b15deca76e133c6a3f16101a25fe10',
    HERE/'root_independent_readback.py':'b5661e26fc7ae69c5f005dbf98191a730e085de26ab66b9c9c4cdd4d45727296',
    HERE/'results.json':'005f1f100da3cd558212d6eceef7a1ca886602e94c64438f08dd009f69e48b61',
    HERE/'PROOF.md':'edf6c75d294a70052cff22302785766e14cea31a036f22958e06169c3a5f4b4f',
    S:'769a65a128aaa79018c050c7bd61d38999166621742946cd0f6a13876854ba83',
    D:'f766b8b41142ce1f82f3e672f423c32485ae436487d9418c42d907b6f27e8c74',
    C:'4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17',
    G:'884e0bc76a27465c86b62a202cab491d764aedf030623592238e09adee86d77d',
    CONVERT:'90148d53d74d8ea5cc8a665e6a3660abf2db72487a1fa59586caba98f9cc38e2',
    HELPER:'50524d01db5003b0839e873efb37d41bb0552ad0f0c162f73bd668ad0e58d31b',
    ASSET/'katex.min.css':'5a0ff74f11887f7260dc16b216906f51352b84ca29f21abb9a14b62a2ca6eb15',
    ASSET/'katex.min.js':'e6bfe5deebd4c7ccd272055bab63bd3ab2c73b907b6e6a22d352740a81381fd4',
    FAVICON:'961e6e85ec16258652bf3c3f322d336fb598b90bccc49821047babfef3925530',
    CHECKPOINT/'CHECKPOINT.html':'38557df31bc3476c3eb558439e252074dc3267e5813b3b9081a863f79569827a',
}
TITLE='Four strict minima of the cached energy'
PILOTS=tuple(f'{a}_{b}' for a in ('max_contrast_change','max_negative_contrast_change')
             for b in ('small_17','large_03'))
LABELS=('Candidate + · small L','Candidate + · large','Candidate − · small L','Candidate − · large')
BG,FG,MUTED,GRID='#0b1320','#eaf3fc','#a8bbd0','#314660'
COLORS=('#81d8ef','#81d8ef','#e2ac86','#e2ac86')


def require(ok,msg):
    if not ok:
        raise ValueError(msg)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def link(path):
    return Path(os.path.relpath(path,HERE)).as_posix()


def pin(path):
    raw=path.read_bytes()
    require(sha256(raw).hexdigest()==PINS[path],f'changed input: {path}')
    return json.loads(raw)


def numbers():
    for path,h in PINS.items():
        require(sha(path)==h,f'changed source, evidence or asset: {path}')
    review,draft=map(pin,(HERE/'independent_readback.json',HERE/'results.json'))
    require(review['status']=='PASS_INDEPENDENT_FOUR_CACHED_MINIMUM_APPLICATION'
            and review['cached_local_minima_certified']==4 and review['exact_spatial_minima_certified']==0
            and review['stable_loading_or_continuum_or_biological_fit_claimed'] is False
            and review['full_five_gate_goal_complete'] is False and all(v is True for v in review['checks'].values())
            and review['review_source_sha256']==PINS[HERE/'root_independent_readback.py']
            and review['applicant_results_sha256']==PINS[HERE/'results.json']
            and tuple(r['pilot_id'] for r in review['rows'])==PILOTS,
            'complete four-case independent promotion receipt required')
    for name,h in draft['source_and_input_sha256'].items():
        path=(ROOT/name).resolve()
        require(path.is_relative_to(ROOT) and path.relative_to(ROOT).as_posix()==name and sha(path)==h,
                'a proof premise changed')
    rows=[]
    for row,old in zip(review['rows'],draft['rows'],strict=True):
        m,gap,margin=(Q(row[k]) for k in ('normalized_cached_H_eigenlower_exact','energy_gap_coefficient_exact',
                                          'root_large_cube_interior_margin_exact'))
        require(row['pilot_id']==old['target_identity']['pilot_id'] and m==Q(old['normalized_cached_H_eigenlower_exact'])
                and gap==Q(old['normalized_cached_energy_quadratic_gap_coefficient_exact'])==m/2
                and 0<margin==Q(old['root_interior_margin_coercive_cube_exact'])<Q(1,2**34),
                'independent and draft exact arithmetic differ')
        rows.append((m,gap,margin))
    return rows


def downward(x,digits=4):
    """Floor a positive exact rational for a conservative displayed lower bound."""
    den=10**digits
    n=(x*den).__floor__()
    return f'{n//den}.{n%den:0{digits}d}'


def figure(rows):
    values=[10**4*m for m,_,_ in rows]
    fig,ax=plt.subplots(figsize=(9.6,4.3),dpi=170,facecolor=BG)
    ax.set_facecolor(BG)
    ax.set_axisbelow(True)
    ax.grid(axis='x',color=GRID,linewidth=.7)
    ax.barh(range(4),list(map(float,values)),color=COLORS,height=.55)
    ax.set_yticks(range(4),LABELS,color=FG,fontsize=11)
    ax.invert_yaxis()
    ax.set_xlim(0,7.9)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.tick_params(axis='x',colors=MUTED)
    ax.tick_params(axis='y',length=0,pad=9)
    ax.set_xlabel('Lower bound on normalized Hessian eigenvalue (×10⁻⁴)',color=FG,labelpad=12)
    for side in ('top','right','left'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color(GRID)
    for i,value in enumerate(values):
        ax.text(float(value)+.1,i,downward(value),va='center',ha='left',color=FG,fontsize=12)
    fig.text(.06,.92,'Positive curvature in all four cached cubes',color=FG,fontsize=19,weight='semibold')
    fig.text(.06,.855,'Independent exact bounds  ·  closed free-coordinate radius 2⁻³⁴',color=MUTED,fontsize=11)
    fig.subplots_adjust(left=.285,right=.965,bottom=.18,top=.78)
    with FIGURE.open('xb') as f:
        fig.savefig(f,format='png',dpi=170,facecolor=BG)
    plt.close(fig)


TEMPLATE=r"""# Four strict minima of the cached energy

The two frozen material candidates each have a verified strict local minimum at small record L and at the matched large input. These four minima belong to the normalized binary \(Q8/S12\) quadrature energy. Exact-spatial energy minima remain a separate proof obligation.

## Four positive curvature bounds

The earlier root certificate locates the stationary point \(u_*\) inside a closed \(2^{-36}\) free-coordinate cube. The new minimum proof gives uniform positive curvature on the concentric \(2^{-34}\) cube with the same center, linear free coordinates, exact gain times dose, and frozen binary preconditioner \(R\).

![Four positive lower bounds on the smallest normalized cached Hessian eigenvalue, ranging from 4.0988 to 6.5349 times 10 to the minus four.](curvature_bounds.png)

| Frozen target | Curvature \(m_Q\) (×10⁻⁴) | Root interior margin (×10⁻¹¹) |
|:--|--:|--:|
[[TABLE]]

Both columns are lower bounds rounded downward. The interior margins refer to the larger cube in dimensionless free coordinates. The curvature is for the normalized cached energy in those coordinates.

## Why the sign follows

On each admissible cube, positive cached-site analytic margins make the scalar energy \(E_Q\) twice differentiable with symmetric Hessian \(H_Q\). The fresh symmetric-part proof gives \(P=\operatorname{sym}R\succ0\). The larger-cube cached derivative bound gives \(\|I-RH_Q\|_\infty \le \bar Z_Q<1\). With \(A=RH_Q\) and \(S=H_Q^{-1}\), \(A\) is invertible and has eigenvalues with positive real part. The Lyapunov identity yields

\[
\begin{aligned}
AS+SA^\mathsf{T}&=2P,\\
S&=\int_0^\infty e^{-tA}\,2P\,e^{-tA^\mathsf{T}}\,dt\succ0.
\end{aligned}
\]

The Neumann bound and the exact binary norm \(N_R=\|R\|_\infty\) give the uniform quantitative result

\[
\lambda_{\min}(H_Q)\ge m_Q=\frac{1-\bar Z_Q}{N_R}>0 .
\]

Invertibility and contraction alone give no sign. For example, \(R=H_Q=-I\) gives \(RH_Q=I\), yet the quadratic stationary point is a strict maximum.

The root lies strictly inside the larger cube. Every segment from \(u_*\) to \(u\) in that cube stays admissible. Taylor integration with vanishing gradient at \(u_*\) gives

\[
E_Q(u)-E_Q(u_*)\ge \tfrac12 m_Q\|u-u_*\|_2^2 .
\]

This proves one strict minimum per specified cached-energy cube. The four certificates concern isolated fixed targets. Continuous loading, exact spatial integration, continuum accuracy and a measured biological fit require further evidence.

## Proof records

The independently verified receipt resolves the review status of the archived application draft. The exact rational bounds and the same-target proof links are retained for replay. [Independent minimum receipt](independent_readback.json) · [Detailed derivation](PROOF.md) · [Original scalar draft](results.json)

[Cached roots]([[C]]) · [Whole-cube geometry]([[G]]) · [Larger cached derivative bounds]([[D]]) · [Actual symmetric-part proofs]([[S]]) · [Four-pilot checkpoint]([[CHECKPOINT]])
"""


def write(path,text):
    with path.open('x',encoding='utf-8',newline='\n') as f:
        f.write(text)


def main():
    rows=numbers()
    table='\n'.join(f'| {label} | {downward(10**4*m)} | {downward(10**11*h)} |'
                    for label,(m,_,h) in zip(LABELS,rows,strict=True))
    tokens={'TABLE':table,**{k:link(v) for k,v in {'S':S,'D':D,'G':G,'C':C,'CHECKPOINT':CHECKPOINT/'CHECKPOINT.html'}.items()}}
    md=re.sub(r'\[\[([A-Z_]+)\]\]',lambda m:tokens[m[1]],TEMPLATE)
    spec=importlib.util.spec_from_file_location('_pinned_checkpoint_math_compiler',HELPER)
    helper=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    sys.path.insert(0,str(ROOT/'scripts'))
    import build_report
    require(Path(build_report.__file__).resolve()==CONVERT,'converter source differs')
    body,count=helper.compile_math(build_report.convert(md))
    page=('<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
          '<meta name="viewport" content="width=device-width,initial-scale=1">'
          f'<title>{TITLE}</title><link rel="icon" type="image/png" href="{link(FAVICON)}">'
          f'<link rel="stylesheet" href="{link(ASSET/"katex.min.css")}"><style>{helper.CSS}</style>'
          '</head><body><main><nav aria-label="Related documents">'
          f'<a href="{link(CHECKPOINT/"CHECKPOINT.html")}">Certification checkpoint</a>'
          f'<a href="{link(RUN/"report/REPORT.html")}">Nonlinear report</a>'
          '<a href="MINIMA.md">Markdown</a></nav>'+body+'</main></body></html>\n')
    paths=helper.LocalLinks()
    paths.feed(page)
    for url in paths.paths:
        path=(HERE/url).resolve()
        require(path.is_relative_to(SITE)
                and (path.is_file() or path in {HERE/'MINIMA.md',FIGURE}),'broken local minima link')
    figure(rows)
    write(HERE/'MINIMA.md',md)
    write(HERE/'MINIMA.html',page)
    receipt={'status':'BUILT_FOUR_VERIFIED_CACHED_MINIMA_NOTE','source_sha256':sha(Path(__file__)),
             'input_sha256':{p.relative_to(SITE).as_posix():h for p,h in PINS.items()},
             'file_sha256':{n:sha(HERE/n) for n in ('MINIMA.md','MINIMA.html',FIGURE.name)},
             'rendered_math_count':count,'proof_premise_hashes_read_back':57,
             'exact_figure_and_table_curvature_lower_bound':[str(r[0]) for r in rows],
             'exact_energy_gap_coefficient':[str(r[1]) for r in rows],
             'exact_large_cube_root_interior_margin':[str(r[2]) for r in rows],
             'certified_cached_minima':4,'exact_spatial_minima_from_this_certificate':0,
             'new_matrix_force_Hessian_AD_or_solves':0}
    path=HERE/'minima_build_receipt.json'
    write(path,json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'receipt_sha256':sha(path),'compiled_math':count,'files':receipt['file_sha256']}))


if __name__=='__main__':
    main()
