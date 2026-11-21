"""Build the standalone checkpoint from frozen independent numeric receipts."""

from fractions import Fraction as Q
from hashlib import sha256
from html import escape, unescape
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import subprocess
import sys


HERE = Path(__file__).resolve().parent
EXACT = HERE.parent
DESIGN = EXACT.parent
RUN = DESIGN.parent
ROOT = HERE.parents[4]
SITE = ROOT.parent
ASSETS = SITE/'docs/assets/katex'
FAVICON = RUN/'report/favicon.png'
DERIV = EXACT/'cached_derivative_extension_v1/independent_actual_preparation'
CACHE = DESIGN/'cached_pilots_v1'
OBS = CACHE/'observer_application_v1/independent_readback'
LEDGER = DESIGN/'full_outcome_coverage_v1'
M_PLUS = LEDGER/'remaining_cached_wave1/independent_actual_preparation/actual_replay'
M_MINUS = LEDGER/'remaining_cached_wave1/negative_M_narrow_v1/actual'
PINS = {
    HERE/'receipt.json':'1227947005359ddd23c27c4edbde368419dfe835ad1f32992ff549f68e0da337',
    HERE/'figure_receipt.json':'88ce72ed63f18feac51046f5e8c0d6b5c59a5626f6c45a11c8545ba4ef775117',
    HERE/'cached_derivative_headroom.png':'a68b0f1c018580716e5b69082ad4eb9c79a89550b4025908b1673b29fae0cb11',
    DERIV/'bridge_admission.json':'269229312793ec101cc301086d634dd4de89c62a5050ccfc60f185cf761dbc7d',
    DERIV/'actual_replay_v2/receipt.json':'f766b8b41142ce1f82f3e672f423c32485ae436487d9418c42d907b6f27e8c74',
    CACHE/'independent_actual_preparation/actual_replay/receipt.json':
        '4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17',
    EXACT/'independent_actual/replay/receipt.json':'884e0bc76a27465c86b62a202cab491d764aedf030623592238e09adee86d77d',
    OBS/'receipt.json':'939212f0690f4ccdf27ede99a5b7e7be6eff4a080adfd2dda78c0352c2e9c2e5',
    OBS/'artifact_manifest.json':'a3171d87420caafdef25737898563ab39be112194e5b7afa6c6421150fdc4c4a',
    LEDGER/'receipt.json':'235e8336c735b1bf1a75889018ae6a18586fd0905a65d38cb683c5ed236df19a',
    M_PLUS/'receipt.json':'9854b5033a055ea33d7af818cf8da27b8eee35b3f4bb57606fcbf9bca96e3b56',
    M_PLUS/'max_contrast_change_small_09.json':'c5c78d779e4537859c6a0aaf4b206a1339ba6cc7d1392a99fdcb4d7255d007b4',
    M_PLUS/'os_exit_observation.json':'42e11a4b7065b07eca318e5127951e884e8bd7703773b37df184bc9b747706ef',
    M_MINUS/'receipt.json':'7aa4c765a6df459e802e2d1e99c7fb007ab533e3b29ef07cee2d92df3f52d4f1',
    M_MINUS/'os_exit_observation.json':'3fa3737870a585f58ae7ca846b9d696c2ee98df87820c44f503529214273a6b8',
    RUN/'observer_robustness/actual/results.json':'53a17296ddae08b4d1fe6b7e9b20ad0a377bf8dc962142bbfb9eb9877a381647',
    RUN/'observer_robustness/independent_actual/receipt.json':'703b574bc51a083b10c7328d9622862e939ddabe791b18ce9e1817f07718f69f',
    RUN/'report/REPORT.html':'d33f3d594406ed498d364c3462450b62cb27054f0146f4bbca3a1c2ae22b41f6',
    FAVICON:'961e6e85ec16258652bf3c3f322d336fb598b90bccc49821047babfef3925530',
    ROOT/'PLAN.md':'993240af5d30c8e5b72ab31a1278a524314d38ede48486b677b021d84bcde500',
    ROOT/'scripts/build_report.py':'90148d53d74d8ea5cc8a665e6a3660abf2db72487a1fa59586caba98f9cc38e2',
    ASSETS/'katex.min.css':'5a0ff74f11887f7260dc16b216906f51352b84ca29f21abb9a14b62a2ca6eb15',
    ASSETS/'katex.min.js':'e6bfe5deebd4c7ccd272055bab63bd3ab2c73b907b6e6a22d352740a81381fd4',
}
TITLE = "Four cached roots and the exact-spatial test"
VARIANTS = ('current_disk__spatial_max','current_disk__regional_mean',
            'current_radial_line__spatial_max','current_radial_line__regional_mean')
LABELS = ('Disk · separate component maxima','Disk · component means',
          'Radial line · separate component maxima','Radial line · component means')


def require(ok,message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def pinned(path):
    require(sha(path)==PINS[path],f'changed scientific or renderer input: {path.name}')
    return json.loads(path.read_bytes())


def link(path):
    return Path(os.path.relpath(path,HERE)).as_posix()


def site_name(path):
    return path.relative_to(SITE).as_posix()


def fixed_integer(n,places):
    a,b = divmod(n,10**places)
    return f'{a}.{b:0{places}d}'


def outward(lo,hi,places):
    factor=10**places
    return fixed_integer((lo*factor).__floor__(),places),fixed_integer((hi*factor).__ceil__(),places)


def reviewed_numbers():
    for path,h in PINS.items():
        require(sha(path)==h,f'changed pinned input: {site_name(path)}')
    checkpoint=pinned(HERE/'receipt.json')
    root_inventory=checkpoint['checked_file_sha256']
    require(checkpoint['status']=='PASS_BOUND_COMPONENT_READBACK_ONLY' and len(root_inventory)==236
            and all(v is True for v in checkpoint['checks'].values())
            and checkpoint['new_exact_spatial_integral_roots_certified']==0
            and checkpoint['full_five_gate_goal_complete'] is False,'wrong checkpoint authority')
    for name,h in root_inventory.items():
        path=(ROOT/name).resolve()
        require(path.is_relative_to(ROOT) and path.relative_to(ROOT).as_posix()==name
                and sha(path)==h,f'changed checkpoint root-readback bytes: {name}')
    cache=pinned(CACHE/'independent_actual_preparation/actual_replay/receipt.json')
    geo=pinned(EXACT/'independent_actual/replay/receipt.json')
    deriv=pinned(DERIV/'actual_replay_v2/receipt.json')
    wrap=pinned(DERIV/'bridge_admission.json')
    fig=pinned(HERE/'figure_receipt.json')
    require(cache['status']=='PASS_INDEPENDENT_FOUR_CACHED_LOCAL_ROOTS'
            and cache['root_certificates_cached_quadrature_only']==4
            and cache['root_certificates_exact_spatial_integral']==0
            and geo['status']=='PASS_INDEPENDENT_ACTUAL_EIGHT_SINGLETON_GEOMETRY_BOXES'
            and deriv['status']=='PASS_INDEPENDENT_FOUR_LARGER_CACHED_DERIVATIVES'
            and wrap['reviewed_cases']==deriv['reviewed_cases']
            and fig['input_sha256']==PINS[DERIV/'bridge_admission.json']
            and fig['figure_sha256']==PINS[HERE/'cached_derivative_headroom.png'],
            'cached proof, whole geometry and figure do not refer to the same reviewed inputs')
    radii=set()
    for row in cache['field_results']:
        path=CACHE/'independent_actual_preparation/actual_replay'/f"{row['pilot_id']}.json"
        require(root_inventory[path.relative_to(ROOT).as_posix()]==row['receipt_sha256']==sha(path),
                'changed proof of individual cached root')
        proof=json.loads(path.read_bytes())
        radii.add(Q(proof['radius_exact']))
        require(proof['status']=='PASS_INDEPENDENT_CACHED_LOCAL_ROOT','cached root missing')
    zq={r['pilot_id']:Q(r['reviewed_ZQ_upper_exact']) for r in deriv['reviewed_cases']}
    require(len(zq)==4 and set(zq)=={q['pilot_id'] for q in cache['field_results']}
            and zq=={k:Q(v) for k,v in checkpoint['reviewed_ZQ_upper_exact_by_pilot'].items()}
            and list(zq.values())==list(map(Q,fig['plotted_bounds_exact']))
            and set(Q(x['radius_exact']) for x in deriv['reviewed_cases'])=={Q(1,2**34)}
            and radii=={Q(1,2**36)} and all(0<z<1 for z in zq.values()),
            'cached derivative bounds or reviewed radii differ')
    obs=pinned(OBS/'receipt.json')
    obs_manifest=pinned(OBS/'artifact_manifest.json')
    ledger=pinned(LEDGER/'receipt.json')
    sensitivity=pinned(RUN/'observer_robustness/actual/results.json')
    independent_sweep=pinned(RUN/'observer_robustness/independent_actual/receipt.json')
    require(obs['status']=='PASS_INDEPENDENT_CACHED_ROOT_OBSERVER_APPLICATION_READBACK'
            and obs['comparison_counts']==checkpoint['cached_observer_comparison_counts']=={'IN':28,'OUT':4}
            and len(obs['all_32_comparisons'])==32 and obs['tracked_small_root_fields']==2
            and obs['exact_spatial_roots_certified']==0
            and sha(OBS/'DERIVATION.md')==obs_manifest['files_sha256']['DERIVATION.md']
            and independent_sweep['status']=='PASS_INDEPENDENT_SECONDARY_STORED_OBSERVER_SWEEP'
            and sensitivity['measured_small_fields']==34
            and ledger['status']=='PASS_READ_ONLY_36_FIELD_COVERAGE_INVENTORY',
            'cached observer arithmetic or separate sensitivity/coverage scope changed')
    low,high=gap=tuple(map(Q,obs['exact_large_gap_um']))
    delta=Q(obs['conditional_delta_um_exact'])
    require(0<low<=high and gap==tuple(map(Q,checkpoint['cached_observer_large_gap_um_exact']))
            and delta==Q(ledger['conditional_allowance_um_exact'])
            ==Q(sensitivity['frozen_primary_source_allowance_um_exact']),
            'inconsistent conditional gap or allowance')
    requirements=sensitivity['common_pair_requirement_by_variant']
    require([requirements[k]['status_at_frozen_allowance'] for k in VARIANTS]
            ==['ALLOWANCE_SATISFIED',*(['ALLOWANCE_EXCEEDED']*3)],'stored-field observer classification changed')
    require((ledger['outcome_state_count'],ledger['independently_cached_root_count'],
             ledger['independent_exact_spatial_root_count'],ledger['remaining_new_cache_and_singleton_cases'],
             ledger['total_primary_small_constraints'],ledger['primary_small_constraints_if_four_exact_pilots_succeed'],
             ledger['still_uncovered_after_four_successes'])==(36,4,0,32,136,8,128),
             'full-outcome certificate and 68 inequalities per endpoint accounting differ')
    return zq,gap,delta,requirements,ledger,obs,len(root_inventory)


def reviewed_m_followup():
    combined=pinned(M_PLUS/'receipt.json')
    plus=pinned(M_PLUS/'max_contrast_change_small_09.json')
    minus=pinned(M_MINUS/'receipt.json')
    for folder,pid in ((M_PLUS,combined['source_and_phase']['pid']),(M_MINUS,minus['pid'])):
        observed=pinned(folder/'os_exit_observation.json')
        path=folder/'terminal.json'
        require(sha(path)==observed['terminal_sha256'],'M replay terminal bytes changed')
        terminal=json.loads(path.read_bytes())
        require(type(observed['observed_exit_code']) is int and observed['observed_exit_code']==0
                and type(terminal['exit_code']) is int and terminal['exit_code']==0 and terminal['status']=='FINISHED'
                and observed['pid']==terminal['pid']==pid
                and observed['receipt_sha256']==terminal['receipt_sha256']==PINS[folder/'receipt.json'],
                'same observed successful M replay required')
    rplus,rminus=Q(plus['radius_exact']),Q(minus['radius_exact'])
    require(combined['status']=='INDEPENDENT_TWO_M_CACHED_REPLAY_WITH_UNRESOLVED_CASE'
            and combined['cached_local_roots']==1 and combined['spatial_root_count']==0
            and plus['status']=='PASS_INDEPENDENT_M_CACHED_LOCAL_ROOT'
            and plus['case_id']=='max_contrast_change_small_09'
            and combined['results'][0]['native_sha256']==PINS[M_PLUS/'max_contrast_change_small_09.json']
            and combined['results'][1]['status']=='UNRESOLVED_INDEPENDENT_M_CACHED_ROOT'
            and minus['status']=='PASS_INDEPENDENT_NEGATIVE_M_SMALLER_CACHED_ROOT'
            and minus['case_id']=='max_negative_contrast_change_small_09'
            and minus['field_sha256']==combined['results'][1]['field_sha256']
            and minus['cache_packet_sha256']==combined['results'][1]['cache_packet_sha256']
            and minus['original_unresolved_test_preserved'] is True
            and minus['spatial_integration_or_root_claim'] is False
            and (rplus,rminus,Q(minus['original_unresolved_radius_exact']))==(Q(1,2**36),Q(1,2**34),Q(1,2**32)),
            'new narrow cached M results differ from the reviewed follow-up')
    return {plus['case_id']:str(rplus),minus['case_id']:str(rminus)}


TEMPLATE = r"""# Four cached roots and the exact-spatial test

The two frozen material candidates retain a strictly positive large-patch contrast difference in the independently verified **cached-quadrature model**. The signed gap is **[[GAP_LO]] to [[GAP_HI]] µm**. The exact-spatial finite-element root test remains pending in this checkpoint.

## What the four proofs cover

Each candidate contributes one small state, record L, and one large state at the common matched input. The proofs certify four local roots of the **611-coordinate binary Q8/S12 quadrature potential**, with the same frozen discrete model and each state's exact requested optical drive. Q8/S12 refers to the quadrature orders for bulk and surface terms. The drive in the exact arithmetic is the candidate gain times its exact dose. The rounded drive that the numerical solver used is retained separately.

The four independent root proofs use closed cubes of radius \(2^{-36}\) in dimensionless free coordinates. A separate checked geometry certificate covers the entire larger cube of radius \(r=2^{-34}\) for each center. A second independent interval replay encloses the cached derivative over that whole larger cube with the same binary preconditioner \(R\). Define \(\overline{Z}_Q\) as the resulting upper bound, with \(H_Q\) the cached-quadrature Hessian.

\[
\begin{aligned}
\|I-RH_Q(u)\|_\infty &\le \overline Z_Q,\\
\|u-u_0\|_\infty &\le 2^{-34}.
\end{aligned}
\]

| Frozen candidate and outcome | \(\overline Z_Q\) | Room below one |
|:--|--:|--:|
[[DERIVATIVE_TABLE]]

The table and figure display rounded dimensionless bounds. The receipts retain exact rationals. The room below one leaves headroom for the spatial Hessian defect. The signed point-force condition is a separate test.

![Four reviewed cached derivative bounds at radius 2 to the minus 34. The bars show the room left below the contraction threshold of one.](cached_derivative_headroom.png)

[Four cached-root proofs]([[CACHE_RECEIPT]]) · [Whole-cube geometry replay]([[GEOMETRY_RECEIPT]]) · [Larger cached derivative replay]([[DERIV_RECEIPT]]) · [Figure provenance](figure_receipt.json)

## What the certified cached roots say about the outputs

The large-patch observable \(Q\) is the current apical height on the axis minus its height at current radius **[[SHOULDER]] µm**. The reported interval encloses \(Q_{+}-Q_{-}\), with the candidate order frozen before the nonlinear comparison. The independent transfer uses each root displacement bound plus the exact binary addition gap between the analytic center and serialized coordinates. A monotone quadratic trace on the common current-radius interval then bounds the apical shoulder value. This gap measures the response of two fixed model candidates. It gives no estimate of measurement error or experimental resolution.

For the **two corresponding L states**, the cached-root observer application independently checked **[[COMPARISONS]]** inequalities across four static conventions. **[[IN]] were IN and [[OUT]] were OUT** at the original allowance \(\delta \simeq [[DELTA]]\,\mu\mathrm{m}\). Every OUT comparison was the regional mean of total shortening under radial-line weighting. IN and OUT refer to whether a rigorous absolute residual interval lies within or above the frozen allowance. The allowance is a model-discrepancy threshold chosen for this diagnostic.

The full saved-field comparison is sensitive to observation convention. These worst-case common required allowances come from the separate, independently reviewed sweep of the **34 saved small fields**, using the same two material candidates and recorded values.

| Static convention | Common required allowance (µm) | At the frozen allowance |
|:--|--:|:--|
[[OBSERVER_TABLE]]

Only the declared primary convention, disk weighting with separate apical and basal spatial maxima, fits the same allowance for both saved-field candidates. The stored-field sweep leaves the source timing, initial profiles, specimen-level pairing and ROI registration unresolved.

[Exact cached observer readback]([[OBS_DERIVATION]]) · [Exact observer intervals]([[OBS_RECEIPT]]) · [Independent saved-field observer sweep]([[SWEEP_REVIEW]]) · [Full nonlinear report]([[FULL_REPORT]])

## What exact spatial integration must still prove

The next test keeps the same closed material singleton and binary \(R\). Let \(F_I\) and \(F_Q\) denote the exact-spatial and cached-quadrature reduced forces. At the exact center \(u_0\), the spatial force defect \(D_0\) and cached anchor must be combined with their signs intact.

\[
\begin{aligned}
D_0 &= F_I(u_0)-F_Q(u_0),\\
R F_I(u_0) &= R F_Q(u_0)+R D_0 .
\end{aligned}
\]

Over the **entire** matching cube, a bound \(b\) on \(\|R(H_I-H_Q)\|_\infty\) gives a spatial derivative bound \(z=\overline Z_Q+b\). Let \(\eta\) bound the infinity norm of the signed full center force \(R F_I(u_0)\). A sufficient strict local root test is

\[
z<1,\qquad \eta+zr<r .
\]

The plan calls for **144 bulk and 24 surface components** per outcome, with independent replay and signed global reduction. The retained source and synthetic checks verify the arithmetic interface. This checkpoint has no admitted exact-spatial root result for these witnesses. A spatial-root certificate would still leave the energy-minimum proof and continuum error estimate open.

[Independent component-replay method]([[COMPONENT_METHOD]]) · [Signed reduction derivation]([[REDUCTION_PROOF]]) · [Pure reduction source review]([[REDUCTION_REVIEW]])

## Coverage of the original five gates

The frozen comparison has **[[OUTCOMES]] outcomes**, comprising **34 small final fields and two large final fields**. The other eight saved checkpoints document numerical loading. At the four-pilot snapshot, **[[MISSING]] outcomes** lay outside the cached proof set. Even if all four exact-spatial pilots succeed, they cover just **[[COVERED]] of [[INEQUALITIES]]** primary small-region inequalities and the two large operands. **[[REMAINING]] primary inequalities** still need exact-root observation transfer.

**Cached coverage update.** A separate follow-up now certifies the two M states in the cached-quadrature model, giving **[[CACHED_TOTAL]] of [[OUTCOMES]]** cached roots and leaving **[[CACHED_REMAINING]]** at that stage. Candidate + at M passed on radius \(2^{-36}\). Candidate − at M was unresolved on radius \(2^{-32}\). Its independent replay passed on the narrower \(2^{-34}\) cube at the same center, drive and binary preconditioner. [Positive M proof]([[M_PLUS]]) · [Negative M narrow-cube proof]([[M_MINUS]])

1. **Originality.** Existing work already establishes related small-to-large transfer and mechanical discrimination strategies. The bounded audit records the prior art. A new general concept has not been established.
2. **Data feasibility.** Public summaries and observation conventions are documented. A matched geometry, reference-profile, timing, transported ROI and measurement-discrepancy contract still needs evidence.
3. **Model and identifiability.** The two fixed models supply a conditional comparison. Full physical parameter and geometry coverage, branch admissibility, and continuum error control remain open.
4. **Frozen prediction.** The pair and its static diagnostic were frozen. Revised uncertainty-propagated predictions and matched larger-region outcome scoring remain open.
5. **Contribution.** The positive cached contrast gap is quantitative evidence under a declared model and observer. It does not establish biological transfer, experimental detectability or novelty. The study's full contribution assessment remains open.

[36-outcome coverage ledger]([[COVERAGE]]) · [Full checkpoint readback](receipt.json) · [Source observation audit]([[SOURCE_AUDIT]]) · [Originality follow-up]([[ORIGINALITY]]) · [Five-gate contract]([[PLAN]])
"""

def markdown_text(zq,gap,delta,requirements,ledger,obs,m_followup):
    bounds=[]
    for k,z in zq.items():
        sign='+' if k.startswith('max_contrast_change') else '−'
        geometry='small L' if k.endswith('small_17') else 'large'
        bounds.append(f'| Candidate {sign} · {geometry} | {float(z):.6f} | {float(1-z):.6f} |')
    observer=[]
    for k,label in zip(VARIANTS,LABELS,strict=True):
        value=requirements[k]
        bound=Q(value['required_allowance_um_exact'][1])
        decision='Within' if value['status_at_frozen_allowance']=='ALLOWANCE_SATISFIED' else 'Exceeds'
        observer.append(f'| {label} | {float(bound):.6f} | {decision} |')
    lo,hi=outward(*gap,14)
    shoulder=Q(ledger["large_Q_definition"]["shoulder_radius_um_exact"])
    values={'GAP_LO':lo,'GAP_HI':hi,'DERIVATIVE_TABLE':'\n'.join(bounds),'OBSERVER_TABLE':'\n'.join(observer),
            'SHOULDER':f'{float(shoulder):.8f}','DELTA':f'{float(delta):.10f}',
            'COMPARISONS':str(len(obs['all_32_comparisons'])),'IN':str(obs['comparison_counts']['IN']),
            'OUT':str(obs['comparison_counts']['OUT']),
            'OUTCOMES':str(ledger['outcome_state_count']),'MISSING':str(ledger['remaining_new_cache_and_singleton_cases']),
            'CACHED_TOTAL':str(len(zq)+len(m_followup)),
            'CACHED_REMAINING':str(ledger['outcome_state_count']-len(zq)-len(m_followup)),
            'COVERED':str(ledger['primary_small_constraints_if_four_exact_pilots_succeed']),
            'INEQUALITIES':str(ledger['total_primary_small_constraints']),'REMAINING':str(ledger['still_uncovered_after_four_successes'])}
    paths={
        'CACHE_RECEIPT':CACHE/'independent_actual_preparation/actual_replay/receipt.json',
        'GEOMETRY_RECEIPT':EXACT/'independent_actual/replay/receipt.json','DERIV_RECEIPT':DERIV/'actual_replay_v2/receipt.json',
        'OBS_DERIVATION':OBS/'DERIVATION.md','OBS_RECEIPT':OBS/'receipt.json',
        'SWEEP_REVIEW':RUN/'observer_robustness/independent_actual/receipt.json','FULL_REPORT':RUN/'report/REPORT.html',
        'COMPONENT_METHOD':EXACT/'bridge_v1/independent_actual_preparation/ROOT-SOURCE-READBACK.md',
        'REDUCTION_PROOF':EXACT/'bridge_reduction_v1/PROOF.md',
        'REDUCTION_REVIEW':EXACT/'bridge_reduction_v1/independent_review/final_binding.json',
        'COVERAGE':LEDGER/'COVERAGE.md','SOURCE_AUDIT':ROOT/'docs/JOINT-OBSERVATION-CONTRACT.html',
        'ORIGINALITY':RUN/'originality_followup/README.md','PLAN':ROOT/'PLAN.md',
        'M_PLUS':M_PLUS/'max_contrast_change_small_09.json','M_MINUS':M_MINUS/'receipt.json',
    }
    tokens={**values,**{key:link(path) for key,path in paths.items()}}
    return re.sub(r'\[\[([A-Z_]+)\]\]',lambda match:tokens[match[1]],TEMPLATE)


CSS = """
*{box-sizing:border-box}
:root{color-scheme:dark;--bg:#0b1320;--text:#eaf3fc;--muted:#b4c5d9;--line:#344963;--cyan:#96e3f5;--gold:#ffca84;--surface:#122034}
body{margin:0;background:var(--bg);color:var(--text);font:17px/1.68 system-ui,sans-serif}
main{max-width:1000px;min-width:0;margin:0 auto;padding:28px clamp(16px,5vw,52px) 68px}
nav{display:flex;gap:5px 22px;flex-wrap:wrap;border-bottom:1px solid var(--line);padding:0 0 15px;font-size:.88rem}
a{color:var(--cyan);text-underline-offset:.2em}a:hover{color:#d6f8ff}
h1,h2{line-height:1.13;text-wrap:balance;letter-spacing:-.026em}
h1{font-size:clamp(2.05rem,5.8vw,3.55rem);max-width:20ch;margin:30px 0 17px}
h2{color:#d4f1fb;font-size:clamp(1.32rem,3.2vw,1.85rem);margin:1.95em 0 .55em}
p{overflow-wrap:break-word;margin:.9em 0}
h1+p{border-left:3px solid var(--gold);background:var(--surface);padding:13px 18px;border-radius:0 7px 7px 0}
strong{color:#f4faff}li{margin:.48em 0}ol{padding-left:1.5em}
img{width:100%;height:auto;display:block;border:1px solid var(--line);border-radius:9px;margin:20px 0}
table{position:relative;display:block;max-width:100%;overflow-x:auto;border-collapse:collapse;margin:1.1em 0;font-size:.89em;font-variant-numeric:tabular-nums}
th,td{padding:8px 11px;border:1px solid var(--line);vertical-align:top}
thead{color:#f3fbff;background:var(--surface)}
tbody tr:nth-child(even){background:#0f1a2a}
.math.display{max-width:100%;position:relative;overflow-x:auto;background:var(--surface);border-radius:7px;padding:7px 10px}
.math.display .katex-display{margin:.48em 0}
@media(max-width:520px){body{font-size:16px}main{padding-top:18px}h1{margin-top:23px}h1+p{padding:10px 13px}th,td{padding:6px 7px}.math .katex{font-size:1.08em}}
"""
MATH_SNIPPET = r'''
const fs = require("node:fs");
const katex = require(process.argv[1]);
const input = JSON.parse(fs.readFileSync(0, "utf8"));
const output = input.map(e => katex.renderToString(e.tex, {
  displayMode:e.display, throwOnError:true, strict:"error", output:"htmlAndMathml"
}));
process.stdout.write(JSON.stringify(output));
'''


def compile_math(body):
    pattern=re.compile(r'<(span|div) class="math( display)?" data-tex="([^"]+)"></\1>')
    matches=list(pattern.finditer(body))
    tasks=[{'tex':unescape(m[3]),'display':m[1]=='div'} for m in matches]
    require(len(tasks)>=12 and sum(t['display'] for t in tasks)==3,'expected all inline and display expressions')
    result=subprocess.run(['node','-e',MATH_SNIPPET,str(ASSETS/'katex.min.js')],
                          input=json.dumps(tasks),encoding='utf-8',capture_output=True,check=True,cwd=HERE)
    compiled=json.loads(result.stdout)
    require(len(compiled)==len(matches) and all('class="katex"' in line for line in compiled),
            'missing compiled KaTeX expressions')
    stream=iter(compiled)
    return pattern.sub(lambda m:f'<{m[1]} class="math{m[2] or ""}">{next(stream)}</{m[1]}>',body),len(matches)


class LocalLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.paths=[]

    def handle_starttag(self,tag,attrs):
        for key,val in attrs:
            if key in ('href','src'):
                self.paths.append(val)


def build_html(contents):
    sys.path.insert(0,str(ROOT/'scripts'))
    import build_report
    expected=ROOT/'scripts/build_report.py'
    require(Path(build_report.__file__).resolve()==expected and sha(expected)==PINS[expected],
            'renderer module differs from the frozen converter')
    body,math_count=compile_math(build_report.convert(contents))
    page=('<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
          '<meta name="viewport" content="width=device-width,initial-scale=1">'
          f'<title>{escape(TITLE)}</title><link rel="icon" type="image/png" href="{link(FAVICON)}">'
          f'<link rel="stylesheet" href="{link(ASSETS/"katex.min.css")}">'
          f'<style>{CSS}</style></head><body><main><nav aria-label="Related documents">'
          f'<a href="{link(ROOT/"docs/index.html")}">Controlled-furrowing dossier</a>'
          f'<a href="{link(RUN/"report/REPORT.html")}">Nonlinear report</a>'
          '<a href="CHECKPOINT.md">Markdown</a></nav>'+body+'</main></body></html>\n')
    return page,math_count


def write(path,text):
    with path.open('x',encoding='utf-8',newline='\n') as handle:
        handle.write(text)


def main():
    zq,gap,delta,requirements,ledger,obs,checked=reviewed_numbers()
    m_followup=reviewed_m_followup()
    contents=markdown_text(zq,gap,delta,requirements,ledger,obs,m_followup)
    page,math_count=build_html(contents)
    p=LocalLinks()
    p.feed(page)
    resolutions={}
    for url in p.paths:
        target=(HERE/url).resolve()
        require(target.is_relative_to(SITE) and (target.is_file() or target.name=='CHECKPOINT.md'),
                f'checkpoint link outside local site or missing: {url}')
        resolutions[url]=site_name(target)
    write(HERE/'CHECKPOINT.md',contents)
    write(HERE/'CHECKPOINT.html',page)
    files={name:sha(HERE/name) for name in ('CHECKPOINT.md','CHECKPOINT.html')}
    receipt={'status':'BUILD_FROM_PINNED_REVIEWED_CACHED_QUADRATURE_RESULTS','file_sha256':files,
             'source_sha256':sha(Path(__file__)),'input_sha256':{site_name(p):h for p,h in PINS.items()},
             'root_readback_inputs_verified':checked,'compiled_math_blocks':math_count,'links':resolutions,
             'derived_exact_values':{'ZQ_upper_bound_by_pilot':{p:str(z) for p,z in zq.items()},
                 'room_below_one_by_pilot':{p:str(1-z) for p,z in zq.items()},
                 'cached_large_gap_um_exact':list(map(str,gap)),'cached_large_gap_um_outward_display':outward(*gap,14),
                 'delta_um_exact':str(delta),'small_cached_comparisons':obs['comparison_counts'],
                 'same_saved_fields_common_allowances_um_exact':{name:requirements[name]['required_allowance_um_exact']
                                                                  for name in VARIANTS}},
             'original_cached_quad_pilot_count':len(zq),'M_followup_cached_cubes_radius_exact':m_followup,
             'updated_cached_quad_roots':len(zq)+len(m_followup),'exact_spatial_roots_in_checkpoint':0,
             'mathematical_operator_verified':'frozen binary Q8/S12 cached quadrature',
             'new_derivative_force_integral_or_solver_work':0,'full_five_gate_complete':False}
    path=HERE/'checkpoint_build_receipt.json'
    write(path,json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'receipt_sha256':sha(path),
                      'compiled_math':math_count,'files':files}),flush=True)


if __name__=='__main__':
    main()
