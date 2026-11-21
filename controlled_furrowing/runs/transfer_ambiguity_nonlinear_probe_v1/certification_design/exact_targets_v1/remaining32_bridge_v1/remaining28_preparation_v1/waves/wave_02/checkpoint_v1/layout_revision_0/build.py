"""Wave 02 checkpoint from separately observed exact root and observer evidence."""

import argparse
from collections import Counter
import csv
from fractions import Fraction as Q
from hashlib import sha256
from html import escape
from importlib.util import module_from_spec, spec_from_file_location
import json
import math
import os
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter, MultipleLocator


HERE = Path(__file__).resolve().parent
WAVE = HERE.parent
ROOT = next(p for p in HERE.parents if (p / 'scripts/build_report.py').is_file())
SITE = ROOT.parent
EX = next(p for p in HERE.parents if p.name == 'exact_targets_v1')
BR = EX / 'remaining32_bridge_v1'
NO = WAVE / 'global_reduction_v1/independent'
OB = WAVE / 'observer_application_v1/independent_actual_readback'
KENERGY = EX / 'negative_KL_active_energy_v1/independent_actual'
FIRST = BR / 'first_wave_reduction_v1/independent_actual_preparation/actual_replay'
FIRST_OB = BR / 'first_wave_observer_v1/independent_readback/actual/receipt.json'
FOUR = EX / 'bridge_reduction_v1/independent_actual_preparation/all_four_readback.json'
FOUR_OB = EX / 'observer_application_v2/independent_actual_readback/receipt.json'
CACHED = EX.parent / 'full_outcome_coverage_v1/observer_application_v1/independent_readback/receipt.json'
HELPER = EX / 'certification_checkpoint_v1/build_checkpoint.py'
CONVERT = ROOT / 'scripts/build_report.py'
ASSETS = SITE / 'docs/assets/katex'
FAV = EX.parents[1] / 'report/favicon.png'
TITLE = 'Wave 02 exact-spatial roots and observer bounds'
VARIANTS = ('current_disk__spatial_max','current_disk__regional_mean',
            'current_radial_line__spatial_max','current_radial_line__regional_mean')
VNAMES = ('Disk, component maxima','Disk, component means','Radial line, component maxima','Radial line, component means')
CASE_NAMES = ('+Q','+M','−Q','−M')
CASE_IDS = tuple('max_' + sign + 'contrast_change_small_' + step
                 for sign in ('','negative_') for step in ('08','09'))
RADIUS = Q(1,2**34)
PINS = {
    NO / 'results.json': '1a2a785c5bceaf4d301295162747a63cfec4c4da235ba709a46c8e383d04b411',
    NO / 'root_observed_exit.json': 'b6f8b9014b397d8da5b962cd757993d1060096cadeaa8a76d7e53c29d461e082',
    OB / 'receipt.json': '88d5950f9b0ec2cc236d57c9b70f686d11fcb4a1879870e2c36e91eecc578fc0',
    OB / 'os_exit_observation.json': 'b9d886ba79bc2e8a9cd59b969fdf07cc959fff78331ab6bf072a93b74cfa4d56',
    WAVE / 'observer_application_v1/actual/results.json':
        '5356604ff5c6314dea0a393348aa9855199a41ebac7343f01b62dfe9a4e1f6e3',
    FIRST / 'receipt.json': 'e73ea05a4703256235a46ce7b914298d3f8258e9fef029f4dbe952b04bb9798a',
    FIRST / 'os_exit_observation.json': 'aa3c90101c64ccefe1403af2c296bcb22cafad3b2e35d73539aa7e3c5227f281',
    FIRST_OB: '304e8d6df8da8ca8c91737421d42c26384c2610a188736d79f327456a9ad7bd6',
    FOUR: '859f507ad3762cfef44d53672bbe15a43592ce54d6e08ba8828bf50a925ec340',
    FOUR_OB: 'e57c6e18e4e585b002ebab8aa0782d063e281c201ab1449cfac484ef325287e7',
    CACHED: 'a06549f048abba4eb2b4d1c1fb398f06e8f524196c92a673bd0493c527c8f40d',
    KENERGY / 'receipt.json': '84bbe79fac39fbf0ce72b45b9e0c18f76310cad981da004c5f0b5e2a7f6d1abf',
    KENERGY / 'os_exit_observation.json': '060f31afa87c5571760ba0d57c4c51052bafbd2994d76a3328951c5dcef89316',
    HELPER: '50524d01db5003b0839e873efb37d41bb0552ad0f0c162f73bd668ad0e58d31b',
    CONVERT: '90148d53d74d8ea5cc8a665e6a3660abf2db72487a1fa59586caba98f9cc38e2',
    ASSETS / 'katex.min.css': '5a0ff74f11887f7260dc16b216906f51352b84ca29f21abb9a14b62a2ca6eb15',
    ASSETS / 'katex.min.js': 'e6bfe5deebd4c7ccd272055bab63bd3ab2c73b907b6e6a22d352740a81381fd4',
    FAV: '961e6e85ec16258652bf3c3f322d336fb598b90bccc49821047babfef3925530',
}
BG, GRID, PALE, INK, CYAN, BLUE, RED = '#122034','#3a516a','#c1d1e4','#f1f6ff','#80d9e8','#9caef1','#ffb6a0'


def demand(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def link(path):
    return Path(os.path.relpath(path,HERE)).as_posix()


def pair(seq):
    a,b=map(Q,seq)
    demand(a<=b,'ordered interval')
    return a,b


def directed(value, places, up):
    x=value*10**places
    q=x.__ceil__() if up else x.__floor__()
    whole,fraction=divmod(abs(q),10**places)
    return f'{"−" if q<0 else ""}{whole}.{fraction:0{places}d}'


def zero(v):
    return type(v) is int and v == 0


def text(path, data):
    with path.open('x', encoding='utf-8', newline='\n') as f:
        f.write(data)


def binary_outward(x, high):
    f = float(x)
    if (Q(f) < x if high else Q(f) > x):
        f = math.nextafter(f, math.inf if high else -math.inf)
    return f


def load(path, pins, expected=None):
    p=path.resolve()
    demand(p.is_relative_to(SITE),'local immutable authority')
    data=p.read_bytes()
    digest=sha256(data).hexdigest()
    demand(digest==pins.setdefault(p,expected or digest) and (expected is None or digest==expected),
           'changed evidence: '+str(p))
    return json.loads(data)


def residue(row, allowance):
    lower,upper=pair(row['exact_spatial_root_observer_interval_um_exact'])
    measurement=Q(row['observed_um_exact'])
    x,y=lower-measurement,upper-measurement
    lo,hi=max(Q(0),x,-y),max(abs(x),abs(y))
    state='IN' if hi<=allowance else 'OUT' if lo>allowance else 'UNRESOLVED'
    demand((lo,hi)==pair(row['absolute_residual_um_exact']) and row['standing']==state,'actual signed interval mismatch')
    return lo,hi,state


def checks(args):
    pins=PINS.copy()
    demand(all(sha(p)==h for p,h in pins.items()),'frozen evidence and renderer changed')
    no,no_os=load(NO/'results.json',pins),load(NO/'root_observed_exit.json',pins)
    ob,ob_os=load(OB/'receipt.json',pins),load(OB/'os_exit_observation.json',pins)
    ob_term=load(OB/'terminal.json',pins,ob_os['terminal_sha256'])
    ob_started=load(OB/'started.json',pins,ob_os['started_sha256'])
    adopted=load(args.observer_acceptance,pins,args.observer_acceptance_sha256)
    demand(adopted['receipt_sha256']==pins[OB/'receipt.json']
           and adopted['agent_true_os_sha256']==pins[OB/'os_exit_observation.json']
           and adopted['status']=='PASS_ROOT_ACCEPTED_WAVE02_INDEPENDENT_EXACT_SPATIAL_OBSERVERS'
           and adopted['all64_scalar_residual_statuses_rechecked'] is True and adopted['full_five_gate_complete'] is False,
           'explicit accepted observer replay and true-process authority')
    demand(no['status']=='COMPLETE_INDEPENDENT_LATER_WAVE_SIGNED_ENDPOINT_REPLAY'
           and no_os['status']=='PASS_ROOT_OBSERVED_LATER_WAVE_INDEPENDENT_GLOBAL_EXIT'
           and no_os['results_sha256']==pins[NO/'results.json'] and zero(no_os['observed_exit_code'])
           and no_os['direct_root_tool_observation'] is True and no_os['all672_local_hull_links_rehashed'] is True
           and no['wave_id']==no_os['wave_id']==ob['wave_id']==ob_os['wave_id']=='wave_02'
           and len(no_os['sha256'])==11 and no['signed_hull_components']==672 and no['full_five_gate_complete'] is False,
           'complete observed global native wave')
    for name,digest in no_os['sha256'].items():
        p=(ROOT/name).resolve()
        demand(p.is_relative_to(NO) and sha(p)==digest and pins.setdefault(p,digest)==digest,'native global closure mismatch')
    demand(ob['status']=='PASS_INDEPENDENT_ACTUAL_LATER_WAVE_SPATIAL_OBSERVER_READBACK'
           and ob_os['status']=='PASS_OBSERVED_INDEPENDENT_LATER_WAVE_SPATIAL_OBSERVER_READBACK_EXIT'
           and ob_os['receipt_sha256']==ob_term['receipt_sha256']==pins[OB/'receipt.json']
           and zero(ob_os['observed_exit_code']) and zero(ob_term['exit_code'])
           and all(type(x['pid']) is int and x['pid']==ob_os['pid']>0 for x in (ob,ob_term,ob_started))
           and ob['native_global_os_sha256']==ob_os['native_global_os_sha256']==pins[NO/'root_observed_exit.json']
           and ob['observer_producer_results_sha256']==pins[WAVE/'observer_application_v1/actual/results.json']
           and ob['all_64_exact_comparisons_recreated'] is True and ob['complete_calibration_or_source_fit_claimed'] is False,
           'complete independent observer proof and process')
    demand(tuple(r['case_id'] for r in no['cases'])==tuple(r['case_id'] for r in ob['cases'])==CASE_IDS,
           'the four fixed wave 02 targets in profile order')
    allowance=Q(ob['allowance_um_exact'])
    roots, comps = [], []
    for batch, observer, name in zip(no['cases'],ob['cases'],CASE_NAMES,strict=True):
        case=batch['case_id']
        native=load(NO/case/'results.json',pins,batch['case_results_sha256'])
        operands=load(NO/case/'operands.json',pins,batch['operands_sha256'])
        a=native['arithmetic']
        t=native['observer_and_optional_minimum']
        demand(batch['exact_spatial_root_certified'] is True and batch['exact_spatial_minimum_certified'] is True
               and native['status']=='PASS_INDEPENDENT_LATER_WAVE_EXACT_SPATIAL_ROOT' and native['root_asserted'] is True
               and native['target']==operands['target'] and Q(native['radius_exact'])==RADIUS
               and observer['root_proof_sha256']==batch['case_results_sha256'] and t['strict_minimum_proved'] is True
               and observer['original_field_sha256']==native['target']['original_observer_field_sha256']
               and observer['record_id'].endswith('_'+name[-1]),'binding exact wave root, same original field, minimum')
        eta,z=Q(a['eta_upper_exact']),Q(a['z_upper_exact'])
        e,sigma,xi=map(Q,(observer['posterior_exact'],observer['archived_sigma_exact'],observer['original_xi_exact']))
        lhs=eta+z*RADIUS
        ep=Q(67,2)*(e+sigma+xi)
        m=Q(t['normalized_H_euclidean_lower_exact'])
        demand(0<=eta and 0<=z<1 and lhs==Q(a['selfmap_lhs_upper_exact'])<RADIUS and a['strict_selfmap'] is True
               and Q(a['ZQ_upper_exact'])+Q(a['b_H_upper_exact'])==z
               and 0<=e==eta/(1-z)==Q(a['posterior_error_upper_exact'])<RADIUS
               and m==(1-z)/Q(operands['R_infinity_norm_exact'])>0
               and sigma==Q(native['target']['archived_center_shift_exact'])
               and xi==Q(native['target']['original_addition_gap_exact'])
               and ep==Q(t['physical_nodal_error_um_upper_exact'])==Q(observer['physical_nodal_error_um_exact']),
               'exact contraction, curvature and original observer posteriors')
        roots.append({'case_id':case,'label':name,'record_id':observer['record_id'],'radius_exact':str(RADIUS),
            'eta_upper_exact':str(eta),'eta_over_radius_exact':str(eta/RADIUS),'z_upper_exact':str(z),
            'selfmap_headroom_fraction_exact':str(1-lhs/RADIUS),'curvature_lower_exact':str(m),
            'posterior_exact':str(e),'archived_sigma_exact':str(sigma),'original_xi_exact':str(xi),
            'physical_nodal_error_um_exact':str(ep),'exact_spatial_proof_sha256':batch['case_results_sha256']})
        demand(len(observer['comparisons'])==16 and len({(r['variant'],r['observable']) for r in observer['comparisons']})==16,
               'complete independent case comparisons')
        for comp in observer['comparisons']:
            lo,hi,state=residue(comp,allowance)
            p0,p1=pair(comp['exact_spatial_root_observer_interval_um_exact'])
            comps.append({'case_id':case,'label':name,'variant':comp['variant'],'observable':comp['observable'],
                'prediction_lo_um_exact':str(p0),'prediction_hi_um_exact':str(p1),
                'observed_um_exact':comp['observed_um_exact'],'residual_lo_um_exact':str(lo),'residual_hi_um_exact':str(hi),
                'allowance_um_exact':str(allowance),'standing':state,
                'strict_OUT_excess_um_exact':str(lo-allowance) if state=='OUT' else ''})
    count=Counter((r['variant'],r['standing']) for r in comps)
    expect=[{'IN':16},{'IN':14,'OUT':2},{'IN':16},{'IN':14,'OUT':2}]
    demand([dict(Counter({s:count[v,s] for s in ('IN','OUT','UNRESOLVED') if count[v,s]})) for v in VARIANTS]
           ==expect==[ob['counts'][v] for v in VARIANTS],'exact wave observer standings')
    out=[r for r in comps if r['standing']=='OUT']
    demand(len(out)==len(ob['strict_out'])==4 and all(r['label'] in ('+M','−M') and
           r['observable']=='apical_max_shortening_um' and r['variant'] in VARIANTS[1::2] for r in out),
           'four alternative M failures and no primary failure')
    for r in out:
        [match]=[x for x in ob['strict_out'] if (x['case_id'],x['variant'],x['observable'])==
                 (r['case_id'],r['variant'],r['observable'])]
        demand(Q(match['strict_OUT_margin_um_exact'])==Q(r['strict_OUT_excess_um_exact']),'exact strict OUT margin')

    a,ao=load(FOUR,pins),load(FOUR_OB,pins)
    b,bo=load(FIRST/'receipt.json',pins),load(FIRST_OB,pins)
    b_os=load(FIRST/'os_exit_observation.json',pins)
    ids_a={r['pilot_id'] for r in a['cases']}
    ids_b={r['case_id'] for r in b['case_results']}
    ids_c={r['case_id'] for r in roots}
    cached=load(CACHED,pins)
    cached_ids={r['case_id'] for r in cached['observer_transfer_terms']}
    demand(len(ids_a)==len(ids_b)==len(ids_c)==4 and ids_a.isdisjoint(ids_b)
           and ids_a.isdisjoint(ids_c) and ids_b.isdisjoint(ids_c)
           and ids_a=={r['pilot_id'] for r in ao['cases']} and ids_b=={r['case_id'] for r in bo['cases']}
           and b_os['receipt_sha256']==pins[FIRST/'receipt.json'] and zero(b_os['observed_exit_code'])
           and a['status']=='PASS_FOUR_INDEPENDENT_FIXED_INPUT_EXACT_SPATIAL_ROOTS_AND_MINIMA'
           and a['independent_exact_spatial_root_count']==b['exact_spatial_roots']==b['strict_local_minima']==4
           and len(cached_ids)==cached['cached_fixed_parameter_roots']==36 and (ids_a|ids_b|ids_c)<=cached_ids
           and Q(ao['conditional_allowance_um_exact'])==Q(bo['allowance_um_exact'])==allowance,
           'all three prior and current populations are genuinely disjoint with common frozen allowance')
    total=count.copy()
    for cases,key in ((ao['cases'],'status'),(bo['cases'],'standing')):
        for case in cases:
            for row in case.get('comparisons',[]):
                total[row['variant'],row[key]]+=1
    total_in=[total[v,'IN'] for v in VARIANTS]
    total_out=[total[v,'OUT'] for v in VARIANTS]
    demand(total_in==[40,37,38,36] and total_out==[0,3,2,4] and sum(total.values())==160,
           'all 160 cumulative rows and verified group coverage')
    union=ids_a|ids_b|ids_c
    demand(len([c for c in union if '_small_' in c])==10 and len([c for c in union if '_large_' in c])==2,
           'ten exact small and two exact large fixed inputs')

    ke,kos=load(KENERGY/'receipt.json',pins),load(KENERGY/'os_exit_observation.json',pins)
    accept_energy=load(args.energy_acceptance,pins,args.energy_acceptance_sha256)
    demand(accept_energy['status']=='PASS_ROOT_ACCEPTED_INDEPENDENT_NEGATIVE_KL_NECESSARY_ORDER'
           and accept_energy['independent_receipt_sha256']==pins[KENERGY/'receipt.json']
           and accept_energy['agent_true_os_sha256']==pins[KENERGY/'os_exit_observation.json']
           and accept_energy['independent_difference_strictly_negative'] is True and accept_energy['certified_loading_path'] is False
           and kos['receipt_sha256']==pins[KENERGY/'receipt.json'] and zero(kos['observed_exit_code'])
           and ke['status']=='PASS_INDEPENDENT_ACTUAL_NEGATIVE_KL_NECESSARY_ENERGY_ORDER'
           and ke['case_count']==2 and ke['active_elements_each']==24
           and ke['positive_common_factor_2pi_over_scale'] is True and ke['successful_necessary_order_only'] is True
           and ke['certified_loading_path'] is False,'accepted necessary active-energy result only')
    kgap=pair(ke['difference']['independent_interval_exact'])
    demand(kgap[1]<0 and [x['case_id'] for x in ke['two_endpoint_proofs']]
           ==['max_negative_contrast_change_small_16','max_negative_contrast_change_small_17'],'correct energy sign and direction')
    cover={'original_four':sorted(ids_a),'first_wave':sorted(ids_b),'wave_02':sorted(ids_c),
           'cumulative_root_count':len(union),'exact_spatial_small_fields':10,'exact_spatial_large_fields':2,
           'counts_by_variant':{v:{'IN':i,'OUT':o} for v,i,o in zip(VARIANTS,total_in,total_out,strict=True)},
           'remaining_small_fields':24,'remaining_primary_inequalities':96}
    return roots,comps,out,expect,allowance,kgap,cover,pins


def csv_file(name, rows):
    path=HERE/name
    with path.open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    with path.open(encoding='utf-8',newline='') as f:
        demand(list(csv.DictReader(f))==rows,'exact rational table writeback')


def chart_axes(size):
    fig,ax=plt.subplots(figsize=size,dpi=180,facecolor=BG)
    ax.set_facecolor(BG)
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['bottom','left']].set_color(GRID)
    ax.tick_params(colors=PALE,labelsize=11)
    ax.grid(axis='x',color=GRID,alpha=.62)
    ax.set_axisbelow(True)
    ax.invert_yaxis()
    return fig,ax


def figures(rows, out, allowance):
    fig,ax=chart_axes((7.8,3.45))
    fig.subplots_adjust(left=.13,right=.96,top=.83,bottom=.22)
    e=[Q(r['eta_over_radius_exact']) for r in rows]
    z=[Q(r['z_upper_exact']) for r in rows]
    positions=range(4)
    for y,(a,b) in enumerate(zip(e,z,strict=True)):
        aa=ax.add_patch(Rectangle((0,y-.28),float(a),.56,color=CYAN,label='Center residual η/r' if y==0 else None))
        bb=ax.add_patch(Rectangle((float(a),y-.28),float(b),.56,color=BLUE,label='Whole-cube derivative z' if y==0 else None))
        demand(aa.get_width()==float(a) and bb.get_width()==float(b) and bb.get_x()==float(a),
               'plotted exact self-map terms do not match')
        ax.text(float(a+b)+.025,y,f'{float(a+b):.3f}',va='center',color=INK,fontsize=11)
    ax.set_xlim(0,1.03)
    ax.set_ylim(3.6,-.58)
    ax.set_yticks(positions,[r['label'] for r in rows])
    ax.xaxis.set_major_locator(MultipleLocator(.2))
    ax.axvline(1,color=RED,linestyle='--',linewidth=1.7)
    ax.set_xlabel('Fraction of cube radius in the self-map bound',color=PALE,fontsize=11)
    ax.legend(loc='upper left',bbox_to_anchor=(0,1.21),ncols=2,frameon=False,labelcolor=INK,fontsize=10)
    ax.text(1.0,-.68,'limit 1',ha='center',color=RED,fontsize=10)
    with (HERE/'selfmap_wave02.png').open('xb') as f:
        fig.savefig(f,format='png')
    plt.close(fig)

    fig,ax=chart_axes((7.8,3.55))
    fig.subplots_adjust(left=.28,right=.95,top=.87,bottom=.215)
    labels=[]
    glyphs=[]
    threshold=binary_outward(allowance,True)
    ax.axvline(threshold,color=CYAN,linestyle='--',linewidth=2)
    span=[]
    for y,row in enumerate(out):
        a,b=map(Q,(row['residual_lo_um_exact'],row['residual_hi_um_exact']))
        x0,x1=binary_outward(a,False),binary_outward(b,True)
        center=float((a+b)/2)
        line,=ax.plot([x0,x1],[y,y],color=RED,linewidth=2.5)
        point,=ax.plot(center,y,color=RED,marker='o',markersize=7)
        demand(list(line.get_xdata())==[x0,x1] and Q(x0)<=a<=b<=Q(x1)
               and point.get_xdata()[0]==center and Q(threshold)>=allowance and Q(threshold)<a,
               'plotted threshold and directed residual interval invalid')
        ax.text(center+.025,y,f'{center:.3f}',va='center',color=INK,fontsize=11)
        labels.append(row['label']+' · '+('disk mean' if 'disk' in row['variant'] else 'line mean'))
        glyphs.append({'case_id':row['case_id'],'variant':row['variant'],
                       'exact_interval':list(map(str,(a,b))),'plot_binary_endpoints':[x0,x1]})
        span.append(x1)
    ax.set_xlim(threshold-.17,max(span)+.19)
    ax.set_yticks(range(4),labels)
    ax.xaxis.set_major_locator(MultipleLocator(.25))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlabel('Absolute residual in the apical regional-mean channel (µm)',color=PALE,fontsize=10)
    ax.text(threshold,-.78,'frozen allowance',ha='left',color=CYAN,fontsize=10)
    with (HERE/'observer_failures_wave02.png').open('xb') as f:
        fig.savefig(f,format='png')
    plt.close(fig)
    return {'exact_mechanics_terms':[[str(a),str(b)] for a,b in zip(e,z,strict=True)],
            'out_residual_captures':glyphs,'threshold_plot_binary_upper':threshold}


def markdown(roots, out, count, allowance, kgap, cover):
    m='\n'.join(f"| {r['label']} | {float(Q(r['z_upper_exact'])):.5f} | "
        f"{100*float(Q(r['selfmap_headroom_fraction_exact'])):.2f}% | "
        f"{float(Q(r['curvature_lower_exact'])):.9f} | {float(Q(r['physical_nodal_error_um_exact'])):.4g} |" for r in roots)
    obs='\n'.join(f"| {name} | {c['IN']} | {c.get('OUT',0)} | {n['IN']} | {n['OUT']} |"
         for name,c,n in zip(VNAMES,count,cover['counts_by_variant'].values(),strict=True))
    low,high=(directed(x,16,upper) for x,upper in zip(kgap,(False,True),strict=True))
    return fr"""# {TITLE}

Four more fixed-input equilibria have independently verified exact-spatial roots and strict local energy minima. The count is now **12 of 36** after matching case IDs across the original four, the first new wave and wave 02. This includes ten small fields and two large fields. All 36 already have cached-quadrature proofs. [Wave 02 signed reduction]({link(NO/'results.json')}) · [First-wave report]({link(BR/'first_wave_checkpoint_v1/REPORT.html')})

## Root and minimum proofs

The four new cases use small records Q and M at both frozen material choices. The + and − labels name the two material choices. For each cube of radius \(r=2^{{-34}}\), \(\eta\) bounds the signed preconditioned center residual and \(z\) bounds the full-cube derivative. The independently verified local self-map requires

\[
\eta/r+z<1.
\]

The normalized curvature lower bound is \((1-z)/\lVert R\rVert_\infty\) using the verified positive symmetric part of the same \(R\). Positive exact body and trace domains apply on each certified cube. Headroom below is \(1-z-\eta/r\). The observer error uses the root posterior \(e\), archived center shift \(\sigma\), and original binary addition gap \(\xi\).

\[
\epsilon=33.5(e+\sigma+\xi)\,\mathrm{{\mu m}}.
\]

| Case | \(z\) | Headroom | Curvature floor | Nodal error (µm) |
|:--|--:|--:|--:|--:|
{m}

![For four wave 02 roots, the stacked center residual and derivative terms occupy between 0.648 and 0.682 of the unit self-map budget.](selfmap_wave02.png)

The table shows rounded displays. Exact proof endpoints determine all decisions. [Exact bounds CSV](exact_bounds.csv) · [Observed independent process]({link(NO/'root_observed_exit.json')})

## Observer comparisons at the same allowance

The frozen maximum residual allowance is approximately **{float(allowance):.9f} µm**. It expresses a model discrepancy threshold. The wave 02 independent scalar replay recovered **16 primary IN**, and **44 IN / 4 OUT** across the three alternatives. All four OUT rows are comparisons of the apical regional mean at M. The plotted intervals are far narrower than one pixel.

| Convention | Wave 02 IN | OUT | Cumulative IN | OUT |
|:--|--:|--:|--:|--:|
{obs}

![The four M apical regional-mean residuals exceed the frozen allowance of about 5.39 micrometres. Their positive-direction values are about 6.58 for disk and 6.12 for line weighting, and the negative-direction values are about 6.42 and 5.95.](observer_failures_wave02.png)

The cumulative exact small-field coverage consists of **40 primary IN** and **111 alternative IN / 9 OUT**. Twenty-four small fields and 96 primary inequalities remain for the specified exactly integrated two-candidate comparison. Timing, a source-matched large shape outcome, and an experimental error model remain open. [All 64 exact wave comparisons](exact_observer_comparisons.csv) · [Independent scalar receipt]({link(OB/'receipt.json')})

## Separate K–L loading check

At the already certified negative-material K and L states, the same finite-element law and fixed 24-element active support give a negative signed difference in the coefficient \(V\) of activation energy. With the common positive normalization factor removed, the independent bound encloses

\[
\begin{{aligned}}
\frac{{\Delta V}}{{2\pi/S}} &\in [{low.replace('−','-')},\\
&\hspace{{1em}}{high.replace('−','-')}].
\end{{aligned}}
\]

On a differentiable branch of strict minima of \(E_0+gIV\), the derivative is \(dV/dI=-g\nabla V^T H^{{-1}}\nabla V\le 0\). The observed endpoint order satisfies this necessary condition. Intermediate roots, branch continuity, and loading from zero still need separate certificates. The two endpoints alone also leave biological dynamics, continuum accuracy, and source compatibility unresolved. [Independent active-energy calculation]({link(KENERGY/'receipt.json')}) · [Loading-tube design]({link(ROOT/'runs/exact_loading_tube_design_v1/METHOD.md')})
"""


STYLE = """
:root{color-scheme:dark;--b:#0b1320;--p:#122034;--l:#3a526e;--ink:#edf6ff;--a:#92e4f0}
*{box-sizing:border-box}body{margin:0;background:var(--b);color:var(--ink);font:17px/1.62 system-ui,sans-serif}
main{max-width:920px;min-width:0;margin:auto;padding:24px clamp(16px,5vw,44px) 56px}
nav{display:flex;gap:7px 18px;flex-wrap:wrap;border-bottom:1px solid var(--l);padding:0 0 12px;font-size:.85rem}
a{color:var(--a);text-underline-offset:.18em;overflow-wrap:anywhere}h1,h2{line-height:1.14;text-wrap:balance}
h1{max-width:21ch;letter-spacing:-.035em;font-size:clamp(2rem,6vw,3.5rem);margin:.55em 0 .32em}
h1+p{border-left:3px solid var(--a);padding:10px 14px;background:var(--p);border-radius:0 7px 7px 0}
h2{color:#cef0fc;font-size:clamp(1.4rem,4vw,1.94rem);margin:1.65em 0 .55em}
img{display:block;width:100%;height:auto;border:1px solid var(--l);border-radius:7px;margin:17px 0}
table{position:relative;max-width:100%;display:block;overflow-x:auto;border-collapse:collapse;font-variant-numeric:tabular-nums;margin:1em 0}
th,td{border:1px solid var(--l);padding:.43em .64em;vertical-align:top;text-align:left}thead{background:var(--p)}
td:nth-child(n+2){white-space:nowrap}.math.display{position:relative;overflow-x:auto;max-width:100%;background:var(--p);padding:6px 9px;border-radius:6px}
@media(max-width:520px){body{font-size:16px}main{padding-top:16px}th,td{font-size:.82em;padding:.42em .3em}.math.display .katex{font-size:1.04em}}
"""


def render(md):
    spec=spec_from_file_location('_pinned_static_katex',HELPER)
    helper=module_from_spec(spec)
    spec.loader.exec_module(helper)
    sys.path.insert(0,str(CONVERT.parent))
    import build_report
    demand(Path(build_report.__file__).resolve()==CONVERT,'wrong local converter')
    body,n=helper.compile_math(build_report.convert(md))
    html=('<!doctype html>\n<html lang="en"><head><meta charset="utf-8"><meta name="viewport" '
          'content="width=device-width,initial-scale=1">'
          f'<title>{escape(TITLE)}</title><link rel="icon" type="image/png" href="{link(FAV)}">'
          f'<link rel="stylesheet" href="{link(ASSETS/"katex.min.css")}"><style>{STYLE}</style></head>'
          '<body><main><nav aria-label="Report resources"><a href="REPORT.md">Markdown</a>'
          '<a href="exact_bounds.csv">Exact bounds</a><a href="exact_observer_comparisons.csv">Exact observers</a>'
          '</nav>'+body+'</main></body></html>\n')
    p=helper.LocalLinks()
    p.feed(html)
    links={}
    for f in p.paths:
        path=(HERE/f).resolve()
        demand(path.is_relative_to(SITE) and path.is_file(),'missing local report resource: '+f)
        links[f]=path.relative_to(SITE).as_posix()
    return html,n,links


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for field in ('observer','energy'):
        p.add_argument('--'+field+'-acceptance',type=Path,required=True)
        p.add_argument('--'+field+'-acceptance-sha256',required=True)
    args=p.parse_args()
    roots,comps,out,counts,allowance,kgap,cover,pins=checks(args)
    csv_file('exact_bounds.csv',roots)
    csv_file('exact_observer_comparisons.csv',comps)
    plot=figures(roots,out,allowance)
    md=markdown(roots,out,counts,allowance,kgap,cover)
    text(HERE/'REPORT.md',md)
    html, n, links=render(md)
    text(HERE/'REPORT.html',html)
    outputs=('REPORT.md','REPORT.html','exact_bounds.csv','exact_observer_comparisons.csv','selfmap_wave02.png','observer_failures_wave02.png')
    receipt={'status':'PASS_BUILD_FROM_OBSERVED_INDEPENDENT_WAVE02_READBACK','source_sha256':sha(Path(__file__)),
             'input_sha256':{p.relative_to(SITE).as_posix():h for p,h in pins.items()},'local_links':links,
             'file_sha256':{name:sha(HERE/name) for name in outputs},'rendered_math_count':n,'exact_coverage':cover,
             'wave_observer_counts':dict(zip(VARIANTS,counts,strict=True)),'independently_checked_exact_CSV_counts':[4,64],
             'energy_difference_normalization_removed_interval_exact':list(map(str,kgap)),
             'verified_plot_inputs':plot,'certified_joint_loading_path':False,'full_five_gate_complete':False}
    text(HERE/'build_receipt.json',json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'receipt_sha256':sha(HERE/'build_receipt.json'),'math':n,'output_sha256':receipt['file_sha256']}))


if __name__=='__main__':
    main()
