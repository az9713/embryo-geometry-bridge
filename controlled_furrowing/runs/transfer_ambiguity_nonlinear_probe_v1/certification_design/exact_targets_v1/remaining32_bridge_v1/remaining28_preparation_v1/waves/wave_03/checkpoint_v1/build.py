"""Read accepted Wave 03 proof inputs and render a local checkpoint report."""

import argparse
from collections import Counter
import csv
from fractions import Fraction as Q
from hashlib import sha256
from html import escape, unescape
from importlib.util import module_from_spec, spec_from_file_location
import json
import os
import re
import subprocess
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle


HERE=Path(__file__).resolve().parent
WAVE=HERE.parent
ROOT=next(p for p in HERE.parents if (p/'scripts/build_report.py').is_file())
SITE=ROOT.parent
EX=next(p for p in HERE.parents if p.name=='exact_targets_v1')
PREV=WAVE.parent/'wave_02/checkpoint_v1'
GLOBAL=WAVE/'global_reduction_v1/independent'
OB=WAVE/'observer_application_v1'
OB_INDEPENDENT=OB/'independent_actual_readback'
TRANS=ROOT/'runs/nonaxisymmetric_surface_sign_preparation_v1'
THEORY=ROOT/'runs/nonaxisymmetric_adequacy_audit_v1'
LOAD=ROOT/'runs/exact_loading_tube_pilot_geometry_actual_readback_v1'
CACHED=EX.parent/'full_outcome_coverage_v1/observer_application_v1/independent_readback/receipt.json'
HELPER=EX/'certification_checkpoint_v1/build_checkpoint.py'
CONVERT=ROOT/'scripts/build_report.py'
ASSETS=SITE/'docs/assets/katex'
FAV=EX.parents[1]/'report/favicon.png'
TITLE='Wave 03 exact-spatial roots and transverse stability'
PRIMARY='https://arxiv.org/html/2101.04006'
VARIANTS=('current_disk__spatial_max','current_disk__regional_mean',
          'current_radial_line__spatial_max','current_radial_line__regional_mean')
VARIANT_NAMES=('Disk, maxima','Disk, means','Radial line, maxima','Radial line, means')
CASE_NAMES=('+I','+G','−I','−G')
CASE_IDS=tuple('max_'+s+'contrast_change_small_'+j for s in ('','negative_') for j in ('05','15'))
KL=('max_negative_contrast_change_small_16','max_negative_contrast_change_small_17')
RADIUS=Q(1,2**34)
PINS={
    GLOBAL/'results.json':'8ddecf0622e57e8dd8d3e5b113c3261ad1c4981d5f6ac231585cf8adc69af2b8',
    GLOBAL/'root_observed_exit.json':'8352f1d4fda4f0d198d802fd79ba3a2c9f3f8278497f62ed6e017b840c6df2b6',
    OB/'actual/results.json':'4c5a601062c307d281a49ae4a0febfd4619ac9064723b1978c2d24e87b6e340b',
    OB/'actual/os_exit_observation.json':'dea0cebe28e89a8afdc36080bddea257782dd9d8001e6921a605c982202ea236',
    PREV/'build_receipt.json':'2f3704013bc9f26c10ba327c88051cc7c5bcf3566528ce8b407a7c5fc716d0d5',
    PREV/'root_readback.json':'e214452153309746dfe19712a5bfdb14798f11ae8349d6f1049559981c1ea8e1',
    PREV/'build.py':'d2c11b9948df0b4aef18cc5a9a6a84198cc438efb13d6c8f037cbdde9e24baff',
    CACHED:'a06549f048abba4eb2b4d1c1fb398f06e8f524196c92a673bd0493c527c8f40d',
    TRANS/'independent_actual_readback.json':'bf20e708826c08ed674c183e8c7364048748340501645885a12aa035f68d794e',
    TRANS/'root_actual_acceptance.json':'ed65fd8cb2e4d052ab465581c5567a488f18423396f6c79e7f8ddd4b3b361f46',
    THEORY/'actual/results.json':'331ee98456e7a8ee83c5e9e085dd88a1a764bd4b2b216901dd56f32ae060b04f',
    LOAD/'actual/receipt.json':'595988a5a2b15b237ad13ec575290d915046dc2d9ae01c8243547923bac27ba9',
    LOAD/'actual/os_exit_observation.json':'d33e301a76decc888c803227a88058fab7aa1714c373536f49acd36750c340a0',
    LOAD/'actual/root_acceptance.json':'c3b02fbc1e121bbc8febe7a84eb787274aee44b29f306ba2a874454f7f92a5a4',
    LOAD/'segment_actual/receipt.json':'97953457778f39011d96cc071070d2d5d342eb362f7e8e904d28d7947fd6e214',
    LOAD/'root_segment_actual_readback/receipt.json':'3e3cbe88a83aa3e8723b4c2ebfbf28809944ee5734850f0206c3b97269b74513',
    HELPER:'50524d01db5003b0839e873efb37d41bb0552ad0f0c162f73bd668ad0e58d31b',
    CONVERT:'90148d53d74d8ea5cc8a665e6a3660abf2db72487a1fa59586caba98f9cc38e2',
    ASSETS/'katex.min.css':'5a0ff74f11887f7260dc16b216906f51352b84ca29f21abb9a14b62a2ca6eb15',
    ASSETS/'katex.min.js':'e6bfe5deebd4c7ccd272055bab63bd3ab2c73b907b6e6a22d352740a81381fd4',
    FAV:'961e6e85ec16258652bf3c3f322d336fb598b90bccc49821047babfef3925530',
}
BG, GRID, DIM, INK, CYAN, PERI, CORAL, GOLD = '#122034','#3a516a','#c1d1e4','#f1f6ff','#80d9e8','#9caef1','#ffb6a0','#f6d16f'


def require(flag,message):
    if not flag:
        raise ValueError(message)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path,pins,h=None):
    path=path.resolve()
    require(path.is_relative_to(SITE),'foreign evidence path')
    expected=pins[path] if h is None else h
    require(sha(path)==expected==pins.setdefault(path,expected),f'bad proof digest for {path.name}')
    return json.loads(path.read_bytes())


def link(path):
    return Path(os.path.relpath(path,HERE)).as_posix()


def qpair(raw):
    require(type(raw) is list and len(raw)==2 and all(type(x) is str for x in raw),'rational endpoints expected')
    a,b=map(Q,raw)
    require(a<=b,'ordered endpoints required')
    return a,b


def zero(x):
    return type(x) is int and x==0


def membership(v,allowance):
    lo,hi=qpair(v['exact_spatial_root_observer_interval_um_exact'])
    observed=Q(v['observed_um_exact'])
    a,b=lo-observed,hi-observed
    res=(Q(0) if a<=0<=b else min(abs(a),abs(b)),max(abs(a),abs(b)))
    state='IN' if res[1]<=allowance else ('OUT' if res[0]>allowance else 'UNRESOLVED')
    require(res==qpair(v['absolute_residual_um_exact']) and state==v['standing'],'independent observer classification mismatch')
    return lo,hi,*res,state


def root_observer(pins,args):
    batch=read(GLOBAL/'results.json',pins)
    host=read(GLOBAL/'root_observed_exit.json',pins)
    accepted=read(OB_INDEPENDENT/'root_acceptance.json',pins,args.observer_acceptance_sha256)
    verify=read(OB_INDEPENDENT/'receipt.json',pins,accepted['receipt_sha256'])
    osdoc=read(OB_INDEPENDENT/'os_exit_observation.json',pins,accepted['root_true_os_sha256'])
    producer=read(OB/'actual/results.json',pins)
    pos=read(OB/'actual/os_exit_observation.json',pins)
    require(host['status']=='PASS_ROOT_OBSERVED_LATER_WAVE_INDEPENDENT_GLOBAL_EXIT'
            and host['root_scalar_and_672_leaf_checks_passed'] is True and zero(host['observed_exit_code'])
            and host['results_sha256']==pins[GLOBAL/'results.json'] and len(host['sha256'])==11
            and batch['status']=='COMPLETE_INDEPENDENT_LATER_WAVE_SIGNED_ENDPOINT_REPLAY'
            and batch['signed_hull_components']==672 and batch['full_five_gate_complete'] is False,
            'complete independent signed wave and genuinely observed exit required')
    for name,h in host['sha256'].items():
        p=(ROOT/name).resolve()
        require(p.is_relative_to(GLOBAL) and sha(p)==h==pins.setdefault(p,h),'native wave input link mismatch')
    require(accepted['status']=='PASS_ROOT_LATER_WAVE_OBSERVER_REHASH_AND_SCALARS'
            and type(accepted['all64_scalar_residual_statuses_rechecked']) is int
            and accepted['all64_scalar_residual_statuses_rechecked']==64 and accepted['full_five_gate_complete'] is False
            and verify['status']=='PASS_INDEPENDENT_ACTUAL_LATER_WAVE_SPATIAL_OBSERVER_READBACK'
            and osdoc['status']=='PASS_ROOT_OBSERVED_INDEPENDENT_LATER_OBSERVER_READBACK_EXIT'
            and osdoc['receipt_sha256']==pins[OB_INDEPENDENT/'receipt.json'] and zero(osdoc['observed_exit_code'])
            and pos['status']=='PASS_ROOT_OBSERVED_LATER_WAVE_OBSERVER_APPLICATION_EXIT'
            and pos['direct_root_tool_observation'] is True and zero(pos['observed_exit_code'])
            and pos['results_sha256']==pins[OB/'actual/results.json']
            and verify['producer_os_sha256']==osdoc['producer_os_sha256']==pins[OB/'actual/os_exit_observation.json']
            and verify['native_global_os_sha256']==osdoc['native_global_os_sha256']==pins[GLOBAL/'root_observed_exit.json']
            and verify['observer_producer_results_sha256']==pins[OB/'actual/results.json']
            and verify['all_64_exact_comparisons_recreated'] is True and verify['complete_calibration_or_source_fit_claimed'] is False,
            'observed independent Wave 03 observer acceptance required')
    terminal=read(OB_INDEPENDENT/'terminal.json',pins,osdoc['terminal_sha256'])
    started=read(OB_INDEPENDENT/'started.json',pins,osdoc['started_sha256'])
    require(zero(terminal['exit_code']) and terminal['receipt_sha256']==pins[OB_INDEPENDENT/'receipt.json']
            and all(type(d['pid']) is int and d['pid']==osdoc['pid']>0 for d in (verify,terminal,started))
            and all(d['wave_id']=='wave_03' for d in (batch,host,verify,producer,pos))
            and tuple(x['case_id'] for x in verify['cases'])==tuple(x['case_id'] for x in batch['cases'])==CASE_IDS,
            'same observers, wave, and successful terminal required')
    allowance=Q(verify['allowance_um_exact'])
    roots,obs=[],[]
    for leaf,o,label in zip(batch['cases'],verify['cases'],CASE_NAMES,strict=True):
        cid=leaf['case_id']
        native=read(GLOBAL/cid/'results.json',pins,leaf['case_results_sha256'])
        operand=read(GLOBAL/cid/'operands.json',pins,leaf['operands_sha256'])
        a,t=native['arithmetic'],native['observer_and_optional_minimum']
        eta,z=Q(a['eta_upper_exact']),Q(a['z_upper_exact'])
        e,sigma,xi=(Q(o[x]) for x in ('posterior_exact','archived_sigma_exact','original_xi_exact'))
        m=(1-z)/Q(operand['R_infinity_norm_exact'])
        eps=Q(67,2)*(sigma+e+xi)
        frac=1-z-eta/RADIUS
        require(native['status']=='PASS_INDEPENDENT_LATER_WAVE_EXACT_SPATIAL_ROOT'
                and leaf['exact_spatial_root_certified'] is True and leaf['exact_spatial_minimum_certified'] is True
                and native['root_asserted'] is True and t['strict_minimum_proved'] is True
                and native['target']==operand['target'] and o['root_proof_sha256']==leaf['case_results_sha256']
                and o['original_field_sha256']==native['target']['original_observer_field_sha256']
                and o['record_id'].endswith('_'+label[-1])
                and Q(native['radius_exact'])==RADIUS and 0<=z<1 and 0<=eta
                and eta+z*RADIUS==Q(a['selfmap_lhs_upper_exact'])<RADIUS
                and z==Q(a['ZQ_upper_exact'])+Q(a['b_H_upper_exact']) and a['strict_selfmap'] is True
                and e==eta/(1-z)==Q(a['posterior_error_upper_exact'])<RADIUS
                and m==Q(t['normalized_H_euclidean_lower_exact'])>0 and frac>0
                and sigma==Q(native['target']['archived_center_shift_exact'])
                and xi==Q(native['target']['original_addition_gap_exact'])
                and eps==Q(t['physical_nodal_error_um_upper_exact'])==Q(o['physical_nodal_error_um_exact']),
                'exact force radius, derivative, curvature and same original observer field required')
        roots.append({'case_id':cid,'label':label,'record_id':o['record_id'],'radius_exact':RADIUS,'eta_upper_exact':eta,
                      'eta_over_radius_exact':eta/RADIUS,'z_upper_exact':z,'selfmap_headroom_fraction_exact':frac,
                      'normalized_curvature_lower_exact':m,'posterior_exact':e,'archive_shift_exact':sigma,
                      'binary_addition_gap_exact':xi,'physical_nodal_error_um_exact':eps,'proof_sha256':leaf['case_results_sha256']})
        require(len(o['comparisons'])==len({(v['variant'],v['observable']) for v in o['comparisons']})==16,
                'complete new observer comparison map')
        for v in o['comparisons']:
            p0,p1,res0,res1,state=membership(v,allowance)
            obs.append({'case_id':cid,'label':label,'variant':v['variant'],'observable':v['observable'],
                        'prediction_lo_um_exact':p0,'prediction_hi_um_exact':p1,'observed_um_exact':v['observed_um_exact'],
                        'residual_lo_um_exact':res0,'residual_hi_um_exact':res1,'allowance_um_exact':allowance,
                        'standing':state,'strict_OUT_excess_um_exact':res0-allowance if state=='OUT' else ''})
    counts=Counter((o['variant'],o['standing']) for o in obs)
    require(len(obs)==64 and [counts[v,'IN'] for v in VARIANTS]==[16,14,16,14]
            and [counts[v,'OUT'] for v in VARIANTS]==[0,2,0,2] and counts.total()==64,
            'four negative alternative decisions and no unresolved comparisons required')
    fail=[x for x in obs if x['standing']=='OUT']
    require(all(x['label'] in ('+G','−G') and x['observable']=='apical_max_shortening_um'
                and x['variant'] in VARIANTS[1::2] for x in fail),'same alternative apical mean failures')
    prev,prev_root=read(PREV/'build_receipt.json',pins),read(PREV/'root_readback.json',pins)
    require(prev_root['status']=='PASS_ROOT_WAVE02_REPORT_SCIENTIFIC_VISUAL_READBACK'
            and prev_root['build_receipt_sha256']==pins[PREV/'build_receipt.json']
            and prev_root['disjoint12_root_coverage_checked_separately'] is True,
            'independently checked cumulative baseline needed')
    old=prev['exact_coverage']
    groups={k:old[k] for k in ('original_four','first_wave','wave_02')}
    groups['wave_03']=list(CASE_IDS)
    ids=[s for vals in groups.values() for s in vals]
    cached=read(CACHED,pins)
    cache_ids={x['case_id'] for x in cached['observer_transfer_terms']}
    require(len(ids)==len(set(ids))==16 and len(cache_ids)==cached['cached_fixed_parameter_roots']==36
            and set(ids)==set(accepted['disjoint_case_ids'])<=cache_ids
            and sum('_small_' in x for x in ids)==14 and sum('_large_' in x for x in ids)==2,
            'disjoint new coverage and subset of cached census required')
    cumulative={v:{s:old['counts_by_variant'][v][s]+counts[v,s] for s in ('IN','OUT')} for v in VARIANTS}
    require([cumulative[v]['IN'] for v in VARIANTS]==[56,51,54,50]
            and [cumulative[v]['OUT'] for v in VARIANTS]==[0,5,2,6],'same cumulative 224 observer comparisons')
    return roots,obs,counts,allowance,{'groups':groups,'root_count':16,'small_fields':14,'large_fields':2,
        'remaining_small_fields':20,'remaining_primary_inequalities':80,'observers':cumulative}


def transverse(pins):
    ra=read(TRANS/'root_actual_acceptance.json',pins)
    ind=read(TRANS/'independent_actual_readback.json',pins)
    original=read(THEORY/'actual/results.json',pins)
    require(ra['status']=='PASS_ROOT_OBSERVED_INDEPENDENT_TRANSVERSE_SCALAR_REPLAY'
            and zero(ra['observed_exit_code']) and ra['windows_checked']==ind['checked_windows']==original['evaluated_windows']==672
            and ra['independent_receipt_sha256']==pins[TRANS/'independent_actual_readback.json']
            and ra['conditional_natural_3D_energy_descent_verified'] is True
            and ra['physical_wavelength_or_continuum_equilibrium_claimed'] is False
            and ind['status']=='PASS_INDEPENDENT_ACTUAL_ALL672_TRANSVERSE_STRESS_WINDOWS'
            and ind['no_nesting_or_intersection_used_as_proof'] is True
            and tuple(x['case_id'] for x in ind['case_readback'])==tuple(x['case_id'] for x in original['cases'])==KL,
            'two accepted independent K/L transverse window sets required')
    for name,h in ind['verified_source_and_actual_sha256'].items():
        f=(ROOT/name).resolve()
        require(f.is_relative_to(ROOT) and sha(f)==h==pins.setdefault(f,h),'transverse proof source closure mismatch')
    flat,witnesses=[],[]
    for case,letter in zip(ind['case_readback'],('K','L'),strict=True):
        cid=case['case_id']
        seen=set()
        for w in case['all_windows']:
            a,b=qpair(w['window_t_exact'])
            r0,r1=qpair(w['reference_R_exact'])
            edge,face=w['edge'],w['face']
            state=w['independent_standing']
            tau=w['independent_tau_negative_upper_exact']
            weight=w['independent_minus_tau_r_R_squared_over_Js_lower_exact']
            require(face in ('apical','basal') and type(edge) is int and 0<=edge<12 and
                    Q(1,16)<=a<b<=Q(15,16) and b-a==Q(1,16) and 0<r0<r1
                    and (face,edge,a) not in seen and w['D_bounds_nondisjoint'] is True,
                    'each independent checked interior trace window is unique')
            seen.add((face,edge,a))
            if state=='STRICT_NEGATIVE_PATCH':
                require(face=='apical' and Q(tau)<0<Q(weight)
                        and qpair(w['independent_r_R_exact'])[0]>0,'negative and geometric strict witness required')
            row={'case_id':cid,'anchor':letter,'face':face,'edge_index':edge,'t_lo_exact':a,'t_hi_exact':b,
                 'reference_R_lo_exact':r0,'reference_R_hi_exact':r1,'standing':state,
                 'tau_strict_upper_exact':tau,'minus_tau_r_R_squared_over_Js_lower_exact':weight,
                 'r_R_lo_exact':qpair(w['independent_r_R_exact'])[0]}
            flat.append(row)
            if face=='apical' and edge==7 and (a,b)==(Q(5,8),Q(11,16)):
                require(state=='STRICT_NEGATIVE_PATCH','illustrative independent strict patch needed')
                witnesses.append(row)
        counts=Counter(x['standing'] for x in flat if x['case_id']==cid)
        require(len(seen)==336 and counts==case['independent_counts']==ra['counts'][cid]
                and counts=={'NONNEGATIVE_ON_WINDOW':276,'UNRESOLVED':2,'STRICT_NEGATIVE_PATCH':58}
                and case['angular_instability_direction_exists_in_stated_enriched_natural_3D_model'] is True,
                '58 complete negative windows at each root')
    require(len(witnesses)==2,'each exact example window available')
    return flat,witnesses


def geometry(pins,args):
    full=read(LOAD/'actual/receipt.json',pins)
    full_host=read(LOAD/'actual/root_acceptance.json',pins)
    before=read(LOAD/'actual/os_exit_observation.json',pins)
    trial=read(LOAD/'segment_actual/receipt.json',pins)
    seg=read(args.segment_acceptance,pins,args.segment_acceptance_sha256)
    independent=read(LOAD/'root_segment_actual_readback/receipt.json',pins)
    require(full['status']=='VERIFIED_UNRESOLVED_KL_WHOLE_ENVELOPE_BODY'
            and full_host['status']=='PASS_ROOT_REHASH_OF_DELEGATED_KL_UNRESOLVED_GEOMETRY_REPLAY'
            and full_host['independent_receipt_sha256']==pins[LOAD/'actual/receipt.json']==before['receipt_sha256']
            and full_host['delegated_os_observation_sha256']==pins[LOAD/'actual/os_exit_observation.json']
            and zero(before['observed_exit_code']) and full['root_or_loading_branch_proved'] is False,
            'accepted full-envelope geometry remains unresolved')
    require(independent['status']=='PASS_INDEPENDENT_FULL_SPAN_KL_FOUR_SLAB_BODY_LOG_COVERAGE'
            and independent['producer_receipt_sha256']==pins[LOAD/'segment_actual/receipt.json']
            and independent['first_complete_tested_dyadic_level']==trial['first_complete_tested_dyadic_level']==4
            and independent['surface_cache_force_H_or_loading_root_proved'] is False
            and independent['exact_original_closed_interval']==trial['exact_original_closed_interval']
            and Q(independent['constant_radius_exact'])==Q(trial['constant_corrector_radius_exact'])==Q(1,2**12)
            and trial['prior_full_scope_geometry_readback_os_sha256']==pins[LOAD/'actual/os_exit_observation.json']
            and trial['surface_and_cached_domains_certified'] is False and trial['root_or_loading_branch_proved'] is False,
            'accepted segmented body-only proof required')
    left,right=qpair(trial['exact_original_closed_interval'])
    b=full['body_summary']
    lo,hi=qpair(b['body_J_bounds_exact'])
    rows=[{'slabs':1,'slab':0,'I_lo_exact':left,'I_hi_exact':right,'J_lo_exact':lo,'J_hi_exact':hi,
           'monotonicity_lower_exact':'','rho_upper_exact':'','radial_margin_exact':Q(full['literal_radial_floor_exact']),
           'insufficient_J':len(b['witnesses_failed_positive_bound']),
           'insufficient_monotonicity':len(b['witnesses_failed_monotonicity_bound']),
           'insufficient_log_rho':'','body_and_log_guards_pass':False}]
    require(rows[0]['insufficient_J']==20 and rows[0]['insufficient_monotonicity']==34 and lo<0,
            'readback of exactly the broad-box guard failures')
    require([entry['slabs'] for entry in trial['levels']]==[2,4]
            and trial['first_complete_tested_dyadic_level']==4,'first complete tested dyadic level')
    for lv,iv in zip(trial['levels'],independent['levels'],strict=True):
        n=lv['slabs']
        require(len(lv['checks'])==len(iv['rows'])==iv['slabs']==n,'whole exact dose coverage in each tested level')
        for i,(leaf,replay) in enumerate(zip(lv['checks'],iv['rows'],strict=True)):
            interval=qpair(leaf['exact_dose_subinterval'])
            jl,jh=qpair(leaf['J_enclosure_exact'])
            passed=leaf['passes_body_J_monotonicity_log_rho_and_radial_bound']
            independent_body=read(LOAD/'root_segment_actual_readback'/f'body_{n}_{i}.json',pins,replay['body_sha256'])
            require(type(leaf['slot']) is int and leaf['slot']==i and passed is (n==4)
                    and type(replay['slot']) is int and replay['slot']==i and replay['passes'] is passed
                    and replay['J_bound']==independent_body['J_bound']
                    and (jl,jh)==(Q(replay['J_bound']['lower_exact']),Q(replay['J_bound']['upper_exact']))
                    and len(independent_body['triangles'])==144
                    and interval==(left+(right-left)*Q(i,n),left+(right-left)*Q(i+1,n)),
                    'closed exact dyadic subinterval and same result required')
            rows.append({'slabs':n,'slab':i,'I_lo_exact':interval[0],'I_hi_exact':interval[1],'J_lo_exact':jl,'J_hi_exact':jh,
                         'monotonicity_lower_exact':leaf['minimum_monotonicity_exact'] or '',
                         'rho_upper_exact':leaf['maximum_bulk_log_rho_exact'] or '',
                         'radial_margin_exact':Q(leaf['radial_nodal_margin_exact']),
                         'insufficient_J':len(leaf['insufficient_J_slots']),
                         'insufficient_monotonicity':len(leaf['insufficient_monotonicity_slots']),
                         'insufficient_log_rho':len(leaf['insufficient_bulk_log_rho_slots']),
                         'body_and_log_guards_pass':passed})
    return rows


def write(path,text):
    with path.open('x',encoding='utf-8',newline='\n') as out:
        out.write(text)


def csv_write(name,rows):
    with (HERE/name).open('x',encoding='utf-8',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def panel(size):
    fig,ax=plt.subplots(figsize=size,facecolor=BG)
    ax.set_facecolor(BG)
    for spine in ax.spines.values():
        spine.set_color(GRID)
    ax.tick_params(colors=DIM,labelsize=11)
    ax.xaxis.label.set_color(DIM)
    ax.yaxis.label.set_color(DIM)
    return fig,ax


def finish(fig,name):
    with (HERE/name).open('xb') as out:
        fig.savefig(out,dpi=178,format='png',facecolor=fig.get_facecolor())
    plt.close(fig)


def charts(roots,windows,slabs):
    fig,ax=panel((8.1,3.25))
    ax.set_xlim(0,1)
    ax.set_axisbelow(True)
    ax.xaxis.grid(color=GRID)
    for i,row in enumerate(roots):
        a,b=(float(row[k]) for k in ('eta_over_radius_exact','z_upper_exact'))
        ax.barh(i,a,height=.56,color=PERI)
        ax.barh(i,b,left=a,height=.56,color=CYAN)
        ax.text(a+b+.017,i,f'{(1-a-b)*100:.1f}% free',va='center',color=INK,fontsize=10)
    ax.set(yticks=range(4),yticklabels=CASE_NAMES,xlabel='Fraction of unit self-map budget')
    ax.invert_yaxis()
    ax.legend([Patch(color=PERI),Patch(color=CYAN)],['Center force / radius','Whole-cube derivative'],
              frameon=False,labelcolor=INK,loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2,fontsize=10)
    fig.subplots_adjust(left=.11,right=.97,top=.84,bottom=.20)
    finish(fig,'selfmap.png')

    fig,ax=panel((8.6,3.38))
    color={'NONNEGATIVE_ON_WINDOW':CYAN,'STRICT_NEGATIVE_PATCH':CORAL,'UNRESOLVED':GOLD}
    labels=('K apical','K basal','L apical','L basal')
    for w in windows:
        y=labels.index(w['anchor']+' '+w['face'])
        start=w['edge_index']+float(w['t_lo_exact'])
        width=float(w['t_hi_exact']-w['t_lo_exact'])
        ax.add_patch(Rectangle((start,y-.23),width,.46,color=color[w['standing']],lw=0))
    ax.set(xlim=(0,12),ylim=(3.61,-.61),xticks=range(13),yticks=range(4),yticklabels=labels,
           xlabel='Reference trace edge index + within-edge coordinate')
    ax.set_axisbelow(True)
    ax.xaxis.grid(color=GRID,alpha=.75)
    ax.legend([Patch(color=CYAN),Patch(color=CORAL),Patch(color=GOLD)],['Nonnegative','Strict negative','Unresolved'],
              frameon=False,labelcolor=INK,loc='lower center',bbox_to_anchor=(.5,1.01),ncol=3,fontsize=10)
    fig.subplots_adjust(left=.16,right=.96,top=.84,bottom=.19)
    finish(fig,'transverse_windows.png')

    fig,ax=panel((8.1,3.85))
    x=list(range(len(slabs)))
    for i,row in enumerate(slabs):
        ax.bar(i,float(row['J_lo_exact']),color=CYAN if row['body_and_log_guards_pass'] else CORAL,width=.65)
    ax.axhline(0,color=INK,lw=1)
    ax.set_axisbelow(True)
    ax.yaxis.grid(color=GRID)
    ax.set(xticks=x,xticklabels=['Whole','Half 1','Half 2','¼ 1','¼ 2','¼ 3','¼ 4'],ylabel='Certified lower bound for J')
    ax.legend([Patch(color=CYAN),Patch(color=CORAL)],['Body + log guards passed','One or more guards unresolved'],
              loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2,frameon=False,labelcolor=INK,fontsize=9.5)
    fig.subplots_adjust(left=.15,right=.97,top=.84,bottom=.15)
    finish(fig,'slab_guards.png')


def decimal_outward(q,places,up):
    q=Q(q)
    units=10**places
    n=-((-q.numerator*units)//q.denominator) if up else (q.numerator*units)//q.denominator
    return f'{"−" if n<0 else ""}{abs(n)//units}.{abs(n)%units:0{places}d}'


def markdown(roots,counts,allowance,coverage,witnesses,slabs):
    root_table='\n'.join(f"| {x['label']} | {float(x['z_upper_exact']):.5f} | "
                         f"{100*float(x['selfmap_headroom_fraction_exact']):.2f}% | "
                         f"{float(x['normalized_curvature_lower_exact']):.8f} | "
                         f"{float(x['physical_nodal_error_um_exact']):.4g} |" for x in roots)
    all_counts=coverage['observers']
    obs_table='\n'.join(f"| {s} | {counts[v,'IN']} | {counts[v,'OUT']} | {all_counts[v]['IN']} | {all_counts[v]['OUT']} |"
                       for v,s in zip(VARIANTS,VARIANT_NAMES,strict=True))
    example_table='\n'.join(f"| {w['anchor']} | {decimal_outward(w['tau_strict_upper_exact'],12,True)} | "
                            f"{decimal_outward(w['minus_tau_r_R_squared_over_Js_lower_exact'],12,False)} |"
                            for w in witnesses)
    min_j=min(x['J_lo_exact'] for x in slabs if x['slabs']==4)
    min_spd=min(Q(x['monotonicity_lower_exact']) for x in slabs if x['slabs']==4)
    max_rho=max(Q(x['rho_upper_exact']) for x in slabs if x['slabs']==4)
    return fr"""# {TITLE}

Four new exact-spatial root and minimum proofs bring the total to **16 of 36** frozen fixed-input outcomes. A separate transverse test gives a stricter limit on interpretation. At the previously certified negative-material K and L fields, the natural three-dimensional extension of the same energy admits lower-energy angular perturbations. Both facts hold because the finite axisymmetric certificate tests only zero-swirl meridional degrees of freedom. [Wave 03 root proofs]({link(GLOBAL/'results.json')}) · [Independent K/L transverse result]({link(TRANS/'root_actual_acceptance.json')})

## Transverse stability at K and L

The surface density on either disk is \(K_s(J_s-1)^2/2+T_s(J_s-1)\). Its derivative with respect to the area ratio \(J_s\) gives the total surface stress coefficient,

\[
\tau_s=K_s(J_s-1)+T_s.
\]

At a compressed patch, the area-elastic term can make \(\tau_s\) negative. The exact P2 traces and certified posterior nodal enclosures give **58 strictly negative, 276 nonnegative, and 2 unresolved windows per anchor**. The independent readback verified all 672 windows, and every strict negative window lies on the apical face. The same apical witness on edge 7, at local coordinate \(t\in[5/8,11/16]\), has the following outward-rounded bounds. Here \(r_R\) is the radial derivative with respect to reference radius, and \(-\tau_a r_R^2/J_s\) is the positive weight associated with the negative angular surface term.

| Anchor | Upper bound on \(\tau_a\) | Lower bound on \(-\tau_a r_R^2/J_s\) |
|:--|--:|--:|
{example_table}

![Each of the K and L apical traces contains a band of red strictly negative checked windows across several interior edges. The basal checked windows are nonnegative. Two apical windows at each anchor remain unresolved.](transverse_windows.png)

Each plotted edge has equal width. Blank pieces are outside the predeclared interior windows. The CSV supplies exact reference radii and the individual bounds used to classify each window. The result rests on full windows away from the axis and the clamped wall. [All 672 exact window records](transverse_windows.csv) · [Independent scalar readback]({link(TRANS/'independent_actual_readback.json')})

For the declared natural three-dimensional extension with no surface bending, take a smooth annular axial displacement \(\psi(R,Z)\cos(m\theta)e_z\) supported near one strict apical patch. The first variation in that angular harmonic averages to zero. The second variation has form

\[
D^2E[v_m,v_m]=\frac{{\pi}}{{S}}(C_0+m^2C_2).
\]

The negative surface part of \(C_2\) stays fixed when the support extends a shorter distance into the bulk. Its positive bulk part decreases with that support depth. After choosing a depth with \(C_2<0\), \(C_0\) is finite and sufficiently large integer \(m\) gives a negative second variation. Small amplitudes preserve the wall clamp and admissible geometry. This proves energy descent in the enlarged variation space at the revolved K and L fields. The axisymmetric exact finite-element minima remain valid in their own space.

The construction probes arbitrarily short wavelengths of the bending-free extension. The thin layer and angular wavelength used in the proof may fall below cell-scale validity. A physically resolved transverse prediction calls for a declared scale cutoff or bending law and tests of finite angular modes. A three-dimensional continuum equilibrium, the first unstable mode, a physical wavelength and tissue behavior remain unestablished. Negative-surface-stress instability and its short-wave regularization by bending are established in prior work. [Derivation and scope]({link(THEORY/'READBACK.md')}) · [Binysh, Wilks and Souslov, *Science Advances*, 2022]({PRIMARY})

## Wave 03 fixed-input proofs

The four new exact-spatial cases are small records I and G at both frozen material choices. The signs label those choices. For each 611-dimensional cube of radius \(r=2^{{-34}}\), the signed preconditioned center residual has upper bound \(\eta\) and the whole-cube derivative has contraction bound \(z\). The strict test is \(\eta/r+z<1\). The positive symmetric part of the same binary \(R\) gives the normalized Hessian lower bound \((1-z)/\|R\|_\infty\). The observer error in the table includes the root posterior, any archived center shift and the original binary addition gap.

| Case | \(z\) | Headroom | Curvature floor | Nodal error (µm) |
|:--|--:|--:|--:|--:|
{root_table}

![Stacked center-force and derivative terms for +I, +G, −I and −G lie below the unit self-map budget. The free percentage equals the certified margin divided by the cube radius.](selfmap.png)

The new four case IDs are distinct from the previous twelve. The exact-spatial total contains fourteen small and two large fields, all within the already certified 36 cached-quadrature outcomes. The plotted marks are rounded displays. The exact signed certificates govern the classifications. [Exact bounds CSV](exact_bounds.csv) · [Observed independent root process]({link(GLOBAL/'root_observed_exit.json')}) · [Wave 02 checkpoint]({link(PREV/'REPORT.html')})

## Observer comparisons at the fixed allowance

The frozen residual allowance is approximately **{float(allowance):.9f} µm**. It is a model discrepancy threshold. An IN comparison has an entire absolute residual interval at or below that allowance. An OUT comparison has its entire interval above it. All four new OUT rows concern the apical shortening channel at G when computed using regional means. The primary disk-max comparisons are all IN.

| Convention | Wave 03 IN | OUT | Cumulative IN | OUT |
|:--|--:|--:|--:|--:|
{obs_table}

Cumulative exact coverage gives **56 primary IN** and **155 alternative IN / 13 OUT**. Twenty small fields and 80 primary inequalities remain for the numerically exact two-candidate comparison across all 34 small fields. Source joining, experimental discrepancy bounds, and a physically matched large shape outcome remain separate open questions. [All 64 Wave 03 exact comparisons](exact_observer_comparisons.csv) · [Accepted independent observer replay]({link(OB_INDEPENDENT/'root_acceptance.json')})

## K–L continuation geometry

The wide rectangular hull of the exact secant predictor failed the sufficient body-enclosure test. Twenty element bounds could not prove positive \(J\), and 34 could not prove positive-definite symmetric meridional gradients. The failure concerns this enclosure. It supplies no evidence of an inverted state or a broken branch. Retaining the same predictor and corrector radius, the first complete tested dyadic partition uses four closed subintervals. Each subinterval passes the body embedding, bulk log expansion and literal nodal radial bound checks. Across those four enclosures the worst lower bound on \(J\) is at least {decimal_outward(min_j,6,False)}, the worst lower bound for strong monotonicity is at least {decimal_outward(min_spd,6,False)}, and the bulk log expansion bound is at most {decimal_outward(max_rho,6,True)} on its dimensionless convergence parameter.

![The single K–L envelope gives a negative lower bound on J. The two half-interval boxes give positive determinant lower bounds but fail other geometry guards. The four quarter-interval boxes give positive determinant lower bounds and pass the tested body and bulk-log guards.](slab_guards.png)

Surface and cached analytic domains, uniform mechanical bounds on the slabs, and same-root junction proofs remain required. The accepted geometry result narrows the continuation calculation. A strict-minimum branch between these endpoint fields, zero-to-full loading, and biological dynamics remain open. [Broad-box independent readback]({link(LOAD/'actual/receipt.json')}) · [Exact geometry bounds CSV](exact_slab_geometry.csv) · [Independent dyadic replay]({link(LOAD/'root_segment_actual_readback/receipt.json')})
"""


def render(md):
    spec=spec_from_file_location('_checkpoint_wave03_static_math',HELPER)
    helper=module_from_spec(spec)
    spec.loader.exec_module(helper)
    sys.path.insert(0,str(CONVERT.parent))
    import build_report
    require(Path(build_report.__file__).resolve()==CONVERT,'wrong pinned converter')
    body=build_report.convert(md)
    pattern=re.compile(r'<(span|div) class="math( display)?" data-tex="([^"]+)"></\1>')
    matches=list(pattern.finditer(body))
    tasks=[{'tex':unescape(m[3]),'display':m[1]=='div'} for m in matches]
    displays=len(re.findall(r'\\\[(.*?)\\\]',md,re.S))
    inlines=len(re.findall(r'\\\((.*?)\\\)',md,re.S))
    require(len(tasks)==displays+inlines>0 and sum(t['display'] for t in tasks)==displays,
            'every source inline and display expression must reach the renderer')
    compiled=subprocess.run(['node','-e',helper.MATH_SNIPPET,str(ASSETS/'katex.min.js')],
                            input=json.dumps(tasks),encoding='utf-8',capture_output=True,check=True,cwd=HERE)
    rendered=json.loads(compiled.stdout)
    require(len(rendered)==len(matches) and all('class="katex"' in x for x in rendered),
            'missing strict compiled KaTeX expression')
    stream=iter(rendered)
    body=pattern.sub(lambda m:f'<{m[1]} class="math{m[2] or ""}">{next(stream)}</{m[1]}>',body)
    count=len(matches)
    oldspec=spec_from_file_location('_checkpoint_wave02_style',PREV/'build.py')
    old=module_from_spec(oldspec)
    oldspec.loader.exec_module(old)
    html=('<!doctype html>\n<html lang="en"><head><meta charset="utf-8"><meta name="viewport" '
          'content="width=device-width,initial-scale=1">'
          f'<title>{escape(TITLE)}</title><link rel="icon" type="image/png" href="{link(FAV)}">'
          f'<link rel="stylesheet" href="{link(ASSETS/"katex.min.css")}"><style>{old.STYLE}</style></head>'
          '<body><main><nav aria-label="Report resources"><a href="REPORT.md">Markdown</a>'
          '<a href="exact_bounds.csv">Exact roots</a><a href="exact_observer_comparisons.csv">Observer comparisons</a>'
          '<a href="transverse_windows.csv">Stress windows</a><a href="exact_slab_geometry.csv">Slab geometry</a>'
          '</nav>'+body+'</main></body></html>\n')
    parser=helper.LocalLinks()
    parser.feed(html)
    local={}
    for v in parser.paths:
        if v==PRIMARY:
            continue
        p=(HERE/v).resolve()
        require(p.is_relative_to(SITE) and p.is_file(),f'broken local link {v}')
        local[v]=p.relative_to(SITE).as_posix()
    return html,count,local


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--observer-acceptance-sha256',required=True)
    parser.add_argument('--segment-acceptance',type=Path,required=True)
    parser.add_argument('--segment-acceptance-sha256',required=True)
    p=parser.parse_args()
    pins=PINS.copy()
    require(all(sha(path)==h for path,h in pins.items()),'frozen proof or render source drift')
    roots,obs,counts,delta,coverage=root_observer(pins,p)
    windows,witnesses=transverse(pins)
    slabs=geometry(pins,p)
    rows={'exact_bounds.csv':roots,'exact_observer_comparisons.csv':obs,'transverse_windows.csv':windows,
          'exact_slab_geometry.csv':slabs}
    for name,data in rows.items():
        csv_write(name,data)
    charts(roots,windows,slabs)
    md=markdown(roots,counts,delta,coverage,witnesses,slabs)
    write(HERE/'REPORT.md',md)
    html,math,local=render(md)
    write(HERE/'REPORT.html',html)
    outputs=('REPORT.md','REPORT.html',*rows,'selfmap.png','transverse_windows.png','slab_guards.png')
    receipt={'status':'PASS_BUILD_FROM_OBSERVED_WAVE03_ROOT_OBSERVER_TRANSVERSE_AND_GEOMETRY_READBACKS',
             'source_sha256':sha(Path(__file__)), 'input_sha256':{f.relative_to(SITE).as_posix():h for f,h in pins.items()},
             'file_sha256':{out:sha(HERE/out) for out in outputs},'local_links':local,'external_primary_source':PRIMARY,
             'rendered_math_count':math,'exact_coverage':coverage,'csv_rows':{k:len(v) for k,v in rows.items()},
             'exact_plot_values':{'root_terms':[[r['eta_over_radius_exact'],r['z_upper_exact']] for r in roots],
                                  'slab_J_lower_bounds':[r['J_lo_exact'] for r in slabs]},
             'axisymmetric_minima_preserved':True,'conditional_transverse_descent_at_K_and_L':True,
             'continued_KL_path_or_continuum_or_physical_wavelength_established':False,'full_five_gate_complete':False}
    write(HERE/'build_receipt.json',json.dumps(receipt,indent=2,default=str)+'\n')
    print(json.dumps({'receipt_sha256':sha(HERE/'build_receipt.json'),'math':math,'rows':receipt['csv_rows']}))


if __name__=='__main__':
    main()
