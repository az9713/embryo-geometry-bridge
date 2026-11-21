"""All-68-entry joint gain candidate and direct nonlinear replay.

One declared mechanical scenario; four observation alternatives retained.
The exact surrogate optimum is NOT an exact nonlinear or continuum optimum.
No held-out response file is read. Every failed nonlinear attempt is retained.
"""
from dataclasses import asdict
from fractions import Fraction as Q
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import time
import numpy as np
import torch
from joint_calibration_contract import OBSERVABLES, load_observations, evaluate
from piecewise_gain_minimax import fit
from quadratic_endpoint_observer import observe_state
from sparse_quadratic import AxisymmetricModel, Parameters

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'runs/joint_gain_calibration'
OBSERVER_SHA = '075dac2148ed26a3bfe83fadf1b97602e4b08536f2be38cb95c2054776c0b926'
VARIANTS = [(w,c) for w in ('current_disk','current_radial_line')
            for c in ('spatial_max','regional_mean')]
TABLE_FILES = {f:ROOT/'runs/activated_quadrature'/f'lambda_5000.0_strength_{float(f)}'/'reequilibrated.json'
               for f in (1,5,10)}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path,data):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def hashes():
    names = ('joint_gain_calibration.py','joint_calibration_contract.py','piecewise_gain_minimax.py',
             'quadratic_endpoint_observer.py','sparse_quadratic.py','quadratic_quadrature_study.py',
             'robust_quadratic.py','axisymmetric_quadratic.py','axisymmetric_mixed.py',
             'axisymmetric_prestress.py','axisymmetric_mechanics.py')
    return {n:sha(Path(__file__).with_name(n)) for n in names}


def inputs():
    with (ROOT/'data/derived/calibration_small_18um.csv').open(newline='',encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    return rows,load_observations()


def endpoints(state):
    out = observe_state(state,weighting='current_disk',component_aggregation='spatial_max',tolerance_um='0.001')
    if not out['maximum_shortening_certificate']['requested_width_achieved']:
        raise RuntimeError('endpoint reporting interval not resolved')
    all_values = {}
    for weight,component in VARIANTS:
        means = out['mean_alternatives'][weight]
        a,b = ((out['spatial_max_components']['apical_downward'],out['spatial_max_components']['basal_upward'])
               if component == 'spatial_max' else
               (means['apical_regional_mean_downward'],means['basal_regional_mean_upward']))
        values = (out['four_endpoints']['maximum_Eulerian_height_shortening'],
                  means['height_regional_mean_shortening'],a,b)
        all_values[weight+'__'+component] = {key:[str(v.lo),str(v.hi)] for key,v in zip(OBSERVABLES,values)}
    return all_values


def prepare():
    OUT.mkdir(parents=True,exist_ok=True)
    protocol_path = OUT/'protocol.json'
    if protocol_path.exists():
        raise RuntimeError('preparation already recorded; retain it rather than overwrite or duplicate')
    if sha(Path(__file__).with_name('quadratic_endpoint_observer.py')) != OBSERVER_SHA:
        raise RuntimeError('observer changed since component verification')
    rows,observed = inputs(); doses = [Q(r['activation_magnitude']) for r in rows]
    states = {f:json.loads(p.read_text()) for f,p in TABLE_FILES.items()}
    p = states[1]['parameters']
    if any(s['parameters'] != p or not s['converged'] for s in states.values()):
        raise RuntimeError('matched verified response-table fields required')
    protocol = {'scope':'Conditional joint gain calibration in ONE fixed mechanical scenario; not full material/geometry identified-set inference',
                'physical_parameters':p, 'length_unit_um':'33.5', 'source_hashes':hashes(),
                'calibration_csv_sha256':sha(ROOT/'data/derived/calibration_small_18um.csv'),
                'input_state_hashes':{str(f):sha(path) for f,path in TABLE_FILES.items()},
                'strength_nodes':['0','1','5','10'],
                'gain_bounds_exact':['0',str(Q(10)/max(doses))],
                'gain_bound_reason':'No extrapolation of available response-table strengths; not a physiological bound. Boundary selection triggers an explicit limitation.',
                'observation_variants':[{'weighting':w,'component_aggregation':c} for w,c in VARIANTS],
                'common_assumptions':['static endpoint used as a peak proxy, not a validated temporal maximum',
                     'flat initial reference with no growth and reference component heights H/0',
                     'fixed current-radius ROI equal to the initial patch radius',
                     'canonical author-provided source-column pairing',
                     'uniform micrometre minimax objective across 68 endpoints, not a measurement-error model'],
                'objective':'min_g max_i,k |piecewise-linear midpoint-table prediction(g I_i,k) - observed_i,k|',
                'certification_scope':'Exact optimum for the DEFINED midpoint response-table surrogate only; nonlinear interpolation/continuum/experimental-peak errors are not enclosed',
                'direct_verification':'Solve all distinct g_variant*I_i in ascending order from the last verified state, using unchanged sparse high-Q mechanics. Retain failures and exact static-curve endpoint intervals.',
                'direct_status':'No claim of nonlinear optimality, full compatible set or biological validation from this candidate replay',
                'held_out_response_files_read':False}
    save(protocol_path,protocol)  # Written before response extraction and fitting.
    model = AxisymmetricModel(Parameters(**p))
    zero = model.advance(0.,output_directory=OUT/'zero')
    if not zero['reached_target']:
        raise RuntimeError('zero reference did not verify')
    save(OUT/'zero/final.json',zero); states[0] = zero
    table = []
    for strength in (0,1,5,10):
        values = endpoints(states[strength])
        table.append({'strength_exact':str(strength),'variant_endpoint_intervals':values})
        print(json.dumps({'stage':'response_table','strength':strength}),flush=True)
    save(OUT/'response_table.json',{'protocol_sha256':sha(protocol_path),'rows':table})
    variants = []
    for weight,component in VARIANTS:
        key = weight+'__'+component
        response = [[(Q(row['variant_endpoint_intervals'][key][k][0])+Q(row['variant_endpoint_intervals'][key][k][1]))/2
                     for k in OBSERVABLES] for row in table]
        observed_matrix = [[observed[r['record_id']][k] for k in OBSERVABLES] for r in rows]
        fitted = fit([0,1,5,10],response,doses,observed_matrix,protocol['gain_bounds_exact'])
        variants.append({'variant':key,'weighting':weight,'component_aggregation':component,'surrogate_fit':fitted})
        print(json.dumps({'stage':'exact_surrogate_fit','variant':key,'gain':fitted['gain_display'],
                          'loss_um':fitted['loss_display'],'cells':len(fitted['certified_cells']),
                          'boundary':fitted['at_declared_gain_boundary']}),flush=True)
        save(OUT/'surrogate_candidates.json',{'status':'COMPLETE' if len(variants)==4 else 'IN_PROGRESS',
             'protocol_sha256':sha(protocol_path),'response_table_sha256':sha(OUT/'response_table.json'),
             'variants':variants,'scope':protocol['certification_scope']})
    if protocol['source_hashes'] != hashes():
        raise RuntimeError('sources changed during preparation')


def direct():
    protocol = json.loads((OUT/'protocol.json').read_text())
    candidates = json.loads((OUT/'surrogate_candidates.json').read_text())
    if candidates['status'] != 'COMPLETE' or protocol['source_hashes'] != hashes():
        raise RuntimeError('complete unchanged-source candidate preparation required')
    if candidates['protocol_sha256'] != sha(OUT/'protocol.json') or candidates['response_table_sha256'] != sha(OUT/'response_table.json'):
        raise RuntimeError('candidate input hashes changed')
    if protocol['calibration_csv_sha256'] != sha(ROOT/'data/derived/calibration_small_18um.csv'):
        raise RuntimeError('small-region observations changed')
    job = OUT/'direct_job.json'
    if job.exists():
        raise RuntimeError('direct job record exists; inspect live/terminal state, do not duplicate')
    save(job,{'status':'RUNNING','pid':os.getpid(),'source_hashes':hashes()})
    rows,observed = inputs(); requests = {}
    for variant in candidates['variants']:
        g = Q(variant['surrogate_fit']['gain_exact'])
        for row in rows:
            strength = g*Q(row['activation_magnitude'])
            requests.setdefault(strength,[]).append((variant['variant'],row['record_id']))
    model = AxisymmetricModel(Parameters(**protocol['physical_parameters']))
    previous = json.loads((OUT/'zero/final.json').read_text()); predictions = {v['variant']:{} for v in candidates['variants']}
    records = []; started = time.time()
    try:
        for index,(strength,uses) in enumerate(sorted(requests.items())):
            folder = OUT/'direct'/f'{index:03d}'
            result = model.advance(float(strength),previous,output_directory=folder)
            save(folder/'final.json',result)
            record = {'index':index,'requested_strength_exact':str(strength),
                      'actual_solver_strength':result['activation_strength'],
                      'uses':[{'variant':v,'record_id':rid} for v,rid in uses],
                      'file':(folder/'final.json').relative_to(ROOT).as_posix(),'sha256':sha(folder/'final.json'),
                      'reached':result['reached_target'],'force_linf':result['scaled_free_force_linf'],
                      'attempts':len(result['robust_continuation']['attempts']),
                      'rejected_attempts':sum(not r['accepted'] for r in result['robust_continuation']['attempts'])}
            records.append(record)
            if not result['reached_target']:
                save(OUT/'direct_progress.json',{'status':'STOPPED_TARGET_UNREACHED','records':records,'predictions':predictions})
                raise RuntimeError('direct target not reached; last verified state and failed trials retained')
            values = endpoints(result)
            for variant,rid in uses:
                predictions[variant][rid] = values[variant]
            previous = result
            save(OUT/'direct_progress.json',{'status':'COMPLETE' if index+1 == len(requests) else 'RUNNING',
                 'records':records,'predictions':predictions,'total_targets':len(requests),'elapsed_seconds':time.time()-started})
            print(json.dumps({'stage':'direct','index':index,'of':len(requests),'strength':float(strength),
                              'force':record['force_linf'],'attempts':record['attempts']}),flush=True)
        summaries = []
        for variant in candidates['variants']:
            key = variant['variant']; pred = predictions[key]
            residuals = evaluate(observed,pred,{k:0 for k in OBSERVABLES})
            lower = max(Q(v) for v in residuals['discrepancy_lower_requirements_exact'].values())
            upper = max(Q(v) for v in residuals['discrepancy_sufficient_requirements_exact'].values())
            summaries.append({'variant':key,'surrogate_gain_exact':variant['surrogate_fit']['gain_exact'],
                 'surrogate_loss_exact':variant['surrogate_fit']['loss_exact'],
                 'actual_static_discrepancy_requirement_um_exact':[str(lower),str(upper)],
                 'actual_static_discrepancy_requirement_um_display':[float(lower),float(upper)],
                 'all_68_residuals':residuals,
                 'interpretation':'Discrepancy requirement at this direct candidate, not estimated noise or an optimum over nonlinear mechanical parameters'})
        if protocol['source_hashes'] != hashes():
            raise RuntimeError('sources changed during direct verification')
        save(OUT/'direct_results.json',{'status':'COMPLETE_CONDITIONAL_GAIN_CANDIDATE_REPLAY',
             'protocol_sha256':sha(OUT/'protocol.json'),'records':records,'variants':summaries,
             'scope':protocol['direct_status'],'original_five_gate_goal_achieved':False})
        save(job,{'status':'FINISHED','pid':os.getpid(),'targets_verified':len(records)})
    except BaseException as error:
        save(job,{'status':'STOPPED_WITH_ERROR','pid':os.getpid(),'error':str(error),'completed_targets':sum(r['reached'] for r in records)})
        raise


if __name__ == '__main__':
    torch.set_num_threads(1)
    parser = argparse.ArgumentParser(); parser.add_argument('stage',choices=['prepare','direct']); args = parser.parse_args()
    {'prepare':prepare,'direct':direct}[args.stage]()
