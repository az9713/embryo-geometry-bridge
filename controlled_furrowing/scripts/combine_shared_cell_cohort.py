"""Join reviewed 60+8 outputs into one conditional all17 stationary-tuple claim.

This never fits a parameter or estimates noise. The two tolerances are explicit
illustrative queries chosen after the earlier M threshold was known.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,csv,hashlib,json
from joint_calibration_contract import OBSERVABLES

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs/shared_cell_cohort_v1'
OUT=ROOT/'evidence/shared_cell_cohort_inference.json'
read=lambda p:json.loads(Path(p).read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
rel=lambda p:Path(p).relative_to(ROOT).as_posix()
alltrue=lambda d:isinstance(d,dict) and bool(d) and all(v is True for v in d.values())
QUANTIFIER='for_all_parameters_in_closed_cell_exists_unique_root_in_declared_branch_domain'

FILES={
    'source':'data/derived/calibration_small_18um.csv',
    'context':'runs/cohort_uniform_roots_v1/protocol.json',
    'contract':'runs/cohort_static_observers_v1/model_observation_contract.json',
    'old_outputs':'runs/cohort_static_observers_v1/composition_final/results.json',
    'old_roots':'evidence/cohort_uniform_roots_independent_review.json',
    'old_application':'evidence/cohort_static_observers_independent_review.json',
    'old_application_meta':'evidence/cohort_static_application_meta_review.json',
    'new_outputs':'runs/cohort_recentered_observers_v1/results_scoped_v2.json',
    'new_outputs_original':'runs/cohort_recentered_observers_v1/results.json',
    'scope_correction':'runs/cohort_recentered_observers_v1/scope_correction.json',
    'new_roots':'evidence/cohort_recentered_roots_independent_review.json',
    'new_application':'evidence/recentered_observers_independent_review.json',
    'original_plan':'PLAN.md'}

def save(p,v):
    if p.exists():raise RuntimeError('preserve frozen inference artifact '+str(p))
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def freeze():
    files={value:sha(ROOT/value) for value in FILES.values()}
    for prefix in ('old','new'):
        output=read(ROOT/FILES[prefix+'_outputs']);roots=read(ROOT/FILES[prefix+'_roots'])
        for entry in output['records']:
            file=ROOT/entry['file']
            if sha(file)!=entry['sha256']:raise RuntimeError('changed observer result')
            files[entry['file']]=sha(file)
        for entry in roots['record_audits']:
            file=ROOT/entry['review_file']
            if sha(file)!=entry['review_sha256']:raise RuntimeError('changed independent root capsule')
            files[entry['review_file']]=sha(file)
        # Bind actual root provider files from the independently reviewed source.
        protocol=read(ROOT/('runs/cohort_uniform_roots_v2/protocol.json' if prefix=='old' else 'runs/cohort_recentered_roots_v1/protocol.json'))
        for row in protocol['records']:
            file=ROOT/row['output_folder']/'result.json';files[rel(file)]=sha(file)
    files[rel(Path(__file__))]=sha(__file__)
    save(RUN/'protocol.json',{'version':1,'file_hashes':files,'file_roles':FILES,'source_sha256':sha(__file__),
        'queries_common_delta_um_exact':['51/10','26/5'],'query_selection':'Illustrative conditional queries after the M threshold was known; neither fitted noise nor a prospective hypothesis test.',
        'target':'All68 uniform outputs on one common parameter box under independent stationary-root-per-record closure',
        'missing_policy':'Reject a complete claim if any record/channel or independent proof binding is missing; no interpolation.',
        'scope':'Exact finite binary-cache potential, fixed geometry/normalization/source doses and static current-radius disk/spatial-max observer. Full physical domain, source temporal validation and continuum accuracy remain incomplete.'})
    print(json.dumps({'status':'FROZEN_COMBINATION_PROTOCOL','inputs':len(files),'protocol_sha256':sha(RUN/'protocol.json')}))

def run():
    p=read(RUN/'protocol.json')
    if any(sha(ROOT/name)!=h for name,h in p['file_hashes'].items()):raise RuntimeError('frozen inputs changed')
    context=read(ROOT/FILES['context']);contract=read(ROOT/FILES['contract']);C=context['shared_parameter_cell'];cellhash=context['shared_parameter_cell_sha256'];contracthash=sha(ROOT/FILES['contract'])
    if contracthash!=context['common_model_observation_contract_sha256'] or contract['parameter_cell_sha256']!=cellhash:raise RuntimeError('context mismatch')
    with (ROOT/FILES['source']).open(encoding='utf-8-sig',newline='') as f:data=list(csv.DictReader(f))
    expected={r['record_id'] for r in data};outputs={};proofs=[];checks={}
    for prefix,allowed in [('old',expected-{'optoGEF2_small_col_F','optoGEF2_small_col_K'}),('new',{'optoGEF2_small_col_F','optoGEF2_small_col_K'})]:
        packet=read(ROOT/FILES[prefix+'_outputs']);rootreview=read(ROOT/FILES[prefix+'_roots']);application=read(ROOT/FILES[prefix+'_application'])
        if not application['status'].startswith('PASS') or not alltrue(application['checks']):raise RuntimeError('passing application review required')
        roots={entry['record_id']:entry['audit'] for entry in rootreview['record_audits']}
        applications={entry['record_id']:entry for entry in application['record_audits']}
        providerp=read(ROOT/('runs/cohort_uniform_roots_v2/protocol.json' if prefix=='old' else 'runs/cohort_recentered_roots_v1/protocol.json'))
        providers={r['record_id']:r for r in providerp['records']}
        for entry in packet['records']:
            rid=entry['record_id']
            if rid not in allowed:
                if prefix=='old' and rid in ('optoGEF2_small_col_F','optoGEF2_small_col_K'):
                    if read(ROOT/entry['file'])['status']!='UNRESOLVED':raise RuntimeError('original unresolved history changed')
                    continue
                raise RuntimeError('unexpected selected record')
            if rid in outputs:raise RuntimeError('duplicate selected record')
            o=read(ROOT/entry['file']);rpath=ROOT/providers[rid]['output_folder']/'result.json';root=read(rpath);audit=roots[rid];app=applications[rid]
            if audit['status']!='PASS' or not alltrue(audit['checks']) or not audit['uniform_quantifier_verified'] or audit['certificate_sha256']!=sha(rpath):raise RuntimeError('uniform root proof not bound')
            if not app['status'].startswith('PASS') or not alltrue(app['checks']) or app['root_certificate_sha256']!=sha(rpath) or app['observer_certificate_sha256']!=sha(ROOT/entry['file']):raise RuntimeError('observer application proof not bound')
            if root['parameter_cell_sha256']!=cellhash or root['parameter_cell']!=C or root['model_observation_contract_sha256']!=contracthash or root['quantifier']!=QUANTIFIER or root['branch_law_semantics']!='independent_stationary_root_per_record':raise RuntimeError('root family mismatch')
            if o['parameter_cell_sha256']!=cellhash or o['common_contract_sha256']!=contracthash or o['branch_domain_id']!=root['branch_domain_id']:raise RuntimeError('observer family mismatch')
            source=next(x for x in data if x['record_id']==rid)
            if Q(root['dose_exact'])!=Q(source['activation_magnitude']):raise RuntimeError('source dose mismatch')
            if set(o['prediction_intervals_exact'])!=set(OBSERVABLES):raise RuntimeError('all four channels required')
            outputs[rid]={k:list(map(Q,v)) for k,v in o['prediction_intervals_exact'].items()}
            proofs.append({'record_id':rid,'root_certificate_file':rel(rpath),'root_certificate_sha256':sha(rpath),
                'observer_certificate_file':entry['file'],'observer_certificate_sha256':sha(ROOT/entry['file']),
                'root_review_file':FILES[prefix+'_roots'],'root_review_sha256':sha(ROOT/FILES[prefix+'_roots']),
                'application_review_file':FILES[prefix+'_application'],'application_review_sha256':sha(ROOT/FILES[prefix+'_application']),
                'branch_domain_id':root['branch_domain_id'],'source_dose_exact':str(Q(source['activation_magnitude']))})
    if set(outputs)!=expected or len(expected)!=17:raise RuntimeError('missing/unexpected source record')
    rows=[]
    for source in data:
        rid=source['record_id']
        for k in OBSERVABLES:
            lo,hi=outputs[rid][k];y=Q(source[k])
            if lo>hi:raise RuntimeError('invalid prediction interval')
            rows.append({'record_id':rid,'observable':k,'observed_exact':str(y),'prediction_interval_exact':[str(lo),str(hi)],
                'minimum_residual_exact':str(max(lo-y,y-hi,Q(0))),'maximum_residual_exact':str(max(abs(lo-y),abs(hi-y)))})
    lower=max(Q(v['minimum_residual_exact']) for v in rows);upper=max(Q(v['maximum_residual_exact']) for v in rows)
    queries=[]
    for value in p['queries_common_delta_um_exact']:
        delta=Q(value);status='OUT_SELECTED_TUPLE' if delta<lower else 'IN_ALL_PARAMETERS_SELECTED_TUPLE' if delta>=upper else 'UNRESOLVED'
        queries.append({'common_delta_um_exact':str(delta),'selected_tuple_status':status,
            'all_branches_existential_parameter_status':'IN_STATIONARY_ROOT_TUPLE_WITNESS' if status.startswith('IN_') else 'UNRESOLVED',
            'source_noise_estimate':False,'common_physical_time_path_certified':False})
    checks={'all17_distinct_records':len(outputs)==17,'all68_channels_present':len(rows)==68,'all17_root_observer_proof_pairs_bound':len(proofs)==17,
        'nonempty_positive_width_common_box':all(Q(a)<Q(b) for a,b in C.values()),
        'correct_nested_FK_global_coverage_flag':read(ROOT/FILES['new_outputs'])['aggregate_F_K_only']['all68_uniform_static_tuple_enclosed'] is False,
        'unknown_source_temporal_semantics_retained':contract['source_temporal_semantics_validated'] is False,
        'all_inputs_unchanged_after_combination':all(sha(ROOT/name)==h for name,h in p['file_hashes'].items())}
    if not all(checks.values()):raise RuntimeError(checks)
    save(OUT,{'status':'CERTIFIED_CONDITIONAL_SHARED_CELL_STATIONARY_ENCLOSURES','checks':checks,
        'protocol_sha256':sha(RUN/'protocol.json'),'source_sha256':sha(__file__),'parameter_cell':C,'parameter_cell_sha256':cellhash,
        'model_observation_contract_sha256':contracthash,'records':17,'constraints':68,'uniform_constraints':rows,'proof_provenance':proofs,
        'maximum_residual_lower_exact':str(lower),'maximum_residual_upper_exact':str(upper),'range_width_exact':str(upper-lower),
        'display_um':{'lower':float(lower),'upper':float(upper),'width':float(upper-lower)},'queries':queries,
        'theorem':'For every theta in this declared common cell, a selected stationary root tuple exists and all68 static outputs obey their intervals. Its maximum absolute residual D(theta) lies in[lower,upper]. Allowances belowlower exclude this selected tuple throughoutC; allowances at leastupper admit it throughoutC. Other branches are not excluded by selected-tuple failure.',
        'branch_law_semantics':'independent_stationary_root_per_record','physical_time_path_or_stability_certified':False,
        'numerical_target':'Exact finite potential using frozen binary quadrature and geometry; not the ideal integral or continuum PDE.',
        'physiological_parameter_prior':False,'estimated_measurement_noise':False,'full_identified_set_computed':False,
        'source_temporal_semantics_validated':False,'full_five_gate_goal_complete':False,'confirmed_novelty':False})
    print(json.dumps({'status':read(OUT)['status'],'records':17,'constraints':68,'bounds_um':[float(lower),float(upper)],'queries':queries,'receipt_sha256':sha(OUT)}))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['freeze','run']);args=a.parse_args()
    freeze() if args.stage=='freeze' else run()
