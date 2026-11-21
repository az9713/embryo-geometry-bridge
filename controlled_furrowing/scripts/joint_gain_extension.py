"""Versioned computational-support extension and direct all-endpoint replay.

One 12x6 mechanical scenario; all four observation variants retained. Parent
preparation, physical kernels and numerical policies remain unchanged. No
held-out response files are read. Exact fits certify a midpoint-table surrogate,
not nonlinear mechanics, biological peak closure or a full identified set.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,bisect,hashlib,json,os,time
import numpy as np
import torch
from joint_gain_calibration import inputs,endpoints,VARIANTS,OBSERVABLES,TABLE_FILES
from piecewise_gain_minimax import fit
from joint_calibration_contract import evaluate,interval_residual
from sparse_quadratic import AxisymmetricModel,Parameters
from quadratic_quadrature_study import IntegrationPolicy,evaluate as mechanical_evaluate,exact_geometry

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/joint_gain_extension_v2'
PARENT=ROOT/'runs/joint_gain_calibration'
SCHEDULE=(20,40,80)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def relative_output(path):return Path(os.path.relpath(path,Path.cwd()))
def hashes():
    names=('joint_gain_extension.py','joint_gain_calibration.py','joint_calibration_contract.py','piecewise_gain_minimax.py',
           'quadratic_endpoint_observer.py','sparse_quadratic.py','quadratic_quadrature_study.py','robust_quadratic.py',
           'axisymmetric_quadratic.py','axisymmetric_mixed.py','axisymmetric_prestress.py','axisymmetric_mechanics.py',
           'mixed_geometry_certificate.py','mixed_geometry_certificate_v2.py')
    return {name:sha(Path(__file__).with_name(name)) for name in names}

def assert_protocol(protocol):
    if protocol['source_hashes']!=hashes():raise RuntimeError('locked source changed')
    for name,value in protocol['parent_artifact_hashes'].items():
        if sha(PARENT/name)!=value:raise RuntimeError('parent artifact changed: '+name)
    if sha(ROOT/'data/derived/calibration_small_18um.csv')!=protocol['calibration_csv_sha256']:raise RuntimeError('small-source table changed')
    for name,value in protocol['seed_field_hashes'].items():
        if sha(ROOT/name)!=value:raise RuntimeError('seed field changed: '+name)

def prepare():
    path=OUT/'protocol.json'
    if path.exists():raise RuntimeError('v2 preparation already exists; retain it rather than overwrite')
    parent=json.loads((PARENT/'protocol.json').read_text());seed=json.loads((PARENT/'response_table.json').read_text())
    candidates=json.loads((PARENT/'surrogate_candidates.json').read_text())
    if candidates['status']!='COMPLETE' or not all(v['surrogate_fit']['at_declared_gain_boundary'] for v in candidates['variants']):raise RuntimeError('expected retained v1 boundary diagnosis')
    rows,observed=inputs();doses=[Q(r['activation_magnitude']) for r in rows]
    fields={p.relative_to(ROOT).as_posix():sha(p) for p in TABLE_FILES.values()}
    fields[(PARENT/'zero/final.json').relative_to(ROOT).as_posix()]=sha(PARENT/'zero/final.json')
    protocol={'version':2,'scope':'One fixed mechanical scenario, all68small endpoints, four conditional observation variants; computational support extension before direct replay',
      'physical_parameters':parent['physical_parameters'],'integration_policy':{'volume_duffy_gauss_order':8,'surface_gauss_order':12,'volume_chunk_triangles':12},
      'source_hashes':hashes(),'parent_artifact_hashes':{name:sha(PARENT/name) for name in ('protocol.json','response_table.json','surrogate_candidates.json')},
      'seed_field_hashes':fields,'calibration_csv_sha256':sha(ROOT/'data/derived/calibration_small_18um.csv'),
      'maximum_dose_exact':str(max(doses)),'initial_strength_nodes':['0','1','5','10'],'computational_extension_schedule_exact':list(map(str,SCHEDULE)),
      'extension_rule':'Always add20; add40then80onlyif any retained variant still chooses a declared gain boundary. These are computational table-support limits, NOT physiological priors.',
      'gain_interval_rule':'0 <= g <= maximum_verified_table_strength/maxdose; no extrapolation',
      'common_assumptions':parent['common_assumptions'],'observation_variants':parent['observation_variants'],
      'surrogate_definition':'Exact piecewise-linear interpolation of interval-midpoint responses; physical nonlinear interpolation error is NOT enclosed',
      'candidate_review_gate':'Independent all-cell coverage/active-line/subgradient audit with candidate_file_sha256 and response_table_sha256 before direct acceptance',
      'numerical_gates':['last-verified-only adaptive continuation','Q8 force<1e-7','Q12 cross-force<1e-7','exact discrete positiveJ and global embedding','observer requested interval resolution'],
      'direct_rule':'Union of all distinct exact g_variant*dose targets, ascending from verifiedzero; all17records/all4endpoints for EACHvariant; every failed attempt and source/field hash retained',
      'not_claimed':['physiological gain bound','nonlinear global gain optimum','full material/geometry identified set','statistical noise or confidence','continuum error bound','validated temporal peak observation'],
      'held_out_response_files_read':False,'direct_v1_run':False}
    save(path,protocol)  # Freeze before extracting or computing any new node.
    save(OUT/'response_table.json',{'protocol_sha256':sha(path),'rows':seed['rows'],'source':'Immutable copiedv1seed; new verifiednodes appended only in this v2directory'})
    save(OUT/'preparation.json',{'status':'PREPARED_FROZEN_BEFORE_EXTENSION','protocol_sha256':sha(path),'source_hashes':hashes()})
    print(json.dumps({'stage':'prepared_v2','source_sha256':hashes()['joint_gain_extension.py'],'protocol_sha256':sha(path)}),flush=True)

def certify_state(model,result):
    if not result.get('reached_target',result.get('converged')) or not result['converged']:raise RuntimeError('requested nonlinear target not reached; no endpoint promoted')
    model._verify_state(result)
    other=AxisymmetricModel(Parameters(**result['parameters']),integration_policy=IntegrationPolicy(12,16))
    cross=mechanical_evaluate(other,result)
    geometry=exact_geometry(result)
    if cross['scaled_force_linf']>=1e-7:raise RuntimeError('Q12 cross-force fails; retained result is not promoted into response table')
    if not geometry['all_exact_J_positive'] or not geometry['global_embedding_certified']:raise RuntimeError('exact discrete geometry not certified; result retained but not promoted')
    values=endpoints(result)
    return {'Q8_force_linf':result['scaled_free_force_linf'],'Q12_cross_evaluation':cross,'geometry':geometry},values

def fit_all(table,protocol):
    rows,observed=inputs();dose=[Q(r['activation_magnitude']) for r in rows]
    strengths=[Q(r['strength_exact']) for r in table['rows']];variants=[]
    target=[[observed[r['record_id']][k] for k in OBSERVABLES] for r in rows]
    for weight,component in VARIANTS:
        name=weight+'__'+component
        response=[[(Q(r['variant_endpoint_intervals'][name][k][0])+Q(r['variant_endpoint_intervals'][name][k][1]))/2 for k in OBSERVABLES] for r in table['rows']]
        result=fit(strengths,response,dose,target,(0,max(strengths)/max(dose)))
        variants.append({'variant':name,'weighting':weight,'component_aggregation':component,'surrogate_fit':result})
        print(json.dumps({'stage':'exact_v2_surrogate','support':float(max(strengths)),'variant':name,'gain':result['gain_display'],'loss_um':result['loss_display'],'boundary':result['at_declared_gain_boundary']}),flush=True)
    return {'status':'ALL_INTERIOR_AWAIT_INDEPENDENT_REVIEW' if all(not v['surrogate_fit']['at_declared_gain_boundary'] for v in variants) else 'COMPUTATIONAL_BOUNDARY_REMAINS',
            'protocol_sha256':sha(OUT/'protocol.json'),'response_table_sha256':sha(OUT/'response_table.json'),'variants':variants,
            'maximum_verified_strength_exact':str(max(strengths)),'scope':protocol['surrogate_definition']}

def extend():
    protocol=json.loads((OUT/'protocol.json').read_text());assert_protocol(protocol)
    job=OUT/'extension_job.json'
    if job.exists():raise RuntimeError('extensionjob record exists; inspectlive/terminal state, do not duplicate')
    save(job,{'status':'RUNNING_EXTENSION','pid':os.getpid(),'source_hashes':hashes()})
    model=AxisymmetricModel(Parameters(**protocol['physical_parameters']))
    previous=json.loads(TABLE_FILES[10].read_text());model._verify_state(previous)
    table=json.loads((OUT/'response_table.json').read_text());records=[]
    try:
        for strength in SCHEDULE:
            assert_protocol(protocol);folder=OUT/'extension'/f'strength_{strength}';started=time.monotonic()
            result=model.advance(float(strength),previous,output_directory=relative_output(folder));save(folder/'final.json',result)
            row={'requested_strength':strength,'previous_verified_strength':previous['activation_strength'],'reached':result['reached_target'],
                 'file':(folder/'final.json').relative_to(ROOT).as_posix(),'sha256':sha(folder/'final.json'),
                 'attempts':len(result['robust_continuation']['attempts']),'rejected_attempts':sum(not a['accepted'] for a in result['robust_continuation']['attempts']),
                 'seconds':time.monotonic()-started}
            records.append(row);save(OUT/'extension_progress.json',{'status':'VERIFYING_NODE','records':records})
            verification,values=certify_state(model,result)
            save(folder/'verification.json',verification);row['verification_sha256']=sha(folder/'verification.json')
            table['rows'].append({'strength_exact':str(strength),'variant_endpoint_intervals':values,'field_file':row['file'],'field_sha256':row['sha256']})
            save(OUT/'response_table.json',table)
            fitted=fit_all(table,protocol);save(OUT/f'surrogate_candidates_strength_{strength}.json',fitted);save(OUT/'surrogate_candidates.json',fitted)
            save(OUT/'extension_progress.json',{'status':fitted['status'],'records':records,'response_table_sha256':sha(OUT/'response_table.json'),'candidate_file_sha256':sha(OUT/'surrogate_candidates.json')})
            assert_protocol(protocol)
            if fitted['status']=='ALL_INTERIOR_AWAIT_INDEPENDENT_REVIEW':
                save(job,{'status':'FINISHED_AWAIT_INDEPENDENT_REVIEW','pid':os.getpid(),'support':strength,'source_hashes':hashes()})
                print(json.dumps({'stage':'extension_ready_for_audit','support':strength,'candidate_file_sha256':sha(OUT/'surrogate_candidates.json'),'response_table_sha256':sha(OUT/'response_table.json')}),flush=True);return
            previous=result
        save(job,{'status':'COMPUTATIONAL_SUPPORT_EXHAUSTED','pid':os.getpid(),'support':SCHEDULE[-1],'source_hashes':hashes()})
        print(json.dumps({'stage':'support_exhausted','gain_not_inferred':True}),flush=True)
    except BaseException as error:
        save(job,{'status':'STOPPED_WITH_ERROR','pid':os.getpid(),'error':str(error),'records':records,'source_hashes':hashes()});raise

def surrogate_at(table,variant,strength):
    nodes=[Q(r['strength_exact']) for r in table['rows']];strength=Q(strength)
    if not nodes[0]<=strength<=nodes[-1]:raise ValueError('no surrogate extrapolation allowed')
    j=min(bisect.bisect_right(nodes,strength)-1,len(nodes)-2);weight=(strength-nodes[j])/(nodes[j+1]-nodes[j]);out={}
    for key in OBSERVABLES:
        values=[sum(map(Q,table['rows'][i]['variant_endpoint_intervals'][variant][key]))/2 for i in (j,j+1)]
        out[key]=values[0]+weight*(values[1]-values[0])
    return out

def direct(review_path):
    protocol=json.loads((OUT/'protocol.json').read_text());assert_protocol(protocol)
    candidates=json.loads((OUT/'surrogate_candidates.json').read_text());table=json.loads((OUT/'response_table.json').read_text())
    review=json.loads(Path(review_path).read_text())
    if candidates['status']!='ALL_INTERIOR_AWAIT_INDEPENDENT_REVIEW':raise RuntimeError('all four candidates must be interior before direct replay')
    if review.get('status')!='PASS' or not review.get('checks') or not all(review['checks'].values()):raise RuntimeError('independent exact surrogate audit has not passed')
    if review.get('candidate_file_sha256')!=sha(OUT/'surrogate_candidates.json') or review.get('response_table_sha256')!=sha(OUT/'response_table.json'):raise RuntimeError('independent audit does not bind these candidates/table')
    if review.get('protocol_sha256')!=sha(OUT/'protocol.json') or review.get('optimizer_source_sha256')!=sha(Path(__file__).with_name('piecewise_gain_minimax.py')):raise RuntimeError('independent audit does not bind this protocol/optimizer')
    if candidates['protocol_sha256']!=sha(OUT/'protocol.json') or candidates['response_table_sha256']!=sha(OUT/'response_table.json'):raise RuntimeError('candidate hashes changed')
    if json.loads((OUT/'extension_job.json').read_text()).get('status')!='FINISHED_AWAIT_INDEPENDENT_REVIEW':raise RuntimeError('extension must be terminal and awaiting independent review')
    job=OUT/'direct_job.json'
    if job.exists():raise RuntimeError('directjob record exists; inspect live/terminal state instead of duplication')
    rows,observed=inputs();requests={};predictions={v['variant']:{} for v in candidates['variants']}
    for variant in candidates['variants']:
        gain=Q(variant['surrogate_fit']['gain_exact'])
        for row in rows:requests.setdefault(gain*Q(row['activation_magnitude']),[]).append((variant['variant'],row['record_id']))
    direct_protocol={'protocol_sha256':sha(OUT/'protocol.json'),'candidate_file_sha256':sha(OUT/'surrogate_candidates.json'),
        'response_table_sha256':sha(OUT/'response_table.json'),'independent_review_sha256':sha(review_path),'source_hashes':hashes(),
        'total_distinct_exact_targets':len(requests),'requests':[{'strength_exact':str(s),'actual_solver_strength_exact':str(Q.from_float(float(s))),
            'uses':[{'variant':v,'record_id':rid} for v,rid in uses]} for s,uses in sorted(requests.items())],
        'all_observations_per_variant':68,'variants':4,'note':'Same68source observations under4declared contracts, NOT4independent data sets.'}
    save(OUT/'direct_protocol.json',direct_protocol);save(job,{'status':'RUNNING_DIRECT','pid':os.getpid(),'source_hashes':hashes()})
    model=AxisymmetricModel(Parameters(**protocol['physical_parameters']));previous=json.loads((PARENT/'zero/final.json').read_text());model._verify_state(previous)
    records=[];started=time.monotonic()
    try:
        for index,(strength,uses) in enumerate(sorted(requests.items())):
            assert_protocol(protocol);folder=OUT/'direct'/f'{index:03d}';tick=time.monotonic()
            result=model.advance(float(strength),previous,output_directory=relative_output(folder));save(folder/'final.json',result)
            record={'index':index,'requested_strength_exact':str(strength),'solver_input_strength_exact':str(Q.from_float(float(strength))),
                'strength_rounding_difference_exact':str(Q.from_float(float(strength))-strength),'actual_verified_strength':result['activation_strength'],
                'uses':[{'variant':v,'record_id':rid} for v,rid in uses],'previous_verified_strength':previous['activation_strength'],
                'file':(folder/'final.json').relative_to(ROOT).as_posix(),'sha256':sha(folder/'final.json'),'reached':result['reached_target'],
                'attempts':len(result['robust_continuation']['attempts']),'rejected_attempts':sum(not a['accepted'] for a in result['robust_continuation']['attempts']),
                'seconds':time.monotonic()-tick}
            records.append(record);save(OUT/'direct_progress.json',{'status':'VERIFYING_TARGET','records':records,'predictions':predictions,'total_targets':len(requests)})
            verification,values=certify_state(model,result);save(folder/'verification.json',verification)
            record['force_Q8']=verification['Q8_force_linf'];record['force_Q12']=verification['Q12_cross_evaluation']['scaled_force_linf'];record['verification_sha256']=sha(folder/'verification.json')
            for variant,rid in uses:predictions[variant][rid]=values[variant]
            previous=result
            save(OUT/'direct_progress.json',{'status':'COMPLETE' if index+1==len(requests) else 'RUNNING','records':records,'predictions':predictions,'total_targets':len(requests),'elapsed_seconds':time.monotonic()-started})
            print(json.dumps({'stage':'direct_v2','index':index+1,'of':len(requests),'strength':float(strength),'force_Q8':record['force_Q8'],'force_Q12':record['force_Q12'],'rejected':record['rejected_attempts']}),flush=True)
        summaries=[]
        byid={r['record_id']:r for r in rows}
        for variant in candidates['variants']:
            key=variant['variant'];pred=predictions[key];residual=evaluate(observed,pred,{k:0 for k in OBSERVABLES})
            lower=max(map(Q,residual['discrepancy_lower_requirements_exact'].values()));upper=max(map(Q,residual['discrepancy_sufficient_requirements_exact'].values()));errors=[]
            gain=Q(variant['surrogate_fit']['gain_exact'])
            for rid,values in pred.items():
                strength=gain*Q(byid[rid]['activation_magnitude']);approx=surrogate_at(table,key,strength)
                for endpoint in OBSERVABLES:
                    lo,hi=interval_residual(approx[endpoint],values[endpoint]);errors.append({'record_id':rid,'observable':endpoint,
                        'surrogate_midpoint_prediction_exact':str(approx[endpoint]),'direct_interval_exact':values[endpoint],
                        'absolute_surrogate_to_direct_error_interval_exact':[str(lo),str(hi)]})
            summaries.append({'variant':key,'selected_exact_SURROGATE_gain':variant['surrogate_fit']['gain_exact'],'surrogate_loss_exact':variant['surrogate_fit']['loss_exact'],
                'direct_candidate_discrepancy_requirement_um_exact':[str(lower),str(upper)],'direct_candidate_discrepancy_requirement_um_display':[float(lower),float(upper)],
                'all68_residuals':residual,'all68_surrogate_direct_errors':errors,
                'maximum_surrogate_to_direct_error_upper_um':float(max(Q(v['absolute_surrogate_to_direct_error_interval_exact'][1]) for v in errors)),
                'scope':'Actual static candidate replay, NOT a nonlinear optimal gain/full physical identified set/statistical noise estimate.'})
        assert_protocol(protocol)
        result={'status':'COMPLETE_ALL_FOUR_CONDITIONAL_CANDIDATE_REPLAYS','protocol_sha256':sha(OUT/'protocol.json'),'direct_protocol_sha256':sha(OUT/'direct_protocol.json'),
            'records':records,'variants':summaries,'original_five_gate_goal_achieved':False,'source_hashes':hashes()}
        save(OUT/'direct_results.json',result);save(ROOT/'evidence/joint_gain_extension_direct.json',result)
        save(job,{'status':'FINISHED','pid':os.getpid(),'targets_verified':len(records),'source_hashes':hashes()})
    except BaseException as error:
        save(job,{'status':'STOPPED_WITH_ERROR','pid':os.getpid(),'error':str(error),'targets_verified':sum('verification_sha256' in r for r in records),'source_hashes':hashes()});raise

if __name__=='__main__':
    torch.set_num_threads(1);parser=argparse.ArgumentParser();parser.add_argument('stage',choices=('prepare','extend','direct'));parser.add_argument('--review',type=Path);args=parser.parse_args()
    if args.stage=='prepare':prepare()
    elif args.stage=='extend':extend()
    elif args.review is None:raise ValueError('direct stage requires an independent audit receipt')
    else:direct(args.review)
