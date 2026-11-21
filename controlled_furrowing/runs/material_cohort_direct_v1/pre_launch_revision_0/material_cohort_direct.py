"""Direct all-17 verification of one predeclared multivariate material proposal.

Primary disk/spatial-max continuation only. No new fitting or held-out outcomes.
All four observation conventions are evaluated at the SAME new material fields.
"""
from pathlib import Path
from fractions import Fraction as Q
from dataclasses import replace,asdict
import argparse,hashlib,json,math,os,time
import torch
from sparse_quadratic import AxisymmetricModel,Parameters
from joint_gain_extension import certify_state,hashes as parent_hashes
from joint_calibration_contract import evaluate,OBSERVABLES

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/material_cohort_direct_v1'
COHORT=ROOT/'runs/material_cohort_descent_v1'
CANDIDATES=COHORT/'exact_extraction_v2/exact_candidates.json'
AFFINE_REVIEW=ROOT/'evidence/bounded_affine_minimax_independent_review.json'
PRIMARY='current_disk__spatial_max'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def save(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def hashes():
    values=parent_hashes()
    for name in ('joint_gain_extension.py','material_cohort_direct.py','material_cohort_descent.py',
                 'material_tangent_inference.py','quadratic_mean_derivatives.py','bounded_affine_minimax.py'):
        values[name]=sha(Path(__file__).with_name(name))
    return values

def proposal(parameters,base_gain,step):
    if len(step)!=5:raise ValueError('exact five-coordinate proposal required')
    x=list(map(Q,step));bounds=[(-Q.from_float(.1),Q.from_float(.1))]*5
    bounds[2]=(Q(0),Q.from_float(.02))
    if any(not a<=v<=b for v,(a,b) in zip(x,bounds)):raise ValueError('proposal exceeds declared numerical trust step')
    p=Parameters(**parameters)
    if p.bulk!=5000 or p.apical_surface_tension!=0:raise ValueError('fixed baseline compression/zero apical tension required')
    # Stored binary64 exp/multiplication values define this numerical proposal;
    # no exact transcendental or physical parameter measurement is claimed.
    changes={'apical_area_modulus':p.apical_area_modulus*math.exp(float(x[0])),
             'basal_area_modulus':p.basal_area_modulus*math.exp(float(x[1])),
             'apical_surface_tension':p.apical_surface_tension+float(x[2]),
             'basal_surface_tension':p.basal_surface_tension*math.exp(float(x[3]))}
    new=replace(p,**changes);gain=Q(base_gain)*Q.from_float(math.exp(float(x[4])))
    if gain<=0 or any(not math.isfinite(v) or v<=0 for k,v in changes.items() if k!='apical_surface_tension'):
        raise ValueError('finite positive transformed material/gain required')
    return asdict(new),gain

def tests():
    checks={};p=asdict(Parameters(bulk=5000,apical_surface_tension=0,basal_surface_tension=2,
        apical_area_modulus=.2,basal_area_modulus=2))
    same,g=proposal(p,'8',['0']*5)
    checks['identity_step_preserves_parameters_gain']=same==p and g==8
    updated,g=proposal(p,'8',['0','0','0','0','-1/20'])
    checks['gain_transform_and_fixed_penalty']=g==8*Q.from_float(math.exp(-.05)) and updated['bulk']==5000
    allowed=set(('apical_area_modulus','basal_area_modulus','apical_surface_tension','basal_surface_tension'))
    checks['geometry_and_other_materials_unchanged']=all(updated[k]==v for k,v in p.items() if k not in allowed)
    for i,call in enumerate([lambda:proposal(p,'8',['0']*4),lambda:proposal(p,'8',['0','0','-1/1000','0','0']),
        lambda:proposal(p,'8',['1','0','0','0','0']),lambda:proposal({**p,'bulk':50},'8',['0']*5)]):
        try:call();checks['input_guard_'+str(i)]=False
        except ValueError:checks['input_guard_'+str(i)]=True
    save(OUT/'synthetic_verification.json',{'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'source_sha256':sha(__file__),
        'scope':'Proposal transform/guards only; no equilibrium computed'})
    print(json.dumps({'stage':'tests','pass':all(checks.values()),'checks':len(checks)}),flush=True)
    if not all(checks.values()):raise RuntimeError('direct proposal tests failed')

def prepare():
    if (OUT/'protocol.json').exists():raise RuntimeError('protocol exists; preserve instead of overwrite')
    test=read(OUT/'synthetic_verification.json')
    if test['status']!='PASS' or test['source_sha256']!=sha(__file__):raise RuntimeError('matching synthetic checks required')
    cohort=read(COHORT/'protocol.json');old=read(COHORT/'results.json');candidates=read(CANDIDATES)
    if old['status']!='COMPLETE_17_RECORD_LOCAL_AFFINE_DESCENT' or candidates['status']!='PASS_ALL_FOUR_EXACT_AFFINE_CANDIDATES':raise RuntimeError('completed tangent/exact affine artifacts required')
    c=next(v for v in candidates['candidates'] if v['variant']==PRIMARY)
    if not all(c['exact_affine_result']['checks'].values()):raise RuntimeError('primary exact affine certificate fails')
    params,gain=proposal(cohort['parameters'],cohort['baseline_gain_exact'],c['exact_affine_result']['step_exact'])
    records=[]
    for r in cohort['records']:
        exact=gain*Q(r['dose_exact']);records.append({**r,'requested_strength_exact':str(exact),
            'solver_input_strength_exact':str(Q.from_float(float(exact))),
            'strength_rounding_difference_exact':str(Q.from_float(float(exact))-exact)})
    records.sort(key=lambda r:Q(r['requested_strength_exact']))
    if len(records)!=17 or len({r['record_id'] for r in records})!=17:raise RuntimeError('all17 distinct source records required')
    protocol={'version':1,'scope':'Direct nonlinear verification of ONE predeclared primary material/gain step; no new optimization or full identified-set calculation',
        'selection_rule':'Continue the primary disk/spatial-max baseline selected by prior small-calibration discrepancy; no selection from larger outcomes or new nonlinear results',
        'primary_variant':PRIMARY,'source_hashes':hashes(),'cohort_protocol_sha256':sha(COHORT/'protocol.json'),
        'cohort_results_sha256':sha(COHORT/'results.json'),'candidate_file_sha256':sha(CANDIDATES),
        'exact_extraction_protocol_sha256':sha(CANDIDATES.parent/'protocol.json'),
        'base_parameters':cohort['parameters'],'new_parameters':params,'base_gain_exact':cohort['baseline_gain_exact'],
        'new_gain_exact':str(gain),'new_gain_display':float(gain),'primary_step_exact':c['exact_affine_result']['step_exact'],
        'transformation_rule':'Binary64 exp(step_display) and stored binary64 material multiplication; gain uses exact oldgain times binary64exp. Requested/actualbinary64strength and rational rounding retained.',
        'primary_affine_baseline_loss_um':c['baseline_loss_um'],'primary_affine_predicted_loss_exact':c['exact_affine_result']['loss_exact'],
        'records':records,'observation_variants':[v['variant'] for v in candidates['candidates']],
        'integration_policy':{'volume_duffy_gauss_order':8,'surface_gauss_order':12,'volume_chunk_triangles':12},
        'higher_rule_crosscheck':{'volume_duffy_gauss_order':12,'surface_gauss_order':16},
        'continuation':'Generate and verify NEW-material zero; advance ascending requested input from last verified state only. Original-material states are never promoted as verified initial states.',
        'required_gates':['source/input hashes unchanged','Q8force<1e-7','Q12force<1e-7','exact discrete positiveJ/globalembedding','all4 observer intervals resolved'],
        'affine_review_required_fields':['status=PASS','checks alltrue','candidate_file_sha256','cohort_results_sha256','cohort_protocol_sha256','helper_source_sha256','exact_extraction_protocol_sha256'],
        'prelaunch_review_required_fields':['status=PASS','checks alltrue','driver_sha256','protocol_sha256','candidate_file_sha256','affine_review_sha256'],
        'common_assumptions':'Same static endpoint-as-source-peak proxy, flat component references, currentROI and source-column pairing as earlier study; no new temporal validation',
        'heldout_response_files_read':False,'no_refitting_in_direct_stage':True,'other_three_affine_proposals_directly_verified':False,
        'not_claimed':['nonlinear optimum','global/stable/unique equilibrium','uniform forward/PDE bound','full physical Theta_delta','new biological mechanism','novelty']}
    save(OUT/'protocol.json',protocol)
    print(json.dumps({'stage':'prepared','source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),
        'new_parameters':params,'new_gain':float(gain),'records':17}),flush=True)

def assert_inputs(p):
    if p['source_hashes']!=hashes():raise RuntimeError('frozen driver/dependency source changed')
    for path,key in [(COHORT/'protocol.json','cohort_protocol_sha256'),(COHORT/'results.json','cohort_results_sha256'),
                     (CANDIDATES,'candidate_file_sha256'),(CANDIDATES.parent/'protocol.json','exact_extraction_protocol_sha256')]:
        if sha(path)!=p[key]:raise RuntimeError('frozen candidate/cohort input changed')
    for r in p['records']:
        if sha(ROOT/r['field'])!=r['field_sha256']:raise RuntimeError('baseline field changed')

def prelaunch(p,path):
    a=read(AFFINE_REVIEW)
    if a['status']!='PASS' or not a.get('checks') or not all(v is True for v in a['checks'].values()):raise RuntimeError('independent actual affine review has not passed')
    bindings={'candidate_file_sha256':p['candidate_file_sha256'],'cohort_results_sha256':p['cohort_results_sha256'],
        'cohort_protocol_sha256':p['cohort_protocol_sha256'],'helper_source_sha256':p['source_hashes']['bounded_affine_minimax.py'],
        'exact_extraction_protocol_sha256':p['exact_extraction_protocol_sha256']}
    if any(a.get(k)!=v for k,v in bindings.items()):raise RuntimeError('independent affine audit does not bind these exact inputs')
    review=read(path)
    if review.get('status')!='PASS' or not review.get('checks') or not all(v is True for v in review['checks'].values()):raise RuntimeError('concrete driver/protocol prelaunch review required')
    wanted={'driver_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),
            'candidate_file_sha256':p['candidate_file_sha256'],'affine_review_sha256':sha(AFFINE_REVIEW)}
    if any(review.get(k)!=v for k,v in wanted.items()):raise RuntimeError('prelaunch review does not bind current source/protocol/affine audit')

def run(review_path):
    p=read(OUT/'protocol.json');assert_inputs(p);prelaunch(p,review_path)
    if (OUT/'job.json').exists():raise RuntimeError('job record exists; inspect live/terminal state, never duplicate')
    save(OUT/'job.json',{'status':'RUNNING','pid':os.getpid(),'protocol_sha256':sha(OUT/'protocol.json'),'prelaunch_review_sha256':sha(review_path)})
    model=AxisymmetricModel(Parameters(**p['new_parameters']));rows=[];predictions={v:{} for v in p['observation_variants']};started=time.monotonic()
    try:
        zero=model.advance(0.,output_directory=(OUT/'zero').relative_to(ROOT.parent))
        save(OUT/'zero/final.json',zero);verification,values=certify_state(model,zero);save(OUT/'zero/verification.json',verification);previous=zero
        for i,r in enumerate(p['records']):
            assert_inputs(p);folder=OUT/'direct'/f'{i:02d}';strength=Q(r['requested_strength_exact']);tick=time.monotonic()
            state=model.advance(float(strength),previous,output_directory=folder.relative_to(ROOT.parent));save(folder/'final.json',state)
            record={'index':i,'record_id':r['record_id'],'dose_exact':r['dose_exact'],'requested_strength_exact':str(strength),
                'solver_input_strength_exact':r['solver_input_strength_exact'],'strength_rounding_difference_exact':r['strength_rounding_difference_exact'],
                'previous_verified_strength':previous['activation_strength'],'actual_strength':state['activation_strength'],
                'file':(folder/'final.json').relative_to(ROOT).as_posix(),'sha256':sha(folder/'final.json'),'seconds':time.monotonic()-tick,
                'attempts':len(state['robust_continuation']['attempts']),'rejected_attempts':sum(not a['accepted'] for a in state['robust_continuation']['attempts'])}
            rows.append(record);save(OUT/'progress.json',{'status':'VERIFYING','records':rows,'predictions':predictions,'total':17})
            verification,values=certify_state(model,state);save(folder/'verification.json',verification)
            record.update({'verification_sha256':sha(folder/'verification.json'),'force_Q8':verification['Q8_force_linf'],
                'force_Q12':verification['Q12_cross_evaluation']['scaled_force_linf']})
            for v in predictions:predictions[v][r['record_id']]=values[v]
            previous=state;save(OUT/'progress.json',{'status':'RUNNING' if i<16 else 'COMPLETE','records':rows,'predictions':predictions,'total':17,'elapsed_seconds':time.monotonic()-started})
            print(json.dumps({'stage':'material_direct','record':i+1,'of':17,'strength':float(strength),'rejected':record['rejected_attempts'],'force_Q12':record['force_Q12']}),flush=True)
        observed={r['record_id']:r['observed_exact'] for r in p['records']};summaries=[]
        for name,pred in predictions.items():
            residual=evaluate(observed,pred,{k:0 for k in OBSERVABLES})
            lo=max(map(Q,residual['discrepancy_lower_requirements_exact'].values()));hi=max(map(Q,residual['discrepancy_sufficient_requirements_exact'].values()))
            summaries.append({'variant':name,'all68_residuals':residual,'candidate_discrepancy_interval_exact':[str(lo),str(hi)],
                'candidate_discrepancy_interval_um':[float(lo),float(hi)],'scope':'All variants at SAME newly equilibrated primary material step; other3variant-specific affine steps not solved'})
        assert_inputs(p)
        result={'status':'COMPLETE_PRIMARY_MATERIAL_STEP_DIRECT_REPLAY','protocol_sha256':sha(OUT/'protocol.json'),'prelaunch_review_sha256':sha(review_path),
            'source_hashes':hashes(),'records':rows,'variants':summaries,'zero_field_sha256':sha(OUT/'zero/final.json'),
            'primary_affine_predicted_loss_exact':p['primary_affine_predicted_loss_exact'],'original_five_gate_goal_achieved':False,
            'scope':p['scope'],'elapsed_seconds':time.monotonic()-started}
        save(OUT/'results.json',result);save(OUT/'job.json',{'status':'FINISHED','pid':os.getpid(),'records_verified':17,'results_sha256':sha(OUT/'results.json')})
        print(json.dumps({'stage':'complete_material_direct','records':17,'results_sha256':sha(OUT/'results.json')}),flush=True)
    except BaseException as error:
        save(OUT/'job.json',{'status':'STOPPED_WITH_ERROR','pid':os.getpid(),'error':str(error),'verified_records':sum('verification_sha256' in r for r in rows),'source_hashes':hashes()});raise

if __name__=='__main__':
    torch.set_num_threads(1);parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['test','prepare','run']);parser.add_argument('--review',type=Path);args=parser.parse_args()
    if args.stage=='test':tests()
    elif args.stage=='prepare':prepare()
    elif args.review is None:raise ValueError('run requires concrete prelaunch review')
    else:run(args.review)
