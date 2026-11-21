"""All-17-record local material/gain descent, with no new equilibrium solves.

Floating implicit tangents and graph finite differences define a local affine
surrogate. Exact rational primal/dual bounds certify that surrogate LP only.
Synthetic displaced graphs are never labelled equilibrium solutions.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,copy,csv,hashlib,json,os,time
import numpy as np
import torch
from scipy.optimize import linprog
from material_tangent_inference import mechanical_tangent,COLUMNS,check_sources
from quadratic_endpoint_observer import _state_controls,observe_graphs
from quadratic_mean_derivatives import mean_state_gradients
from joint_gain_calibration import OBSERVABLES,VARIANTS

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/material_cohort_descent_v1'
BASE=ROOT/'runs/joint_gain_extension_v2'
SELECTED='current_disk__spatial_max'
ACTIVE=(1,2,3,4,5)
COORDS=tuple(COLUMNS[i] for i in ACTIVE)
BOUNDS=((-0.1,0.1),(-0.1,0.1),(0.,0.02),(-0.1,0.1),(-0.1,0.1))
L0=Q('33.5')
FD_H=(1e-3,5e-4)
MAX_SHRINK=4
OBSERVER_TOL='0.0000001'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def save(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def hashes():
    names=list(check_sources())+['material_tangent_inference.py','quadratic_mean_derivatives.py',
        'quadratic_endpoint_observer.py','joint_gain_calibration.py','joint_calibration_contract.py',
        'annulus_observer.py','material_cohort_descent.py']
    return {n:sha(Path(__file__).with_name(n)) for n in names}
def qfloat(v):return Q.from_float(float(v))
def extract(out):
    if not out['maximum_shortening_certificate']['requested_width_achieved']:
        raise ValueError('graph maximum interval resolution not achieved')
    values={}
    for w,c in VARIANTS:
        means=out['mean_alternatives'][w]
        a,b=((out['spatial_max_components']['apical_downward'],out['spatial_max_components']['basal_upward'])
             if c=='spatial_max' else (means['apical_regional_mean_downward'],means['basal_regional_mean_upward']))
        seq=(out['four_endpoints']['maximum_Eulerian_height_shortening'],means['height_regional_mean_shortening'],a,b)
        values[w+'__'+c]={k:[str(v.lo),str(v.hi)] for k,v in zip(OBSERVABLES,seq)}
    return values
def graphs(state,current=None):
    trial=copy.deepcopy(state)
    if current is not None:trial['current_nodes']=np.asarray(current).tolist()
    for surface in ('top','bottom'):trial.pop('boundary_'+surface+'_quadratic_controls',None)
    controls,H=_state_controls(trial,L0)
    # Direct graph observer: this perturbed coefficient vector is NOT a root.
    return extract(observe_graphs(controls['top'],controls['bottom'],
        roi_um=(0,qfloat(state['parameters']['patch_radius'])*L0),reference_apical_um=H,
        reference_basal_um=0,tolerance_um=OBSERVER_TOL,maximum_cells=32768))

def affine_minimax(base,jacobian,observed,bounds=BOUNDS):
    """Exact weak-dual bounds for the explicitly rounded affine LP.

    r and D are exact rationals for the saved affine model. Any nonnegative
    multiplier with sum<=1 gives lower=min_box (D_signed^T lambda).d
    +lambda.r_signed. A feasible step gives upper=max|r+D.d|.
    """
    r=[Q(a)-Q(b) for a,b in zip(base,observed)]
    D=[[qfloat(x) for x in row] for row in jacobian]
    if not r or len(D)!=len(r) or any(len(row)!=len(bounds) for row in D):raise ValueError('matching nonempty affine dimensions required')
    box=[tuple(qfloat(x) for x in pair) for pair in bounds]
    if any(a>b for a,b in box):raise ValueError('ordered trust box required')
    signedD=D+[[-x for x in row] for row in D];signedr=r+[-x for x in r]
    A=np.array([[float(x) for x in row]+[-1.] for row in signedD]);b=-np.array(list(map(float,signedr)))
    result=linprog([0.]*len(bounds)+[1.],A_ub=A,b_ub=b,bounds=list(bounds)+[(0,None)],method='highs')
    if not result.success:raise RuntimeError('affine minimax LP failed: '+result.message)
    step=[min(max(qfloat(v),lo),hi) for v,(lo,hi) in zip(result.x[:-1],box)]
    residual=[rr+sum(a*x for a,x in zip(row,step)) for rr,row in zip(r,D)]
    upper=max(map(abs,residual))
    lam=[max(Q(0),qfloat(-x)) for x in result.ineqlin.marginals]
    mass=sum(lam)
    if mass>1:lam=[v/mass for v in lam]
    c=[sum(lam[i]*signedD[i][j] for i in range(len(lam))) for j in range(len(bounds))]
    lower=sum(lam[i]*signedr[i] for i in range(len(lam)))+sum(min(v*lo,v*hi) for v,(lo,hi) in zip(c,box))
    gap=upper-lower
    checks={'step_within_exact_trust_box':all(lo<=x<=hi for x,(lo,hi) in zip(step,box)),
        'dual_nonnegative_mass_at_most_one':all(x>=0 for x in lam) and sum(lam)<=1,
        'exact_weak_duality':gap>=0,'exact_affine_optimality_gap_below_1e_minus8_um':gap<=Q('0.00000001')}
    if not all(checks.values()):raise RuntimeError('exact affine primal/dual certificate failed')
    return {'step_exact':list(map(str,step)),'step_display':list(map(float,step)),
        'affine_optimum_interval_exact':[str(lower),str(upper)],'affine_optimum_interval_um':[float(lower),float(upper)],
        'exact_gap_um':str(gap),'baseline_loss_um':float(max(map(abs,r))),
        'affine_predicted_residuals_exact':list(map(str,residual)),
        'dual_multipliers_exact':list(map(str,lam)),
        'affine_residual_intercepts_exact':list(map(str,r)),
        'affine_derivative_coefficients_exact':[[str(x) for x in row] for row in D],
        'trust_box_exact':[[str(x) for x in pair] for pair in box],
        'checks':checks,'scipy_status':result.message,
        'scope':'Exact primal/weak-dual certificate for this rounded finite affine LP ONLY; no derivative/root/nonlinear/physical identified-set certificate.'}

def tests():
    OUT.mkdir(parents=True,exist_ok=True);checks={}
    a=affine_minimax(['0','2'],[[1],[1]],['0','0'],bounds=((-2,2),))
    checks['analytic_two_constraint_minimax']=abs(a['step_display'][0]+1)<1e-12 and a['affine_optimum_interval_um'][0]<=1<=a['affine_optimum_interval_um'][1]
    b=affine_minimax(['-2'],[[1]],['0'],bounds=((0,.02),))
    checks['one_sided_trust_boundary']=b['step_exact'][0]==str(qfloat(.02))
    z=affine_minimax(['3'],[[0]],['0'],bounds=((-1,1),))
    checks['zero_tangent_does_not_invent_improvement']=z['affine_optimum_interval_exact']==['3','3']
    controls=lambda z:[[(Q(0),Q(z)),(Q(1,2),Q(z)),(Q(1),Q(z))]]
    obs=lambda a,b:extract(observe_graphs(controls(a),controls(b),roi_um=(0,1),reference_apical_um=10,tolerance_um=OBSERVER_TOL))
    plus=obs(Q(8)-Q(1,1000),Q(1)+Q(2,1000));minus=obs(Q(8)+Q(1,1000),Q(1)-Q(2,1000))
    for variant in plus:
        expected=(3,3,1,2)
        checks['constant_graph_tangent_'+variant]=all((Q(plus[variant][k][0])-Q(minus[variant][k][1]))/Q(2,1000)==v for k,v in zip(OBSERVABLES,expected))
    checks['nonsmooth_max_warning_example']=(max(.001,-.001)-max(-.001,.001))/(.002)==0 and (max(.001,-.001)-0)/.001==1
    for i,fn in enumerate([lambda:affine_minimax([],[],[]),lambda:affine_minimax(['0'],[[1]],['0'],bounds=((1,-1),))]):
        try:fn();checks['input_guard_'+str(i)]=False
        except ValueError:checks['input_guard_'+str(i)]=True
    save(OUT/'synthetic_verification.json',{'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'source_sha256':sha(__file__),'examples':[a,b,z],
         'scope':'Synthetic affine optimization and graph derivative plumbing only; no new mechanical equilibrium.'})
    print(json.dumps({'stage':'synthetics','checks':len(checks),'pass':all(checks.values())}),flush=True)
    if not all(checks.values()):raise RuntimeError('cohort descent synthetics failed')

def prepare():
    if (OUT/'protocol.json').exists():raise RuntimeError('protocol exists; do not overwrite or duplicate')
    test=read(OUT/'synthetic_verification.json')
    if test['status']!='PASS' or test['source_sha256']!=sha(__file__):raise RuntimeError('matching synthetic checks required')
    final=read(BASE/'direct_results.json');review=read(ROOT/'evidence/joint_direct_independent_review.json')
    if final['status']!='COMPLETE_ALL_FOUR_CONDITIONAL_CANDIDATE_REPLAYS' or review['status']!='PASS' or not all(review['final_checks'].values()):raise RuntimeError('completed independently reviewed direct fields required')
    best=next(v for v in final['variants'] if v['variant']==SELECTED)
    if any(Q(best['direct_candidate_discrepancy_requirement_um_exact'][1])>=Q(v['direct_candidate_discrepancy_requirement_um_exact'][0]) for v in final['variants'] if v['variant']!=SELECTED):raise RuntimeError('declared baseline not robustly best calibration discrepancy')
    with (ROOT/'data/derived/calibration_small_18um.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    records=[]
    for row in rows:
        matches=[r for r in final['records'] if {'variant':SELECTED,'record_id':row['record_id']} in r['uses']]
        if len(matches)!=1:raise RuntimeError('one baseline field per source record required')
        r=matches[0]
        if sha(ROOT/r['file'])!=r['sha256']:raise RuntimeError('baseline field hash changed')
        records.append({'record_id':row['record_id'],'dose_exact':row['activation_magnitude'],'field':r['file'],'field_sha256':r['sha256'],
            'observed_exact':{k:row[k] for k in OBSERVABLES}})
    p=read(ROOT/records[0]['field'])['parameters']
    if len(records)!=17 or p['bulk']!=5000:raise RuntimeError('full17 fixedpenalty scenario required')
    protocol={'version':1,'scope':'One actual local multivariate tangent calibration STEP, not full material/geometry identified-set inference',
        'sources':hashes(),'baseline_direct_results_sha256':sha(BASE/'direct_results.json'),
        'independent_direct_review_sha256':sha(ROOT/'evidence/joint_direct_independent_review.json'),
        'baseline_selection':'Lowest already verified calibration-only common discrepancy across four prior candidates; no larger outcomes read',
        'baseline_variant':SELECTED,'baseline_gain_exact':best['selected_exact_SURROGATE_gain'],
        'same17_base_states_for_all4_variants':True,'records':records,'parameters':p,'fixed_length_um':str(L0),
        'all_mechanical_tangent_columns':list(COLUMNS),'active_LP_columns':list(COORDS),'active_column_indices':list(ACTIVE),
        'fixed_compression_coefficient':5000,'fixed_geometry':True,'trust_box':BOUNDS,
        'trust_box_scope':'Numerical one-step trust bounds, NOT measured material priors or an exhaustive physical domain',
        'observer_difference_steps':FD_H,'observer_step_shrink_levels':MAX_SHRINK,
        'observer_maximum_interval_tolerance_um':OBSERVER_TOL,'affine_derivative_rule':'Use smaller-step central graph-difference midpoint; retain both steps and interval widths; all derivatives floating diagnostic except separate exact trace-mean formula',
        'mean_crosscheck':'Exact rational physical-node mean derivative intervals contracted with serialized floating displacement tangents',
        'LP_certificate':'Exact rational primal upper / weak-dual lower for defined rounded affine surrogate, gap<=1e-8um',
        'no_new_equilibrium_solves':True,'heldout_response_files_read':False,
        'not_claimed':['valid derivative of nonsmooth maximum','validated root/IFT hypothesis','uniform nonlinear forward bound','new equilibrated material fit','physical prior','full Theta_delta','experimental validation','novelty']}
    save(OUT/'protocol.json',protocol)
    print(json.dumps({'stage':'prepared','source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'records':17}),flush=True)

def assert_protocol(p):
    if p['sources']!=hashes():raise RuntimeError('frozen dependency/source changed')
    if p['baseline_direct_results_sha256']!=sha(BASE/'direct_results.json'):raise RuntimeError('baseline result changed')

def perturb_observers(state,model,u,S,column,h):
    values={};geometry={}
    for sign in (-1,1):
        vector=u+sign*h*S[:,column]
        current=model.coordinates(torch.tensor(vector,dtype=torch.float64))
        cert=model.orientation_certificate(current)
        if not cert['certified_positive']:raise ValueError('perturbed coefficient vector fails local orientation screen')
        values[sign]=graphs(state,current.detach().numpy())
        geometry[str(sign)]={'floating_local_orientation_positive':True,'scope':'Graph perturbation diagnostic, not a stationary or exact-root field'}
    derivatives={}
    for variant in values[1]:
        derivatives[variant]={}
        for key in OBSERVABLES:
            a,b=map(Q,values[1][variant][key]);c,d=map(Q,values[-1][variant][key]);den=2*qfloat(h)
            derivatives[variant][key]=[str((a-d)/den),str((b-c)/den)]
    return {'step':h,'endpoint_intervals':{str(k):v for k,v in values.items()},'central_difference_intervals':derivatives,'geometry':geometry}

def smooth_gradients(state,model,S):
    result={}
    for weight in ('current_disk','current_radial_line'):
        g=mean_state_gradients(state,weighting=weight);out={}
        for name,gradient in g['node_gradients'].items():
            flat=[v for node in gradient for v in node]
            intervals=[]
            for col in ACTIVE:
                lo=hi=Q(0)
                for idx,coefficient in zip(model.free,S[:,col]):
                    v=flat[idx];a=qfloat(coefficient)*L0
                    lo+=min(a*v.lo,a*v.hi);hi+=max(a*v.lo,a*v.hi)
                intervals.append([str(lo),str(hi)])
            out[name]=intervals
        result[weight]=out
    return result

def run():
    protocol=read(OUT/'protocol.json');assert_protocol(protocol)
    if (OUT/'job.json').exists():raise RuntimeError('job record exists; inspect live/terminal state instead of restarting')
    save(OUT/'job.json',{'status':'RUNNING','pid':os.getpid(),'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json')})
    rows=[];started=time.monotonic()
    try:
        for i,rec in enumerate(protocol['records']):
            assert_protocol(protocol);tick=time.monotonic();state=read(ROOT/rec['field'])
            if sha(ROOT/rec['field'])!=rec['field_sha256'] or state['parameters']!=protocol['parameters']:raise RuntimeError('baseline field/parameters changed')
            model,u,H,B,S,verification=mechanical_tangent(state)
            folder=OUT/'records'/f'{i:02d}';folder.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(folder/'mechanical_tangents.npz',free_vector=u,free_indices=model.free,
                H_data=H.data,H_indices=H.indices,H_indptr=H.indptr,H_shape=H.shape,
                parameter_force_derivatives=B,displacement_tangents=S)
            save(folder/'mechanical_verification.json',verification)
            if not all(verification['checks'].values()):raise RuntimeError('mechanical tangent audit failed')
            baseline=graphs(state);directions=[]
            for column in ACTIVE:
                trials=[];pair=None
                for level in range(MAX_SHRINK+1):
                    try:
                        pair=[perturb_observers(state,model,u,S,column,h/(2**level)) for h in FD_H]
                        trials.append({'level':level,'accepted_graph_difference_pair':True});break
                    except (ValueError,RuntimeError) as error:
                        trials.append({'level':level,'accepted_graph_difference_pair':False,'error':str(error)})
                if pair is None:raise RuntimeError('no admissible graph-difference pair for '+COLUMNS[column])
                discrepancy={v:{k:abs(float(sum(map(Q,pair[0]['central_difference_intervals'][v][k]))/2)-float(sum(map(Q,pair[1]['central_difference_intervals'][v][k]))/2)) for k in OBSERVABLES} for v in baseline}
                directions.append({'coordinate':COLUMNS[column],'trials':trials,'two_steps':pair,'step_midpoint_discrepancy_um_per_coordinate':discrepancy})
            analytic=smooth_gradients(state,model,S);mean_checks=[]
            for weight in ('current_disk','current_radial_line'):
                variant=weight+'__regional_mean'
                for key,mean_key in [('region_averaged_shortening_um','height_regional_mean_shortening'),
                    ('apical_max_shortening_um','apical_regional_mean_downward'),
                    ('basal_max_shortening_um','basal_regional_mean_upward')]:
                    for j,direction in enumerate(directions):
                        a,b=map(Q,analytic[weight][mean_key][j]);differences=[]
                        for pair in direction['two_steps']:
                            c,d=map(Q,pair['central_difference_intervals'][variant][key]);mid=(c+d)/2
                            differences.append({'step':pair['step'],'distance_of_FD_midpoint_to_analytic_tangent_interval':float(max(a-mid,mid-b,Q(0))),
                                'FD_interval_width_um_per_coordinate':float(d-c)})
                        mean_checks.append({'weighting':weight,'mean':mean_key,'coordinate':COORDS[j],
                            'analytic_tangent_interval_exact':[str(a),str(b)],'two_step_checks':differences,
                            'scope':'Floating implicit tangent contracted with exact trace derivative; central differences are numerical diagnostics, not derivative error bounds.'})
            row={'index':i,'record_id':rec['record_id'],'input_field':rec['field'],'input_field_sha256':rec['field_sha256'],
                'activation_strength':state['activation_strength'],'baseline_endpoint_intervals':baseline,
                'coordinate_diagnostics':directions,'analytic_smooth_mean_tangent_intervals':analytic,'analytic_mean_FD_crosschecks':mean_checks,
                'mechanical_archive_sha256':sha(folder/'mechanical_tangents.npz'),'mechanical_verification_sha256':sha(folder/'mechanical_verification.json'),
                'seconds':time.monotonic()-tick,'scope':'Floating equilibrium tangent and non-equilibrated graph perturbations only; no material update equilibrium has been solved.'}
            save(folder/'observer_tangents.json',row);rows.append({'record_id':rec['record_id'],'file':(folder/'observer_tangents.json').relative_to(ROOT).as_posix(),'sha256':sha(folder/'observer_tangents.json'),'seconds':row['seconds']})
            save(OUT/'progress.json',{'status':'RUNNING','rows':rows,'completed':len(rows),'total':17,'elapsed_seconds':time.monotonic()-started})
            print(json.dumps({'stage':'cohort_tangent','record':i+1,'of':17,'record_id':rec['record_id'],'seconds':row['seconds']}),flush=True)
        candidates=[]
        for weight,component in VARIANTS:
            name=weight+'__'+component;base=[];observed=[];J=[];labels=[]
            for rec,rowref in zip(protocol['records'],rows):
                row=read(ROOT/rowref['file'])
                for key in OBSERVABLES:
                    a,b=map(Q,row['baseline_endpoint_intervals'][name][key]);base.append(str((a+b)/2));observed.append(rec['observed_exact'][key]);labels.append({'record_id':rec['record_id'],'observable':key})
                    J.append([float(sum(map(Q,d['two_steps'][1]['central_difference_intervals'][name][key]))/2) for d in row['coordinate_diagnostics']])
            lp=affine_minimax(base,J,observed)
            candidates.append({'variant':name,'all68_constraint_order':labels,'affine_LP':lp,
                'proposed_parameter_multipliers':{k:float(np.exp(v)) for k,v in zip(COORDS,lp['step_display']) if k!='apical_tension_additive'},
                'proposed_apical_tension':protocol['parameters']['apical_surface_tension']+lp['step_display'][2],
                'no_direct_new_equilibrium_verification':True})
        assert_protocol(protocol)
        result={'status':'COMPLETE_17_RECORD_LOCAL_AFFINE_DESCENT','protocol_sha256':sha(OUT/'protocol.json'),'source_hashes':hashes(),
            'records':rows,'candidates':candidates,'elapsed_seconds':time.monotonic()-started,
            'scope':protocol['scope'],'original_five_gate_goal_achieved':False,'new_equilibrium_solves':0}
        save(OUT/'results.json',result);save(OUT/'job.json',{'status':'FINISHED','pid':os.getpid(),'records_verified':17,'candidates':4,'result_sha256':sha(OUT/'results.json')})
        print(json.dumps({'stage':'complete','records':17,'candidates':4,'results_sha256':sha(OUT/'results.json')}),flush=True)
    except BaseException as error:
        save(OUT/'job.json',{'status':'STOPPED_WITH_ERROR','pid':os.getpid(),'records_verified':len(rows),'error':str(error),'source_sha256':sha(__file__)})
        raise

if __name__=='__main__':
    torch.set_num_threads(1);parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['test','prepare','run']);args=parser.parse_args()
    {'test':tests,'prepare':prepare,'run':run}[args.stage]()
