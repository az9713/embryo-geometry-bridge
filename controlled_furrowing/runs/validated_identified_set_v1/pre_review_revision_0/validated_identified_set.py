"""Conditional identified-set cell decisions with exact cover bookkeeping.

Point intervals are never treated as uniform forward bounds. IN requires every
record/channel; OUT on a selected branch does not exclude other branches.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,csv,hashlib,json,math
from joint_calibration_contract import OBSERVABLES,load_observations,interval_residual

ROOT=Path(__file__).resolve().parents[1]
QUANTIFIER='for_all_parameters_in_closed_cell_exists_unique_root_in_declared_branch_domain'
STATIONARY_CLOSURE='independent_stationary_root_per_record'
VARIANTS=('static_current_radius_spatial_max_components','static_current_radius_regional_mean_components')

def rational(v):
    if isinstance(v,bool):raise ValueError('Boolean is not an exact physical value')
    if isinstance(v,float) and not math.isfinite(v):raise ValueError('finite values required')
    try:return v if isinstance(v,Q) else Q(str(v))
    except (ValueError,TypeError,ZeroDivisionError) as e:raise ValueError('finite rational required')from e
def interval(v):
    if not isinstance(v,(list,tuple)) or len(v)!=2:raise ValueError('two interval endpoints required')
    a,b=map(rational,v)
    if a>b:raise ValueError('ordered interval required')
    return a,b
def cell(v):
    if not isinstance(v,dict)or not v or any(not isinstance(k,str)or not k for k in v):raise ValueError('named nonempty parameter box required')
    return {k:interval(v[k])for k in sorted(v)}
def serial_cell(v):return {k:list(map(str,pair))for k,pair in cell(v).items()}
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def cell_hash(v):return digest(serial_cell(v))
def subset(inner,outer):
    a,b=cell(inner),cell(outer)
    return set(a)==set(b)and all(b[k][0]<=a[k][0]<=a[k][1]<=b[k][1]for k in a)
def split(v,coordinate,at=None):
    p=cell(v)
    if coordinate not in p:raise ValueError('split coordinate absent')
    a,b=p[coordinate];mid=(a+b)/2 if at is None else rational(at)
    if not a<mid<b:raise ValueError('strictly interior split required')
    left=dict(p);right=dict(p);left[coordinate]=(a,mid);right[coordinate]=(mid,b)
    return serial_cell(left),serial_cell(right)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def read(p):return json.loads(Path(p).read_text())

def validate_observations(observations,discrepancies,doses):
    if not observations or set(doses)!=set(observations):raise ValueError('every source record requires its declared dose')
    if set(discrepancies)!=set(OBSERVABLES):raise ValueError('explicit discrepancy for every endpoint required')
    delta={k:rational(discrepancies[k])for k in OBSERVABLES}
    if min(delta.values())<0:raise ValueError('nonnegative discrepancies required')
    for rid,row in observations.items():
        if set(row)!=set(OBSERVABLES):raise ValueError('all four endpoint channels required for every record')
        if rational(doses[rid])<0:raise ValueError('nonnegative source dose required')
        for value in row.values():rational(value)
    return delta

def static_operator_floor(observations,discrepancies,operator):
    """Necessary constraints for declared STATIC signed graph functionals.

    No mechanical monotonicity or source temporal-peak identity is assumed.
    The component-mean equality is conditional on SAME static state/ROI/reference.
    """
    if operator not in VARIANTS:raise ValueError('explicit supported static operator required')
    delta={k:rational(discrepancies[k])for k in OBSERVABLES};rows=[];floor=Q(0)
    M,S,A,B=OBSERVABLES
    for rid,row in observations.items():
        y={k:rational(row[k])for k in OBSERVABLES}
        candidates=[('mean_not_above_max',y[S]-y[M],delta[S]+delta[M],2)]
        if operator==VARIANTS[0]:
            candidates += [('max_not_above_sum_component_maxima',y[M]-y[A]-y[B],delta[M]+delta[A]+delta[B],3),
                ('mean_not_above_sum_component_maxima',y[S]-y[A]-y[B],delta[S]+delta[A]+delta[B],3)]
        else:
            candidates += [('same_state_mean_identity',abs(y[S]-y[A]-y[B]),delta[S]+delta[A]+delta[B],3),
                ('sum_component_means_not_above_max',y[A]+y[B]-y[M],delta[A]+delta[B]+delta[M],3)]
        for name,lhs,rhs,norm in candidates:
            floor=max(floor,lhs/norm,Q(0));rows.append({'record_id':rid,'necessary_relation':name,'lhs_exact':str(lhs),
                'allowed_error_sum_exact':str(rhs),'violated':lhs>rhs,'common_delta_lower_bound_exact':str(max(Q(0),lhs/norm))})
    return {'operator':operator,'witnesses':rows,'any_violation':any(r['violated']for r in rows),
        'necessary_common_delta_floor_exact':str(floor),'necessary_common_delta_floor_display':float(floor),
        'quantifier':'All parameters and all branches whose outputs obey this declared static graph operator',
        'scope':'Necessary only; not measured noise, a fitted error allowance, source temporal validation, or a full physical identified set.'}

def artifact(reference,base):
    if set(reference)!=set(('file','sha256')):raise ValueError('artifact requires exact file/hash fields')
    p=Path(reference['file']);p=p if p.is_absolute()else Path(base)/p
    if not p.is_file()or sha(p)!=reference['sha256']:raise ValueError('missing/changed proof artifact')
    return read(p),reference['sha256']
def audit_matches(a,certificate_hash,outer_hash,contract_hash):
    return a.get('status')=='PASS' and bool(a.get('checks')) and all(v is True for v in a['checks'].values()) and a.get('certificate_sha256')==certificate_hash and a.get('parameter_cell_sha256')==outer_hash and a.get('model_observation_contract_sha256')==contract_hash and a.get('uniform_quantifier_verified') is True

def uniform_bounds(package,requested_cell,doses,contract_hash,base=ROOT):
    """Accept only hash-bound uniformly reviewed root+observer packages.

    This checks proof provenance/quantifier/data association. The referenced
    independent audits are the trusted proof checkers; JSON status alone is not
    a new mathematical proof. Missing/point-only packages contribute NO bounds.
    """
    bounds={};provenance=[];issues=[];seen=set()
    if package is None:return bounds,provenance,['No uniform forward package supplied']
    if package.get('kind')!='AUDITED_UNIFORM_CELL_FORWARD_V1':return bounds,provenance,['Point/root-only/unrecognized artifact is not a uniform cohort observation package']
    if package.get('model_observation_contract_sha256')!=contract_hash:return bounds,provenance,['Model/geometry/observation contract mismatch']
    if package.get('branch_law_semantics')!=STATIONARY_CLOSURE:return bounds,provenance,['Only independent stationary-root-per-record closure is supported; a common path/stability law needs its own audit adapter']
    for entry in package.get('records',[]):
        rid=entry.get('record_id')
        if not isinstance(rid,str):issues.append({'record_id':None,'reason':'String source record identifier required'});continue
        if rid in seen:
            bounds.pop(rid,None);provenance=[p for p in provenance if p['record_id']!=rid]
            issues.append({'record_id':rid,'reason':'Duplicate source record invalidates that entire record package'});continue
        seen.add(rid)
        try:
            if rid not in doses:raise ValueError('unknown source record')
            root,rh=artifact(entry['root_certificate'],base);rr,rrh=artifact(entry['root_review'],base)
            obs,oh=artifact(entry['observer_certificate'],base);orev,orh=artifact(entry['observer_review'],base)
            outer=root['parameter_cell'];outer_hash=cell_hash(outer)
            if not subset(requested_cell,outer):raise ValueError('requested box not contained in certified parameter cell')
            if root.get('quantifier')!=QUANTIFIER or root.get('status')!='CERTIFIED_UNIFORM_LOCAL_ROOT':raise ValueError('root is not certified uniformly on cell')
            if root.get('branch_scope')!='selected_unique_root_in_declared_domain':raise ValueError('unsupported root branch coverage claim')
            if root.get('input_kind')!='source_record_dose' or root.get('strength_relation')!='common_gain_times_source_dose':raise ValueError('diagnostic strength is not a source-dose binding')
            if root.get('record_id')!=rid or rational(root['dose_exact'])!=rational(doses[rid]):raise ValueError('root source/dose mismatch')
            if not audit_matches(rr,rh,outer_hash,contract_hash):raise ValueError('uniform root audit not bound or not passed')
            if rr.get('branch_law_semantics')!=STATIONARY_CLOSURE:raise ValueError('root review does not bind the declared stationary closure')
            if obs.get('status')!='CERTIFIED_UNIFORM_CELL_OBSERVABLES' or obs.get('root_certificate_sha256')!=rh or cell_hash(obs['parameter_cell'])!=outer_hash:raise ValueError('observer lacks same-cell uniform root binding')
            if obs.get('record_id')!=rid or rational(obs['dose_exact'])!=rational(doses[rid]):raise ValueError('observer source/dose mismatch')
            if obs.get('branch_domain_id')!=root.get('branch_domain_id') or not root.get('branch_domain_id'):raise ValueError('branch-domain mismatch')
            if not audit_matches(orev,oh,outer_hash,contract_hash) or orev.get('forward_family_admissible') is not True:raise ValueError('uniform observer/admissibility audit not bound or not passed')
            if orev.get('branch_law_semantics')!=STATIONARY_CLOSURE:raise ValueError('observer review does not bind the declared stationary closure')
            if not obs.get('endpoint_intervals') or not set(obs['endpoint_intervals'])<=set(OBSERVABLES):raise ValueError('observer must supply a nonempty subset of declared channels')
            if rid in bounds:raise ValueError('duplicate source record package; no silent selection')
            bounds[rid]={k:interval(v)for k,v in obs['endpoint_intervals'].items()}
            provenance.append({'record_id':rid,'root_sha256':rh,'root_review_sha256':rrh,'observer_sha256':oh,'observer_review_sha256':orh,
                'outer_cell_sha256':outer_hash,'requested_subset_exactly_checked':True,'branch_domain_id':root['branch_domain_id'],
                'channels_covered':sorted(bounds[rid])})
        except (KeyError,ValueError,TypeError)as error:
            # Ambiguous duplicate/invalid entries never select an earlier value.
            if rid in bounds:bounds.pop(rid);provenance=[p for p in provenance if p['record_id']!=rid]
            issues.append({'record_id':rid,'reason':str(error)})
    return bounds,provenance,issues

def classify_cell(observations,discrepancies,doses,parameter_cell,*,contract_hash,package=None,operator=None,base=ROOT):
    delta=validate_observations(observations,discrepancies,doses);parameter_cell=serial_cell(parameter_cell)
    bounds,proofs,issues=uniform_bounds(package,parameter_cell,doses,contract_hash,base);rows=[]
    for rid,obs in observations.items():
        for key in OBSERVABLES:
            y=rational(obs[key]);row={'record_id':rid,'observable':key,'observed_exact':str(y),
                'allowed_band_exact':[str(y-delta[key]),str(y+delta[key])],'status':'UNRESOLVED','reason':'No audited uniform forward interval for this source row'}
            if rid in bounds and key in bounds[rid]:
                lo,hi=bounds[rid][key];minimum,maximum=interval_residual(y,(lo,hi))
                row.update({'prediction_interval_exact':[str(lo),str(hi)],'minimum_residual_exact':str(minimum),
                    'maximum_residual_exact':str(maximum),'status':'OUT'if minimum>delta[key]else'IN'if maximum<=delta[key]else'UNRESOLVED',
                    'reason':'Interval comparison for every parameter in requested cell on specified certified root branch'})
            rows.append(row)
    selected='OUT'if any(r['status']=='OUT'for r in rows)else'IN'if all(r['status']=='IN'for r in rows)else'UNRESOLVED'
    necessary=static_operator_floor(observations,delta,operator)if operator is not None else None
    universal_out=necessary is not None and necessary['any_violation']
    # OUT on one selected root family cannot exclude a different admissible root.
    all_branches='OUT'if universal_out else'IN_STATIONARY_ROOT_TUPLE_WITNESS'if selected=='IN'else'UNRESOLVED'
    context={'observations':{rid:{k:str(rational(v))for k,v in row.items()}for rid,row in observations.items()},
        'doses':{rid:str(rational(v))for rid,v in doses.items()},'discrepancies':{k:str(v)for k,v in delta.items()},
        'operator':operator,'model_observation_contract_sha256':contract_hash}
    return {'status':'OUT'if universal_out else selected,'selected_branch_status':selected,
        'all_branches_existential_parameter_status':all_branches,'parameter_cell':parameter_cell,'parameter_cell_sha256':cell_hash(parameter_cell),
        'model_observation_contract_sha256':contract_hash,'record_count':len(observations),'constraint_count':len(rows),'constraints':rows,
        'uniform_records_covered':sorted(bounds),'proof_provenance':proofs,'unresolved_proof_issues':issues,'conditional_static_operator_necessary_test':necessary,
        'uniform_constraints_covered':[{'record_id':rid,'observable':k}for rid in sorted(bounds)for k in sorted(bounds[rid])],
        'inference_context':context,'inference_context_sha256':digest(context),'uniform_forward_package':package,
        'branch_law_semantics':STATIONARY_CLOSURE,'common_path_or_energetic_stability_certified':False,
        'discrepancy_vector_exact':{k:str(v)for k,v in delta.items()},
        'scope':'Selected-branch cell decisions only unless an explicit universal static-operator contradiction applies. IN requires every row; OUT may use one uniform contradiction. Neither incomplete coverage nor a finite box exhausts the physical domain.'}

class CoverLedger:
    """Exact finite-box subdivision with explicit ambient-domain remainder."""
    def __init__(self,root_cell,ambient_domain,*,uncovered_physics):
        self.root=serial_cell(root_cell);self.domain=ambient_domain
        if set(self.root)!=set(ambient_domain):raise ValueError('ambient coordinate names must match box')
        for k,(a,b)in cell(self.root).items():
            d=ambient_domain[k];lo=None if d.get('lower')is None else rational(d['lower']);hi=None if d.get('upper')is None else rational(d['upper'])
            if lo is not None and(a<lo or(a==lo and d.get('lower_strict',False))):raise ValueError('root box outside ambient lower domain')
            if hi is not None and(b>hi or(b==hi and d.get('upper_strict',False))):raise ValueError('root box outside ambient upper domain')
        self.nodes={'0':{'cell':self.root,'parent':None,'children':[],'classification':{'status':'UNRESOLVED','reason':'No uniform cell proof attached'}}}
        equal=all(d.get('lower')is not None and d.get('upper')is not None and not d.get('lower_strict',False)and not d.get('upper_strict',False)and interval((d['lower'],d['upper']))==cell(self.root)[k]for k,d in ambient_domain.items())
        self.remainder={'status':'EMPTY'if equal else'UNCOVERED','set_definition':'ambient_domain minus closed root_cell','ambient_domain':ambient_domain,'root_cell':self.root,'unbounded_tail_proof':None}
        self.uncovered_physics=list(uncovered_physics);self.inference_context_sha256=None
    def subdivide(self,node,coordinate,at=None):
        old=self.nodes[node]
        if old['children']:raise ValueError('only a current leaf can be subdivided')
        children=split(old['cell'],coordinate,at);ids=[node+'0',node+'1']
        for name,box in zip(ids,children):self.nodes[name]={'cell':box,'parent':node,'children':[],'classification':{'status':'UNRESOLVED','reason':'New leaf retained pending proof'}}
        old['children']=ids;old['split_coordinate']=coordinate;return ids
    def attach(self,node,classification):
        if self.nodes[node]['children']:raise ValueError('classify current leaves only')
        if cell_hash(classification['parameter_cell'])!=cell_hash(self.nodes[node]['cell']):raise ValueError('decision bound to a different parameter cell')
        if classification.get('status')not in ('IN','OUT','UNRESOLVED') or classification.get('inference_context_sha256')!=digest(classification.get('inference_context')):raise ValueError('unbound/invalid cell decision')
        if self.inference_context_sha256 is None:self.inference_context_sha256=classification['inference_context_sha256']
        if self.inference_context_sha256!=classification['inference_context_sha256']:raise ValueError('cannot mix data, discrepancy vectors or observation contracts in one cover')
        self.nodes[node]['classification']=classification
    def verify(self):
        visited=set()
        def walk(name):
            if name in visited:raise ValueError('cycle/shared subtree in cover')
            visited.add(name);n=self.nodes[name]
            if n['children']:
                if len(n['children'])!=2:raise ValueError('binary subdivision required')
                left,right=(self.nodes[x]for x in n['children']);coordinate=n['split_coordinate']
                a,b=cell(n['cell'])[coordinate];mid=cell(left['cell'])[coordinate][1]
                expected=split(n['cell'],coordinate,mid)
                if expected!=(serial_cell(left['cell']),serial_cell(right['cell'])):raise ValueError('children fail exact parent-cover identity')
                for child in n['children']:
                    if self.nodes[child]['parent']!=name:raise ValueError('child provenance mismatch')
                    walk(child)
        walk('0')
        if visited!=set(self.nodes):raise ValueError('orphan/dropped cells in cover ledger')
        leaves=[name for name,n in self.nodes.items()if not n['children']]
        for name in leaves:
            decision=self.nodes[name]['classification']
            if decision['status']not in ('IN','OUT','UNRESOLVED'):raise ValueError('unknown leaf classification')
            if 'inference_context'in decision:
                c=decision['inference_context']
                if decision['inference_context_sha256']!=self.inference_context_sha256:raise ValueError('mixed inference contexts')
                replay=classify_cell(c['observations'],c['discrepancies'],c['doses'],self.nodes[name]['cell'],
                    contract_hash=c['model_observation_contract_sha256'],package=decision['uniform_forward_package'],operator=c['operator'])
                if digest(replay)!=digest(decision):raise ValueError('leaf proof/decision replay changed')
        unresolved=[n for n in leaves if self.nodes[n]['classification']['status']=='UNRESOLVED']
        return {'exact_closed_leaf_union_equals_root_box':True,'leaves':leaves,
            'unresolved_leaves':unresolved,'finite_selected_branch_cover_complete':not unresolved,
            'ambient_remainder_status':self.remainder['status'],'uncovered_physics':self.uncovered_physics,
            'global_branch_path_certified':False,'full_physical_identified_set_complete':False,
            'scope':'Cover geometry is proved by exact binary union induction. Scientific cell decisions require their own proof receipts; unbounded tails and excluded physics are retained.'}
    def export(self):return {'root_cell':self.root,'ambient_domain':self.domain,'nodes':self.nodes,'outside_remainder':self.remainder,'uncovered_physics':self.uncovered_physics,'inference_context_sha256':self.inference_context_sha256,'verification':self.verify()}

def synthetic_package(directory,parameter_cell,observations,doses,affine_maps,contract):
    """Exact affine synthetic forward theorem; NEVER labelled mechanical evidence."""
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);entries=[];C=cell(parameter_cell);ch=cell_hash(C)
    def put(name,value):
        p=directory/name;save(p,value);return {'file':p.resolve().relative_to(ROOT.resolve()).as_posix(),'sha256':sha(p)}
    for rid in observations:
        root={'status':'CERTIFIED_UNIFORM_LOCAL_ROOT','quantifier':QUANTIFIER,'parameter_cell':serial_cell(C),
            'branch_scope':'selected_unique_root_in_declared_domain','branch_domain_id':'synthetic_affine_root',
            'input_kind':'source_record_dose','record_id':rid,'dose_exact':str(doses[rid]),'strength_relation':'common_gain_times_source_dose',
            'scope':'SYNTHETIC theorem F(u,theta)=u-affine(theta)=0; not a mechanical root'}
        rr=put(rid+'_root.json',root)
        review=lambda h:{'status':'PASS','checks':{'exact_affine_root_identity':True},'certificate_sha256':h,'parameter_cell_sha256':ch,'model_observation_contract_sha256':contract,'uniform_quantifier_verified':True,'forward_family_admissible':True,'branch_law_semantics':STATIONARY_CLOSURE,'scope':'SYNTHETIC exact affine theorem only'}
        rootreview=put(rid+'_root_review.json',review(rr['sha256']));values={}
        for key,(constant,coeffs)in affine_maps[rid].items():
            lo=hi=rational(constant)
            for coordinate,a in coeffs.items():
                a=rational(a);x,y=C[coordinate];lo+=min(a*x,a*y);hi+=max(a*x,a*y)
            values[key]=[str(lo),str(hi)]
        obs={'status':'CERTIFIED_UNIFORM_CELL_OBSERVABLES','root_certificate_sha256':rr['sha256'],'parameter_cell':serial_cell(C),
            'record_id':rid,'dose_exact':str(doses[rid]),'branch_domain_id':'synthetic_affine_root','endpoint_intervals':values,
            'scope':'SYNTHETIC exact affine extrema on box only'}
        oo=put(rid+'_observer.json',obs);orv=put(rid+'_observer_review.json',review(oo['sha256']))
        entries.append({'record_id':rid,'root_certificate':rr,'root_review':rootreview,'observer_certificate':oo,'observer_review':orv})
    return {'kind':'AUDITED_UNIFORM_CELL_FORWARD_V1','model_observation_contract_sha256':contract,'branch_law_semantics':STATIONARY_CLOSURE,'records':entries,'scope':'SYNTHETIC exact affine fixture, not biological/mechanical validation'}

def tests():
    out=ROOT/'runs/validated_identified_set_v1/synthetic';out.mkdir(parents=True,exist_ok=True)
    obs={'r':{k:Q(0)for k in OBSERVABLES}};doses={'r':Q(1)};delta={k:Q(1,2)for k in OBSERVABLES};contract=digest({'synthetic':'u=theta'})
    maps={'r':{k:(0,{'x':1})for k in OBSERVABLES}};box={'x':[-1,1]}
    package=synthetic_package(out/'proofs',box,obs,doses,maps,contract)
    checks={};whole=classify_cell(obs,delta,doses,box,contract_hash=contract,package=package)
    checks['whole_cell_not_mistaken_for_point_IN']=whole['status']=='UNRESOLVED'
    inside=classify_cell(obs,delta,doses,{'x':['-1/2','1/2']},contract_hash=contract,
        package=synthetic_package(out/'inside',{'x':['-1/2','1/2']},obs,doses,maps,contract))
    checks['uniform_all_rows_IN']=inside['status']=='IN' and inside['all_branches_existential_parameter_status']=='IN_STATIONARY_ROOT_TUPLE_WITNESS'
    outside=classify_cell(obs,delta,doses,{'x':['3/4',1]},contract_hash=contract,
        package=synthetic_package(out/'outside',{'x':['3/4',1]},obs,doses,maps,contract))
    checks['selected_branch_OUT_not_all_branches_OUT']=outside['status']=='OUT'and outside['all_branches_existential_parameter_status']=='UNRESOLVED'
    missing=classify_cell(obs,delta,doses,box,contract_hash=contract,package={'kind':'CERTIFIED_LOCAL_ROOT_BINARY_QUADRATURE'})
    checks['actual_point_kind_not_promoted_to_cell']=missing['status']=='UNRESOLVED' and not missing['uniform_records_covered']
    bad=json.loads(json.dumps(package));bad['records'][0]['root_certificate']['sha256']='0'*64
    checks['mutated_certificate_hash_rejected']=classify_cell(obs,delta,doses,box,contract_hash=contract,package=bad)['status']=='UNRESOLVED'
    partial=dict(obs);partial['missing']={k:Q(0)for k in OBSERVABLES};pd={**doses,'missing':1}
    checks['missing_record_retained_not_dropped']=classify_cell(partial,delta,pd,box,contract_hash=contract,package=package)['constraint_count']==8
    one=synthetic_package(out/'onechannel',{'x':['3/4',1]},obs,doses,{'r':{OBSERVABLES[0]:(0,{'x':1})}},contract)
    one_result=classify_cell(obs,delta,doses,{'x':['3/4',1]},contract_hash=contract,package=one)
    checks['one_uniform_exclusion_keeps_other3_rows_missing']=one_result['status']=='OUT' and sum(r['status']=='UNRESOLVED'for r in one_result['constraints'])==3
    checks['wrong_contract_not_accepted']=classify_cell(obs,delta,doses,box,contract_hash='wrong',package=package)['status']=='UNRESOLVED'
    duplicate=json.loads(json.dumps(package));duplicate['records']*=2
    checks['duplicate_record_not_silently_selected']=not classify_cell(obs,delta,doses,box,contract_hash=contract,package=duplicate)['uniform_records_covered']
    # Exact full cover with two OUT leaves and an IN middle box, no gap.
    ledger=CoverLedger(box,{'x':{'lower':-1,'upper':1}},uncovered_physics=['Synthetic only; not physical material domain'])
    a,b=ledger.subdivide('0','x','-1/2');c,d=ledger.subdivide(b,'x','1/2')
    for name in (a,c,d):
        C=ledger.nodes[name]['cell'];pkg=synthetic_package(out/name,C,obs,doses,maps,contract)
        decision=classify_cell(obs,delta,doses,C,contract_hash=contract,package=pkg);ledger.attach(name,decision)
    checks['exact_closed_cover_retains_boundary_UNRESOLVED']=ledger.verify()['exact_closed_leaf_union_equals_root_box'] and len(ledger.verify()['unresolved_leaves'])==2
    # Closed boxes touching the allowed boundary contain an exact compatible point;
    # a single mixed box must remain unresolved rather than falsely OUT.
    checks['boundary_cells_not_falsely_excluded']=all(ledger.nodes[n]['classification']['status']=='UNRESOLVED'for n in (a,d))
    full=CoverLedger(box,{'x':{'lower':-1,'upper':1}},uncovered_physics=['Synthetic model only'])
    children=full.subdivide('0','x')
    for name in children:
        C=full.nodes[name]['cell'];pkg=synthetic_package(out/('full'+name),C,obs,doses,maps,contract)
        full.attach(name,classify_cell(obs,{k:1 for k in OBSERVABLES},doses,C,contract_hash=contract,package=pkg))
    checks['fully_resolved_analytic_finite_cover']=full.verify()['finite_selected_branch_cover_complete'] and full.verify()['ambient_remainder_status']=='EMPTY' and all(full.nodes[n]['classification']['status']=='IN'for n in children)
    changed=classify_cell(obs,{k:2 for k in OBSERVABLES},doses,full.nodes[children[0]]['cell'],contract_hash=contract)
    try:full.attach(children[0],changed);checks['mixed_discrepancy_context_rejected']=False
    except ValueError:checks['mixed_discrepancy_context_rejected']=True
    old_decision=ledger.nodes[a]['classification'];ledger.nodes[a]['classification']=json.loads(json.dumps(old_decision));ledger.nodes[a]['classification']['status']='IN'
    try:ledger.verify();checks['forged_leaf_status_fails_proof_replay']=False
    except ValueError:checks['forged_leaf_status_fails_proof_replay']=True
    ledger.nodes[a]['classification']=old_decision
    unbounded=CoverLedger(box,{'x':{'lower':None,'upper':None}},uncovered_physics=['All other material/geometry coordinates'])
    unbounded.subdivide('0','x');checks['unbounded_tails_preserved']=unbounded.verify()['ambient_remainder_status']=='UNCOVERED'
    damaged=json.loads(json.dumps(ledger.export()));old=ledger.nodes[a]['cell'];ledger.nodes[a]['cell']={'x':[-1,'-3/4']}
    try:ledger.verify();checks['gap_in_cover_rejected']=False
    except ValueError:checks['gap_in_cover_rejected']=True
    ledger.nodes[a]['cell']=old
    relation_obs={'r':dict(zip(OBSERVABLES,[5,1,4,3]))};floor=static_operator_floor(relation_obs,{k:0 for k in OBSERVABLES},VARIANTS[1])
    checks['exact_mean_identity_necessary_floor']=floor['any_violation'] and Q(floor['necessary_common_delta_floor_exact'])==2
    result={'status':'PASS'if all(checks.values())else'FAIL','checks':checks,'source_sha256':sha(__file__),
        'finite_cover':ledger.export(),'fully_resolved_synthetic_cover':full.export(),'unbounded_cover':unbounded.export(),'uniform_IN_example':inside,'selected_OUT_example':outside,
        'scope':'Analytic synthetic proof/plumbing only. Not an actual mechanical forward cell or completed identified set.'}
    result['dependency_source_hashes']={'joint_calibration_contract.py':sha(Path(__file__).with_name('joint_calibration_contract.py'))}
    save(ROOT/'evidence/validated_identified_set_synthetic.json',result)
    print(json.dumps({'status':result['status'],'checks':len(checks)}),flush=True)
    if not all(checks.values()):raise RuntimeError('identified-set synthetic checks failed')

def actual():
    test=read(ROOT/'evidence/validated_identified_set_synthetic.json')
    if test['status']!='PASS'or test['source_sha256']!=sha(__file__):raise RuntimeError('matching synthetic checks required')
    observations=load_observations()
    with (ROOT/'data/derived/calibration_small_18um.csv').open(newline='',encoding='utf-8-sig')as f:source=list(csv.DictReader(f))
    doses={r['record_id']:r['activation_magnitude']for r in source}
    C={'apical_area_modulus':['1/10','3/10'],'basal_area_modulus':[1,3],'apical_surface_tension':[0,'1/50'],
        'basal_surface_tension':[1,3],'common_gain':[1,10]}
    # The physical quotient permits zero surface/gain coefficients. Log-step
    # coordinates used locally do not remove those boundary faces globally.
    domain={k:{'lower':0,'lower_strict':False,'upper':None}for k in C}
    contract=digest({'scope':'Conditional static flat current-radius operator/source-column pairing; not source-validated temporal peaks','fixed_mu':1,'fixed_kappa':5000,'fixed_geometry':'33.5um flat small-region reference'})
    profiles=[]
    for operator in VARIANTS:
        for value in (0,1,2,3,4,5,6):
            delta={k:value for k in OBSERVABLES}
            profiles.append({'operator':operator,'common_delta_exact':str(value),'decision':classify_cell(observations,delta,doses,C,contract_hash=contract,operator=operator)})
    ledger=CoverLedger(C,domain,uncovered_physics=['Compression coefficient/shear ratio beyond fixed5000','Geometry, source/reference/temporal operator uncertainty','Other root branches'])
    ledger.subdivide('0','common_gain')
    # Genuine existing root certificate has no all-parameter/all-source envelope.
    point=read(ROOT/'evidence/algebraic_root_certification_feasibility.json')
    zero_delta={k:0 for k in OBSERVABLES};point_rejection=classify_cell(observations,zero_delta,doses,C,contract_hash=contract,package=point)
    checks={'all17_records_68_entries_retained':all(v['decision']['record_count']==17 and v['decision']['constraint_count']==68 for v in profiles),
        'point_root_not_promoted_to_uniform_cohort':point_rejection['status']=='UNRESOLVED'and not point_rejection['uniform_records_covered'],
        'all_actual_numerical_forward_cells_still_uncovered':all(not v['decision']['uniform_records_covered']for v in profiles),
        'unbounded_tails_retained':ledger.verify()['ambient_remainder_status']=='UNCOVERED',
        'physical_identified_set_not_complete':not ledger.verify()['full_physical_identified_set_complete']}
    receipt={'status':'PASS'if all(checks.values())else'FAIL','checks':checks,'source_sha256':sha(__file__),
        'dependency_source_hashes':{'joint_calibration_contract.py':sha(Path(__file__).with_name('joint_calibration_contract.py'))},
        'calibration_csv_sha256':sha(ROOT/'data/derived/calibration_small_18um.csv'),'actual_existing_point_certificate_sha256':sha(ROOT/'evidence/algebraic_root_certification_feasibility.json'),
        'conditional_tolerance_family':profiles,'cover_ledger':ledger.export(),'point_as_cell_rejection':point_rejection,
        'operator_floor_summary':[static_operator_floor(observations,zero_delta,op)for op in VARIANTS],
        'no_new_equilibrium_solves':True,'full_physical_identified_set_complete':False,
        'scope':'Actual-data necessary exclusions are conditional static graph identities only. No actual all17 uniform mechanical/observer cell has been certified; unresolved source assumptions and physical tails remain.'}
    save(ROOT/'evidence/validated_identified_set_actual.json',receipt)
    print(json.dumps({'status':receipt['status'],'records':17,'constraints':68,
        'conditional_common_delta_floors':{x['operator']:x['necessary_common_delta_floor_display']for x in receipt['operator_floor_summary']},
        'full_identified_set_complete':False}),flush=True)
    if not all(checks.values()):raise RuntimeError('actual inference-layer verification failed')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['test','actual']);args=parser.parse_args()
    {'test':tests,'actual':actual}[args.stage]()
