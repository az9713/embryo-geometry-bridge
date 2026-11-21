"""Independent analytic/adversarial review of the identified-set layer.

Own fixtures only; never calls producer test/actual stages or physical solvers.
Exact affine output extrema are enumerated at vertices independently.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,copy,csv,itertools,json,hashlib,time
import validated_identified_set as engine

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/review/identified_set'
RECEIPT=ROOT/'evidence/validated_identified_set_independent_review.json'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def write(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf8')
def ref(p):return {'file':p.relative_to(ROOT).as_posix(),'sha256':sha(p)}


def fixture(name,box,observations,doses,affine,contract,channels=None):
    """F(u,theta)=u-a(theta); exact unique root, Cartesian channels per record."""
    folder=OUT/name;entries=[]
    coords=sorted(box);vertices=list(itertools.product(*[(Q(box[k][0]),Q(box[k][1])) for k in coords]))
    def put(filename,value):
        p=folder/filename;write(p,value);return ref(p)
    for rid in observations:
        root={'status':'CERTIFIED_UNIFORM_LOCAL_ROOT','quantifier':engine.QUANTIFIER,'parameter_cell':engine.serial_cell(box),
            'branch_scope':'selected_unique_root_in_declared_domain','branch_domain_id':'independent_affine',
            'input_kind':'source_record_dose','record_id':rid,'dose_exact':str(doses[rid]),'strength_relation':'common_gain_times_source_dose',
            'scope':'Independent SYNTHETIC exact equation u-affine(theta)=0, never biological evidence'}
        rp=put(rid+'_root.json',root)
        def audit(h):return {'status':'PASS','checks':{'exact_affine_equation':True,'vertex_extrema_exact':True},'certificate_sha256':h,
            'parameter_cell_sha256':engine.cell_hash(box),'model_observation_contract_sha256':contract,'uniform_quantifier_verified':True,
            'forward_family_admissible':True,'branch_law_semantics':engine.STATIONARY_CLOSURE}
        rr=put(rid+'_root_review.json',audit(rp['sha256']));values={}
        for key in observations[rid]:
            if channels is not None and (rid,key) not in channels:continue
            constant,coeff=affine[rid][key]
            vals=[Q(constant)+sum((Q(coeff.get(k,0))*v for k,v in zip(coords,vertex)),Q(0)) for vertex in vertices]
            values[key]=list(map(str,(min(vals),max(vals))))
        op=put(rid+'_observer.json',{'status':'CERTIFIED_UNIFORM_CELL_OBSERVABLES','root_certificate_sha256':rp['sha256'],'parameter_cell':engine.serial_cell(box),
            'record_id':rid,'dose_exact':str(doses[rid]),'branch_domain_id':'independent_affine','endpoint_intervals':values})
        ora=put(rid+'_observer_review.json',audit(op['sha256']))
        entries.append({'record_id':rid,'root_certificate':rp,'root_review':rr,'observer_certificate':op,'observer_review':ora})
    return {'kind':'AUDITED_UNIFORM_CELL_FORWARD_V1','model_observation_contract_sha256':contract,'branch_law_semantics':engine.STATIONARY_CLOSURE,'records':entries}


def change_artifact(package,section,mutator,name):
    out=copy.deepcopy(package);old=ROOT/out['records'][0][section]['file'];value=read(old);mutator(value)
    new=OUT/'mutations'/(name+'.json');write(new,value);out['records'][0][section]=ref(new);return out


def must_reject(call):
    try:call();return False
    except (ValueError,KeyError,TypeError):return True


def main():
    if RECEIPT.exists() or OUT.exists():raise RuntimeError('retain prior independent fixtures')
    started=time.time();checks={};rows={};K=engine.OBSERVABLES;contract=engine.digest({'synthetic':'independent affine root and Cartesian endpoint extrema'})
    obs={f'r{i:02}':{k:Q(0) for k in K} for i in range(17)};doses={k:Q(i+1,19) for i,k in enumerate(obs)}
    box={'x':['-1/10','1/10']};affine={rid:{k:(0,{'x':1}) for k in K} for rid in obs};delta={k:Q(1,10) for k in K}
    package=fixture('full68',box,obs,doses,affine,contract)
    call=lambda pkg=package,C=box,ys=obs,ds=doses,dd=delta,op=None:engine.classify_cell(ys,dd,ds,C,contract_hash=contract,package=pkg,operator=op)
    full=call();rows['full68']=full
    checks['all68_IN_only_with_all68_uniform_constraints']=full['status']=='IN' and full['constraint_count']==68 and len(full['uniform_constraints_covered'])==68
    checks['IN_is_stationary_tuple_not_common_stable_path']=full['all_branches_existential_parameter_status']=='IN_STATIONARY_ROOT_TUPLE_WITNESS' and full['common_path_or_energetic_stability_certified'] is False
    noone=copy.deepcopy(package);noone['records'].pop()
    missing=call(noone);checks['missing_record_blocks_IN_and_keeps68rows']=missing['status']=='UNRESOLVED' and missing['constraint_count']==68 and sum(r['status']=='UNRESOLVED' for r in missing['constraints'])==4
    almost=fixture('only67',box,obs,doses,affine,contract,channels={(rid,k) for rid in obs for k in K}-{('r16',K[-1])})
    missing=call(almost);checks['one_missing_channel_blocks_IN']=missing['status']=='UNRESOLVED' and len(missing['uniform_constraints_covered'])==67
    outbox={'x':[1,2]};one=fixture('oneOUT',outbox,obs,doses,affine,contract,channels={('r00',K[2])})
    oneout=call(one,outbox);rows['oneOUT']=oneout
    checks['one_uniform_channel_OUT_with_other67missing']=oneout['status']=='OUT' and sum(r['status']=='UNRESOLVED' for r in oneout['constraints'])==67
    checks['selected_OUT_does_not_exclude_other_branches']=oneout['all_branches_existential_parameter_status']=='UNRESOLVED'
    checks['point_only_package_never_uniform']=not call({'kind':'CERTIFIED_LOCAL_ROOT_BINARY_QUADRATURE'})['uniform_records_covered']
    touch={'x':['1/10','1/5']};tp=fixture('touching',touch,obs,doses,affine,contract);touchr=call(tp,touch)
    checks['closed_band_touch_is_not_false_OUT']=touchr['status']=='UNRESOLVED'
    point={'x':['1/10','1/10']};pp=fixture('boundarypoint',point,obs,doses,affine,contract)
    checks['inclusive_threshold_exact_point_IN']=call(pp,point)['status']=='IN'
    perm=copy.deepcopy(package);perm['records'].reverse();pr=call(perm,ys=dict(reversed(list(obs.items()))),ds=dict(reversed(list(doses.items()))))
    normalize=lambda x:sorted((r['record_id'],r['observable'],r['status']) for r in x['constraints'])
    checks['record_and_package_permutation_invariant']=normalize(pr)==normalize(full) and pr['status']==full['status'] and pr['inference_context_sha256']==full['inference_context_sha256']
    duplicate=copy.deepcopy(package);duplicate['records']=[package['records'][0]]*3+package['records'][1:]
    dr=call(duplicate);checks['triple_duplicate_invalidates_entire_record']=dr['status']=='UNRESOLVED' and 'r00' not in dr['uniform_records_covered'] and len(dr['uniform_records_covered'])==16
    badhash=copy.deepcopy(package);badhash['records'][0]['root_certificate']['sha256']='0'*64
    checks['certificate_corruption_blocks_IN']=call(badhash)['status']=='UNRESOLVED'
    mutations=[('root_review','emptychecks',lambda x:x.update(checks={})),('root_review','nonboolcheck',lambda x:x.update(checks={'bad':1})),
        ('root_review','failedreview',lambda x:x.update(status='FAIL')),('root_review','wrongcell',lambda x:x.update(parameter_cell_sha256='0'*64)),
        ('observer_review','wrongcontract',lambda x:x.update(model_observation_contract_sha256='0'*64)),('observer_review','noadmissibility',lambda x:x.update(forward_family_admissible=False)),
        ('root_review','nouniformaudit',lambda x:x.update(uniform_quantifier_verified=False)),('root_certificate','wrongrootquantifier',lambda x:x.update(quantifier='exists_one_parameter')),
        ('root_certificate','wrongdose',lambda x:x.update(dose_exact='3')),('root_certificate','alldiagnostic',lambda x:x.update(input_kind='diagnostic')),
        ('root_certificate','unsupportedbranch',lambda x:x.update(branch_scope='all_branches')),('observer_certificate','wrongbranch',lambda x:x.update(branch_domain_id='other')),
        ('observer_certificate','wrongobserverdose',lambda x:x.update(dose_exact='3')),('observer_certificate','unknownchannel',lambda x:x['endpoint_intervals'].update(unknown=[0,1]))]
    for section,name,mut in mutations:
        bad=change_artifact(package,section,mut,name);decision=call(bad);checks['reject_'+name]=decision['status']=='UNRESOLVED' and 'r00' not in decision['uniform_records_covered']
    badlaw=copy.deepcopy(package);badlaw['branch_law_semantics']='common_stable_path';checks['reject_unproved_common_path_law']=not call(badlaw)['uniform_records_covered']
    checks['requested_outside_certified_cell_rejected']=not call(package,{'x':[-1,1]})['uniform_records_covered']
    for name,bad in [('bool',True),('nan',float('nan')),('inf',float('inf')),('invalid','NaN')]:checks['reject_invalid_rational_'+name]=must_reject(lambda b=bad:engine.rational(b))
    checks['reverse_interval_rejected']=must_reject(lambda:engine.interval([1,0]))
    # Independent two-dimensional vertex extrema and sharp residual calculation.
    two={'a':['-1/5','3/10'],'b':['1/10','2/5']};amap={rid:{k:(Q(j-2,7),{'a':Q(i+1,13),'b':Q((-1)**j*(j+1),11)}) for j,k in enumerate(K)} for i,rid in enumerate(obs)}
    twop=fixture('two_dimensional',two,obs,doses,amap,contract);twor=call(twop,two,dd={k:Q(3) for k in K});exact=True
    for row in twor['constraints']:
        c,co=amap[row['record_id']][row['observable']];vals=[c+co['a']*a+co['b']*b for a,b in itertools.product(map(Q,two['a']),map(Q,two['b']))];lo,hi=min(vals),max(vals)
        exact &= tuple(map(Q,row['prediction_interval_exact']))==(lo,hi) and Q(row['minimum_residual_exact'])==max(lo,-hi,Q(0)) and Q(row['maximum_residual_exact'])==max(abs(lo),abs(hi))
    checks['all68_two_dimensional_vertex_and_residual_arithmetic_exact']=exact
    # Declared STATIC operator identities, including signed components.
    relobs={'r':dict(zip(K,[5,1,4,3]))};zero={k:0 for k in K};rd={'r':1}
    spatial=engine.static_operator_floor(relobs,zero,engine.VARIANTS[0]);mean=engine.static_operator_floor(relobs,zero,engine.VARIANTS[1])
    checks['mean_identity_is_conditional_not_spatialmax_identity']=not spatial['any_violation'] and mean['any_violation'] and Q(mean['necessary_common_delta_floor_exact'])==2
    unknown=engine.classify_cell(relobs,zero,rd,box,contract_hash=contract,operator=engine.VARIANTS[1]);checks['static_universal_OUT_distinct_from_selected_roots']=unknown['status']=='OUT' and unknown['selected_branch_status']=='UNRESOLVED' and unknown['all_branches_existential_parameter_status']=='OUT'
    checks['unsupported_temporal_operator_rejected']=must_reject(lambda:engine.static_operator_floor(relobs,zero,'mean_of_separate_temporal_peaks'))
    # Cover and standalone replay attacks, using our own proof files only.
    def ledger():
        q=engine.CoverLedger(box,{'x':{'lower':'-1/10','upper':'1/10'}},uncovered_physics=['Synthetic only'])
        q.attach('0',full);return q
    good=ledger();export=good.export();write(OUT/'good_export.json',export)
    checks['standalone_valid_export_replays']=engine.verify_cover_snapshot(read(OUT/'good_export.json'))==export['verification']
    for status in ('IN','OUT'):
        q=ledger();q.nodes['0']['classification']={'status':status};checks['proofless_'+status+'_rejected']=must_reject(q.verify)
    q=ledger();q.nodes['0']['classification']={'status':'IN','inference_context':{}};checks['partial_context_resolved_leaf_rejected']=must_reject(q.verify)
    q=ledger();q.nodes['0']['cell']={'x':[0,'1/10']};checks['root_node_mutation_rejected']=must_reject(q.verify)
    q=ledger();q.root={'x':[0,'1/10']};checks['live_root_declaration_mutation_rejected']=must_reject(q.verify)
    q=ledger();q.domain['x']['upper']='1';checks['live_ambient_declaration_mutation_rejected']=must_reject(q.verify)
    q=ledger();q.uncovered_physics=[];checks['live_uncovered_physics_mutation_rejected']=must_reject(q.verify)
    unbounded=engine.CoverLedger(box,{'x':{'lower':0,'upper':None}},uncovered_physics=[]) if False else engine.CoverLedger(box,{'x':{'lower':None,'upper':None}},uncovered_physics=['zero faces, tails'])
    checks['unbounded_tails_not_erased']=unbounded.verify()['ambient_remainder_status']=='UNCOVERED'
    unbounded.remainder['status']='EMPTY';checks['forged_EMPTY_tail_rejected']=must_reject(unbounded.verify)
    for name,mut in [('root',lambda s:s['nodes']['0'].update(cell={'x':[0,'1/10']})),('proofless',lambda s:s['nodes']['0'].update(classification={'status':'OUT'})),('remainder',lambda s:s['outside_remainder'].update(status='UNCOVERED'))]:
        s=copy.deepcopy(export);mut(s);checks['standalone_'+name+'_mutation_rejected']=must_reject(lambda:engine.verify_cover_snapshot(s))
    q=engine.CoverLedger(box,{'x':{'lower':None,'upper':None}},uncovered_physics=[]);a,b=q.subdivide('0','x')
    q.nodes[a]['cell']={'x':['-1/10','-1/20']};checks['child_gap_rejected']=must_reject(q.verify)
    q=engine.CoverLedger(box,{'x':{'lower':None,'upper':None}},uncovered_physics=[]);a,b=q.subdivide('0','x');q.nodes[a]['parent']='wrong';checks['wrong_child_parent_rejected']=must_reject(q.verify)
    q=ledger();q.nodes['orphan']=copy.deepcopy(q.nodes['0']);checks['orphan_node_rejected']=must_reject(q.verify)
    q=ledger();changed=call(dd={k:2 for k in K});checks['mixed_discrepancy_context_rejected']=must_reject(lambda:q.attach('0',changed))
    q=ledger();q.nodes['0']['classification']=copy.deepcopy(full);q.nodes['0']['classification']['status']='OUT';checks['forged_leaf_status_rejected_by_proof_replay']=must_reject(q.verify)
    # API semantic premise: contract hash is opaque, not a parsed operator schema.
    limitation={'operator_contract_semantics':'Caller must link static operator enum to reviewed model/observation contract. Generic engine does not parse an opaque contract hash; mismatched semantic premises are not validated internally.',
        'audit_trust_boundary':'Independent audit files are trusted mathematical attestations bound by hashes; arbitrary self-authored PASS JSON is not proof or an authenticated reviewer identity.',
        'snapshot_scope':'Standalone verification proves coverage of the declared snapshot root/context. Binding that declaration to an earlier request requires an external recorded hash.'}
    write(OUT/'analytic_results.json',rows)
    source=ROOT/'scripts/validated_identified_set.py';deps=ROOT/'scripts/joint_calibration_contract.py'
    result={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'engine_source_sha256':sha(source),'dependency_source_sha256':sha(deps),'review_source_sha256':sha(__file__),
        'own_fixture_directory':OUT.relative_to(ROOT).as_posix(),'own_fixture_hashes':{p.relative_to(ROOT).as_posix():sha(p) for p in OUT.rglob('*.json')},
        'semantic_preconditions_and_trust_boundary':limitation,'actual_cover_review':'PENDING','no_physical_or_producer_fixture_changes':True,'elapsed_seconds':time.time()-started}
    write(RECEIPT,result);print(json.dumps({'status':result['status'],'checks':len(checks),'failed':[k for k,v in checks.items() if not v],'elapsed':result['elapsed_seconds']},indent=2))
    if not all(checks.values()):raise RuntimeError('Independent inference-layer finding retained')


def actual():
    record=read(RECEIPT)
    if record['status']!='PASS' or record['actual_cover_review']!='PENDING':raise RuntimeError('fresh independent engine review required')
    started=time.time();checks={};base=ROOT/'runs/validated_identified_set_v1';single=read(ROOT/'evidence/validated_identified_set_actual_M.json');broad=read(ROOT/'evidence/validated_identified_set_actual_M_cover.json')
    protocol=read(base/'actual_M/protocol.json');bprotocol=read(base/'actual_M_cover_v2/protocol.json');package=read(base/'actual_M/forward_package.json')
    combined=read(ROOT/'evidence/parametric_actual_M_independent_review.json');obsfile=ROOT/'runs/parametric_root_cells_v1/actual_M/observer_certificate.json';obsproof=read(obsfile);contractfile=obsfile.with_name('model_observation_contract.json');contract=read(contractfile)
    observations=engine.load_observations()
    with (ROOT/'data/derived/calibration_small_18um.csv').open(encoding='utf-8-sig',newline='') as f:doses={r['record_id']:r['activation_magnitude'] for r in csv.DictReader(f)}
    checks['actual_source_protocol_receipt_hash_closure']=single['status']==broad['status']=='PASS' and all(v is True for v in single['checks'].values()) and all(v is True for v in broad['checks'].values()) and single['protocol_sha256']==sha(base/'actual_M/protocol.json') and broad['protocol_sha256']==sha(base/'actual_M_cover_v2/protocol.json') and sha(base/'apply_actual_M.py')==protocol['generator_sha256'] and sha(base/'expand_actual_M_cover.py')==bprotocol['generator_sha256'] and sha(ROOT/'scripts/validated_identified_set.py')==protocol['engine_source_sha256']==bprotocol['engine_source_sha256']
    checks['combined_actual_M_independent_review_bound']=combined['status']=='PASS' and sha(ROOT/'evidence/parametric_actual_M_independent_review.json')==protocol['combined_independent_review_sha256'] and combined['root_audit']['status']==combined['observer_audit']['status']=='PASS'
    for section,name in [('root_audit','root_review_adapter.json'),('observer_audit','observer_review_adapter.json')]:
        adapter=read(base/'actual_M'/name);der=adapter.pop('transparent_derivation')
        checks[section+'_adapter_exact_frozen_review_section']=adapter==combined[section] and der['section']==section and der['new_mathematical_claims'] is False and der['combined_review']['sha256']==sha(ROOT/'evidence/parametric_actual_M_independent_review.json')
    checks['actual_static_spatialmax_contract_and_no_mean_identity_shortcut']=contract['primary_variant']=='current_disk__spatial_max' and contract['component_aggregation'].startswith('spatial_max') and contract['time_observer'].startswith('single static') and single['primary_decision']['inference_context']['operator'] is None and single['primary_decision']['conditional_static_operator_necessary_test'] is None and sha(contractfile)==protocol['contract_sha256']==single['primary_decision']['model_observation_contract_sha256']
    exact_package=sha(base/'actual_M/forward_package.json')==single['forward_package_sha256']
    bounds,prov,issues=engine.uniform_bounds(package,single['primary_decision']['parameter_cell'],doses,sha(contractfile))
    checks['actual_provider_accepts_exactly_reviewed_one_channel']=exact_package and len(bounds)==1 and set(bounds['optoGEF2_small_col_M'])=={'apical_max_shortening_um'} and not issues and len(prov)==1
    lo,hi=map(Q,obsproof['endpoint_intervals']['apical_max_shortening_um']);y=observations['optoGEF2_small_col_M']['apical_max_shortening_um'];cutlo=max(lo-y,y-hi,Q(0));cuthi=max(abs(lo-y),abs(hi-y))
    cold=True;exact=True
    for profile in single['conditional_delta_family']:
        d=Q(profile['common_delta_exact']);saved=profile['decision'];fresh=engine.classify_cell(observations,{k:d for k in engine.OBSERVABLES},doses,saved['parameter_cell'],contract_hash=sha(contractfile),package=package)
        cold &= engine.digest(fresh)==engine.digest(saved)
        row=next(r for r in saved['constraints'] if r['record_id']=='optoGEF2_small_col_M' and r['observable']=='apical_max_shortening_um')
        exact &= saved['status']==('OUT' if d<cutlo else 'UNRESOLVED') and row['status']==('OUT' if d<cutlo else 'IN' if d>=cuthi else 'UNRESOLVED') and sum(r['status']=='UNRESOLVED' for r in saved['constraints'] if r is not row)==67 and saved['all_branches_existential_parameter_status']=='UNRESOLVED' and saved['constraint_count']==68 and saved['conditional_static_operator_necessary_test'] is None
    checks['all_actual_tolerance_decisions_cold_replayed']=cold
    checks['exact_delta_cutoffs_touching_and_67missing_all_profiles']=exact
    checks['declared_5p1_tolerance_is_strict_selected_family_OUT']=cutlo>Q('5.1') and single['primary_decision']['status']=='OUT' and single['primary_decision']['common_path_or_energetic_stability_certified'] is False
    snap=broad['cover_ledger'];nodes=snap['nodes'];visited=set();leaves=[];valid=True
    def walk(name):
        nonlocal valid
        if name in visited:valid=False;return
        visited.add(name);node=nodes[name];C={k:tuple(map(Q,v)) for k,v in node['cell'].items()}
        if not node['children']:leaves.append(name);return
        if len(node['children'])!=2:valid=False;return
        left,right=(nodes[k] for k in node['children']);coord=node['split_coordinate'];CL={k:tuple(map(Q,v)) for k,v in left['cell'].items()};CR={k:tuple(map(Q,v)) for k,v in right['cell'].items()}
        valid &= C[coord][0]==CL[coord][0]<CL[coord][1]==CR[coord][0]<CR[coord][1]==C[coord][1] and all(CL[k]==CR[k]==v for k,v in C.items() if k!=coord)
        valid &= all(nodes[ch]['parent']==name for ch in node['children'])
        for ch in node['children']:walk(ch)
    walk('0')
    def volume(C):
        v=Q(1)
        for a,b in C.values():v*=Q(b)-Q(a)
        return v
    checks['independent_closed_union_tree_induction_nine_cuts']=valid and visited==set(nodes) and len(nodes)==19 and len(leaves)==10 and len(broad['exact_splits'])==9
    checks['exact_5D_leaf_volumes_sum_to_seed']=sum((volume(nodes[n]['cell']) for n in leaves),Q(0))==volume(snap['root_cell'])
    target=broad['target_leaf'];expectedhash=combined['root_audit']['parameter_cell_sha256']
    checks['only_exact_target_leaf_covered_and_selected_OUT']=engine.cell_hash(nodes[target]['cell'])==expectedhash and nodes[target]['classification']['status']=='OUT' and sum(nodes[n]['classification']['status']=='OUT' for n in leaves)==1
    checks['other_nine_cells_not_extrapolated_from_tiny_certificate']=all(nodes[n]['classification']['status']=='UNRESOLVED' and not nodes[n]['classification']['uniform_constraints_covered'] and all(r['status']=='UNRESOLVED' for r in nodes[n]['classification']['constraints']) for n in leaves if n!=target)
    checks['all_ten_leaves_retain68_and_otherbranches_unresolved']=all(nodes[n]['classification']['constraint_count']==68 and nodes[n]['classification']['record_count']==17 and nodes[n]['classification']['all_branches_existential_parameter_status']=='UNRESOLVED' and nodes[n]['classification']['inference_context']['operator'] is None for n in leaves)
    checks['target_keeps_other67missing']=sum(r['status']=='UNRESOLVED' for r in nodes[target]['classification']['constraints'])==67
    replay=engine.verify_cover_snapshot(snap)
    checks['actual_standalone_export_replays_exactly']=replay==snap['verification'] and replay['exact_closed_leaf_union_equals_root_box'] and len(replay['unresolved_leaves'])==9
    checks['nonnegative_zero_faces_and_unbounded_tails_explicit']=all(Q(v['lower'])==0 and v['lower_strict'] is False and v['upper'] is None for v in snap['ambient_domain'].values()) and snap['outside_remainder']['status']=='UNCOVERED' and snap['outside_remainder']['unbounded_tail_proof'] is None
    checks['finite_and_global_completion_not_overstated']=replay['finite_selected_branch_cover_complete'] is False and replay['global_branch_path_certified'] is False and replay['full_physical_identified_set_complete'] is False
    for name,mut in [('proofless_target_OUT',lambda s:s['nodes'][target].update(classification={'status':'OUT'})),('root_node',lambda s:s['nodes']['0'].update(cell=nodes[target]['cell'])),('fakeEMPTYtail',lambda s:s['outside_remainder'].update(status='EMPTY'))]:
        changed=copy.deepcopy(snap);mut(changed);write(OUT/(name+'.json'),changed);checks['actual_export_'+name+'_attack_rejected']=must_reject(lambda:engine.verify_cover_snapshot(changed))
    checks['source_dependency_matches_independent_analytic_review']=sha(ROOT/'scripts/validated_identified_set.py')==record['engine_source_sha256'] and sha(ROOT/'scripts/joint_calibration_contract.py')==record['dependency_source_sha256']
    files=[ROOT/'evidence/validated_identified_set_actual_M.json',ROOT/'evidence/validated_identified_set_actual_M_cover.json',base/'actual_M/protocol.json',base/'actual_M_cover_v2/protocol.json',base/'actual_M/forward_package.json',base/'actual_M/root_review_adapter.json',base/'actual_M/observer_review_adapter.json',base/'apply_actual_M.py',base/'expand_actual_M_cover.py',ROOT/'evidence/parametric_actual_M_independent_review.json',obsfile,contractfile,ROOT/'data/derived/calibration_small_18um.csv']
    record['actual_cover_review']='PASS' if all(checks.values()) else 'FAIL';record['actual_checks']=checks;record['actual_input_hashes']={p.relative_to(ROOT).as_posix():sha(p) for p in files}
    record['actual_cover_summary']={'cuts':9,'finite_leaves':len(leaves),'selected_OUT':1,'unresolved_leaves':9,'target_leaf':target,'target_missing_constraints':67,'cutoffs_exact':[str(cutlo),str(cuthi)],'operator':None,'contract_primary_variant':contract['primary_variant'],'otherbranches':'UNRESOLVED','ambient_remainder':'UNCOVERED','full_physical_identified_set_complete':False}
    record['current_review_source_sha256']=sha(__file__);record['status']='PASS' if all(record['checks'].values()) and all(checks.values()) else 'FAIL';record['actual_review_elapsed_seconds']=time.time()-started
    write(RECEIPT,record);print(json.dumps({'status':record['status'],'actual_checks':len(checks),'failed':[k for k,v in checks.items() if not v],'elapsed':record['actual_review_elapsed_seconds'],'summary':record['actual_cover_summary']},indent=2))
    if record['status']!='PASS':raise RuntimeError('Independent actual inference-layer finding retained')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--actual',action='store_true');args=p.parse_args()
    actual() if args.actual else main()
