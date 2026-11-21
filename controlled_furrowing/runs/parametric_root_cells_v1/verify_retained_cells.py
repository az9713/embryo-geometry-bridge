"""Replay retained exact inequalities, cover accounting and artifact hashes."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,re

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
checks={};counts={'cases':0,'certified':0,'unresolved':0,'radius_tests':0};inputs={}
for folder in(OUT,OUT/'actual_M'):
    label='pilot'if folder==OUT else'actual_M';p=read(folder/'protocol.json');r=read(folder/'results.json');tests=read(folder/'tests.json')
    checks[label+'_protocol_binding']=r['protocol_sha256']==sha(folder/'protocol.json')==tests['protocol_sha256']
    checks[label+'_tests_pass']=tests['status']=='PASS'and all(tests['checks'].values())and all(tests['kernel_checks']['checks'].values())
    source=OUT/'pilot_source_snapshot/parametric_root_cells.py'if folder==OUT else ROOT/'scripts/parametric_root_cells.py'
    checks[label+'_frozen_source_binding']=sha(source)==p['source_hashes']['parametric_root_cells.py']
    checks[label+'_kernel_binding']=sha(ROOT/'scripts/algebraic_root_feasibility.py')==p['source_hashes']['algebraic_root_feasibility.py']
    checks[label+'_predictor_binding']=sha(folder/p['predictor_file'])==p['predictor_sha256']
    checks[label+'_no_false_global_completion']=r['full_identified_set_computed']is False and r['no_OUT_classification']is True and r['cover_ledger']['full_physical_domain_covered']is False
    for row in r['cell_attempts']:
        counts['cases']+=1;counts['certified'if row['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT'else'unresolved']+=1
        assert row['status']in('CERTIFIED_UNIFORM_LOCAL_ROOT','UNRESOLVED')
        assert all(F(pair[0])<F(pair[1])for pair in row['cell_bounds_exact'].values())
        passed=[]
        for a in row['radius_attempts']:
            counts['radius_tests']+=1
            if'Z_uniform'not in a:assert'failure'in a;continue
            z=F(a['Z_uniform']['upper_exact']);eta=F(row['eta_uniform']['upper_exact']);radius=F(a['radius_exact']);lhs=F(a['eta_plus_Zr']['upper_exact'])
            assert lhs>=eta+z*radius
            valid=z<1 and lhs<radius
            assert valid==a['strict_selfmap_contraction']
            if valid:passed.append(a)
        assert bool(passed)==(row['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT')
        if passed:
            a=passed[0];z=F(a['Z_uniform']['upper_exact']);eta=F(row['eta_uniform']['upper_exact'])
            assert F(row['uniform_error_from_predictor']['upper_exact'])>=eta/(1-z)
            assert F(row['uniform_outer_nodal_radius_about_base']['upper_exact'])>=max(map(F,row['full_nodal_predictor_bounds_exact']))+F(row['uniform_error_from_predictor']['upper_exact'])
    ledger=r['cover_ledger'];assert len(ledger['split_records'])==1 and len(ledger['leaves'])==2
    split=ledger['split_records'][0];name=split['split_coordinate'];parent=ledger['root_cell'];children=split['children']
    assert F(children[0][name][0])==F(parent[name][0])and F(children[1][name][1])==F(parent[name][1])and F(children[0][name][1])==F(children[1][name][0])
    assert all(child[k]==parent[k]for child in children for k in parent if k!=name)
    checks[label+'_closed_child_union_and_all_leaves_retained']=all(l['status']=='UNRESOLVED'for l in ledger['leaves'])
    for f in(folder/'protocol.json',folder/'results.json',folder/'tests.json'):inputs[f.relative_to(ROOT).as_posix()]=sha(f)
checks['all28_cases_retained']=counts['cases']==28 and counts['certified']==10 and counts['unresolved']==18
p=OUT/'actual_M';root=read(p/'root_certificate.json');obs=read(p/'observer_certificate.json');diagnostic=read(p/'contradiction_diagnostic.json');body=read(p/'body_geometry.json');contract=read(p/'model_observation_contract.json')
checks['root_observer_model_cell_bindings']=obs['root_certificate_sha256']==sha(p/'root_certificate.json')and root['model_observation_contract_sha256']==obs['model_observation_contract_sha256']==sha(p/'model_observation_contract.json')and root['parameter_cell_sha256']==obs['parameter_cell_sha256']
checks['stationary_branch_not_common_path']=root['branch_law_semantics']==obs['branch_law_semantics']=='independent_stationary_root_per_record'and not root['common_path_or_energetic_stability_certified']and not obs['common_path_or_energetic_stability_certified']
checks['correct_actual_record_dose_strength_law']=root['record_id']=='optoGEF2_small_col_M'and root['dose_exact']=='0.78208913'and root['strength_relation']=='common_gain_times_source_dose'
checks['conservative_full_corrector_composition']=obs['corrector_radius_exact']==root['uniform_cell_result']['selected_radius_exact']=='1/100000000000'
checks['every_body_field_geometrically_admissible']=body['status']=='CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING'and body['global_injectivity_for_entire_box']and body['all_element_J_positive']
lo,hi=map(F,obs['endpoint_intervals']['apical_max_shortening_um']);y=F(12);res=(y-hi,y-lo)
checks['exact_entire_delta_family_cuts']=list(map(F,diagnostic['residual_interval_um_exact']))==list(res)and F(diagnostic['excluded_for_delta_apical_below_exact'])==res[0]and F(diagnostic['one_constraint_uniformly_satisfied_for_delta_at_least_exact'])==res[1]
checks['declared_5point1_example_not_noise']=F(diagnostic['declared_delta_um_exact'])==F('5.1')<res[0]and diagnostic['delta_is_measured_noise_estimate']is False
checks['only_one_of68_observations_supplied']=set(obs['endpoint_intervals'])=={'apical_max_shortening_um'}and obs['other67constraints_not_supplied']is True
checks['operator_is_primary_static_disk_spatialmax']=contract['primary_variant']=='current_disk__spatial_max'and not contract['source_temporal_semantics_validated']
doc=ROOT/'docs/PARAMETRIC-ROOT-CELLS.md';text=doc.read_text();checks['documentation_math_balanced']=text.count('\\[')==text.count('\\]')and text.count('\\(')==text.count('\\)')
for target in re.findall(r'\]\(([^)]+)\)',text):
    if not target.startswith(('http:','https:','#')):assert(doc.parent/target).resolve().exists(),target
checks['documentation_local_links_exist']=True
for f in(p/'root_certificate.json',p/'observer_certificate.json',p/'body_geometry.json',p/'model_observation_contract.json',p/'contradiction_diagnostic.json',doc):inputs[f.relative_to(ROOT).as_posix()]=sha(f)
assert all(checks.values()),checks
receipt={'status':'AUTHOR_EXACT_ARITHMETIC_AND_PROVENANCE_PASS','checks':checks,'counts':counts,'input_hashes':inputs,'source_sha256':sha(__file__),
 'scope':'All retained inequality/ledger/hash checks, not an independent rerun of all28 interval assemblies. SelectedactualMcell has its separate independentroot/observer review.'}
(OUT/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'checks':len(checks),'counts':counts}))
