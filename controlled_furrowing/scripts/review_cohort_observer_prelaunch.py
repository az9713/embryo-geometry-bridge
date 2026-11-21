"""Review static-cohort composition before any root outputs are consumed."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,importlib,itertools,json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs/cohort_static_observers_v1'
read=lambda p:json.loads(Path(p).read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    out=RUN/'prelaunch_independent_review.json'
    if out.exists():raise RuntimeError('preserve concrete launch review')
    m=importlib.import_module('cohort_static_observers')
    p=read(RUN/'protocol.json');pre=read(RUN/'preflight.json');tests=read(RUN/'synthetic_checks.json')
    contract=read(RUN/'model_observation_contract.json')
    checks={
        'observer_source_bound':p['source_hashes']['scripts/cohort_static_observers.py']==sha(ROOT/'scripts/cohort_static_observers.py')==pre['source_sha256']==tests['source_sha256'],
        'observer_protocol_preflight_bound':pre['protocol_sha256']==sha(RUN/'protocol.json'),
        'preflight_all_boolean_true':pre['status']=='PASS' and bool(pre['checks']) and all(v is True for v in pre['checks'].values()),
        'synthetics_all_boolean_true':tests['status']=='PASS' and bool(tests['checks']) and all(v is True for v in tests['checks'].values()),
        'all_recursive_sources_and_inputs_current':all(sha(ROOT/k)==v for k,v in {**p['source_hashes'],**p['input_hashes']}.items()),
        'root_prelaunch_review_required_and_bound':'runs/cohort_uniform_roots_v1/prelaunch_independent_review.json' in p['input_hashes'],
        'same_fixed_common_contract':sha(RUN/'model_observation_contract.json')==p['common_contract_sha256']=='2498493bf50049295eb2ad0d2360b1defaf387ee8db8828ceb576f3d7e942f82',
        'same_shared_parameter_cell':contract['parameter_cell_sha256']==p['parameter_cell_sha256']=='64f93fefa034f4c004b93451d98b2fd04e0072855be8c5e3bba9babb8f394dbf',
        'full_corrector_not_posterior_tightening':p['full_corrector_cube_required'] is True,
        'no_solve_or_noise_assumption':p['no_new_physical_solve'] is True and p['no_biological_tolerance_assumed'] is True,
        'independent_review_still_required':p['producer_results_require_independent_review'] is True,
        'stationary_closure_only':p['branch_law_semantics']=='independent_stationary_root_per_record' and contract['common_path_or_energetic_stability_certified'] is False and contract['source_temporal_semantics_validated'] is False,
        'no_composition_already_written':not (RUN/'results.json').exists() and not (RUN/'composition_final/results.json').exists(),
        'generic_body_and_transport_have_prior_independent_review':all(read(ROOT/'evidence'/name)['status'].startswith('PASS') for name in ('root_body_geometry_independent_review.json','nodal_observer_error_independent_review.json'))}
    # Independent exact affine-box arithmetic, including nonzero bases.
    names=m.ORDER;C={k:[str(-i-2),str(i+3)] for i,k in enumerate(names)};base={k:str(Q(i,3)) for i,k in enumerate(names)}
    S=np.array([[1.,-.25,2.,-.5,.125],[-2.,.5,0.,1.,-.25]])
    got=m.predictor_hulls(S,C,base);want=[]
    for row in S:
        values=[]
        for vertex in itertools.product(*[tuple(map(Q,C[k])) for k in names]):
            values.append(sum((Q.from_float(float(s))*(v-Q(base[k])) for s,k,v in zip(row,names,vertex)),Q(0)))
        want.append((min(values),max(values)))
    checks['independent_two_rows_all32_vertices_match']=got==want
    expected={'dose_exact':'78208913/100000000','base_strength_exact':'1/2','quantifier':'q','branch_scope':'b'}
    m.validate_root_context({'dose_exact':'0.78208913','base_strength_exact':'0.5','quantifier':'q','branch_scope':'b','parameter_cell':C},expected,C)
    checks['rational_numeric_aliases_validated_without_id_normalization']=True
    observations={'a':{k:'1' for k in m.OBSERVABLES},'b':{k:'3' for k in m.OBSERVABLES}}
    predicted={'a':{k:['-2','0'] for k in m.OBSERVABLES},'b':{k:['4','5'] for k in m.OBSERVABLES}}
    a=m.aggregate(observations,predicted)
    checks['two_record_exact_residual_range']=a['maximum_residual_lower_exact']==1 and a['maximum_residual_upper_exact']==3 and a['covered']==8
    del predicted['b'][m.OBSERVABLES[-1]];a=m.aggregate(observations,predicted)
    checks['any_missing_row_forbids_sufficient_upper']=a['maximum_residual_upper_exact'] is None and a['common_delta_family']['IN_static_stationary_root_tuple_if'] is None and a['covered']==7
    checks['wider_valid_interval_policy_explicit']=bool(p.get('maximum_reporting_policy')) and tests['checks'].get('valid_wider_maximum_enclosure_accepted') is True
    source=(ROOT/'scripts/cohort_static_observers.py').read_text(encoding='utf-8')
    checks['per_surface_radius_and_height_sum_present']="epsilon=L*max(radii[int(n)][j] for edge in ids[name] for n in edge for j in range(2))" in source and 'errors=dict(zip(OBSERVABLES,(ea+eb,ea+eb,ea,eb)))' in source
    checks['exact_serialization_gap_in_full_nodal_box']='abs(exact[i][j]-centers[i][j])' in source and 'max(map(abs,hulls[i]))+radius' in source
    checks['each_body_and_both_fixed_ROI_endpoints_checked']='body=certify_box(' in source and "radii[node][0]!=0 or gaps[node][0]!=0" in source and 'fixed_radial_endpoints=True,roi_um=roi' in source
    checks['provider_manifest_and_protocol_hashes_locked_before_use']='locks[selection_file]=sha(path(selection_file))' in source and "providers[provider['record_id']]=provider" in source
    checks['actual_base_width_and_reporting_flag_retained']="'base_maximum_reporting':maximum_report" in source and "observed.get('status')!='CERTIFIED_STATIC_INTERVALS'" in source
    result={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,
        'observer_source_sha256':sha(ROOT/'scripts/cohort_static_observers.py'),
        'observer_protocol_sha256':sha(RUN/'protocol.json'),'observer_preflight_sha256':sha(RUN/'preflight.json'),
        'review_source_sha256':sha(__file__),
        'mathematical_review':'Exact affine predictor hull plus full corrector and serialization gap gives a nodal box. Independently reviewed body/curve lemmas apply separately to each boundary. With common fixed current ROI, pointwise Eulerian errors ea,eb imply height max and positive normalized disk mean errors at most ea+eb; component maxima have errors ea and eb. Exact interval residual extrema then bound max of all68 residuals. All68 are required for a sufficient common upper. This is conditional on each future root/provider packet and still requires completed independent reconstruction.',
        'scope':'Prelaunch source/contract/math-composition review; no completed cohort output, all17 root proof, source temporal validation, biological fit, continuum or novelty is certified.'}
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'checks':len(checks),'failed':[k for k,v in checks.items() if not v],'receipt_sha256':sha(out)}))
    if result['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
