"""Exact scalar composition of reviewed B integral-root transfer premises.

No new equilibrium solve, integration, material fit, PDE or biological claim.
All predicates are exact rational sufficient conditions; failed tests remain
UNRESOLVED. Preparation does not assert that a root has been certified.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,datetime,base64,io
import numpy as np
ROOT=Path(__file__).resolve().parents[1];SOURCE=Path(__file__)
BASE=ROOT/'runs/exact_integration_bridge_v1'
OUT=BASE/'integral_root_B_safe_v2'
SOURCES={
 'point_review':ROOT/'evidence/exact_complete_point_B_independent_review_safe_v2.json',
 'point_exact':BASE/'bulk_independent_replay_safe_v2/exact_168_reduction.json',
 'retube':BASE/'retube_B_preparation/protocol.json',
 'vectors':BASE/'retube_B_preparation/candidate_vectors.npz',
 'geometry':BASE/'retube_B_geometry/protocol.json',
 'geometry_review':ROOT/'evidence/exact_retube_B_geometry_independent_review.json',
 'cached':BASE/'cached_Q_common_B_v1/results.json',
 'cached_review':ROOT/'evidence/exact_cached_Q_common_B_reduction_review.json',
 'first_bound':BASE/'uniform_full_B_independent_reduction_v1/results.json',
 'physical_bound':BASE/'uniform_full_B_physical_axis_independent_reduction_v2/results.json',
 'parameter':ROOT/'evidence/exact_parameter_defect_B_independent_reduction.json',
 'parameter_protocol':BASE/'parameter_defect_B_safe_v2/protocol.json',
 'parameter_producer':BASE/'parameter_defect_B_safe_v2/results.json',
 'cached_protocol':BASE/'cached_Q_common_B_v1/protocol.json',
 'energy_review':ROOT/'evidence/cohort_energy_stability_independent_review.json',
 'energy_B':ROOT/'runs/cohort_energy_stability_review_v1/records/optoGEF2_small_col_B.json',
 'lemma':BASE/'retube_B_safe_relink/ROOT-TRANSFER-LEMMA.md',
}
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 if p.exists():raise ValueError('preserve existing output '+str(p))
 p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n',encoding='utf8')
def passed(d):return isinstance(d,dict)and bool(d)and all(v is True for v in d.values())
def gate(r,s,p,f):
 if not isinstance(r,dict)or r.get('status')!='PASS'or not passed(r.get('checks'))or r.get('source_sha256')!=s or r.get('protocol_sha256')!=p or r.get('preflight_sha256')!=f:raise ValueError('bound nonempty literal PASS review required')
def sufficient(etaQ,m,zQ,b,dnorm,s,parameter,r):
 vals=list(map(Q,(etaQ,m,zQ,b,dnorm,s,parameter,r)))
 if min(vals)<0 or vals[-1]==0:raise ValueError('nonnegative bounds and positive radius required')
 etaQ,m,zQ,b,dnorm,s,parameter,r=vals;z=zQ+b;V=b*s+parameter;eta=etaQ+m+z*dnorm+V
 return {'z':z,'V':V,'eta':eta,'selfmap_lhs':eta+z*r,'strict_sufficient_conditions':z<1 and eta+z*r<r}
def tests():
 checks={}
 checks['exact_sufficient_pass']=sufficient(Q(1,100),0,Q(1,10),Q(1,100),Q(1,100),Q(1,10),0,Q(1,10))['strict_sufficient_conditions']
 checks['large_predictor_variation_not_omitted']=not sufficient(Q(1,100),0,Q(1,10),Q(1,100),Q(1,100),10,0,Q(1,10))['strict_sufficient_conditions']
 checks['contraction_only_not_selfmap']=not sufficient(2,0,Q(1,10),0,0,0,0,1)['strict_sufficient_conditions']
 checks['z_equal_one_rejected']=not sufficient(0,0,Q(1,2),Q(1,2),0,0,0,1)['strict_sufficient_conditions']
 checks['strict_boundary_not_pass']=not sufficient(Q(1,2),0,Q(1,2),0,0,0,0,1)['strict_sufficient_conditions']
 checks['map_center_escape_does_not_exclude_root']=Q(6,5)-Q(1,2)*Q(4,5)==Q(4,5)and Q(4,5)<1<Q(6,5)
 # Exact nonsymmetric R / symmetric SPD target: the homotopy argument does
 # not require R to be symmetric.
 R=[[Q(1),Q(1,10)],[Q(0),Q(1)]];H=[[Q(1,2),Q(1,20)],[Q(1,20),Q(1)]]
 z=max(sum(abs(Q(int(i==j))-sum(R[i][k]*H[k][j]for k in range(2)))for j in range(2))for i in range(2))
 checks['nonsymmetric_R_SPD_fixture_contractive']=z<1 and H[0][0]>0 and H[0][0]*H[1][1]-H[0][1]**2>0
 checks['root_not_asserted_by_preparation']=True
 if not passed(checks):raise AssertionError(checks)
 return checks
def prepare():
 if (OUT/'protocol.json').exists():raise ValueError('preserve preparation')
 checks=tests();inputs={p.relative_to(ROOT).as_posix():sha(p)for p in SOURCES.values()}
 inputs[SOURCE.relative_to(ROOT).as_posix()]=sha(SOURCE)
 # Full current reviewer dependency bindings are retained, rather than using
 # an isolated PASS label as proof of a different file snapshot.
 for name in ['point_review','physical_bound','first_bound','parameter','parameter_protocol','cached_review','cached_protocol']:
  doc=read(SOURCES[name]);binds=doc.get('input_hashes',doc.get('source_and_input_hashes',{}))
  for path,h in binds.items():
   if sha(ROOT/path)!=h:raise ValueError('changed premise dependency '+path)
   if path in inputs and inputs[path]!=h:raise ValueError('conflicting premise binding '+path)
   inputs[path]=h
 record=read(SOURCES['retube'])['selected_record']
 for name in ['point_cache_file','predictor_file','root_certificate_file','independent_root_capsule_file']:
  path=record[name];expected=record[name.replace('_file','_sha256')]
  if sha(ROOT/path)!=expected:raise ValueError('changed authoritative selected-record input '+path)
  if path in inputs and inputs[path]!=expected:raise ValueError('conflicting authoritative selected-record binding '+path)
  inputs[path]=expected
 capsule=read(ROOT/record['independent_root_capsule_file'])
 entry=next(v for v in capsule['record_audits']if v['record_id']==record['record_id'])
 path=entry['review_file'];expected=entry['review_sha256']
 if sha(ROOT/path)!=expected or(path in inputs and inputs[path]!=expected):raise ValueError('changed selected independent root capsule')
 inputs[path]=expected
 p={'source_and_input_hashes':inputs,'source_files':{k:v.relative_to(ROOT).as_posix()for k,v in SOURCES.items()},
    'record_id':record['record_id'],'parameter_cell':record['parameter_cell'],'parameter_cell_sha256':record['parameter_cell_sha256'],
    'dose_exact':record['dose_exact'],'root_cube_radius_exact':record['selected_corrector_radius_exact'],
    'candidate_rule':'p1(theta)=u0+S(theta-theta0)-d, exact d=actual binary u0 minus u1.',
    'equations':'z=zQ+b; V=b*s+parameter; etaI=etaQ+m+z*||d||+V; require z<1 and etaI+z*r<r.',
    'local_stability_rule':'Same cached base H0 has independently proved SPD. Symmetric Ht=(1-t)H0+t*HI stays nonsingular under the convex ||I-RHt||<1 bound; inertia is constant. lambda_min(HI)>=(1-z)/||R||infinity.',
    'scope':'ONE B source dose, all theta in SAME small C, fixed 611-coordinate axisymmetric EXACT-SPATIAL-INTEGRAL finite-element target. Not all17 integral roots, PDE, global minimum/path, physical identified set, observational fit or confirmed novelty.',
    'root_asserted_at_preparation':False}
 write(OUT/'protocol.json',p);write(OUT/'preflight.json',{'status':'PASS_PREPARATION_ONLY','checks':checks,'source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json'),'root_certified':False})
 print(json.dumps({'phase':'prepared','source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json')}),flush=True)
def run():
 p=read(OUT/'protocol.json');gate(read(OUT/'prelaunch_independent_review.json'),sha(SOURCE),sha(OUT/'protocol.json'),sha(OUT/'preflight.json'))
 for path,h in p['source_and_input_hashes'].items():
  if sha(ROOT/path)!=h:raise ValueError('changed frozen premise '+path)
 docs={k:read(ROOT/v)for k,v in p['source_files'].items()if k not in ['vectors','lemma']}
 checks={}
 def ck(k,v):
  checks[k]=bool(v)
  if not v:raise AssertionError(k)
 for name in ['point_review','cached','cached_review','physical_bound','first_bound','parameter','energy_review','energy_B','geometry_review']:
  ck('nonempty_pass_'+name,docs[name]['status'].startswith(('PASS','CERTIFIED_'))and passed(docs[name]['checks']))
 point=docs['point_exact'];audit=docs['point_review'];record=docs['retube']['selected_record'];geometry=docs['geometry'];cached=docs['cached'];physical=docs['physical_bound'];parameter=docs['parameter'];spd=docs['energy_B']
 ck('audited_independent_exact_point',sha(SOURCES['point_exact'])==audit['exact_reduction_sha256'])
 ck('same_original_closed_C',all(x['parameter_cell_sha256']==p['parameter_cell_sha256']for x in [geometry,cached,physical,spd]))
 ck('same_explicit_C_base_and_dose',all(x['parameter_cell']==p['parameter_cell']and x['base_parameters_exact']==record['base_parameters_exact']for x in [geometry,docs['cached_protocol'],docs['parameter_protocol']])and cached['dose_exact']==p['dose_exact']==record['dose_exact'])
 box_path=ROOT/docs['parameter_protocol']['common_geometry_boxes_file'];box_hash=sha(box_path);box=read(box_path)
 ck('same_canonical_box_and_domain',docs['parameter_protocol']['common_geometry_boxes_sha256']==box_hash==docs['cached_protocol']['geometry_boxes_sha256']==geometry['prepared_box_sha256']and physical['common_geometry_domain_id']==box['domain_id']==docs['cached_protocol']['geometry_domain_id'])
 oldroot=read(ROOT/record['root_certificate_file']);cachefile=ROOT/record['point_cache_file'];cache=read(cachefile)
 ck('authoritative_same_pointcache_and_predictor',all(oldroot[k]==record[k]for k in ['point_cache_file','predictor_file'])and all(oldroot['artifact_hashes'][record[k]]==record[k.replace('_file','_sha256')]==sha(ROOT/record[k])for k in ['point_cache_file','predictor_file']))
 ck('same_cached_protocol_target_and_predictor',docs['cached_protocol']['build_point_parameters']['point_certificate_file']==record['point_cache_file']and docs['cached_protocol']['build_point_parameters']['predictor_file']==record['predictor_file'])
 ck('physical_bound_same_R_and_S_provenance',all(physical['input_hashes'][record[k]]==record[k.replace('_file','_sha256')]for k in ['point_cache_file','predictor_file']))
 raw=base64.b64decode(cache['frozen_binary_cache']['data'],validate=True)
 ck('same_actual_binary_R_cache',hashlib.sha256(raw).hexdigest()==cache['frozen_binary_cache']['sha256']==record['frozen_binary_cache_sha256'])
 with np.load(io.BytesIO(raw),allow_pickle=False)as z:cacheu=z['u'].copy()
 ck('parameter_review_binds_producer_and_protocol',parameter['producer_results_sha256']==sha(SOURCES['parameter_producer'])and parameter['protocol_sha256']==sha(SOURCES['parameter_protocol']))
 ck('parameter_component_same_C_predictor_and_dose',docs['parameter_protocol']['parameter_cell']==p['parameter_cell']and docs['parameter_protocol']['dose_exact']==p['dose_exact']and docs['parameter_protocol']['predictor_file']==record['predictor_file'])
 ck('same_cached_base_H0_energy_proof',spd['root_certificate_sha256']==sha(ROOT/record['root_certificate_file'])and spd['record_id']==record['record_id']and Q(spd['uniform_Hessian_lower_bound_exact'])>0)
 ck('linked_base_SPD_record',next(v for v in docs['energy_review']['record_audits']if v['record_id']==record['record_id'])['review_sha256']==sha(SOURCES['energy_B']))
 with np.load(ROOT/p['source_files']['vectors'],allow_pickle=False)as z:
  ck('actual_base_vector_matches_cached_H0',np.array_equal(z['u0'],cacheu))
  delta=[Q.from_float(float(a))-Q.from_float(float(b))for a,b in zip(z['u0'],z['u1'])]
 ck('all611_candidate_coordinates',len(delta)==len(point['RF_endpoint_sums'])==611)
 dnorm=max(map(abs,delta));m=max(max(abs(Q(a)-v),abs(Q(b)-v))for(a,b),v in zip(point['RF_endpoint_sums'],delta))
 ck('same_common_radius_covers_sweep',Q(geometry['new_common_enclosing_radius_exact'])==Q(p['root_cube_radius_exact'])+dnorm==Q(cached['common_radius_exact'])==Q(physical['common_radius_exact']))
 etaQ=Q(point['eta_used_exact']);zQ=Q(cached['ZQ_upper_exact']);s=Q(parameter['predictor_displacement_norm_upper_exact']);parameter_bound=Q(parameter['parameter_only_variation_upper_exact']);r=Q(p['root_cube_radius_exact'])
 ck('same_cached_ZQ_review',cached['ZQ_upper_exact']==docs['cached_review']['ZQ_upper_exact'])
 ck('eta_uses_larger_author_and_independent_bound',etaQ==max(Q(point['eta_author_exact']),Q(point['eta_independent_exact'])))
 oldcapsule=read(ROOT/record['independent_root_capsule_file']);entry=next(v for v in oldcapsule['record_audits']if v['record_id']==record['record_id']);oldreview=read(ROOT/entry['review_file'])
 ck('same_author_eta_and_independent_capsule_eta',Q(oldroot['uniform_cell_result']['eta_uniform']['upper_exact'])==Q(point['eta_author_exact'])and Q(oldreview['bounds']['eta']['upper_exact'])==Q(point['eta_independent_exact']))
 ck('old_normal_radius_not_silently_changed',Q(point['old_radius_exact'])==r==Q(geometry['old_radius_exact'])==Q(spd['normal_root_radius_exact']))
 ck('old_root_normal_radius_and_C',Q(oldroot['uniform_cell_result']['selected_radius_exact'])==r and oldroot['uniform_cell_result']['cell_bounds_exact']==p['parameter_cell'])
 b=Q(physical['norm_bounds']['b_absR_H_upper_exact']);firstb=Q(docs['first_bound']['norm_bounds']['b_absR_H_upper_exact'])
 first=sufficient(etaQ,m,zQ,firstb,dnorm,s,parameter_bound,r);chosen=sufficient(etaQ,m,zQ,b,dnorm,s,parameter_bound,r)
 # R's exact norm is independently bound by the original point reduction and
 # old energy audit; the new cached family uses precisely this unchanged R.
 Rnorm=Q(point['Rnorm_exact']);ck('same_exact_Rnorm',Rnorm==Q(spd['exact_R_infinity_norm'])and Rnorm>0)
 ck('base_cached_Z0_less_than_one',Q(cached['Z0_upper_exact'])<1)
 ck('new_full_bound_not_greater_than_first',b<=firstb)
 accepted=chosen['strict_sufficient_conditions'];z=chosen['z'];error=chosen['eta']/(1-z)if accepted else None
 curvature=(1-z)/Rnorm if accepted else None
 location_lower=max(Q(0),abs(delta[588])-error)if accepted else None
 oldQ_error=etaQ/(1-zQ)if accepted else None
 branch_gap=max(Q(0),abs(delta[588])-error-oldQ_error)if accepted else None
 old_cube_excluded=bool(accepted and location_lower>r)
 status='CERTIFIED_UNIFORM_EXACT_FE_ROOT_AND_STRICT_LOCAL_MINIMUM'if accepted else'UNRESOLVED_SUFFICIENT_ROOT_TRANSFER_TEST'
 def serial(d):return {k:(str(v)if isinstance(v,Q)else v)for k,v in d.items()}
 out={'status':status,'checks':checks,'protocol_sha256':sha(OUT/'protocol.json'),'source_sha256':sha(SOURCE),
      'record_id':p['record_id'],'parameter_cell_sha256':p['parameter_cell_sha256'],'dose_exact':p['dose_exact'],
      'r_exact':str(r),'etaQ_exact':str(etaQ),'m_exact':str(m),'zQ_exact':str(zQ),'b_exact':str(b),'d_norm_exact':str(dnorm),
      'predictor_norm_exact':str(s),'parameter_only_exact':str(parameter_bound),'b_policy':'Fully exact-rational |R| times absolute Hessian-row majorant; tighter signed-Arb norm not required for this root proof.','first_representation_test':serial(first),
      'physical_axis_test':serial(chosen),'root_error_from_translated_predictor_upper_exact':str(error)if error is not None else None,
      'uniform_integral_energy_curvature_lower_exact':str(curvature)if curvature is not None else None,
      'location_witness_free_index':588,'integral_root_distance_from_old_predictor_lower_exact':str(location_lower)if location_lower is not None else None,
      'old_Q_root_error_upper_exact':str(oldQ_error)if oldQ_error is not None else None,
      'I_Q_branch_component_separation_lower_exact':str(branch_gap)if branch_gap is not None else None,
      'no_integral_root_in_old_cube_certified':old_cube_excluded,
      'location_proof_scope':'Only if the new candidate root exists, its certified error places coordinate588 outside the old cube, and the common cube contains both regions with derivative norm below1. At most one integral root in the common convex cube, so no integral root in the old cube. No exclusion outside the common cube or C.',
      'same_R_infinity_norm_exact':str(Rnorm),'root_certified':accepted,'strict_local_minimum_certified':accepted,
      'all17_integral_or_PDE_or_biological_validation':False,'scope':p['scope']}
 write(OUT/'results.json',out)
 print(json.dumps({'status':status,'b':float(b),'z':float(z),'eta':float(chosen['eta']),'lhs':float(chosen['selfmap_lhs']),'r':float(r),'error':float(error)if error else None,'curvature':float(curvature)if curvature else None,'results_sha256':sha(OUT/'results.json')}),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('command',choices=['prepare','run']);a=ap.parse_args();prepare()if a.command=='prepare'else run()
