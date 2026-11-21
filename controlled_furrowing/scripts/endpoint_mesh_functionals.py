"""Target-free endpoint and ROI curve mesh sensitivity on saved equilibria only.

No forward model, optimization, fitting, or source endpoint workbook is used.
Exact rational intervals refer to supplied stored P2 coordinates only.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
from datetime import datetime,timezone
import argparse,hashlib,json,sys,time
from annulus_observer import Interval,serial
from run_activated_mesh_study import QuadraticGraph,compare_graphs
from quadratic_endpoint_observer import observe_state,observe_graphs

ROOT=Path(__file__).resolve().parents[1]
OBSERVER_SHA='075dac2148ed26a3bfe83fadf1b97602e4b08536f2be38cb95c2054776c0b926'
HELPER_SHA='78295655c236d465a4c4ac9d324c9d810d23359855e7d54761496535ffe0d8a2'
LENGTH=Q('33.5')
VARIANTS=list(product(('current_radial_line','current_disk'),('spatial_max','regional_mean')))
FILES={12:'runs/activated_quadrature/lambda_5000.0_strength_10.0/reequilibrated.json',16:'runs/joint_refinement/final.json',24:'runs/sparse_quadratic/24x12_lambda5000/gain_10.0/final.json'}

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,obj):Path(path).write_text(json.dumps(serial(obj),indent=2)+'\n',encoding='utf-8')
def signed_difference(fine,coarse):return Interval(fine.lo-coarse.hi,fine.hi-coarse.lo)
def absolute_difference(fine,coarse):
 d=signed_difference(fine,coarse)
 return Interval(Q(0) if d.lo<=0<=d.hi else min(abs(d.lo),abs(d.hi)),max(abs(d.lo),abs(d.hi)))

class RestrictedQuadraticGraph(QuadraticGraph):
 """Narrow only the declared domain; inherited exact range/inversion unchanged."""
 def __init__(self,controls,roi_um,length_unit_um=LENGTH):
  super().__init__(controls,length_unit_um)
  lo,hi=map(Q,roi_um)
  if not 0<=lo<hi or lo<self.domain[0] or hi>self.domain[1]:raise ValueError('positive-width ROI must be covered by current graph')
  self.original_domain=self.domain;self.domain=(lo,hi)


def self_test():
 checks={};graph=lambda fn:[[(t,fn(t)) for t in (Q(0),Q(1,2),Q(1))]]
 a=graph(lambda r:Q(1));b=graph(lambda r:Q(6,5));roi=(Q(1,4),Q(3,4))
 z=compare_graphs(RestrictedQuadraticGraph(a,roi,Q(1)),RestrictedQuadraticGraph(b,roi,Q(1)))
 checks['restriction_preserves_exact_constant_difference']=z['supremum_absolute_discrete_height_difference_um']==Interval(Q(1,5),Q(1,5))
 checks['restriction_reports_exact_declared_domain']=z['common_physical_radius_domain_um']==roi
 flat=[[(0,0),(Q(1,4),0),(Q(1,2),0)],[(Q(1,2),0),(Q(3,4),0),(1,0)]]
 bump=[[(0,0),(Q(1,4),0),(Q(1,2),0)],[(Q(1,2),0),(Q(3,4),1),(1,0)]]
 inside=compare_graphs(RestrictedQuadraticGraph(flat,(0,Q(2,5)),Q(1)),RestrictedQuadraticGraph(bump,(0,Q(2,5)),Q(1)))
 full=compare_graphs(QuadraticGraph(flat,Q(1)),QuadraticGraph(bump,Q(1)))
 checks['outside_ROI_shape_change_excluded']=inside['supremum_absolute_discrete_height_difference_um']==Interval(Q(0),Q(0)) and full['supremum_absolute_discrete_height_difference_um'].lo>=1
 same=RestrictedQuadraticGraph([[(0,0),(Q(1,4),Q(1,4)),(1,1)]],roi,Q(1))
 direct=RestrictedQuadraticGraph(graph(lambda r:r),roi,Q(1))
 equivalent=compare_graphs(same,direct,tolerance_um=Q('0.001'))
 checks['restriction_retains_current_radius_not_material_parameter']=equivalent['supremum_absolute_discrete_height_difference_um'].lo==0 and equivalent['supremum_absolute_discrete_height_difference_um'].hi<=Q('0.001')
 for i,r in enumerate(((0,2),(-1,Q(1,2)),(Q(1,2),Q(1,2)),(Q(3,4),Q(1,4)))):
  try:RestrictedQuadraticGraph(flat,r,Q(1));checks[f'invalid_ROI_guard_{i}']=False
  except ValueError:checks[f'invalid_ROI_guard_{i}']=True
 ea,eb=Q(1,5),Q(3,10)
 for weight,aggregation in VARIANTS:
  old=observe_graphs(graph(lambda r:10-r),graph(lambda r:0),roi_um=(0,1),reference_apical_um=10,weighting=weight,component_aggregation=aggregation)
  new=observe_graphs(graph(lambda r:10-r-ea),graph(lambda r:eb),roi_um=(0,1),reference_apical_um=10,weighting=weight,component_aggregation=aggregation)
  differences={k:signed_difference(new['four_endpoints'][k],v) for k,v in old['four_endpoints'].items()}
  for key,expected in zip(old['four_endpoints'],(ea,eb,ea+eb,ea+eb)):
   v=differences[key];checks[f'sharp_constant_Lipschitz_{weight}_{aggregation}_{key}']=v.lo<=expected<=v.hi and v.hi-v.lo<Q(1,10**15)
 checks['TV_line_disk_exact_quarter']=Q(1,2)*((Q(1,2)-Q(1,4))+(Q(1)-Q(1,4)-Q(1,2)))==Q(1,4)
 checks['TV_bound_sharp_step']=abs(Q(3,4)-Q(1,2))==Q(1,4)
 checks['TV_bound_linear_example']=abs(Q(2,3)-Q(1,2))<=Q(1,4)
 receipt={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'source_sha256':sha(__file__),'scope':'Synthetic exact ROI/inversion and sharp functional inequalities; no forward mechanics or data fit.'}
 write(ROOT/'evidence/endpoint_mesh_functionals_synthetic.json',receipt)
 print(json.dumps({'status':receipt['status'],'checks':len(checks),'failed':[k for k,v in checks.items() if not v]}),flush=True)
 if not all(checks.values()):raise RuntimeError('Retained synthetic check failure')
 return receipt


def run():
 synthetic=self_test();out=ROOT/'runs/endpoint_mesh_functionals';out.mkdir(parents=True,exist_ok=True)
 if sha(ROOT/'scripts/quadratic_endpoint_observer.py')!=OBSERVER_SHA or sha(ROOT/'scripts/run_activated_mesh_study.py')!=HELPER_SHA:raise RuntimeError('frozen observation/comparison source changed')
 states={k:json.loads((ROOT/v).read_text(encoding='utf-8')) for k,v in FILES.items()}
 first=states[12];physical=lambda s:{k:v for k,v in s['parameters'].items() if k not in ('nr','nz')}
 if any(physical(s)!=physical(first) or s['activation_strength']!=10. or s['integration_policy']!=first['integration_policy'] or not s.get('converged') for s in states.values()):raise RuntimeError('Saved models must have matching physical parameters, gain10, actual integration policy and verified status')
 roi=(Q(0),Q(first['parameters']['patch_radius'])*LENGTH)
 protocol={'UTC':datetime.now(timezone.utc).isoformat(),'scope':'No new optimization, target data or gain fit. All observer variants and all mesh comparisons retained.','input_fields':{str(k):{'file':v,'sha256':sha(ROOT/v)} for k,v in FILES.items()},'length_unit_um':LENGTH,'current_radius_ROI_um':roi,'ROI_rule':'Same fixed registered current-radius [0,initial_patch_radius*length_unit], not material transport or measured PIV registration.','references_um':{'apical':Q(first['parameters']['height'])*LENGTH,'basal':Q(0)},'variants':VARIANTS,'observer_maximum_reporting_width_um':Q('0.001'),'curve_sup_reporting_width_um':Q('0.001'),'curve_sup_maximum_cells':16384,'not_acceptance_tolerance':'Reporting resolution is not an experimental or continuum accuracy threshold.','integration_policy':first['integration_policy'],'source_sha256':sha(__file__),'observer_sha256':OBSERVER_SHA,'comparison_helper_sha256':HELPER_SHA}
 write(ROOT/'evidence/endpoint_mesh_functionals_protocol.json',protocol)
 observations={};start=time.perf_counter()
 for mesh,state in states.items():
  variants={}
  for weight,aggregation in VARIANTS:
   key=weight+'__'+aggregation
   variants[key]=observe_state(state,roi_um=roi,length_unit_um=LENGTH,weighting=weight,component_aggregation=aggregation,tolerance_um='0.001',maximum_cells=16384)
   print(f'OBSERVED {mesh} {key}',flush=True)
  observations[mesh]=variants
  write(out/f'mesh_{mesh}_observations.json',{'mesh':mesh,'source_field_sha256':sha(ROOT/FILES[mesh]),'variants':variants,'scope':protocol['scope']})
 comparisons=[]
 for coarse,fine in ((12,16),(16,24),(12,24)):
  row={'coarse_mesh':[coarse,coarse//2],'fine_mesh':[fine,fine//2],'same_current_radius_ROI_um':roi,'surfaces':{},'variants':{}}
  for surface in ('top','bottom'):
   controls='boundary_'+surface+'_quadratic_controls'
   row['surfaces'][surface]=compare_graphs(RestrictedQuadraticGraph(states[fine][controls],roi),RestrictedQuadraticGraph(states[coarse][controls],roi),tolerance_um=Q('0.001'),maximum_cells=16384)
   print(f'CURVES {coarse}->{fine} {surface}',flush=True)
  ea=row['surfaces']['top']['supremum_absolute_discrete_height_difference_um'].hi;eb=row['surfaces']['bottom']['supremum_absolute_discrete_height_difference_um'].hi
  for variant in observations[fine]:
   new=observations[fine][variant]['four_endpoints'];old=observations[coarse][variant]['four_endpoints'];diffs={}
   for key in new:
    bound=ea if key=='apical_downward' else eb if key=='basal_upward' else ea+eb
    actual=absolute_difference(new[key],old[key]);diffs[key]={'signed_fine_minus_coarse_um':signed_difference(new[key],old[key]),'absolute_difference_um':actual,'Lipschitz_upper_um':bound,'exact_interval_consistent_with_bound':actual.hi<=bound,'note':'This is a discrete-pair sensitivity bound, not PDE error.'}
   row['variants'][variant]=diffs
  comparisons.append(row)
  write(ROOT/'evidence/endpoint_mesh_functionals_results.json',{'status':'COMPLETE' if len(comparisons)==3 else 'RUNNING','protocol':protocol,'comparisons':comparisons,'observation_files':[f'runs/endpoint_mesh_functionals/mesh_{k}_observations.json' for k in states],'seconds':time.perf_counter()-start})
 full={'12_to16':'evidence/joint_refinement_comparison.json','16_to24':'evidence/sparse_quadratic_24_comparison.json'}
 checks={'all12_observer_variants_processed':sum(len(v) for v in observations.values())==12,'all_maximum_widths_achieved':all(o['maximum_shortening_certificate']['requested_width_achieved'] for variants in observations.values() for o in variants.values()),'all6_ROI_sup_widths_achieved':all(v['requested_width_achieved'] for c in comparisons for v in c['surfaces'].values()),'all48_actual_differences_below_Lipschitz_bounds':all(v['exact_interval_consistent_with_bound'] for c in comparisons for d in c['variants'].values() for v in d.values()),'all6_comparisons_use_exact_same_current_ROI':all(v['common_physical_radius_domain_um']==roi for c in comparisons for v in c['surfaces'].values()),'input_hashes_unchanged':all(sha(ROOT/v)==protocol['input_fields'][str(k)]['sha256'] for k,v in FILES.items()),'frozen_observer_and_helper_unchanged':sha(ROOT/'scripts/quadratic_endpoint_observer.py')==OBSERVER_SHA and sha(ROOT/'scripts/run_activated_mesh_study.py')==HELPER_SHA}
 write(ROOT/'evidence/endpoint_mesh_functionals_verification.json',{'status':'PASS_SCOPED_SAVED_FIELD_SENSITIVITY' if all(checks.values()) else 'FAIL','checks':checks,'synthetic_receipt_sha256':sha(ROOT/'evidence/endpoint_mesh_functionals_synthetic.json'),'results_sha256':sha(ROOT/'evidence/endpoint_mesh_functionals_results.json'),'source_sha256':sha(__file__),'full_domain_sensitivity_receipts':full,'scope':'No PDE error bound, experimental measurement uncertainty, statistical inference adequacy, original observation resolution or novel concept claim.'})
 print(json.dumps({'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'seconds':time.perf_counter()-start}),flush=True)
 if not all(checks.values()):raise RuntimeError('Retained functional-sensitivity verification failed')

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--test',action='store_true');parser.add_argument('--run',action='store_true');args=parser.parse_args()
 if args.test:self_test()
 if args.run:run()
