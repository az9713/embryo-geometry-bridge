"""Exact interval node derivatives of smooth fixed-current-ROI P2 means.

No forward optimization or maximum derivative. Full arrays represent the linear
functional on admissible trace variations; forbidden domain-endpoint radial
motions are not asserted to have independent two-sided physical partials.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json
from annulus_observer import Interval,rational,serial,intersection,_quadratic_controls,_root_bracket
from quadratic_endpoint_observer import _state_controls
ROOT=Path(__file__).resolve().parents[1]
OBSERVER_SHA='075dac2148ed26a3bfe83fadf1b97602e4b08536f2be38cb95c2054776c0b926'
BASIS=((Q(1),Q(-3),Q(2)),(Q(0),Q(4),Q(-4)),(Q(0),Q(-1),Q(2)))

def zero():return Interval(Q(0),Q(0))
def plus(a,b):return Interval(a.lo+b.lo,a.hi+b.hi)
def times(a,c):
 c=rational(c);return Interval(min(c*a.lo,c*a.hi),max(c*a.lo,c*a.hi))
def product(a,b):
 out=[Q(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return out

def horner(coeff,ab):
 value=zero()
 for c in reversed(coeff):
  terms=[value.lo*ab[0],value.lo*ab[1],value.hi*ab[0],value.hi*ab[1]]
  value=Interval(min(terms)+c,max(terms)+c)
 return value

def integral(coeff,left,right):
 primitive=[Q(0)]+[v/Q(i+1) for i,v in enumerate(coeff)];a=horner(primitive,left);b=horner(primitive,right)
 return Interval(b.lo-a.hi,b.hi-a.lo)

def mean_curve_gradient(controls,*,roi_um,weighting='current_disk',root_steps=80):
 """Physical mean Z and its per-edge [node][R,Z] gradient intervals.

 Shared-node contributions must be added; mean_state_gradients does this.
 The ROI and probability normalizer remain FIXED under variation.
 """
 if not isinstance(root_steps,int) or root_steps<1:raise ValueError('positive integer root budget required')
 lo,hi=map(rational,roi_um)
 if not 0<=lo<hi:raise ValueError('positive-width nonnegative current ROI required')
 if weighting not in ('current_radial_line','current_disk'):raise ValueError('explicit line/disk measure required')
 curves=_quadratic_controls(controls);domain=(curves[0][0][0][0],curves[-1][0][2][0])
 if lo<domain[0] or hi>domain[1]:raise ValueError('fixed physical ROI not covered by current graph')
 mass=hi-lo if weighting=='current_radial_line' else (hi*hi-lo*lo)/2
 gradient=[[[zero(),zero()] for _ in range(3)] for _ in curves];total=zero();slopes=[]
 for i,(q,c0,c1,c2) in enumerate(curves):
  cut=intersection(Interval(lo,hi),Interval(q[0][0],q[2][0]))
  if cut is None or cut.lo==cut.hi:continue
  left=_root_bracket(c0[0],c1[0],c2[0],cut.lo,root_steps);right=_root_bracket(c0[0],c1[0],c2[0],cut.hi,root_steps)
  R=[c0[0],c1[0],c2[0]];Z=[c0[1],c1[1],c2[1]];dR=[c1[0],2*c2[0]];dZ=[c1[1],2*c2[1]]
  slope=min(dR[0]+dR[1]*left[0],dR[0]+dR[1]*right[1]);slopes.append(slope)
  if slope<=0:raise ValueError('strict current radial derivative required on visible segment for smooth inverse')
  w=[Q(1)] if weighting=='current_radial_line' else R
  total=plus(total,integral(product(product(w,Z),dR),left,right))
  for j,N in enumerate(BASIS):
   gradient[i][j][0]=times(integral(product(product(w,dZ),N),left,right),-1/mass)
   gradient[i][j][1]=times(integral(product(product(w,dR),N),left,right),1/mass)
 return {'mean_height_um':times(total,1/mass),'edge_node_gradients':gradient,'current_radius_ROI_um':(lo,hi),'normalizing_mass':mass,'minimum_visible_radial_derivative':min(slopes),'root_steps':root_steps,
         'scope':'Exact rational derivative enclosure of supplied graph and fixed physical ROI/reference/weight. Boundary-domain radial directions must preserve coverage; no maximum derivatives.'}


def _boundary_node_ids(state):
 ref=[[rational(v) for v in p] for p in state['reference_nodes']];H=rational(state['parameters']['height']);edges={}
 for a,b,c,d,e,f in state['triangles']:
  for i,j,m in ((a,b,d),(b,c,e),(c,a,f)):edges[tuple(sorted((i,j)))]=m
 result={}
 for name,z in (('top',H),('bottom',Q(0))):
  selected=[]
  for (i,j),m in edges.items():
   if ref[i][1]==ref[j][1]==z:
    if ref[i][0]>ref[j][0]:i,j=j,i
    selected.append((ref[i][0],[i,m,j]))
  result[name]=[ids for r,ids in sorted(selected)]
 return result


def mean_state_gradients(state,*,roi_um=None,length_unit_um='33.5',weighting='current_disk',reference_apical_um=None,reference_basal_um=0,root_steps=80):
 """Full node gradient arrays w.r.t PHYSICAL current coordinates; no solve.

 Returned arrays [node][R,Z] have dimensionless entries. Multiply by physical
 length unit and select mechanics free DOFs for dimensionless-node derivatives.
 """
 if not state.get('converged',False):raise ValueError('verified saved state required')
 if hashlib.sha256((ROOT/'scripts/quadratic_endpoint_observer.py').read_bytes()).hexdigest()!=OBSERVER_SHA:raise RuntimeError('frozen source curve extractor changed')
 controls,H=_state_controls(state,length_unit_um);L=rational(length_unit_um)
 if roi_um is None:roi_um=(0,rational(state['parameters']['patch_radius'])*L)
 a0=H if reference_apical_um is None else rational(reference_apical_um);b0=rational(reference_basal_um)
 if a0<=b0:raise ValueError('positive reference height required')
 ids=_boundary_node_ids(state);count=len(state['current_nodes']);grads={};means={};details={}
 for surface,label,sign,reference in [('top','apical_regional_mean_downward',-1,a0),('bottom','basal_regional_mean_upward',1,-b0)]:
  observed=mean_curve_gradient(controls[surface],roi_um=roi_um,weighting=weighting,root_steps=root_steps);array=[[zero(),zero()] for _ in range(count)]
  if len(ids[surface])!=len(observed['edge_node_gradients']):raise ValueError('boundary trace node mapping mismatch')
  for edge,local in zip(ids[surface],observed['edge_node_gradients']):
   for node,gradient in zip(edge,local):
    for coordinate in (0,1):array[node][coordinate]=plus(array[node][coordinate],times(gradient[coordinate],sign))
  grads[label]=array;means[label]=plus(Interval(reference,reference),times(observed['mean_height_um'],sign));details[surface]=observed
 A='apical_regional_mean_downward';B='basal_regional_mean_upward';S='height_regional_mean_shortening'
 grads[S]=[[plus(a,b) for a,b in zip(aa,bb)] for aa,bb in zip(grads[A],grads[B])];means[S]=plus(means[A],means[B])
 return {'status':'CERTIFIED_FIXED_ROI_MEAN_DERIVATIVES','means_um':means,'node_gradients':grads,'boundary_node_ids':ids,'curve_derivative_details':details,
         'contract':{'length_unit_um':L,'current_radius_ROI_um':tuple(map(rational,roi_um)),'weighting':weighting,'reference_apical_um':a0,'reference_basal_um':b0,'derivative_coordinates':'Physical node [R,Z], so mean gradients dimensionless; multiply L0 before selecting dimensionless mechanics free coordinates.','admissibility':'Formula is a trace directional derivative on variations preserving fixed-ROI graph coverage and continuity. At a domain endpoint coinciding with an ROI endpoint, forbidden radial boundary motion is not asserted to admit an independent two-sided physical partial; select actual free DOFs.','state_scope':'One verified saved static state. No maximum, time/path, material or implicit-equilibrium derivative claimed.'}}


def midpoint_node_gradient(result,key):
 """Convenience float centres only, not a replacement for exact intervals."""
 return [[float((v.lo+v.hi)/2) for v in node] for node in result['node_gradients'][key]]


def _write(name,obj):
 (ROOT/'evidence'/name).write_text(json.dumps(serial(obj),indent=2)+'\n',encoding='utf-8')


def _numeric_mean(controls,roi,weighting):
 # Independent numeric path used ONLY for finite-difference tests: float
 # polynomial inverse and scalar SciPy quadrature, no interval primitive reuse.
 import numpy as np
 from scipy.optimize import brentq
 from scipy.integrate import quad
 numerator=0.;lo,hi=map(float,roi);disk=weighting=='current_disk'
 for q in np.array(controls,dtype=float):
  c0=q[0];c1=4*q[1]-3*q[0]-q[2];c2=2*(q[0]+q[2]-2*q[1]);left=max(lo,q[0,0]);right=min(hi,q[2,0])
  if left>=right:continue
  def value(r):
   if abs(r-q[0,0])<1e-13:t=0.
   elif abs(r-q[2,0])<1e-13:t=1.
   else:t=brentq(lambda t:c0[0]+c1[0]*t+c2[0]*t*t-r,0,1,xtol=5e-15)
   z=c0[1]+c1[1]*t+c2[1]*t*t
   return z*(r if disk else 1.)
  numerator+=quad(value,left,right,epsabs=2e-11,epsrel=2e-13)[0]
 return numerator/((hi*hi-lo*lo)/2 if disk else hi-lo)


def tests():
 import numpy as np
 checks={};rows=[];roi=(Q(5,16),Q(21,16));controls=[[(Q(0),Q(0)),(Q(3,4),Q(1,4)),(Q(2),Q(1))]]
 expected={'current_radial_line':[[Q(0),Q(1,48)],[Q(-11,24),Q(11,12)],[Q(-1,24),Q(1,16)]],
           'current_disk':[[Q(17,780),Q(-17,624)],[Q(-149,312),Q(237,260)],[Q(-109,1560),Q(361,3120)]]}
 for weight in expected:
  out=mean_curve_gradient(controls,roi_um=roi,weighting=weight)
  for j in range(3):
   for k in range(2):
    interval=out['edge_node_gradients'][0][j][k];target=expected[weight][j][k]
    checks[f'exact_R_t_plus_t2_{weight}_{j}_{k}']=interval.contains(target) and interval.hi-interval.lo<Q(1,10**20)
    epsilon=1e-5;plusc=np.array(controls,dtype=float);minusc=plusc.copy();plusc[0,j,k]+=epsilon;minusc[0,j,k]-=epsilon
    fd=(_numeric_mean(plusc,roi,weight)-_numeric_mean(minusc,roi,weight))/(2*epsilon);error=abs(fd-float(target));checks[f'FD_synthetic_{weight}_{j}_{k}']=error<2e-8;rows.append({'weight':weight,'node':j,'coordinate':k,'FD':fd,'analytic':target,'absolute_error':error})
  zsum=zero()
  for node in out['edge_node_gradients'][0]:zsum=plus(zsum,node[1])
  checks[f'vertical_translation_identity_{weight}']=zsum.contains(1)
 checks['radial_transport_term_nonzero']=expected['current_radial_line'][1][0]==Q(-11,24)
 checks['fixed_t_endpoint_error_exposed']=expected['current_radial_line'][1][0]!=Q(-1,12) and expected['current_radial_line'][1][0]+Q(3,8)==Q(-1,12)
 # Constant Z on a fixed current ROI has zero derivative under any allowed
 # radial reparameterization, even though fixed-t endpoint terms are nonzero.
 constant=[[(0,3),(Q(3,4),3),(2,3)]]
 for weight in expected:
  out=mean_curve_gradient(constant,roi_um=roi,weighting=weight)
  checks[f'constant_graph_all_radial_derivatives_zero_{weight}']=all(v[0].lo==v[0].hi==0 for v in out['edge_node_gradients'][0])
 guards=[lambda:mean_curve_gradient(controls,roi_um=(0,3)),lambda:mean_curve_gradient(controls,roi_um=(1,1)),lambda:mean_curve_gradient(controls,roi_um=roi,weighting='unspecified'),lambda:mean_curve_gradient(controls,roi_um=roi,root_steps=0),lambda:mean_curve_gradient([[[(0),(0)],[(0),(1)],[(1),(2)]]],roi_um=(0,1))]
 for i,fn in enumerate(guards):
  try:fn();checks[f'input_guard_{i}']=False
  except ValueError:checks[f'input_guard_{i}']=True
 receipt={'status':'PASS_SYNTHETIC_MEAN_DERIVATIVES' if all(checks.values()) else 'FAIL','checks':checks,'FD_rows':rows,'FD_scope':'Independent numerical check, not a proof; analytic exact coefficients and variation derivation supply the mathematical claim.','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'observer_sha256':OBSERVER_SHA}
 _write('quadratic_mean_derivatives_synthetic.json',receipt);print(json.dumps({'status':receipt['status'],'checks':len(checks),'failed':[k for k,v in checks.items() if not v]}),flush=True)
 if not all(checks.values()):raise RuntimeError('Retained derivative synthetic failure')


def saved_check():
 import numpy as np
 from quadratic_endpoint_observer import observe_state
 path=ROOT/'runs/activated_quadrature/lambda_5000.0_strength_10.0/reequilibrated.json';before=hashlib.sha256(path.read_bytes()).hexdigest();state=json.loads(path.read_text());L=Q('33.5');controls,H=_state_controls(state,L);roi=(0,rational(state['parameters']['patch_radius'])*L);ids=_boundary_node_ids(state);checks={};rows=[];outputs={};rng=np.random.default_rng(90319);epsilon=2e-5
 for weight in ('current_radial_line','current_disk'):
  result=mean_state_gradients(state,weighting=weight);outputs[weight]=result;reference=observe_state(state,weighting=weight,component_aggregation='regional_mean')
  for our,their in [('apical_regional_mean_downward','apical_downward'),('basal_regional_mean_upward','basal_upward'),('height_regional_mean_shortening','regional_mean_Eulerian_height_shortening')]:
   a=result['means_um'][our];b=reference['four_endpoints'][their];checks[f'mean_matches_frozen_observer_{weight}_{our}']=max(a.lo,b.lo)<=min(a.hi,b.hi)
  for key,gradient in result['node_gradients'].items():
   allowed=set(node for surface in ('top','bottom') for edge in ids[surface] for node in edge) if key=='height_regional_mean_shortening' else set(node for edge in ids['top' if key.startswith('apical') else 'bottom'] for node in edge)
   checks[f'nonboundary_zero_{weight}_{key}']=all(v.lo==v.hi==0 for i,node in enumerate(gradient) if i not in allowed for v in node)
  # Six independent directions, constrained at physical radial domain endpoints.
  for j in range(6):
   direction=rng.normal(size=(len(state['current_nodes']),2));direction[:,0]*=.2
   for i,p in enumerate(state['reference_nodes']):
    if p[0]==0 or p[0]==state['parameters']['outer_radius']:direction[i,0]=0.
   direction/=np.linalg.norm(direction)
   numeric={}
   for sign in (-1,1):
    vals={}
    for surface in ('top','bottom'):
     nodes=np.array(controls[surface],dtype=float);perturbed=nodes+sign*epsilon*direction[np.array(ids[surface])];vals[surface]=_numeric_mean(perturbed,roi,weight)
    numeric[sign]={'apical_regional_mean_downward':float(H)-vals['top'],'basal_regional_mean_upward':vals['bottom'],'height_regional_mean_shortening':float(H)-vals['top']+vals['bottom']}
   for key in result['node_gradients']:
    cert=zero()
    for node,d in zip(result['node_gradients'][key],direction):
     for value,coefficient in zip(node,d):cert=plus(cert,times(value,rational(float(coefficient))))
    fd=(numeric[1][key]-numeric[-1][key])/(2*epsilon);distance=max(float(cert.lo)-fd,fd-float(cert.hi),0.);checks[f'FD_saved_{weight}_{j}_{key}']=distance<2e-8
    rows.append({'weight':weight,'direction':j,'mean':key,'FD':fd,'certified_directional_derivative':cert,'distance_to_interval':distance,'step_um':epsilon})
 checks['input_source_unchanged']=hashlib.sha256(path.read_bytes()).hexdigest()==before
 receipt={'status':'PASS_SAVED12_MEAN_DERIVATIVES' if all(checks.values()) else 'FAIL','checks':checks,'input':str(path.relative_to(ROOT)),'input_sha256':before,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'observer_sha256':OBSERVER_SHA,'FD_absolute_tolerance':2e-8,'FD_rows':rows,'outputs':outputs,'scope':'One saved12HQgain10 field; no optimization, material derivative, maxima derivative or rank claim.'}
 _write('quadratic_mean_derivatives_saved.json',receipt);print(json.dumps({'status':receipt['status'],'checks':len(checks),'failed':[k for k,v in checks.items() if not v]}),flush=True)
 if not all(checks.values()):raise RuntimeError('Retained saved derivative check failure')

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--test',action='store_true');parser.add_argument('--saved-check',action='store_true');args=parser.parse_args()
 if args.test:tests()
 if args.saved_check:saved_check()
