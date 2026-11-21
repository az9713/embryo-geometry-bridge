"""Independent float arithmetic crosscheck, not a second rigorous certificate."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from scipy.integrate import quad
root=Path('controlled_furrowing');L=33.5
files={12:'runs/activated_quadrature/lambda_5000.0_strength_10.0/reequilibrated.json',16:'runs/joint_refinement/final.json',24:'runs/sparse_quadratic/24x12_lambda5000/gain_10.0/final.json'}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
class NumericCurve:
 def __init__(self,nodal):
  self.nodes=np.asarray(nodal)*L;self.coefs=np.stack((self.nodes[:,0],4*self.nodes[:,1]-3*self.nodes[:,0]-self.nodes[:,2],2*(self.nodes[:,0]+self.nodes[:,2]-2*self.nodes[:,1])),axis=1)
 def val(self,R):
  i=min(int(np.searchsorted(self.nodes[:,2,0],R,side='left')),len(self.nodes)-1);c=self.coefs[i]
  if abs(R-self.nodes[i,0,0])<1e-13:t=0.
  elif abs(R-self.nodes[i,2,0])<1e-13:t=1.
  else:t=brentq(lambda t:c[0,0]+c[1,0]*t+c[2,0]*t*t-R,0,1,xtol=1e-14)
  return c[0,1]+c[1,1]*t+c[2,1]*t*t
 def knots(self,A):return [float(v) for v in self.nodes[:,2,0] if 0<v<A]
 def extrema(self,A):
  values=[self.val(0),self.val(A)]
  for c,node in zip(self.coefs,self.nodes):
   for t in [0.,1.]+([-c[1,1]/(2*c[2,1])] if c[2,1]!=0 else []):
    if 0<=t<=1:
     r=c[0,0]+c[1,0]*t+c[2,0]*t*t
     if 0<=r<=A:values.append(c[0,1]+c[1,1]*t+c[2,1]*t*t)
  return min(values),max(values)
def numerical_max(fun,A,knots):
 points=sorted(set([0.,A]+list(knots)));best=max(fun(r) for r in points)
 for x,y in zip(points[:-1],points[1:]):
  sub=np.linspace(x,y,9)
  for a,b in zip(sub[:-1],sub[1:]):
   opt=minimize_scalar(lambda r:-fun(r),bounds=(a,b),method='bounded',options={'xatol':1e-12});best=max(best,-opt.fun)
 return float(best)
def integral(curve,A,disk):
 points=[0.]+curve.knots(A)+[A]
 return sum(quad(lambda r:curve.val(r)*(r if disk else 1),x,y,epsabs=1e-10,epsrel=2e-13)[0] for x,y in zip(points[:-1],points[1:]))/(A*A/2 if disk else A)
def checkinterval(value,bound,tolerance=2e-8):return bound['lower']['display']-tolerance<=value<=bound['upper']['display']+tolerance
checks={};numeric={};curves={};rows=[]
for mesh,rel in files.items():
 state=json.loads((root/rel).read_text());A=state['parameters']['patch_radius']*L;a=NumericCurve(state['boundary_top_quadratic_controls']);b=NumericCurve(state['boundary_bottom_quadratic_controls']);curves[mesh]=(a,b,A)
 observation=json.loads((root/f'runs/endpoint_mesh_functionals/mesh_{mesh}_observations.json').read_text())['variants'];M=numerical_max(lambda r:L-a.val(r)+b.val(r),A,a.knots(A)+b.knots(A));maxA=L-a.extrema(A)[0];maxB=b.extrema(A)[1]
 for variant,result in observation.items():
  weight,aggregation=variant.split('__');disk=weight=='current_disk';meanA=L-integral(a,A,disk);meanB=integral(b,A,disk);values=[maxA,maxB,M,meanA+meanB] if aggregation=='spatial_max' else [meanA,meanB,M,meanA+meanB]
  for (key,bound),value in zip(result['four_endpoints'].items(),values):
   check=f'{mesh}_{variant}_{key}';checks[check]=checkinterval(value,bound);rows.append({'check':check,'independent_value_um':value,'production_interval_um':[bound['lower']['display'],bound['upper']['display']],'within_interval_up_to_float_tolerance':checks[check]})
results=json.loads((root/'evidence/endpoint_mesh_functionals_results.json').read_text());surfacechecks=[]
for c in results['comparisons']:
 coarse=c['coarse_mesh'][0];fine=c['fine_mesh'][0]
 for i,surface in enumerate(('top','bottom')):
  a=curves[fine][i];b=curves[coarse][i];A=curves[fine][2];maximum=numerical_max(lambda r:abs(a.val(r)-b.val(r)),A,a.knots(A)+b.knots(A));bound=c['surfaces'][surface]['supremum_absolute_discrete_height_difference_um'];name=f'ROI_sup_{coarse}_{fine}_{surface}';checks[name]=checkinterval(maximum,bound);surfacechecks.append({'check':name,'independent_numerical_maximum_um':maximum,'production_interval_um':[bound['lower']['display'],bound['upper']['display']],'within_enclosure_with_float_tolerance':checks[name]})
receipt={'status':'PASS_INDEPENDENT_NUMERICAL_CROSSCHECK' if all(checks.values()) else 'FAIL','checks':checks,'endpoint_checks':rows,'ROI_curve_checks':surfacechecks,'independent_methods':'Separately coded physical nodal polynomials; Brent inverse, SciPy quadrature and bounded float optimization. Does not call exact observer or graph helper.','float_comparison_tolerance_um':2e-8,'source_sha256':sha(Path(__file__)),'inputs_sha256':{rel:sha(root/rel) for rel in files.values()},'scope':'Independent numerical sanity check on exact discrete certificates, not an independent rigorous error enclosure or experiment.'};(root/'evidence/endpoint_mesh_functionals_independent.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(json.dumps({'status':receipt['status'],'checks':len(checks),'failed':[k for k,v in checks.items() if not v]},indent=2))
if not all(checks.values()):raise RuntimeError('Retained independent numerical discrepancy')
