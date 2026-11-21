"""Independent generic second-order polynomial AD and conditioning diagnostics."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import probe as P
T=P.T;ROOT=P.ROOT;OUT=Path(__file__).parent
class Jet:
    def __init__(self,v,g,h=None):self.v=v;self.g=g;self.n=len(g);self.h=h or {(i,j):{}for i in range(self.n)for j in range(i,self.n)}
    def __mul__(self,b):
        return Jet(T.mul(self.v,b.v),[T.add(T.mul(a,b.v),T.mul(self.v,c))for a,c in zip(self.g,b.g)],{(i,j):T.add(T.add(T.mul(self.h[i,j],b.v),T.mul(self.v,b.h[i,j])),T.add(T.mul(self.g[i],b.g[j]),T.mul(self.g[j],b.g[i])))for i,j in self.h})
    def __add__(self,b):return Jet(T.add(self.v,b.v),[T.add(a,c)for a,c in zip(self.g,b.g)],{k:T.add(self.h[k],b.h[k])for k in self.h})
    def __neg__(self):return Jet(T.scale(self.v,-1),[T.scale(a,-1)for a in self.g],{k:T.scale(v,-1)for k,v in self.h.items()})
    def __sub__(self,b):return self+-b

def check_function(q,B,constructor,manual):
    jets=[Jet(v,row)for v,row in zip(q,B)];jet=constructor(jets);val,g,h=manual(q,B)
    return jet.v==val and jet.g==g and list(jet.h.values())==h

def main():
    dest=OUT/'derivative_verification.json'
    if dest.exists():raise ValueError('write-once output exists')
    if P.sha(ROOT/'scripts/exact_fe_target.py')!=P.EXPECTED:raise ValueError('frozen target changed')
    protocol=T.read(ROOT/'runs/exact_integration_bridge_v1/target_preparation/protocol.json');row=protocol['records'][0];target=T.Target(T.read(ROOT/row['point_cache_file']));base=T.read(OUT/'results.json');checks={};conditions=[]
    for b in base['bulk']:
        f=target.element(b['element']);checks['generic_polynomial_AD_'+b['kind']]=check_function(f['q'],f['B'],lambda q:(q[0]*q[3]-q[1]*q[2])*q[4],P.jacobian_derivatives)
        if not b['kind'].startswith('off_axis'):
            J,_,_=P.jacobian_derivatives(f['q'],f['B']);z=T.sub(J,T.ONE);centered=T.compose(z,{(0,0):F(1,2),(1,0):F(1,2)},{(0,0):F(1,2),(0,1):F(1,2)})
            conditions.append({'kind':b['kind'],'unshifted_unit_chart_L1_exact':str(P.norm(z)),'centered_unit_chart_L1_exact':str(P.norm(centered)),'ratio':float(P.norm(z)/P.norm(centered))})
    for s in base['surfaces']:
        f=target.surface(s['surface'],s['edge']);checks[f"generic_polynomial_AD_surface{s['surface']}_{s['edge']}"]=check_function(f['q_physical_R'],f['B_physical_R'],lambda q:q[0]*q[0]*(q[1]*q[1]+q[2]*q[2]),P.surface_derivatives)
        G,_,_=P.surface_derivatives(f['q_physical_R'],f['B_physical_R']);z=T.sub(G,T.ONE);lo,hi=f['Rlo'],f['Rhi'];rc=(lo+hi)/2;hr=(hi-lo)/2
        physical=sum(abs(c)*(hi**i if i>=0 else lo**i)for (i,j),c in z.items())
        if lo:
            numerator=T.mul(z,T.power(T.X,2));local=P.norm(T.compose(numerator,{(0,0):rc,(1,0):hr},{}))/lo**2
        else:local=P.norm(T.compose(z,{(0,0):rc,(1,0):hr},{}))
        conditions.append({'kind':f"surface{s['surface']}_{s['edge']}",'physical_R_L1_bound_exact':str(physical),'centered_numerator_over_Rmin_power_bound_exact':str(local),'ratio':float(physical/local)if local else None})
    checks['all_10_polynomial_value_gradient_Hessian_AD_checks']=all(checks.values())
    checks['frozen_target_unchanged']=P.sha(ROOT/'scripts/exact_fe_target.py')==P.EXPECTED
    if not all(checks.values()):raise AssertionError(checks)
    result={'status':'PASS_POLYNOMIAL_DERIVATIVE_AND_CONDITIONING_PROBE','checks':checks,'conditioning':conditions,'probe_results_sha256':P.sha(OUT/'results.json'),'probe_source_sha256':P.sha(OUT/'probe.py'),'verification_source_sha256':P.sha(__file__),'frozen_target_sha256':P.EXPECTED,'scope':'Generic truncated polynomial second-order automatic differentiation reproduces all representative J/Ju/Juv and G/Gu/Guv polynomial coefficients. Conditioning bounds concern the actual B center, not its full parameter/root tube; no integration, power16 coefficient expansion or solve.'}
    dest.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8');print(json.dumps({'status':result['status'],'checks':len(checks),'sha256':P.sha(dest)}))
if __name__=='__main__':main()
