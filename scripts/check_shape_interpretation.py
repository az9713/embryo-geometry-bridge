"""Bind the interpretation to saved results and check its algebraic illustration."""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import argparse,json
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'controlled_furrowing/research_directions/relaxed_slab_execution_v1'
def read(p):return json.loads(p.read_text(encoding='utf8'))
def digest(p):return sha256(p.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    old=read(RUN/'virtual_work_radial_2048.json')
    for row,(lo,hi) in zip(old['rows'],[('-0.001129','-0.000978'),('-0.001180','-0.001028')],strict=True):
        actual_lo,actual_hi=map(Q,row['reduced_virtual_work_exact'])
        assert Q(lo)<=actual_lo<=actual_hi<=Q(hi)<0
        assert row['component']=='radial' and row['node']==323
    r,z,c,alpha=s.symbols('r z c a',real=True)
    energy=(r-alpha)**2+(z-c*r)**2
    assert s.diff(energy,r).subs({r:0,z:0})==-2*alpha
    assert s.diff(energy,z).subs({r:0,z:0})==0
    assert energy.subs({r:alpha,z:c*alpha})==0
    hessian=s.hessian(energy,(r,z));assert s.det(hessian)==4
    assert s.expand(hessian[0,0]-2*(1+c**2))==0
    # The sum-of-squares expression and positive Hessian establish the stated unique minimizer.
    assert energy.subs(c,0).subs({r:alpha,z:0})==0
    t=s.symbols('t',real=True)
    denom=1+c**2*t*(1-t)
    rt=t*alpha/denom;zt=t**2*c*alpha/denom
    interpolated=(1-t)*(r*r+z*z)+t*energy
    assert s.simplify(s.diff(interpolated,r).subs({r:rt,z:zt}))==0
    assert s.simplify(s.diff(interpolated,z).subs({r:rt,z:zt}))==0
    assert s.diff(zt,t).subs(t,0)==0 and zt.subs(t,1)==c*alpha
    # The exclusion margin is strict in this exact arithmetic illustration.
    d,eta=Q(1,2),Q(1,16);assert d*d>eta
    audit=read(ROOT/'docs/RESULTS-LITERATURE-AUDIT.json')
    assert not audit['novelty_confirmed'] and not audit['exhaustive_search'] and not audit['author_contact']
    assert len(audit['sources'])==6
    bind=[ROOT/'docs/RESULTS-AND-NEXT-THEOREM.md',ROOT/'docs/RESULTS-AND-NEXT-THEOREM.html',
          ROOT/'docs/RESULTS-LITERATURE-AUDIT.json',ROOT/'docs/assets/shape-result-logic.svg',
          RUN/'virtual_work_radial_2048.json',RUN/'virtual_work_integer_recentered_4096.json',
          RUN/'THEOREM.md',RUN/'slab_verification.json',Path(__file__)]
    out={'status':'PASS_INTERPRETATION_IDENTITIES_AND_ILLUSTRATIVE_ALGEBRA',
      'sources':{str(p.relative_to(ROOT)).replace('\\','/'):digest(p) for p in bind},
      'strict_saved_work_cases_checked':2,'toy_radial_derivative':'-2*a','toy_height_derivative_at_old_root':'0',
      'toy_unique_minimizer':['a','c*a'],'toy_hessian_determinant':4,
      'toy_continuation_stationarity_checked':True,'toy_zero_initial_height_slope_but_nonzero_endpoint_possible':True,
      'literature_sources':6,'novelty_confirmed':False,
      'new_tissue_equilibrium_or_observable_bound_computed':False,
      'scope':'Saved-result identities, outward readable enclosures and exact algebra of the explanatory toy. The conditional exclusion proposition has a written analytic proof, not a numerical proof. Does not establish its cylinder-specific hypotheses or an external originality review.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps({k:v for k,v in out.items() if k!='sources'}))

if __name__=='__main__':main()
