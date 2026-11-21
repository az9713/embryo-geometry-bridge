"""Exact primal/dual certificates for a DEFINED bounded affine minimax LP.

Floating HiGHS proposes a basis. Exact rational arithmetic decides feasibility,
dual nonnegativity, stationarity and equality of objective values. This says
nothing about nonlinear/derivative error when the affine coefficients are
approximate mechanical sensitivities. Failure to recover a certificate returns
UNRESOLVED, not an accepted floating optimum.
"""
from fractions import Fraction as Q
from pathlib import Path
from itertools import combinations
import hashlib,json
import numpy as np
from scipy.optimize import linprog
from scipy.linalg import qr


def rational(value):
    if isinstance(value,Q):return value
    if isinstance(value,(float,np.floating)):
        if not np.isfinite(value):raise ValueError('finite coefficient required')
        return Q(str(float(value)))
    return Q(str(value))


def exact_solve(A,b):
    """Solve a consistent full-column-rank rational rectangular system."""
    if not A or len(A)!=len(b):raise ValueError('nonempty matching system required')
    n=len(A[0]);M=[[rational(x) for x in row]+[rational(y)] for row,y in zip(A,b)]
    if any(len(row)!=n+1 for row in M):raise ValueError('rectangular matrix required')
    pivotrow=0;pivots=[]
    for col in range(n):
        pivot=next((i for i in range(pivotrow,len(M)) if M[i][col]),None)
        if pivot is None:continue
        M[pivotrow],M[pivot]=M[pivot],M[pivotrow]
        scale=M[pivotrow][col];M[pivotrow]=[x/scale for x in M[pivotrow]]
        for i in range(len(M)):
            if i!=pivotrow and M[i][col]:
                scale=M[i][col];M[i]=[x-scale*y for x,y in zip(M[i],M[pivotrow])]
        pivots.append(col);pivotrow+=1
    if any(not any(row[:n]) and row[-1] for row in M):raise ValueError('inconsistent exact system')
    if len(pivots)!=n:raise ValueError('nonunique exact system')
    answer=[Q(0)]*n
    for i,col in enumerate(pivots):answer[col]=M[i][-1]
    return answer


def problem(residuals,jacobian,bounds):
    r=list(map(rational,residuals));J=[[rational(v) for v in row] for row in jacobian]
    bound=[tuple(map(rational,pair)) for pair in bounds];d=len(bound)
    if not r or len(r)!=len(J) or not d or any(len(row)!=d for row in J):raise ValueError('consistent nonempty dimensions required')
    if any(len(pair)!=2 or pair[0]>pair[1] for pair in bound):raise ValueError('ordered finite step bounds required')
    A=[];b=[];labels=[]
    for i,(value,row) in enumerate(zip(r,J)):
        A.extend([row+[Q(-1)],[-v for v in row]+[Q(-1)]])
        b.extend([-value,value]);labels.extend([f'residual_{i}_positive',f'residual_{i}_negative'])
    for j,(lo,hi) in enumerate(bound):
        row=[Q(0)]*(d+1);row[j]=1;A.append(row);b.append(hi);labels.append(f'coordinate_{j}_upper')
        row=[Q(0)]*(d+1);row[j]=-1;A.append(row);b.append(-lo);labels.append(f'coordinate_{j}_lower')
    A.append([Q(0)]*d+[Q(-1)]);b.append(Q(0));labels.append('nonnegative_loss')
    return r,J,bound,A,b,labels


def verify_certificate(residuals,jacobian,bounds,certificate):
    r,J,bound,A,b,labels=problem(residuals,jacobian,bounds);n=len(bound)+1
    x=list(map(rational,certificate['primal_exact']));lam=list(map(rational,certificate['dual_exact']))
    if len(x)!=n or len(lam)!=len(A):raise ValueError('certificate dimensions mismatch')
    dot=lambda a,z:sum((aa*zz for aa,zz in zip(a,z)),Q(0))
    primal=all(dot(row,x)<=rhs for row,rhs in zip(A,b))
    dual=all(v>=0 for v in lam)
    stationary=all(sum((lam[i]*A[i][j] for i in range(len(A))),Q(0))==(-1 if j==n-1 else 0) for j in range(n))
    lower=-dot(lam,b);equal=lower==x[-1]
    return {'exact_primal_feasibility':primal,'exact_dual_nonnegative':dual,
            'exact_dual_stationarity':stationary,'exact_primal_dual_objective_equality':equal}


def solve_minimax(residuals,jacobian,bounds):
    r,J,bound,A,b,labels=problem(residuals,jacobian,bounds);n=len(bound)+1
    fa=np.array(A,dtype=float);fb=np.array(b,dtype=float);cost=np.zeros(n);cost[-1]=1
    proposed=linprog(cost,A_ub=fa,b_ub=fb,bounds=[(None,None)]*n,method='highs')
    if not proposed.success:return {'status':'UNRESOLVED_FLOATING_PROPOSAL','message':proposed.message}
    slack=fb-fa@proposed.x;duals=-np.asarray(proposed.ineqlin.marginals)
    tol=1e-7*max(1.,float(np.max(abs(fb))),float(np.max(abs(fa@proposed.x))))
    active=np.flatnonzero(abs(slack)<=tol).tolist()
    support=np.flatnonzero(duals>1e-9).tolist()
    # Any full-rank active basis can propose the same optimal vertex. Prefer
    # nonzero-dual rows; fallback search is explicitly bounded and may decline.
    ranked=list(dict.fromkeys(support+active));bases=[]
    if len(ranked)>=n:
        _,_,pivots=qr(fa[ranked].T,pivoting=True,mode='economic')
        bases.append(tuple(ranked[k] for k in pivots[:n]))
    if len(support)==n:bases.insert(0,tuple(support))
    if len(active)<=18:
        for candidate in combinations(active,n):
            if len(bases)>=2000:break
            bases.append(candidate)
    seen=set();rejections=[]
    for basis in bases:
        if basis in seen:continue
        seen.add(basis)
        try:
            x=exact_solve([A[i] for i in basis],[b[i] for i in basis])
            if any(sum((a*v for a,v in zip(row,x)),Q(0))>rhs for row,rhs in zip(A,b)):
                raise ValueError('exact primal inequality fails')
            tight=[i for i,(row,rhs) in enumerate(zip(A,b)) if sum((a*v for a,v in zip(row,x)),Q(0))==rhs]
            supports=[support,list(basis)]
            supports.extend([[i] for i in tight])
            for subset in supports:
                if not subset:continue
                try:values=exact_solve([[A[i][j] for i in subset] for j in range(n)],[-int(j==n-1) for j in range(n)])
                except ValueError:continue
                if any(v<0 for v in values):continue
                lam=[Q(0)]*len(A)
                for i,value in zip(subset,values):lam[i]=value
                certificate={'primal_exact':list(map(str,x)),'dual_exact':list(map(str,lam))}
                checks=verify_certificate(r,J,bound,certificate)
                if not all(checks.values()):continue
                return {'status':'CERTIFIED_DEFINED_AFFINE_LP','step_exact':list(map(str,x[:-1])),
                    'step_display':list(map(float,x[:-1])),'loss_exact':str(x[-1]),'loss_display':float(x[-1]),
                    'certificate':certificate,'checks':checks,'active_constraints':[labels[i] for i in tight],
                    'dual_support':[{'constraint':labels[i],'multiplier_exact':str(v)} for i,v in enumerate(lam) if v],
                    'floating_proposal':{'step':proposed.x[:-1].tolist(),'loss':float(proposed.fun)},
                    'basis_trials':len(seen),'coefficient_convention':'Decimal str(float) for floating inputs; exact rationals otherwise.',
                    'scope':'Global optimum of this explicitly bounded AFFINE problem only. No rigorous mechanical derivative, trust-region approximation or nonlinear inference certificate.'}
            raise ValueError('no exact nonnegative dual recovered')
        except ValueError as error:rejections.append(str(error))
    return {'status':'UNRESOLVED_EXACT_CERTIFICATE','floating_proposal':proposed.x.tolist(),
            'basis_trials':len(seen),'rejection_reasons':sorted(set(rejections)),
            'scope':'No optimum accepted without exact certificate; numerical proposal is not promoted.'}


def tests():
    cases=[('positive_bound',[3],[[1]],[(-1,1)],Q(2)),
           ('negative_bound',[-3],[[1]],[(-1,1)],Q(2)),
           ('opposing_residuals',[2,-2],[[1],[1]],[(-3,3)],Q(2)),
           ('one_sided_tension',[1],[[-10]],[(0,'.02')],Q('0.8')),
           ('exact_zero',[1],[[1]],[(-2,2)],Q(0)),
           ('unused_coordinate',[2,-2],[[1,0],[1,0]],[(-3,3),(-4,4)],Q(2)),
           ('two_independent_channels',[3,1],[[1,0],[0,1]],[(-1,1),(-1,1)],Q(2)),
           ('degenerate_fixed_box',[3],[[1]],[(0,0)],Q(3))]
    checks={};rows=[]
    for name,r,J,bounds,known in cases:
        result=solve_minimax(r,J,bounds);checks[name+'_exact_optimum']=result.get('status')=='CERTIFIED_DEFINED_AFFINE_LP' and Q(result['loss_exact'])==known
        checks[name+'_independent_inequalities']=all(verify_certificate(r,J,bounds,result['certificate']).values()) if 'certificate' in result else False
        if 'certificate' in result:
            changed=json.loads(json.dumps(result['certificate']));changed['primal_exact'][-1]=str(known-Q(1,100))
            checks[name+'_reject_false_lower_loss']=not all(verify_certificate(r,J,bounds,changed).values())
        rows.append({'name':name,'result':result,'analytic_expected_exact':str(known)})
    for name,r,J,bounds in [('nan',[np.nan],[[1]],[(-1,1)]),('bad_bounds',[1],[[1]],[(1,0)]),('dimension',[1],[[1,2]],[(-1,1)])]:
        try:solve_minimax(r,J,bounds);checks[name+'_guard']=False
        except ValueError:checks[name+'_guard']=True
    out={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'rows':rows,
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'proof':'For Ax<=b and lambda>=0 with A^T lambda=-c, c^T x>=-lambda^T b. Equality plus primal feasibility proves global optimality.',
         'scope':'Analytically solvable affine LPs plus adversarial false-loss checks; no biological or nonlinear-model claim.'}
    path=Path(__file__).resolve().parents[1]/'evidence/bounded_affine_minimax_tests.json'
    path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'checks':len(checks)}))
    if not all(checks.values()):raise RuntimeError('affine LP checks failed')


if __name__=='__main__':tests()
