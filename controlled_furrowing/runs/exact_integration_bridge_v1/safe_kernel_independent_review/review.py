"""Independent safe polynomial-policy review; no integrals or solves."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import ast,datetime,hashlib,importlib.util,json,sys
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent
sys.path.insert(0,str(ROOT/'scripts'))
import exact_fe_target_safe as S
import exact_fe_target as OLD
TARGET_SHA='889b61cfa794d432021fbd4670b781169033d4b8451554b18e7f7a569d9d3954'
HELPER_SHA='1439e20ec127d39565467d0788173b1148b5ab47f659b3b39486b81aaf1ee373'
PROTOCOL='runs/exact_integration_bridge_v1/target_preparation_safe_v2/protocol.json'
PREFLIGHT='runs/exact_integration_bridge_v1/target_preparation_safe_v2/preflight.json'
EXPECTED_PROTOCOL='68164b3a922854b05d55c57d6872e92451e176484d7f5bf2ba3127101819850b'
EXPECTED_PREFLIGHT='0d5302e157166172f218609e008622b6238cb2d867a81add65ca9205402507cf'
HELPER='runs/exact_integration_bridge_v1/safe_kernel_v2/polynomial_jet_safe.py'
GEOMETRY='runs/exact_integration_bridge_v1/retube_B_geometry/coordinate_boxes.json'
checks={};counts={};input_hashes={}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(f):
    input_hashes[f]=sha(ROOT/f);return json.loads((ROOT/f).read_text(encoding='utf8'))
def ck(name,value):
    checks[name]=bool(value)
    if not value:raise AssertionError(name)
def reject(name,fn):
    try:fn();ck(name,False)
    except (ValueError,ZeroDivisionError,TypeError):checks[name]=True

def inside(ball,value):
    if not isinstance(ball,arb):return Q(ball)==Q(value)
    if not ball.is_finite():return False
    return Q(str(ball.lower().fmpq()))<=value<=Q(str(ball.upper().fmpq()))
def exact_zero_cross(ball):return isinstance(ball,arb)and ball.is_finite()and not ball.is_zero()and Q(str(ball.lower().fmpq()))<=0<=Q(str(ball.upper().fmpq()))
def load_module(name,file):
    spec=importlib.util.spec_from_file_location(name,ROOT/file);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def qc(p):
    # Independent polynomial oracle: RATIONAL coefficients ONLY, no Arb cleanup.
    result={}
    for k,v in p.items():
        if isinstance(v,arb):raise ValueError('exact oracle cannot accept interval coefficients')
        v=Q(v)
        if v==0:continue
        result[k]=v
    return result

def qadd(a,b):
    p=dict(a)
    for k,v in b.items():p[k]=p.get(k,Q(0))+v
    return qc(p)
def qscale(a,c):return qc({k:v*Q(c)for k,v in a.items()})
def qmul(a,b):
    p={}
    for (i,j),v in a.items():
        for (k,l),w in b.items():p[i+k,j+l]=p.get((i+k,j+l),Q(0))+v*w
    return qc(p)
def qpower(a,n):
    if n<0:
        (i,j),v=next(iter(a.items()));return {(i*n,j*n):v**n}
    p={(0,0):Q(1)}
    for _ in range(n):p=qmul(p,a)
    return p

def qdiff(a,axis):
    out={}
    for k,v in a.items():
        if k[axis]:r=list(k);r[axis]-=1;out[tuple(r)]=v*k[axis]
    return qc(out)
def qcompose(a,x,y):
    p={}
    for(i,j),v in a.items():p=qadd(p,qscale(qmul(qpower(x,i),qpower(y,j)),v))
    return p

def qeval(p,x,y):return sum((v*x**i*y**j for(i,j),v in p.items()),Q(0))

class PolynomialAD:
    def __init__(self,value,gradient,H=None):self.v=value;self.g=gradient;self.n=len(gradient);self.H=H if H is not None else{(i,j):{}for i in range(self.n)for j in range(i,self.n)}
    def __add__(self,b):return PolynomialAD(qadd(self.v,b.v),[qadd(a,c)for a,c in zip(self.g,b.g)],{k:qadd(v,b.H[k])for k,v in self.H.items()})
    def __neg__(self):return PolynomialAD(qscale(self.v,-1),[qscale(a,-1)for a in self.g],{k:qscale(v,-1)for k,v in self.H.items()})
    def __sub__(self,b):return self+-b
    def __mul__(self,b):
        return PolynomialAD(qmul(self.v,b.v),[qadd(qmul(a,b.v),qmul(self.v,c))for a,c in zip(self.g,b.g)],{(i,j):qadd(qadd(qmul(v,b.v),qmul(self.v,b.H[i,j])),qadd(qmul(self.g[i],b.g[j]),qmul(self.g[j],b.g[i])))for(i,j),v in self.H.items()})

class ScalarAD:
    def __init__(self,value,g,H=None):self.v=Q(value);self.g=list(map(Q,g));self.n=len(g);self.H=H if H is not None else{(i,j):Q(0)for i in range(self.n)for j in range(i,self.n)}
    def __add__(self,b):return ScalarAD(self.v+b.v,[a+c for a,c in zip(self.g,b.g)],{k:v+b.H[k]for k,v in self.H.items()})
    def __neg__(self):return ScalarAD(-self.v,[-a for a in self.g],{k:-v for k,v in self.H.items()})
    def __sub__(self,b):return self+-b
    def __mul__(self,b):return ScalarAD(self.v*b.v,[a*b.v+self.v*c for a,c in zip(self.g,b.g)],{(i,j):v*b.v+self.v*b.H[i,j]+self.g[i]*b.g[j]+self.g[j]*b.g[i]for(i,j),v in self.H.items()})

def semantic_tests():
    z=arb(0,S.core.A(Q(1)));small=arb(0,S.core.A(Q(1,10)))
    ck('old_bug_negative_control_reproduced',not(z!=0)and OLD.constant(z)=={})
    ck('safe_zero_crossing_is_not_exact_zero',not S.is_exact_zero(z)and(0,0)in S.constant(z))
    ck('safe_exact_zero_only_removed',S.constant(arb(0))=={}and S.constant(Q(0))=={})
    for label,value in [('nan',arb('nan')),('positive_inf',arb('inf')),('negative_inf',-arb('inf'))]:reject('nonfinite_'+label,lambda value=value:S.constant(value))
    for b in(True,False):reject('Boolean_power_rejected_'+str(b),lambda b=b:S.power(S.X,b))
    reject('fractional_power_rejected',lambda:S.power(S.X,Q(1,2)))
    for label,call in [('negative_power_singular',lambda:S.power(S.constant(z),-1)),('Laurent_evaluation_singular',lambda:S.evaluate({(-1,0):Q(1)},z,Q(0)))]:
        try:r=call();vals=list(r.values())if isinstance(r,dict)else[r];ck(label+'_never_finite_certificate',any(isinstance(v,arb)and not v.is_finite()for v in vals))
        except (ValueError,ZeroDivisionError):checks[label+'_never_finite_certificate']=True
    ai={(0,0):z,(1,0):S.core.A(Q(2)),(0,2):S.core.A(Q(-3))};bi={(0,0):S.core.A(Q(1)),(1,1):small}
    interval_ops={'constant':S.constant(z),'add':S.add(ai,bi),'multiply':S.mul(ai,bi),'power3':S.power(ai,3),'differentiate':S.diff(ai,0),'compose':S.compose(ai,S.add(S.ONE,S.X),S.sub(S.Y,S.ONE)),'shift':S.shift(ai,z,small)}
    for name,poly in interval_ops.items():ck('finite_interval_coefficients_'+name,all(not isinstance(v,arb)or v.is_finite()for v in poly.values()))
    for c,d in product((Q(-1),Q(1)),(Q(-1,10),Q(1,10))):
        ap={(0,0):c,(1,0):Q(2),(0,2):Q(-3)};bp={(0,0):Q(1),(1,1):d}
        exact_ops={'constant':qc({(0,0):c}),'add':qadd(ap,bp),'multiply':qmul(ap,bp),'power3':qpower(ap,3),'differentiate':qdiff(ap,0),'compose':qcompose(ap,{(0,0):Q(1),(1,0):Q(1)},{(0,1):Q(1),(0,0):Q(-1)}),'shift':qcompose(ap,{(1,0):Q(1),(0,0):c},{(0,1):Q(1),(0,0):d})}
        for name,p in exact_ops.items():
            outer=interval_ops[name];ck(f'interval_operation_coefficient_inclusion_{name}_{c}_{d}',all(inside(outer.get(k,arb(0)),p.get(k,Q(0)))for k in set(outer)|set(p)))
            ck(f'interval_operation_value_inclusion_{name}_{c}_{d}',inside(S.evaluate(outer,Q(2,5),Q(-1,3)),qeval(p,Q(2,5),Q(-1,3))))
    counts['semantic_interval_operations']=len(interval_ops)

def direct_geometry(corners,lam):
    x=[[Q(v)for v in p]for p in corners];det=(x[1][0]-x[0][0])*(x[2][1]-x[0][1])-(x[2][0]-x[0][0])*(x[1][1]-x[0][1]);R=sum(l*p[0]for l,p in zip(lam,x));Z=sum(l*p[1]for l,p in zip(lam,x));dl=[((x[1][1]-x[2][1])/det,(x[2][0]-x[1][0])/det),((x[2][1]-x[0][1])/det,(x[0][0]-x[2][0])/det),((x[0][1]-x[1][1])/det,(x[1][0]-x[0][0])/det)]
    N=[l*(2*l-1)for l in lam]+[4*lam[i]*lam[j]for i,j in((0,1),(1,2),(2,0))];grad=[[(4*lam[k]-1)*dl[k][j]for j in(0,1)]for k in range(3)]+[[4*(lam[k]*dl[l][j]+lam[l]*dl[k][j])for j in(0,1)]for k,l in((0,1),(1,2),(2,0))]
    B=[[Q(0)]*12 for _ in range(5)]
    for k in range(6):B[0][2*k]=grad[k][0];B[1][2*k]=grad[k][1];B[2][2*k+1]=grad[k][0];B[3][2*k+1]=grad[k][1];B[4][2*k]=N[k]/R
    return R,Z,B

def chart_point(corners,number):
    zero=[i for i,x in enumerate(corners)if x[0]==0];off=[i for i in range(3)if i not in zero];lam=[Q(0)]*3
    if not zero:lam=[Q(1,3)]*3 if number==0 else[Q(1,5),Q(2,5),Q(2,5)]
    else:
        s,t=(Q(1,3),Q(2,5))if number==0 else(Q(3,4),Q(1,4))
        if len(zero)==2:lam[off[0]]=s;lam[zero[0]]=(1-s)*(1-t);lam[zero[1]]=(1-s)*t
        else:lam[zero[0]]=1-s;lam[off[0]]=s*(1-t);lam[off[1]]=s*t
    R,Z,B=direct_geometry(corners,lam);point=(R,Z-sum(p[1]for p in corners)/3)if not zero else(s,t)
    return point,B

def interval_bulk_case(name,corners,interval_nodes,witnesses,JET,field=None):
    field=S.element_fields(corners,interval_nodes)if field is None else field;jet=JET.bulk_jacobian_jet(field['q'],field['B']);free=field['free_local'];n=len(free);cross=sum(exact_zero_cross(v)for p in field['q']+[jet['value']]for v in p.values());counts['zero_crossing_field_coefficients']=counts.get('zero_crossing_field_coefficients',0)+cross
    for pi in range(2):
        point,Bfull=chart_point(corners,pi);B=[[row[k]for k in free]for row in Bfull];qballs=[S.evaluate(p,*point)for p in field['q']];jball=S.evaluate(jet['value'],*point);gballs=[S.evaluate(p,*point)for p in jet['gradient']];Hballs={k:S.evaluate(v,*point)for k,v in jet['hessian_upper'].items()}
        for wi,nodes in enumerate(witnesses):
            Y=[v for row in nodes for v in row];values=[sum(v*y for v,y in zip(row,Y))for row in Bfull];xs=[ScalarAD(v,row)for v,row in zip(values,B)];oracle=(xs[0]*xs[3]-xs[1]*xs[2])*xs[4]
            ck(f'{name}_q_{pi}_{wi}',all(inside(a,b)for a,b in zip(qballs,values)));ck(f'{name}_J_{pi}_{wi}',inside(jball,oracle.v));ck(f'{name}_Ju_{pi}_{wi}',all(inside(v,g)for v,g in zip(gballs,oracle.g)));ck(f'{name}_Juv_{pi}_{wi}',all(inside(Hballs[k],v)for k,v in oracle.H.items()))
            counts['scalar_derivative_values_checked']=counts.get('scalar_derivative_values_checked',0)+5+1+n+n*(n+1)//2
    return cross


def point_and_helper_tests(t,old_t,JET,old_jet):
    ck('all325_exact_current_nodes_unchanged',t.nodes()==old_t.nodes())
    for i in range(144):
        f=t.element(i);old=old_t.element(i);ck('Fraction_full_element_'+str(i),f==old)
    for face in(0,1):
        for edge in range(12):ck(f'Fraction_full_surface_{face}_{edge}',t.surface(face,edge)==old_t.surface(face,edge))
    ck('all_auxiliary_points_and_weights_unchanged',t.auxiliary_rule()==old_t.auxiliary_rule())
    cases=[('bulk',i,None)for i in(0,6,72,66)]+[('surface',face,e)for face in(0,1)for e in(0,1,11)]
    for kind,i,j in cases:
        if kind=='bulk':
            f=t.element(i);q,B=f['q'],f['B'];got=JET.bulk_jacobian_jet(q,B);old=old_jet.bulk_jacobian_jet(q,B);jets=[PolynomialAD(v,row)for v,row in zip(q,B)];ad=(jets[0]*jets[3]-jets[1]*jets[2])*jets[4]
        else:
            f=t.surface(i,j);q,B=f['q_physical_R'],f['B_physical_R'];got=JET.surface_squared_jacobian_jet(q,B);old=old_jet.surface_squared_jacobian_jet(q,B);jets=[PolynomialAD(v,row)for v,row in zip(q,B)];ad=jets[0]*jets[0]*(jets[1]*jets[1]+jets[2]*jets[2])
        label=f'{kind}_{i}_{j}';ck('Fraction_helper_matches_legacy_'+label,got==old);ck('independent_exact_polynomial_AD_'+label,got['value']==ad.v and got['gradient']==ad.g and got['hessian_upper']==ad.H)
        # Uncertain zero-crossing constitutive coefficients, exact signed functional.
        A={(0,0):arb(0,S.core.A(Q(1))),(1,0):S.core.A(Q(1,2))};Bcoef={(0,0):arb(0,S.core.A(Q(1,5)))}
        def delta(a,b):return Q(2)**a*Q(1,3)**b-Q(3)**a*Q(-2,5)**b/2
        result=JET.contract_volumetric_defect(got,A,Bcoef,delta)
        def functional(poly):return qeval(poly,Q(2),Q(1,3))-qeval(poly,Q(3),Q(-2,5))/2
        for ac,bc in product((Q(-1),Q(1)),(Q(-1,5),Q(1,5))):
            ap={(0,0):ac,(1,0):Q(1,2)};bp={(0,0):bc}
            ck(f'interval_dual_gradient_{label}_{ac}_{bc}',all(inside(v,functional(qmul(ap,p)))for v,p in zip(result['gradient'],ad.g)))
            ck(f'interval_dual_H_{label}_{ac}_{bc}',all(inside(result['hessian_upper'][k],functional(qadd(qmul(ap,v),qmul(bp,qmul(ad.g[k[0]],ad.g[k[1]])))))for k,v in ad.H.items()))
    counts['Fraction_elements_compared']=144;counts['Fraction_surfaces_compared']=24;counts['independent_poly_AD_helper_cases']=10

def synthetic_field_tests(JET):
    triangles=[[[Q(1),Q(0)],[Q(2),Q(0)],[Q(1),Q(1)]],[[Q(0),Q(0)],[Q(1),Q(0)],[Q(0),Q(1)]],[[Q(0),Q(0)],[Q(1),Q(0)],[Q(1),Q(1)]]]
    for ci,corners in enumerate(triangles):
        ref=corners+[[sum(corners[j][c]for j in pair)/2 for c in(0,1)]for pair in((0,1),(1,2),(2,0))];rad=Q(1,100);bounds=[]
        for r,z in ref:bounds.append([((r-rad,r+rad)if r else(Q(0),Q(0))),(z-rad,z+rad)])
        nodes=[[S.core.A(lo).union(S.core.A(hi))for lo,hi in row]for row in bounds];witnesses=[]
        for pattern in range(8):witnesses.append([[pair[(pattern+node*3+coord*5)%2]for coord,pair in enumerate(row)]for node,row in enumerate(bounds)])
        cross=interval_bulk_case('synthetic_chart'+str(ci),corners,nodes,witnesses,JET);ck('synthetic_zero_cross_coefficients_present_'+str(ci),cross>0)
        if ci:
            identity=S.element_fields(corners,ref);ck('axis_exact_identity_endpoint_'+str(ci),S.evaluate_bulk(identity,Q(0),Q(1,3))==[1,0,0,1,1])
            # A non-axis radial node whose CURRENT interval crosses zero must not vanish.
            ri=next(i for i,p in enumerate(ref)if p[0]>0);alt=[row[:]for row in ref];alt[ri][0]=arb(0,S.core.A(rad));altpoint=[]
            for sign in(-1,1):v=[row[:]for row in ref];v[ri][0]=sign*rad;altpoint.append(v)
            interval_bulk_case('synthetic_radial_zero_cross_'+str(ci),corners,alt,altpoint,JET)
            bad=[row[:]for row in ref];bad[0][0]=arb(0,S.core.A(rad));reject('uncertain_axis_radial_constraint_rejected_'+str(ci),lambda:S.element_fields(corners,bad))
    counts['synthetic_chart_classes']=3

def actual_field_tests(t,box,JET):
    bounds=[[(Q(c)-Q(r),Q(c)+Q(r))for c,r in zip(cc,rr)]for cc,rr in zip(box['current_centers_exact'],box['coordinate_radii_exact'])]
    nodes=[[S.core.A(lo).union(S.core.A(hi))for lo,hi in row]for row in bounds]
    for element in(0,2,6,8):
        field=t.element(element,nodes);ids=t.a['triangles'][element];corners=[[S.rational(v)for v in row]for row in t.a['reference'][ids[:3]]];witnesses=[]
        for pattern in range(8):witnesses.append([[bounds[int(node)][coord][(pattern+int(node)*3+coord*5)%2]for coord in(0,1)]for node in ids])
        interval_bulk_case('actual_B_common_element'+str(element),corners,[nodes[int(node)]for node in ids],witnesses,JET,field)
    counts['actual_common_elements']=4
    # Entire-field boxes can produce zero-overlapping polynomial coefficients;
    # record actual count rather than assuming that nodal centers cross zero.
    ck('actual_and_synthetic_zero_cross_terms_retained',counts.get('zero_crossing_field_coefficients',0)>0)

def main():
    dest=ROOT/'evidence/exact_safe_polynomial_kernel_independent_review.json'
    if dest.exists()or(OUT/'results.json').exists():raise ValueError('write-once review exists')
    ctx.prec=384
    ck('safe_target_exact_frozen_hash',sha(S.__file__)==TARGET_SHA);ck('safe_helper_exact_frozen_hash',sha(ROOT/HELPER)==HELPER_SHA);ck('legacy_target_unchanged',sha(OLD.__file__)=='4fd89ddf23c59050a92f18f7b4b6e69ef8e5d08e57ba0c1949cf0b61b2c869a4')
    prep=read(PROTOCOL);preflight=read(PREFLIGHT);ck('safe_preparation_protocol_hash',input_hashes[PROTOCOL]==EXPECTED_PROTOCOL);ck('safe_preparation_preflight_hash',input_hashes[PREFLIGHT]==EXPECTED_PREFLIGHT);ck('safe_preparation66_checks_literal_true',len(preflight['checks'])==66 and all(v is True for v in preflight['checks'].values()))
    for group in(prep['source_hashes'],prep['input_hashes']):
        for f,h in group.items():ck('frozen_input_'+f,sha(ROOT/f)==h);input_hashes[f]=h
    oldprep=read('runs/exact_integration_bridge_v1/target_preparation/protocol.json');ck('same_physical_target_definition',prep['target_definition']==oldprep['target_definition']);ck('same17_parameter_cache_records',prep['records']==oldprep['records'])
    # Normalized AST comparison confirms the helper numerical formulas are unchanged.
    JET=load_module('independent_safe_JET',HELPER);oldhelper='runs/exact_integration_bridge_v1/bridge_design_probe/polynomial_jet.py';oldJET=load_module('independent_old_JET',oldhelper)
    def functions(p):return {n.name:ast.dump(n,include_attributes=False)for n in ast.parse(Path(p).read_text(encoding='utf8')).body if isinstance(n,ast.FunctionDef)}
    ck('helper_function_AST_identical_except_import_policy',functions(ROOT/HELPER)==functions(ROOT/oldhelper))
    semantic_tests();row=prep['records'][0];packet=read(row['point_cache_file']);t=S.Target(packet);old_t=OLD.Target(packet);point_and_helper_tests(t,old_t,JET,oldJET);synthetic_field_tests(JET);box=read(GEOMETRY);actual_field_tests(t,box,JET)
    ck('safe_source_still_frozen',sha(S.__file__)==TARGET_SHA and sha(ROOT/HELPER)==HELPER_SHA);ck('all_inputs_stable_at_completion',all(sha(ROOT/f)==h for f,h in input_hashes.items()))
    result={'status':'PASS','checks':checks,'source_sha256':TARGET_SHA,'helper_source_sha256':HELPER_SHA,'safe_preparation_protocol_sha256':EXPECTED_PROTOCOL,'safe_preparation_preflight_sha256':EXPECTED_PREFLIGHT,'review_source_sha256':sha(__file__),'input_hashes':input_hashes,'counts':counts,'mathematical_review':'For finite Arb coefficients, exact-zero-only deletion preserves every coefficient enclosure; Arb arithmetic and finite polynomial/Laurent operations preserve enclosure on nonsingular coordinate domains. A zero-crossing coefficient is never treated as zero. Source inspection and adversarial regression tests cover both cleanup and radial-node conditionals. Fraction point algebra and target/parameter/geometry definitions remain identical.','scope':'Safe sparse interval-polynomial component and helper ONLY. Generic rational polynomial AD/dual linearity,3synthetic charts and4actualBcommon-field samples supplement code/interval arithmetic review. No spatial integrations, solves, uniform root proofs or biological validation. Prior unsafe integration claims remain withdrawn until separately re-executed/reviewed.','limitations':['Singular inverse/Laurent operations may return nonfinite Arb rather than raising; callers must require finite results and a nonsingular domain before certification.','Field/node vertex samples are diagnostics; universal enclosure rests on interval arithmetic plus exact-zero-only coefficient preservation, not sampling alone.'],'old_integration_claims_restored':False,'new_point_integrals_computed':False,'new_root_certified':False}
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8');result['retained_results_file']=str((OUT/'results.json').relative_to(ROOT)).replace(chr(92),'/');result['retained_results_sha256']=sha(OUT/'results.json');dest.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'status':'PASS','checks':len(checks),'counts':counts,'receipt_sha256':sha(dest),'retained_results_sha256':sha(OUT/'results.json')}),flush=True)
if __name__=='__main__':main()
