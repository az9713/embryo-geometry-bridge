"""Uniform LOCAL surface Qstar/cache force and Hessian consistency.
The old raw cache is normalized before exact coefficient subtraction. A finite
nilpotent third-order Arb jet supplies analytic derivatives, not a spatial
Taylor approximation. No whole-surface sum, integration, solve or root proof.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import factorial
import argparse,datetime,hashlib,json
import numpy as np
from flint import arb,arb_mat,ctx
import exact_fe_target as target
import algebraic_root_feasibility as core
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/exact_integration_bridge_v1/surface_cache_consistency_preparation'
TARGET_PREP='runs/exact_integration_bridge_v1/target_preparation/protocol.json'
GEOMETRY='runs/exact_integration_bridge_v1/retube_B_geometry/protocol.json'
ZERO=(0,0,0)
POWERS=tuple((i,j,k)for i in range(4)for j in range(4-i)for k in range(4-i-j))

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def A(v):
    if isinstance(v,arb):
        if not v.is_finite():raise ValueError('finite interval required')
        return v
    return core.A(target.rational(v))
def q(v):return target.rational(v)
def mat(rows):return arb_mat(len(rows),len(rows[0]),[A(v)for row in rows for v in row])
def maxabs(values):return max((core.upper(v)for v in values),default=arb(0))
def serial(v):
    if isinstance(v,F):return str(v)
    if isinstance(v,dict):return {str(k):serial(x)for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [serial(x)for x in v]
    return v

def write(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():raise ValueError('write-once output exists: '+str(p))
    p.write_text(json.dumps(serial(obj),indent=2,allow_nan=False)+'\n',encoding='utf8')
def canonical_path(path):return str(path).replace(chr(92),'/')
def bind(bindings,path,h):
    path=canonical_path(path);resolved=(ROOT/path).resolve();resolved.relative_to(ROOT.resolve());path=resolved.relative_to(ROOT.resolve()).as_posix();target.bind(bindings,path,h)
    if sha(resolved)!=h:raise ValueError('changed frozen input: '+path)
def contains_exact(ball,value):
    value=q(value);return F(str(ball.lower().fmpq()))<=value<=F(str(ball.upper().fmpq()))
def pack_matrix(m):return [[core.pack(m[i,j])for j in range(m.ncols())]for i in range(m.nrows())]

def jclean(jet):return {p:A(v)for p,v in jet.items()if not(v==0)}
def jconstant(value):return jclean({ZERO:A(value)})
def jadd(a,b):
    out=dict(a)
    for p,v in b.items():out[p]=out.get(p,arb(0))+v
    return jclean(out)
def jscale(a,s):return jclean({p:v*A(s)for p,v in a.items()})
def jmul(a,b):
    out={}
    for p,v in a.items():
        for r,w in b.items():
            exponent=tuple(i+j for i,j in zip(p,r))
            if sum(exponent)<=3:out[exponent]=out.get(exponent,arb(0))+v*w
    return jclean(out)
def jvariable(value,axis):
    if isinstance(axis,bool)or not isinstance(axis,int)or not 0<=axis<3:raise ValueError('three jet axes required')
    power=tuple(int(i==axis)for i in range(3));return jadd(jconstant(value),{power:arb(1)})
def jsqrt(jet):
    constant=jet.get(ZERO,arb(0))
    if not constant>0:raise ValueError('positive constant coefficient required by formal sqrt jet')
    # REMOVE the constant coefficient exactly; interval C-C would not be zero.
    nilpotent={p:v/constant for p,v in jet.items()if p!=ZERO};square=jmul(nilpotent,nilpotent);cube=jmul(square,nilpotent)
    normalized=jadd(jconstant(1),jadd(jscale(nilpotent,F(1,2)),jadd(jscale(square,F(-1,8)),jscale(cube,F(1,16)))))
    return jscale(normalized,constant.sqrt())

def surface_derivatives(strain,K,tension):
    if len(strain)!=3:raise ValueError('three normalized surface strains required')
    a,b,c=map(A,strain);K,tension=A(K),A(tension);norm2=b*b+c*c
    if not a>0 or not norm2>0:raise ValueError('entire surface strain box must have positive hoop and nondegenerate tangent')
    x,y,z=[jvariable(v,i)for i,v in enumerate((a,b,c))];n2=jadd(jmul(y,y),jmul(z,z))
    W=jadd(jscale(jmul(jmul(x,x),n2),K/2),jadd(jscale(jmul(x,jsqrt(n2)),tension-K),jconstant(K/2)))
    def derivative(indices):
        ex=tuple(indices.count(i)for i in range(3));return W.get(ex,arb(0))*int(np.prod([factorial(n)for n in ex]))
    g=[derivative([i])for i in range(3)];H=[[derivative([i,j])for j in range(3)]for i in range(3)];third=[[[derivative([i,j,k])for k in range(3)]for j in range(3)]for i in range(3)]
    return {'value':W.get(ZERO,arb(0)),'gradient':g,'Hessian':H,'third':third,'domain':{'hoop':a,'tangent_norm_squared':norm2,'Js':a*norm2.sqrt()},'jet_coefficients':W}

def consistency_bound(Bstar_full,Bcache_normalized_full,same_Y,selected_local,K,tension):
    Bs=list(Bstar_full);Bc=list(Bcache_normalized_full)
    if len(Bs)!=3 or len(Bc)!=3 or len(same_Y)!=6 or any(len(row)!=6 for row in Bs+Bc):raise ValueError('normalized three-by-six maps and SAME full six-coordinate Y required')
    selected=list(selected_local)
    if not selected or len(set(selected))!=len(selected)or any(isinstance(i,bool)or not isinstance(i,(int,np.integer))or not 0<=i<6 for i in selected):raise ValueError('valid unique local free columns required')
    Bs=[[q(v)for v in row]for row in Bs];Bc=[[q(v)for v in row]for row in Bc];delta=[[s-c for s,c in zip(sr,cr)]for sr,cr in zip(Bs,Bc)];Y=list(map(A,same_Y))
    qc=[sum((A(v)*y for v,y in zip(row,Y)),arb(0))for row in Bc];dq=[sum((A(v)*y for v,y in zip(row,Y)),arb(0))for row in delta];line=[v+arb(0).union(w)for v,w in zip(qc,dq)]
    cache=surface_derivatives(qc,K,tension);tube=surface_derivatives(line,K,tension);Hc=mat(cache['Hessian']);gc=arb_mat(3,1,cache['gradient']);dgradient=[sum((tube['Hessian'][i][k]*dq[k]for k in range(3)),arb(0))for i in range(3)];dH=arb_mat(3,3,[sum((tube['third'][i][j][k]*dq[k]for k in range(3)),arb(0))for i in range(3)for j in range(3)])
    star=mat([[row[k]for k in selected]for row in Bs]);old=mat([[row[k]for k in selected]for row in Bc]);dB=mat([[row[k]for k in selected]for row in delta])
    F1=dB.transpose()*gc;F2=star.transpose()*arb_mat(3,1,dgradient);H1=dB.transpose()*Hc*old;H2=star.transpose()*Hc*dB;H3=star.transpose()*dH*star
    return {'force_difference':F1+F2,'H_difference':H1+H2+H3,'force_deltaBT_gc':F1,'force_BstarT_dg':F2,'H_deltaBT_Hc_Bc':H1,'H_BstarT_Hc_deltaB':H2,'H_BstarT_dH_Bstar':H3,'delta_q':dq,'cache_q':qc,'convex_line_box':line,'delta_B_full':delta,'delta_gradient':dgradient,'delta_local_H':dH,'domain':tube['domain'],'scope':'UNWEIGHTED local point, all SAME-Y/shared-K/T values in supplied boxes. No surface sum or integral/root transfer.'}

def canonical_maps(reference,KcN,KcdN,r0c,L0c,s):
    """Pure exact map recipe; reference/current node order[start,mid,end]."""
    lo,hi=q(reference[0][0]),q(reference[2][0]);L=hi-lo;r0c,L0c,s=map(q,(r0c,L0c,s))
    if L<=0 or lo<0 or r0c<=0 or L0c<=0 or not 0<s<1:raise ValueError('positive canonical/cached geometry and strictly interior point required')
    R=lo+L*s;N=[(1-s)*(1-2*s),4*s*(1-s),s*(2*s-1)];dN=[-3+4*s,4-8*s,-1+4*s];alpha=1/(r0c*L0c);D=[alpha*L,1/L,1/L]
    star=[[F(0)]*6 for _ in range(3)];cache=[[F(0)]*6 for _ in range(3)];raw=[[F(0)]*6 for _ in range(3)]
    for j in range(3):
        star[0][2*j]=N[j]/R;star[1][2*j]=dN[j]/L;star[2][2*j+1]=dN[j]/L
        raw[0][2*j]=q(KcN[j]);raw[1][2*j]=q(KcdN[j]);raw[2][2*j+1]=q(KcdN[j])
    cache=[[v*D[i]for v in row]for i,row in enumerate(raw)]
    return {'star_full':star,'cache_normalized_full':cache,'cache_raw_full':raw,'normalization_D_exact':D,'canonical_L_exact':L,'canonical_R_at_Qstar_exact':R,'old_alpha_exact':alpha}

def actual_maps(t,face,edge,iq):
    if any(isinstance(v,bool)or not isinstance(v,int)for v in(face,edge,iq))or face not in(0,1)or not 0<=edge<12 or not 0<=iq<12:raise ValueError('valid explicit integer face/edge/quadrature indices required')
    a=t.a;ids=a[f'surface_{face}_edges'][edge];N=a[f'surface_{face}_N'][iq];dN=a[f'surface_{face}_dN'][iq];s=q(N[2])+q(N[1])/2
    maps=canonical_maps(a['reference'][ids],N,dN,a[f'surface_{face}_r0'][edge,iq],a[f'surface_{face}_l0'][edge],s);selected=[k for k in range(6)if 2*int(ids[k//2])+k%2 in t.index]
    return {**maps,'nodes':list(map(int,ids)),'selected_local':selected,'global_indices':[t.index[2*int(ids[k//2])+k%2]for k in selected],'auxiliary_normalized_weight_exact':q(a[f'surface_{face}_weights'][edge,iq])/t.energy_scale}


def synthetic_tests():
    checks={}
    def ck(k,v):
        checks[k]=bool(v)
        if not v:raise AssertionError(k)
    import sympy as sy
    a,b,c,K,tension=sy.symbols('a b c K T',real=True);ea,eb,ec=sy.symbols('ea eb ec');n2=b*b+c*c;nil=2*b*eb+2*c*ec+eb*eb+ec*ec
    sqrtjet=sy.sqrt(n2)*(1+nil/(2*n2)-nil**2/(8*n2**2)+nil**3/(16*n2**3));W=K/2*a*a*n2+(tension-K)*a*sy.sqrt(n2)+K/2;Wjet=K/2*(a+ea)**2*((b+eb)**2+(c+ec)**2)+(tension-K)*(a+ea)*sqrtjet+K/2;xs=[a,b,c];es=[ea,eb,ec];zero={e:0 for e in es}
    from itertools import combinations_with_replacement
    for order in(1,2,3):
        for indices in combinations_with_replacement(range(3),order):ck('Sympy_formal_jet_identity_'+''.join(map(str,indices)),sy.simplify(sy.diff(Wjet,*[es[i]for i in indices]).subs(zero)-sy.diff(W,*[xs[i]for i in indices]))==0)
    points=[([F(6,5),F(4,5),F(3,5)],F(3),F(7)),([F(9,10),F(-3,5),F(4,5)],F(2),F(1)),([F(1),F(1,100),F(0)],F(2),F(3))]
    for pn,(values,kk,tt)in enumerate(points):
        result=surface_derivatives(values,kk,tt);subs={x:sy.Rational(v.numerator,v.denominator)for x,v in zip(xs,values)};subs.update({K:sy.Rational(kk.numerator,kk.denominator),tension:sy.Rational(tt.numerator,tt.denominator)})
        ck('point_value_'+str(pn),contains_exact(result['value'],F(str(W.subs(subs)))))
        for i in range(3):
            ck(f'point_gradient_{pn}_{i}',contains_exact(result['gradient'][i],F(str(sy.diff(W,xs[i]).subs(subs)))))
            for j in range(3):
                ck(f'point_Hessian_{pn}_{i}{j}',contains_exact(result['Hessian'][i][j],F(str(sy.diff(W,xs[i],xs[j]).subs(subs)))))
                for k in range(3):ck(f'point_third_{pn}_{i}{j}{k}',contains_exact(result['third'][i][j][k],F(str(sy.diff(W,xs[i],xs[j],xs[k]).subs(subs)))))
    ck('zero_crossing_Arb_coefficient_NOT_dropped',ZERO in jclean({ZERO:arb(0,A(F(1,10)))}))
    v=jadd(jmul(jvariable(F(1,100),1),jvariable(F(1,100),1)),jmul(jvariable(0,2),jvariable(0,2)));sqrt=jsqrt(v);squared=jmul(sqrt,sqrt);difference=jadd(squared,jscale(v,-1))
    ck('formal_sqrt_jet_square_identity',all(x.contains(0)for x in difference.values()));ck('exactly20_degree3_monomials',len(POWERS)==20)
    ck('nilpotent_formula_requires_no_small_formal_coefficient',any(core.upper(x)>1 for p,x in v.items()if p!=ZERO)or A(2*F(1,100))/A(F(1,10000))>1)
    for name,qbox in [('negative_hoop',[-1,1,0]),('collapsed_tangent',[1,0,0]),('hoop_crosses_zero',[A(-1).union(A(1)),1,0])]:
        try:surface_derivatives(qbox,1,1);ck('analytic_guard_'+name,False)
        except ValueError:checks['analytic_guard_'+name]=True
    for bad in(True,-1,2):
        try:actual_maps(None,bad,0,0);ck('bad_face_rejected_'+str(bad),False)
        except ValueError:checks['bad_face_rejected_'+str(bad)]=True
    ref=[[F(1),F(0)],[F(2),F(0)],[F(3),F(0)]];N=[F(0),F(1),F(0)];dN=[F(-1),F(0),F(1)];matched=canonical_maps(ref,N,dN,F(2),F(2),F(1,2));center=[F(6,5),F(1,10),F(12,5),F(1,10),F(18,5),F(1,10)]
    # Old raw and normalized energies have identical SAME-Y pullback derivatives.
    raw=mat(matched['cache_raw_full']);normalized=mat(matched['cache_normalized_full']);yy=arb_mat(6,1,list(map(A,center)));qr=list((raw*yy).entries());qn=list((normalized*yy).entries());gr,Hr,_=core.surface_local(qr,A(3),A(7),A(matched['old_alpha_exact']));gn,Hn,_=core.surface_local(qn,A(3),A(7),A(1))
    Fr=raw.transpose()*arb_mat(3,1,gr);Fn=normalized.transpose()*arb_mat(3,1,gn);Hr=raw.transpose()*mat(Hr)*raw;Hn=normalized.transpose()*mat(Hn)*normalized
    ck('raw_normalized_force_pullback_identical',all(a.overlaps(b)for a,b in zip(Fr.entries(),Fn.entries())));ck('raw_normalized_H_pullback_identical',all(a.overlaps(b)for a,b in zip(Hr.entries(),Hn.entries())))
    radius=F(1,100000);Y=[A(v-radius).union(A(v+radius))for v in center];Kbox=A(1).union(A(2));Tbox=A(F(1,4)).union(A(F(1,2)));zero_bound=consistency_bound(matched['star_full'],matched['cache_normalized_full'],Y,list(range(6)),Kbox,Tbox)
    ck('matched_maps_zero_deltaq_force_H',all(v==0 for v in zero_bound['delta_q']+list(zero_bound['force_difference'].entries())+list(zero_bound['H_difference'].entries())))
    eps=F(1,2**42);changed=canonical_maps(ref,[eps,1-2*eps,eps],[-1+eps,-2*eps,1+eps],F(2)+eps,F(2)+eps/4,F(1,2));bound=consistency_bound(changed['star_full'],changed['cache_normalized_full'],Y,list(range(6)),Kbox,Tbox)
    count=0
    for signs in product((-1,1),repeat=6):
        vertex=[A(v+s*radius)for v,s in zip(center,signs)]
        for kk,tt in product((F(1),F(2)),(F(1,4),F(1,2))):
            Bs=mat(changed['star_full']);Bc=mat(changed['cache_normalized_full']);vv=arb_mat(6,1,vertex);qs=list((Bs*vv).entries());qc=list((Bc*vv).entries());gs,hs,_=core.surface_local(qs,A(kk),A(tt),A(1));gc,hc,_=core.surface_local(qc,A(kk),A(tt),A(1));directF=Bs.transpose()*arb_mat(3,1,gs)-Bc.transpose()*arb_mat(3,1,gc);directH=Bs.transpose()*mat(hs)*Bs-Bc.transpose()*mat(hc)*Bc
            ck('interval_fixture_vertex_'+str(count),all(v.contains(w)for v,w in zip(bound['force_difference'].entries(),directF.entries()))and all(v.contains(w)for v,w in zip(bound['H_difference'].entries(),directH.entries())));count+=1
    qs=[sum((A(v)*y for v,y in zip(row,Y)),arb(0))for row in changed['star_full']];qc=[sum((A(v)*y for v,y in zip(row,Y)),arb(0))for row in changed['cache_normalized_full']];naive=[s-c for s,c in zip(qs,qc)]
    ck('coefficientwise_deltaq_preserves_small_difference',maxabs(bound['delta_q'])<A(F(1,10**8))and maxabs(naive)>A(F(1,100000)))
    fixture={'node_halfwidth_exact':radius,'shared_K_exact':[1,2],'shared_T_exact':[F(1,4),F(1,2)],'binary_geometry_perturbation_exact':eps,'endpoint_witnesses':count,'coefficientwise_deltaq_max_upper':core.pack(maxabs(bound['delta_q'])),'naive_deltaq_max_upper':core.pack(maxabs(naive)),'force_difference_max_upper':core.pack(maxabs(bound['force_difference'].entries())),'H_difference_max_upper':core.pack(maxabs(bound['H_difference'].entries()))}
    return checks,fixture

def actual_fixture(t,row,box,face,edge,iq):
    maps=actual_maps(t,face,edge,iq);centers=[[F(v)for v in rr]for rr in box['current_centers_exact']];radii=[[F(v)for v in rr]for rr in box['coordinate_radii_exact']];Y=[A(centers[node][c]-radii[node][c]).union(A(centers[node][c]+radii[node][c]))for node in maps['nodes']for c in(0,1)]
    kname='basal_area_modulus'if face==0 else'apical_area_modulus';tname='basal_surface_tension'if face==0 else'apical_surface_tension';Klo,Khi=map(F,row['parameter_cell'][kname]);Tlo,Thi=map(F,row['parameter_cell'][tname])
    if Klo<=0 or Tlo<0:raise ValueError('same declared physical K/T cell required')
    Kbox=A(Klo).union(A(Khi));Tbox=A(Tlo).union(A(Thi));result=consistency_bound(maps['star_full'],maps['cache_normalized_full'],Y,maps['selected_local'],Kbox,Tbox)
    return {'status':'UNIFORM_LOCAL_SURFACE_CACHE_CONSISTENCY_UNREVIEWED','face':face,'edge':edge,'quadrature_index':iq,'global_indices':maps['global_indices'],'parameter_cell_sha256':row['parameter_cell_sha256'],'input_nodal_geometry_domain':box['domain_id'],'normalization_D_exact':maps['normalization_D_exact'],'canonical_L_exact':maps['canonical_L_exact'],'canonical_R_at_Qstar_exact':maps['canonical_R_at_Qstar_exact'],'old_alpha_exact':maps['old_alpha_exact'],'shared_K_exact':[Klo,Khi],'shared_T_exact':[Tlo,Thi],'delta_B_max_abs_exact':max(abs(v)for r in result['delta_B_full']for v in r),'delta_q':[core.pack(v)for v in result['delta_q']],'convex_line_domain':{k:core.pack(v)for k,v in result['domain'].items()},'unweighted_force_difference':pack_matrix(result['force_difference']),'unweighted_H_difference':pack_matrix(result['H_difference']),'signed_term_maxima':{k:core.pack(maxabs(result[k].entries()))for k in('force_deltaBT_gc','force_BstarT_dg','H_deltaBT_Hc_Bc','H_BstarT_Hc_deltaB','H_BstarT_dH_Bstar')},'auxiliary_normalized_weight_exact_NOT_applied':maps['auxiliary_normalized_weight_exact'],'scope':'One LOCAL auxiliary site over the NEW common geometry nodal box and SAME physical K/T cell. No whole24surface sum, exact-integral error, H assembly or root transfer.'}

def prepare():
    if(OUT/'protocol.json').exists():raise ValueError('write-once preparation exists')
    ctx.prec=256;checks,fixture=synthetic_tests();bindings={};prep=read(ROOT/TARGET_PREP);row=next(r for r in prep['records']if r['record_id']=='optoGEF2_small_col_B');geometry=read(ROOT/GEOMETRY)
    for group in(prep['source_hashes'],prep['input_hashes'],geometry['source_and_input_hashes']):
        for f,h in group.items():bind(bindings,f,h)
    for f in[TARGET_PREP,GEOMETRY,'scripts/exact_surface_cache_consistency.py','runs/exact_integration_bridge_v1/retube_B_geometry/coordinate_boxes.json','runs/exact_integration_bridge_v1/retube_B_geometry/preflight.json']:
        bind(bindings,f,sha(ROOT/f))
    box=read(ROOT/'runs/exact_integration_bridge_v1/retube_B_geometry/coordinate_boxes.json')
    if sha(ROOT/'runs/exact_integration_bridge_v1/retube_B_geometry/coordinate_boxes.json')!=geometry['prepared_box_sha256']or box['parameter_cell_sha256']!=row['parameter_cell_sha256']or geometry['parameter_cell_sha256']!=row['parameter_cell_sha256']:raise ValueError('unchanged NEWcommon geometry and SAME C required')
    packet=read(ROOT/row['point_cache_file']);t=target.Target(packet)
    if sha(ROOT/row['point_cache_file'])!=row['point_cache_sha256']:raise ValueError('old B cache changed')
    protocol={'version':1,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_and_input_hashes':bindings,'selected_record':row,'precision_bits':256,'new_common_geometry_domain_id':box['domain_id'],'new_common_geometry_box_sha256':geometry['prepared_box_sha256'],'actual_local_fixture_sites':[[0,0,0],[1,11,11]],'old_normalization':'Dc=diag(alpha_c*L_I,1/L_I,1/L_I), alpha_c=1/(r0c*L0c). qcn=Dc*qraw; old Js=qcn0*sqrt(qcn1²+qcn2²).','potential':'W=K/2*a²*(b²+c²)+(T-K)*a*sqrt(b²+c²)+K/2.','derivative_method':'Three-variable formal jet modulo totaldegree4,20coefficients. sqrt(C+N)=sqrt(C)*(1+N/(2C)-N²/(8C²)+N³/(16C³)); constant N coefficient REMOVED exactly. Finite nilpotent identity, not spatial-series convergence.','analytic_guards':'Entire q-line box: hoop>0,tangent_norm_squared>0; finite shared K/T intervals. Exact declared C enforces physical K>0,T>=0 before outward Arb enclosures.','small_difference':'deltaq=(Bstar_full-Bcn_full)*SAME full6nodeY coefficientwise. dg=int H(qline)*deltaq, dH=int Wthird(qline)*deltaq.','signed_F':'deltaBT*gc + BstarT*dg','signed_H':'deltaBT*Hc*Bc + BstarT*Hc*deltaB + BstarT*dH*Bstar','scope':'LOCAL BULK-INDEPENDENT SURFACE CACHE consistency only. Two site fixtures after protocol freeze, not whole24surface or integral-H/root job.','physical_parameter_or_geometry_change':False,'whole24_job':False,'root_transfer_claim':False}
    write(OUT/'protocol.json',protocol);write(OUT/'synthetic_tests.json',{'status':'PASS','checks':checks,'matched_and_mismatched_interval_fixture':fixture,'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json')})
    actual=[]
    for face,edge,iq in protocol['actual_local_fixture_sites']:
        result=actual_fixture(t,row,box,face,edge,iq);name=f'local_face{face}_edge{edge}_q{iq}.json';write(OUT/name,result);actual.append({'file':name,'sha256':sha(OUT/name),'status':result['status']})
    for f,h in bindings.items():
        if sha(ROOT/f)!=h:raise ValueError('bound input changed during local preparation: '+f)
    preflight={'status':'PASS_LOCAL_SURFACE_CACHE_PREPARATION_UNREVIEWED','source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'synthetic_tests_sha256':sha(OUT/'synthetic_tests.json'),'synthetic_checks':len(checks),'all_synthetic_checks_literal_true':all(v is True for v in checks.values()),'two_actual_LOCAL_newcommonbox_fixtures':actual,'all_frozen_bindings_current':True,'whole24_surface_job_run':False,'uniform_integral_H_or_root_transfer_claim':False};write(OUT/'preflight.json',preflight)
    print(json.dumps({'status':preflight['status'],'checks':len(checks),'source_sha256':sha(__file__),'protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json'),'local_actual_sites':len(actual)}),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['test','prepare']);args=ap.parse_args();ctx.prec=256
    if args.command=='prepare':prepare()
    else:
        checks,fixture=synthetic_tests();print(json.dumps({'status':'PASS','checks':len(checks),'fixture':serial(fixture)}),flush=True)
