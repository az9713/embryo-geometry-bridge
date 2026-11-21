from fractions import Fraction as F
from flint import arb
import time
def sample_element(I,target,index,rule,record):
    start=time.monotonic();f=target.element(index);jet=I.independent_jet(f['q'],f['B']);n=jet['local_dofs'];fun=I.Functional(target,f,index,rule);a,b=I.series(I.T.sub(jet['value'],I.T.ONE),18);vol=I.DUAL.contract_volumetric_defect(jet,a,b,fun.defect);direct_g=fun.apply(I.T.mul(a,jet['gradient'][0]));direct_h=fun.apply(I.T.add(I.T.mul(a,jet['hessian_upper'][0,0]),I.T.mul(b,I.T.mul(jet['gradient'][0],jet['gradient'][0]))));direct_checks={'gradient0_direct_moment_product_overlaps_dual':direct_g.overlaps(vol['gradient'][0]),'H00_direct_moment_product_overlaps_dual':direct_h.overlaps(vol['hessian_upper'][0,0])};assert all(direct_checks.values());drive=F(record['base_strength_exact'])*I.rational(target.a['active'][index,0]);diagonal=[F(1),1+drive,F(1),1+drive,F(1)];t1,t2=I.M.volumetric_derivative_tail_bounds(18,F(21,500),1,5000);mass=I.C.upper(fun.massI+fun.massQ);jg=[I.majorant(x,f)for x in jet['gradient']]
    force=[];hessian={}
    for i in range(n):
        p={}
        for k in range(5):p=I.T.add(p,I.T.scale(I.T.mul(f['q'][k],f['B'][k][i]),diagonal[k]))
        tail=mass*I.C.A(t1*jg[i]);force.append(fun.apply(p)+vol['gradient'][i]+(-tail).union(tail))
    for i,j in jet['hessian_upper']:
        p={}
        for k in range(5):p=I.T.add(p,I.T.scale(I.T.mul(f['B'][k][i],f['B'][k][j]),diagonal[k]))
        tail=mass*I.C.A(t2*jg[i]*jg[j]+t1*I.majorant(jet['hessian_upper'][i,j],f));hessian[i,j]=fun.apply(p)+vol['hessian_upper'][i,j]+(-tail).union(tail)
    cg,ch,count=I.correction(target,f,index,fun,drive);force=[v+cg[i]for i,v in enumerate(force)];hessian={k:v+ch[k]for k,v in hessian.items()}
    return {'element':index,'checks':direct_checks,'order':18,'precision_bits':640,'global_indices':f['global_indices'],'chart':f['chart']['kind'],'canonical_basis_entry_checks':count,'force':[I.C.pack(v)for v in force],'H_upper':[{'i':i,'j':j,'interval':I.C.pack(v)}for(i,j),v in hessian.items()],'seconds':time.monotonic()-start,'scope':'Only preselected POINT element; no global/tube/root claim.'}
