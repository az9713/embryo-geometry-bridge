"""Assemble uniform bulk integration errors on the unchanged common B box.

Reuses the seven completed local packets only when their exact selector and
domain bindings match; computes the remaining 137. No surface/root claim.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, json, time, os
from flint import arb, arb_mat, ctx
import exact_uniform_bulk_defect_B_safe as U
ROOT=U.ROOT
SOURCE=Path(__file__)
OUT=ROOT/'runs/exact_integration_bridge_v1/uniform_bulk_assembly_B_safe_v2'
PILOT=U.OUT
def read(p):return U.read(p)
def sha(p):return U.sha(p)
def write(p,d):return U.write(p,d)
def ball(p):return U.A(F(p['lower_exact'])).union(U.A(F(p['upper_exact'])))

def prepare():
    if (OUT/'protocol.json').exists():raise RuntimeError('preserve existing frozen protocol')
    ctx.prec=512;U.core.A.cache_clear()
    p=read(PILOT/'protocol.json');terminal=read(PILOT/'terminal.json');results=read(PILOT/'results.json')
    if terminal['exit_code']!=0 or terminal['results_sha256']!=sha(PILOT/'results.json'):raise ValueError('terminal pilot required')
    binds=dict(p['source_and_input_hashes'])
    for path,h in binds.items():
        if sha(ROOT/path)!=h:raise ValueError('pilot binding changed '+path)
    paths=[SOURCE,PILOT/'protocol.json',PILOT/'preflight.json',PILOT/'prelaunch_independent_review.json',PILOT/'results.json',PILOT/'terminal.json']
    for q in paths:binds[q.relative_to(ROOT).as_posix()]=sha(q)
    target=U.T.Target(read(ROOT/read(ROOT/U.P.RECORD)['point_cache_file']));rule=target.auxiliary_rule();Y=U.field_nodes(read(U.GEOM/'coordinate_boxes.json'))
    body=read(U.GEOM/'certification_run/whole_body_geometry.json')
    selectors=[]
    for i in range(144):
        jb=body['triangles'][i]['J_bound'];raw=max(abs(F(jb['lower_exact'])-1),abs(F(jb['upper_exact'])-1))
        v=raw*2**48;rho=F((v.numerator+v.denominator-1)//v.denominator,2**48)
        if not raw<=rho<1 or F(jb['lower_exact'])<=F(1,10**6):raise ValueError('fresh analytic domain required')
        selectors.append(U.choose(target,target.element(i,Y),i,rule,rho))
    reuse={}
    for item in results['elements']:
        q=PILOT/item['file'];packet=read(q);i=item['element']
        if sha(q)!=item['sha256'] or selectors[i]!=next(v for v in p['selectors']if v['element']==i):raise ValueError('pilot selector changed')
        if not packet['checks'] or not all(v is True for v in packet['checks'].values()):raise ValueError('pilot checks required')
        reuse[str(i)]={'file':q.relative_to(ROOT).as_posix(),'sha256':sha(q)};binds[q.relative_to(ROOT).as_posix()]=sha(q)
    checks={'all144_unique_selectors':len(selectors)==144 and [s['element']for s in selectors]==list(range(144)),
            'exact_seven_reused':set(map(int,reuse))=={0,2,6,8,72,74,66},
            'same_common_geometry':p['geometry_boxes_sha256']==sha(U.GEOM/'coordinate_boxes.json'),
            'same_precision':p['precision_bits']==512,
            'corrected_arithmetic_only':U.T.__name__=='exact_fe_target_safe'}
    if not all(checks.values()):raise AssertionError(checks)
    protocol={'source_and_input_hashes':binds,'precision_bits':512,'selectors':selectors,'reuse_completed_local_packets':reuse,
              'uniform_local_protocol':p,'matrix_norm_policy':'Assemble signed local interval endpoints in same 611 free coordinates; bound ||R deltaH||infinity by signed interval product and also by |R| times rowwise absolute deltaH bounds. No interval coefficient discarded.',
              'scope':'All144 bulk elements, same common nodal box, original five-dimensional C, original dose and physical target; no surfaces or root-transfer theorem.',
              'no_new_equilibrium_optimization':True,'root_transfer_proved':False}
    write(OUT/'protocol.json',protocol)
    write(OUT/'preflight.json',{'status':'PASS_PREPARATION_ONLY','checks':checks,'source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json'),
                              'orders':{str(n):sum(v['degree']==n for v in selectors)for n in sorted({v['degree']for v in selectors})},'actual_new_integrations_run':False})
    print(json.dumps({'phase':'prepared','source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json')}),flush=True)

def run():
    p=read(OUT/'protocol.json');review=read(OUT/'prelaunch_independent_review.json')
    if review.get('status')!='PASS' or not review.get('checks') or not all(v is True for v in review['checks'].values()) or review['source_sha256']!=sha(SOURCE) or review['protocol_sha256']!=sha(OUT/'protocol.json') or review['preflight_sha256']!=sha(OUT/'preflight.json'):raise ValueError('bound independent prelaunch required')
    if any((OUT/n).exists()for n in ['started.json','results.json','terminal.json']):raise RuntimeError('single launch only')
    for path,h in p['source_and_input_hashes'].items():
        if sha(ROOT/path)!=h:raise ValueError('changed frozen input '+path)
    ctx.prec=p['precision_bits'];U.core.A.cache_clear();start=time.monotonic()
    write(OUT/'started.json',{'pid':os.getpid(),'utc':U.P.now(),'protocol_sha256':sha(OUT/'protocol.json')})
    target=U.T.Target(read(ROOT/read(ROOT/U.P.RECORD)['point_cache_file']));rule=target.auxiliary_rule();Y=U.field_nodes(read(U.GEOM/'coordinate_boxes.json'))
    g=[arb(0)for _ in range(611)];H={};packets=[]
    for selector in p['selectors']:
        i=selector['element'];reused=str(i)in p['reuse_completed_local_packets']
        if reused:
            info=p['reuse_completed_local_packets'][str(i)];dest=ROOT/info['file'];packet=read(dest)
        else:
            packet=U.evaluate(target,i,Y,rule,selector,p['uniform_local_protocol']);dest=OUT/'elements'/f'{i:03d}.json';write(dest,packet)
        ids=packet['global_indices']
        for j,v in zip(ids,packet['force_integral_minus_cache']):g[j]+=ball(v)
        for v in packet['H_integral_minus_cache_upper']:
            key=tuple(sorted((ids[v['i']],ids[v['j']])));H[key]=H.get(key,arb(0))+ball(v['enclosure'])
        packets.append({'element':i,'file':dest.relative_to(ROOT).as_posix(),'sha256':sha(dest),'reused_completed_pilot':reused})
        print(json.dumps({'phase':'element_complete','element':i,'reused':reused,'seconds':time.monotonic()-start}),flush=True)
    RM=arb_mat(611,611,[U.A(float(v))for v in target.a['preconditioner'].ravel()]);HM=arb_mat(611,611);rows=[arb(0)for _ in g]
    for (i,j),v in H.items():
        HM[i,j]=v;HM[j,i]=v;rows[i]+=U.core.upper(v)
        if i!=j:rows[j]+=U.core.upper(v)
    RH=RM*HM;RF=RM*arb_mat(611,1,g)
    signed_bound=max(U.core.upper(sum((U.core.upper(RH[i,j])for j in range(611)),arb(0)))for i in range(611))
    AR=arb_mat(611,611,[U.core.upper(v)for v in RM.entries()]);rowproduct=AR*arb_mat(611,1,rows)
    rowwise_bound=max(U.core.upper(rowproduct[i,0])for i in range(611))
    for path,h in p['source_and_input_hashes'].items():
        if sha(ROOT/path)!=h:raise ValueError('input changed during run '+path)
    result={'status':'ALL144_UNIFORM_BULK_ASSEMBLED_UNREVIEWED','protocol_sha256':sha(OUT/'protocol.json'),'elements':packets,'free_dofs':611,'H_upper_slots':len(H),
            'force_integral_minus_cache':[U.core.pack(v)for v in g],
            'H_integral_minus_cache_upper':[{'i':i,'j':j,'enclosure':U.core.pack(v)}for (i,j),v in sorted(H.items())],
            'preconditioned_force_defect':[U.core.pack(RF[i,0])for i in range(611)],
            'H_absolute_row_bounds':[U.core.pack(v)for v in rows],
            'signed_RH_norm_upper':U.core.pack(signed_bound),'rowwise_absR_H_norm_upper':U.core.pack(rowwise_bound),
            'best_preconditioned_H_norm_upper':U.core.pack(min(signed_bound,rowwise_bound)),
            'seconds':time.monotonic()-start,'scope':p['scope'],'root_transfer_proved':False}
    write(OUT/'results.json',result)
    write(OUT/'terminal.json',{'status':'FINISHED','exit_code':0,'pid':os.getpid(),'utc':U.P.now(),'seconds':time.monotonic()-start,'results_sha256':sha(OUT/'results.json')})
    print(json.dumps({'phase':'finished','b_bulk':str(result['best_preconditioned_H_norm_upper']['ball']),'seconds':result['seconds']}),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['prepare','run']);a=ap.parse_args();prepare()if a.command=='prepare'else run()
