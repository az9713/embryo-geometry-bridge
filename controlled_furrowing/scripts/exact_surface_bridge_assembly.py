"""All24 B-center surface discrepancies; point model only, no root transfer."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,argparse,time
from flint import arb,ctx
import exact_surface_bridge_pilot as P
ROOT=Path(__file__).resolve().parents[1]
SOURCE=Path(__file__)
OUT=ROOT/'runs/exact_integration_bridge_v1/surface_assembly_B'
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):P.write(p,d)
def unpack(d):return P.A(F(d['lower_exact'])).union(P.A(F(d['upper_exact'])))
def prepare():
    ctx.prec=512;pilot=read(P.OUT/'protocol.json');P.assert_bindings(pilot)
    review=read(P.OUT/'prelaunch_independent_review.json')
    P.validate_prelaunch(review,sha(P.__file__),sha(P.OUT/'protocol.json'),sha(P.OUT/'preflight.json'))
    bindings=dict(pilot['source_and_input_hashes'])
    for path in [SOURCE,P.OUT/'protocol.json',P.OUT/'preflight.json',P.OUT/'prelaunch_independent_review.json']:
        bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    target=P.T.Target(read(ROOT/pilot['record']['point_cache_file']));cases=[]
    for face in (0,1):
        for edge in range(12):
            field=target.surface(face,edge);domain=P.positive_G_bounds(field)
            cases.append({'face':face,'edge':edge,'domain':domain,'global_indices':field['global_indices']})
    protocol={'utc':P.stamp(),'source_and_input_hashes':bindings,'point_cache_file':pilot['record']['point_cache_file'],'record':pilot['record'],'precision_bits':512,'truncation_order':28,'cases':cases,'scope':'ALL24 surface edges at ONE saved B-center point. No parameter/root tube or equilibrium transfer.','root_transfer_proved':False}
    write(OUT/'protocol.json',protocol);write(OUT/'preflight.json',{'status':'PASS_PREPARATION_ONLY','source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json'),'checks':{'all24_declared_once':[(c['face'],c['edge']) for c in cases]==[(k,e) for k in (0,1) for e in range(12)],'all_positive_hoop_and_G':all(F(c['domain']['hoop_lower_exact'])>0 and F(c['domain']['G_lower_exact'])>0 for c in cases),'all_relative_series_domains':all(0<=F(c['domain']['rho_exact'])<1 for c in cases),'point_only_scope_explicit':True},'actual_integrals_run':False})
    print(json.dumps({'status':'PASS_PREPARATION_ONLY','source_sha256':sha(SOURCE),'protocol_sha256':sha(OUT/'protocol.json'),'preflight_sha256':sha(OUT/'preflight.json')}),flush=True)
def run():
    p=read(OUT/'protocol.json');ctx.prec=p['precision_bits']
    P.validate_prelaunch(read(OUT/'prelaunch_independent_review.json'),sha(SOURCE),sha(OUT/'protocol.json'),sha(OUT/'preflight.json'))
    for path,h in p['source_and_input_hashes'].items():
        if sha(ROOT/path)!=h:raise ValueError('changed frozen input '+path)
    target=P.T.Target(read(ROOT/p['point_cache_file']));rule=target.auxiliary_rule();g=[arb(0) for _ in range(611)];H={};files=[];start=time.monotonic()
    write(OUT/'started.json',{'utc':P.stamp(),'protocol_sha256':sha(OUT/'protocol.json'),'source_sha256':sha(SOURCE)})
    for case in p['cases']:
        face,edge=case['face'],case['edge'];result=P.run_edge(target,face,edge,rule,case['domain'],p['truncation_order']);indices=result['global_indices']
        if indices!=case['global_indices']:raise ValueError('changed frozen free map')
        for i,v in zip(indices,result['force_integral_minus_cache']):g[i]+=unpack(v)
        for item in result['H_integral_minus_cache_upper']:
            key=tuple(sorted((indices[item['i']],indices[item['j']])));H[key]=H.get(key,arb(0))+unpack(item['enclosure'])
        dest=OUT/'edges'/f'face{face}_edge{edge:02d}.json';write(dest,result);files.append({'face':face,'edge':edge,'file':dest.relative_to(ROOT).as_posix(),'sha256':sha(dest)})
        print(json.dumps({'phase':'edge_complete','face':face,'edge':edge}),flush=True)
    for path,h in p['source_and_input_hashes'].items():
        if sha(ROOT/path)!=h:raise ValueError('input changed during integration '+path)
    write(OUT/'results.json',{'status':'ALL24_POINT_SURFACE_ASSEMBLED_UNREVIEWED','protocol_sha256':sha(OUT/'protocol.json'),'source_sha256':sha(SOURCE),'edges':files,'force_integral_minus_cache':[P.core.pack(v) for v in g],'H_upper':[{'i':i,'j':j,'enclosure':P.core.pack(v)} for (i,j),v in sorted(H.items())],'all_frozen_bindings_unchanged':True,'seconds':time.monotonic()-start,'root_transfer_proved':False,'scope':p['scope']})
    print(json.dumps({'status':'ALL24_POINT_SURFACE_ASSEMBLED_UNREVIEWED','seconds':time.monotonic()-start}),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['prepare','run']);args=ap.parse_args();prepare() if args.command=='prepare' else run()
