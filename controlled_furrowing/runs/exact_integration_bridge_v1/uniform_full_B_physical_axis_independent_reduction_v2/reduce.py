"""Exact uniform168-packet reduction; no local reintegration or root solve."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,base64,datetime,gzip,hashlib,io,json,os,time,traceback
import numpy as np
from flint import arb,arb_mat,fmpq,ctx
import algebra as A
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
BO=ROOT/'runs/exact_integration_bridge_v1/uniform_bulk_assembly_B_safe_v2'
SO=ROOT/'runs/exact_integration_bridge_v1/uniform_surface_B_safe_v2'
POB=ROOT/'runs/exact_integration_bridge_v1/bulk_assembly_B_safe_v2'
POS=ROOT/'runs/exact_integration_bridge_v1/surface_assembly_B_safe_v2'
AX=ROOT/'runs/exact_integration_bridge_v1/uniform_axis_edge_B_safe_v2'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):
    if p.exists():raise ValueError('preserve '+str(p))
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def exact_float(v):return Q.from_float(float(v))
def assert_inputs(p):
    for k,h in p['input_hashes'].items():
        if sha(ROOT/k)!=h:raise ValueError('frozen input changed '+k)
def bound(inputs,p,expected=None):
    p=Path(p);h=sha(p)
    if expected is not None and h!=expected:raise ValueError('authoritative hash mismatch '+str(p))
    key=p.resolve().relative_to(ROOT.resolve()).as_posix()
    if key in inputs and inputs[key]!=h:raise ValueError('conflicting binding '+key)
    inputs[key]=h;return read(p)if p.suffix=='.json'else None
def prepare(bulksha,surfsha):
    if not all(x and len(x)==64 and all(c in '0123456789abcdef'for c in x)for x in [bulksha,surfsha]):raise ValueError('two full authoritative terminal hashes required')
    inputs={};bulk=bound(inputs,BO/'results.json',bulksha);surface=bound(inputs,SO/'results.json',surfsha);bp=bound(inputs,BO/'protocol.json');sp=bound(inputs,SO/'protocol.json')
    if bulk['status']!='ALL144_UNIFORM_BULK_ASSEMBLED_UNREVIEWED'or surface['status']!='UNIFORM_ALL24_SURFACE_COMPONENT_ASSEMBLED_UNREVIEWED'or len(bulk['elements'])!=144 or len(surface['edges'])!=24:raise ValueError('true complete producer packets required')
    for folder,h in [(BO,bulksha),(SO,surfsha)]:
        term=bound(inputs,folder/'launch_terminal.json')
        if term['status']!='FINISHED'or term['exit_code']!=0 or term['results_sha256']!=h or term['all_launch_bindings_current']is not True:raise ValueError('actual launcher terminal evidence required')
        review=bound(inputs,folder/'prelaunch_independent_review.json');pre=bound(inputs,folder/'preflight.json')
        if review['status']!='PASS'or not review['checks']or not all(v is True for v in review['checks'].values())or review['protocol_sha256']!=sha(folder/'protocol.json')or review['preflight_sha256']!=sha(folder/'preflight.json'):raise ValueError('bound complete prelaunch review required')
        if review['source_sha256']!=pre['source_sha256']:raise ValueError('prelaunch source mismatch')
    for p in [bp,sp]:
        for k,h in p['source_and_input_hashes'].items():
            path=(ROOT/k).resolve();key=path.relative_to(ROOT.resolve()).as_posix()
            if sha(path)!=h or(key in inputs and inputs[key]!=h):raise ValueError('changed/conflicting source binding '+key)
            inputs[key]=h
    if any(sha(ROOT/k)!=h for k,h in inputs.items()):raise ValueError('source inputs changed')
    for manifest in [bulk['elements'],surface['edges']]:
        for x in manifest:bound(inputs,ROOT/x['file'],x['sha256'])
    for folder in [POB,POS]:
        m=bound(inputs,folder/'results.json')
        for x in m['elements'if folder==POB else'edges']:bound(inputs,ROOT/x['file'],x['sha256'])
    for path in [Path(__file__),HERE/'algebra.py',HERE/'selfcheck.py',HERE/'selfchecks.json',ROOT/'evidence/exact_complete_point_B_independent_review_safe_v2.json',ROOT/'runs/exact_integration_bridge_v1/cached_Q_common_B_v1/results.json',ROOT/'runs/exact_integration_bridge_v1/retube_B_geometry/protocol.json',ROOT/'runs/exact_integration_bridge_v1/retube_B_geometry/coordinate_boxes.json']:
        bound(inputs,path)
    if sha(ROOT/'evidence/exact_complete_point_B_independent_review_safe_v2.json')!='433125c73e1ee6856f425ab7dc1d62658d9505a7c7ff8968c5ad2d4f7363dfc2':raise ValueError('corrected SAFE point review required')
    row=sp['selected_record'];ulp=bp['uniform_local_protocol'];box=read(ROOT/sp['common_geometry_box_file']);geom=read(ROOT/'runs/exact_integration_bridge_v1/retube_B_geometry/protocol.json')
    if ulp['parameter_cell_sha256']!=sp['parameter_cell_sha256']or ulp['parameter_cell']!=row['parameter_cell']or ulp['geometry_domain_id']!=sp['common_geometry_domain_id']or ulp['geometry_boxes_sha256']!=sp['common_geometry_box_sha256']or sha(ROOT/sp['common_geometry_box_file'])!=sp['common_geometry_box_sha256']or box['root_membership_assumed']is not False:raise ValueError('same explicit common C/geometry family required')
    ax=bound(inputs,AX/'results.json','06b6bed9eedcb17ae2f3abf1e2097489a8f8f74a6486deb9ebd9152c75f6ea9c');ap=bound(inputs,AX/'protocol.json',ax['protocol_sha256']);execution=bound(inputs,AX/'execution_receipt.json');ar=bound(inputs,AX/'prelaunch_independent_review.json');representation=bound(inputs,ROOT/'evidence/exact_axis_edge_physical_independent_review.json','8eaa76bd6ee3a6980815e303aba619c3f21243191c44035e0b1c3dabb26edffd')
    if execution['actual_unified_exec_exit_code']!=0 or execution['results_sha256']!=sha(AX/'results.json')or execution['independent_prelaunch_review_sha256']!=sha(AX/'prelaunch_independent_review.json')or ar['status']!='PASS'or not all(v is True for v in ar['checks'].values())or ar['source_sha256']!=execution['source_sha256']or ar['protocol_sha256']!=sha(AX/'protocol.json')or representation['status']!='PASS_INDEPENDENT_AXIS_EDGE_PHYSICAL_REPRESENTATION'or not all(v is True for v in representation['checks'].values()):raise ValueError('actual reviewed physical-six execution and representation required')
    if ax['status']!='SIX_SAME_TARGET_PHYSICAL_AXIS_EDGE_COMPONENTS_UNREVIEWED'or not ax['all_point_diagnostics_pass']or {x['element']for x in ax['elements']}!=set(range(72,78))or ap['parameter_cell']!=row['parameter_cell']or ap['parameter_cell_sha256']!=sp['parameter_cell_sha256']or ap['geometry_domain_id']!=sp['common_geometry_domain_id']or ap['geometry_boxes_sha256']!=sp['common_geometry_box_sha256']or ap['dose_exact']!=row['dose_exact']:raise ValueError('same complete physical six/common geometry/C/dose required')
    for name,h in ap['source_and_input_hashes'].items():
        path=(ROOT/name).resolve();key=path.relative_to(ROOT.resolve()).as_posix()
        if sha(path)!=h or(key in inputs and inputs[key]!=h):raise ValueError('changed axis source/input '+key)
        inputs[key]=h
    axis_items={x['element']:x for x in ax['elements']};effective=[]
    for item in bulk['elements']:
        i=item['element']
        if i in axis_items:
            new=axis_items[i];file=AX/new['file'];bound(inputs,file,new['sha256']);effective.append({'element':i,'file':file.relative_to(ROOT).as_posix(),'sha256':new['sha256'],'representation_override':'axis_edge_physical','original_packet_file':item['file'],'original_packet_sha256':item['sha256']})
        else:effective.append(dict(item))
    bound(inputs,ROOT/'evidence/exact_uniform_full_B_independent_review.json','5e51264b47bf9ae407021dc429555ca02fceb54efc023dba18f71f022f8b5b76')
    bound(inputs,ROOT/row['point_cache_file'],row['point_cache_sha256']);tests=read(HERE/'selfchecks.json')
    if tests['status']!='PASS_SYNTHETIC_REDUCER_ONLY'or not all(v is True for v in tests['checks'].values())or tests['algebra_sha256']!=sha(HERE/'algebra.py'):raise ValueError('current synthetic reducer tests required')
    p={'status':'FROZEN_TRUE_TERMINAL_UNIFORM_REDUCTION','utc':now(),'input_hashes':inputs,'bulk_results_sha256':bulksha,'surface_results_sha256':surfsha,'source_sha256':sha(__file__),'selected_record':row,'parameter_cell_sha256':sp['parameter_cell_sha256'],'common_geometry_domain_id':sp['common_geometry_domain_id'],'common_radius_exact':geom['new_common_enclosing_radius_exact'],'signed_product_precision_bits':256,'methods':['exact_rational_absR_times_H_absolute_rows','independently_assembled_signed_Arb256_RH_with_exact_endpoint_row_sums'],'reused_pilot_policy':'Exactly seven hash-bound completed safe pilot packets, no substitution by new physical-axis representations in this version.','scope':'Uniform finite-P2 exact-integral-minus-cached-Q component, original B/common explicit nodal box and SAME C. No independent local reintegration here, no root/selfmap/PDE/biology or full-domain inference claim.'}
    p['effective_bulk_packets']=effective;p['axis_protocol_file']=(AX/'protocol.json').relative_to(ROOT).as_posix();p['axis_results_sha256']=sha(AX/'results.json');p['packet_selection']='138 original bulk +6 same-target physical-axis packets72..77 +24 original surface; original first bound preserved';p['reused_pilot_policy']='Exactly the remaining original safe pilot records outside72..77 remain reused; all six full-axis-edge records replaced as prescribed, without outcome-dependent selection.'
    assert_inputs(p);write(HERE/'protocol.json',p);print(json.dumps({'status':p['status'],'protocol_sha256':sha(HERE/'protocol.json'),'inputs':len(inputs),'source_sha256':sha(__file__)}),flush=True)
def run():
    p=read(HERE/'protocol.json');assert_inputs(p)
    if p['source_sha256']!=sha(__file__)or(HERE/'started.json').exists():raise ValueError('frozen single reduction only')
    start=time.monotonic();write(HERE/'started.json',{'pid':os.getpid(),'utc':now(),'protocol_sha256':sha(HERE/'protocol.json')});checks={}
    bulk=read(BO/'results.json');surface=read(SO/'results.json');bp=read(BO/'protocol.json');sp=read(SO/'protocol.json');ulp=bp['uniform_local_protocol'];physical_selectors={x['element']:x for x in read(ROOT/p['axis_protocol_file'])['selectors']};row=p['selected_record'];cache=read(ROOT/row['point_cache_file']);raw=base64.b64decode(cache['frozen_binary_cache']['data'],validate=True)
    if hashlib.sha256(raw).hexdigest()!=cache['frozen_binary_cache']['sha256']:raise ValueError('embedded binary cache checksum mismatch')
    with np.load(io.BytesIO(raw),allow_pickle=False)as z:arrays={k:z[k].copy()for k in z.files}
    R=[[exact_float(v)for v in rr]for rr in arrays['preconditioner']];n=611;assert len(R)==n and all(len(x)==n for x in R)
    F=[(Q(0),Q(0))for _ in range(n)];H={};details=[];pointoverlaps=[];freeindex={int(k):i for i,k in enumerate(arrays['free'])};body=read(ROOT/'runs/exact_integration_bridge_v1/retube_B_geometry/certification_run/whole_body_geometry.json')
    for kind,items,output in [('bulk',p['effective_bulk_packets'],bulk),('surface',surface['edges'],surface)]:
        partF=[(Q(0),Q(0))for _ in range(n)];partH={};seen=set();selectors={x['element']:x for x in bp['selectors']}if kind=='bulk'else{(x['face'],x['edge']):x for x in sp['cases']}
        for item in items:
            local=read(ROOT/item['file']);tag=local['element']if kind=='bulk'else(local['face'],local['edge']);assert tag not in seen;seen.add(tag);ids=local['global_indices'];selected=physical_selectors[tag]if kind=='bulk'and tag in physical_selectors else selectors[tag];nn=len(ids);assert ids==selected.get('global_indices',ids);hkey='H_integral_minus_cache_upper';hh={(x['i'],x['j']):A.interval(x['enclosure'])for x in local[hkey]}
            nodes=arrays['triangles'][tag]if kind=='bulk'else arrays[f'surface_{tag[0]}_edges'][tag[1]];dofs=[2*int(node)+c for node in nodes for c in [0,1]];assert ids==[freeindex[k]for k in dofs if k in freeindex]
            if len(hh)!=len(local[hkey]):raise ValueError('duplicate local Hessian entries')
            ff=list(map(A.interval,local['force_integral_minus_cache']));A.scatter(ids,ff,hh,F,H);A.scatter(ids,ff,hh,partF,partH)
            if kind=='bulk':
                good=local['degree']==selected['degree']and Q(local['rho_exact'])==Q(selected['rho_exact'])and bool(local.get('checks',local.get('point_overlap_diagnostics_only',{})))and all(v is True for v in local.get('checks',local.get('point_overlap_diagnostics_only',{})).values())and local['cache_consistency']['all64_convex_strain_lines_positive']is True and Q(local['cache_consistency']['minimum_cache_convex_line_J_lower']['lower_exact'])>Q(1,10**6)and all(A.interval(v)[1]<=Q(ulp['force_tail_target_exact'])for v in local['force_two_mass_tail_upper'])and all(A.interval(v['two_mass_tail_upper'])[1]<=Q(ulp['H_entry_tail_target_exact'])for v in local[hkey])
                jb=body['triangles'][tag]['J_bound'];rawrho=max(abs(Q(jb['lower_exact'])-1),abs(Q(jb['upper_exact'])-1));scaled=rawrho*2**48;dyadic=Q((scaled.numerator+scaled.denominator-1)//scaled.denominator,2**48);good=good and Q(local['rho_exact'])==dyadic and Q(jb['lower_exact'])>Q(1,10**6)and local['active']==str(exact_float(arrays['active'][tag,0]))
                if tag in physical_selectors:good=good and item['representation_override']=='axis_edge_physical'and local['chart']=='axis_edge_physical'
                else:
                    reused=str(tag)in bp['reuse_completed_local_packets'];good=good and item['reused_completed_pilot']is reused
                    if reused:good=good and item['file']==bp['reuse_completed_local_packets'][str(tag)]['file']and item['sha256']==bp['reuse_completed_local_packets'][str(tag)]['sha256']
                point=read(POB/'elements'/f'{tag:03d}.json');phkey='H_integral_minus_cache_upper_triangle'
            else:
                dom=selected['domain'];sites=local['actual_cache_sites'];good=local['order']==28 and local['precision_bits']==640 and local['q0_exact']==dom['q0_exact']and local['rho_exact']==dom['rho_exact']and Q(dom['G_bounds_exact'][0])>0 and Q(dom['hoop_lower_exact'])>0 and 0<=Q(dom['rho_exact'])<1 and len(sites)==12 and [x['quadrature_index']for x in sites]==list(range(12))and all(Q(x['weight_applied_exactly_once'])==Q(y['auxiliary_normalized_weight_exact'])>0 for x,y in zip(sites,selected['cache_site_domains']));gl,gh=map(Q,dom['G_bounds_exact']);good=good and Q(dom['q0_exact'])==(gl+gh)/2 and Q(dom['rho_exact'])==(gh-gl)/(gh+gl)and dom['status']=='CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN'
                point=read(POS/'edges'/f'face{tag[0]}_edge{tag[1]:02d}.json');phkey='H_integral_minus_cache_upper'
            checks[f'{kind}_{tag}_local_map_domain_tail_reuse_policy']=bool(good)
            ph={(x['i'],x['j']):A.interval(x['enclosure'])for x in point[phkey]};pf=list(map(A.interval,point['force_integral_minus_cache']));meeting=ids==point['global_indices']and len(pf)==nn and all(A.overlap(a,b)for a,b in zip(ff,pf))and set(hh)==set(ph)and all(A.overlap(v,ph[k])for k,v in hh.items());checks[f'{kind}_{tag}_safe_point_packet_overlap']=meeting
            pointoverlaps.append({'kind':kind,'tag':str(tag),'all_point_intervals_contained':all(A.contains(a,b)for a,b in zip(ff,pf))and all(A.contains(v,ph[k])for k,v in hh.items()),'all_point_intervals_overlap':meeting})
            details.append({'kind':kind,'tag':str(tag),'file':item['file'],'sha256':item['sha256']})
        checks[kind+'_complete_unique_case_set']=seen==set(range(144))if kind=='bulk'else seen=={(f,e)for f in [0,1]for e in range(12)}
        if kind=='surface':
            checks[kind+'_all611_force_exact_reduction']=len(output['force_integral_minus_cache'])==n and all(A.contains(A.interval(x),v)for x,v in zip(output['force_integral_minus_cache'],partF))
            oh={(x['i'],x['j']):A.interval(x['enclosure'])for x in output['H_upper']};checks[kind+'_all_symmetric_H_exact_reduction']=len(oh)==len(output['H_upper'])and set(oh)==set(partH)and all(A.contains(oh[k],v)for k,v in partH.items())
        else:
            original={x['element']:x for x in bulk['elements']};checks['exactly138_original_packets_unchanged']=sum(x['element']not in physical_selectors for x in p['effective_bulk_packets'])==138 and all(x==original[x['element']]for x in p['effective_bulk_packets']if x['element']not in physical_selectors)
            checks['exactly6_prescribed_physical_axis_overrides']={x['element']for x in p['effective_bulk_packets']if x['element']in physical_selectors}==set(range(72,78))
    checks['all168_contributions_and6676_Hslots']=len(details)==168 and len(H)==6676
    box=read(ROOT/sp['common_geometry_box_file']);free={int(k):i for i,k in enumerate(arrays['free'])};membership=[]
    for dof in range(650):
        node,c=divmod(dof,2);value=exact_float(arrays['reference'][node,c])+(exact_float(arrays['u'][free[dof]])if dof in free else 0);center=Q(box['current_centers_exact'][node][c]);radius=Q(box['coordinate_radii_exact'][node][c]);membership.append(center-radius<=value<=center+radius)
    theta={k:Q(v)for k,v in row['base_parameters_exact'].items()};checks['all650_original_point_coordinates_in_common_box']=all(membership);checks['base_theta_in_same_closed_C']=all(Q(row['parameter_cell'][k][0])<=v<=Q(row['parameter_cell'][k][1])for k,v in theta.items())
    hrows=A.h_rows(H,n);ARrows=A.absolute_row_product(R,hrows);brow=max(ARrows);RF=A.signed_matvec(R,F)
    ctx.prec=p['signed_product_precision_bits']
    def ball(q):return arb(fmpq(q.numerator,q.denominator))
    RM=arb_mat(n,n,[ball(x)for rr in R for x in rr]);HM=arb_mat(n,n)
    for(i,j),v in H.items():HM[i,j]=HM[j,i]=ball(v[0]).union(ball(v[1]))
    RH=RM*HM;signedrows=[]
    with gzip.open(HERE/'signed_RH_endpoints.ndjson.gz','wt',encoding='utf8')as out:
        for i in range(n):
            values=[(Q(str(RH[i,j].lower().fmpq())),Q(str(RH[i,j].upper().fmpq())))for j in range(n)];signedrows.append(sum(map(A.abs_upper,values),Q(0)));out.write(json.dumps({'i':i,'intervals':[[str(a),str(b)]for a,b in values]},separators=(',',':'))+'\n')
    bsigned=max(signedrows);b=min(brow,bsigned);qbound=read(ROOT/'runs/exact_integration_bridge_v1/cached_Q_common_B_v1/results.json');ZQ=Q(qbound['ZQ_upper_exact']);checks['fresh_ZQ_same_C_radius']=qbound['parameter_cell_sha256']==p['parameter_cell_sha256']and qbound['common_radius_exact']==p['common_radius_exact']and all(v is True for v in qbound['checks'].values());checks['all_frozen_inputs_still_current']=all(sha(ROOT/k)==v for k,v in p['input_hashes'].items())
    exact={'force_endpoint_sums':[[str(a),str(b)]for a,b in F],'H_upper_endpoint_sums':[{'i':i,'j':j,'interval':[str(a),str(b)]}for(i,j),(a,b)in sorted(H.items())],'H_absolute_rows_exact':list(map(str,hrows)),'absR_H_rows_exact':list(map(str,ARrows)),'signed_RH_rows_exact':list(map(str,signedrows)),'signed_RF_endpoint_sums':[[str(a),str(b)]for a,b in RF],'b_absR_H_upper_exact':str(brow),'b_signed_Arb_upper_exact':str(bsigned),'b_chosen_upper_exact':str(b),'ZQ_upper_exact':str(ZQ),'ZQ_plus_b_upper_exact':str(ZQ+b),'force_preconditioned_uniform_upper_exact':str(max(map(A.abs_upper,RF)))}
    write(HERE/'exact_reduction.json',exact)
    result={'status':'PASS_COMBINED_UNIFORM_PHYSICAL_AXIS_V2_COMPONENT_BOUND'if all(v is True for v in checks.values())else'FAIL_COMBINED_UNIFORM_COMPONENT_REVIEW','checks':checks,'utc':now(),'protocol_sha256':sha(HERE/'protocol.json'),'source_sha256':sha(__file__),'input_hashes':p['input_hashes'],'axis_selection':p['packet_selection'],'first_bound_retained_sha256':'5e51264b47bf9ae407021dc429555ca02fceb54efc023dba18f71f022f8b5b76','component_packets':details,'safe_point_consistency':pointoverlaps,'exact_reduction_sha256':sha(HERE/'exact_reduction.json'),'signed_RH_matrix_sha256':sha(HERE/'signed_RH_endpoints.ndjson.gz'),'norm_bounds':{k:v for k,v in exact.items()if k.endswith('_upper_exact')},'parameter_cell_sha256':p['parameter_cell_sha256'],'common_geometry_domain_id':p['common_geometry_domain_id'],'common_radius_exact':p['common_radius_exact'],'ZQ_plus_b_less_than_one':ZQ+b<1,'root_selfmap_or_integral_root_certified':False,'local_reintegration_by_reviewer':False,'reuse_disclosure':'Generic safe local engines and frozen independently reviewed body/cache/domain/prelaunch proofs reused; all168 packets independently reduced and checked, no local integration rerun in this stage. All138 unchanged original bulk packets plus6 prescribed physical-axis packets and24 unchanged surface; first bound preserved. Same-target representation proof and separately reviewed adapter reused.','seconds':time.monotonic()-start,'scope':p['scope']}
    write(HERE/'results.json',result);write(ROOT/'evidence/exact_uniform_full_B_physical_axis_independent_review_v2.json',result);write(HERE/'terminal.json',{'status':'FINISHED'if all(v is True for v in checks.values())else'FAILED_REVIEW','exit_code':0 if all(v is True for v in checks.values())else 1,'pid':os.getpid(),'utc':now(),'results_sha256':sha(HERE/'results.json')});print(json.dumps({'status':result['status'],'checks':len(checks),'failed':[k for k,v in checks.items()if not v],'b_row':float(brow),'b_signed':float(bsigned),'b':float(b),'ZQ_plus_b':float(ZQ+b),'seconds':result['seconds'],'receipt_sha256':sha(ROOT/'evidence/exact_uniform_full_B_physical_axis_independent_review_v2.json')}),flush=True)
    assert all(v is True for v in checks.values()),checks
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('command',choices=['prepare','run']);a.add_argument('--bulk-results-sha');a.add_argument('--surface-results-sha');x=a.parse_args();prepare(x.bulk_results_sha,x.surface_results_sha)if x.command=='prepare'else run()
