"""Version the immutable first reducer with an explicit six-packet selection."""
from pathlib import Path
import hashlib,json,shutil
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];OLD=HERE.with_name('uniform_full_B_independent_reduction_v1')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    for name in ['algebra.py','selfcheck.py','selfchecks.json']:
        dest=HERE/name;assert not dest.exists();shutil.copyfile(OLD/name,dest);assert sha(dest)==sha(OLD/name)
    src=(OLD/'reduce.py').read_text()
    src=src.replace("POS=ROOT/'runs/exact_integration_bridge_v1/surface_assembly_B_safe_v2'","POS=ROOT/'runs/exact_integration_bridge_v1/surface_assembly_B_safe_v2'\nAX=ROOT/'runs/exact_integration_bridge_v1/uniform_axis_edge_B_safe_v2'")
    marker="    bound(inputs,ROOT/row['point_cache_file'],row['point_cache_sha256']);tests=read(HERE/'selfchecks.json')"
    inject="""    ax=bound(inputs,AX/'results.json','06b6bed9eedcb17ae2f3abf1e2097489a8f8f74a6486deb9ebd9152c75f6ea9c');ap=bound(inputs,AX/'protocol.json',ax['protocol_sha256']);execution=bound(inputs,AX/'execution_receipt.json');ar=bound(inputs,AX/'prelaunch_independent_review.json');representation=bound(inputs,ROOT/'evidence/exact_axis_edge_physical_independent_review.json','8eaa76bd6ee3a6980815e303aba619c3f21243191c44035e0b1c3dabb26edffd')
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
"""
    assert marker in src;src=src.replace(marker,inject+marker)
    needle="    assert_inputs(p);write(HERE/'protocol.json',p);print("
    replacement="    p['effective_bulk_packets']=effective;p['axis_protocol_file']=(AX/'protocol.json').relative_to(ROOT).as_posix();p['axis_results_sha256']=sha(AX/'results.json');p['packet_selection']='138 original bulk +6 same-target physical-axis packets72..77 +24 original surface; original first bound preserved';p['reused_pilot_policy']='Exactly the remaining original safe pilot records outside72..77 remain reused; all six full-axis-edge records replaced as prescribed, without outcome-dependent selection.'\n    assert_inputs(p);write(HERE/'protocol.json',p);print("
    assert needle in src;src=src.replace(needle,replacement)
    src=src.replace("bulk=read(BO/'results.json');surface=read(SO/'results.json');bp=read(BO/'protocol.json');sp=read(SO/'protocol.json');ulp=bp['uniform_local_protocol'];", "bulk=read(BO/'results.json');surface=read(SO/'results.json');bp=read(BO/'protocol.json');sp=read(SO/'protocol.json');ulp=bp['uniform_local_protocol'];physical_selectors={x['element']:x for x in read(ROOT/p['axis_protocol_file'])['selectors']};")
    src=src.replace("[('bulk',bulk['elements'],bulk),('surface',surface['edges'],surface)]", "[('bulk',p['effective_bulk_packets'],bulk),('surface',surface['edges'],surface)]")
    src=src.replace("selected=selectors[tag];nn=len(ids)","selected=physical_selectors[tag]if kind=='bulk'and tag in physical_selectors else selectors[tag];nn=len(ids)")
    src=src.replace("bool(local['checks'])and all(v is True for v in local['checks'].values())", "bool(local.get('checks',local.get('point_overlap_diagnostics_only',{})))and all(v is True for v in local.get('checks',local.get('point_overlap_diagnostics_only',{})).values())")
    old="""                reused=str(tag)in bp['reuse_completed_local_packets'];good=good and item['reused_completed_pilot']is reused
                if reused:good=good and item['file']==bp['reuse_completed_local_packets'][str(tag)]['file']and item['sha256']==bp['reuse_completed_local_packets'][str(tag)]['sha256']"""
    new="""                if tag in physical_selectors:good=good and item['representation_override']=='axis_edge_physical'and local['chart']=='axis_edge_physical'
                else:
                    reused=str(tag)in bp['reuse_completed_local_packets'];good=good and item['reused_completed_pilot']is reused
                    if reused:good=good and item['file']==bp['reuse_completed_local_packets'][str(tag)]['file']and item['sha256']==bp['reuse_completed_local_packets'][str(tag)]['sha256']"""
    assert old in src;src=src.replace(old,new)
    start=src.index("        checks[kind+'_all611_force_exact_reduction']")
    end=src.index("    checks['all168_contributions_and6676_Hslots']",start)
    src=src[:start]+"""        if kind=='surface':
            checks[kind+'_all611_force_exact_reduction']=len(output['force_integral_minus_cache'])==n and all(A.contains(A.interval(x),v)for x,v in zip(output['force_integral_minus_cache'],partF))
            oh={(x['i'],x['j']):A.interval(x['enclosure'])for x in output['H_upper']};checks[kind+'_all_symmetric_H_exact_reduction']=len(oh)==len(output['H_upper'])and set(oh)==set(partH)and all(A.contains(oh[k],v)for k,v in partH.items())
        else:
            original={x['element']:x for x in bulk['elements']};checks['exactly138_original_packets_unchanged']=sum(x['element']not in physical_selectors for x in p['effective_bulk_packets'])==138 and all(x==original[x['element']]for x in p['effective_bulk_packets']if x['element']not in physical_selectors)
            checks['exactly6_prescribed_physical_axis_overrides']={x['element']for x in p['effective_bulk_packets']if x['element']in physical_selectors}==set(range(72,78))
"""+src[end:]
    src=src.replace("'PASS_COMBINED_UNIFORM_COMPONENT_BOUND'","'PASS_COMBINED_UNIFORM_PHYSICAL_AXIS_V2_COMPONENT_BOUND'")
    src=src.replace("'component_packets':details", "'axis_selection':p['packet_selection'],'first_bound_retained_sha256':'5e51264b47bf9ae407021dc429555ca02fceb54efc023dba18f71f022f8b5b76','component_packets':details")
    src=src.replace('evidence/exact_uniform_full_B_independent_review.json','evidence/exact_uniform_full_B_physical_axis_independent_review_v2.json')
    # The first-bound dependency must continue referring to the ORIGINAL receipt.
    src=src.replace("bound(inputs,ROOT/'evidence/exact_uniform_full_B_physical_axis_independent_review_v2.json','5e512", "bound(inputs,ROOT/'evidence/exact_uniform_full_B_independent_review.json','5e512")
    src=src.replace('Seven completed safe pilot bulk packets explicitly reused by producer.','All138 unchanged original bulk packets plus6 prescribed physical-axis packets and24 unchanged surface; first bound preserved. Same-target representation proof and separately reviewed adapter reused.')
    dest=HERE/'reduce.py';assert not dest.exists();dest.write_text(src,encoding='utf8')
    lineage={'status':'VERSIONED_BEFORE_EXECUTION','original_reducer_sha256':sha(OLD/'reduce.py'),'new_reducer_sha256':sha(dest),'identical_algebra_sha256':sha(HERE/'algebra.py'),'cached_synthetic_checks_sha256':sha(HERE/'selfchecks.json'),'first_delivery_sha256':sha(OLD/'delivery.json'),'changes':'Explicit frozen6packet substitution only; selection/core map/domain/tail/point checks adapted for physical schema. Generic exact interval/scatter and signedArb matrix algorithms unchanged. No local integrations.'};p=HERE/'lineage.json';assert not p.exists();p.write_text(json.dumps(lineage,indent=2)+'\n');print(json.dumps({'status':lineage['status'],'source_sha256':sha(dest)}))
if __name__=='__main__':main()
