"""Read-only final shared-cell component delivery evidence audit.
No equilibrium, force/Hessian, numerical proof, full browser or document edits.
Live Space source is the retained independent native readback from this task.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
ROOT=Path(__file__).resolve().parents[2]
STUDY=ROOT.parent
RUN=Path(__file__).resolve().parent
OUT=ROOT/'evidence/shared_cell_delivery_independent_review.json'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
yes=lambda d:isinstance(d,dict) and bool(d) and all(v is True for v in d.values())


def main():
    if OUT.exists():raise RuntimeError('preserve completed delivery review')
    names=['readiness_snapshot_verification.json','report_browser_verification.json','raw_prose_math_scan.json','rendering_guard_independent_review.json','shared_cell_cohort_independent_review.json','shared_cell_cohort_inference.json','parent_index_checkpoint_verification.json','space_shared_cell_update.json']
    files={f'evidence/{n}':sha(ROOT/'evidence'/n) for n in names};snapshot=read(ROOT/'evidence'/names[0]);browser=read(ROOT/'evidence'/names[1]);raw=read(ROOT/'evidence'/names[2]);render=read(ROOT/'evidence'/names[3]);old180=read(ROOT/'evidence'/names[4]);math=read(ROOT/'evidence'/names[5]);parent=read(ROOT/'evidence'/names[6]);space=read(ROOT/'evidence'/names[7])
    pages={p.name for p in (ROOT/'docs').glob('*.html')};pairs={(c['page'],c['width']) for c in browser['checks']};expected={(p,w) for p in pages for w in (390,1440)}
    checks={'immutable_snapshot_reports_PASS':snapshot['status']=='PASS' and yes(snapshot['checks']),
            'all214_current_snapshot_files_hashmatch':len(snapshot['snapshot_sha256'])==214 and all(sha(ROOT/p)==v for p,v in snapshot['snapshot_sha256'].items()),
            'exact102_pages_204_viewports_no_duplicates_missing':len(pages)==snapshot['html_pages']==102 and len(browser['checks'])==snapshot['page_viewport_checks']==len(pairs)==204 and pairs==expected,
            'browser_receipt_exact_hash_binding':snapshot['browser_receipt_sha256']==sha(ROOT/'evidence/report_browser_verification.json'),
            'all204_browser_rows_pass_literalTrue':browser['status']=='PASS' and all(c['pass'] is True for c in browser['checks']),
            'all204_math_counts_raw_prose_assets_layout_links_clean':all(c['math']==c['rendered'] and not c['rawProseTeX'] and c['images'] is True and c['overflow'] is False and not c['errors'] and not c['browser_errors'] and not c['missing_links'] for c in browser['checks']),
            'raw_scan_exact102_unique_current_pages':raw['status']=='PASS' and raw['html_pages']==len(raw['rows'])==len({c['page'] for c in raw['rows']})==102 and {c['page'] for c in raw['rows']}==pages,
            'raw_scan_currentHTMLhashes_no_commands':all(not c['raw_TeX_commands_in_prose'] and c['HTML_sha256']==sha(ROOT/'docs'/c['page']) for c in raw['rows']),
            'raw_scan_own_verifier_source_current':raw['verifier_source_sha256']==sha(ROOT/'scripts/scan_raw_prose_math.py'),
            'independent_rendering38checks_PASS':render['status']=='PASS' and len(render['checks'])==38 and yes(render['checks']),
            'rendering_source_and_fixture_hashes_current':all(sha(ROOT/p)==render['source_hashes'][k] for k,p in render['source_files'].items()) and all(sha(ROOT/p)==v for p,v in render['fixture_file_hashes'].items()),
            'old180_math_receipt_unchanged':sha(ROOT/'evidence/shared_cell_cohort_independent_review.json')==sha(ROOT/'runs/rendering_guard_revision_0/shared_cell_cohort_independent_review.json')=='2f7da9f2bee4c1b7e4fdd7f3a88ed90dd914fb218269b3d65bf12bb8a6ab9114'}
    changed=[p for p,v in old180['file_hashes'].items() if sha(STUDY/p)!=v];delta=render['old180_control_source_supersession']
    checks['old180_ONLY_declared_builder_delta']=changed==[delta['path']] and delta['old']==old180['file_hashes'][delta['path']] and delta['new']==sha(STUDY/delta['path'])
    checks['current_render_builder_verifier_in_immutable_browser_snapshot']=snapshot['snapshot_sha256']['scripts/build_report.py']==render['source_hashes']['new_builder'] and snapshot['snapshot_sha256']['scripts/verify_report.py']==render['source_hashes']['new_verifier']
    science=read(ROOT/'runs/shared_cell_cohort_v1/protocol.json');contract=read(ROOT/'runs/cohort_static_observers_v1/model_observation_contract.json')
    checks['all71_scientific_inputs_still_current']=len(science['file_hashes'])==71 and all(sha(ROOT/p)==v for p,v in science['file_hashes'].items())
    checks['all18_physical_observer_kernel_sources_unchanged']=len(contract['source_kernel_hashes'])==18 and all(sha(ROOT/'scripts'/p)==v for p,v in contract['source_kernel_hashes'].items())
    checks['component_math_records17_outputs68_originalgoal_notcomplete']=math['records']==17 and math['constraints']==68 and math['full_five_gate_goal_complete'] is False and math['full_identified_set_computed'] is False and math['confirmed_novelty'] is False
    checks['parent_index_current_2width_PASS']=parent['status']=='PASS' and len(parent['checks'])==2 and {c['width'] for c in parent['checks']}=={390,1440} and all(c['pass'] is True and not c['errors'] and not c['mathError'] and c['allImages'] is True and c['overflow'] is False and c['jointLink'] is True and c['currentText'] is True for c in parent['checks'])
    checks['parent_index_current_html_and_builder_hashes']=parent['page_sha256']==sha(STUDY/'docs/index.html') and parent['builder_sha256']==sha(STUDY/'scripts/build_dossier.py')
    correction=ROOT/'runs/documentation_rendering_corrections_v1';edits=read(correction/'edit_log.json');parsed=read(correction/'math_render_checks.json')
    checks['three_doc_rendering_only_delta_current_and_old_archived']=edits['status']=='RENDERING_ONLY' and len(edits['documents'])==3 and all(sha(ROOT/c['file'])==c['new_sha256'] and sha(ROOT/c['archive_file'])==c['old_sha256'] and c['existing_display_math_untouched'] is True and c['wording_math_and_links_preserved_after_declared_typographical_normalization'] is True for c in edits['documents'])
    checks['three_doc_source_math_render_checks_current']=parsed['status']=='PASS_SOURCE_MATH_RENDERING' and yes(parsed['checks']) and len(parsed['documents'])==3 and all(sha(ROOT/c['file'])==c['sha256'] and yes(c['checks']) and not c['errors'] and not c['raw_TeX_commands_in_source_prose'] and not c['table_hazards'] for c in parsed['documents'])
    live=read(RUN/'space_live_readback.json');page=live['page'];blocks=live['blocks'];text='\n'.join(b['markdown'] for b in blocks);latest='\n'.join(b['markdown'] for b in blocks[-10:])
    checks['space_author_receipt_PASS_and_math_hashes_current']=space['status']=='PASS' and yes(space['checks']) and space['component_inference_sha256']==sha(ROOT/'evidence/shared_cell_cohort_inference.json') and space['independent_combination_review_sha256']==sha(ROOT/'evidence/shared_cell_cohort_independent_review.json')
    checks['fresh_independent_Space91blocks_sequence13']=page['page_id']=='page_c50551506ff881919407cc8f9f3477c2' and page['current_sequence']==page['projection_through_sequence']==space['through_sequence']==13 and len(blocks)==space['after_blocks']==len({b['id'] for b in blocks})==91 and live['references_truncated'] is False
    checks['fresh_Space_latest_section_exact_math_hashes_scope']=sha(ROOT/'evidence/shared_cell_cohort_inference.json') in latest and sha(ROOT/'evidence/shared_cell_cohort_independent_review.json') in latest and all(s in latest for s in ('All 17 records and 68','5.1','5.2')) and all(s in latest.lower() for s in ('noise','continuum','novel','stationary'))
    visu=read(ROOT/'runs/shared_cell_cohort_v1/visual_receipt.json');inspected=read(ROOT/'runs/shared_cell_cohort_v1/visual_inspection.json');native={v['reference']:v for v in live['native_visual_reads']}
    checks['visual_generation_followed_by_bound_inspection']=inspected['status']=='VISUALLY_INSPECTED' and yes(inspected['checks']) and inspected['visual_receipt_sha256']==sha(ROOT/'runs/shared_cell_cohort_v1/visual_receipt.json') and visu['inference_sha256']==sha(ROOT/'evidence/shared_cell_cohort_inference.json')
    checks['all_four_local_visual_files_match_receipt']=all(sha(ROOT/'runs/shared_cell_cohort_v1'/p)==v for p,v in visu['output_hashes'].items())
    imagefiles=('shared-cell-68-observations.png','shared-cell-residual-threshold.png')
    checks['fresh_native_both_image_reads_match_upload_localbytes']=len(native)==2 and all(v['reference'] in text and v['reference'] in native and native[v['reference']]['mime_type']=='image/png' and native[v['reference']]['byte_size']==v['byte_size']==(ROOT/'runs/shared_cell_cohort_v1'/name).stat().st_size and native[v['reference']]['sha256']==v['sha256']==sha(ROOT/'runs/shared_cell_cohort_v1'/name)==inspected['images'][name] for name,v in zip(imagefiles,space['visuals']))
    gates=read(ROOT/'evidence/gates.json');checks['full_five_gates_still_incomplete']=gates['goal_achieved'] is False and gates['remaining_original_requirements']['full_physical_identified_set']=='INCOMPLETE'
    historical=[];h=read(ROOT/'evidence/final_component_binding_check.json');bad={}
    for group in h['groups']:
        for row in group.get('rows',[]):
            path=ROOT/row['path']
            if path.is_file() and sha(path)!=row['expected_sha256']:bad[row['path']]={'historical_expected':row['expected_sha256'],'current':sha(path)}
    historical.append({'receipt':'evidence/final_component_binding_check.json','interpretation':'Point-root-phase historical PASS, not a current delivery verifier; excluded from acceptance.','stale_bindings':bad})
    focus=read(ROOT/'evidence/shared_cell_render_focus.json');historical.append({'receipt':'evidence/shared_cell_render_focus.json','interpretation':'Earlier focused page test superseded by full immutable102/204 browser snapshot; excluded from acceptance.','historical_HTML_sha256':focus['HTML_sha256'],'current_HTML_sha256':sha(ROOT/'docs/SHARED-CELL-COHORT.html')})
    for extra in ['runs/shared_cell_cohort_v1/protocol.json','runs/shared_cell_cohort_v1/visual_receipt.json','runs/shared_cell_cohort_v1/visual_inspection.json','runs/documentation_rendering_corrections_v1/edit_log.json','runs/documentation_rendering_corrections_v1/math_render_checks.json','runs/shared_cell_delivery_review_v1/space_live_readback.json','evidence/gates.json','scripts/scan_raw_prose_math.py']:
        files[extra]=sha(ROOT/extra)
    files['parent:docs/index.html']=sha(STUDY/'docs/index.html');files['parent:scripts/build_dossier.py']=sha(STUDY/'scripts/build_dossier.py')
    result={'status':'PASS_INDEPENDENT_COMPONENT_DELIVERY' if yes(checks) else 'FAIL_CURRENT_DELIVERY_BINDINGS','checks':checks,'failures':[k for k,v in checks.items() if v is not True],'created_utc':datetime.now(timezone.utc).isoformat(),'review_source_sha256':sha(__file__),'file_hashes_relative_to_child_except_parent_prefix':files,'historical_stale_receipts_NOT_used_as_current':historical,'browser_scope':{'current_child_pages':102,'current_child_page_viewports':204,'parent_landing_viewports':2,'fresh_browser_rerun_by_this_reviewer':False,'method':'Current exacthashbinding+inspection of existing immutable browser rows; own independent renderingfixtures already executed separately.'},'Space_scope':'Fresh independent native Page read seq13/91blocks and actual PNGreference reads with exactbytehashes matching local/authorreceipts; no Space edits or contact. Account/ACL identities omitted from retained readback.','root_final_delivery_receipt_intentionally_not_bound':'Root finalizer writes shared_cell_final_delivery.json after consuming thisreview; bindingpendingfilewouldcausecircularstaleness.','component_complete':yes(checks),'original_five_gate_goal_complete':False,'full_physical_identified_set_computed':False,'confirmed_novelty':False,'scope':'Componentdeliveryverification only; no numerical/root/Hessian jobs, productiondocument/HTML edits, publicdeployment or newscientificproof. Backslash-command scan remainsheuristic, notuniversalsemanticmathdetection.'}
    OUT.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8');print(json.dumps({'status':result['status'],'checks':len(checks),'failures':result['failures'],'receipt_sha256':sha(OUT),'historical_stale_binding_count':len(bad)}))
    if not yes(checks):raise SystemExit(1)
if __name__=='__main__':main()
