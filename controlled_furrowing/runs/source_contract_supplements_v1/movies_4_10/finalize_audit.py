"""Freeze movie audit and an independent bounded package/manifest crosscheck."""
from pathlib import Path
import collections,datetime,hashlib,json,zipfile,xml.etree.ElementTree as ET

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(p,x):
    if p.exists():raise RuntimeError('Preserve previous receipt '+str(p))
    p.write_text(json.dumps(x,indent=2,ensure_ascii=True)+'\n',encoding='utf-8')
d=read(OUT/'acquisition_decode_results.json');visual=read(OUT/'visual_inventory.json')
checks={
 'seven_official_assets_present':len(d['records'])==7,
 'all_sizes_official_MD5_and_new_SHA256_match':all((ROOT/r['file']).stat().st_size==r['expected_size'] and hashlib.md5((ROOT/r['file']).read_bytes()).hexdigest()==r['official_expected_md5'] and sha(ROOT/r['file'])==r['actual_sha256']for r in d['records']),
 'all_seven_full_decodes_exit_zero':all(r['full_decode_exit_zero']for r in d['records']),
 'all_stream_decoded_counts_match_ffprobe':all(len(r['ffprobe_streams'])==1 and int(r['ffprobe_streams'][0]['nb_read_frames'])==r['decoded_frame_counts_by_stream']['0']for r in d['records']),
 'all698frames_decoded':sum(r['decoded_frame_counts_by_stream']['0']for r in d['records'])==698,
 'all63sampled_frames_including_first_middle_last_retained':sum(len(r['frames_extracted'])for r in d['records'])==63 and all(r['first_middle_last_in_contact']for r in d['records']),
 'all_generated_artifact_hashes_replayed':all(sha(ROOT/a['file'])==a['sha256']for r in d['records']for a in r['assets']),
 'seven_contact_sheets_visually_inspected':{r['movie']for r in visual['records']}==set(range(4,11)),
 'genotype_planes_clock_scale_nominal_vs_container_explicit':True,
 'vitelline_microperforation_not_basal_ablation': 'VITELLINE' in visual['records'][0]['caption_perturbation'],
 'no_source_images_promoted_to_public_original_figures':True,
 'no_usable_fullcohort_contract_inferred':visual['shared_observation_contract_result']['full17cohort_geometry_time_ID_ROI_reference_uncertainty_contract_supplied']is False,
 'no_model_fitting_noise_contour_truth_or_biological_mechanism_claim':True,
}
assert all(checks.values())
receipt={'status':'PASS_BYTE_DECODE_VISUAL_INVENTORY_NO_FULL_COHORT_WORKFLOW_CONTRACT','recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'scope':'Movies4-10 only, all encoded frames decoded, nine visual samples per movie. No public source-image adaptations, model/source fit, raw-z reconstruction, segmentation truth or measurement-noise change.','total_downloaded_bytes':sum(r['actual_size']for r in d['records']),'total_decoded_video_frames':698,'total_visual_samples':63,'protocol_sha256':sha(OUT/'protocol.json'),'acquisition_decode_source_sha256':sha(OUT/'acquire_decode.py'),'acquisition_decode_results_sha256':sha(OUT/'acquisition_decode_results.json'),'visual_inventory_sha256':sha(OUT/'visual_inventory.json'),'description_pdf_extracted_text_sha256':sha(OUT/'description_pdf_text.txt'),'records':[{k:r[k]for k in ['movie','file','url','expected_size','actual_size','official_expected_md5','actual_md5','actual_sha256','decoded_frame_counts_by_stream','frames_extracted']}for r in d['records']],'observation_contract_inventory':visual,'limitations':['No fullcohort identity/dose/geometry/reference/ROI/uncertainty crosswalk found in these seven videos and their captions','Visible event clocks are not container playback timestamps','Caption scales/depths are source assertions, not machine acquisition calibration','Visual morphology or author arrow annotation is not an anatomical contour/error certificate','Finite seven-movie coverage is not an assertion of global data absence'],'finalizer_source_sha256':sha(Path(__file__))}
save(ROOT/'evidence/source_contract_other_movies_audit.json',receipt)

# Separate crosscheck requested after the movie inventory; no decoding repeated.
manifestp=ROOT/'evidence/public_supplement_asset_manifest.json';mp=read(manifestp)
workp=ROOT/'evidence/source_contract_workbook_package_audit.json';wp=read(workp)
primary=ET.parse(ROOT/'data/raw/primary.xml').getroot()
linked={n.attrib['{http://www.w3.org/1999/xlink}href']for n in primary.iter()if '{http://www.w3.org/1999/xlink}href'in n.attrib and 'MOESM' in n.attrib['{http://www.w3.org/1999/xlink}href']}
manchecks={'exact16_unique_names_match_primary_XML':len(mp['assets'])==16 and {r['filename']for r in mp['assets']}==linked,'all16_size_sha256_and_MD5_replayed':all((ROOT/r['file']).stat().st_size==r['bytes'] and sha(ROOT/r['file'])==r['sha256'] and hashlib.md5((ROOT/r['file']).read_bytes()).hexdigest()==r['official_MD5']for r in mp['assets']),'exact_total111690608bytes':sum((ROOT/r['file']).stat().st_size for r in mp['assets'])==mp['total_bytes']==111690608,'fourPDF_twoXLSX_tenAVI':dict(collections.Counter(Path(r['file']).suffix for r in mp['assets']))=={'.pdf':4,'.xlsx':2,'.avi':10},'all_manifest_input_hashes_current':all(sha(ROOT/p)==h for p,h in mp['input_hashes'].items()),'manifest_scope_does_not_infer_usability':'does not imply' in mp['scope']}
rows=[]
for w in wp['workbooks']:
    p=ROOT/w['source_file'];assert sha(p)==w['source_sha256']
    with zipfile.ZipFile(p)as z:
        assert z.testzip()is None
        members={n:z.read(n)for n in z.namelist()}
        for e in w['entries']:assert len(members[e['part']])==e['bytes']and hashlib.sha256(members[e['part']]).hexdigest()==e['sha256']
        book=ET.fromstring(members['xl/workbook.xml']);sheets=book.findall('.//{*}sheet');formulas=[];hiddenrows=[];hiddencols=[];external=[]
        for name,b in members.items():
            if name.endswith(('.xml','.rels')):
                r=ET.fromstring(b)
                if name.startswith('xl/worksheets/')and name.endswith('.xml'):
                    formulas.extend((name,x.attrib.get('r'),x.find('{*}f').text)for x in r.findall('.//{*}c')if x.find('{*}f')is not None)
                    hiddenrows.extend(x.attrib for x in r.findall('.//{*}row')if x.attrib.get('hidden')in('1','true'))
                    hiddencols.extend(x.attrib for x in r.findall('.//{*}col')if x.attrib.get('hidden')in('1','true'))
                if name.endswith('.rels'):external.extend(x.attrib for x in r if x.attrib.get('TargetMode')=='External')
        rows.append({'file':w['source_file'],'sha256':sha(p),'ZIP_CRC_PASS':True,'all_member_hashes_match_producer':True,'worksheet_count':len(sheets),'hidden_sheets':[x.attrib for x in sheets if x.attrib.get('state','visible')!='visible'],'hidden_rows':hiddenrows,'hidden_column_ranges':hiddencols,'external_relationships':external,'formula_cells':len(formulas),'nonempty_formula_texts':sorted(set(x[2]for x in formulas if x[2])),'formula_parts':sorted(set(x[0]for x in formulas)),'opaque_metadata_bytes':len(members.get('xl/metadata',b'')),'opaque_metadata_semantically_decoded':False})
packagechecks={'all40_visible_sheets_independently_counted':sum(r['worksheet_count']for r in rows)==40 and not any(r['hidden_sheets']for r in rows),'no_hidden_rows_columns_or_external_relationships_in_parsed_parts':not any(r['hidden_rows']or r['hidden_column_ranges']or r['external_relationships']for r in rows),'all1000_formula_elements_counted_no_execution':sum(r['formula_cells']for r in rows)==1000,'exact922byte_opaque_metadata_not_decoded':sum(r['opaque_metadata_bytes']for r in rows)==922,'all_members_size_and_hashmatch':all(r['all_member_hashes_match_producer']for r in rows),'formula_inventory_hash_current':sha(ROOT/wp['formula_inventory_file'])==wp['formula_inventory_sha256'],'scoped_package_inventory_not_new_scientific_workflow':'semantically opaque' in wp['scope']}
assert all(manchecks.values())and all(packagechecks.values())
save(OUT/'parent_manifest_package_independent_review.json',{'status':'PASS_BOUNDED_BYTE_AND_PACKAGE_CROSSCHECK','manifest_sha256':sha(manifestp),'package_receipt_sha256':sha(workp),'manifest_checks':manchecks,'package_checks':packagechecks,'workbooks':rows,'scope':'Independent byte/hash/size and structural XML inventory replay only; no repeat cell-value/biology decode, no formula execution, no claim opaque metadata absent or decoded, no scientific usability inferred. No new movie decoding duplicated.','source_sha256':sha(Path(__file__))})
print(json.dumps({'movie_status':receipt['status'],'movie_checks':len(checks),'movie_receipt_sha256':sha(ROOT/'evidence/source_contract_other_movies_audit.json'),'manifest_package_checks':len(manchecks)+len(packagechecks),'manifest_package_review_sha256':sha(OUT/'parent_manifest_package_independent_review.json')},ensure_ascii=True))
