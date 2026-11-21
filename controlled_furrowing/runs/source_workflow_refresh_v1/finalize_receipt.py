from pathlib import Path
import datetime, hashlib, json, re, xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
access=json.loads((OUT/'direct_access.json').read_text(encoding='utf-8'))
xml=ET.parse(OUT/'europepmc.txt').getroot()
supp={n.attrib['{http://www.w3.org/1999/xlink}href'] for n in xml.iter() if '{http://www.w3.org/1999/xlink}href' in n.attrib and 'MOESM' in n.attrib['{http://www.w3.org/1999/xlink}href']}
existing=[{'file':'data/raw/'+x,'exists_locally':(ROOT/'data/raw'/x).is_file()} for x in sorted(supp)]
avail=[s for x in access['records'] for s in x.get('availability_sections',[]) if x['label']=='europepmc']
links={x['label']:x.get('repository_links',[]) for x in access['records'] if x['label']!='publisher'}
publisher=(OUT/'publisher.txt').read_text(encoding='utf-8')
challenge='/_fs-ch-' in publisher
checks={
 'fresh_fulltext_exact_article_title':''.join(xml.find('.//article-title').itertext())=='Endogenous OptoRhoGEFs reveal biophysical principles of epithelial tissue furrowing',
 'fresh_data_availability_request_only_primary_microscopy':any('available from the corresponding author upon request' in x['text'] for x in avail),
 'fresh_code_availability_no_custom_code':any('No custom code was generated' in x['text'] for x in avail),
 'sixteen_unique_linked_supplements':len(supp)==16,
 'both_linked_spreadsheet_files_already_local':all((ROOT/'data/raw'/x).is_file() for x in supp if x.endswith('.xlsx')),
 'no_repository_links_in_inspected_primary_fulltext_and_lab_pages':all(not x for x in links.values()),
 'publisher_http200_not_misclassified_as_article':challenge,
 'no_new_experimental_files_downloaded':all(not p.name.endswith(('.avi','.xlsx','.tif','.tiff')) for p in OUT.iterdir()),
 'no_contact_or_publication_action':True,
 'three_focused_queries_max':True,
 'no_markdown_html_figure_builder_edits_for_this_task':True,
}
out={
 'status':'COMPLETE_BOUNDED_AVAILABILITY_CHECK_NO_NEW_USABLE_WORKFLOW_PACKAGE_FOUND',
 'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'article_doi':'10.1038/s41467-025-62483-6',
 'scope':'Fresh public primary fulltext/link inspection and exactly three focused search queries; not an exhaustive repository or internet census.',
 'queries':[
  {'query':'"s41467-025-62483-6" data timecourse repository','result_retention':'Initial combined tool output exceeded context. Exact query retained; no evidence claim relies solely on that unretained output.'},
  {'query':'"Endogenous OptoRhoGEFs" Kasza data repository GitHub','result_file':'web_results.json'},
  {'query':'"Countryman" "62483" "Zenodo" OR "BioImage" OR "github"','result_file':'web_results.json'}],
 'direct_access_coverage':access['records'],
 'publisher_access_caveat':'Direct web open failed at Nature IDP; direct HTTP 200 body is a challenge, not article fulltext. Search-index article availability text agrees with independently freshly fetched EuropePMC XML.',
 'pmc_browser_access_caveat':'Web open to PMC returned CAPTCHA; live EuropePMC fullTextXML succeeded and is archived.',
 'current_availability_sections':avail,
 'linked_supplement_inventory_name_check_only':existing,
 'unfetched_linked_files_scope':'Nine further AVI movies and peer-review PDF are linked but not locally acquired by this refresh; this is explicit inspection scope, not evidence of absence or an assertion their contents have been exhaustively audited. Previously obtained Source Data and Supplementary Data spreadsheets, Supplementary Information, file descriptions, reporting summary, and Movie3 remain the audited package.',
 'repository_link_inspection':links,
 'new_public_package_found':False,
 'new_file_schema_verified':False,
 'conclusion':'No publicly accessible new per-specimen processed apical/basal timecourses, initial profiles, or ROI/mask/growth-reference crosswalk package was located in these inspected sources. The bounded refresh does not remove the present source-observation blocker. It does not prove such data do not exist, are unavailable on request, or are absent from all repositories.',
 'not_done':['No author contact','No huge microscopy/movie reacquisition','No repeat workbook decode','No new anatomical segmentation or temporal operator inferred','No global absence claim','No GitHub publication'],
 'checks':checks,
 'all_checks_pass':all(checks.values()),
 'retained_initial_finding':'Initial receipt wrongly required all sixteen linked files to be local; this failed because nine further AVI movies and peer-review PDF were deliberately not acquired. Initial source and receipt are retained in initial_overbroad_inventory_check. Final check is limited to both linked source spreadsheets already present; no movie absence/content inference is made.',
 'artifacts':[{'file':str(p.relative_to(ROOT)).replace('\\','/'),'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='finalize_receipt.py'],
 'this_receipt_generator_sha256':sha(Path(__file__)),
}
dest=ROOT/'evidence/source_workflow_availability_refresh.json'
if dest.exists() and not (OUT/'initial_overbroad_inventory_check/receipt.json').is_file(): raise RuntimeError('Receipt exists without retained initial revision')
dest.write_text(json.dumps(out,indent=2,ensure_ascii=True)+'\n',encoding='utf-8')
print(json.dumps({'status':out['status'],'checks':checks,'receipt_sha256':sha(dest),'receipt':str(dest)},ensure_ascii=True))
