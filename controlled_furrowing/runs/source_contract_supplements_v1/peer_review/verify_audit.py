"""Read-only source verification; writes only this audit's replay/evidence files."""
from pathlib import Path
import hashlib, json, re, xml.etree.ElementTree as ET
from datetime import datetime, timezone
import fitz
from pypdf import PdfReader

HERE = Path(__file__).resolve().parent
CF = HERE.parents[2]
RAW = CF / 'data/raw'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p, x): p.write_text(json.dumps(x, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
def norm(s): return re.sub(r'\s+', ' ', s).casefold()

def main():
    original = json.loads((HERE/'download_and_extraction.json').read_text(encoding='utf-8'))
    pdf = RAW/'41467_2025_62483_MOESM2_ESM.pdf'
    b = pdf.read_bytes()
    doc = fitz.open(pdf)
    other = PdfReader(pdf)
    checks = {}
    def check(k,v): checks[k] = bool(v)
    check('official_HTTP_200', original['http']['status']==200)
    check('size_matches_official_listing',len(b)==original['expected_size_from_official_listing']==9250570)
    check('sha256_matches_download',sha(pdf)==original['sha256'])
    check('MD5_matches_official_metadata',hashlib.md5(b).hexdigest()==original['expected_md5_from_official_json'])
    check('metadata_source_unchanged',sha(RAW/'PMC12361397.1.json')==original['primary_metadata_sha256'])
    check('object_listing_unchanged',sha(RAW/'oa_object_listing.xml')==original['listing_sha256'])
    check('both_parsers_36_pages',len(doc)==len(other.pages)==36)
    # Initial extraction used spatially sorted text; default PDF object order is different.
    # The failed default-order equality replay is retained in verification_revision_0.
    pages = [p.get_text(sort=True) for p in doc]
    check('every_page_nonempty',all(t.strip() for t in pages))
    check('all_36_saved_page_hashes_match',all(sha(HERE/x['text_file'])==x['text_sha256'] for x in original['page_inventory']))
    check('all_36_reextracted_page_texts_match',all((HERE/f'page_{i+1:03}.txt').read_text(encoding='utf-8')==t for i,t in enumerate(pages)))
    check('full_primary_parser_extract_hash',sha(HERE/'peer_review.full.txt')==original['fitz_fulltext_sha256'])
    check('full_secondary_parser_extract_hash',sha(HERE/'peer_review.pypdf.txt')==original['pypdf_fulltext_sha256'])
    anchors=[(17,'32-35'),(17,'3-minute developmental window'),(31,'finish scanning our lateral activation volume'),(35,'z-sampling'),(36,'acquisition protocol')]
    anchor_results=[]
    for pg, text in anchors:
        a = norm(text) in norm(pages[pg-1])
        c = norm(text) in norm(other.pages[pg-1].extract_text())
        anchor_results.append({'PDF_page_1based':pg,'anchor':text,'both_parsers':a and c})
    check('key_anchors_in_both_independent_parsers',all(x['both_parsers'] for x in anchor_results))
    render = json.loads((HERE/'local_page_renders.json').read_text(encoding='utf-8'))
    check('five_inspected_peer_source_renders_unchanged',all(sha(HERE/x['file'])==x['sha256'] for x in render['pages']))
    ann=[{'PDF_page_1based':i+1,'uri':l['uri']} for i,p in enumerate(doc) for l in p.get_links() if l.get('uri')]
    check('nine_URI_annotations_all_Zotero',len(ann)==9 and all(x['uri'].startswith('https://www.zotero.org/google-docs/') for x in ann))
    check('no_embedded_files',doc.embfile_count()==0)
    xml=ET.parse(RAW/'primary.xml').getroot()
    ptexts={p.attrib.get('id'):''.join(p.itertext()) for p in xml.iter('p') if p.attrib.get('id')}
    check('distinct_sample_unit_explicit_in_final_primary','distinct samples' in ptexts['Par67'])
    check('PIV_activation_ROI_tracking_explicit','PIV' in ptexts['Par59'])
    check('manual_z_drift_correction_explicit','drift' in ptexts['Par66'])
    fig2=' '.join(xml.find('.//fig[@id="Fig2"]').itertext())
    fig6=' '.join(xml.find('.//fig[@id="Fig6"]').itertext())
    check('Fig2_40_second_activation_and_preimage_sequence','40-s activation' in fig2 and 'unperturbed tissue' in fig2)
    check('explicit_initiation_sentence_located_in_Fig6_not_transplanted_to_Fig2','initiation of photoactivation' in fig6 and 'initiation of photoactivation' not in fig2)
    terms=['maximum','average','apical','basal','reference','cellularization','activation','sampling','resolution','segmentation','reconstruction','drift','PIV','repository','raw data','data availability','embryo','32-35','40']
    hits={term:[i+1 for i,t in enumerate(pages) if term.casefold() in t.casefold()] for term in terms}
    urls=[]
    for i,t in enumerate(pages):
        urls.extend({'PDF_page_1based':i+1,'url_text':m.group(0)} for m in re.finditer(r'https?://[^\s<>]+',t))
    save(HERE/'search_ledger.json',{'coverage':'All 36 peer-review pages read in full; keyword search supplements, not substitutes for reading. Relevant final Methods, Supplementary Information and Reporting Summary cross-checked.','keywords_to_pages':hits,'plaintext_URL_candidates':urls,'URI_annotations':ann,'embedded_file_count':doc.embfile_count(),'anchors':anchor_results})
    bound=[pdf,RAW/'PMC12361397.1.json',RAW/'oa_object_listing.xml',RAW/'primary.xml',RAW/'41467_2025_62483_MOESM1_ESM.pdf',RAW/'41467_2025_62483_MOESM15_ESM.pdf',HERE/'download_and_extraction.json',HERE/'link_annotations.json',HERE/'local_page_renders.json',HERE/'primary_crosscheck_locators.json',HERE/'primary_selected_sections.json',HERE/'search_ledger.json',HERE/'PEER-REVIEW-SOURCE-CONTRACT-AUDIT.md',Path(__file__).resolve()]
    bound += [HERE/'verification_revision_0/verify_audit.py',HERE/'verification_revision_0/verification.json']
    artifact_hashes={str(p.relative_to(CF)).replace('\\','/'):sha(p) for p in bound}
    result={'status':'PASS_BOUNDED_SOURCE_AUDIT' if all(checks.values()) else 'FAIL_SOURCE_VERIFICATION','created_utc':datetime.now(timezone.utc).isoformat(),'checks':checks,'artifact_sha256':artifact_hashes,'source_pdf_sha256':sha(pdf),'source_pdf_bytes':len(b),'source_pages_read':list(range(1,37)),'peer_figure_pages_visually_inspected':[11,12,21,24,27],'supplement_figure_pages_visually_inspected':[10,12],'reporting_summary_pages_visually_inspected':[2],'finding':'PARTIALLY_CONSTRAINED_NOT_RESOLVED','new_source_constraints':[{'locator':'Peer-review PDF p17, first author rebuttal, Reviewer1 minor11 response; final Sec14 Par56; ReportingSummary p2 Data exclusions','finding':'Measured 32–35 micrometre selection, approximately three-minute developmental window; not a certified true specimen geometry interval.'},{'locator':'Peer-review PDF p31, second author rebuttal, Reviewer1 point2 response; final Fig2 b/b-prime/b-double-prime','finding':'Shortening already occurs during the approximately40-second lateral scan; scan-end cannot generally be assumed an undeformed reference; this does not assert response in every specimen. Final Fig2 establishes acquisition sequence and activation duration; explicit time-zero initiation sentence appears in final Fig6 caption.'},{'locator':'Peer-review PDF pp35–36, second author rebuttal, Reviewer2 point3d response; final Sec14 Par56','finding':'Authors compare roughly2-micrometre basal displacement with2-micrometre axial sampling. Sampling pitch is not a segmentation uncertainty bound or a5.2-micrometre discrepancy justification.'}],'contract_gaps_remaining':['Exact regional-average implementation: literal final Methods favors mean of local temporal maxima, but no executable aggregation or individual timecourses are supplied in this source.','Component spatial/temporal aggregation and separate growth-reference allocation not established by this supplement.','Per-acquisition clock/slice mapping, registered initial surfaces and segmentation/ROI weights not reconstructed here.','No new individual embryo crosswalk, per-embryo3D initial geometry, measured noise law or continuum material-parameter ranges recovered.'],'availability_scope':'No new public raw-microscopy repository/accession found in these36pages or9URI annotations; final paper says raw microscopy available by author request. This is not proof data do not exist or cannot be obtained. No contact performed.','source_timing_caption_issue':{'peer_review_PDF_page':27,'panel_time_labels_minutes':[0,2,6,10],'peer_caption_minutes':[0,2,4,6],'final_supplement_PDF_page':12,'final_panel_labels_minutes':[0,2,6,10],'scope':'Different cephalic-furrow condition; do not use inconsistent draft caption as calibration timing.'},'quote_policy':{'report_peer_review_verbatim_words':0,'report_is_paraphrase':True},'limitations':['Read-only source-contract audit, not new biology or novelty.','No model changes, fits, noise estimates, workbook-value redecoding, author contact or publication changes.','Does not change static-proxy conditional proof contracts or complete the original five gates.','All original-source images retained for local inspection only.']}
    result['audit_corrections_retained'] = [
        {'issue':'Initial verifier compared default PDF object-order extraction with the original sorted extraction.','resolution':'Archived initial failed verifier/receipt; get_text(sort=True) reproduces every saved page exactly. No source text or original PDF changed.'},
        {'issue':'A proposed source attribution placed an explicit t=0 initiation sentence in Fig2.','resolution':'Ancestor-specific XML inspection locates the sentence in Fig6. Parent acknowledged and corrected the unpublished request draft. Fig2 duration and acquisition sequence remain source-supported.'}
    ]
    save(HERE/'verification.json',result)
    save(CF/'evidence/source_contract_peer_review_audit.json',result)
    print(json.dumps({'status':result['status'],'checks':len(checks),'failed':[k for k,v in checks.items() if not v],'receipt_sha256':sha(CF/'evidence/source_contract_peer_review_audit.json')},indent=2))
    if not all(checks.values()): raise SystemExit(1)

if __name__=='__main__': main()
