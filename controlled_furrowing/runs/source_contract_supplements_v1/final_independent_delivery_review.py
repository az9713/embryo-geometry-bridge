"""Read-only byte/schema audit of settled delivery v2; no scientific reruns."""
import hashlib,json,re
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlparse,parse_qs
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).with_suffix('.json')
FINAL='evidence/source_contract_final_delivery_v2.json'
EXPECTED='7e6cd17aa1183dc848f6a47d16615a2243c1d97578ed8426357ff4a84eb021bf'
checks={}; inspected={}
def digest(path,algorithm='sha256'):
    h=hashlib.new(algorithm)
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''): h.update(block)
    return h.hexdigest()
def sha(f):
    p=ROOT/f; h=digest(p); inspected[f]=h; return h
def load(f):
    sha(f); return json.loads((ROOT/f).read_text(encoding='utf-8'))
def ck(name,value):
    checks[name]=bool(value)
    if not value: raise AssertionError(name)
def mapping(d,label):
    ck(label,all(sha(k)==v for k,v in d.items()))
def truechecks(j,label):
    d=j['checks']; ck(label,isinstance(d,dict) and bool(d) and all(v is True for v in d.values()))
assert not OUT.exists(),'Write-once review already exists'
final=load(FINAL); ck('authoritative_v2_bytes',sha(FINAL)==EXPECTED)
truechecks(final,'producer20_checks_literal_true'); ck('producer20_check_count',len(final['checks'])==20)
mapping(final['file_hashes_relative_to_child'],'all24_current_delivery_bindings')
mapping(final['rendered_snapshot_hashes'],'entire_current_rendered_snapshot')
asset=load('evidence/public_supplement_asset_manifest.json'); truechecks(asset,'asset_manifest_checks')
meta=load('data/raw/PMC12361397.1.json')
official={Path(urlparse(u).path).name:parse_qs(urlparse(u).query)['md5'][0] for u in meta['media_urls']}
assets=asset['assets']; ck('16_unique_supplement_ids',len(assets)==16 and {x['source_id'] for x in assets}=={f'MOESM{i}' for i in range(1,17)})
asset_results=[]
for a in assets:
    p=ROOT/a['file']; actualsha=sha(a['file']); actualmd5=digest(p,'md5')
    ck('asset_'+a['source_id'],p.stat().st_size==a['bytes'] and actualsha==a['sha256'] and actualmd5==a['official_MD5']==official[a['filename']] and p.name==a['filename'])
    asset_results.append({'file':a['file'],'bytes':p.stat().st_size,'sha256':actualsha,'md5':actualmd5})
ck('package_total_bytes_formats',sum(x['bytes'] for x in assets)==111690608 and Counter(Path(x['file']).suffix for x in assets)=={'.pdf':4,'.xlsx':2,'.avi':10})
audits={}
for key in ['workbook_package','small_movies','other_movies','peer_review']:
    audits[key]=load('evidence/source_contract_'+key+'_audit.json');truechecks(audits[key],key+'_checks_literal_true')
ck('current_peer_review_version',sha('evidence/source_contract_peer_review_audit.json')=='bb4303755c92831cf64c955d68c1345473501aedbc46b2c212080a250cab938f')
small=audits['small_movies']; other=audits['other_movies']; m3=load('evidence/movie_geometry_audit.json')
ck('small_movie_frame_evidence',len(small['movies'])==2 and [x['decoded_frames'] for x in small['movies']]==[159,20] and all(sha(x['source_file'])==x['source_sha256'] for x in small['movies']))
ck('Movie3_all26_frame_hashes',m3['container']['decoded_frames']==26 and [x['index'] for x in m3['frames']]==list(range(26)) and all(re.fullmatch('[0-9a-f]{64}',x['decoded_bgr_sha256']) for x in m3['frames']))
mapping({x['path']:x['sha256'] for x in m3['source_hashes']},'Movie3_source_bindings')
logfile=m3['container']['ffmpeg_log'];sha(logfile);log=(ROOT/logfile).read_text(encoding='utf-8',errors='replace')
ck('Movie3_retained_complete_decode_log',bool(re.search(r'frame=\s*26\b',log)))
other_counts=[sum(x['decoded_frame_counts_by_stream'].values()) for x in other['records']]
ck('other_seven_full_decode_counts',len(other_counts)==7 and sum(other_counts)==other['total_decoded_video_frames']==698)
ck('903_frames_reconstructed',sum(x['decoded_frames'] for x in small['movies'])+len(m3['frames'])+sum(other_counts)==903)
books=audits['workbook_package']['workbooks']
ck('two_workbooks40_sheets1000_formulas',len(books)==2 and sum(x['worksheet_count'] for x in books)==40 and sum(x['formula_cells'] for x in books)==1000)
base='runs/source_contract_supplements_v1/'
binding=load(base+'public_source_audit_binding.json');truechecks(binding,'current_document_binding_checks')
mapping(binding['bindings'],'draft_and_visual_bindings'); mapping(binding['source_audit_bindings'],'source_audit_and_math_bindings')
ck('chapter_is_exact_corrected_draft',sha('docs/PUBLIC-SUPPLEMENT-AUDIT.md')==sha(base+'public-source-audit-draft.md'))
ck('request_is_exact_v2_draft',sha('docs/DRAFT-DATA-REQUEST.md')==sha(base+'data-request-v2.md'))
ck('fullMovie3_bound_not_only_three_snapshots',binding['source_audit_bindings']['evidence/movie_geometry_audit.json']==sha('evidence/movie_geometry_audit.json'))
chapter=(ROOT/'docs/PUBLIC-SUPPLEMENT-AUDIT.md').read_text(encoding='utf-8')
ck('corrected_chapter_direct_full_decode_citation','movie_geometry_audit.json' in chapter)
correction=load(base+'published_citation_correction.json')
ck('citation_correction_versions_explicit',correction['new_published_MD_sha256']==sha('docs/PUBLIC-SUPPLEMENT-AUDIT.md') and correction['binding_sha256']==sha(base+'public_source_audit_binding.json') and correction['correct_full26frame_receipt_sha256']==sha('evidence/movie_geometry_audit.json') and correction['old_published_MD_sha256']!=correction['new_published_MD_sha256'])
pub=load(base+'publication_source_receipt.json')
ck('four_publication_controls_unchanged',len(pub['after_controls_relative_to_repository'])==4 and all(digest(ROOT.parent/k)==v for k,v in pub['after_controls_relative_to_repository'].items()))
ck('historical_initial_receipt_preserved',sha('evidence/source_contract_final_delivery.json')==correction['prior_initial104delivery_sha256'])
snap=load('evidence/readiness_snapshot_verification.json');truechecks(snap,'fresh_snapshot_checks')
ck('snapshot_exactly_current_final_v2',snap['status']=='PASS' and snap['snapshot_sha256']==final['rendered_snapshot_hashes'] and snap['html_pages']==104 and snap['page_viewport_checks']==208)
browser=load('evidence/report_browser_verification.json');rows=browser['checks']
expected={(p.name,w) for p in (ROOT/'docs').glob('*.html') for w in [390,1440]}
ck('104_pages208_exact_viewport_coverage',len(expected)==208 and len(rows)==208 and {(x['page'],x['width']) for x in rows}==expected)
ck('all_existing_browser_rows_clean',browser['status']=='PASS' and all(x['pass'] is True and not x['browser_errors'] and not x['missing_links'] and not x['errors'] and x['math']==x['rendered'] and not x['rawProseTeX'] and x['images'] is True and x['overflow'] is False for x in rows))
ck('fresh_snapshot_browser_receipt_bound',snap['browser_receipt_sha256']==sha('evidence/report_browser_verification.json') and json.loads(snap['browser_stdout'])=={'status':'PASS','page_viewport_checks':208,'failures':[]} and not snap['browser_stderr'])
raw=load('evidence/raw_prose_math_scan.json')
ck('raw_scan104_current_html_clean',raw['status']=='PASS' and raw['html_pages']==104 and len(raw['rows'])==104 and {x['page'] for x in raw['rows']}=={p.name for p in (ROOT/'docs').glob('*.html')} and all(not x['raw_TeX_commands_in_prose'] and sha('docs/'+x['page'])==x['HTML_sha256'] for x in raw['rows']))
space=load('evidence/source_contract_space_update.json');truechecks(space,'retained_Space_readback_checks')
image=base+'public-source-contract.png'; visual=load(base+'root_visual_inspection.json');truechecks(visual,'retained_root_visual_inspection_checks')
ck('retained_Space106_seq15_original_visual',space['after_blocks']==106 and space['through_sequence']==15 and space['visual']['sha256']==sha(image)==visual['image_sha256'] and (ROOT/image).stat().st_size==space['visual']['byte_size'] and space['peer_review_audit_sha256']==sha('evidence/source_contract_peer_review_audit.json'))
goldens={
'evidence/shared_cell_cohort_inference.json':'1f8c4dc1ad8924f8fda13e6748d5e73aac426b5771b02404803d27c1955aa677',
'evidence/shared_cell_cohort_independent_review.json':'2f7da9f2bee4c1b7e4fdd7f3a88ed90dd914fb218269b3d65bf12bb8a6ab9114',
'runs/cohort_energy_stability_v1/results.json':'5f590ed1bea078e2e88d17bb39f8030d144a95fb73bf96ab10b5823951883164',
'evidence/cohort_energy_stability_independent_review.json':'d4f4c4d30d16bbbf2141cabf1cee18b47fcc35ab958a2a8a203000f0a7acaf95',
'runs/cohort_energy_stability_v1/protocol.json':'804633c41fafc57b4254204ae389a71cc8ada09b2bb69c1c2f36cb031473a82e',
'scripts/cohort_energy_stability.py':'a3fe42c9d79907103c0e02af35db3bdc4c2f89aee59755c006e9ab09e322ee7d',
'runs/shared_cell_cohort_v1/protocol.json':'9a513b33eb8ef589a8ea82561e476ffc9ca8440dcef5ee4488ff11710165f7ea'}
mapping(goldens,'seven_mathematical_golden_packets_unchanged')
energy=load('runs/cohort_energy_stability_v1/protocol.json'); mapping(energy['source_hashes'],'14_frozen_energy_sources');mapping(energy['input_hashes'],'181_frozen_energy_inputs')
shared=load('runs/shared_cell_cohort_v1/protocol.json');mapping(shared['file_hashes'],'71_frozen_shared_cell_inputs')
checkpoint=load('runs/cohort_static_observers_v1/checkpoint_documentation.json')
ck('original_PLAN_unchanged',sha('PLAN.md')==checkpoint['bindings']['PLAN.md']=='993240af5d30c8e5b72ab31a1278a524314d38ede48486b677b021d84bcde500')
gates=load('evidence/gates.json');contact=load(base+'public_contact_verification.json')
ck('full_goal_and_novelty_remain_unestablished',gates['goal_achieved'] is False and final['full_original_five_gate_goal_complete'] is False and final['confirmed_novelty'] is False)
ck('no_author_contact_or_invented_sender',final['author_request_sent'] is False and final['external_contact_authorization_question_pending'] is True and contact['message_sent'] is False and contact['sender_identity_supplied'] is False and contact['external_contact_authorized'] is False)
ck('final_v2_and_snapshot_stable_during_readonly_review',digest(ROOT/FINAL)==EXPECTED and all(digest(ROOT/k)==v for k,v in final['rendered_snapshot_hashes'].items()))
result={'status':'PASS_COMPONENT_DELIVERY_REVIEW','recorded_utc':datetime.now(timezone.utc).isoformat(),'checks':checks,'authoritative_delivery_file':FINAL,'authoritative_delivery_sha256':EXPECTED,'source_verifier_sha256':digest(Path(__file__)),'assets_independently_byte_hashed':asset_results,'decoded_frame_evidence_reconstruction':{'movies_1_2':179,'movie_3':26,'movies_4_10':698,'total':903},'inspected_file_sha256':inspected,'resolved_finding':'Initial chapter linked a three-snapshot receipt for complete Movie3 decoding. Corrected current chapter and v2 binding explicitly cite the existing full26-frame audit; old delivery snapshots remain historical. No scientific bytes were changed.','verification_scope':'Independent read-only byte, schema, count, current-rendered-snapshot and retained-receipt audit. Existing browser/movie/Space receipts inspected; no new browser runs, movie decodes, PDF content reads, remote Space calls, simulations or author contact. This is component delivery review, not independent biological validation or re-proving numerical theorems.','full_original_goal_complete':False,'confirmed_novelty':False,'author_contact_sent':False}
OUT.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':result['status'],'checks':len(checks),'files_hashed':len(inspected),'receipt':str(OUT),'receipt_sha256':digest(OUT),'source_sha256':digest(Path(__file__))}))
