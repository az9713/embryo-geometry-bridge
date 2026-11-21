"""Audit the completed component and current publication, not all five gates."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];PROJECT=ROOT.parent
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
checks={}
def ck(k,v):
 checks[k]=bool(v)
 if not v:raise AssertionError(k)
base=ROOT/'runs/exact_integration_bridge_v1';cert=read(base/'integral_root_B_safe_v2/results.json');protocol=read(base/'integral_root_B_safe_v2/protocol.json')
for path,h in protocol['source_and_input_hashes'].items():ck('current_proof_input_'+path,sha(ROOT/path)==h)
ck('actual_uniform_root_and_minimum',cert['root_certified']is True and cert['strict_local_minimum_certified']is True)
ck('new_location_not_from_point_escape_alone',cert['no_integral_root_in_old_cube_certified']is True and Q(cert['integral_root_distance_from_old_predictor_lower_exact'])>Q(cert['r_exact']))
tests=cert['physical_axis_test'];z=Q(tests['z']);eta=Q(tests['eta']);r=Q(cert['r_exact'])
ck('exact_selfmap',z<1 and eta+z*r==Q(tests['selfmap_lhs'])<r)
ck('exact_error',eta/(1-z)==Q(cert['root_error_from_translated_predictor_upper_exact'])<Q('2.557e-12'))
ck('exact_curvature',(1-z)/Q(cert['same_R_infinity_norm_exact'])==Q(cert['uniform_integral_energy_curvature_lower_exact'])>Q('0.0005305'))
ck('first_insufficient_bound_retained',cert['first_representation_test']['strict_sufficient_conditions']is False)
mathreview=ROOT/'evidence/exact_integral_root_B_safe_independent_review.json'
ck('independent_review_literal',sha(mathreview)=='4500e5b34039eafe2e8af742cad739b1133ab414521fd6ed604cf4787d80d561')
ck('independent_review_all_checks',read(mathreview)['checks']and all(v is True for v in read(mathreview)['checks'].values()))
direct=ROOT/'evidence/exact_uniform_direct_root_B_crosscheck.json'
ck('second_residual_route_literal',sha(direct)=='f28b458b66d718c0fa7b4ef35a29ee463fc59612d887dee52f7e661d0d41ecfe')
browser=read(ROOT/'evidence/readiness_snapshot_verification.json')
ck('107_pages_214_actual_viewports',browser['status']=='PASS'and browser['html_pages']==107 and browser['page_viewport_checks']==214 and all(browser['checks'].values()))
ck('browser_document_snapshot_current',all(sha(ROOT/name)==h for name,h in browser['snapshot_sha256'].items()))
ck('browser_receipt_current',browser['browser_receipt_sha256']==sha(ROOT/'evidence/report_browser_verification.json'))
parent=read(ROOT/'evidence/parent_index_checkpoint_verification.json')
ck('parent_two_viewports_current',parent['status']=='PASS'and len(parent['checks'])==2 and all(v['pass']for v in parent['checks'])and parent['page_sha256']==sha(PROJECT/'docs/index.html')and parent['builder_sha256']==sha(PROJECT/'scripts/build_dossier.py'))
raw=read(ROOT/'evidence/raw_prose_math_scan.json')
ck('raw_math_scan_current',raw['status']=='PASS'and raw['html_pages']==107 and all(sha(ROOT/'docs'/v['page'])==v['HTML_sha256']for v in raw['rows']))
sci=ROOT/'evidence/safe_uniform_publication_scientific_readback.json'
ck('targeted_scientific_DOM_literal',sha(sci)=='865d32c8854a36c2bcb1dec275a80931af87ecf5ef09e07dc62e52abe216f711')
ck('targeted_readback_all_checks',read(sci)['status'].startswith('PASS')and read(sci)['checks']and all(v is True for v in read(sci)['checks'].values()))
space=ROOT/'evidence/exact_integral_root_B_space_update.json';cloud=read(space)
ck('space_exact_readback',cloud['status']=='PASS_SPACE_EXACT_READBACK'and all(cloud['checks'].values())and cloud['sequence']==19 and cloud['current_block_count']==130)
ck('staged_root_chapter_published_exactly',sha(ROOT/'docs/EXACT-INTEGRAL-ROOT.md')==sha(HERE/'preparation/EXACT-INTEGRAL-ROOT.md')=='f18745b87096baa52531111ad1dd0576c3d0c6c1015e28195d0dc597a1a3ed56')
ck('margins_figure_bound',sha(HERE/'preparation/root_margins.png')=='e3f5ec5ec719a2bb2f55b77149250697992b462322fee23c258f203ef90bd369')
ck('root_location_figure_bound',sha(HERE/'preparation/root_location.png')=='c480670511205b6c8f36a3f47bc5022e2089afdf0be34e2c93c4d005f14d9c42')
for name in ['manifest.json','additional_receipts_manifest.json']:
 archive=HERE/'prepublication_snapshot';manifest=read(archive/name)
 ck('historical_bytes_preserved_'+name,all(sha(archive/path)==h for path,h in manifest['files'].items()))
ck('original_objective_still_incomplete',read(ROOT/'evidence/gates.json')['goal_achieved']is False and cert['all17_integral_or_PDE_or_biological_validation']is False)
files=[base/'integral_root_B_safe_v2/results.json',base/'integral_root_B_safe_v2/protocol.json',base/'integral_root_B_safe_v2/prelaunch_independent_review.json',base/'integral_root_B_safe_v2/terminal.json',mathreview,direct,sci,space,ROOT/'docs/EXACT-INTEGRAL-ROOT.md',ROOT/'docs/EXACT-INTEGRAL-ROOT.html',ROOT/'docs/SAFE-INTEGRATION-REPLAY.md',ROOT/'docs/SAFE-INTEGRATION-REPLAY.html',ROOT/'docs/index.html',PROJECT/'docs/index.html',PROJECT/'README.md',ROOT/'evidence/gates.json',ROOT/'evidence/readiness_snapshot_verification.json',ROOT/'evidence/parent_index_checkpoint_verification.json',ROOT/'evidence/raw_prose_math_scan.json',HERE/'preparation/root_margins.png',HERE/'preparation/root_location.png',HERE/'preparation/authoring_binding.json',HERE/'local_publication.json',Path(__file__)]
out={'status':'PASS_SAFE_UNIFORM_INTEGRAL_ROOT_COMPONENT_AND_PUBLICATION','checks':checks,
     'artifact_hashes':{p.relative_to(PROJECT).as_posix():sha(p)for p in files},
     'executed_verification_scope':{'mathematics':'Source-bound proofs;1500 independent final scalar/identity checks and33 rejection cases. Full integration kernels and reductions have separately retained reviews.','browser':'107 child pages,214 desktop/mobile viewports; separate parent2 views and targeted new2page4view scientific DOM readback.','space':'130blocks,sequence19,exact status/tail/native-image hashes and preservation of untouched old content.'},
     'full_five_gate_goal_complete':False,'confirmed_novelty':False,'github_pushed':False,'github_pages_deployed':False,
     'scope':'One B-dose exact-spatial-integral finite-element equilibrium for all parameters in the same small C, strict local energy minimum, old-cube exclusion and same-parameter branch gap. No other16 integral roots, PDE/nonaxisymmetric/dynamic/global/biological/physical identified-set claim.'}
dest=ROOT/'evidence/exact_integral_root_B_final_delivery.json';assert not dest.exists()
dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps({'status':out['status'],'checks':len(checks),'artifact_bindings':len(files),'receipt_sha256':sha(dest)}),flush=True)
