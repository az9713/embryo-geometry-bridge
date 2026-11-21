"""Read-only current scientific publication and actual Edge DOM readback."""
from pathlib import Path
from fractions import Fraction as Q
from functools import partial
from http.server import ThreadingHTTPServer
import json,hashlib,sys,threading,datetime
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[3];PROJECT=ROOT.parent;OUT=Path(__file__).parent
sys.path.insert(0,str(ROOT/'scripts'))
import build_report as B
import verify_report as V
PUB=ROOT/'runs/exact_integration_bridge_v1/safe_uniform_publication_v1'
RECEIPT=ROOT/'evidence/safe_uniform_publication_scientific_readback.json'
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
checks={}
def ck(k,v):
 checks[k]=bool(v)
 if not v:raise AssertionError(k)
def write(p,d):
 if p.exists():raise ValueError('write-once output exists '+str(p))
 p.write_text(json.dumps(d,indent=2)+'\n')
if RECEIPT.exists():raise ValueError('preserve existing final readback')
pub=read(PUB/'local_publication.json')
ck('publication_current',sha(PUB/'local_publication.json')=='0ce76385458e7ca5036004bf80a6bbbf610a63c8110e95b6b82acef51fa2ae4a')
ck('local_not_deployment_or_fullgoal',pub['github_pushed']is False and pub['public_site_deployed']is False and pub['full_five_gate_goal_complete']is False)
for name,h in pub['source_and_output_hashes'].items():ck('publication_binding_'+name,sha(PROJECT/name)==h)
rootmd=ROOT/'docs/EXACT-INTEGRAL-ROOT.md';pointmd=ROOT/'docs/SAFE-INTEGRATION-REPLAY.md';staged=PUB/'preparation/EXACT-INTEGRAL-ROOT.md'
ck('staged_final_MD_exact',rootmd.read_bytes()==staged.read_bytes()and sha(staged)=='f18745b87096baa52531111ad1dd0576c3d0c6c1015e28195d0dc597a1a3ed56')
for path in [rootmd,pointmd]:ck('HTML_matches_exact_MD_'+path.stem,(ROOT/f'docs/{path.stem}.html').read_text()==B.shell(path.stem,B.convert(path.read_text())))
cert=read(ROOT/'runs/exact_integration_bridge_v1/integral_root_B_safe_v2/results.json');review=read(ROOT/'evidence/exact_integral_root_B_safe_independent_review.json');direct=read(ROOT/'evidence/exact_uniform_direct_root_B_crosscheck.json')
ck('root_review1500_actual',review['status']=='PASS'and len(review['checks'])==1500 and all(v is True for v in review['checks'].values())and review['certificate_sha256']==sha(ROOT/'runs/exact_integration_bridge_v1/integral_root_B_safe_v2/results.json'))
ck('direct_route_shared_kernels_disclosed',direct['status']=='PASS_DIRECT_UNIFORM_FORCE_ROOT_CROSSCHECK'and direct['strict_selfmap_pass']is True and direct['whole_physical_kernel_independently_reimplemented']is False and all(v is True for v in direct['checks'].values()))
t=rootmd.read_text();pt=pointmd.read_text();old=(ROOT/'docs/EXACT-INTEGRATION-BRIDGE.md').read_text();rm=(PROJECT/'README.md').read_text();ci=(ROOT/'docs/index.html').read_text();pi=(PROJECT/'docs/index.html').read_text()
ck('root_scope_oneB_C_611_PDElimits',all(s in t for s in ['one source-dose record B','same small closed cell C','611-coordinate','does not yet concern all 17 records','continuum elasticity PDE']))
ck('new_existence_exclusion_and_shared_direct_limits',all(s in t for s in ['existence plus common-cube uniqueness','point-map escape alone could not','not a wholly independent physics implementation','confirmed novelty','completion of the original five gates','experimentally detectable tissue displacement']))
ck('point_scope_and_new_successor',all(s in pt for s in ['It does **not** prove','It is not the un-discretized elasticity PDE','subsequent [uniform integral-root theorem](EXACT-INTEGRAL-ROOT.md)','Historical unsafe receipts remain withdrawn']))
ck('historical_warning_current_links',old.splitlines()[2].startswith('**HISTORICAL UNSAFE RECEIPTS REMAIN WITHDRAWN.**')and'(EXACT-INTEGRAL-ROOT.md)'in old.splitlines()[2]and'(SAFE-INTEGRATION-REPLAY.md)'in old.splitlines()[2])
ck('README_current_no_blanket_pending',all(s in rm for s in ['Arithmetic correction completed','Exactly integrated finite-element equilibrium certified','fixed 611-coordinate model','full five-gate study remains incomplete'])and not any(s in rm for s in ['Correction: integration certificates withdrawn pending replay','WITHDRAWN pending replay — proposed integration']))
ck('child_index_current',all(s in ci for s in ['Corrected arithmetic and an exactly integrated equilibrium','EXACT-INTEGRAL-ROOT.html','historical withdrawn receipts','original five gates remain incomplete']))
ck('parent_index_current',all(s in pi for s in ['Corrected arithmetic and an exact finite-element equilibrium','controlled_furrowing/docs/EXACT-INTEGRAL-ROOT.html','Historical unsafe receipts remain withdrawn']))
for name,k,a in [('error','root_error_from_translated_predictor_upper_exact','2.556053177e-12'),('curvature','uniform_integral_energy_curvature_lower_exact','5.305244640e-4'),('location','integral_root_distance_from_old_predictor_lower_exact','1.418087117e-11'),('IQgap','I_Q_branch_component_separation_lower_exact','1.359006738e-11')]:ck('advertised_'+name,abs(Q(cert[k])-Q(a))<=abs(Q(a))*Q(1,10**8))
for name,k,a in [('z','z','.1050417212'),('selfmap','selfmap_lhs','3.337978164e-12')]:ck('advertised_'+name,abs(Q(cert['physical_axis_test'][k])-Q(a))<=abs(Q(a))*Q(1,10**8))
ck('README_outward_error_bound',Q(cert['root_error_from_translated_predictor_upper_exact'])<Q('2.557e-12'))
pages=[ROOT/'docs/EXACT-INTEGRAL-ROOT.html',ROOT/'docs/SAFE-INTEGRATION-REPLAY.html']
figs=[PUB/'preparation/root_margins.png',PUB/'preparation/root_location.png',ROOT/'runs/exact_integration_bridge_v1/safe_v2_publication_preparation/safe_point_restoration.png']
files=[rootmd,pointmd,staged,ROOT/'docs/EXACT-INTEGRATION-BRIDGE.md',PROJECT/'README.md',ROOT/'docs/index.html',PROJECT/'docs/index.html',ROOT/'scripts/build_report.py',ROOT/'scripts/verify_report.py',PUB/'local_publication.json']+pages+figs
before={p.relative_to(PROJECT).as_posix():sha(p)for p in files}
for source in pages:
 parser=V.Links();parser.feed(source.read_text());missing=[]
 for link in parser.links:
  u=V.urlsplit(link)
  if not u.scheme and not u.netloc and u.path and not(source.parent/V.unquote(u.path)).resolve().exists():missing.append(link)
 ck('all_local_links_'+source.stem,not missing)
server=ThreadingHTTPServer(('127.0.0.1',0),partial(V.Quiet,directory=str(PROJECT)));threading.Thread(target=server.serve_forever,daemon=True).start();rows=[]
try:
 with sync_playwright()as pw:
  browser=pw.chromium.launch(channel='msedge',headless=True)
  for width in [390,1440]:
   for source in pages:
    page=browser.new_page(viewport={'width':width,'height':1000});errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
    response=page.goto(f'http://127.0.0.1:{server.server_port}/'+source.relative_to(PROJECT).as_posix(),wait_until='networkidle')
    data=page.evaluate(r"""() => {const prose=document.body.cloneNode(true);prose.querySelectorAll('.math,script,style,pre,code').forEach(x=>x.remove());let strict=[];document.querySelectorAll('[data-tex]').forEach((e,i)=>{try{katex.renderToString(e.dataset.tex,{displayMode:e.classList.contains('display'),throwOnError:true,strict:'error'})}catch(x){strict.push({index:i,error:String(x)})}});return {math:document.querySelectorAll('[data-tex]').length,rendered:document.querySelectorAll('.katex').length,mathErrors:document.getElementById('math-errors')?.textContent||'',strictErrors:strict,rawProseTeX:[...new Set(prose.textContent.match(/\\[A-Za-z]+/g)||[])],overflow:document.documentElement.scrollWidth>innerWidth+1,images:[...document.images].map(i=>({src:i.getAttribute('src'),complete:i.complete,width:i.naturalWidth,height:i.naturalHeight})),bodyText:document.body.innerText}}""")
    for label,value in [('HTTP200',response.status==200),('no_JS_errors',not errors),('math_all_rendered',data['math']>0 and data['math']==data['rendered']),('strict_KaTeX',not data['mathErrors']and not data['strictErrors']),('no_raw_TeX',not data['rawProseTeX']),('no_overflow',not data['overflow']),('all_images_loaded',bool(data['images'])and all(i['complete']and i['width']>0 and i['height']>0 for i in data['images']))]:ck(f'DOM_{source.stem}_{width}_{label}',value)
    body=data.pop('bodyText');ck(f'DOM_{source.stem}_{width}_claim_visible','finite-element'in body and('1,500'in body if source.stem=='EXACT-INTEGRAL-ROOT'else'outside the old certified cube'in body))
    if source.stem=='EXACT-INTEGRAL-ROOT':ck(f'DOM_root_{width}_two_new_figures',len(data['images'])==2 and{Path(i['src']).name for i in data['images']}=={'root_margins.png','root_location.png'})
    screen=OUT/f'{source.stem}_{width}.png';page.screenshot(path=str(screen),full_page=True);data.update({'page':source.name,'width':width,'browser_errors':errors,'screenshot_file':screen.relative_to(ROOT).as_posix(),'screenshot_sha256':sha(screen),'body_text_sha256':hashlib.sha256(body.encode()).hexdigest()});rows.append(data);page.close()
  browser.close()
finally:server.shutdown();server.server_close()
after={p.relative_to(PROJECT).as_posix():sha(p)for p in files};ck('snapshot_unchanged_readonly',before==after)
result={'status':'PASS','checks':checks,'review_source_sha256':sha(__file__),'publication_receipt_sha256':sha(PUB/'local_publication.json'),'staged_MD_sha256':sha(staged),'snapshot':before,'DOM':rows,'browser':'Microsoft Edge Playwright headlessHTTP390/1440','scope':'Read-only local scientific-claim consistency,staged/source/rendered identity,math/assets/layout. No repeated physics/proofs or deployment/fullgoal claim.','point_chapter_scope_note':'Restricted point-proof stage retains its own next-operation account; finalparagraph explicitly links subsequent uniform-root success. Historical unsafe chapter is labelled historical at top.','production_files_modified':False,'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
write(OUT/'results.json',result);write(RECEIPT,result)
print(json.dumps({'status':'PASS','checks':len(checks),'receipt_sha256':sha(RECEIPT),'math_counts':[(r['page'],r['width'],r['math'])for r in rows]}),flush=True)
