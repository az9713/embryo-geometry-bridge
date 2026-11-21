"""Independent rendering fixtures; original sources and scientific results read-only."""
from pathlib import Path
from http.server import SimpleHTTPRequestHandler,ThreadingHTTPServer
from functools import partial
import ast,difflib,hashlib,importlib.util,importlib.metadata,json,re,threading
from playwright.sync_api import sync_playwright
RUN=Path(__file__).resolve().parent
ROOT=RUN.parents[1]
STUDY=ROOT.parent
OUT=ROOT/'evidence/rendering_guard_independent_review.json'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
rel=lambda p:Path(p).resolve().relative_to(ROOT.resolve()).as_posix()

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def probe(path):
    v=[n.value for n in ast.walk(ast.parse(path.read_text())) if isinstance(n,ast.Constant) and isinstance(n.value,str) and 'document.querySelectorAll' in n.value and 'rendered:' in n.value]
    if len(v)!=1:raise ValueError('unique exactbrowserprobe required')
    return v[0]

def passes(v,raw):return not v['errors'] and v['math']==v['rendered'] and v['images'] and not v['overflow'] and (not v.get('rawProseTeX') if raw else True)
class Quiet(SimpleHTTPRequestHandler):
    def log_message(self,*args):pass

def main():
    if OUT.exists():raise RuntimeError('preserve completed rendering review')
    arc=ROOT/'runs/rendering_guard_revision_0';paths={'old_builder':arc/'build_report.py','new_builder':ROOT/'scripts/build_report.py','old_verifier':arc/'verify_report.py','new_verifier':ROOT/'scripts/verify_report.py','old_review':arc/'shared_cell_cohort_independent_review.json'};sources={k:sha(v) for k,v in paths.items()}
    snap=RUN/'source_snapshots';snap.mkdir()
    for k,v in paths.items():(snap/(k+v.suffix)).write_bytes(v.read_bytes())
    old=module('old_builder_fixture',paths['old_builder']);new=module('new_builder_fixture',paths['new_builder']);before=probe(paths['old_verifier']);after=probe(paths['new_verifier'])
    (RUN/'old_probe.js').write_text(before,encoding='utf8');(RUN/'new_probe.js').write_text(after,encoding='utf8')
    ast_without_convert=lambda p:ast.dump(ast.Module(body=[n for n in ast.parse(p.read_text()).body if not isinstance(n,ast.FunctionDef) or n.name!='convert'],type_ignores=[]),include_attributes=False)
    checks={'builder_AST_unchanged_except_convert':ast_without_convert(paths['old_builder'])==ast_without_convert(paths['new_builder'])}
    for kind in ('builder','verifier'):(RUN/(kind+'_delta.diff')).write_text(''.join(difflib.unified_diff(paths['old_'+kind].read_text().splitlines(True),paths['new_'+kind].read_text().splitlines(True),fromfile='old_'+kind,tofile='new_'+kind)),encoding='utf8')
    md='| Name | Condition | Note |\n|---|---|---|\n| increment | \\(|s|\\le0.001\\) | kept |\n';expected='|s|\\le0.001';(RUN/'table_fixture.md').write_text(md,encoding='utf8');newbody=new.convert(md)
    checks['pipe_entity_present_before_tableparse']='&#124;s&#124;' in newbody
    standalone=r'Inline \(|x|=1\).'+ '\n\n'+r'\[\|x\|_2=1\]'
    code=new.convert('```text\n\\rho \\frac{1}{2}\n```\n\n`\\alpha`\n')+r'<style>/* \sigma */</style><script type="application/json">{"sample":"\\beta"}</script>'
    bodies={'new_table':newbody,'old_table':old.convert(md),'raw_prose':r'<p>Unrendered \rho and \frac{1}{2}.</p>','code_exempt':code,'standalone':new.convert(standalone),'guard_limit':r'<p>Plain x_1+y, and \(x+y\) without a control word.</p>'}
    render=re.findall(r'<script>(.*?)</script>',new.shell('fixture',''),re.S)
    if len(render)!=1:raise ValueError('unique exactrenderer inline script required')
    css='body{background:#0b111b;color:#e7eef7;font:20px/1.6 system-ui;padding:30px;margin:0}table{border-collapse:collapse;width:100%;font-size:18px}td,th{padding:15px;border:1px solid #657286}pre{white-space:pre-wrap}#math-errors{color:#ff9999}h1{font-size:26px}'
    for name,body in bodies.items():(RUN/(name+'.html')).write_text('<!doctype html><html><head><meta charset="utf-8"><link rel="stylesheet" href="/docs/assets/katex/katex.min.css"><style>'+css+'</style></head><body><h1>'+name+'</h1>'+body+'<p id="math-errors"></p><script src="/docs/assets/katex/katex.min.js"></script><script>'+render[0]+'</script></body></html>',encoding='utf8')
    server=ThreadingHTTPServer(('127.0.0.1',0),partial(Quiet,directory=str(STUDY)));threading.Thread(target=server.serve_forever,daemon=True).start();results=[]
    try:
        with sync_playwright() as pw:
            browser=pw.chromium.launch(channel='msedge',headless=True);version=browser.version
            for width in (1440,390):
                page=browser.new_page(viewport={'width':width,'height':700})
                for name in bodies:
                    errors=[];badhttp=[];onerror=lambda x:errors.append(str(x));onresponse=lambda x:badhttp.append(x.url) if x.status>=400 else None;page.on('pageerror',onerror);page.on('response',onresponse)
                    response=page.goto(f'http://127.0.0.1:{server.server_port}/'+(RUN/(name+'.html')).relative_to(STUDY).as_posix(),wait_until='networkidle');b=page.evaluate(before);a=page.evaluate(after)
                    extra=page.evaluate(r"""() => ({dataTex:[...document.querySelectorAll('[data-tex]')].map(x=>x.dataset.tex),headCells:[...document.querySelectorAll('thead th')].map(x=>x.textContent),bodyRows:[...document.querySelectorAll('tbody tr')].map(x=>[...x.querySelectorAll('td')].map(y=>y.textContent)),codeText:[...document.querySelectorAll('pre,code')].map(x=>x.textContent),annotations:[...document.querySelectorAll('annotation[encoding="application/x-tex"]')].map(x=>x.textContent),katexVersion:katex.version})""")
                    results.append({'fixture':name,'width':width,'HTTP_status':response.status,'browser_errors':errors,'HTTP_errors':badhttp,'old_guard':b,'new_guard':a,'DOM':extra});tag=name+'_'+str(width)
                    checks[tag+'_local_assets_no_browser_errors']=response.status==200 and not errors and not badhttp
                    if name=='new_table':
                        checks[tag+'_one_mathnode_exact_dataset']=a['math']==1 and extra['dataTex']==[expected]
                        checks[tag+'_exact_three_header_bodycells']=extra['headCells']==['Name','Condition','Note'] and len(extra['bodyRows'])==1 and len(extra['bodyRows'][0])==3 and extra['bodyRows'][0][0]=='increment' and extra['bodyRows'][0][2]=='kept'
                        checks[tag+'_KaTeX_success_exact_annotation']=a['rendered']==1 and extra['annotations']==[expected] and passes(a,True)
                    elif name=='old_table':checks[tag+'_old_format_failure_reproduced']=extra['dataTex']!=[expected] or a['math']!=1
                    elif name=='raw_prose':
                        checks[tag+'_old_count_guard_misses_raw']=passes(b,False)
                        checks[tag+'_new_guard_rejects_raw']=set(a['rawProseTeX'])=={'\\rho','\\frac'} and not passes(a,True)
                    elif name=='code_exempt':checks[tag+'_code_style_script_exempt']=not a['rawProseTeX'] and passes(a,True) and any('\\rho \\frac{1}{2}' in s for s in extra['codeText']) and any('\\alpha' in s for s in extra['codeText'])
                    elif name=='standalone':checks[tag+'_inline_display_regression']=extra['dataTex']==['|x|=1','\\|x\\|_2=1'] and a['rendered']==2 and passes(a,True)
                    else:checks[tag+'_heuristic_limit_explicit']=not a['rawProseTeX']
                    if name in ('new_table','old_table'):page.screenshot(path=str(RUN/(tag+'.png')),full_page=True)
                    page.remove_listener('pageerror',onerror);page.remove_listener('response',onresponse)
                page.close()
            browser.close()
    finally:server.shutdown();server.server_close()
    historical=read(paths['old_review']);changed=[k for k,v in historical['file_hashes'].items() if sha(STUDY/k)!=v]
    checks['old180_receipt_bytes_unchanged']=sha(ROOT/'evidence/shared_cell_cohort_independent_review.json')==sources['old_review']=='2f7da9f2bee4c1b7e4fdd7f3a88ed90dd914fb218269b3d65bf12bb8a6ab9114'
    checks['ONLY_old_builder_binding_superseded']=changed==['controlled_furrowing/scripts/build_report.py']
    checks['archived_builder_equals_old180_binding']=historical['file_hashes'].get('controlled_furrowing/scripts/build_report.py')==sources['old_builder']
    science=read(ROOT/'runs/shared_cell_cohort_v1/protocol.json');checks['all71_science_inputs_unchanged']=len(science['file_hashes'])==71 and all(sha(ROOT/n)==v for n,v in science['file_hashes'].items())
    con=read(ROOT/'runs/cohort_static_observers_v1/model_observation_contract.json');checks['all18_physical_observer_sources_unchanged']=len(con['source_kernel_hashes'])==18 and all(sha(ROOT/'scripts'/n)==v for n,v in con['source_kernel_hashes'].items());checks['tested_render_sources_unchanged']=all(sha(v)==sources[k] for k,v in paths.items())
    result={'status':'PASS' if all(v is True for v in checks.values()) else 'FAIL','checks':checks,'review_source_sha256':sha(__file__),'source_hashes':sources,'source_files':{k:rel(v) for k,v in paths.items()},'browser':'Microsoft Edge, Playwright headless','browser_version':version,'playwright_version':importlib.metadata.version('playwright'),'markdown_version':importlib.metadata.version('Markdown'),'fixture_results':results,'fixture_file_hashes':{rel(q):sha(q) for q in RUN.rglob('*') if q.is_file()},'old180_control_source_supersession':{'path':'controlled_furrowing/scripts/build_report.py','old':sources['old_builder'],'new':sources['new_builder'],'reason':'Escape attribute pipes before Markdown; browser restores exact TeX.'},'scope':'Independent localfixtures only, not final102pages/204viewports, livedeployment or newsciencecomputation. Detector matches backslash-letter control words; .math/script/style/pre/code exempt. Plain algebra or control-symbol-only text is not universally detected.','scientific_proofs_or_claims_changed':False,'development_note':'Initial sourcegenerator nested-string syntax failed before writing any review files or executing browser/science; corrected generator. Original production sources were never edited.'}
    OUT.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8');print(json.dumps({'status':result['status'],'checks':len(checks),'failed':[k for k,v in checks.items() if not v],'receipt_sha256':sha(OUT)}))
    if result['status']!='PASS':raise SystemExit(1)
if __name__=='__main__':main()
