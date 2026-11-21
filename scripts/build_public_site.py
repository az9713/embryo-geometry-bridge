"""Build the public journey, catalogue and static mathematical presentation."""
from pathlib import Path
from hashlib import sha256
import html
import importlib.util
import json
import os
import re
import subprocess
from urllib.parse import quote

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/'docs'
BASE='https://az9713.github.io/embryo-geometry-bridge/'
WAVE='controlled_furrowing/runs/transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1/remaining32_bridge_v1/remaining28_preparation_v1/waves/wave_03/checkpoint_v1/'
PATTERN=re.compile(r'(<(span|div) class="math(?: display)?" data-tex="([^"]*)">)\s*</\2>',re.S)
CSS='''*{box-sizing:border-box}html{color-scheme:dark}body{margin:0;background:#0d1828;color:#e9f0fa;font:17px/1.7 system-ui,sans-serif;overflow-wrap:anywhere}main{max-width:1040px;margin:auto;padding:32px clamp(18px,5vw,58px) 70px}nav{display:flex;gap:12px 25px;flex-wrap:wrap;border-bottom:1px solid #3d536c;padding-bottom:20px}a{color:#92e0e4;text-underline-offset:.2em}h1{font-size:clamp(2.1rem,5vw,3.5rem);line-height:1.1;letter-spacing:-.035em}h2{margin-top:2.1em;color:#b7e8f3;line-height:1.25}h3{line-height:1.3}p{max-width:82ch}article,.notice,blockquote{background:#172a40;border:1px solid #425a75;border-radius:10px;padding:18px 23px;margin:24px 0}.notice{border-left:4px solid #f2c76e}blockquote{border-left:4px solid #8cdbdf}img{display:block;max-width:100%;height:auto;border-radius:8px}.math.display,pre{max-width:100%;overflow:auto;background:#16283c;padding:10px 14px;border-radius:7px}table{display:block;max-width:100%;overflow:auto;border-collapse:collapse;font-size:.92em}th,td{border:1px solid #3b526b;padding:11px 13px;text-align:left;vertical-align:top}.cards{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}.cards article{margin:0}.cards strong{font-size:1.5rem}.buttons{display:flex;gap:12px;flex-wrap:wrap}.button{display:inline-block;padding:10px 15px;border:1px solid #75c8d1;border-radius:7px;text-decoration:none}li{margin:.45em 0}code{overflow-wrap:anywhere}.small{font-size:.9rem;color:#becbdb}@media(max-width:560px){body{font-size:16px}.cards{grid-template-columns:1fr}main{padding-top:20px}article{padding:15px}.katex{font-size:1.05em}}'''
NODE=r'''const fs=require('node:fs');const k=require(process.argv[1]);const xs=JSON.parse(fs.readFileSync(0,'utf8'));process.stdout.write(JSON.stringify(xs.map(x=>{try{return {html:k.renderToString(x.tex,{displayMode:x.display,throwOnError:true,strict:'error',output:'htmlAndMathml'})}}catch(e){return {error:String(e),id:x.id}}})));'''

def rel(p,base): return Path(os.path.relpath(p,base)).as_posix()
def pages(): return sorted(p for p in ROOT.rglob('*.html') if '.git' not in p.relative_to(ROOT).parts)
def url(p): return BASE+quote(p.relative_to(ROOT).as_posix(),safe='/')
def shell(title,body):
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
            '<title>'+html.escape(title)+'</title><link rel="stylesheet" href="assets/katex/katex.min.css"><style>'+CSS+'</style></head><body><main>'
            '<nav><a href="index.html">Scientific dossier</a><a href="DEVELOPMENT-JOURNEY.html">Development journey</a>'
            '<a href="ALL-PAGES.html">All web pages</a><a href="../README.md">README source</a></nav>'+body+'</main></body></html>')

def convert(md):
    p=ROOT/'controlled_furrowing/scripts/build_report.py'
    spec=importlib.util.spec_from_file_location('_publication_markdown_converter',p)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    body=module.convert(md)
    # Stable anchors for browsing the journal.
    def heading(m):
        text=re.sub('<[^>]+>','',m[2]);slug=re.sub(r'[^a-z0-9 -]','',html.unescape(text).lower()).strip().replace(' ','-')
        return '<h'+m[1]+' id="'+slug+'">'+m[2]+'</h'+m[1]+'>'
    return re.sub(r'<h([123])>(.*?)</h\1>',heading,body,flags=re.S)

def build():
    md=(DOC/'DEVELOPMENT-JOURNEY.md').read_text(encoding='utf8')
    (DOC/'DEVELOPMENT-JOURNEY.html').write_text(shell('Development journey: Tsimerman to tissue geometry',convert(md)),encoding='utf8',newline='\n')
    archived=DOC/'ARCHIVED-AVIAN-FINDINGS.md'
    if archived.exists():
        body='<div class="notice">Historical avian branch. These archived findings are not the current controlled-folding coverage counts. The complete older simulation record is not included in this curated release.</div>'+convert(archived.read_text(encoding='utf8'))
        (DOC/'ARCHIVED-AVIAN-FINDINGS.html').write_text(shell('Archived findings: earlier avian branch',body),encoding='utf8',newline='\n')
    body='''<h1>Embryo Geometry Bridge</h1><p>From Tsimerman's challenge about scientific ambition to a testable question: what can a specified tissue model tell us about embryonic geometry?</p><p class="buttons"><a class="button" href="DEVELOPMENT-JOURNEY.html">Read the development journey</a><a class="button" href="../'''+WAVE+'''REPORT.html">Proofs and verification</a><a class="button" href="ALL-PAGES.html">Browse every HTML page</a></p><div class="notice"><strong>Research in progress.</strong> No validated biological shape forecast or confirmed new general scientific concept is claimed.</div><div class="cards"><article><strong>16 / 36</strong><p>Selected fields have verified exact-spatial finite-element roots and minima within the axisymmetric displacement space.</p></article><article><strong>56 primary comparisons</strong><p>All pass under the declared discrepancy allowance. Alternative observation conventions produce 13 violations among 168 comparisons.</p></article><article><strong>Four loading slabs</strong><p>Body geometry guards pass across the original K-L interval. The mechanical loading branch remains uncertified.</p></article><article><strong>A symmetry limitation</strong><p>Compressive surface patches permit energy descent in the stated bending-free 3D extension. Physically resolved wavelengths and continuum stability remain open.</p></article></div><h2>What the project connects</h2><img src="assets/research-chain.svg" alt="Molecular activity, active stress, mechanics, tissue deformation and geometry each require a measured or justified connection."><p>The current calculations test portions of a mechanics-to-geometry chain. Missing biological observations and constitutive calibration remain separate requirements.</p><p><a href="../PUBLICATION_MANIFEST.json">File manifest</a> · <a href="../PUBLICATION_OMISSIONS.json">Omitted artifacts and replay limits</a> · <a href="../verification/browser/receipt.json">Earlier publication rendering receipt</a></p>'''
    (DOC/'index.html').write_text(shell('Embryo Geometry Bridge: scientific dossier',body),encoding='utf8',newline='\n')
    (ROOT/'index.html').write_text('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta http-equiv="refresh" content="0; url=docs/index.html"><title>Embryo Geometry Bridge</title></head><body><a href="docs/index.html">Open the scientific dossier</a></body></html>',encoding='utf8',newline='\n')
    # Include the catalogue itself on its first build.
    catalog=DOC/'ALL-PAGES.html'; allpages=sorted(set(pages())|{catalog})
    items=[]
    for p in allpages:
        text=p.read_text(encoding='utf8') if p.exists() else '<title>All HTML pages</title>'
        title=re.search(r'<title>(.*?)</title>',text,re.S|re.I)
        title=html.unescape(title[1]).strip() if title else p.stem
        items.append('<li><a href="'+html.escape(rel(p,DOC),quote=True)+'">'+html.escape(title)+'</a><br><span class="small">'+html.escape(p.relative_to(ROOT).as_posix())+'</span></li>')
    catalog.write_text(shell('All HTML pages: Embryo Geometry Bridge','<h1>All '+str(len(allpages))+' HTML pages</h1><p>Every link opens a rendered web page on this static site. Older pages are dated checkpoints; use the latest dossier for current scope. Omission notices identify artifacts outside the public bundle.</p><ol>'+''.join(items)+'</ol>'),encoding='utf8',newline='\n')
    (DOC/'PAGE-CATALOG.json').write_text(json.dumps({'site':BASE,'html_pages':[{'path':p.relative_to(ROOT).as_posix(),'url':url(p)} for p in allpages]},indent=2)+'\n',encoding='utf8')
    readme=ROOT/'README.md';text=readme.read_text(encoding='utf8')
    text=re.sub(r'<!-- live-pages:start -->.*?<!-- live-pages:end -->\s*','',text,flags=re.S)
    text=text.replace('[Open the scientific dossier](docs/index.html)', '[Open the live scientific dossier]('+BASE+'docs/index.html)')
    text=text.replace('](controlled_furrowing/runs/transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1/remaining32_bridge_v1/remaining28_preparation_v1/waves/wave_03/checkpoint_v1/REPORT.html)', ']('+BASE+WAVE+'REPORT.html)')
    block='<!-- live-pages:start -->\n[**Development journey: beginning with Tsimerman\'s video**]('+BASE+'docs/DEVELOPMENT-JOURNEY.html) · [**Every HTML page**]('+BASE+'docs/ALL-PAGES.html)\n\n[![Rendered scientific dossier](docs/assets/dossier-preview.png)]('+BASE+'docs/index.html)\n\n[![Rendered development journey](docs/assets/journey-preview.png)]('+BASE+'docs/DEVELOPMENT-JOURNEY.html)\n\nGitHub displays clickable image previews here. Follow them to the fully rendered GitHub Pages documents.\n\n<details>\n<summary>Direct webpage links to all '+str(len(allpages))+' HTML files</summary>\n\n'
    for p in allpages:
        block+='- ['+p.relative_to(ROOT).as_posix()+']('+url(p)+')\n'
    block+='\n</details>\n<!-- live-pages:end -->\n\n'
    marker='## Findings and significance';assert marker in text;text=text.replace(marker,block+marker,1)
    readme.write_text(text,encoding='utf8',newline='\n')
    return allpages

def compile_all():
    jobs=[]; inputs={};matchsets={};normalizations=[]
    for p in pages():
        text=p.read_text(encoding='utf8');matches=list(PATTERN.finditer(text))
        if not matches:continue
        inputs[p]=text;matchsets[p]=matches
        for i,m in enumerate(matches):
            tex=html.unescape(m[3]);fixed=tex
            if p.name=='ROOT-OBSERVER-APPLICATION-REVIEW.html':
                fixed=fixed.replace('{\nm \u00b5m}',r'\mathrm{\mu m}')
            fixed=fixed.replace('\u00b5',r'\mu ')
            label=p.relative_to(ROOT).as_posix()+':'+str(i)
            if fixed!=tex:normalizations.append({'expression':label,'before':tex,'after':fixed,'reason':'Legacy unit typography normalized to standard TeX; numerical values unchanged.'})
            jobs.append({'id':label,'tex':fixed,'original_tex':tex,'display':m[2]=='div'})
    if jobs:
        run=subprocess.run(['node','-e',NODE,str(DOC/'assets/katex/katex.min.js')],input=json.dumps(jobs),text=True,encoding='utf8',capture_output=True,check=True,cwd=ROOT)
        rendered=json.loads(run.stdout);errors=[x for x in rendered if 'error' in x]
        if errors:raise ValueError(json.dumps(errors,ensure_ascii=True))
        assert len(rendered)==len(jobs);cursor=iter(zip(rendered,jobs,strict=True))
        for p,text in inputs.items():
            # Retain source TeX; existing runtime renderers may harmlessly rerender.
            def replace(m):
                result,job=next(cursor)
                opening=m[1]
                if job['tex']!=job['original_tex']:
                    opening=opening.replace('data-tex="'+m[3]+'"','data-tex="'+html.escape(job['tex'],quote=True)+'"')
                return opening+result['html']+'</'+m[2]+'>'
            updated=PATTERN.sub(replace,text)
            p.write_text(updated,encoding='utf8',newline='\n')
    receipt={'status':'PASS_STRICT_STATIC_MATH_COMPILATION','newly_compiled_expressions':len(jobs),'pages_compiled':len(inputs),
             'math_source_retained_except_logged_unit_typography':True,'unit_typography_normalizations':normalizations,
             'throw_on_error':True,'strict':'error','scientific_numerical_files_changed':False}
    (DOC/'STATIC-MATH-VERIFICATION.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    return receipt

if __name__=='__main__':
    pp=build();result=compile_all();print(json.dumps({'html_pages':len(pp),**result}))
