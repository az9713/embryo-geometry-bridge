"""Strict static rendering for the executed three-step research record."""
from pathlib import Path
from hashlib import sha256
import importlib.util
import html,json,re,subprocess

HERE=Path(__file__).resolve().parent
MECH=HERE.parents[1]
PROJECT=MECH.parent


def main():
    spec=importlib.util.spec_from_file_location('_unchanged_markdown_renderer',MECH/'scripts/build_report.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    pattern=re.compile(r'(<(span|div) class="math(?: display)?" data-tex="([^"]*)">)\s*</\2>',re.S)
    css='''*{box-sizing:border-box}html{color-scheme:dark}body{margin:0;background:#0e1929;color:#e8f0fb;font:17px/1.72 system-ui,sans-serif;overflow-wrap:anywhere}main{max-width:1100px;margin:auto;padding:32px clamp(16px,4vw,46px) 75px}a{color:#86dfdf}h1{font-size:clamp(2rem,4.5vw,3.2rem);line-height:1.15}h2{color:#b8e7f4;margin-top:2.1em;line-height:1.3}nav{display:flex;gap:12px 24px;flex-wrap:wrap;border-bottom:1px solid #425976;padding-bottom:18px}.notice{background:#1c2e43;border-left:4px solid #e8bd72;padding:16px}.math.display{overflow:auto;max-width:100%;padding:14px;background:#162a40;border-radius:8px}img{display:block;max-width:100%;height:auto}table{display:block;max-width:100%;overflow:auto;border-collapse:collapse;font-size:.92em}td,th{border:1px solid #405771;padding:12px;text-align:left;vertical-align:top}li{margin:.4em 0}code{overflow-wrap:anywhere}@media(max-width:600px){body{font-size:16px}.katex{font-size:1.04em}}'''
    js=r'''const fs=require('node:fs');const k=require(process.argv[1]);const xs=JSON.parse(fs.readFileSync(0,'utf8'));process.stdout.write(JSON.stringify(xs.map(x=>k.renderToString(x.tex,{displayMode:x.display,throwOnError:true,strict:'error',output:'htmlAndMathml'}))));'''
    rows=[]
    for name,title in [('REPORT','Relaxed force balance and finite-height surface recovery'),('THEOREM','Constructive fixed-height periodic slab theorem')]:
        source=HERE/(name+'.md');text=source.read_text(encoding='utf8');body=module.convert(text)
        page='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><link rel="stylesheet" href="../../../docs/assets/katex/katex.min.css"><style>'+css+'</style></head><body><main><nav><a href="../../../docs/SHAPE-DESCRIPTION.html">Shape convention</a><a href="REPORT.html">Executed results</a><a href="THEOREM.html">Full slab proof</a><a href="'+name+'.md">Source</a></nav><p class="notice">Scoped mathematical results. No full cylinder relaxation theorem, biological prediction, sharp-rate or novelty claim.</p>'+body+'</main></body></html>'
        matches=list(pattern.finditer(page));jobs=[{'tex':html.unescape(x[3]),'display':x[2]=='div'} for x in matches]
        r=subprocess.run(['node','-e',js,str(PROJECT/'docs/assets/katex/katex.min.js')],input=json.dumps(jobs),encoding='utf8',capture_output=True,check=True)
        rendered=iter(json.loads(r.stdout));page=pattern.sub(lambda m:m[1]+next(rendered)+'</'+m[2]+'>',page)
        assert len(jobs)==len(re.findall(r'\\\((.*?)\\\)|\\\[(.*?)\\\]',text,re.S))
        target=HERE/(name+'.html');target.write_text(page,encoding='utf8',newline='\n')
        rows.append({'source':source.name,'source_sha256':sha256(source.read_bytes()).hexdigest(),'html':target.name,'html_sha256':sha256(target.read_bytes()).hexdigest(),'strict_math_expressions':len(jobs)})
    receipt={'status':'PASS_STRICT_EXECUTED_REPORT_MATH','rows':rows,'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'document_build.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps({'status':receipt['status'],'math_counts':{x['html']:x['strict_math_expressions'] for x in rows}}))


if __name__=='__main__':main()
