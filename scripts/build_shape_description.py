"""Compile the canonical shape guide to static accessible mathematics."""
from pathlib import Path
from hashlib import sha256
import html
import importlib.util
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT/'docs'


def main():
    renderer = ROOT/'controlled_furrowing/scripts/build_report.py'
    spec = importlib.util.spec_from_file_location('_shape_markdown', renderer)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source = DOC/'SHAPE-DESCRIPTION.md'
    text = source.read_text(encoding='utf8')
    body = module.convert(text)
    css = '''*{box-sizing:border-box}html{color-scheme:dark}body{margin:0;background:#0d1828;color:#eaf1fb;font:17px/1.72 system-ui,sans-serif;overflow-wrap:anywhere}main{max-width:1050px;margin:auto;padding:34px clamp(16px,4vw,48px) 70px}a{color:#8ce2e3}h1{font-size:clamp(2rem,4.7vw,3.2rem);line-height:1.15}h2{color:#b3e7f3;margin-top:2em;line-height:1.3}.math.display{max-width:100%;overflow:auto;background:#15293f;padding:14px;border-radius:8px}img{display:block;max-width:100%;height:auto}table{display:block;max-width:100%;overflow:auto;border-collapse:collapse}th,td{border:1px solid #455a73;padding:12px;text-align:left;vertical-align:top}nav{display:flex;flex-wrap:wrap;gap:12px 24px;border-bottom:1px solid #455a73;padding-bottom:18px}li{margin:.5em 0}code{overflow-wrap:anywhere}@media(max-width:600px){body{font-size:16px}main{padding-top:22px}.katex{font-size:1.04em}}'''
    page = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            '<title>How we describe tissue shape mathematically</title>'
            '<link rel="stylesheet" href="assets/katex/katex.min.css"><style>'+css+'</style>'
            '</head><body><main><nav><a href="index.html">Scientific dossier</a>'
            '<a href="SHAPE-DESCRIPTION.md">Guide source</a><a href="../AGENTS.md">Documentation rule</a>'
            '</nav>'+body+'</main></body></html>')
    pattern = re.compile(r'(<(span|div) class="math(?: display)?" data-tex="([^"]*)">)\s*</\2>', re.S)
    matches = list(pattern.finditer(page))
    jobs = [{'tex': html.unescape(m[3]), 'display': m[2] == 'div'} for m in matches]
    js = r'''const fs=require('node:fs');const k=require(process.argv[1]);const xs=JSON.parse(fs.readFileSync(0,'utf8'));process.stdout.write(JSON.stringify(xs.map(x=>k.renderToString(x.tex,{displayMode:x.display,throwOnError:true,strict:'error',output:'htmlAndMathml'}))));'''
    result = subprocess.run(['node', '-e', js, str(DOC/'assets/katex/katex.min.js')],
                            input=json.dumps(jobs), encoding='utf8', capture_output=True, check=True)
    rendered = iter(json.loads(result.stdout))
    page = pattern.sub(lambda m: m[1]+next(rendered)+'</'+m[2]+'>', page)
    expected = len(re.findall(r'\\\((.*?)\\\)|\\\[(.*?)\\\]', text, re.S))
    assert expected == len(jobs)
    target = DOC/'SHAPE-DESCRIPTION.html'
    target.write_text(page, encoding='utf8', newline='\n')
    inputs = [source, target, DOC/'assets/shape-description-chain.svg',
              ROOT/'AGENTS.md', Path(__file__), ROOT/'scripts/build_shape_visuals.py',
              DOC/'SHAPE-VISUALS-RECEIPT.json']
    for stem in ['shape-maps-concrete', 'shape-features-concrete', 'shape-metric-versus-curvature']:
        inputs.extend([DOC/'assets'/(stem+'.png'), DOC/'assets'/(stem+'.svg')])
    inputs.extend([DOC/'assets/shape-features-concrete.csv',
                   DOC/'assets/shape-maps-concrete-data.npz'])
    receipt = {'status': 'PASS_STRICT_SHAPE_GUIDE_BUILD',
               'scope': 'Source expression coverage and strict static TeX compilation; explanatory mathematics, not a research certificate.',
               'math_expressions': len(jobs), 'javascript_required': False,
               'sha256': {str(p.relative_to(ROOT)).replace('\\', '/'): sha256(p.read_bytes()).hexdigest()
                          for p in inputs}}
    (DOC/'SHAPE-DESCRIPTION-BUILD.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf8')
    print(json.dumps({'status': receipt['status'], 'math_expressions': len(jobs)}))


if __name__ == '__main__':
    main()
