"""Verify the one updated parent landing page; does not revalidate parent science."""
from pathlib import Path
import hashlib
import json
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[2]
page_file=ROOT/'docs/index.html'
checks=[]
with sync_playwright() as p:
    browser=p.chromium.launch(channel='msedge',headless=True)
    for width in (390,1440):
        page=browser.new_page(viewport={'width':width,'height':1000});errors=[]
        page.on('pageerror',lambda error:errors.append(str(error)))
        page.goto(page_file.as_uri(),wait_until='networkidle')
        result=page.evaluate("""() => ({mathError:document.querySelector('#latex-error').textContent,
            allImages:[...document.images].every(x=>x.complete&&x.naturalWidth>0),
            overflow:document.documentElement.scrollWidth>innerWidth+1,
            jointLink:[...document.links].some(x=>x.getAttribute('href')==='../controlled_furrowing/docs/JOINT-INFERENCE-PROGRESS.html'),
            currentText:document.body.innerText.includes('all 68 distinct nonlinear candidate targets have completed independent verification')})""")
        passed=not errors and not result['mathError'] and result['allImages'] and not result['overflow'] and result['jointLink'] and result['currentText']
        checks.append({'width':width,'pass':passed,'errors':errors,**result});page.close()
    browser.close()
receipt={'status':'PASS' if all(c['pass'] for c in checks) else 'FAIL','checks':checks,
         'page_sha256':hashlib.sha256(page_file.read_bytes()).hexdigest(),
         'builder_sha256':hashlib.sha256((ROOT/'scripts/build_dossier.py').read_bytes()).hexdigest(),
         'scope':'Updated parent landing page only in headless Edge/file view; not the unchanged parent scientific pages, public deployment or research completion.'}
(ROOT/'controlled_furrowing/evidence/parent_index_checkpoint_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
if receipt['status']!='PASS':raise SystemExit(1)
