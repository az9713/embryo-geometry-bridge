"""Test the new parent landing link and layout; parent science is out of scope."""
from pathlib import Path
import hashlib
import json
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[3]
page_file = ROOT/'docs/index.html'
rows = []
with sync_playwright() as p:
    browser = p.chromium.launch(channel='msedge', headless=True)
    for width in (390, 1440):
        page = browser.new_page(viewport={'width': width, 'height': 1000})
        errors = []
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.goto(page_file.as_uri(), wait_until='networkidle')
        value = page.evaluate("""() => ({
          mathError:document.querySelector('#latex-error')?.textContent||'',
          images:[...document.images].every(x=>x.complete&&x.naturalWidth>0),
          overflow:document.documentElement.scrollWidth>innerWidth+1,
          energyLink:[...document.links].some(x=>x.getAttribute('href')===
            '../controlled_furrowing/docs/COHORT-ENERGY-STABILITY.html'),
          energyText:document.body.innerText.includes('strict local energy minima'),
          scope:document.body.innerText.includes('continuum stability, biological fit and novelty remain unproved')
        })""")
        passed = not errors and not value['mathError'] and value['images'] and not value['overflow'] and value['energyLink'] and value['energyText'] and value['scope']
        rows.append({'width': width, 'pass': passed, 'browser_errors': errors, **value})
        page.close()
    browser.close()
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
out = ROOT/'controlled_furrowing/evidence/cohort_energy_parent_landing_verification.json'
result = {'status': 'PASS' if all(r['pass'] for r in rows) else 'FAIL', 'checks': rows,
          'page_sha256': sha(page_file), 'builder_sha256': sha(ROOT/'scripts/build_dossier.py'),
          'verifier_sha256': sha(Path(__file__)),
          'scope': 'Parent landing only, headless Edge at two widths. No parent science-page or public deployment claim.'}
out.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
print(json.dumps({'status': result['status'], 'viewports': len(rows)}))
if result['status'] != 'PASS':
    raise SystemExit(1)
