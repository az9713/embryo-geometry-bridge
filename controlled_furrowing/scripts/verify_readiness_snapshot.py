"""Verify report coverage and unchanged source bytes across the browser check.

This supplements, rather than edits, the browser verifier. It cannot establish
scientific correctness; its receipt records exactly which document snapshot
was rendered and detects concurrent document edits or omitted viewports.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
FROZEN = {
    'axisymmetric_mechanics.py':'ca431e30a2574a02f4a21c98935c50217aa79bc5512cf7279143f22ce1a580e0',
    'axisymmetric_prestress.py':'141b4b2ea92f4e04f4076b18edc9bafd6218f75a7c4c9633e66c5d2a9fc03c41',
    'axisymmetric_mixed.py':'daf9b4f52181b18b1a75cc69909444bf2a85725c3896b32b125b2f20f10b1f14',
    'axisymmetric_quadratic.py':'5bb64d6b441f6db8ada9da0f1bcb44716c86ed5d816e44f75bd4e9e021fdf410',
    'robust_quadratic.py':'be6f17e620ddf04533d456fff9d48cd48517f4a0709d921ab510616926b4176b',
    'quadratic_quadrature_study.py':'0eac4bdcfd2bb0b725940d15874fc762c16065d2a29f2a10d0da70c5725ac5dc',
    'sparse_quadratic.py':'d4833318151c7b94b4c56d40a96158d55a06258619ef080548d489e2be0394f9',
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot():
    files = list((ROOT/'docs').glob('*.md')) + list((ROOT/'docs').glob('*.html'))
    files += [ROOT/'scripts'/n for n in FROZEN]
    files += [ROOT/'scripts'/n for n in ('build_report.py','verify_report.py',
              'verify_readiness_snapshot.py','physical_parameter_domain.py')]
    return {p.relative_to(ROOT).as_posix():digest(p) for p in sorted(files)}


def main():
    before = snapshot()
    start = time.time()
    proc = subprocess.run([sys.executable,'-B',str(ROOT/'scripts/verify_report.py')],
                          cwd=ROOT.parent, capture_output=True, text=True)
    after = snapshot()
    browser = json.loads((ROOT/'evidence/report_browser_verification.json').read_text())
    pages = {p.name for p in (ROOT/'docs').glob('*.html')}
    expected = {(p,w) for p in pages for w in (390,1440)}
    actual = [(r['page'],r['width']) for r in browser['checks']]
    checks = {
        'browser_process_exit_zero':proc.returncode == 0,
        'browser_report_pass':browser['status'] == 'PASS',
        'all_browser_rows_pass':all(r['pass'] for r in browser['checks']),
        'exact_page_viewport_coverage':set(actual) == expected and len(actual) == len(expected),
        'documents_and_verifiers_unchanged_during_browser_check':before == after,
        'all_frozen_physical_and_policy_sources_match':all(
            before['scripts/'+name] == value for name,value in FROZEN.items()),
    }
    receipt = {'status':'PASS' if all(checks.values()) else 'FAIL', 'checks':checks,
               'html_pages':len(pages), 'page_viewport_checks':len(actual),
               'elapsed_seconds':time.time()-start, 'snapshot_sha256':before,
               'browser_stdout':proc.stdout, 'browser_stderr':proc.stderr,
               'browser_receipt_sha256':digest(ROOT/'evidence/report_browser_verification.json'),
               'scope':'Immutable document snapshot, exact browser coverage and preserved model sources; not scientific proof, continuum accuracy or completion of five gates.'}
    path = ROOT/'evidence/readiness_snapshot_verification.json'
    path.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('status','checks','html_pages','page_viewport_checks')}),flush=True)
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
