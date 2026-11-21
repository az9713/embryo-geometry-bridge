"""Preserve the prior local publication and verify the scoped navigation update."""
import argparse
import ast
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import urlsplit


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CHILD = 'controlled_furrowing'
RUN = f'{CHILD}/runs/transfer_ambiguity_nonlinear_probe_v1'
L_PUBLICATION = (f'{CHILD}/runs/exact_integration_cohort_v1/L_defect_diagnostic/'
                 'signed_predictor_L_v1/publication')
ARCHIVE = HERE / 'archive'
PLAN = f'{CHILD}/PLAN.md'
CHILD_README = f'{CHILD}/README.md'
READMES = ['README.md'] + ([CHILD_README] if (ROOT / CHILD_README).is_file() else [])
LANDINGS = ('docs/index.html', f'{CHILD}/docs/index.html')
BUILDERS = ('scripts/build_dossier.py', f'{CHILD}/scripts/build_report.py')
TARGETS = (*READMES, *LANDINGS, *BUILDERS)
L_DELIVERY = f'{CHILD}/evidence/L_signed_root_final_delivery.json'
L_READBACK = f'{L_PUBLICATION}/final_delivery_readback.json'
L_BROWSER = f'{L_PUBLICATION}/browser_readback.json'
NONLINEAR_BROWSER = f'{RUN}/report/browser_readback_v2/receipt.json'
REPORT_BUILD = f'{RUN}/report/report_build_receipt.json'
REPORT_SCIENCE = f'{RUN}/report/scientific_readback/receipt.json'
ORIGINALITY = f'{RUN}/originality_followup/receipt.json'
LINKS = (f'{RUN}/report/REPORT.html', f'{RUN}/originality_followup/README.md',
         f'{RUN}/independent_actual_preparation/actual_replay/receipt.json',
         f'{RUN}/observer_robustness/independent_actual/receipt.json',
         REPORT_SCIENCE)
PINS = {
    PLAN: '993240af5d30c8e5b72ab31a1278a524314d38ede48486b677b021d84bcde500',
    L_DELIVERY: '3753330a0b669ef364ff0f0dd3d34f0e0ea3c661a535840ff2b673909bf73162',
    L_READBACK: 'e384b459ed93c539ab26455bb11ea8ed16200b2e3584a3539f0e33dd715eee95',
    L_BROWSER: '56103ea64cdb209de5b218d02b814d74dbeefd232da529588b48b03e324e5be6',
    NONLINEAR_BROWSER: 'c0236dd4a464cd29b2b4925a3232fb4017ca0ddc65952cefb839118005a52d8d',
    REPORT_BUILD: '79cae6792ea8fbfef89bd9d69a3c87f14f419fb1a6c03384f781caf92eeff7c9',
    LINKS[0]: 'd33f3d594406ed498d364c3462450b62cb27054f0146f4bbca3a1c2ae22b41f6',
    LINKS[1]: 'be18025c853273ebed98db4cbb47dac48006236631bee671090d22f279410c43',
    LINKS[2]: 'c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd',
    LINKS[3]: '703b574bc51a083b10c7328d9622862e939ddabe791b18ce9e1817f07718f69f',
    LINKS[4]: '42a40beabd178385bc14afa5d15106f43d390883427413c90474cd34eafc33e3',
    ORIGINALITY: '6797ad5d05f71535f13edc519276e57bc5176d7d5b0e6273b1aa91eba48653d0',
}


def require(check, message):
    if not check:
        raise ValueError(message)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def save(path, data):
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write('\n')
    return sha(path)


def check_bytes(snapshot):
    for name, expected in snapshot.items():
        require(sha(ROOT / name) == expected, f'protected bytes changed: {name}')


def page_snapshot():
    pages = [ROOT / 'index.html']
    for folder in ('docs', f'{CHILD}/docs'):
        pages.extend(sorted((ROOT / folder).glob('*.html')))
    return {p.relative_to(ROOT).as_posix(): sha(p) for p in pages}


def scientific_authority():
    check_bytes(PINS)
    primary = read(ROOT / LINKS[2])
    secondary = read(ROOT / LINKS[3])
    science = read(ROOT / LINKS[4])
    original = read(ROOT / ORIGINALITY)
    delta = Q(science['primary_allowance_um_exact'])
    lo, hi = map(Q, science['signed_large_gap_um_exact'])
    pair_thresholds = list(science['observer_pair_allowances_um_exact'].values())
    require(science['status'] == 'PASS_INDEPENDENT_SCIENTIFIC_REPORT_READBACK'
            and primary['status'] == 'PASS_INDEPENDENT_ACTUAL_NONLINEAR_STORED_FIELD_REPLAY'
            and primary['both_endpoints_satisfy_conditional_stored_allowance'] is True
            and secondary['status'] == 'PASS_INDEPENDENT_SECONDARY_STORED_OBSERVER_SWEEP'
            and secondary['robust_across_conventions_at_frozen_delta'] is False
            and len(pair_thresholds) == 4 and sum(Q(p[1]) <= delta for p in pair_thresholds) == 1
            and all(f'{float(q):.6f}' == '0.154755' for q in (lo, hi))
            and f'{float(delta):.4f}' == '5.3898'
            and original['novelty_confirmed'] is False,
            'declared navigation statement drifted from the independently reviewed numbers and scope')
    report = read(ROOT / REPORT_BUILD)
    check_bytes({f'{CHILD}/{k}': v for k, v in report['report_and_figure_hashes'].items()})
    return report


def archive():
    report = scientific_authority()
    delivery = read(ROOT / L_DELIVERY)
    review = read(ROOT / L_READBACK)
    old = delivery['source_and_artifact_hashes']
    require(delivery['status'] == 'PASS_L_ROOT_COMPONENT_DELIVERY_NOT_FULL_GOAL'
            and review['status'] == 'PASS_INDEPENDENT_L_COMPONENT_FINAL_DELIVERY_READBACK'
            and review['delivery_sha256'] == PINS[L_DELIVERY]
            and len(old) == review['verified_current_delivery_bindings'] == 335,
            'prior L delivery must retain 335 exact bindings')
    check_bytes(old)
    pages = page_snapshot()
    require(sum(k.startswith(f'{CHILD}/docs/') for k in pages) == 110,
            'existing controlled-furrowing dossier expected 110 rendered pages')
    ARCHIVE.mkdir(exist_ok=False)
    copied = {}
    for name in (*TARGETS, L_DELIVERY, L_READBACK, L_BROWSER, NONLINEAR_BROWSER, REPORT_BUILD, REPORT_SCIENCE):
        origin, dest = ROOT / name, ARCHIVE / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        payload = origin.read_bytes()
        with dest.open('xb') as handle:
            handle.write(payload)
        require(dest.read_bytes() == payload, f'archive copy failed: {name}')
        copied[name] = {'archive_path': dest.relative_to(ROOT).as_posix(), 'sha256': sha(dest)}
    snapshots = {k: v for k, v in old.items() if k in TARGETS}
    require(all(copied[p]['sha256'] == s for p, s in snapshots.items())
            and copied[BUILDERS[1]]['sha256'] == report['input_hashes']['scripts/build_report.py'],
            'old report converter and L shared navigation bindings must match archived bytes')
    path = ARCHIVE / 'manifest.json'
    digest = save(path, {'status': 'PRESERVED_BOTH_LANDINGS_AND_PRIOR_L_DELIVERY_BEFORE_NAVIGATION',
                         'archived_exact_bytes': copied, 'affected_existing_paths': list(TARGETS),
                         'child_readme_exists': CHILD_README in READMES,
                         'old_L_final_delivery_sha256': PINS[L_DELIVERY],
                         'all_335_prior_delivery_bindings_verified_before_edit': True,
                         'L_delivery_navigation_bindings_archived_by_path': snapshots,
                         'all_existing_landing_and_nonlanding_html_sha256': pages,
                         'PLAN_sha256': PINS[PLAN], 'nonlinear_input_pins': PINS})
    print(json.dumps({'archive_manifest_sha256': digest, 'archived_files': len(copied),
                      'L_bindings_rehashed': len(old), 'pages_rehashed': len(pages)}, indent=2))


class Anchors(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = Counter()

    def handle_starttag(self, tag, attrs):
        attr = dict(attrs)
        if tag == 'a' and 'href' in attr:
            self.hrefs[attr['href']] += 1


def hrefs(path, kind):
    text = path.read_text(encoding='utf-8')
    if kind == '.md':
        return Counter(re.findall(r'\]\(([^)\n]+)\)', text))
    parser = Anchors()
    parser.feed(text)
    return parser.hrefs


def local_from(page, href):
    uri = urlsplit(href)
    require(not uri.scheme and not uri.netloc and uri.path and not uri.query,
            f'new navigation link must be a local path: {href}')
    to = (page.parent / uri.path).resolve()
    require(to.is_relative_to(ROOT) and to.is_file(), f'broken or out-of-dossier link: {href}')
    return to.relative_to(ROOT).as_posix()


def named_functions(path, names):
    text = path.read_text(encoding='utf-8')
    return {node.name: ast.get_source_segment(text, node) for node in ast.parse(text).body
            if isinstance(node, ast.FunctionDef) and node.name in names}


def landing_section(payload):
    marker = b'<h2>Nonlinear transfer depends on the observer</h2>'
    require(payload.count(marker) == 1, 'one nonlinear navigation section required')
    start = payload.index(marker)
    end = payload.index(b'<h2>', start + len(marker))
    return payload[start:end], payload[:start] + payload[end:], payload[end:payload.index(b'</h2>', end) + 5]


def isolate_navigation():
    # The full child landing generator can also refresh its evidence inventory.
    # Keep that independent inventory at its prior publication snapshot here.
    path = ROOT / LANDINGS[1]
    baseline = ARCHIVE / LANDINGS[1]
    generated = path.read_bytes()
    section, _, following = landing_section(generated)
    original = baseline.read_bytes()
    require(original.count(following) == 1, 'landing section insertion anchor drift')
    extras = {h for h in hrefs(path, '.html') - hrefs(baseline, '.html')
              if 'transfer_ambiguity_nonlinear_probe_v1/' not in h}
    require(extras == {'../evidence/L_signed_root_final_delivery.json',
                       '../evidence/L_signed_root_space_update.json'},
            'unexpected landing changes in initial generator pass')
    first_pass = HERE / 'initial_landing_regeneration'
    first_pass.mkdir(exist_ok=False)
    snapshot = first_pass / 'controlled_furrowing_index.html'
    with snapshot.open('xb') as handle:
        handle.write(generated)
    updated = original.replace(following, section + following, 1)
    temp = path.with_suffix('.navigation_tmp')
    with temp.open('xb') as handle:
        handle.write(updated)
    temp.replace(path)
    claim = save(first_pass / 'receipt.json', {
        'status': 'ISOLATED_INTENDED_NAVIGATION_FROM_DYNAMIC_EVIDENCE_INVENTORY',
        'initial_generated_child_landing_sha256': sha(snapshot),
        'archived_prior_child_landing_sha256': sha(baseline),
        'current_child_landing_sha256': sha(path), 'generator_sha256': sha(ROOT / BUILDERS[1]),
        'incidental_evidence_inventory_links': sorted(extras),
        'scope': 'The linked report section is generated. All other landing bytes follow the original publication.'})
    print(json.dumps({'isolated_landing_sha256': sha(path), 'first_pass_receipt_sha256': claim}))


def verify():
    scientific_authority()
    archive_manifest = ARCHIVE / 'manifest.json'
    prior = read(archive_manifest)
    require(prior['affected_existing_paths'] == list(TARGETS), 'navigation paths drifted')
    copied = prior['archived_exact_bytes']
    for target, item in copied.items():
        require(sha(ROOT / item['archive_path']) == item['sha256'], f'archived bytes drifted: {target}')
    refreshed = {name: sha(ROOT / name) for name in TARGETS}
    require(all(refreshed[name] != copied[name]['sha256'] for name in TARGETS),
            'every intended navigation target should have the new bytes')
    l_bindings = read(ROOT / L_DELIVERY)['source_and_artifact_hashes']
    superseded = prior['L_delivery_navigation_bindings_archived_by_path']
    require(superseded == {k: h for k, h in l_bindings.items() if k in TARGETS},
            'L delivery transition must explicitly enumerate all superseded navigation paths')
    retained = {k: h for k, h in l_bindings.items() if k not in superseded}
    check_bytes(retained)
    require(len(retained) + len(superseded) == 335, 'all L delivery bindings required')
    before = prior['all_existing_landing_and_nonlanding_html_sha256']
    after = page_snapshot()
    require(set(before) == set(after) and {k for k in before if before[k] != after[k]} == set(LANDINGS),
            'only the two affected landing pages may be re-rendered')
    readback = {}
    for name in (*READMES, *LANDINGS):
        origin, frozen = ROOT / name, ROOT / copied[name]['archive_path']
        old, new = hrefs(frozen, origin.suffix), hrefs(origin, origin.suffix)
        increments = new - old
        require(not (old - new), f'original navigation link removed in {name}')
        added = {ref: local_from(origin, ref) for ref in increments}
        require(sorted(added.values()) == sorted(LINKS) and all(v == 1 for v in increments.values()),
                f'five scoped links per navigation path expected: {name}')
        if name in LANDINGS:
            _, remaining, _ = landing_section(origin.read_bytes())
            require(remaining == frozen.read_bytes(), f'landing mutated outside the new section: {name}')
        readback[name] = added
    for builder in BUILDERS:
        old = named_functions(ROOT / copied[builder]['archive_path'], {'convert', 'shell'})
        new = named_functions(ROOT / builder, {'convert', 'shell'})
        require(set(old) == {'convert', 'shell'} and len(new) == 2 and old == new,
                f'proof-page or standalone report rendering function changed: {builder}')
    result = {
        'status': 'PASS_LOCAL_NONLINEAR_NAVIGATION_STATIC_READBACK',
        'source_file': Path(__file__).relative_to(ROOT).as_posix(), 'source_sha256': sha(Path(__file__)),
        'archive_manifest': archive_manifest.relative_to(ROOT).as_posix(),
        'archive_manifest_sha256': sha(archive_manifest), 'L_delivery_sha256': PINS[L_DELIVERY],
        'original_L_delivery_paths_checked': 335,
        'L_delivery_bindings_that_still_match_current_bytes': len(retained),
        'superseded_shared_navigation_binding_archives': {k: copied[k] for k in superseded},
        'before_after_navigation_files': {k: {'before': copied[k]['sha256'], 'after': refreshed[k]} for k in TARGETS},
        'existing_child_html_pages_verified': sum(k.startswith(f'{CHILD}/docs/') for k in after),
        'all_existing_html_including_parent_redirect': len(after), 'rerendered_landing_pages': list(LANDINGS),
        'unchanged_chapter_and_redirect_count': len(after) - len(LANDINGS),
        'prior_nonlinear_report_converter_bytes_and_browser_landings_archived': {
            k: copied[k] for k in (BUILDERS[1], *LANDINGS, NONLINEAR_BROWSER)},
        'converter_and_shell_function_sources_match_previous_versions': True,
        'both_landings_match_prior_bytes_outside_new_navigation_section': True,
        'new_verified_local_links_by_navigation_document': readback,
        'canonical_PLAN_sha256': PINS[PLAN], 'current_science_and_report_receipts_sha256': PINS,
        'scope': 'Historical publication receipts retain their original navigation snapshots. '
                 'Current proof files and scientific bindings are verified byte-for-byte. '
                 'The new browser delta should exercise both current landings.',
        'current_navigation_browser_readback_pending': True, 'remote_publication_or_contact': False,
    }
    path = HERE / 'receipt.json'
    digest = save(path, result)
    print(json.dumps({'status': result['status'], 'receipt_sha256': digest,
                      'L_delivery_current': len(retained), 'L_navigation_archived': len(superseded),
                      'local_links_checked': len(LINKS) * len(readback)}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('archive', 'isolate', 'verify'))
    args = parser.parse_args()
    {'archive': archive, 'isolate': isolate_navigation, 'verify': verify}[args.mode]()

