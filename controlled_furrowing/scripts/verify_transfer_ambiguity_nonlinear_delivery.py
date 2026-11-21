"""Bind the delivered numerical component, historical navigation and open gates."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SITE=ROOT.parent
BASE=ROOT/'runs/transfer_ambiguity_nonlinear_probe_v1'
PIN={
 'actual/results.json':'2ddc2fb3f0da426551bf74b6e2aac7bb489dfe848f8d9e73cfc42fddecbe99d6',
 'independent_actual_preparation/actual_replay/receipt.json':'c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd',
 'observer_robustness/actual/results.json':'53a17296ddae08b4d1fe6b7e9b20ad0a377bf8dc962142bbfb9eb9877a381647',
 'observer_robustness/independent_actual/receipt.json':'703b574bc51a083b10c7328d9622862e939ddabe791b18ce9e1817f07718f69f',
 'certification_design/budgets_v1/actual/results.json':'689948766db950ee4ddf399e6731aea08890cc7a1f84c5615d19e7ded68b0872',
 'certification_design/budgets_v1/independent_actual/receipt.json':'978a2a4095bb8735103a5ee2c11ad266a6b0e052e22af2b7a36ef29f1f38265c',
 'certification_design/cached_pilots_v1/independent_actual_preparation/actual_replay/receipt.json':'4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17',
 'report/report_build_receipt.json':'79cae6792ea8fbfef89bd9d69a3c87f14f419fb1a6c03384f781caf92eeff7c9',
 'report/scientific_readback/receipt.json':'42a40beabd178385bc14afa5d15106f43d390883427413c90474cd34eafc33e3',
 'report/browser_readback_v2/receipt.json':'c0236dd4a464cd29b2b4925a3232fb4017ca0ddc65952cefb839118005a52d8d',
 'nonlinear_local_navigation_v1/receipt.json':'e336d501aa2c37061616b11d8bef175b4b6c968bce13650f9cd9fd46e1b157ab',
 'nonlinear_local_navigation_v1/browser_readback_v2/receipt.json':'a80f2f22cc3d1037e72f34e3fc7fab95e518cbc9d0dea90193c4d4f2dd3dc61d',
 'originality_followup/receipt.json':'6797ad5d05f71535f13edc519276e57bc5176d7d5b0e6273b1aa91eba48653d0',
}


def sha(p):
    return sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    return json.loads(Path(p).read_bytes())


def main():
    current={}; historical={}
    nav=read(BASE/'nonlinear_local_navigation_v1/receipt.json')
    archive=nav['superseded_shared_navigation_binding_archives']
    def bound(p,h):
        p=Path(p).resolve()
        assert p.is_relative_to(SITE.resolve())
        rel=p.relative_to(SITE).as_posix()
        if sha(p)==h:
            current[rel]=h
        else:
            old=archive[rel]
            assert old['sha256']==h and sha(SITE/old['archive_path'])==h
            historical[rel]={'previous_sha256':h,'archive_path':old['archive_path'],
                             'current_sha256':sha(p),'reason':'Preserved navigation update; scientific converter functions unchanged.'}
            current[old['archive_path']]=h
    for name,h in PIN.items():
        bound(BASE/name,h)
    primary=read(BASE/'actual/results.json')
    audit=read(BASE/'independent_actual_preparation/actual_replay/receipt.json')
    assert primary['both_endpoints_satisfy_declared_allowance'] is True
    assert audit['status']=='PASS_INDEPENDENT_ACTUAL_NONLINEAR_STORED_FIELD_REPLAY'
    assert audit['checkpoints']==44 and audit['independently_recomputed_small_residuals']==136
    assert audit['source_compatible_witness'] is False and audit['exact_equilibrium_position_bounded'] is False
    gap=list(map(Q,audit['stored_large_Q_signed_difference_um_exact']))
    assert Q('0.15475515')<gap[0]<=gap[1]<Q('0.15475517')
    for name,h in audit['all_hashes'].items():
        bound(ROOT/name,h)
    secondary=read(BASE/'observer_robustness/actual/results.json')
    assert secondary['small_constraint_evaluations']==544 and secondary['primary_replay_exactly_matches'] is True
    for name,v in secondary['common_pair_requirement_by_variant'].items():
        assert v['status_at_frozen_allowance']==('ALLOWANCE_SATISFIED' if name=='current_disk__spatial_max' else 'ALLOWANCE_EXCEEDED')
    budgets=read(BASE/'certification_design/budgets_v1/actual/results.json')
    assert len(budgets['field_budgets'])==36 and len(budgets['four_pilot_selection'])==4
    assert budgets['exact_root_error_bounds_proved']==0
    cache=read(BASE/'certification_design/cached_pilots_v1/independent_actual_preparation/actual_replay/receipt.json')
    assert cache['root_certificates_cached_quadrature_only']==4
    assert cache['root_certificates_exact_spatial_integral']==0
    assert cache['stability_or_whole_body_embedding_certified'] is False and cache['source_compatibility_or_PDE_certified'] is False
    assert all(v is True for v in cache['checks'].values())
    for name,h in cache['producer_and_payload_sha256'].items():
        bound(ROOT/name,h)
    build=read(BASE/'report/report_build_receipt.json')
    for name,h in {**build['input_hashes'],**build['report_and_figure_hashes']}.items():
        bound(ROOT/name,h)
    science=read(BASE/'report/scientific_readback/receipt.json')
    assert science['status'].startswith('PASS')
    for name,h in science['rehashed_files'].items():
        bound(ROOT/name,h)
    browser=read(BASE/'report/browser_readback_v2/receipt.json')
    assert browser['status']=='PASS_BROWSER_READBACK' and len(browser['checks'])==2
    for row in browser['checks']:
        assert row['status']=='PASS' and row['math']==row['rendered']==34
        assert len(row['images'])==3 and all(x['status']==200 for x in row['local_links'])
        for image in row['screenshots']:
            bound(BASE/'report/browser_readback_v2'/image['file'],image['sha256'])
    nav_browser=read(BASE/'nonlinear_local_navigation_v1/browser_readback_v2/receipt.json')
    assert nav_browser['status']=='PASS_NEW_LOCAL_NAVIGATION_BROWSER' and len(nav_browser['checks'])==4
    for name,h in nav_browser['protected_hashes'].items():
        bound(SITE/name,h)
    for name,v in nav['before_after_navigation_files'].items():
        bound(SITE/name,v['after'])
    l=read(ROOT/'evidence/L_signed_root_final_delivery.json')
    for name,h in l['source_and_artifact_hashes'].items():
        bound(SITE/name,h)
    assert len(historical)==4 and nav['L_delivery_bindings_that_still_match_current_bytes']==331
    bound(ROOT/'PLAN.md','993240af5d30c8e5b72ab31a1278a524314d38ede48486b677b021d84bcde500')
    space=read(BASE/'space_update_receipt.json')
    assert space['sequence']==25 and space['blocks']==168 and all(space['checks'].values())
    for image in space['images']:
        bound(BASE/'report'/image['name'],image['sha256'])
    for name in ('space_update_receipt.json','field_manifest.json','terminal_admission.json','FIVE-GATE-COMPLETION-AUDIT.md'):
        bound(BASE/name,sha(BASE/name))
    out={'status':'PASS_NONLINEAR_NUMERICAL_COMPONENT_AND_LOCAL_DOCUMENTATION',
         'verifier_source_sha256':sha(__file__),'current_hashes':current,'historical_navigation_bindings':historical,
         'primary_checkpoints_verified':44,'small_primary_residuals':136,'secondary_residuals':544,
         'stored_large_gap_um_exact':list(map(str,gap)),'cached_quadrature_roots_independently_certified':4,
         'new_exact_spatial_integral_roots_certified':0,'full_compatible_set_computed':False,
         'source_compatibility_established':False,'novelty_established':False,'full_five_gate_goal_complete':False,
         'scope':'Delivery binding and requirement-scope audit. Relies on retained independent mathematical, browser and Space readbacks; no new proof engine, browser session or remote read here.'}
    p=BASE/'delivery_receipt.json'
    with p.open('x',encoding='utf-8') as f:
        json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps({'status':out['status'],'current_hashes':len(current),'historical_bindings':len(historical),'delivery_sha256':sha(p)}))


if __name__=='__main__':
    main()
