"""Independent source-to-CSV and artifact readback for the wave03 report."""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import csv
import json

HERE = Path(__file__).resolve().parent
WAVE = HERE.parent
ROOT = next(p for p in HERE.parents if (p/'scripts/build_report.py').is_file())
SITE = ROOT.parent

def sha(p): return sha256(p.read_bytes()).hexdigest()
def load(p): return json.loads(p.read_bytes())
def rows(name):
    with (HERE/name).open(encoding='utf8',newline='') as f:
        return list(csv.DictReader(f))

def main():
    build = load(HERE/'build_receipt.json')
    assert sha(HERE/'build_receipt.json') == 'd5c0418863062609955e0735852ecb889e559d557bbf98d04164249d6d4e6370'
    assert sha(HERE/'build.py') == build['source_sha256']
    for p,h in build['input_sha256'].items(): assert sha(SITE/p) == h
    for p,h in build['file_sha256'].items(): assert sha(HERE/p) == h
    # The plot terms and each root error transfer are rebuilt from native output.
    rootrows = rows('exact_bounds.csv')
    assert len(rootrows) == 4
    for r in rootrows:
        proof = WAVE/'global_reduction_v1/independent'/r['case_id']/'results.json'
        d = load(proof)
        a,v,t = d['arithmetic'], d['observer_and_optional_minimum'], d['target']
        radius, eta, z = Q(d['radius_exact']), Q(a['eta_upper_exact']), Q(a['z_upper_exact'])
        expected = {
            'radius_exact':radius, 'eta_upper_exact':eta, 'eta_over_radius_exact':eta/radius,
            'z_upper_exact':z, 'selfmap_headroom_fraction_exact':1-z-eta/radius,
            'normalized_curvature_lower_exact':Q(v['normalized_H_euclidean_lower_exact']),
            'posterior_exact':Q(a['posterior_error_upper_exact']),
            'archive_shift_exact':Q(t['archived_center_shift_exact']),
            'binary_addition_gap_exact':Q(t['original_addition_gap_exact']),
            'physical_nodal_error_um_exact':Q(v['physical_nodal_error_um_upper_exact'])}
        assert all(Q(r[k]) == x for k,x in expected.items()) and r['proof_sha256'] == sha(proof)
        assert eta+z*radius < radius and z < 1 and v['strict_minimum_proved'] is True
    obs = load(WAVE/'observer_application_v1/independent_actual_readback/receipt.json')
    indexed = {(c['case_id'],x['variant'],x['observable']):x
               for c in obs['cases'] for x in c['comparisons']}
    obsrows = rows('exact_observer_comparisons.csv')
    assert len(obsrows) == len(indexed) == 64
    for r in obsrows:
        x = indexed[(r['case_id'],r['variant'],r['observable'])]
        assert list(map(Q,(r['prediction_lo_um_exact'],r['prediction_hi_um_exact']))) == list(map(Q,x['exact_spatial_root_observer_interval_um_exact']))
        assert list(map(Q,(r['residual_lo_um_exact'],r['residual_hi_um_exact']))) == list(map(Q,x['absolute_residual_um_exact']))
        assert Q(r['observed_um_exact']) == Q(x['observed_um_exact']) and r['standing'] == x['standing']
        assert Q(r['allowance_um_exact']) == Q(obs['allowance_um_exact'])
    stress = load(ROOT/'runs/nonaxisymmetric_surface_sign_preparation_v1/independent_actual_readback.json')
    indexed = {(c['case_id'],x['face'],x['edge'],tuple(map(Q,x['window_t_exact']))):x
               for c in stress['case_readback'] for x in c['all_windows']}
    stressrows = rows('transverse_windows.csv')
    assert len(stressrows) == len(indexed) == 672
    for r in stressrows:
        x = indexed[(r['case_id'],r['face'],int(r['edge_index']),tuple(map(Q,(r['t_lo_exact'],r['t_hi_exact']))))]
        assert r['standing'] == x['independent_standing']
        assert list(map(Q,(r['reference_R_lo_exact'],r['reference_R_hi_exact']))) == list(map(Q,x['reference_R_exact']))
        assert Q(r['r_R_lo_exact']) == Q(x['independent_r_R_exact'][0])
        if r['standing'] == 'STRICT_NEGATIVE_PATCH':
            assert Q(r['tau_strict_upper_exact']) == Q(x['independent_tau_negative_upper_exact']) < 0
            assert Q(r['minus_tau_r_R_squared_over_Js_lower_exact']) == Q(x['independent_minus_tau_r_R_squared_over_Js_lower_exact']) > 0
    slabrows = rows('exact_slab_geometry.csv')
    assert len(slabrows) == 7
    own = load(ROOT/'runs/exact_loading_tube_pilot_geometry_actual_readback_v1/root_segment_actual_readback/receipt.json')
    native = load(ROOT/'runs/exact_loading_tube_pilot_geometry_actual_readback_v1/segment_actual/receipt.json')
    mapping = {(v['slabs'],r['slot']):(r,n) for v,w in zip(own['levels'],native['levels'],strict=True)
               for r,n in zip(v['rows'],w['checks'],strict=True)}
    for r in slabrows:
        if int(r['slabs']) == 1:
            d = load(ROOT/'runs/exact_loading_tube_pilot_geometry_actual_readback_v1/actual/receipt.json')['body_summary']
            assert list(map(Q,(r['J_lo_exact'],r['J_hi_exact']))) == list(map(Q,d['body_J_bounds_exact']))
            assert r['body_and_log_guards_pass'] == 'False'
            continue
        a,n = mapping[(int(r['slabs']),int(r['slab']))]
        assert list(map(Q,(r['J_lo_exact'],r['J_hi_exact']))) == list(map(Q,n['J_enclosure_exact']))
        assert r['body_and_log_guards_pass'] == str(a['passes'])
        assert list(map(Q,(r['I_lo_exact'],r['I_hi_exact']))) == list(map(Q,n['exact_dose_subinterval']))
    browser = load(HERE/'browser_readback/receipt.json')
    assert sha(HERE/'browser_readback/receipt.json') == 'b4d8a9c225cb3250803808762806e57ba4dc0e0a09674d7a1d931ba666149446'
    assert browser['status'] == 'PASS_EDGE_WAVE03_DESKTOP_MOBILE' and browser['javascript_disabled'] is True
    receipt = {'status':'PASS_ROOT_WAVE03_REPORT_HASHES_AND_EXACT_CSVS',
               'source_sha256':sha(Path(__file__)), 'build_receipt_sha256':sha(HERE/'build_receipt.json'),
               'browser_receipt_sha256':sha(HERE/'browser_readback/receipt.json'),
               'browser_root_observed_exit_code':0,'browser_root_observed_chunk':'8f94e7',
               'build_root_observed_exit_code':0,'build_root_observed_chunk':'008c0f',
               'input_hashes_rechecked':len(build['input_sha256']), 'outputs_rechecked':len(build['file_sha256']),
               'root_csv_rows':4,'observer_csv_rows':64,'stress_csv_rows':672,'slab_csv_rows':7,
               'three_scientific_figures_visually_inspected':True,'full_five_gate_complete':False}
    with (HERE/'root_readback.json').open('x',encoding='utf8',newline='\n') as f:
        json.dump(receipt,f,indent=2); f.write('\n')
    print(json.dumps({'receipt_sha256':sha(HERE/'root_readback.json'),'status':receipt['status']}))

if __name__ == '__main__': main()
