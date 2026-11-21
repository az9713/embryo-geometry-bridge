"""Independent exact scalar and premise readback; no applicant functions called."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]


def main():
    def sha(p):return sha256(p.read_bytes()).hexdigest()
    result=HERE/'results.json'
    assert sha(result)=='005f1f100da3cd558212d6eceef7a1ca886602e94c64438f08dd009f69e48b61'
    data=json.loads(result.read_text())
    inputs=data['source_and_input_sha256']
    assert len(inputs)==57 and all(sha(ROOT/p)==h for p,h in inputs.items())
    assert sha(HERE/'PROOF.md')==data['proof_sha256']=='edf6c75d294a70052cff22302785766e14cea31a036f22958e06169c3a5f4b4f'
    def load(h):
        paths=[ROOT/p for p,v in inputs.items() if v==h and p.endswith('.json')]
        assert len(paths)==1,(h,len(paths))
        return json.loads(paths[0].read_text())
    for h in ('4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17',
              'f766b8b41142ce1f82f3e672f423c32485ae436487d9418c42d907b6f27e8c74',
              '769a65a128aaa79018c050c7bd61d38999166621742946cd0f6a13876854ba83'):
        parent=load(h)
        assert parent['status'].startswith('PASS_INDEPENDENT')
    rows=[]
    assert len(data['rows'])==data['candidate_count']==4
    for row in data['rows']:
        ids=row['target_identity']; hashes=row['premise_sha256']; pid=ids['pilot_id']
        point=load(hashes['cached_point_sha256']); spd=load(hashes['native_SPD_case_sha256'])
        large=load(hashes['native_large_cached_ZQ_case_sha256']); producer=load(hashes['producer_large_cached_ZQ_case_sha256'])
        smallbox=load(hashes['small_box_sha256']); largebox=load(hashes['large_box_sha256'])
        sg=load(hashes['small_geometry_sha256']); lg=load(hashes['large_geometry_sha256'])
        assert point['status']=='PASS_INDEPENDENT_CACHED_LOCAL_ROOT'
        assert spd['status']=='PASS_INDEPENDENT_EXACT_SYMMETRIC_PART_SPD' and spd['target_identity']==ids
        assert all(x['pilot_id']==pid for x in (point,spd,large,sg,lg))
        assert point['field_sha256']==ids['field_sha256']==sg['proof_center_field_sha256']==lg['proof_center_field_sha256']
        assert point['cache_packet_sha256']==spd['cache_packet_sha256']==large['cache_packet_sha256']==ids['cache_packet_sha256']
        assert spd['cached_point_proof_sha256']==hashes['cached_point_sha256']
        assert spd['derivative_case_sha256']==hashes['native_large_cached_ZQ_case_sha256']
        assert large['producer_result_sha256']==hashes['producer_large_cached_ZQ_case_sha256']
        assert Q(point['exact_strength'])==Q(large['proof_drive_exact'])==Q(ids['drive_exact'])==Q(ids['gain_exact'])*Q(ids['dose_exact'])
        assert Q(ids['energy_scale_binary_exact'])>0
        rs,rl=Q(1,2**36),Q(1,2**34)
        assert Q(point['radius_exact'])==Q(smallbox['corrector_radius_exact'])==Q(sg['radius_exact'])==rs
        assert Q(large['radius_exact'])==Q(largebox['corrector_radius_exact'])==Q(lg['radius_exact'])==rl
        assert spd['larger_radius_box_sha256']==large['box_sha256']==lg['box_sha256']==hashes['large_box_sha256']
        assert sg['box_sha256']==hashes['small_box_sha256']
        for key in ('reference_nodes_exact','current_centers_exact','triangles','free'):
            assert smallbox[key]==largebox[key]
        assert len(smallbox['current_centers_exact'])==325 and len(smallbox['free'])==611
        assert all(Q(ll)<=Q(sl)<=Q(sh)<=Q(lh) for s,l in zip(smallbox['coordinate_hulls_exact'],largebox['coordinate_hulls_exact'],strict=True)
                   for (sl,sh),(ll,lh) in zip(s,l,strict=True))
        norm=Q(spd['exact_R_infinity_norm'])
        assert norm==Q(point['preconditioner_infinity_norm_exact'])==Q(large['preconditioner_infinity_norm_exact'])>0
        assert spd['exact_target_symmetric_part_spd'] is True and 0<=Q(spd['max_reviewed_q_upper_exact'])<1
        z=max(Q(large['reviewed_ZQ_upper_exact']),Q(producer['ZQ_upper_exact']))
        assert 0<=z<1 and z==Q(row['cached_ZQ_on_coercive_cube_upper_exact'])
        zsmall=Q(point['independent_z_upper_exact']); eta=Q(point['independent_eta_upper_exact'])
        assert 0<=zsmall<1 and eta+zsmall*rs<rs
        error=eta/(1-zsmall)
        assert error==Q(point['independent_posterior_error_upper_exact'])==Q(row['root_distance_upper_exact'])<rs<rl
        assert rs-error==Q(row['root_interior_margin_small_cube_exact'])>0
        assert rl-error==Q(row['root_interior_margin_coercive_cube_exact'])>0
        domain=large['whole_box_domains']
        assert large['certified_site_counts']=={'bulk':9216,'surface':288}
        assert all(Q(v['lower_exact'])>0 for v in domain['strict_lower_margins'].values())
        assert Q(domain['strict_lower_margins']['bulk_J_min']['lower_exact'])>Q(1,10**6)
        m=(1-z)/norm
        assert m==Q(row['normalized_cached_H_eigenlower_exact'])>0
        assert m/2==Q(row['normalized_cached_energy_quadratic_gap_coefficient_exact'])>0
        rows.append({'pilot_id':pid,'normalized_cached_H_eigenlower_exact':str(m),
                     'energy_gap_coefficient_exact':str(m/2),'root_large_cube_interior_margin_exact':str(rl-error)})
    out={'status':'PASS_INDEPENDENT_FOUR_CACHED_MINIMUM_APPLICATION','checks':{
        'all_57_current_premise_hashes_and_parent_reviews_checked':True,
        'same_field_exact_gain_dose_cache_R_and_energy_normalization':True,
        'all_four_nested_closed_coordinate_cubes_verified':True,
        'all_four_exact_curvature_and_strict_root_interior_bounds_reproduced':True,
        'classical_Lyapunov_and_Taylor_arguments_reviewed_by_root':True},
        'applicant_results_sha256':sha(result),'review_source_sha256':sha(Path(__file__)),
        'reused_mathematical_premises':'Existing independently reviewed C2 cached potential/site derivative, geometry, root and exact Sym(R) certificates. No new primitive proof or matrix/AD evaluation.',
        'rows':rows,'cached_local_minima_certified':4,'exact_spatial_minima_certified':0,
        'stable_loading_or_continuum_or_biological_fit_claimed':False,'full_five_gate_goal_complete':False}
    f=HERE/'independent_readback.json'
    with f.open('x',encoding='utf-8') as stream:json.dump(out,stream,indent=2);stream.write('\n')
    print(json.dumps({'review_sha256':sha(f),'cached_minima':4,'curvature_display':[float(Q(r['normalized_cached_H_eigenlower_exact'])) for r in rows]}))


if __name__=='__main__':main()
