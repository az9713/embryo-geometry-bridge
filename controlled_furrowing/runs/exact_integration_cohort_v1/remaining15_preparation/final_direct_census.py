"""Cross-reference all fifteen completed direct reviews; do not reclassify roots."""
from hashlib import sha256
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
STAGE=ROOT/'runs/exact_integration_cohort_v1/remaining15_reduction_preparation'


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_bytes())


def main():
    rows=[]
    for letter in 'CDEFGHIJKMNOPQR':
        name='optoGEF2_small_col_'+letter
        path=STAGE/'independent_actual'/f'{name}.json'
        doc=read(path)
        assert doc['record_id']==name and doc['status']=='PASS_DIRECT_REDUCTION_UNRESOLVED'
        assert doc['F_entries']==doc['RF_entries_containing_zero']==611 and doc['H_upper_entries']==6676
        assert doc['root_certified'] is False and doc['strict_local_minimum_certified'] is False
        assert doc['packet_route_counts']=={'bulk':144,'surface':24,'axis_vertex':6,'axis_edge_physical':6,'scaled':132}
        pins=doc['producer_pins']
        for filename,key in ((STAGE/name/'results.json','results_sha256'),
                             (STAGE/name/'terminal.json','terminal_sha256'),
                             (STAGE/f'{name}_input_admission.json','admission_sha256'),
                             (ROOT/pins['observation_file'],'observation_sha256')):
            assert sha(filename)==pins[key]
        terminal=read(STAGE/name/'terminal.json')
        assert terminal['status']=='FINISHED' and terminal['exit_code']==0 and terminal['results_sha256']==pins['results_sha256']
        rows.append({'record_id':name,'independent_review_file':path.relative_to(ROOT).as_posix(),
                     'independent_review_sha256':sha(path),'producer_pins':pins,
                     'direct_status':doc['status'],'same_closed_C_sha256':doc['same_closed_C_sha256']})
    assert len(rows)==15 and len({x['same_closed_C_sha256'] for x in rows})==1
    result={'status':'ALL_15_DIRECT_REVIEWS_BOUND_TESTS_UNRESOLVED','source_sha256':sha(Path(__file__)),
            'records':rows,'bulk_packets':2160,'surface_packets':360,'total_packets':2520,
            'signed_RF_rows_containing_zero':15*611,'new_roots_certified_by_direct_tests':0,
            'scope':'Receipt and terminal cross-reference only; earlier independent mathematical replays are the premises. '
                    'All-zero-containing direct intervals do not exclude roots. Later centered C exclusion and D unresolved work are separate stronger analyses.',
            'five_gate_objective_complete':False}
    dest=HERE/'FINAL-DIRECT-CENSUS.json'
    with dest.open('x',encoding='utf-8') as f:
        json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'status':result['status'],'records':len(rows),'census_sha256':sha(dest)}))


if __name__=='__main__':
    main()
