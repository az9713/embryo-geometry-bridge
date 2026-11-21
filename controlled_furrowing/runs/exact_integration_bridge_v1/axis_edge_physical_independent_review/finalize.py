from pathlib import Path
import hashlib,json,datetime
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def main():
    r=read(HERE/'review.json');checks={'all188_review_checks_true':len(r['checks'])==188 and all(v is True for v in r['checks'].values()),'all6_original_axis_edges':len(r['cases'])==6 and [x['element']for x in r['cases']]==list(range(72,78)),'all_point_derivative_and_uniform_coefficient_proofs':all(x['point_and_exact_free_derivative_equivalence']is True and x['uniform_coefficients_proved_by_signed_endpoint_extrema']is True for x in r['cases']),'all_input_hashes_current':all(sha(ROOT/k)==v for k,v in r['input_hashes'].items()),'owned_source_current':r['review_source_sha256']==sha(HERE/'review.py'),'evidence_matches_retained_bytes':(HERE/'review.json').read_bytes()==(ROOT/'evidence/exact_axis_edge_physical_independent_review.json').read_bytes()};assert all(v is True for v in checks.values()),checks
    p={'status':'PASS_AXIS_EDGE_REPRESENTATION_REVIEW_DELIVERY','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'review_sha256':sha(HERE/'review.json'),'report_sha256':sha(HERE/'REVIEW.md'),'finalizer_sha256':sha(__file__),'session_id':35474,'observed_terminal_exit_code':0,'actual_moment_integrations_or_cache_sums_run':False,'new_error_bound_or_root_certified':False};o=HERE/'delivery.json';assert not o.exists();o.write_text(json.dumps(p,indent=2)+'\n');print(json.dumps({'status':p['status'],'sha256':sha(o)}))
if __name__=='__main__':main()
