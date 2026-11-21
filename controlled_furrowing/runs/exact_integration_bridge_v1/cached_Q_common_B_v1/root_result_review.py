"""Independent exact endpoint reduction and domain/output binding audit.

Reuses previously reviewed dense constitutive/variation kernels; this is NOT
an independent 9504-site Hessian reintegration or new root certificate.
"""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def read(p):return json.loads(p.read_text(encoding='utf8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
checks={}
def ck(k,v):
    checks[k]=bool(v)
    assert v,k
p=read(HERE/'protocol.json');d=read(HERE/'results.json');terminal=read(HERE/'terminal.json')
ck('expected_results_hash',sha(HERE/'results.json')=='d96ef672ea40ce807acdecd454ed7ec7a8f342b044dd0c802f4b9e57e57b4613')
ck('true_terminal_bound',terminal['status']=='FINISHED'and terminal['results_sha256']==sha(HERE/'results.json'))
ck('same_launch_review',d['prelaunch_review_sha256']==sha(HERE/'prelaunch_independent_review.json')=='ae413c6e11f28d2489bd382a18223267adefa75f5384215f7d2a95f7d6169410')
ck('current_source_protocol',d['source_sha256']==sha(ROOT/'scripts/exact_cached_contraction_B_common.py')and d['protocol_sha256']==sha(HERE/'protocol.json'))
for name,h in p['source_and_input_hashes'].items():ck('binding_'+name,sha(ROOT/name)==h)
row=read(HERE/'fresh_row_bounds.json');z0=Q(d['Z0_upper_exact']);v=list(map(Q,d['variation_row_upper_exact']));uniform=list(map(Q,d['uniform_row_upper_exact']))
ck('611_rows',len(v)==len(uniform)==611)
ck('exact_independent_endpoint_row_sums',uniform==[z0+x for x in v])
ck('exact_independent_maximum',Q(d['ZQ_upper_exact'])==max(uniform)==z0+max(v))
ck('all_nonnegative',z0>=0 and min(v)>=0)
ck('row_artifact_matches',row['variation_row_upper_exact']==d['variation_row_upper_exact']and row['uniform_row_upper_exact']==d['uniform_row_upper_exact'])
ck('strict_contraction_outward_readable_bound',Q(d['ZQ_upper_exact'])<=Q(105040426,10**9)<1)
ck('same_common_radius',d['common_radius_exact']==p['common_radius_exact']==row['common_radius_exact'])
ck('same_C_dose_and_drive',d['parameter_cell']==p['parameter_cell']and d['dose_exact']==p['dose_exact']and d['base_strength_exact']==p['base_strength_exact'])
domfile=HERE/'cached_domain_bounds.json';domains=read(domfile)
ck('domain_artifact_hash',d['cached_domain_bounds_sha256']==sha(domfile))
bulk=[r for r in domains['rows']if r['kind']=='bulk'];surf=[r for r in domains['rows']if r['kind']=='surface']
ck('all9504_site_domains',len(bulk)==9216 and len(surf)==288 and len(domains['rows'])==9504)
ck('whole_bulk_convex_boxes_strict_positive',all(Q(r['det_lower_exact'])>0 and Q(r['hoop_lower_exact'])>0 and Q(r['J_lower_exact'])>Q(1,10**6)for r in bulk))
ck('whole_surface_convex_boxes_strict_positive',all(Q(r['cached_radial_lower_exact'])>0 and Q(r['cached_tangent_norm2_lower_exact'])>0 for r in surf))
ck('independent_bulk_domain_minimum',min(Q(r['J_lower_exact'])for r in bulk)==Q(domains['minimum_bulk_J_lower_exact']))
ck('independent_surface_domain_minima',min(Q(r['cached_radial_lower_exact'])for r in surf)==Q(domains['minimum_surface_radial_lower_exact'])and min(Q(r['cached_tangent_norm2_lower_exact'])for r in surf)==Q(domains['minimum_surface_tangent_norm2_lower_exact']))
ck('root_integral_PDE_claims_excluded',d['root_selfmap_or_integral_transfer_certified']is False and d['integral_H_error_bound_computed']is False and d['full_17_or_Theta_or_PDE_computed']is False)
out={'status':'PASS_EXACT_REDUCTION_AND_DOMAIN_AUDIT','checks':checks,'ZQ_upper_exact':d['ZQ_upper_exact'],
     'input_hashes':{p.relative_to(ROOT).as_posix():sha(p)for p in [HERE/'protocol.json',HERE/'prelaunch_independent_review.json',HERE/'results.json',HERE/'fresh_row_bounds.json',domfile,HERE/'terminal.json',Path(__file__)]},
     'scope':'Independent exact endpoint/max/domain and binding audit of the new cached-Q common-tube derivative result. Previously reviewed dense analytic kernels reused; not independent full-site assembly or integral/root/PDE verification.'}
dest=ROOT/'evidence/exact_cached_Q_common_B_reduction_review.json'
assert not dest.exists()
dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps({'status':out['status'],'checks':len(checks),'ZQ':float(Q(d['ZQ_upper_exact'])),'review_sha256':sha(dest)}))
