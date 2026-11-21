"""Apply old independent N18/640 intervals to the new adaptive144 outputs."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
H=Path(__file__).resolve().parent;C=H.parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def bounds(p):return F(p['lower_exact']),F(p['upper_exact'])
def overlap(p,q):a,b=bounds(p);c,d=bounds(q);return a<=d and c<=b
def main():
    prior=C/'evidence/exact_bulk_pilot_independent_review.json';new=C/'runs/exact_integration_bridge_v1/bulk_assembly_B/results.json'
    if sha(prior)!='fe497971bb0c140e292194ebf1735cf2a1fb490ed8e973e75c392432ea13624b'or sha(new)!='c7938f431321d5b189efcd057fbf6635fbf0ea0d583cf0431f633ceeb7c23dfd':raise ValueError('authoritative receipt differs')
    old=read(prior);root=read(new);oldhome=C/'runs/exact_integration_bridge_v1/bulk_pilot_independent_replay';p=read(oldhome/'protocol.json');checks={};rows=[]
    checks['prior_independent_review_PASS']=old['status']=='PASS_INDEPENDENT_POINT_BULK_PILOT'and all(v is True for v in old['checks'].values())
    checks['prior_independent_protocol_matches']=sha(oldhome/'protocol.json')==old['protocol_sha256']
    checks['all_prior_scientific_and_review_inputs_still_match']=all(sha(C/k)==d for k,d in p['input_hashes'].items())
    mapping={x['element']:x for x in root['elements']}
    for item in old['elements']:
        file=oldhome/item['file'];x=read(file);index=x['element'];fresh=C/mapping[index]['file'];y=read(fresh)
        if sha(file)!=item['sha256']or sha(fresh)!=mapping[index]['sha256']:raise ValueError('element changed')
        a={(v['i'],v['j']):v['interval']for v in x['independent_Hessian_upper']};b={(v['i'],v['j']):v['enclosure']for v in y['H_integral_minus_cache_upper_triangle']}
        good=x['global_indices']==y['global_indices']and len(x['independent_force'])==len(y['force_integral_minus_cache'])and all(overlap(u,v)for u,v in zip(x['independent_force'],y['force_integral_minus_cache']))and set(a)==set(b)and all(overlap(a[k],v)for k,v in b.items())
        checks[f'prior_independent_element{index}_overlaps_new_adaptive_output']=good
        rows.append({'element':index,'force_entries':len(x['independent_force']),'H_entries':len(a),'old_independent_sha256':sha(file),'new_adaptive_sha256':sha(fresh),'PASS':good})
    checks['exactly_original_four_reused']={x['element']for x in rows}=={0,6,72,66}
    out={'status':'PASS_REUSE_PRIOR_FOUR_FOR_ADAPTIVE_ASSEMBLY'if all(checks.values())else'FAIL','checks':checks,'rows':rows,'source_sha256':sha(Path(__file__)),'prior_receipt_sha256':sha(prior),'new_assembly_sha256':sha(new),'scope':'Exact endpoint-overlap and immutable-context comparison only. No new integration. Together with two new48/143 replays, six of144local elements have independent numerical replay coverage.'}
    (H/'prior_four_reuse.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'checks':len(checks),'failed':[k for k,v in checks.items()if not v],'sha256':sha(H/'prior_four_reuse.json')}))
    if not all(checks.values()):raise SystemExit(1)
if __name__=='__main__':main()
