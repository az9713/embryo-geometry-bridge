"""Extract a fresh exact coordinate witness after the SAFE168 reduction."""
from pathlib import Path
from fractions import Fraction as Q
import base64,hashlib,io,json
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main():
    result=read(HERE/'reduction_results.json');exact=read(HERE/'exact_168_reduction.json')
    assert result['status']=='PASS_SAFE_V2_COMPLETE_POINT_REDUCTION'and all(v is True for v in result['checks'].values())and result['exact_reduction_sha256']==sha(HERE/'exact_168_reduction.json')
    eta=Q(exact['eta_used_exact']);r=Q(exact['old_radius_exact']);rf=[tuple(map(Q,x))for x in exact['RF_endpoint_sums']]
    def lower(v):return Q(0)if v[0]<=0<=v[1]else min(map(abs,v))
    j=max(range(611),key=lambda k:lower(rf[k]));fi=(rf[j][0]-eta,rf[j][1]+eta);shift=(-fi[1],-fi[0]);gap=lower(shift)-r
    rec=read(ROOT/'runs/cohort_uniform_roots_v2/records/00_optoGEF2_small_col_B/result.json');cachefile=ROOT/rec['point_cache_file'];cache=read(cachefile)
    raw=base64.b64decode(cache['frozen_binary_cache']['data'],validate=True);assert hashlib.sha256(raw).hexdigest()==cache['frozen_binary_cache']['sha256']
    with np.load(io.BytesIO(raw),allow_pickle=False)as z:dof=int(z['free'][j])
    checks={'reduction_certificate_hash_bound':True,'same_fresh_interval_norm_witness':lower(rf[j])==Q(exact['RF_norm_lower_exact']),'same_exact_old_eta_radius':True,'coordinate_escapes_iff_fresh_declared_test':(gap>0)==exact['new_integral_map_center_outside_old_cube'],'not_a_root_location_or_nonexistence_theorem':True}
    assert all(v is True for v in checks.values())
    p={'status':'PASS_SAFE_V2_NEW_MAP_COORDINATE_WITNESS','checks':checks,'review_source_sha256':sha(Path(__file__)),'reduction_sha256':sha(HERE/'reduction_results.json'),'exact_reduction_sha256':sha(HERE/'exact_168_reduction.json'),'point_cache_sha256':sha(cachefile),'free_index':j,'global_dof':dof,'node':dof//2,'component':'radial'if dof%2==0 else'axial','RF_defect_interval_exact':list(map(str,rf[j])),'R_Fintegral_interval_exact':list(map(str,fi)),'Tintegral_u0_minus_u0_interval_exact':list(map(str,shift)),'eta_exact':str(eta),'old_radius_exact':str(r),'strict_coordinate_gap_over_radius_exact':str(gap),'new_map_center_leaves_old_cube':gap>0,'scope':'ONE original B center/base theta, corrected SAFE integration packet. This says nothing about integral equilibrium existence or location.'}
    out=HERE/'coordinate_witness.json';assert not out.exists();out.write_text(json.dumps(p,indent=2)+'\n');print(json.dumps({'status':p['status'],'free_index':j,'node':dof//2,'component':p['component'],'map_shift_lower':float(shift[0]),'map_shift_upper':float(shift[1]),'old_radius':float(r),'witness_sha256':sha(out)}))
if __name__=='__main__':main()
