"""Extract one exact coordinate witness from the frozen full point reduction."""
from pathlib import Path
from fractions import Fraction as F
import base64,hashlib,io,json
import numpy as np
H=Path(__file__).resolve().parent;C=H.parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    receipt=C/'evidence/exact_complete_point_B_independent_review.json'
    if sha(receipt)!='bc22a23e7af96d0d1a99836401112378600c14be61d5cba957e6f15a12d699c5':raise ValueError('final review differs')
    result=json.loads(receipt.read_text());data=json.loads((H/'exact_reduction.json').read_text());rf=[tuple(map(F,x))for x in data['RF_endpoint_sums']]
    def lower(x):return F(0)if x[0]<=0<=x[1]else min(map(abs,x))
    index=max(range(611),key=lambda i:lower(rf[i]));lo,hi=rf[index];eta=F(data['eta_used_exact']);r=F(data['radius_exact']);step=(-hi-eta,-lo+eta)
    root=json.loads((C/'runs/cohort_uniform_roots_v2/records/00_optoGEF2_small_col_B/result.json').read_text());cache=json.loads((C/root['point_cache_file']).read_text());raw=base64.b64decode(cache['frozen_binary_cache']['data'])
    if hashlib.sha256(raw).hexdigest()!=cache['frozen_binary_cache']['sha256']:raise ValueError('cache corrupt')
    z=np.load(io.BytesIO(raw),allow_pickle=False);full=int(z['free'][index]);node,coordinate=divmod(full,2)
    out={'status':'CERTIFIED_EXPLICIT_NEW_MAP_CENTER_COORDINATE_OUTSIDE_OLD_CUBE','free_index_0based':index,'full_nodal_coordinate_index_0based':full,'node_0based':node,'coordinate':'radial'if coordinate==0 else'axial','reference_position_normalized':[str(F.from_float(float(x)))for x in z['reference'][node]],'R_deltaF_coordinate_interval_exact':[str(lo),str(hi)],'old_RFQ_absolute_component_bound_exact':str(eta),'new_T_I_u0_minus_u0_coordinate_interval_exact':list(map(str,step)),'old_cube_coordinate_interval_exact':[str(-r),str(r)],'outside_verified':step[0]>r or step[1]<-r,'distance_lower_over_old_radius_exact':str(lower(rf[index])/r-eta/r),'receipt_sha256':sha(receipt),'exact_reduction_sha256':sha(H/'exact_reduction.json'),'source_sha256':sha(Path(__file__)),'interpretation':'NEW exactintegralmap with oldR/cube sendsoldcenteroutside; oldQmap remainsvalid. This doesnotlocate orruleoutan integralroot.'}
    if not out['outside_verified']:raise ValueError('coordinate not outside')
    (H/'coordinate_witness.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'index':index,'node':node,'coordinate':out['coordinate'],'new_map_step_interval_display':list(map(float,step)),'radius':float(r),'sha256':sha(H/'coordinate_witness.json')}))
if __name__=='__main__':main()
