"""Assemble the complete B-center integral/cache discrepancy; POINT only.
The old-map-center test is not a nonexistence theorem for equilibria.
"""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
from flint import arb,arb_mat,ctx
import exact_fe_target as T
import algebraic_root_feasibility as core
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/exact_integration_bridge_v1/complete_point_B'
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def unpack(d):return core.A(F(d['lower_exact'])).union(core.A(F(d['upper_exact'])))
def abslo(d):
    a,b=F(d['lower_exact']),F(d['upper_exact'])
    return F(0) if a<=0<=b else min(abs(a),abs(b))
def main():
    ctx.prec=512
    bulkfile=ROOT/'runs/exact_integration_bridge_v1/bulk_assembly_B/results.json'
    surffile=ROOT/'runs/exact_integration_bridge_v1/surface_assembly_B/results.json'
    recordfile=ROOT/'runs/cohort_uniform_roots_v2/records/00_optoGEF2_small_col_B/result.json'
    bulk=read(bulkfile);surface=read(surffile);record=read(recordfile)
    bp=read(bulkfile.parent/'protocol_v2.json');sp=read(surffile.parent/'protocol.json')
    if bp['record_id']!=record['record_id'] or sp['point_cache_file']!=record['point_cache_file']:raise ValueError('identical saved B point required')
    if bulk['status']!='ALL144_POINT_BULK_ASSEMBLED_UNREVIEWED' or surface['status']!='ALL24_POINT_SURFACE_ASSEMBLED_UNREVIEWED' or {v['element'] for v in bulk['elements']}!=set(range(144)) or {(v['face'],v['edge']) for v in surface['edges']}!={(k,e) for k in (0,1) for e in range(12)}:raise ValueError('complete declared point contributions required')
    if not bulk['all_frozen_bindings_unchanged'] or not surface['all_frozen_bindings_unchanged']:raise ValueError('immutable point inputs required')
    if len(bulk['force_integral_minus_cache'])!=611 or len(surface['force_integral_minus_cache'])!=611:raise ValueError('matching complete free-coordinate vectors required')
    target=T.Target(read(ROOT/record['point_cache_file']));g=[unpack(a)+unpack(b) for a,b in zip(bulk['force_integral_minus_cache'],surface['force_integral_minus_cache'])]
    if len(g)!=611:raise ValueError('full611force vector required')
    H={(v['i'],v['j']):unpack(v['enclosure']) for v in bulk['bulk_H_upper']}
    for v in surface['H_upper']:
        key=(v['i'],v['j']);H[key]=H.get(key,arb(0))+unpack(v['enclosure'])
    R=arb_mat(611,611,[core.A(float(v)) for v in target.a['preconditioner'].ravel()]);Rg=R*arb_mat(611,1,g)
    packed=[core.pack(Rg[i,0]) for i in range(611)]
    force_norm_upper=max(core.upper(Rg[i,0]) for i in range(611));force_norm_lower=max(abslo(v) for v in packed)
    rows=[arb(0) for _ in range(611)]
    for (i,j),v in H.items():
        rows[i]+=core.upper(v)
        if i!=j:rows[j]+=core.upper(v)
    Hnorm=max(core.upper(v) for v in rows);Rnorm=max(core.upper(sum((core.upper(R[i,j]) for j in range(611)),arb(0))) for i in range(611))
    old=record['uniform_cell_result'];eta=F(old['eta_uniform']['upper_exact']);radius=F(old['selected_radius_exact'])
    parameters=read(ROOT/record['point_cache_file'])['parameters'];theta={k:T.rational(parameters[k]) for k in T.ORDER if k!='common_gain'}
    theta['common_gain']=F(record['base_strength_exact'])/F(record['dose_exact'])
    base_in_cell=all(F(record['parameter_cell'][k][0])<=v<=F(record['parameter_cell'][k][1]) for k,v in theta.items())
    if not base_in_cell:raise ValueError('old uniform residual must actually apply to the selected base parameter')
    map_lower=max(F(0),force_norm_lower-eta)
    paths=[Path(__file__),bulkfile,surffile,bulkfile.parent/'protocol_v2.json',surffile.parent/'protocol.json',recordfile,ROOT/record['point_cache_file'],ROOT/'evidence/exact_bulk_assembly_B_independent_review.json',ROOT/'scripts/exact_bulk_bridge_assembly.py',ROOT/'scripts/exact_surface_bridge_assembly.py']
    result={'status':'COMPLETE_POINT_DISCREPANCY_UNREVIEWED','record_id':record['record_id'],'source_and_input_hashes':{str(p.relative_to(ROOT)):sha(p) for p in paths},'all144bulk_and24surface_contributions':True,'same_base_parameter_in_old_cell':base_in_cell,'base_parameters_exact':{k:str(v) for k,v in theta.items()},'free_dofs':611,'H_upper_slots':len(H),'force_integral_minus_cache':[core.pack(v) for v in g],'preconditioned_force_defect':packed,'preconditioned_force_defect_norm_lower_exact':str(force_norm_lower),'preconditioned_force_defect_norm_upper':core.pack(force_norm_upper),'point_H_defect_row_norm_upper':core.pack(Hnorm),'point_conservative_preconditioned_H_defect_norm_upper':core.pack(Rnorm*Hnorm),'old_uniform_eta_upper_exact':str(eta),'old_radius_exact':str(radius),'new_map_center_distance_lower_exact':str(map_lower),'old_map_center_test':{'center_leaves_old_cube_certified':map_lower>radius,'proof':'At base theta in old C, ||R FI(p0)|| >= ||R(FI-FQ)(p0)|| - ||R FQ(p0)|| >= certified point-defect lower - old uniform eta. If this exceeds r, p0-R FI(p0) lies outside the old cube.','not_equilibrium_nonexistence':True},'uniform_a_or_b_established':False,'new_integral_root_certified':False,'scope':'ONE exact saved B center, shared base parameters, full integral/cache force/H difference. Neither parameter/tube-uniform nor continuum/PDE or biological fit.','missing_before_root_transfer':['Uniform full-tube derivative defects and parameter/predictor dependence.','New numerical center/proof domain if old map fails, with fresh geometry and observer verification.']}
    OUT.mkdir(parents=True,exist_ok=True);dest=OUT/'results.json'
    if dest.exists():raise ValueError('write-once point result exists')
    dest.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'status':result['status'],'RFdefect_lower':float(force_norm_lower),'RFdefect_upper':str(force_norm_upper),'new_map_center_lower':float(map_lower),'old_radius':float(radius),'center_leaves_old_cube_certified':map_lower>radius,'result_sha256':sha(dest)}),flush=True)
if __name__=='__main__':main()
