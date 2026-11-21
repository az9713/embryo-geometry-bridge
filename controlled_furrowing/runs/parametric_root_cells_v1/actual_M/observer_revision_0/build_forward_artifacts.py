"""Bind one uniform actual-M root cell to exact geometry/static apical bounds."""
from pathlib import Path
from fractions import Fraction as F
import base64,csv,hashlib,io,json,sys
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from quadratic_endpoint_observer import _state_controls
from quadratic_mean_derivatives import _boundary_node_ids
from nodal_observer_error_transfer import endpoint_bounds
from root_body_geometry import certify_box
from joint_gain_calibration import endpoints

sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
def save(p,v):
    if Path(p).exists():raise RuntimeError('preserve existing artifact:'+str(p))
    Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf8')
def rational(x):return F.from_float(x)if isinstance(x,float)else F(x)

protocol=read(OUT/'protocol.json');results=read(OUT/'results.json');point=read(OUT/'point_cache.json');prep=read(OUT/'preparation_protocol.json')
assert results['status']=='COMPLETE_DECLARED_CELL_ATTEMPTS'and protocol['calibration_record_ids_covered']==['optoGEF2_small_col_M']
affine=[r for r in results['cell_attempts']if r['method']=='affine'and r['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT'and r['cell_id'].startswith('scale_')]
scale={r['cell_id']:F(r['relative_scale_exact'])for r in protocol['initial_cells']}
cell=max(affine,key=lambda r:scale[r['cell_id']]);assert cell['cell_id']=='scale_3'
parameter_cell=cell['cell_bounds_exact'];cellhash=hashlib.sha256(json.dumps(parameter_cell,sort_keys=True,separators=(',',':')).encode()).hexdigest()
state=read(ROOT/point['input_state']);assert sha(ROOT/point['input_state'])==point['input_state_sha256']
names=['parametric_root_cells.py','algebraic_root_feasibility.py','quadratic_endpoint_observer.py','quadratic_mean_derivatives.py','nodal_observer_error_transfer.py','root_body_geometry.py','annulus_observer.py','joint_gain_calibration.py','mixed_geometry_certificate.py','mixed_geometry_certificate_v2.py','certify_discrete_geometry.py']
sources={name:sha(ROOT/'scripts'/name)for name in names}
contract={'version':1,'model':'finite-height axisymmetric full-pointwise P2 solid with two area-elastic/tension surfaces',
 'potential_target':'exact real smooth potential with frozen binarycache geometry/quadrature/coefficients, exactX+Pu; common_gain*dose is exactrational inputlaw',
 'source_kernel_hashes':{**prep['physical_sources'],**sources},'point_cache_sha256':sha(OUT/'point_cache.json'),
 'base_full_parameters':state['parameters'],'fixed_mu':'1','fixed_kappa':'5000','kappa_role':'fixed method/constitutive coefficient, not inferred physiological compressibility',
 'physical_length_um':'33.5','volume_quadrature':'Duffy8 frozen binary nodes/weights','surface_quadrature':'order12 frozen binary nodes/weights',
 'geometry_uncertainty':'none in this CONDITIONAL cell; fixed flat axisymmetric reference, activation support and boundary constraints',
 'record_id':prep['record_id'],'record_identity':'source-column locator, not an independently verified raw embryoID',
 'dose_exact':prep['dose_exact'],'strength_relation':'common_gain_times_source_dose','dose_error_assumption':'canonical workbookextract dose fixed; not a measured zeroerror statement',
 'parameter_coordinate_order':protocol['coordinate_order'],'parameter_cell':parameter_cell,
 'reference_profiles_um':{'apical':'33.5','basal':'0'},'time_observer':'single static endpoint proxy; NOT experimental temporal maximum/localpeak-average reconstruction',
 'roi':'same fixed registered current physical radius [0,patch_radius*33.5] for both surfaces',
 'height_definition':'reference_apical-current_apical + current_basal-reference_basal atsamecurrentradius',
 'spatial_weighting':'projected_disk: 2*r/R^2 dr for regionalmeans; this differs from line weighting dr/R',
 'component_aggregation':'spatial_max ofsignedapicaldown/basalup; not regionalmeans, no localpositivepartclipping',
 'primary_variant':'current_disk__spatial_max','source_component_pairing':'canonical B:R crosssheetorder assumed; no independent biologicalID crosswalk',
 'branch_law_semantics':'independent_stationary_root_per_record','branch_scope':'selected_unique_root_in_declared_domain','branch_domain_id':protocol['branch_domain_id'],
 'common_path_or_energetic_stability_certified':False,'source_temporal_semantics_validated':False,
 'calibration_source_file':'data/derived/calibration_small_18um.csv','calibration_source_sha256':sha(ROOT/'data/derived/calibration_small_18um.csv'),
 'global_domain_warning':'this is one small positive coefficient box; zero faces, other parameters/geometries/branches and unbounded tails remain outsidecoverage'}
save(OUT/'model_observation_contract.json',contract);contracthash=sha(OUT/'model_observation_contract.json')
observer_protocol={'version':1,'stage':'Frozenbeforegeometry/staticobserverevaluation','source_sha256':sha(__file__),
 'uniform_results_sha256':sha(OUT/'results.json'),'uniform_protocol_sha256':sha(OUT/'protocol.json'),
 'selected_cell_id':cell['cell_id'],'selection_rule':'largest prefrozen relative-width affine cell with successful uniformroot test; no observation-based widthselection',
 'parameter_cell_sha256':cellhash,'model_observation_contract_sha256':contracthash,'source_hashes':sources,
 'scope':'one actualM apical static endpoint only; no directdecision/full68IN/globalparametercoverage'}
save(OUT/'observer_protocol.json',observer_protocol)
common={'parameter_cell':parameter_cell,'parameter_cell_sha256':cellhash,'model_observation_contract_sha256':contracthash,
 'quantifier':'for_all_parameters_in_closed_cell_exists_unique_root_in_declared_branch_domain',
 'branch_scope':'selected_unique_root_in_declared_domain','branch_domain_id':protocol['branch_domain_id'],
 'branch_law_semantics':'independent_stationary_root_per_record','common_path_or_energetic_stability_certified':False,
 'input_kind':'source_record_dose','record_id':prep['record_id'],'dose_exact':prep['dose_exact'],'strength_relation':'common_gain_times_source_dose'}
root={'status':'CERTIFIED_UNIFORM_LOCAL_ROOT',**common,'uniform_cell_result':cell,
 'uniform_results_sha256':sha(OUT/'results.json'),'uniform_protocol_sha256':sha(OUT/'protocol.json'),
 'point_cache_sha256':sha(OUT/'point_cache.json'),'source_hashes':protocol['source_hashes'],
 'base_input_state_sha256':point['input_state_sha256'],'observer_enclosures_supplied_separately':True,
 'scope':'selected nearby stationary branch at actualM input; NOT all17, allbranches, physicalprior or fullidentifiedset'}
save(OUT/'root_certificate.json',root)
blob=base64.b64decode(point['frozen_binary_cache']['data']);assert hashlib.sha256(blob).hexdigest()==point['frozen_binary_cache']['sha256']
with np.load(io.BytesIO(blob),allow_pickle=False)as z:reference=z['reference'];free=z['free'];u=z['u']
selected={int(index):i for i,index in enumerate(free)};pred=list(map(F,cell['full_nodal_predictor_bounds_exact']));error=F(cell['uniform_error_from_predictor']['upper_exact'])
exact=[];radii=[];gap=F(0)
for j,x in enumerate(reference.ravel()):
    v=F.from_float(float(x))+(F.from_float(float(u[selected[j]]))if j in selected else F(0));exact.append(v)
    rounded=F.from_float(float(np.asarray(state['current_nodes']).ravel()[j]));d=abs(v-rounded);gap=max(gap,d)
    radii.append(d+(pred[selected[j]]+error if j in selected else F(0)))
current=[[F.from_float(float(v))for v in row]for row in state['current_nodes']];refs=[[F.from_float(float(v))for v in row]for row in reference]
body=certify_box(refs,current,state['triangles'],[radii[j:j+2]for j in range(0,len(radii),2)])
body.update({'root_certificate_sha256':sha(OUT/'root_certificate.json'),'source_sha256':sources['root_body_geometry.py'],'coordinate_radii_exact':[str(v)for v in radii]})
save(OUT/'body_geometry.json',body)
if body['status']!='CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING':raise RuntimeError('uniformrootfamilybodygeometry remainsUNRESOLVED')
ids=_boundary_node_ids(state);endpoint_dofs=[]
for name in('top','bottom'):endpoint_dofs.extend([2*ids[name][0][0],2*ids[name][-1][-1]])
assert all(i not in selected and radii[i]==0 for i in endpoint_dofs)
L=F('33.5');eps=L*max(radii);controls,H=_state_controls(state,L);roi=(F(0),F.from_float(float(state['parameters']['patch_radius']))*L)
transfer=endpoint_bounds(controls['top'],controls['bottom'],eps,fixed_radial_endpoints=True,roi_um=roi)
assert transfer['status']=='CONDITIONAL_EXACT_ERROR_TRANSFER'
base=endpoints(state)['current_disk__spatial_max']['apical_max_shortening_um'];ea=F(transfer['endpoint_error_um_exact']['apical_max_shortening_um']);bounds=[str(F(base[0])-ea),str(F(base[1])+ea)]
observer={'status':'CERTIFIED_UNIFORM_CELL_OBSERVABLES',**common,'root_certificate_sha256':sha(OUT/'root_certificate.json'),
 'body_geometry_sha256':sha(OUT/'body_geometry.json'),'observer_protocol_sha256':sha(OUT/'observer_protocol.json'),
 'endpoint_intervals':{'apical_max_shortening_um':bounds},'other67constraints_not_supplied':True,
 'coordinate_rounding_gap_max_exact':str(gap),'physical_nodal_coordinate_error_um_exact':str(eps),
 'transfer_certificate':transfer,'base_stored_apical_interval':base,'source_hashes':sources,
 'scope':'allparameterroots inthisselectednearbybranch; apical static endpoint only. No experimentaltemporaloperator/fullTheta/commonpath/stabilityclaim'}
save(OUT/'observer_certificate.json',observer)
row=next(r for r in csv.DictReader((ROOT/contract['calibration_source_file']).open())if r['record_id']==prep['record_id']);y=F(row['apical_max_shortening_um']);delta=F('5.1');minimum=max(F(bounds[0])-y,y-F(bounds[1]),F(0))
diagnostic={'status':'PRODUCER_BRANCH_CONDITIONAL_CONTRADICTION'if minimum>delta else'UNRESOLVED','root_certificate_sha256':sha(OUT/'root_certificate.json'),
 'observer_certificate_sha256':sha(OUT/'observer_certificate.json'),'declared_delta_um_exact':str(delta),'delta_is_measured_noise_estimate':False,
 'observed_apical_um_exact':str(y),'uniform_apical_interval_exact':bounds,'minimum_residual_um_exact':str(minimum),
 'independent_review_required_before_cover_engine_use':True,'other67constraints':'UNRESOLVED/not supplied','otherbranches':'UNRESOLVED/not excluded'}
save(OUT/'contradiction_diagnostic.json',diagnostic)
print(json.dumps({'root_sha256':sha(OUT/'root_certificate.json'),'observer_sha256':sha(OUT/'observer_certificate.json'),
 'model_contract_sha256':contracthash,'cell_sha256':cellhash,'body_status':body['status'],
 'apical_bounds_um':list(map(float,map(F,bounds))),'minimum_residual_um':float(minimum),'declared_delta_um':float(delta),
 'only_selected_branch_excluded_by_producer':minimum>delta}),flush=True)
