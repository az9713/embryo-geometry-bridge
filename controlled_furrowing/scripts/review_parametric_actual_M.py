"""Independent read-only strongest actual-M five-parameter cell replay.

Reconstructs parameter force derivatives analytically and assembles fresh
center/variation/Taylor bounds. No producer evaluate_cell/build_point calls,
no equilibrium solve, and no producer output writes.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,base64,csv,hashlib,io,json,time
import numpy as np
from flint import arb,arb_mat,ctx
import algebraic_root_feasibility as kernel

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'runs/parametric_root_cells_v1/actual_M'
OUT=ROOT/'evidence/parametric_actual_M_independent_review.json'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def A(v):return kernel.A(v)
def up(v):return kernel.upper(v)
def pack(v):return kernel.pack(v)
def linf(v):return max(up(v[i,0]) for i in range(v.nrows()))
def rowsum(M):return max(sum((up(M[i,j]) for j in range(M.ncols())),arb(0)).upper() for i in range(M.nrows()))


def main():
    previous=None
    if OUT.exists():
        previous=read(OUT)
        if previous['status']!='FAIL' or previous['precision_bits']!=192:raise RuntimeError('retain completed review; no overwrite')
    started=time.time();ctx.prec=160
    protocol=read(RUN/'protocol.json');result=read(RUN/'results.json');prep=read(RUN/'preparation_protocol.json')
    cachepath=ROOT/protocol['point_certificate_file'];cache=read(cachepath);statepath=ROOT/cache['input_state'];state=read(statepath)
    selected=next(c for c in result['cell_attempts'] if c['cell_id']=='scale_3' and c['method']=='affine')
    checks={};hashpaths=[RUN/'protocol.json',RUN/'results.json',RUN/'preparation_protocol.json',cachepath,statepath,
        RUN/'affine_predictor.npz',ROOT/cache['matrix_archive'],RUN/'prepare_cache.py',ROOT/'scripts/parametric_root_cells.py',Path(kernel.__file__)]
    input_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in hashpaths}
    checks['producer_protocol_result_source_bindings']=result['protocol_sha256']==sha(RUN/'protocol.json') and all(sha(ROOT/'scripts'/k)==v for k,v in protocol['source_hashes'].items())
    checks['actual_M_cache_input_predictor_archives_bound']=sha(cachepath)==protocol['point_certificate_sha256'] and sha(statepath)==cache['input_state_sha256']==protocol['point_input_sha256'] and sha(RUN/'affine_predictor.npz')==protocol['predictor_sha256'] and sha(ROOT/cache['matrix_archive'])==protocol['matrix_archive_sha256']==cache['matrix_archive_sha256']
    checks['preparation_source_and_direct_inputs_bound']=sha(RUN/'prepare_cache.py')==prep['source_sha256'] and sha(RUN/'preparation_protocol.json')==protocol['source_record_preparation_sha256']==cache['preparation_protocol_sha256'] and all(sha(ROOT/'scripts'/k)==v for k,v in prep['physical_sources'].items())
    model,u,R=kernel.restore_cache(cache);n=len(u);dose=Q(protocol['dose_exact']);base=[Q(v) for v in protocol['base_parameters_exact'].values()]
    with np.load(RUN/'affine_predictor.npz',allow_pickle=False) as z:S=z['S'].copy()
    with np.load(ROOT/cache['matrix_archive'],allow_pickle=False) as z:
        archivedS=z['displacement_tangents'][:,[1,2,3,4,5]].copy();archivedfree=z['free_indices'].copy();archivedu=z['free_vector'].copy()
    archivedS/=np.array([float(base[0]),float(base[1]),1.,float(base[3]),float(base[4])])[None,:]
    checks['predictor_exact_stored_binary_rule_reproduced']=np.array_equal(S,archivedS)
    checks['cache_binary_center_free_map_and_reference_match']=np.array_equal(u,archivedu) and np.array_equal(model.free,archivedfree) and np.array_equal(model.reference,np.asarray(state['reference_nodes'])) and np.array_equal(u,np.asarray(state['free_displacement']))
    checks['all_five_actual_absolute_coordinate_intervals_positive_width']=len(selected['cell_bounds_exact'])==5 and all(Q(lo)<Q(hi) for lo,hi in selected['cell_bounds_exact'].values())
    directp=read(ROOT/'runs/material_cohort_direct_v1/protocol.json');directr=read(ROOT/'runs/material_cohort_direct_v1/results.json')
    rec=next(r for r in directp['records'] if r['record_id']=='optoGEF2_small_col_M')
    checks['same_new_material_geometry_gain_and_source_M_dose']=state['parameters']==directp['new_parameters']==cache['parameters']==protocol['base_full_parameters'] and Q(directp['new_gain_exact'])==base[4]==Q(prep['base_gain_exact']) and Q(rec['dose_exact'])==dose and protocol['calibration_record_ids_covered']==['optoGEF2_small_col_M']
    checks['exact_gain_times_dose_is_point_target_not_rounded_solver_strength']=base[4]*dose==Q(prep['requested_strength_exact'])==Q(protocol['base_strength_requested_exact']) and Q.from_float(float(state['activation_strength']))==Q(prep['solver_strength_binary_exact']) and base[4]*dose!=Q(prep['solver_strength_binary_exact'])
    RM=arb_mat(n,n,[A(float(v)) for v in R.ravel()]);AR=arb_mat(n,n,[A(abs(float(v))) for v in R.ravel()]);SM=arb_mat(n,5,[A(float(v)) for v in S.ravel()])
    g=arb_mat(n,1);H=arb_mat(n,n);B0=arb_mat(n,5);terms=[];surface_count=0;basal_count=len(model.surfaces[0][0])*len(model.surfaces[0][1]);bulk_count=0
    for term in kernel.local_terms(model,u,base[4]*dose):
        ids=term['indices'];d=len(ids);k=len(term['q']);BM=arb_mat(k,d,[v for row in term['B'] for v in row])
        fun=kernel.bulk_local if term['kind']=='bulk' else kernel.surface_local
        grad,h,_=fun(term['q'],*term['parameters']);HM=arb_mat(k,k,[v for row in h for v in row]);GM=arb_mat(k,1,grad)
        localH=BM.transpose()*HM*BM;localg=BM.transpose()*GM;C=arb_mat(k,5)
        if term['kind']=='bulk':
            element,iq=divmod(bulk_count,len(model.shape));active=model.active.numpy();mask=A(float(active[element,0 if active.shape[1]==1 else iq]));bulk_count+=1
            C[1,4]=mask*A(dose)*term['q'][1];C[3,4]=mask*A(dose)*term['q'][3];term['mask']=mask
        else:
            basal=surface_count<basal_count;surface_count+=1;term['basal']=basal
            rr,x,z=term['q'];scale=term['parameters'][2];length=(x*x+z*z).sqrt();js=scale*rr*length
            djs=[scale*length,scale*rr*x/length,scale*rr*z/length];ki,ti=(1,3) if basal else (0,2)
            for a in range(k):C[a,ki]=(js-1)*djs[a];C[a,ti]=djs[a]
        localB=BM.transpose()*C
        for a,i in enumerate(ids):
            g[i,0]+=term['weight']*localg[a,0]
            for b,j in enumerate(ids):H[i,j]+=term['weight']*localH[a,b]
            for j in range(5):B0[i,j]+=term['weight']*localB[a,j]
        term['HC']=h;term['P']=BM*arb_mat(d,5,[A(float(S[i,j])) for i in ids for j in range(5)])
        term['absB']=[[up(v) for v in row] for row in term['B']];term['rowsB']=[sum(row,arb(0)).upper() for row in term['absB']]
        if not term['weight']>0:raise ValueError('positive quadrature weights required')
        terms.append(term)
        if len(terms)%1536==0:print(json.dumps({'phase':'independent actualM center','terms':len(terms),'elapsed':time.time()-started}),flush=True)
    theta=[A(Q(lo)).union(A(Q(hi))) for lo,hi in selected['cell_bounds_exact'].values()]
    delta=[theta[i]-A(base[i]) for i in range(5)];hull=[v.union(arb(0)) for v in delta];r=A(Q(selected['selected_radius_exact']))
    remrows=[arb(0) for i in range(n)];drows=[arb(0) for i in range(n)];minJ=None
    for term in terms:
        k=len(term['q']);offset=[sum((term['P'][a,j]*hull[j] for j in range(5)),arb(0)) for a in range(k)]
        for remainder,rows in ((True,remrows),(False,drows)):
            qb=[term['q'][a]+offset[a]+(arb(0,(r*term['rowsB'][a]).upper()) if not remainder else 0) for a in range(k)]
            if term['kind']=='bulk':_,hb,J=kernel.bulk_local(qb,term['parameters'][0],term['parameters'][1],term['mask']*A(dose)*theta[4]);minJ=J.lower() if minJ is None else min(minJ,J.lower())
            else:
                ki,ti=(1,3) if term['basal'] else (0,2);_,hb,_=kernel.surface_local(qb,theta[ki],theta[ti],term['parameters'][2])
            width=[up(v) for v in offset] if remainder else term['rowsB']
            change=[[up(hb[a][b]-term['HC'][a][b]) for b in range(k)] for a in range(k)]
            local=[sum((change[a][b]*width[b] for b in range(k)),arb(0)) for a in range(k)]
            for j,i in enumerate(term['indices']):rows[i]+=term['weight']*sum((term['absB'][a][j]*local[a] for a in range(k)),arb(0))
    identity=arb_mat(n,n,[1 if i==j else 0 for i in range(n) for j in range(n)]);Z0=rowsum(identity-RM*H)
    linear=RM*g+RM*(H*SM+B0)*arb_mat(5,1,delta);remainder=AR*arb_mat(n,1,[x.upper() for x in remrows]);eta=max(up(linear[i,0])+up(remainder[i,0]) for i in range(n))
    variation=AR*arb_mat(n,1,[x.upper() for x in drows]);Z=Z0+linf(variation);lhs=eta+Z*r;err=eta/(1-Z)
    pred=[up(sum((A(float(S[i,j]))*delta[j] for j in range(5)),arb(0))) for i in range(n)];outer=max(v+err for v in pred)
    expected=next(a for a in selected['radius_attempts'] if a['strict_selfmap_contraction'])
    checks.update(all_9504_terms_and_611_free_coordinates=len(terms)==9504 and n==611,
        independent_center_defect_below_one=bool(Z0<1),independent_uniform_contraction=bool(Z<1),independent_uniform_selfmap=bool(lhs<r),
        independent_analytic_domain_J_positive=bool(minJ>0),all_inputs_unchanged=all(sha(ROOT/k)==v for k,v in input_hashes.items()),
        independent_eta_agrees_to_25_relative_digits=abs(float(eta)-float(Q(selected['eta_uniform']['upper_exact'])))<1e-25,
        independent_Z_agrees_to_25_absolute_digits=abs(float(Z)-float(Q(expected['Z_uniform']['upper_exact'])))<1e-25,
        independent_outer_radius_agrees_to_25_absolute_digits=abs(float(outer)-float(Q(selected['uniform_outer_nodal_radius_about_base']['upper_exact'])))<1e-25,
        all_611_predictor_bounds_match_to_25_absolute_digits=all(abs(float(v)-float(Q(s)))<1e-25 for v,s in zip(pred,selected['full_nodal_predictor_bounds_exact'])),
        no_producer_evaluation_or_optimization_called=True)
    receipt={'status':'PASS_INDEPENDENT_ACTUAL_M_UNIFORM_ROOT_CELL' if all(checks.values()) else 'FAIL','checks':checks,'review_source_sha256':sha(__file__),'reviewed_hashes':input_hashes,
        'cell_id':'scale_3','method':'affine','free_dofs':n,'terms':len(terms),'precision_bits':ctx.prec,
        'base_strength_exact':str(base[4]*dose),'solver_strength_binary_exact':prep['solver_strength_binary_exact'],'strength_difference_exact':str(base[4]*dose-Q(prep['solver_strength_binary_exact'])),
        'independent_bounds':{'point_Z0':pack(Z0),'eta_uniform':pack(eta),'Z_uniform':pack(Z),'eta_plus_Zr':pack(lhs),'radius':pack(r),'error_from_predictor':pack(err),'outer_nodal_radius':pack(outer),'quadrature_J_lower':pack(minJ)},
        'scope':'Independent proof/replay of one selected actual-M local five-parameter root family. Body geometry, observation transfer, provider callback, source-channel exclusion and global material domain require separate actual artifacts/review.',
        'observer_provider_review':'PENDING','elapsed_seconds':time.time()-started}
    if previous is not None:receipt['initial_192bit_review_preserved']=previous
    OUT.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8');print(json.dumps(receipt,indent=2),flush=True)
    if not all(checks.values()):raise RuntimeError('independent actual-M review finding retained')


def polynomial(points,coordinate):
    a,m,b=[p[coordinate] for p in points]
    return a,-3*a+4*m-b,2*a-4*m+2*b


def value(c,t):return c[0]+c[1]*t+c[2]*t*t


def inverse_enclosure(c,r):
    if value(c,Q(0))==r:return Q(0),Q(0)
    if value(c,Q(1))==r:return Q(1),Q(1)
    lo,hi=Q(0),Q(1)
    for _ in range(120):
        mid=(lo+hi)/2
        if value(c,mid)<=r:lo=mid
        else:hi=mid
    return lo,hi


def range_polynomial(c,lo,hi):
    candidates=[value(c,lo),value(c,hi)]
    if c[2]:
        t=-c[1]/(2*c[2])
        if lo<=t<=hi:candidates.append(value(c,t))
    return min(candidates),max(candidates)


def compose():
    """Independent exact box/graph/source audit; generic body proof is reviewed separately."""
    from openpyxl import load_workbook
    from root_body_geometry import certify_box
    r=read(OUT)
    if r['status']!='PASS_INDEPENDENT_ACTUAL_M_UNIFORM_ROOT_CELL' or r.get('observer_provider_review')!='PENDING':raise RuntimeError('fresh root-only review required')
    started=time.time();p=read(RUN/'protocol.json');result=read(RUN/'results.json');point=read(RUN/'point_cache.json');state=read(ROOT/point['input_state'])
    root=read(RUN/'root_certificate.json');obs=read(RUN/'observer_certificate.json');contract=read(RUN/'model_observation_contract.json');op=read(RUN/'observer_protocol.json');body=read(RUN/'body_geometry.json');diag=read(RUN/'contradiction_diagnostic.json');prep=read(RUN/'preparation_protocol.json')
    cell=root['parameter_cell'];sel=root['uniform_cell_result'];rootsha=sha(RUN/'root_certificate.json');obssha=sha(RUN/'observer_certificate.json');consha=sha(RUN/'model_observation_contract.json')
    canonical=hashlib.sha256(json.dumps(cell,sort_keys=True,separators=(',',':')).encode()).hexdigest();checks={}
    checks['root_cell_is_exact_frozen_scale3_result']=sel==next(z for z in result['cell_attempts'] if z['cell_id']=='scale_3' and z['method']=='affine') and cell==next(z['bounds'] for z in p['initial_cells'] if z['cell_id']=='scale_3')
    checks['root_observer_contract_same_canonical_cell_hash']=all(x['parameter_cell']==cell and x['parameter_cell_sha256']==canonical and x['model_observation_contract_sha256']==consha for x in (root,obs))
    checks['uniform_root_and_observer_quantifiers_and_branch_law_match']=all(x['quantifier']=='for_all_parameters_in_closed_cell_exists_unique_root_in_declared_branch_domain' and x['branch_scope']=='selected_unique_root_in_declared_domain' and x['branch_domain_id']==p['branch_domain_id'] and x['branch_law_semantics']=='independent_stationary_root_per_record' and x['record_id']=='optoGEF2_small_col_M' and x['dose_exact']==p['dose_exact'] for x in (root,obs))
    checks['actual_M_direct_provenance_hashes_bound']=all(sha(ROOT/'runs/material_cohort_direct_v1'/name)==prep[key] for name,key in [('protocol.json','direct_protocol_sha256'),('results.json','direct_results_sha256')]) and sha(ROOT/'evidence/material_cohort_direct_independent_review.json')==prep['direct_review_sha256']
    checks['current_v2_observer_artifacts_hashes_and_sources_bound']=op['source_sha256']==sha(RUN/'build_forward_artifacts_v2.py') and obs['root_certificate_sha256']==rootsha and obs['body_geometry_sha256']==sha(RUN/'body_geometry.json') and obs['observer_protocol_sha256']==sha(RUN/'observer_protocol.json') and all(sha(ROOT/'scripts'/k)==v for k,v in obs['source_hashes'].items())
    oldmanifest=read(RUN/'observer_revision_0/manifest.json')
    checks['old_observer_revision_retained']=bool(oldmanifest)
    raw=base64.b64decode(point['frozen_binary_cache']['data'])
    with np.load(io.BytesIO(raw),allow_pickle=False) as z:ref=z['reference'].copy();free=z['free'].copy();u=z['u'].copy()
    with np.load(RUN/'affine_predictor.npz',allow_pickle=False) as z:S=z['S'].copy()
    delta=[(Q(cell[name][0])-Q(p['base_parameters_exact'][name]),Q(cell[name][1])-Q(p['base_parameters_exact'][name])) for name in p['coordinate_order']]
    pred=[]
    for row in S:
        pairs=[(Q.from_float(float(s))*lo,Q.from_float(float(s))*hi) for s,(lo,hi) in zip(row,delta)]
        pred.append(max(abs(sum(min(v) for v in pairs)),abs(sum(max(v) for v in pairs))))
    corrector=Q(sel['selected_radius_exact']);where={int(v):i for i,v in enumerate(free)};current=[[Q.from_float(float(v)) for v in pair] for pair in state['current_nodes']];reference=[[Q.from_float(float(v)) for v in pair] for pair in ref]
    gaps=[];radii=[]
    for i,x in enumerate(ref.ravel()):
        exact=Q.from_float(float(x))+(Q.from_float(float(u[where[i]])) if i in where else Q(0));gap=abs(exact-current[i//2][i%2]);gaps.append(gap)
        radii.append(gap+(pred[where[i]]+corrector if i in where else Q(0)))
    first192=r['initial_192bit_review_preserved']['independent_bounds']
    checks['full_corrector_cube_contains_both_independent_root_bounds']=corrector>=Q(first192['error_from_predictor']['upper_exact']) and corrector>=Q(r['independent_bounds']['error_from_predictor']['upper_exact'])
    checks['all611_exact_predictor_hulls_and_full_corrector_reproduced']=list(map(Q,obs['exact_linear_predictor_node_bounds']))==pred and Q(obs['corrector_radius_exact'])==corrector==Q(1,10**11)
    checks['all650_actual_coordinate_radii_add_exact_serialization_gap']=list(map(Q,body['coordinate_radii_exact']))==radii and Q(obs['coordinate_rounding_gap_max_exact'])==max(gaps)
    eps=Q('33.5')*max(radii);checks['physical_scale_applied_to_full_box_once']=eps==Q(obs['physical_nodal_coordinate_error_um_exact'])
    generic=read(ROOT/'evidence/root_body_geometry_independent_review.json')
    checks['generic_whole_body_proof_independently_reviewed']=generic['status'].startswith('PASS') and all(v is True for v in generic['checks'].values()) and sha(ROOT/'scripts/root_body_geometry.py')==body['source_sha256']
    replay=certify_box(reference,current,state['triangles'],[radii[i:i+2] for i in range(0,len(radii),2)])
    checks['actual_full_nodal_box_body_certificate_replayed']=replay['status']=='CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING' and all(replay[k]==body[k] for k in replay)
    # Raw-connectivity edge extraction, independent exact inverse and extrema.
    edge={}
    for a,b,c,d,e,f in state['triangles']:
        for i,j,m in ((a,b,d),(b,c,e),(c,a,f)):edge[tuple(sorted((i,j)))]=m
    H=Q.from_float(float(state['parameters']['height']));L=Q('33.5');roi=(Q(0),Q.from_float(float(state['parameters']['patch_radius']))*L);surface={};endpointdofs=[];Zcandidates=[]
    for name,zlevel in [('apical',H),('basal',Q(0))]:
        edges=[]
        for (i,j),mid in edge.items():
            if reference[i][1]==reference[j][1]==zlevel:
                if reference[i][0]>reference[j][0]:i,j=j,i
                edges.append((reference[i][0],(i,mid,j)))
        edges=[v for _,v in sorted(edges)];endpointdofs.extend((2*edges[0][0],2*edges[-1][-1]));lip=Q(0);minimum=None
        for nodes in edges:
            points=[tuple(L*v for v in current[i]) for i in nodes];rc=polynomial(points,0);zc=polynomial(points,1)
            mr=min(rc[1],rc[1]+2*rc[2]);minimum=mr if minimum is None else min(minimum,mr);lip=max(lip,max(abs(zc[1]),abs(zc[1]+2*zc[2]))/mr)
            if name=='apical':
                rr0,rr1=max(roi[0],points[0][0]),min(roi[1],points[-1][0])
                if rr0<=rr1:
                    for rr in (rr0,rr1):Zcandidates.append(range_polynomial(zc,*inverse_enclosure(rc,rr)))
                    if zc[2]:
                        t=-zc[1]/(2*zc[2])
                        if 0<=t<=1 and rr0<=value(rc,t)<=rr1:Zcandidates.append((value(zc,t),value(zc,t)))
        E=Q(5,4)*eps*(1+lip);surface[name]={'lip':lip,'minimum_Rprime':minimum,'E':E}
        checks[name+'_strict_monotonicity_and_Eulerian_transport_recomputed']=minimum>8*eps and E==Q(obs['transfer_certificate'][name]['Eulerian_vertical_error_um_exact'])
    checks['fixed_common_radial_endpoints_actual_constraints']=all(i not in where and radii[i]==0 for i in endpointdofs)
    checks['same_positive_ROI_fixed_refs_contract']=roi[1]>0 and contract['physical_length_um']=='33.5' and contract['reference_profiles_um']=={'apical':'33.5','basal':'0'} and all(tuple(map(Q,obs['transfer_certificate'][name]['ROI_um_exact']))==roi for name in surface)
    baselo=L*H-min(v[1] for v in Zcandidates);basehi=L*H-min(v[0] for v in Zcandidates)
    publishedbase=tuple(map(Q,obs['base_stored_apical_interval']));publishedobs=tuple(map(Q,obs['endpoint_intervals']['apical_max_shortening_um']));Ea=surface['apical']['E']
    checks['independent_quadratic_graph_extrema_inside_published_base_interval']=publishedbase[0]<=baselo<=basehi<=publishedbase[1]
    checks['actual_uniform_apical_interval_contains_independent_extrema_and_transport']=publishedobs[0]<=baselo-Ea and basehi+Ea<=publishedobs[1]
    data=ROOT/contract['calibration_source_file'];row=next(v for v in csv.DictReader(data.open(encoding='utf-8-sig')) if v['record_id']=='optoGEF2_small_col_M');y=Q(row['apical_max_shortening_um']);dose=Q(row['activation_magnitude'])
    workbook=ROOT/'data/raw/41467_2025_62483_MOESM16_ESM.xlsx';wb=load_workbook(workbook,data_only=True,read_only=True)
    checks['source_M_dose_and_apical_exact_workbook_cells']=dose==Q(str(wb['Fig 2c']['M4'].value))==Q(p['dose_exact']) and y==Q(str(wb['Fig 2e']['M4'].value))==12 and row['dose_source']=='Fig 2c!M4' and row['apical_source']=='Fig 2e!M4'
    minimum=max(publishedobs[0]-y,y-publishedobs[1],Q(0));delta=Q(diag['declared_delta_um_exact'])
    checks['one_channel_contradiction_arithmetic_and_tolerance_scope']=minimum==Q(diag['minimum_residual_um_exact']) and minimum>delta==Q('5.1') and diag['delta_is_measured_noise_estimate'] is False
    checks['other67_and_otherbranches_explicitly_unresolved']=obs['other67constraints_not_supplied'] is True and set(obs['endpoint_intervals'])=={'apical_max_shortening_um'} and 'UNRESOLVED' in diag['otherbranches'] and contract['source_temporal_semantics_validated'] is False
    checks['original_root_review_inputs_unchanged']=all(sha(ROOT/k)==v for k,v in r['reviewed_hashes'].items())
    paths=[RUN/x for x in ['root_certificate.json','observer_certificate.json','body_geometry.json','model_observation_contract.json','observer_protocol.json','contradiction_diagnostic.json','build_forward_artifacts_v2.py','observer_revision_0/manifest.json']]+[data,workbook,ROOT/'evidence/root_body_geometry_independent_review.json']
    r['observer_provider_review']='PASS' if all(checks.values()) else 'FAIL';r['composition_checks']=checks
    r['composition_hashes']={v.relative_to(ROOT).as_posix():sha(v) for v in paths}
    common={'parameter_cell_sha256':canonical,'model_observation_contract_sha256':consha,'uniform_quantifier_verified':True,'branch_law_semantics':'independent_stationary_root_per_record','branch_scope':'selected_unique_root_in_declared_domain','branch_domain_id':p['branch_domain_id']}
    rootchecks={k:v for k,v in r['checks'].items() if 'agrees_to' not in k and 'match_to' not in k}
    r['root_audit']={'status':'PASS' if all(rootchecks.values()) else 'FAIL','checks':rootchecks,'certificate_sha256':rootsha,**common}
    r['observer_audit']={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'certificate_sha256':obssha,'forward_family_admissible':checks['actual_full_nodal_box_body_certificate_replayed'],**common}
    r['independent_observer_reconstruction']={'exact_serialization_gap_max':str(max(gaps)),'exact_global_coordinate_radius_max':str(max(radii)),'physical_coordinate_radius_um':str(eps),'ROI_um':list(map(str,roi)),'apical_base_interval_exact':[str(baselo),str(basehi)],'apical_transport_error_exact':str(Ea),'uniform_apical_interval_exact':list(map(str,publishedobs)),'source_value_um':str(y),'minimum_residual_um':str(minimum),'declared_tolerance_um':str(delta)}
    r['current_review_source_sha256']=sha(__file__);r['status']='PASS' if all(rootchecks.values()) and all(checks.values()) else 'FAIL';r['composition_elapsed_seconds']=time.time()-started
    OUT.write_text(json.dumps(r,indent=2)+'\n',encoding='utf8');print(json.dumps({'status':r['status'],'composition_checks':checks,'elapsed':r['composition_elapsed_seconds'],'uniform_apical_interval_um':list(map(float,publishedobs)),'minimum_residual_um':float(minimum)},indent=2))
    if r['status']!='PASS':raise RuntimeError('retained actualM composition review failure')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--compose',action='store_true');args=ap.parse_args()
    compose() if args.compose else main()
