from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,datetime
P=Path(__file__).resolve().parent
ROOT=P.parents[2]
def rd(p): return json.loads(p.read_text(encoding='utf-8'))
def hs(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def wr(p,x): p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')
r=rd(P/'results.json');p=rd(P/'protocol.json')
assert all(v is True for v in r['checks'].values())
assert r['source_sha256']==hs(P/'crosscheck.py') and r['protocol_sha256']==hs(P/'protocol.json')
assert all(hs(ROOT/k)==v for k,v in p['input_hashes'].items())
passed=r['strict_selfmap_pass'];error=None if not passed else Q(r['root_error_exact'])
led=rd(ROOT/'runs/exact_integration_bridge_v1/retube_B_preparation/rounding_ledger.json')
d588=Q(led['applied_decrement_exact'][588]);radius=Q(r['old_radius_exact']);zq=Q(r['ZQ_exact'])
eq=Q(r['old_eta_exact'])/(1-zq)
location=bool(passed and abs(d588)-error>radius)
sep=None if not passed else max(Q(0),abs(d588)-error-eq)
consequence={'uniform_for_all_theta_in_same_C':True,'same_domain_uniqueness_by_contraction':Q(r['Z_integral_exact'])<1,'coordinate_index':588,'applied_shift_coordinate_exact':str(d588),'integral_root_error_about_new_predictor_exact':None if error is None else str(error),'old_cached_Q_root_error_about_old_predictor_exact':str(eq),'integral_root_outside_old_closed_cube':location,'no_exact_integral_root_in_old_cube':location,'no_integral_root_in_old_cube_reason':'Existence in new cube plus strict contraction on convex enclosing common cube gives uniqueness throughout that common cube; new root is outside old cube by coordinate588. This conclusion does NOT follow from old-map-center escape alone.','Q_vs_integral_root_distance_lower_exact':None if sep is None else str(sep),'scope':'Uniform same B/C and finite-P2 exact-integration potential only; no global uniqueness, PDE, nonaxisymmetric stability or biology.'}
wr(P/'location_consequences.json',consequence)
proof=r'''# Direct uniform-force crosscheck of the integral-root family

This is a separately executed exact-rational scalar composition. It shares the reviewed integration/component kernels and their enclosures with the principal proof. It is not an independent physical implementation. No integration, Hessian assembly, matrix product, optimization or equilibrium solve was repeated.

Write p0(theta)=u0+S(theta-theta0), p1(theta)=p0(theta)-d, and let R be the fixed nonsingular preconditioner. The old finite-quadrature result supplies ||R FQ(p0(theta),theta)||∞ ≤ etaQ for every theta in the same closed C. The new uniform component enclosure supplies Dbox containing R(FI-FQ)(u,theta) on the entire reviewed common nodal box and C. The independent endpoint reduction retains signs of all 611 preconditioned force coordinates. Use b from the fully exact rational |R|-times-Hessian-row majorant, rather than the tighter signed-Arb matrix bound.

Define m=max_i sup |Dbox_i-d_i|, z=zQ+b, and dn=||d||∞. All bounds are exact rational numbers. The common radius is rc=r+dn. Thus both old/new radius-r cubes and the complete line between their centers are in the common domain.

The fundamental theorem of calculus gives

    R FI(p1)=R FQ(p0)+[R(FI-FQ)(p0)-d]
                + integral_0^1 [I-R HI(p0-t*d)] d dt.

Consequently eta_direct=etaQ+m+z*dn bounds ||R FI(p1)||∞ uniformly in theta. On the whole new radius-r cube, T(u)=u-R FI(u) has derivative norm at most z. Strict eta_direct+z*r<r and z<1 prove existence and uniqueness in that cube by Banach, with error ≤eta_direct/(1-z). R is nonsingular by the already verified center defect ||I-R HQ||<1. Thus a fixed point is an FI zero. The geometry and analytic-domain certificates ensure differentiability throughout these cubes; they do not assert a continuum PDE solution.

Actual scalar bounds (decimal displays only): m=1.284615006e-14, z=0.1050417212, eta_direct=2.299666996e-12. The self-map left side is below 3.350084208e-12, strictly less than r=1e-11. The resulting error is below 2.569580114e-12. Acceptance used the complete rational fractions in results.json, not these displays.

There is also at most one FI zero anywhere in the convex common cube: two zeros are fixed points and their distance would be at most z times itself. This uniqueness argument does not require a self-map of the common cube. Coordinate588 satisfies |d588|-error>r, so the newly certified root lies outside the old closed cube. Since the old cube is inside the common cube, there is no exact-integral FI zero in the old cube. This is a new consequence of existence plus common-domain uniqueness plus location, not a consequence of the earlier map-center escape alone. The known old Q root has error at most etaQ/(1-zQ) about p0; the distance between the Q and integral roots is at least max(0,|d588|-error-etaQ/(1-zQ)). Exact values are retained in location_consequences.json.

The result concerns one record B, the declared tiny five-parameter C, fixed geometry/mesh and axisymmetric finite-P2 energy with exact spatial integration. It establishes no PDE convergence, other-record integral roots, global parameter-domain cover, dynamic path, biological observation alignment, energetic stability, or novelty. The historical unsafe point receipts remain withdrawn. This route uses corrected safe v2 uniform packets only. The arbitrary candidate d is read from exact differences of the two stored binary vectors; its historical proposal recipe is not a validated residual premise here.

Nine synthetic checks preceded the frozen protocol and twelve actual binding/arithmetic checks passed. One initially mistaken synthetic subtraction example was corrected before the protocol: .3 and .2 lie within the exact-subtraction range of Sterbenz's lemma. The original source and failure are preserved in engineering_revision_0. There was no failed or repeated physical calculation.
'''
(P/'PROOF.md').write_text(proof,encoding='utf-8')
wr(P/'terminal.json',{'status':'FINISHED','exit_code':0,'actual_command':'python -B -X utf8 crosscheck.py run','stdout_sha256':hs(P/'stdout.log'),'results_sha256':hs(P/'results.json'),'no_integration_matrix_or_solve':True,'tool_observed_wall_seconds':3.6312368})
checks={'source_current':r['source_sha256']==hs(P/'crosscheck.py'),'protocol_current':r['protocol_sha256']==hs(P/'protocol.json'),'all_inputs_current':all(hs(ROOT/k)==v for k,v in p['input_hashes'].items()),'all_actual_checks_true':all(v is True for v in r['checks'].values()),'all_synthetic_checks_true':all(v is True for v in rd(P/'selfchecks.json')['checks'].values()),'strict_selfmap_exact':Q(r['selfmap_lhs_exact'])<radius,'strict_contraction_exact':Q(r['Z_integral_exact'])<1,'error_exact':error==Q(r['eta_direct_exact'])/(1-Q(r['Z_integral_exact'])),'common_radius_exact':radius+Q(r['actual_shift_norm_exact'])==Q(r['common_radius_exact']),'new_root_outside_old_cube_exact':location,'Q_integral_roots_separation_positive':sep>0,'no_whole_goal_completion_claim':True}
assert all(v is True for v in checks.values())
out={'status':'PASS_DIRECT_UNIFORM_FORCE_ROOT_CROSSCHECK','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'source_sha256':hs(P/'crosscheck.py'),'protocol_sha256':hs(P/'protocol.json'),'results_sha256':hs(P/'results.json'),'proof_sha256':hs(P/'PROOF.md'),'location_consequences_sha256':hs(P/'location_consequences.json'),'terminal_sha256':hs(P/'terminal.json'),'finalizer_sha256':hs(Path(__file__)),'actual_status':r['status'],'strict_selfmap_pass':passed,'parameter_cell_sha256':r['parameter_cell_sha256'],'common_geometry_domain_id':r['common_geometry_domain_id'],'display':r['display'],'whole_physical_kernel_independently_reimplemented':False,'integration_and_local_proof_kernels_reused':True,'scope':r['scope']}
e=ROOT/'evidence/exact_uniform_direct_root_B_crosscheck.json';assert not e.exists();wr(e,out)
wr(P/'delivery.json',{'status':out['status'],'evidence_file':str(e.relative_to(ROOT)).replace('\\','/'),'evidence_sha256':hs(e),'checks':checks,'files_sha256':{x.name:hs(x)for x in P.iterdir()if x.is_file()and x.name!='delivery.json'}})
print(out['status'],len(checks),'checks','evidence',hs(e),'delivery',hs(P/'delivery.json'),'separation',float(sep))
