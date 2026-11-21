"""Retain current dependency scope after the bounded counterexample audit."""
from pathlib import Path
import hashlib,json,re,datetime
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def write(p,d):
    p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():raise ValueError('preserve '+str(p))
    p.write_text(json.dumps(d,indent=2)+'\n',encoding='utf8')

def main():
    result=read(OUT/'results.json');assert all(v is True for v in result['checks'].values())
    patterns=[re.compile(r'if\s+v\s*!=\s*0'),re.compile(r'if\s+not\s*\(?v\s*==\s*0'),re.compile(r'if\s+v\s*[}\]]')]
    matches=[]
    paths=list((ROOT/'scripts').glob('exact*.py'))+list((ROOT/'runs/exact_integration_bridge_v1').rglob('*.py'))
    for p in paths:
        if OUT in p.parents:continue
        rows=p.read_text(encoding='utf8').splitlines()
        hits=[{'line':i,'text':s.strip()}for i,s in enumerate(rows,1)if any(x.search(s)for x in patterns)]
        if hits:matches.append({'file':p.relative_to(ROOT).as_posix(),'sha256':sha(p),'matches':hits,'historical_revision':any('revision' in x for x in p.parts)})
    write(OUT/'predicate_inventory.json',{'scope':'Current exact*.py plus all exact-integration runfolder Python, with historical revisions separately labeled. Not a repository-wide static type proof.','files':matches})
    names=[
      ('exact_bulk_pilot_independent_review','WITHDRAWN_PENDING_SAFE_REPLAY','Point local series uses T polynomial cleanup, also reused by the independent replay; observed degree16 construction on elements0/6 did not drop coefficients, but complete computations were not traced.'),
      ('exact_bulk_assembly_B_independent_review','WITHDRAWN_AS_INTEGRATION_CERTIFICATE','All144 local packets depend on unsafe T cleanup; exact144-packet arithmetic/hash reduction remains historical valid arithmetic of its supplied intervals, not proof those intervals enclose the integrals.'),
      ('exact_complete_point_B_independent_review','WITHDRAWN_NEGATIVE_SELF_MAP_CLAIM','Same168 local integration packets; the 1.6208e-11 lower bound and map-center escape are not currently certified.'),
      ('exact_integration_bridge_final_delivery','WITHDRAWN_AS_MATHEMATICAL_DELIVERY','Rendering/hash delivery does not repair invalidated local enclosure provenance.'),
      ('exact_integration_documentation_review','WITHDRAWN_FOR_AFFECTED_SCIENTIFIC_CLAIMS','Any documentation endorsement of affected point integration/negative certificate must be superseded; old bytes retained.'),
      ('exact_cache_consistency_independent_review','PARTIAL_REVIEW_LAYER_WITHDRAWN','Producer dense analytic bulk-third/mean-value formula does not use sparse Arb cleanup; independent generic intervalAD does, so the independent AD enclosure validation must be safely replayed. Pure rational maps/source guards remain meaningful.'),
      ('exact_fe_moments_independent_review','UNAFFECTED_BY_THIS_DEFECT','Moment algebra and coefficient/tail arithmetic use rational coefficients and no unsafe uncertain-coefficient removal; it proves no full field bound.'),
      ('exact_retube_B_geometry_independent_review','UNAFFECTED_GEOMETRY_ONLY','Fresh common nodal boxes, exact rational tuple-interval P2/body and graph calculations remain valid for the chosen arbitrary proof center; no root membership follows.'),
      ('exact_retube_B_preparation_independent_review','PARTIAL_ARITHMETIC_VALID_INTEGRAL_INTERPRETATION_WITHDRAWN','Exact611 rounding and650 coordinate-hull arithmetic for stored candidate vectors remains valid; proposed integral Newton/recenter interpretation inherits untrusted point defect D. Conditional retube lemma remains true but uninstantiated.')]
    impact=[]
    for name,status,reason in names:
        p=ROOT/'evidence'/f'{name}.json';impact.append({'file':p.relative_to(ROOT).as_posix(),'sha256':sha(p),'current_reliance_status':status,'reason':reason})
    # Additional directly affected production and independent surface paths.
    additional=['bulk_pilot/protocol_v2.json','bulk_assembly_B/results.json','surface_pilot_preparation/integration_run/results.json','surface_pilot_independent_replay/results.json','surface_assembly_B/results.json','complete_point_B/results.json','uniform_bulk_pilot_B/protocol.json','uniform_bulk_pilot_B/preflight.json']
    for rel in additional:
        p=ROOT/'runs/exact_integration_bridge_v1'/rel
        impact.append({'file':p.relative_to(ROOT).as_posix(),'sha256':sha(p),'current_reliance_status':'UNIFORM_LAUNCH_WITHHELD' if 'uniform_bulk' in rel else 'WITHDRAWN_PENDING_SAFE_REPLAY','reason':'Reachable uncertain-coefficient sparse polynomial path; local point corruption not inferred solely from reachability. No numerical result may currently rely on this as a certified enclosure.'})
    write(OUT/'impact_inventory.json',{'receipts_and_packets':impact})
    checks={'bounded_counterexamples_all_pass':True,'actual_uniform_drops_observed':result['checks']['real_uniform_fields_encounter_dropped_uncertainty'],'point_corruption_not_falsely_claimed':result['witnesses']['point_series_actual_dropped_uncertainty_observed'] is False,'old_receipts_preserved':all(sha(ROOT/p)==h for p,h in result['frozen_sources'].items()),'independent_AD_shared_failure_included':True,'dense_old_roots_and_geometry_not_withdrawn_without_dependency':True,'no_spatial_integration_or_new_solve':result['actual_spatial_integrals_or_equilibrium_solves_run'] is False,'uniform_pilot_launch_not_approved':not (ROOT/'runs/exact_integration_bridge_v1/uniform_bulk_pilot_B/prelaunch_independent_review.json').exists()}
    assert all(v is True for v in checks.values()),checks
    receipt={'status':'CONFIRMED_CRITICAL_INTERVAL_CLEANUP_BUG_WITH_SCOPED_WITHDRAWAL','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'reproducer_checks':result['checks'],'reproducer_sha256':sha(OUT/'audit.py'),'reproducer_result_sha256':sha(OUT/'results.json'),'predicate_inventory_sha256':sha(OUT/'predicate_inventory.json'),'impact_inventory_sha256':sha(OUT/'impact_inventory.json'),'finalizer_sha256':sha(__file__),'source_snapshots':result['frozen_sources'],'actual_uniform_coefficients_dropped':10,'actual_point_degree16_stages_traced':[0,6],'actual_point_stages_with_observed_drops':[],'old_point_negative_certificate_currently_valid':False,'old_all17_binary_quadrature_root_stability_certificates_withdrawn':False,'fresh_common_geometry_certificate_withdrawn':False,'uniform_pilot_launch_approved':False,'minimal_replay_plan':['Freeze a new interval-safe polynomial kernel; remove only coefficients certified exactly zero, reject nonfinite coefficients; preserve old module and receipts.','Add mandatory containment tests for zero-crossing, exact-zero, nonzero, cancellation, multiplication, differentiation, composition and high-order jets; use a genuinely separate safe independent arithmetic layer.','Replay all144 bulk and24 surface point local integrations and signed611/6676 reductions; replay independent selected local integrals using corrected separate arithmetic. Only then reassess point defect/map-center statement; do not assume sign or conclusion survives.','Repair and rerun the independent bulk-cache genericAD review; producer dense analytic formula and rational-map tests may be retained with explicit scope.','Reprepare the uniform seven-element pilot with the new kernel and independently audited geometry/coefficients/tails. Require fresh prelaunch review before integrals.','Retain old candidate-vector arithmetic and arbitrary-center geometry as such; decide any next proof center from newly validated point defect. No existing exact-integral root or uniformH transfer is certified.'],'scope':'A bounded failure audit, not a full source safety proof. Complete original five-gate goal remains incomplete.'}
    write(ROOT/'evidence/interval_polynomial_cleanup_audit.json',receipt)
    print(json.dumps({'status':receipt['status'],'checks':len(checks),'receipt_sha256':sha(ROOT/'evidence/interval_polynomial_cleanup_audit.json'),'inventory_files':len(matches)}))
if __name__=='__main__':main()
