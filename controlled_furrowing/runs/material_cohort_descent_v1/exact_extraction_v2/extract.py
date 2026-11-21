"""Separate exact active-system extraction for the frozen cohort affine LPs."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[3];BASE=Path(__file__).resolve().parent.parent;OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
from bounded_affine_minimax import solve_minimax,verify_certificate
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
def save(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
if read(BASE/'job.json')['status']!='FINISHED':raise RuntimeError('Complete cohort job required; no partial candidate extraction')
path=OUT/'protocol.json'
if path.exists():raise RuntimeError('Extraction protocol exists; preserve previous computation instead of overwrite')
helper=ROOT/'scripts/bounded_affine_minimax.py'
protocol={'version':2,'stage':'Hardened-helper exact active-system extraction for SAME rounded affine problems',
    'helper_sha256':sha(helper),'generator_sha256':sha(__file__),
    'cohort_protocol_sha256':sha(BASE/'protocol.json'),'cohort_results_sha256':sha(BASE/'results.json'),
    'coefficient_rule':'Use saved exact Fraction coefficients/intercepts/bounds verbatim; no decimal/binary re-rounding',
    'new_equilibrium_solves':0,'scope':'Global exact affine LP certificate, not nonlinear model fit or mechanical derivative validation'}
save(path,protocol)
result=read(BASE/'results.json');cohort=read(BASE/'protocol.json');candidates=[]
for c in result['candidates']:
    p=c['affine_LP'];r=list(map(Q,p['affine_residual_intercepts_exact']))
    J=[list(map(Q,row)) for row in p['affine_derivative_coefficients_exact']]
    box=[tuple(map(Q,pair)) for pair in p['trust_box_exact']]
    exact=solve_minimax(r,J,box)
    if exact['status']!='CERTIFIED_DEFINED_AFFINE_LP':raise RuntimeError('Root helper could not certify '+c['variant'])
    checks=verify_certificate(r,J,box,exact['certificate'])
    lo,hi=map(Q,p['affine_optimum_interval_exact'])
    if not all(checks.values()) or not lo<=Q(exact['loss_exact'])<=hi:raise RuntimeError('Exact extraction inconsistent with producer weak-dual interval')
    step=exact['step_display'];base=cohort['parameters'];import math
    parameters={'apical_area_modulus':base['apical_area_modulus']*math.exp(step[0]),
        'basal_area_modulus':base['basal_area_modulus']*math.exp(step[1]),
        'apical_surface_tension':base['apical_surface_tension']+step[2],
        'basal_surface_tension':base['basal_surface_tension']*math.exp(step[3]),
        'common_gain':float(Q(cohort['baseline_gain_exact']))*math.exp(step[4]),
        'bulk_fixed':5000}
    candidates.append({'variant':c['variant'],'exact_affine_result':exact,
        'producer_optimum_interval_contains_exact_loss':True,
        'baseline_loss_um':p['baseline_loss_um'],'proposed_material_parameters_display':parameters,
        'no_new_equilibrium_verification':True})
if sha(helper)!=protocol['helper_sha256'] or sha(BASE/'results.json')!=protocol['cohort_results_sha256']:raise RuntimeError('Inputs changed during extraction')
final={'status':'PASS_ALL_FOUR_EXACT_AFFINE_CANDIDATES','protocol_sha256':sha(path),'candidates':candidates,
    'scope':protocol['scope'],'new_equilibrium_solves':0,'original_five_gate_goal_achieved':False}
save(OUT/'exact_candidates.json',final)
print(json.dumps({'status':final['status'],'sha256':sha(OUT/'exact_candidates.json'),
    'candidates':[{'variant':c['variant'],'baseline_loss':c['baseline_loss_um'],'affine_loss':c['exact_affine_result']['loss_display'],'step':c['exact_affine_result']['step_display']} for c in candidates]}),flush=True)
