"""Create a separate calibration-driven revision; preserve all v1 source bytes."""
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]
old=ROOT/'scripts/spatial_transfer.py'
text=old.read_text(encoding='utf-8')
old_hash=hashlib.sha256(old.read_bytes()).hexdigest()


def replace(old,new,count=None):
    global text
    actual=text.count(old)
    if not actual or (count is not None and actual!=count):
        raise RuntimeError(f'Unexpected source replacement count {actual}: {old}')
    text=text.replace(old,new)


replace('from axisymmetric_mechanics import Parameters,build_model','from axisymmetric_prestress import Parameters,build_model',1)
replace('SPATIAL-PROTOCOL.json','SPATIAL-TENSION-PROTOCOL.json')
replace('scripts/spatial_transfer.py','scripts/spatial_tension_transfer.py')
replace("'scripts/axisymmetric_mechanics.py'","'scripts/axisymmetric_mechanics.py','scripts/axisymmetric_prestress.py'",1)
replace('runs/spatial_transfer/','runs/spatial_tension/')
replace('evidence/spatial_transfer_','evidence/spatial_tension_')
replace("return f'basal_{kb:g}'","return f'tension_{kb:g}'",1)
replace("protocol['material_scenarios']['basal_area_modulus_over_mu_height']","protocol['material_scenarios']['basal_surface_tension_over_mu_height']",1)
replace('basal_area_modulus=kb,apical_area_modulus=',"basal_area_modulus=material['basal_area_modulus_over_mu_height'],basal_surface_tension=kb,apical_surface_tension=0.,apical_area_modulus=",1)
replace("'basal_modulus':kb","'basal_surface_tension':kb",1)
replace("'--basal'","'--tension'",1)
replace('args.basal','args.tension',1)
text='"""Separate v2 tension revision generated from retained v1; see generation receipt."""\n'+text
# Keep the future import at the beginning after a single module docstring.
text=text.replace('from __future__ import annotations\n','',1)
destination=ROOT/'scripts/spatial_tension_transfer.py'
if destination.exists():raise RuntimeError('Do not overwrite an existing experiment revision')
destination.write_text(text,encoding='utf-8')
protocol=json.loads((ROOT/'SPATIAL-PROTOCOL.json').read_text(encoding='utf-8'))
protocol['revision_reason']='Small-region-only audit: all three v1 area-stiffness scenarios predict basal-dominant motion in 17/17 records, source apical-dominant 16/17. Baseline membrane tension supplies missing quadratic normal stiffness. No v2 larger-region outcomes have been calculated or used to choose this family.'
protocol['status']='PRE-CALIBRATION separate prestress revision; v1 retained unchanged'
protocol['material_scenarios']['basal_area_modulus_over_mu_height']=2.
protocol['material_scenarios']['basal_surface_tension_over_mu_height']=[0.,.2,2.,20.]
protocol['material_scenarios']['apical_surface_tension_over_mu_height']=0.
protocol['material_scenarios']['scope']='Exploratory fixed area moduli and basal prestress scenarios; no measured/exhaustive parameter bounds. All retained. Strong compression and mesh sensitivity remain unpassed adequacy requirements.'
protocol['calibration']['gain_grid']={'minimum':0.,'maximum':40.,'step':.05}
protocol['calibration']['initial_response_strength_nodes']=[i/2 for i in range(33)]+[float(i) for i in range(18,65,2)]
protocol['calibration']['secondary_check']='Report separate small-region apical/basal component maxima and directionality; never substitute their sum for maximum height change. Static-versus-transient observation mismatch remains explicit.'
protocol['original_scope_gaps'].append('V1 local J reached about0.8 under fitted loads; prestress does not itself fix compressibility, locking or continuum error. V2 must report them.')
path=ROOT/'SPATIAL-TENSION-PROTOCOL.json'
path.write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8')
receipt={'status':'CREATED_BEFORE_V2_CALIBRATION','v1_driver_sha256':old_hash,
         'v2_driver_sha256':hashlib.sha256(destination.read_bytes()).hexdigest(),
         'v2_protocol_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
         'original_model_sha256':hashlib.sha256((ROOT/'scripts/axisymmetric_mechanics.py').read_bytes()).hexdigest(),
         'prestress_model_sha256':hashlib.sha256((ROOT/'scripts/axisymmetric_prestress.py').read_bytes()).hexdigest(),
         'scope':'Source and protocol revision driven by small-region failure; no v1 overwrite or target-conditioned parameter selection'}
(ROOT/'evidence/tension_experiment_generation.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
print(json.dumps(receipt))
