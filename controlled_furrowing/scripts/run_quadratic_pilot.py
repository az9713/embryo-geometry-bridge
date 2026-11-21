"""Target-free accuracy pilot on the actual small-patch protocol geometry."""
from dataclasses import asdict
from pathlib import Path
import hashlib
import json
import time
import torch
from axisymmetric_quadratic import Parameters,build_model

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/quadratic_pilot';OUT.mkdir(parents=True,exist_ok=True)
torch.set_num_threads(1)
H=33.5
protocol={'scope':'Target-free finite-height numerical pilot, not fit or held-out evaluation',
          'height_um':H,'nr':12,'nz':6,'penalties':[500.,5000.],
          'basal_tension':2.,'apical_tension':0.,'strengths':[0.,1.,5.,10.],
          'activation_radius_um':18.71024074/2,'activation_depth_below_top_um':[5.,25.],
          'source_sha256':{name:hashlib.sha256((ROOT/'scripts'/name).read_bytes()).hexdigest() for name in ['axisymmetric_mechanics.py','axisymmetric_prestress.py','axisymmetric_mixed.py','axisymmetric_quadratic.py']},
          'interpretation':'Report stationarity, sampled J and sufficient Bernstein enclosure separately; no continuum or incompressibility claim from mean volume.'}
(OUT/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8')
summaries=[];started=time.monotonic()
for penalty in protocol['penalties']:
    p=Parameters(nr=12,nz=6,bulk=penalty,patch_radius=protocol['activation_radius_um']/H,
                 activation_z_min=(H-25)/H,activation_z_max=(H-5)/H,
                 basal_surface_tension=2.,apical_surface_tension=0.,
                 basal_area_modulus=2.,apical_area_modulus=.2,quadrature=4,surface_quadrature=8)
    model=build_model(p);previous=None;previous_strength=0.
    for strength in protocol['strengths']:
        steps=max(1,round(strength-previous_strength))
        result=model.solve(strength,initial=previous,continuation_start_strength=previous_strength,continuation_steps=steps,maxiter=100)
        path=OUT/f'penalty_{penalty:g}_strength_{strength:g}.json'
        path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        summary={key:result.get(key) for key in ['converged','stationary_discrete_equations','scaled_free_force_linf','minimum_J','maximum_J','sampled_minimum_J','sampled_maximum_J','relative_volume_change','maximum_activated_column_shortening','central_apical_displacement','central_basal_displacement']}
        summary.update(penalty=penalty,strength=strength,path=path.relative_to(ROOT).as_posix(),seconds=time.monotonic()-started)
        summaries.append(summary);print(json.dumps(summary),flush=True)
        previous=result['free_displacement'];previous_strength=strength
(OUT/'summary.json').write_text(json.dumps({'protocol':protocol,'cases':summaries},indent=2)+'\n',encoding='utf-8')
