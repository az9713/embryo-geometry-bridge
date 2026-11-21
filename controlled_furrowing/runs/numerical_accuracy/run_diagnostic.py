"""Target-free fixed-strength compressibility diagnostic of the frozen v1 model."""
from pathlib import Path
import sys, json, hashlib, time
from dataclasses import asdict
import torch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from axisymmetric_mechanics import Parameters, build_model

torch.set_num_threads(1)
out = Path(__file__).resolve().parent
source = ROOT / 'scripts/axisymmetric_mechanics.py'
source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
base = Parameters(nr=12, nz=6, height=1., outer_radius=3.,
                  patch_radius=9.35512037/33.5,
                  activation_z_min=1-25/33.5, activation_z_max=1-5/33.5,
                  mu=1., basal_area_modulus=2., apical_area_modulus=.2)
protocol = {'scope': 'Target-free conditional numerical diagnostic, no data imports or fits',
            'frozen_source_sha256': source_hash, 'base_parameters': asdict(base),
            'fixed_activation_strength': 10., 'bulk_values': [5., 50., 500.],
            'continuation_steps': 10, 'maxiter': 200,
            'interpretation': 'Three compressibilities on one mesh do not establish volumetric locking; the full cross-mesh test and stable mixed reference remain required.'}
(out/'protocol.json').write_text(json.dumps(protocol, indent=2)+'\n')
results=[]
for bulk in protocol['bulk_values']:
    from dataclasses import replace
    started=time.monotonic()
    m=build_model(replace(base, bulk=bulk))
    result=m.solve(protocol['fixed_activation_strength'], continuation_steps=10, maxiter=200)
    record={'source_sha256': source_hash, 'protocol_sha256': hashlib.sha256((out/'protocol.json').read_bytes()).hexdigest(),
            'seconds':time.monotonic()-started, 'result':result}
    target=out/f'bulk_{bulk:g}.json'
    target.write_text(json.dumps(record,indent=2)+'\n')
    summary={k:result[k] for k in ('converged','minimum_J','maximum_J','relative_volume_change','scaled_free_force_linf','maximum_activated_column_shortening','central_apical_displacement','central_basal_displacement')}
    summary.update({'bulk':bulk,'seconds':record['seconds'],'artifact':target.name,'sha256':hashlib.sha256(target.read_bytes()).hexdigest()})
    results.append(summary)
    print(json.dumps(summary), flush=True)
(out/'summary.json').write_text(json.dumps({'protocol':protocol,'results':results},indent=2)+'\n')
