"""Plot actual bound receipts without fabricating a parameter response curve."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
ep=ROOT/'evidence/energy_gain_cell.json';gp=ROOT/'evidence/energy_gain_geometry.json'
e=json.loads(ep.read_text(encoding='utf-8'));g=json.loads(gp.read_text(encoding='utf-8'))
if e['status']!='CERTIFIED_UNIFORM_ENERGY_GAIN_CELL' or g['status']!='CERTIFIED_ENERGY_GAIN_FAMILY_GEOMETRY':raise RuntimeError('actual certificates required')
plans=e['completed_plan_bounds'];x=[0.]+[float(Q(p['time_interval'][1])) for p in plans[:-1]]
y=[float(Q(p['total_transformed_Hessian_variation_bound']['upper_exact'])) for p in plans[:-1]]
rho=float(Q(e['metric_error_radius_exact']));eta=float(Q(e['uniform_predictor_eta']['upper_exact']));Z=float(Q(e['uniform_contraction_Z']['upper_exact']))
lo=float(Q(g['geometry']['J_bound']['lower_exact']));hi=float(Q(g['geometry']['J_bound']['upper_exact']))
plt.style.use('dark_background');fig,axes=plt.subplots(1,3,figsize=(15,4.7))
axes[0].stairs(y,x,color='#81d8b5',lw=2,fill=True,alpha=.5)
axes[0].set(xlabel='Taylor path fraction t (not physical time)',ylabel='Certified Hessian-variation upper',title='Eight exact path subintervals');axes[0].grid(alpha=.2)
axes[1].bar([0],[eta/rho],color='#72c8ed',label='Uniform residual η / ρ')
axes[1].bar([0],[Z],bottom=[eta/rho],color='#efbb71',label='Contraction Z')
axes[1].axhline(1,color='white',ls=':',label='Self-map threshold')
axes[1].set(xticks=[0],xticklabels=['η/ρ + Z'],ylabel='Fraction of allowed radius',ylim=(0,1.15),title='Strict self-map condition passes');axes[1].legend(fontsize=9,loc='upper left');axes[1].grid(axis='y',alpha=.2)
axes[2].axhspan(.99,1.01,color='white',alpha=.07,label='One-percent band')
axes[2].fill_between([.999,1.001],[lo]*2,[hi]*2,color='#a997e8',alpha=.5,label='Common whole-body enclosure')
axes[2].axhline(1,color='white',ls=':');axes[2].set(xlabel='Active strength at diagnostic unit input',ylabel='Jacobian J',ylim=(.988,1.012),title='All144 triangles, all admitted gains');axes[2].legend(fontsize=9,loc='lower left');axes[2].grid(alpha=.2)
axes[2].set_xticks([.999,.9995,1.,1.0005,1.001])
axes[2].set_xticklabels(['0.9990','0.9995','1.0000','1.0005','1.0010'])
fig.suptitle('One fixed weak-load model: a uniform gain family with validated local roots and whole-body geometry\nMathematical error enclosures, not measured tissue responses or full material inference',fontsize=12)
fig.tight_layout();path=ROOT/'docs/assets/energy_gain_cell_validation.png';fig.savefig(path,dpi=170,bbox_inches='tight');plt.close(fig)
receipt={'status':'GENERATED_FROM_ACTUAL_UNIFORM_CERTIFICATES','energy_sha256':hashlib.sha256(ep.read_bytes()).hexdigest(),
         'geometry_sha256':hashlib.sha256(gp.read_bytes()).hexdigest(),'image_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
         'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'scope':'Taylor fraction is not physical time; J band is a common enclosure, not a sampled physical curve; fixed weak benchmark only.'}
(ROOT/'evidence/energy_gain_cell_visual.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':receipt['status'],'image':str(path)}))
