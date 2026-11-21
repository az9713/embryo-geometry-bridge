from pathlib import Path
from fractions import Fraction as Q
import sys,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from run_activated_mesh_study import QuadraticGraph,compare_graphs
from annulus_observer import serial
root=Path(__file__).resolve().parents[2];finepath=root/'runs/sparse_quadratic/24x12_lambda5000/gain_10.0/final.json';coarsepath=root/'runs/joint_refinement/final.json'
fine=json.loads(finepath.read_text());coarse=json.loads(coarsepath.read_text())
if not fine['converged'] or not fine['reached_target'] or not coarse['converged']:raise RuntimeError('both comparison states must verify')
a=dict(fine['parameters']);b=dict(coarse['parameters'])
for name in ['nr','nz']:a.pop(name);b.pop(name)
if a!=b or fine['integration_policy']!=coarse['integration_policy'] or fine['activation_strength']!=coarse['activation_strength']:raise RuntimeError('physical/integration/gain mismatch')
result={'scope':'MatchedhighQ24x12versus16x8discretefields; noPDEerrorbound, localbranchidentitynotproved','gain':10.,'lambda':5000.,'physical_height_um':33.5,'integration_policy':fine['integration_policy'],'source_fields':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [finepath,coarsepath]},'shortening_fine_um':fine['maximum_activated_column_shortening']*33.5,'shortening_coarse_um':coarse['maximum_activated_column_shortening']*33.5,'central_apical_difference_um':(fine['central_apical_displacement']-coarse['central_apical_displacement'])*33.5,'central_basal_difference_um':(fine['central_basal_displacement']-coarse['central_basal_displacement'])*33.5}
for surface in ['top','bottom']:
 key='boundary_'+surface+'_quadratic_controls';result[surface]=compare_graphs(QuadraticGraph(fine[key]),QuadraticGraph(coarse[key]))
(root/'evidence/sparse_quadratic_24_comparison.json').write_text(json.dumps(serial(result),indent=2))
plt.style.use('dark_background');fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for j,(surface,label) in enumerate([('top','Apical'),('bottom','Basal')]):
 q=result[surface]['sampled_eulerian_curves'];r=[x['radius_um'] for x in q]
 axs[0,j].plot(r,[x['coarse_z_um'] for x in q],'--',color='#dbae66',label='16 x 8, high Q');axs[0,j].plot(r,[x['fine_z_um'] for x in q],'-',color='#7bc7b4',label='24 x 12, high Q');axs[0,j].set(xlabel='Current physical radius (micrometres)',ylabel='Height (micrometres)',title=label+' full Eulerian graph');axs[0,j].legend(fontsize=8)
 axs[1,j].plot(r,[(x['difference_lower_um']+x['difference_upper_um'])/2 for x in q],color='#b3a0df');axs[1,j].axhline(0,color='grey',lw=.5);bound=result[surface]['supremum_absolute_discrete_height_difference_um'];axs[1,j].set(xlabel='Current physical radius (micrometres)',ylabel='24 minus16 height (micrometres)',title=f'Stored-curve sup: [{float(bound.lo):.4f}, {float(bound.hi):.4f}] micrometres')
fig.suptitle('Same high-quadrature constitutive model and gain10; mesh/branch sensitivity, not a PDE error estimate',fontsize=10)
image=root/'runs/sparse_quadratic/24_vs16_highQ.png';fig.savefig(image,dpi=170,bbox_inches='tight');plt.close(fig)
(root/'evidence/sparse_quadratic_24_comparison_visual.json').write_text(json.dumps({'image':str(image.relative_to(root)),'sha256':hashlib.sha256(image.read_bytes()).hexdigest(),'comparison_sha256':hashlib.sha256((root/'evidence/sparse_quadratic_24_comparison.json').read_bytes()).hexdigest()},indent=2))
print(json.dumps({'shortening_fine_um':result['shortening_fine_um'],'shortening_coarse_um':result['shortening_coarse_um'],'apical_sup_um':[float(result['top']['supremum_absolute_discrete_height_difference_um'].lo),float(result['top']['supremum_absolute_discrete_height_difference_um'].hi)],'basal_sup_um':[float(result['bottom']['supremum_absolute_discrete_height_difference_um'].lo),float(result['bottom']['supremum_absolute_discrete_height_difference_um'].hi)]},indent=2))
