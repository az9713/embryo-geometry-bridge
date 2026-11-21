from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
root=Path('controlled_furrowing');out=root/'runs/axisymmetric'
records=[json.loads((out/f'warm_gain_{g}.json').read_text()) for g in [.2,.5,1.0,2.0,4.0]]
r=records[-1];x=np.asarray(r['current_nodes']);x0=np.asarray(r['reference_nodes']);t=np.asarray(r['triangles']);active=np.asarray(r['activation_elements'])
plt.style.use('dark_background');fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
ax=axes[0,0];poly=PolyCollection(x[t],facecolors=np.where(active,'#d88c57','#315b7e'),edgecolors='#7898aa',linewidths=.4,alpha=.8);ax.add_collection(poly)
for surface in ['initial_top_profile','initial_bottom_profile']:
 q=np.asarray(r[surface]);ax.plot(q[:,0],q[:,1],'--',color='#c0c0c0',lw=1)
ax.set(xlim=(0,1.25),ylim=(-.1,1.1),xlabel='R / H',ylabel='Z / H',title='Finite-height coupled body; g / mu = 4 (no exaggeration)');ax.set_aspect('equal')
ax.text(.45,.38,'Nonactivated tissue',fontsize=10);ax.text(.025,.65,'Active\nfibers',fontsize=9)
ax=axes[0,1];g=np.array([q['activation_strength'] for q in records]);ax.plot(g,[q['central_apical_displacement'] for q in records],'o-',label='Apical signed displacement');ax.plot(g,[q['central_basal_displacement'] for q in records],'s-',label='Basal signed displacement');ax.plot(g,[q['maximum_activated_column_shortening'] for q in records],'^-',label='Static material-column shortening');ax.axhline(0,color='grey',lw=.6);ax.set(xlabel='Activation gain / mu',ylabel='Displacement / H',title='Synthetic small circular patch; Kb / (mu H) = 0.2');ax.legend(fontsize=8)
ax=axes[1,0];mesh=[json.loads((out/q).read_text()) for q in ['base.json','mesh_12x6.json','mesh_16x8.json']];ax.plot([8,12,16],[q['maximum_activated_column_shortening'] for q in mesh],'o-',color='#edb65a');ax.set(xticks=[8,12,16],xticklabels=['8 x 4','12 x 6','16 x 8'],xlabel='Meridional mesh cells (radial x depth)',ylabel='Static shortening / H',title='Unresolved mesh sensitivity at g / mu = 0.2');ax.text(.04,.88,'Small force residual is not a continuum error bound',transform=ax.transAxes,fontsize=9,color='#edb65a')
ax=axes[1,1];base=mesh[0];stages=base['continuation'];gain=[q['activation_strength'] for q in stages];ax.plot(gain,[q['external_activation_work'] for q in stages],'o-',label='Gain-increment potential work');ax.plot(gain,[q['relaxation_energy_drop'] for q in stages],'s-',label='Relaxation potential drop');ax.set(xlabel='Activation gain / mu',ylabel='Energy / (mu H^3)',title='Declared external-potential bookkeeping');ax.legend(fontsize=8);ax.text(.04,.77,'No optical / ATP calibration; gauge dependent',transform=ax.transAxes,fontsize=9)
fig.suptitle('Idealized axisymmetric mechanics: verified discrete equations, unresolved biological and continuum gates',fontsize=12)
path=out/'axisymmetric_summary.png';fig.savefig(path,dpi=170);plt.close(fig)
receipt={'figure':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'inputs':{str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in list(out.glob('warm_gain_*.json'))+[out/'base.json',out/'mesh_12x6.json',out/'mesh_16x8.json']},'scope':'Synthetic data only; actual deformation no exaggeration in mesh panel; no observed large-patch outcomes read.'}
(root/'evidence/axisymmetric_visual.json').write_text(json.dumps(receipt,indent=2));print(json.dumps({'figure':str(path),'sha256':receipt['sha256']},indent=2))
