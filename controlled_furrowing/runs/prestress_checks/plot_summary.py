from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('controlled_furrowing');out=root/'runs/prestress_checks'
plt.style.use('dark_background');fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
colors=['#d9a65b','#7bc7bb','#b49fdf'];tensions=[.2,2.,20.];inputs={}
for tension,color in zip(tensions,colors):
 records=[]
 for gain in [.2,1.,4.]:
  path=out/f'Tb_{tension}_gain_{gain}.json';records.append(json.loads(path.read_text()));inputs[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
 r=records[1]
 for ax,key,key0,label in [(axes[0,0],'top_profile','initial_top_profile','Apical'),(axes[0,1],'bottom_profile','initial_bottom_profile','Basal')]:
  q=np.asarray(r[key]);q0=np.asarray(r[key0]);ax.plot(q[:,0],q[:,1]-q0[:,1],'o-',color=color,label=f'Tb / (mu H) = {tension}')
  ax.set(xlim=(0,1.5),xlabel='Current radial coordinate R / H (zoom)',ylabel='Vertical displacement / H',title=f'{label} full FE surface sampled at boundary nodes; gain 1')
  ax.axhline(0,color='grey',lw=.5)
 axes[1,1].plot([q['activation_strength'] for q in records],[q['maximum_activated_column_shortening'] for q in records],'o-',color=color,label=f'Tb = {tension}')
 axes[1,0].bar(np.log10(tension)-.10,abs(r['central_apical_displacement']),width=.18,color=color,alpha=1)
 axes[1,0].bar(np.log10(tension)+.10,abs(r['central_basal_displacement']),width=.18,color=color,alpha=.4)
axes[0,0].legend(fontsize=8);axes[1,0].set(xticks=np.log10(tensions),xticklabels=['0.2','2','20'],xlabel='Basal surface tension / (mu H)',ylabel='Absolute central vertical displacement / H',title='Gain 1: left opaque = apical; right translucent = basal')
axes[1,1].set(xlabel='Activation gain / mu',ylabel='Static maximum column shortening / H',title='Changing basal tension changes the full nonlinear response');axes[1,1].legend(fontsize=8)
fig.suptitle('Synthetic flat-reference tension revision: Kb=2, Ka=0.2, Ta=0.2; no biological outcomes fitted',fontsize=12)
path=out/'prestress_summary.png';fig.savefig(path,dpi=170);plt.close(fig)
receipt={'image':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'inputs':inputs,'scope':'Synthetic coupled equilibria only, no measured target reads or deformation magnification; radial surface panels zoomed.'}
(root/'evidence/prestress_visual.json').write_text(json.dumps(receipt,indent=2));print(json.dumps({'image':str(path),'sha256':receipt['sha256']},indent=2))
