"""Original analytic visualization: material labels versus current radius."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
plt.style.use('dark_background')
t=np.linspace(0,1,501);radius=np.linspace(0,2,501)
ta=(np.sqrt(1+4*radius)-1)/2
eulerian=1+ta-radius/2
fig,axes=plt.subplots(1,3,figsize=(16,4.7))
axes[0].plot(t,t+t*t,label='Apical R(t)',color='#55c4ef')
axes[0].plot(t,2*t,label='Basal R(t)',color='#ffbd69')
axes[0].set(xlabel='Material parameter t',ylabel='Current radius R',title='Equal material labels reach different radii')
axes[0].legend();axes[0].grid(alpha=.2)
axes[1].plot(radius,eulerian,label='Same current radius: S(R)',color='#80d4bc')
axes[1].axhline(1,color='#b0b7c3',linestyle='--',label='Material-paired shortening = 1')
axes[1].scatter([.75],[9/8],color='#ffbd69',zorder=5)
axes[1].annotate('Exact maximum 9/8',(.75,9/8),xytext=(.92,1.115),arrowprops={'arrowstyle':'->','color':'white'})
axes[1].set(xlabel='Current radius R',ylabel='Vertical shortening',title='Pairing changes the measured quantity',ylim=(.985,1.145))
axes[1].legend(loc='lower right');axes[1].grid(alpha=.2)
values=[1,43/40,13/12,9/8,2]
axes[2].bar(range(5),values,color=['#b0b7c3','#55c4ef','#80d4bc','#ffbd69','#ef9191'])
axes[2].set_xticks(range(5),['Material\npair','Disk\nmean','Line\nmean','Eulerian\nmax','Component\nmax sum'],fontsize=10)
for i,(v,label) in enumerate(zip(values,['1','43/40','13/12','9/8','2'])):
    axes[2].text(i,v+.035,label,ha='center')
axes[2].set(ylabel='Shortening',title='Five different summaries of one geometry',ylim=(0,2.25))
axes[2].grid(axis='y',alpha=.2)
fig.suptitle('Synthetic exact example: apical (R,Z)=(t+t²,10−t); basal (R,Z)=(2t,1−t)\nDeclared flat reference height 10; no biological measurement or fitted parameters',fontsize=13)
fig.tight_layout()
target=ROOT/'docs/assets/endpoint_counterexample.png'
fig.savefig(target,dpi=170,bbox_inches='tight');plt.close(fig)
receipt={'status':'GENERATED_ANALYTIC_VISUAL','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'image_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
         'analytic_values_exact':['1','43/40','13/12','9/8','2'],
         'scope':'Original synthetic analytic geometry only; plot sampling does not supply the exact proof.'}
(ROOT/'evidence/endpoint_counterexample_visual.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'image':str(target)}))
