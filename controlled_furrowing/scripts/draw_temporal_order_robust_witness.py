"""Visualize exact kinematic receipts; no measured spatial traces invented."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
source=ROOT/'evidence/temporal_mean_order_robust_witness.json'
data=json.loads(source.read_text())
if not all(data['checks'].values()):raise RuntimeError('classification failed')
row=next(r for r in data['rows'] if r['record_id'].endswith('_G'))
boundary=data['nonattainment_boundary_rows'][0]
x=np.linspace(0,1,1001)
plt.style.use('dark_background')
fig,axes=plt.subplots(1,3,figsize=(15,4.5))
for weight,color in [('line','#efb56a'),('disk','#80d4bc')]:
    item=row['profiles'][weight];M=float(Q(item['maximum_exact']));C=float(Q(item['nonnegative_baseline_exact']))
    values=C+(M-C)*(1-x*x)**item['degree_parameter']
    axes[0].plot(x,values,label=f"{weight}, p={item['degree_parameter']}",color=color,lw=2)
    entries=boundary['approximating_sequences'][weight]
    axes[2].loglog([v['p'] for v in entries],[v['mean_error_display'] for v in entries],'o-',color=color,label=weight)
for p in (8,32,128,512):
    axes[1].plot(x,float(Q(boundary['maximum_exact']))*(1-x*x)**p,label=f'p={p}')
axes[0].set(title='Record G: same maximum and mean',xlabel='Current radius / ROI radius',ylabel='Constructed shortening (µm)')
axes[1].set(title='Record O: fixed positive maximum',xlabel='Current radius / ROI radius',ylabel='Constructed shortening (µm)',xlim=(0,.65))
axes[2].set(title='O: mean tends to zero, never equals zero',xlabel='Concentration parameter p',ylabel='Positive mean error (µm)')
for ax in axes:ax.legend(fontsize=9);ax.grid(alpha=.2)
fig.suptitle('16 exact matches for both height orderings; record O: mean-of-local-peaks is nonattainable\nSynthetic volume-preserving kinematics; no common force law or measured spatial profiles',fontsize=12)
fig.tight_layout();path=ROOT/'docs/assets/temporal_order_robust_witness.png'
fig.savefig(path,dpi=170,bbox_inches='tight');plt.close(fig)
receipt={'status':'GENERATED_FROM_EXACT_CLASSIFICATION','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'source_receipt_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
         'image_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
         'scope':'Synthetic kinematics; height-summary ordering only. Components remain separately peaked regional means. No common mechanics or biological length cutoff.'}
(ROOT/'evidence/temporal_order_robust_visual.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'image':str(path)}))
