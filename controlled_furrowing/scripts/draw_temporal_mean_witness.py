"""Original exact-summary path illustration, arbitrary time units."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
data=json.loads((ROOT/'evidence/temporal_mean_witness.json').read_text())
row=next(r for r in data['records'] if r['record_id'].endswith('_G'))
M,V,A,B=(float(Q(row['recorded'][k])) for k in ('M','V','A','B'))
t=np.linspace(0,2,401);a=np.where(t<=1,t*A,(2-t)*A+(t-1)*(V-B))
b=np.where(t<=1,t*(V-A),V-a);v=np.minimum(t,1)*V;m=np.minimum(t,1)*M
plt.style.use('dark_background');fig,axes=plt.subplots(1,2,figsize=(13,4.7))
for values,label,color in [(a,'Apical regional mean','#ffbc69'),(b,'Basal regional mean','#80d4bc'),
                           (v,'Height regional mean','#c89df0'),(m,'Spatial height maximum','#67c7ef')]:
    axes[0].plot(t,values,label=label,color=color,linewidth=2)
axes[0].axvspan(1,2,color='white',alpha=.05);axes[0].axhline(0,color='white',alpha=.3)
axes[0].set(xlabel='Arbitrary path parameter t (not seconds)',ylabel='Signed response (µm)',title='Instantaneous apical + basal = height mean')
axes[0].legend(fontsize=9,loc='lower right');axes[0].grid(alpha=.2)
values=[A,B,A+B,V,M]
axes[1].bar(range(5),values,color=['#ffbc69','#80d4bc','#ef9494','#c89df0','#67c7ef'])
axes[1].set_xticks(range(5),['Peak\napical mean','Peak\nbasal mean','Sum of\ncomponent peaks','Peak\nheight mean','Peak spatial\nheight maximum'],fontsize=9)
for i,value in enumerate(values):axes[1].text(i,value+.2,f'{value:.3f}',ha='center',fontsize=10)
axes[1].set(ylabel='Summary response (µm)',title='Four recorded summaries reproduced exactly',ylim=(0,19))
axes[1].grid(axis='y',alpha=.2)
fig.suptitle('Canonical record G: a volume-preserving kinematic path, not a dynamical tissue model\nV is peak-of-mean; mean-of-local-peaks is a different observation contract',fontsize=12)
fig.tight_layout();path=ROOT/'docs/assets/temporal_mean_witness.png'
fig.savefig(path,dpi=170,bbox_inches='tight');plt.close(fig)
receipt={'status':'GENERATED_TEMPORAL_KINEMATIC_VISUAL','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'image_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
         'source_receipt_sha256':hashlib.sha256((ROOT/'evidence/temporal_mean_witness.json').read_bytes()).hexdigest(),
         'scope':'Original synthetic path using source summary numbers; not measured time traces, force-balanced dynamics or local-peak-then-mean matching.'}
(ROOT/'evidence/temporal_mean_witness_visual.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'image':str(path)}))
