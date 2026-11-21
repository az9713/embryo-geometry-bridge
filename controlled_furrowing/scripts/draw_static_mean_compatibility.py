"""Original source-summary and kinematic-witness visualization."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
data=json.loads((ROOT/'evidence/static_mean_compatibility.json').read_text())
iso=json.loads((ROOT/'evidence/static_isochoric_construction.json').read_text())
rows=data['canonical_rows'];sums=[];means=[];errors=[];labels=[]
for row in rows:
    raw=row['original'];means.append(float(Q(raw['region_averaged_shortening_um'])))
    sums.append(float(Q(raw['apical_max_shortening_um'])+Q(raw['basal_max_shortening_um'])))
    errors.append(float(Q(row['correction']['error'])));labels.append(row['record_id'].split('_')[-1])
selected=next(r for r in iso['constructions'] if r['pairing']=='canonical' and r['record_id'].endswith('_G') and r['weighting']=='disk')
H,a,Ro=(float(Q(selected[k])) for k in ('height_um','ROI_radius_um','outer_radius_um'))
M,_,negC=map(lambda x:float(Q(x)),selected['core_S_coefficients_normalized_radius']);C=-negC
A=float(Q(selected['core_apical_downward_constant']));e,b,c=map(lambda x:float(Q(x)),selected['annular_S_coefficients_normalized_annulus_radius'])
r=np.linspace(0,Ro,1501);t=np.maximum(0,(r-a)/(Ro-a))
S=np.where(r<=a,M-C*(r/a)**2,e+b*t+c*t*t)
apical=np.where(r<=a,A,A*(1-t));basal=S-apical;height=H-S;top=H-apical
plt.style.use('dark_background');fig,axes=plt.subplots(1,3,figsize=(17,5))
axes[0].scatter(means,sums,color='#67c7ef',s=35)
axes[0].plot([0,18],[0,18],'--',color='#aab5c7',label='Static mean identity')
i=labels.index('G');axes[0].scatter([means[i]],[sums[i]],color='#ffbc69',s=70)
axes[0].annotate('G: additive defect 10.4253 µm',(means[i],sums[i]),xytext=(.5,17),arrowprops={'arrowstyle':'->','color':'white'},fontsize=9)
axes[0].set(xlabel='Reported regional mean (µm)',ylabel='Reported apical + basal entries (µm)',title='Identity applies only to one static interpretation',xlim=(0,18),ylim=(0,19))
axes[0].legend(fontsize=9);axes[0].grid(alpha=.2)
colors=['#ffbc69' if label=='G' else '#67c7ef' for label in labels]
axes[1].bar(labels,errors,color=colors)
axes[1].axhline(float(Q(data['static_regional_mean_canonical_floor_exact'])),color='#ffbc69',linestyle='--')
axes[1].set(xlabel='Canonical source-column record',ylabel='Necessary common channel discrepancy (µm)',title='Sharp kinematic floor: 3.47510056 µm',ylim=(0,4.1))
axes[1].grid(axis='y',alpha=.2)
axes[2].fill_between(r,basal,top,color='#24485d',alpha=.9)
axes[2].plot(r,top,color='#67c7ef',label='Apical surface')
axes[2].plot(r,basal,color='#e7aa65',label='Basal surface')
for fraction in (.25,.5,.75):axes[2].plot(r,basal+fraction*height,color='white',alpha=.25,linewidth=.7)
axes[2].axvline(a,color='white',linestyle=':',label='ROI boundary')
axes[2].axhline(0,color='white',alpha=.2);axes[2].axhline(H,color='white',alpha=.2)
axes[2].set(xlabel='Current radius (µm)',ylabel='Vertical coordinate (µm)',title='Volume-preserving witness for corrected G',ylim=(-7,38))
axes[2].legend(fontsize=8,loc='lower right');axes[2].grid(alpha=.2)
fig.suptitle('Conditional static regional means: exact consistency bound and geometric sharpness\nNot a noise estimate, temporal-peak identity or common force-balanced model',fontsize=13)
fig.tight_layout();target=ROOT/'docs/assets/static_mean_compatibility.png'
fig.savefig(target,dpi=170,bbox_inches='tight');plt.close(fig)
receipt={'status':'GENERATED_SOURCE_AND_KINEMATIC_VISUAL','image_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
         'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'input_hashes':{name:hashlib.sha256((ROOT/'evidence'/name).read_bytes()).hexdigest() for name in ('static_mean_compatibility.json','static_isochoric_construction.json')},
         'scope':'Published numeric summaries and an original synthetic volume-preserving witness; no source-image pixels or force-balanced biological shape.'}
(ROOT/'evidence/static_mean_compatibility_visual.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'image':str(target)}))
