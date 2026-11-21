"""Receipt-derived local-affine diagnostics, no new mechanical solves."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parent
read=lambda p:json.loads(p.read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
if read(OUT/'job.json')['status']!='FINISHED':raise RuntimeError('Complete cohort required')
path=OUT/'exact_extraction_v2/exact_candidates.json';result=read(path)
if result['status']!='PASS_ALL_FOUR_EXACT_AFFINE_CANDIDATES':raise RuntimeError('Exact affine candidates required')
labels=['Disk / spatial max','Disk / regional mean','Line / spatial max','Line / regional mean']
base=[v['baseline_loss_um'] for v in result['candidates']]
affine=[v['exact_affine_result']['loss_display'] for v in result['candidates']]
steps=np.array([v['exact_affine_result']['step_display'] for v in result['candidates']])
plt.style.use('dark_background');fig,(left,right)=plt.subplots(1,2,figsize=(16,6),layout='constrained')
x=np.arange(4);left.bar(x-.18,base,.36,color='#5cc8ff',label='directly verified common baseline')
left.bar(x+.18,affine,.36,color='#ffb86b',label='exact optimum of local AFFINE model')
for i,(a,b) in enumerate(zip(base,affine)):
    left.text(i-.18,a+.04,f'{a:.4f}',ha='center',fontsize=9)
    left.text(i+.18,b+.04,f'{b:.4f}',ha='center',fontsize=9)
left.set_xticks(x,labels,rotation=18,ha='right',fontsize=9)
left.set_ylim(0,7.4);left.set_ylabel('uniform endpoint discrepancy (µm)');left.legend(fontsize=8)
left.set_title('Same 17 base states and gain for all four conventions');left.grid(axis='y',alpha=.15)
scaled=steps/np.array([.1,.1,.02,.1,.1]);im=right.imshow(scaled,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto')
right.set_xticks(range(5),['Δlog Ka','Δlog Kb','ΔTa','Δlog Tb','Δlog gain'],rotation=20,ha='right')
right.set_yticks(range(4),labels,fontsize=9)
for i in range(4):
    for j in range(5):right.text(j,i,f'{steps[i,j]:.5f}',ha='center',va='center',fontsize=10,color='black' if abs(scaled[i,j])<.45 else 'white')
right.set_title('Proposed local steps often reach numerical trust bounds')
fig.colorbar(im,ax=right,label='step / declared positive bound (.1 for logs; .02 for Ta)')
fig.suptitle('All-record multivariate tangent calibration: proposals, not nonlinear improvements\n17 records × 4 endpoints; no new equilibria, physical error bars, confidence region or full identified set',fontsize=13)
figure=OUT/'material-cohort-affine-descent.png';fig.savefig(figure,dpi=170);plt.close(fig)
receipt={'scope':'Display of completed exact affine LPs and already verified common baseline; no new nonlinear states',
    'generator_sha256':sha(__file__),'exact_candidate_sha256':sha(path),'figure_sha256':sha(figure),
    'all4_use_same17_baseline':True,'no_measurement_or_model_errorbars':True}
save=OUT/'visual_verification.json';save.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
