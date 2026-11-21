"""Original plots of declared cell-test statuses and the conditional tolerance cut."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

OUT=Path(__file__).resolve().parent
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
cases=[('Pilot / constant',OUT,'constant'),('Pilot / affine',OUT,'affine'),
       ('Actual M / constant',OUT/'actual_M','constant'),('Actual M / affine',OUT/'actual_M','affine')]
widths=[F(1,10**k)for k in (12,10,8,6,4,2)];data=[]
for label,folder,method in cases:
    p=read(folder/'protocol.json');r=read(folder/'results.json');scales={F(c['relative_scale_exact']):c['cell_id']for c in p['initial_cells']}
    data.append([int(next(q for q in r['cell_attempts']if q['cell_id']==scales[w]and q['method']==method)['status']=='CERTIFIED_UNIFORM_LOCAL_ROOT')for w in widths])
plt.rcParams.update({'font.size':11,'text.color':'#e4edf7','axes.labelcolor':'#dce6f2','xtick.color':'#bac9db','ytick.color':'#bac9db','axes.edgecolor':'#687d94','svg.hashsalt':'uniform-parameter-cell-evidence'})
fig,ax=plt.subplots(figsize=(9,4.6),facecolor='#0b1421');ax.set_facecolor('#111f2f');ax.imshow(data,cmap=ListedColormap(['#786146','#297d67']),vmin=0,vmax=1,aspect='auto')
ax.set_xticks(range(6),[r'$10^{-'+str(k)+r'}$'for k in(12,10,8,6,4,2)]);ax.set_yticks(range(4),[r[0]for r in cases]);ax.set_xlabel('Declared relative half-width in Ka, Kb, Tb and gain; Ta width = 0.02 times this scale')
for i,row in enumerate(data):
    for j,value in enumerate(row):ax.text(j,i,'CERTIFIED'if value else'UNRESOLVED',ha='center',va='center',color='#f4f7fb',fontsize=9)
ax.set_title('Whole five-parameter cells: tested sufficient conditions',loc='left',color='#f1f5fb',pad=15);fig.text(.03,.03,'Nested cells and methods overlap. Counts are not coverage probabilities. Two additional unresolved split leaves per target are retained.',fontsize=9,color='#aebed2');fig.subplots_adjust(left=.23,right=.98,bottom=.23,top=.86)
for suffix in('svg','png'):fig.savefig(OUT/('cell_status.'+suffix),facecolor=fig.get_facecolor(),dpi=160,metadata={'Date':None}if suffix=='svg'else None)
plt.close(fig)
d=read(OUT/'actual_M/contradiction_diagnostic.json');lo,hi=map(F,d['residual_interval_um_exact']);delta=F(d['declared_delta_um_exact'])
fig,ax=plt.subplots(figsize=(8.8,3.7),facecolor='#0b1421');ax.set_facecolor('#111f2f');ax.axvspan(5.05,float(lo),color='#765743',alpha=.45);ax.axvspan(float(hi),5.25,color='#256853',alpha=.35);ax.axvline(float(lo),color='#81d5ba',lw=2);ax.axvline(float(delta),color='#efc06d',ls='--',lw=2)
ax.annotate('Declared example: delta = 5.1',xy=(5.1,.55),xytext=(5.065,.84),color='#efc06d',arrowprops={'arrowstyle':'->','color':'#efc06d'});ax.text(5.065,.29,'This local branch cell is excluded',color='#f1c0a0');ax.text(5.196,.65,'This one row holds\nthroughout the cell',color='#98dfbf',fontsize=10);ax.annotate('Residual bracket near 5.186387\n(width about 2.18e-7 micrometres)',xy=(float(lo),.48),xytext=(5.12,.08),color='#d7e3f2',fontsize=10,arrowprops={'arrowstyle':'->','color':'#d7e3f2'})
ax.set_xlim(5.05,5.25);ax.set_ylim(0,1);ax.set_yticks([]);ax.set_xlabel('Declared apical discrepancy allowance delta (micrometres)');ax.set_title('A tolerance family from one actual-M uniform prediction',loc='left',color='#f1f5fb',pad=15);fig.text(.04,.03,'Other 67 constraints and other branches remain unresolved. Delta is not an estimated measurement-noise scale.',fontsize=9,color='#aebed2');fig.subplots_adjust(left=.09,right=.98,bottom=.24,top=.83)
for suffix in('svg','png'):fig.savefig(OUT/('delta_profile.'+suffix),facecolor=fig.get_facecolor(),dpi=160,metadata={'Date':None}if suffix=='svg'else None)
plt.close(fig)
receipt={'status':'PLOTS_GENERATED_FROM_FROZEN_RESULTS','generator_sha256':sha(__file__),'inputs':{str(f.relative_to(OUT)):sha(f)for f in [OUT/'results.json',OUT/'protocol.json',OUT/'actual_M/results.json',OUT/'actual_M/protocol.json',OUT/'actual_M/contradiction_diagnostic.json']},'assets':{name:sha(OUT/name)for name in ['cell_status.svg','cell_status.png','delta_profile.svg','delta_profile.png']},'cell_status_matrix':data,'visualization_only_not_new_certificate':True}
(OUT/'visual_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'assets':len(receipt['assets'])}))
