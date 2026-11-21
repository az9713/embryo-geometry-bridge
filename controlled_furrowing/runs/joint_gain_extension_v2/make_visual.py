"""Final receipt-derived calibration diagnostics; never plot a partial replay."""
from pathlib import Path
from fractions import Fraction as Q
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[2];out=Path(__file__).resolve().parent
job=json.loads((out/'direct_job.json').read_text())
if job.get('status')!='FINISHED':raise RuntimeError('Full direct replay must finish before final figures')
result_path=out/'direct_results.json';result=json.loads(result_path.read_text())
if result['status']!='COMPLETE_ALL_FOUR_CONDITIONAL_CANDIDATE_REPLAYS':raise RuntimeError('No complete direct result')
with (root/'data/derived/calibration_small_18um.csv').open(newline='',encoding='utf-8-sig') as stream:rows=list(csv.DictReader(stream))
ids=[r['record_id'] for r in rows];dose={r['record_id']:float(r['activation_magnitude']) for r in rows}
keys=['maximum_shortening_um','region_averaged_shortening_um','apical_max_shortening_um','basal_max_shortening_um']
titles=['Maximum current-radius shortening','Regional mean shortening','Apical downward endpoint','Basal upward endpoint']
labels={'current_disk__spatial_max':'Disk / spatial-max components','current_disk__regional_mean':'Disk / regional-mean components',
        'current_radial_line__spatial_max':'Line / spatial-max components','current_radial_line__regional_mean':'Line / regional-mean components'}
colors=['#5cc8ff','#ffb86b','#bd93f9','#50fa7b'];plt.style.use('dark_background')
fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
for k,(key,title,ax) in enumerate(zip(keys,titles,axes.ravel())):
    ax.scatter([dose[r] for r in ids],[float(rows[i][key]) for i,r in enumerate(ids)],s=35,color='white',label='source-column observations',zorder=5)
    for variant,color in zip(result['variants'],colors):
        channel={r['record_id']:r for r in variant['all68_residuals']['constraint_rows'] if r['observable']==key}
        ordered=sorted(ids,key=lambda r:dose[r]);mid=[];err=[]
        for rid in ordered:
            lo,hi=map(Q,channel[rid]['prediction_interval_exact']);mid.append(float((lo+hi)/2));err.append(float((hi-lo)/2))
        gain=float(Q(variant['selected_exact_SURROGATE_gain']))
        ax.errorbar([dose[r] for r in ordered],mid,yerr=err,fmt='.-',lw=1,ms=4,color=color,label=f"{labels[variant['variant']]}, g={gain:.4f}")
    ax.set(title=title,xlabel='normalized activation input',ylabel='µm');ax.grid(alpha=.15)
axes[0,0].legend(fontsize=7)
fig.suptitle('Direct nonlinear replay at four surrogate-selected gains\nOne fixed 12×6 high-Q scenario; static peak proxy; line joins are display guides, not nonlinear-optimality certificates',fontsize=12)
fig.savefig(out/'joint-direct-endpoints.png',dpi=170);plt.close(fig)

matrix=[]
for variant in result['variants']:
    lookup={(r['record_id'],r['observable']):r for r in variant['all68_residuals']['constraint_rows']}
    values=[]
    for key in keys:
        for rid in ids:
            row=lookup[rid,key];lo,hi=map(Q,row['prediction_interval_exact']);values.append(float((lo+hi)/2-Q(row['observed_exact'])))
    matrix.append(values)
matrix=np.array(matrix);limit=max(abs(matrix.min()),abs(matrix.max()))
fig,ax=plt.subplots(figsize=(16,4.8),layout='constrained')
im=ax.imshow(matrix,aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit)
ax.set_yticks(range(4),[labels[v['variant']] for v in result['variants']],fontsize=9)
ax.set_xticks(range(68),[r.rsplit('_',1)[-1] for key in keys for r in ids],fontsize=6)
ax.set_xlabel('Source-column records B–R repeated by endpoint; canonical pairing assumption, not verified embryo IDs')
top=ax.secondary_xaxis('top');top.set_xticks([8,25,42,59],['Maximum','Regional mean','Apical','Basal'])
for x in (16.5,33.5,50.5):ax.axvline(x,color='white',lw=1)
fig.colorbar(im,ax=ax,label='direct midpoint prediction − recorded endpoint (µm)')
ax.set_title('All 68 residuals under each retained observation convention\nSame 68 observations viewed 4 ways; no discarded high-residual records or inferred noise law',fontsize=12)
fig.savefig(out/'joint-direct-residuals.png',dpi=170);plt.close(fig)

fig,ax=plt.subplots(figsize=(10,4.6),layout='constrained');x=np.arange(4)
sur=[float(Q(v['surrogate_loss_exact'])) for v in result['variants']]
lo=[v['direct_candidate_discrepancy_requirement_um_display'][0] for v in result['variants']]
hi=[v['direct_candidate_discrepancy_requirement_um_display'][1] for v in result['variants']]
ax.bar(x-.18,sur,.36,label='defined exact midpoint-table surrogate')
ax.bar(x+.18,lo,.36,yerr=[np.zeros(4),np.array(hi)-lo],label='direct candidate requirement interval')
ax.set_xticks(x,[labels[v['variant']] for v in result['variants']],fontsize=8);ax.set_ylabel('uniform discrepancy requirement (µm)');ax.grid(axis='y',alpha=.15);ax.legend(fontsize=8)
ax.set_title('Surrogate optimum versus actual candidate replay\nNeither bar estimates measurement noise or proves a nonlinear/physical global optimum',fontsize=12)
fig.savefig(out/'joint-direct-surrogate-comparison.png',dpi=170);plt.close(fig)
files=['joint-direct-endpoints.png','joint-direct-residuals.png','joint-direct-surrogate-comparison.png']
receipt={'scope':'Final complete direct-candidate figures, no new solves or heldout responses',
    'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'direct_results_sha256':hashlib.sha256(result_path.read_bytes()).hexdigest(),
    'calibration_csv_sha256':hashlib.sha256((root/'data/derived/calibration_small_18um.csv').read_bytes()).hexdigest(),
    'figures':{f:hashlib.sha256((out/f).read_bytes()).hexdigest() for f in files}}
(root/'evidence/joint_gain_extension_visuals.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))

