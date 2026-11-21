from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('controlled_furrowing');out=root/'runs/quadratic_robust';sources=[];records=[]
for target in [5.,10.]:
 p=out/f'lambda_5000.0_target_{target}'/'final.json';r=json.loads(p.read_text());sources.append(p);records+=r['robust_continuation']['attempts']
plt.style.use('dark_background');fig,axs=plt.subplots(2,1,figsize=(12,8),sharex=True,layout='constrained')
last=1.;accepted=[]
for i,a in enumerate(records):
 color='#78c6b1' if a['accepted'] else '#e58b70';axs[0].plot([i,i],[a['from_strength'],a['to_strength']],color=color,lw=2,alpha=.65);axs[0].scatter(i,a['to_strength'],color=color,marker='o' if a['accepted'] else 'x',s=65)
 if a['accepted']:last=a['to_strength']
 accepted.append(last);axs[1].scatter(i,a['final_force_linf'],color=color,marker='o' if a['accepted'] else 'x',s=65)
axs[0].step(np.arange(len(records)),accepted,where='post',color='#c8c8c8',ls='--',label='Last verified gain');axs[0].set(ylabel='Activation gain / mu',title='Every attempt starts from the last verified state; rejected vectors are discarded');axs[0].legend(fontsize=9)
axs[1].axhline(1e-7,color='#dfbd68',ls='--',label='Unchanged force threshold 1e-7');axs[1].set(yscale='log',xlabel='Numerical attempt (targets 5 then 10; lambda / mu = 5000)',ylabel='Normalized free-force infinity norm');axs[1].legend(fontsize=9)
fig.suptitle('Robust continuation of the frozen full-P2 energy: actual accepted and rejected computations',fontsize=12)
path=out/'solver_policy.png';fig.savefig(path,dpi=170);plt.close(fig)
(root/'evidence/quadratic_robust_visual.json').write_text(json.dumps({'file':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},'scope':'Actualtargetfreeattempts, notfailure-injectiontests; no biologicalmeasurementortrajectoryinterpretation'},indent=2));print(str(path))
