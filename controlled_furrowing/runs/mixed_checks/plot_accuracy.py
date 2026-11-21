from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('controlled_furrowing');plt.style.use('dark_background')
fixed=json.loads((root/'evidence/mixed_fixed10.json').read_text());kernel=json.loads((root/'evidence/mixed_kernel_audit.json').read_text());manu=json.loads((root/'evidence/quadratic_manufactured.json').read_text())
fig,axs=plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
rows=sorted([r for r in fixed['rows'] if r['lambda']==500],key=lambda r:r['nr']);n=[r['nr'] for r in rows]
axs[0].plot(n,[r['sampled_minimum_J'] for r in rows],'o-',label='Sampled minimum J');axs[0].plot(n,[r['sampled_maximum_J'] for r in rows],'s-',label='Sampled maximum J');axs[0].axhline(1,color='grey',lw=.7);axs[0].set(xlabel='Radial cells (depth = radial / 2)',ylabel='Actual quadrature J',title='Projected constraint: local volume error persists');axs[0].legend(fontsize=8)
axs[1].plot(n,[abs(r['relative_volume_change']) for r in rows],'o-');axs[1].set(yscale='log',xlabel='Radial cells',ylabel='Absolute total volume change / V0',title='Small global volume error is insufficient')
k=sorted([r for r in kernel['rows'] if r['nr']==8],key=lambda r:r['lambda']);l=[r['lambda'] for r in k]
axs[2].plot(l,[r['projected_logJ_RMS'] for r in k],'o-',label='Projected logJ RMS');axs[2].plot(l,[r['unresolved_logJ_RMS'] for r in k],'s-',label='Unresolved logJ RMS');axs[2].set(xscale='log',yscale='log',xlabel='Volume coefficient lambda / mu',ylabel='RMS logJ',title='8 x 4: increasing penalty misses its kernel');axs[2].legend(fontsize=8)
fig.suptitle('Negative numerical result: P2 / P1 projected log-volume control does not enforce local incompressibility',fontsize=12)
path=root/'runs/mixed_checks/mixed_summary.png';fig.savefig(path,dpi=170);plt.close(fig)
(root/'evidence/mixed_visual.json').write_text(json.dumps({'file':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'inputs':{str(root/'evidence'/p):hashlib.sha256((root/'evidence'/p).read_bytes()).hexdigest() for p in ['mixed_fixed10.json','mixed_kernel_audit.json']},'scope':'Synthetic fixedgain10; plottedJvaluesareactualquadrature, notBernsteinbounds; overallcontinuationfailuresretained'},indent=2))
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained');colors=['#dba85e','#79c5b5','#b09bdb']
for lam,color in zip([5.,50.,500.],colors):
 r=sorted([q for q in manu['rows'] if q['lambda']==lam],key=lambda q:q['nr']);labels=[f"{q['nr']} x {q['nz']}" for q in r]
 for ax,key,title in [(axs[0,0],'displacement_L2_RMS','Displacement error against exact continuum map'),(axs[0,1],'gradient_H1_RMS','Full 3D gradient error, including hoop'),(axs[1,0],'J_max_sampled_error','Sampled local volume error')]:
  ax.plot(labels,[q['solved_errors'][key] for q in r],'o-',label=f'lambda / mu = {lam}',color=color);ax.set(yscale='log',xlabel='Meridional cell grid',ylabel=key,title=title)
axs[0,0].legend(fontsize=8)
lo=sorted([r for r in manu['rows'] if r['lambda']==5],key=lambda r:r['nr']);hi=sorted([r for r in manu['rows'] if r['lambda']==500],key=lambda r:r['nr'])
axs[1,1].plot(labels,[a['solved_errors']['gradient_H1_RMS']/b['solved_errors']['gradient_H1_RMS'] for a,b in zip(hi,lo)],'o-',color='#da9d6d');axs[1,1].set(xlabel='Meridional cell grid',ylabel='H1 error at lambda500 / H1 error at lambda5',title='Penalty sensitivity remains after refinement');axs[1,1].axhline(1,color='grey',lw=.7)
fig.suptitle('Forced manufactured equilibrium: same nonlinear exact solution for every volume coefficient',fontsize=12)
path=root/'runs/quadratic_checks/quadratic_accuracy.png';fig.savefig(path,dpi=170);plt.close(fig)
(root/'evidence/quadratic_visual.json').write_text(json.dumps({'file':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'input':str(root/'evidence/quadratic_manufactured.json'),'input_sha256':hashlib.sha256((root/'evidence/quadratic_manufactured.json').read_bytes()).hexdigest(),'scope':'Solvedtargetindependentforcedbenchmark; empiricalquadratureerrorsnotrigorousPDEerrorbounds'},indent=2))
print('Both figures generated and hashed')
