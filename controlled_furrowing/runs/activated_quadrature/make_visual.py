"""Receipt-derived scientific plots; no model solves or target data."""
from pathlib import Path
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[2]
out=Path(__file__).resolve().parent
receipt=json.loads((root/'evidence/activated_quadrature_verification.json').read_text())
plt.style.use('dark_background')
fig,axes=plt.subplots(1,3,figsize=(14,4.3),layout='constrained')
colors={500.:'#5cc8ff',5000.:'#ffb86b'}
for lam,color in colors.items():
    rows=[r for r in receipt['records'] if r['lambda']==lam];f=[r['strength'] for r in rows]
    axes[0].semilogy(f,[r['old_higher_force'] for r in rows],'o--',color=color,label=f'old field, λ={lam:g}, Q8')
    axes[0].semilogy(f,[r['new_force12'] for r in rows],'s-',color=color,label=f'new field, λ={lam:g}, Q12')
    axes[1].plot(f,[r['shortening_change_um'] for r in rows],'o-',color=color,label=f'λ={lam:g}')
axes[0].axhline(1e-7,color='white',ls=':',lw=1,label='acceptance threshold')
axes[0].set(xlabel='activation f/µ',ylabel='scaled force residual',title='Quadrature stationarity corrected')
axes[0].legend(fontsize=7)
axes[1].axhline(0,color='white',ls=':',lw=1)
axes[1].set(xlabel='activation f/µ',ylabel='new − old shortening (µm)',title='Same-space observable correction')
axes[1].legend(fontsize=8)
x=np.arange(6);width=.36
for offset,name,color in [(-width/2,'apical','#5cc8ff'),(width/2,'basal','#bd93f9')]:
    axes[2].bar(x+offset,[r['exact_same_parameter_curve_difference_bounds'][name]['upper_bound_display_um'] for r in receipt['records']],width,label=name,color=color)
axes[2].set_xticks(x,[f"{r['lambda']:g}\nf={r['strength']:g}" for r in receipt['records']],fontsize=8)
axes[2].set(xlabel='penalty λ/µ and activation',ylabel='curve displacement bound (µm)',title='Exact discrete curve-difference bounds')
axes[2].legend(fontsize=8)
for a in axes:a.grid(alpha=.15)
fig.suptitle('Same 12×6 P2 space, same physical law: default rule → Duffy8/surface12\nHigher-rule cross-check passes; no continuum/PDE error certificate',fontsize=12)
fig.savefig(out/'quadrature-overview.png',dpi=170);plt.close(fig)

old=json.loads((root/'runs/quadratic_robust/lambda_5000.0_target_10.0/final.json').read_text())
new=json.loads((out/'lambda_5000.0_strength_10.0/reequilibrated.json').read_text())
reference=np.array(old['reference_nodes']);radii=np.unique(reference[reference[:,1]==old['parameters']['height'],0])[::2]
fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
for index,name in enumerate(('top','bottom')):
    key='boundary_'+name+'_quadratic_controls';before=np.array(old[key]);after=np.array(new[key]);radius=[];a=[];b=[]
    for i,(v,w) in enumerate(zip(before,after)):
        t=np.linspace(0,1,41)
        basis=np.stack(((1-t)*(1-2*t),4*t*(1-t),t*(2*t-1)),axis=1)
        radius.extend((radii[i]+t*(radii[i+1]-radii[i]))*33.5)
        a.extend((basis@v)[:,1]);b.extend((basis@w)[:,1])
    a=np.array(a);b=np.array(b);radius=np.array(radius);baseline=1. if name=='top' else 0.
    axes[index].plot(radius,(a-baseline)*33.5,'--',color='#ffb86b',label='default-rule equilibrium')
    axes[index].plot(radius,(b-baseline)*33.5,color='#5cc8ff',label='higher-rule equilibrium')
    axes[index].set(title='Apical' if name=='top' else 'Basal',xlabel='reference material radius (µm)',ylabel='vertical displacement (µm)')
    axes[2].plot(radius,(b-a)*33.5*1000,label='apical' if name=='top' else 'basal')
for a in axes:a.set_xlim(0,22);a.grid(alpha=.15);a.legend(fontsize=7)
axes[2].set(title='Small correction, resolved separately',xlabel='reference material radius (µm)',ylabel='new − old vertical displacement (nm)')
fig.suptitle('Actual small-region model, λ/µ=5000, f/µ=10\nExact quadratic curves sampled for display; comparison uses common material radius',fontsize=12)
fig.savefig(out/'quadrature-profiles.png',dpi=170);plt.close(fig)
record={'scope':'Visuals generated entirely from retained numerical receipts and fields','generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'input_receipt_sha256':hashlib.sha256((root/'evidence/activated_quadrature_verification.json').read_bytes()).hexdigest(),
        'figures':{name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ['quadrature-overview.png','quadrature-profiles.png']}}
(root/'evidence/activated_quadrature_visuals.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record))
