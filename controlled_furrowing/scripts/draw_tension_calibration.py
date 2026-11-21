"""Plot independently checked small-region adequacy, not held-out performance."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
receipt=json.loads((ROOT/'evidence/tension_calibration_adequacy.json').read_text(encoding='utf-8'))
branches=receipt['verification_results']['branches']
plt.style.use('dark_background')
x=np.arange(len(branches));labels=[f'{b["basal_tension"]:g}' for b in branches]
fig,axes=plt.subplots(1,3,figsize=(14.2,4.8),constrained_layout=True)
for field,title,color in [('height_proxy_RMSE_um','Height proxy','#8ad4ef'),('apical_component_RMSE_um','Apical component','#f0be7c'),('basal_component_RMSE_um','Basal component','#c8aae8')]:
    axes[0].plot(x,[b[field] for b in branches],'o-',label=title,color=color)
axes[0].set(ylabel='Small-region RMSE (µm)',title='Height fit alone conceals motion partition')
axes[0].legend(fontsize=8)
axes[1].bar(x,[b['apical_dominant_predicted'] for b in branches],color='#8ad4ef')
axes[1].axhline(16,color='#f0be7c',ls='--',label='Source proxy: 16 of 17')
axes[1].set(ylabel='Predicted apical-dominant cases / 17',ylim=(0,18),title='Baseline tension changes the motion direction')
axes[1].legend(fontsize=8)
axes[2].bar(x-.17,[100*(1-b['minimum_J']) for b in branches],width=.34,color='#e5a5a5',label='Worst local compression')
axes[2].bar(x+.17,[100*b['maximum_absolute_global_volume_change'] for b in branches],width=.34,color='#a4d7b8',label='Global volume change')
axes[2].set(ylabel='Volume change (%)',yscale='log',title='Improved direction does not fix volume loss')
axes[2].legend(fontsize=8)
for ax in axes:
    ax.set_xticks(x,labels);ax.set_xlabel('Declared basal tension / µh')
fig.suptitle('Calibration-only prestress revision: every scenario retained; no continuum or biological validation',fontsize=12)
fig.savefig(ROOT/'docs/assets/tension_calibration.png',dpi=170)
plt.close(fig)
