"""Visualize the retained exact discrete-curve comparison; no solves/data."""
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[2];out=Path(__file__).resolve().parent
path=root/'evidence/joint_refinement_comparison.json';d=json.loads(path.read_text())['comparison']
plt.style.use('dark_background');fig,axes=plt.subplots(1,3,figsize=(14,4.2),layout='constrained')
for i,side in enumerate(('top','bottom')):
    curve=d[side+'_Eulerian_vertical'];rows=curve['sampled_eulerian_curves'];x=np.array([p['radius_um'] for p in rows]);f=np.array([p['fine_z_um'] for p in rows]);c=np.array([p['coarse_z_um'] for p in rows]);baseline=33.5 if side=='top' else 0.
    axes[i].plot(x,c-baseline,color='#ffb86b',label='12×6, higher rule')
    axes[i].plot(x,f-baseline,color='#5cc8ff',label='16×8, higher rule')
    axes[i].set(title='Apical curve' if side=='top' else 'Basal curve',xlabel='deformed physical radius (µm)',ylabel='vertical displacement (µm)')
    axes[2].plot(x,f-c,label='apical' if side=='top' else 'basal')
for a in axes:a.set_xlim(0,25);a.grid(alpha=.15);a.legend(fontsize=8)
axes[2].axhline(0,color='white',ls=':',lw=1)
axes[2].set(title='Fine − coarse at the same physical radius',xlabel='deformed physical radius (µm)',ylabel='height difference (µm)')
sup=d['top_Eulerian_vertical']['supremum_absolute_discrete_height_difference_um']
fig.suptitle(f"After integration correction: λ/µ=5000, f/µ=10\nCertified discrete apical vertical sup ∈ [{sup['lower']['display']:.6f}, {sup['upper']['display']:.6f}] µm; not a PDE error bound",fontsize=12)
fig.savefig(out/'joint-refinement.png',dpi=170);plt.close(fig)
record={'scope':'Plots of retained common-Eulerian discrete curves, sampled for display; exact interval remains in comparison receipt',
        'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_comparison_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'figure_sha256':hashlib.sha256((out/'joint-refinement.png').read_bytes()).hexdigest()}
(root/'evidence/joint_refinement_visuals.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))
