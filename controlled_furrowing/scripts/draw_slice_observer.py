"""Original schematic of the missing-slice observation ambiguity."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('dark_background')
r=np.linspace(0,2,501)
top=np.interp(r,[0,1,2],[1,.8,1])
fig,axes=plt.subplots(1,2,figsize=(11.8,4.3),constrained_layout=True)
for ax,slice_z,label in zip(axes,[.9,.7],['A slice intersects a plug and an empty annulus','A lower slice contains continuous tissue']):
    ax.fill_between(r,0,top,color='#467c75',alpha=.7,label='Tissue body')
    ax.plot(r,top,color='#9cdec1',lw=2,label='Apical surface')
    ax.axhline(slice_z,color='#f6bd73',ls='--',lw=2,label='Image slice')
    absent=top<slice_z
    if absent.any():
        ax.fill_between(r,top,slice_z,where=absent,color='#d87988',alpha=.5,label='Empty annulus in the slice')
    ax.set(xlabel='Radius (arbitrary units)',ylabel='Height',ylim=(-.03,1.1),title=label)
    ax.legend(fontsize=8,loc='lower left')
fig.suptitle('Identical tissue surfaces can receive different image labels when slice height is unknown',fontsize=12)
root=Path(__file__).resolve().parents[1]
fig.savefig(root/'docs/assets/slice_observer.png',dpi=170)
plt.close(fig)
