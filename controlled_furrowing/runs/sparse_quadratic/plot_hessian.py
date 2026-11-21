from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[2];out=root/'runs/sparse_quadratic'
z=np.load(out/'small_hessian_comparison.npz');dense=z['dense'];colors=z['colors'];V=np.eye(int(colors.max())+1)[colors];compressed=dense@V
plt.style.use('dark_background');fig,axs=plt.subplots(1,3,figsize=(14,4.8),layout='constrained')
limit=max(np.log10(np.max(abs(dense))),np.log10(np.max(abs(compressed))))
for ax,matrix,title in [(axs[0],dense,'Dense AD Hessian: 161 x 161'),(axs[1],compressed,'Compressed AD products: 161 x 42')]:
 mag=np.ma.masked_where(matrix==0,np.log10(np.maximum(abs(matrix),1e-30)));im=ax.imshow(mag,cmap='magma',vmin=-4,vmax=limit,aspect='auto',interpolation='nearest');ax.set(title=title,xlabel='Coordinate / color column',ylabel='Free-coordinate row')
fig.colorbar(im,ax=axs[:2],label='log10 absolute derivative magnitude',shrink=.8)
bench=[json.loads((root/f'evidence/sparse_quadratic_benchmark_{k}.json').read_text()) for k in ['dense','sparse']]
axs[2].bar(['Dense AD','Colored sparse AD'],[q['sampled_incremental_RSS_bytes']/2**20 for q in bench],color=['#d9a765','#7bc8b6']);axs[2].set(ylabel='Sampled incremental RSS (MiB)',title='Separate-process small-mesh measurement')
for i,q in enumerate(bench):axs[2].text(i,q['sampled_incremental_RSS_bytes']/2**20+2,f"{q['seconds']:.3f} s",ha='center')
fig.suptitle('Same energy and gradient; full recovered-H comparison relative error 1.06e-16; no physical-accuracy claim',fontsize=12)
p=out/'sparse_hessian_verification.png';fig.savefig(p,dpi=170);plt.close(fig)
inputs=[out/'small_hessian_comparison.npz',root/'evidence/sparse_quadratic_verification.json']+[root/f'evidence/sparse_quadratic_benchmark_{k}.json' for k in ['dense','sparse']]
receipt={'image':str(p.relative_to(root)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'inputs':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},'scope':'Small6x3verification/benchmarkonly; no24meshperformanceextrapolationorbiologyinference.'}
(root/'evidence/sparse_quadratic_visual.json').write_text(json.dumps(receipt,indent=2));print(str(p))
