from pathlib import Path
from fractions import Fraction as Q
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;d=json.loads((P/'results.json').read_text());rows=d['rows'];x=np.arange(len(rows))
plt.style.use('dark_background');fig,axes=plt.subplots(1,2,figsize=(13,5.4),constrained_layout=True)
for ax,un,pt,tail,title in [(axes[0],'force','point_force','force_tail_upper_exact','Local force discrepancy'),(axes[1],'H_entries','point_H_entries','H_tail_upper_exact','Local Hessian-entry discrepancy')]:
 for offset,key,label,color in [(-.18,pt,'Corrected point norm upper','#61dafb'),(0,un,'Uniform norm upper','#ffca66')]:
  ax.scatter(x+offset,[float(Q(r[key]['norm_upper_exact']))for r in rows],label=label,color=color,s=48,zorder=4)
 ax.scatter(x+.18,[float(Q(r[un]['max_halfwidth_exact']))for r in rows],label='Uniform max interval radius',marker='x',color='#ff7683',s=48)
 ax.scatter(x+.3,[float(Q(r[tail]))for r in rows],label='Two-mass tail upper',marker='_',color='#a4ef95',s=65)
 ax.set_yscale('log');ax.set_xticks(x,[str(r['element'])for r in rows]);ax.set_xlabel('Predeclared bulk element');ax.set_ylabel('Dimensionless normalized local contribution');ax.set_title(title);ax.grid(alpha=.18,which='both')
 ax.axvspan(3.55,5.45,color='#ff7683',alpha=.08)
axes[0].legend(fontsize=8,loc='lower left');fig.suptitle('Seven-element uniform pilot: two axis-edge enclosures are broad\nInterval radii are proof uncertainty, not experimental error; no assembled or root bound',fontsize=12)
fig.savefig(P/'local_bounds.png',dpi=180);fig.savefig(P/'local_bounds.svg');plt.close(fig)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
(P/'visual_generation.json').write_text(json.dumps({'status':'GENERATED_PENDING_NATIVE_INSPECTION','source_sha256':sha(__file__),'input_sha256':sha(P/'results.json'),'outputs':{s:sha(P/s)for s in ['local_bounds.png','local_bounds.svg']},'scope':'Local seven-element norms and exact interval halfwidths only; approximate plot values from exact rational packet endpoints.'},indent=2)+'\n')
