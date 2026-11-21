"""Actual exact whole-root geometry bounds; no measured data or ideal mesh."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

ROOT=Path(__file__).resolve().parents[1]
receipt_path=ROOT/'evidence/root_body_geometry_verification.json'
d=json.loads(receipt_path.read_text(encoding='utf-8'))
if d['status']!='PASS_ACTUAL_ROOT_ENTIRE_BODY_GEOMETRY':raise RuntimeError('actual root-box geometry not certified')
c=json.loads((ROOT/'evidence/algebraic_root_certification_feasibility.json').read_text(encoding='utf-8'))
state=json.loads((ROOT/c['input_state']).read_text(encoding='utf-8'))
records=d['geometry']['triangles'];low=np.array([float(Q(t['J_bound']['lower_exact'])) for t in records]);high=np.array([float(Q(t['J_bound']['upper_exact'])) for t in records])
plt.style.use('dark_background');fig,axes=plt.subplots(1,2,figsize=(13,4.7))
polygons=[np.array([state['reference_nodes'][i] for i in face[:3]])*33.5 for face in state['triangles']]
collection=PolyCollection(polygons,array=low,cmap='viridis',edgecolors='#536274',linewidths=.35)
axes[0].add_collection(collection);axes[0].autoscale();axes[0].set_aspect('equal');axes[0].set(xlabel='Reference radius (µm)',ylabel='Reference height (µm)',title='Certified lower J on every triangle')
fig.colorbar(collection,ax=axes[0],label='Whole-triangle lower bound')
indices=np.arange(len(records));axes[1].fill_between(indices,low,high,color='#68c4ec',alpha=.35,label='Exact interval enclosure')
axes[1].plot(indices,low,color='#81d8b5',lw=1);axes[1].plot(indices,high,color='#efb96c',lw=1);axes[1].axhline(1,color='white',ls=':',lw=1)
axes[1].set(xlabel='Triangle index (not spatial distance)',ylabel='Axisymmetric Jacobian J',title='Every field in the actual root box');axes[1].legend();axes[1].grid(alpha=.2)
lo_q=Q(d['geometry']['J_bound']['lower_exact']);hi_q=Q(d['geometry']['J_bound']['upper_exact']);m_q=Q(d['geometry']['strong_monotonicity_lower_exact'])
lo=float(lo_q);hi=float(hi_q);m=float(m_q)
def safe_decimal(value,digits,upper=False):
    power=10**digits
    integer=-(-value.numerator*power//value.denominator) if upper else value.numerator*power//value.denominator
    return f'{integer//power}.{integer%power:0{digits}d}'
display=[safe_decimal(lo_q,9),safe_decimal(hi_q,9,True)];m_display=safe_decimal(m_q,6)
fig.suptitle(f'One archived weak-load root: 144 whole-element certificates, J in [{display[0]}, {display[1]}]\nUniform meridional monotonicity lower bound {m_display}; finite-model geometry, not continuum accuracy',fontsize=12)
fig.tight_layout();path=ROOT/'docs/assets/root_body_geometry.png';fig.savefig(path,dpi=170,bbox_inches='tight');plt.close(fig)
out={'status':'GENERATED_FROM_EXACT_ROOT_GEOMETRY','safe_J_display':display,'safe_monotonicity_lower_display':m_display,'image_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
     'source_receipt_sha256':hashlib.sha256(receipt_path.read_bytes()).hexdigest(),'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'scope':'Whole-box mathematical geometry enclosures of one root; not measured tissue, continuum accuracy, energetic stability or parameter-set inference.'}
(ROOT/'evidence/root_body_geometry_visual.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8');print(json.dumps({'status':out['status'],'image':str(path),'J':[lo,hi],'m':m}))
