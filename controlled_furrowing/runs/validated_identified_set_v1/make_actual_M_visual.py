"""Receipt-derived cover/constraint scope picture; not a parameter-volume plot."""
from pathlib import Path
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent/'actual_M_cover_v2'
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
path=ROOT/'evidence/validated_identified_set_actual_M_cover.json';r=read(path)
if r['status']!='PASS'or not all(r['checks'].values()):raise RuntimeError('Complete verified actual cover required')
target=r['target_leaf'];nodes=r['cover_ledger']['nodes'];decision=nodes[target]['classification'];splits=r['exact_splits']
plt.style.use('dark_background');fig,(left,right)=plt.subplots(1,2,figsize=(15,8),layout='constrained',gridspec_kw={'width_ratios':[1.2,1]})
short={'apical_area_modulus':'Ka','basal_area_modulus':'Kb','apical_surface_tension':'Ta','basal_surface_tension':'Tb','common_gain':'gain'}
for i,s in enumerate(splits):
    node=s['parent'];children=nodes[node]['children'];other=next(n for n in children if n!=s['target_child'])
    left.plot([0,0],[i,i+1],color='#5cc8ff',lw=2)
    left.plot([0,1],[i,i+.2],color='#a0a0a0',lw=1.5);left.scatter([1],[i+.2],color='#a0a0a0',s=40)
    side='lower'if s['target_child'].endswith('1')else'upper'
    left.text(-.12,i+.35,f'{short[s["coordinate"]]} {side} cut',ha='right',fontsize=10)
    left.text(1.08,i+.2,f'{other}: unresolved',va='center',fontsize=8)
left.scatter([0],[len(splits)],color='#ff6b6b',s=110,zorder=4)
left.text(.15,len(splits),f'{target}: selected local family OUT',va='center',fontsize=10,color='#ff8c8c')
left.text(-.12,-.3,'Broader finite seed',ha='right',fontsize=10)
left.set_xlim(-1.9,3.0);left.set_ylim(len(splits)+.8,-.8);left.axis('off')
left.set_title('9 exact cuts retain 10 finite cells\nDiagram shows proof bookkeeping, not 5D cell volumes',fontsize=12)
ids=list(dict.fromkeys(row['record_id']for row in decision['constraints']))
keys=['maximum_shortening_um','region_averaged_shortening_um','apical_max_shortening_um','basal_max_shortening_um']
matrix=np.zeros((17,4))
for row in decision['constraints']:
    if row['status']=='OUT':matrix[ids.index(row['record_id']),keys.index(row['observable'])]=1
right.imshow(matrix,cmap=ListedColormap(['#4a5260','#e66b6b']),vmin=0,vmax=1,aspect='auto')
right.set_xticks(range(4),['Maximum','Mean','Apical','Basal']);right.set_yticks(range(17),[rid.rsplit('_',1)[-1]for rid in ids])
right.set_ylabel('Published source columns B–R, not verified embryo IDs')
right.set_title('All 68 source constraints remain visible\nOnly M/apical has a uniform forward proof',fontsize=12)
right.text(2,ids.index('optoGEF2_small_col_M'),'OUT',ha='center',va='center',color='white',fontsize=10)
right.set_xticks(np.arange(-.5,4,1),minor=True);right.set_yticks(np.arange(-.5,17,1),minor=True);right.grid(which='minor',color='black',lw=.6);right.tick_params(which='minor',bottom=False,left=False)
fig.suptitle('Actual M: one selected stationary-family cell excluded at DECLARED δ = 5.1 µm\nOther 67 constraints, other 9 finite cells, other branches and unbounded/zero-face remainder stay unresolved',fontsize=13)
figure=OUT/'actual-M-cover-and-constraints.png';fig.savefig(figure,dpi=170);plt.close(fig)
receipt={'scope':'Exact cover topology and coverage status from retained receipt; not a parameter-volume, confidence or biological-rejection plot',
    'source_receipt_sha256':sha(path),'generator_sha256':sha(__file__),'figure_sha256':sha(figure)}
(OUT/'visual_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
