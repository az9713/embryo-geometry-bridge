from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[2];comparison=root/'evidence/activated_mesh_comparison.json';results=root/'evidence/activated_mesh_results.json'
r=json.loads(comparison.read_text());rows=r['comparisons'];fine=json.loads(results.read_text())['rows']
if len(rows)!=6:raise RuntimeError('Plot requires all six matched comparisons; no selection')
plt.style.use('dark_background');fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
colors=['#ddb16b','#79c5b3','#b4a0df']
for row in rows:
 color=colors[[1.,5.,10.].index(row['strength'])];style='-' if row['lambda']==500 else '--';label=f"λ = {row['lambda']:g}, gain = {row['strength']:g}"
 for ax,surface,title in [(axes[0,0],'top','Apical Eulerian height change'),(axes[0,1],'bottom','Basal Eulerian height change')]:
  points=row[surface+'_Eulerian_comparison']['sampled_eulerian_curves'];x=[p['radius_um'] for p in points];y=[(p['difference_lower_um']+p['difference_upper_um'])/2 for p in points]
  ax.plot(x,y,style,color=color,label=label);ax.axhline(0,color='grey',lw=.4);ax.set(xlabel='Current physical radius (micrometres); full 0-100.5 domain',ylabel='16×8 height minus 12×6 height (µm)',title=title)
axes[0,0].legend(fontsize=7,ncol=2)
for lam,marker in [(500.,'o'),(5000.,'s')]:
 rr=sorted([q for q in rows if q['lambda']==lam],key=lambda q:q['strength']);g=[q['strength'] for q in rr]
 for key,style in [('material_column_max_shortening_coarse_um','--'),('material_column_max_shortening_fine_um','-')]:
  axes[1,0].plot(g,[q[key] for q in rr],style,marker=marker,label=f"λ = {lam:g}, {'12×6' if 'coarse' in key else '16×8'}")
axes[1,0].set(xlabel='Unchanged activation gain',ylabel='Maximum material-column shortening (µm)',title='No refit to conceal mesh changes');axes[1,0].legend(fontsize=8)
labels=[f"{q['lambda']:g}\n{q['strength']:g}" for q in rows];positions=np.arange(6)
for offset,surface,color in [(-.13,'top','#79c5b3'),(.13,'bottom','#b4a0df')]:
 lower=[q[surface+'_Eulerian_comparison']['supremum_absolute_discrete_height_difference_um']['lower']['display'] for q in rows];upper=[q[surface+'_Eulerian_comparison']['supremum_absolute_discrete_height_difference_um']['upper']['display'] for q in rows];mid=(np.array(lower)+upper)/2
 axes[1,1].errorbar(positions+offset,mid,yerr=np.array([mid-lower,np.array(upper)-mid]),fmt='o',color=color,capsize=4,label='Apical' if surface=='top' else 'Basal')
axes[1,1].set(xticks=positions,xticklabels=labels,xlabel='Volume coefficient / activation gain',ylabel='Full-curve supremum discrepancy (µm)',title='Rational enclosure between discrete curves; not PDE error');axes[1,1].legend(fontsize=8)
fig.suptitle('Matched activated-mesh refinement: same energy, physical inputs, gain and integration rule',fontsize=12)

for ax in axes.ravel(): ax.tick_params(labelsize=11)
path=root/'runs/activated_mesh/activated_mesh_comparison.png';fig.savefig(path,dpi=170);plt.close(fig)
(root/'evidence/activated_mesh_visual.json').write_text(json.dumps({'file':str(path.relative_to(root)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'inputs':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [comparison,results]},'scope':'All6matchedtargetfreenumericalcases; Eulerianfull-domaincurvedifferences; errorbarsarediscretegeometryenclosures, notPDEerrorboundsorconfidenceintervals.'},indent=2));print(str(path))
