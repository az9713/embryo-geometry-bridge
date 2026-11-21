from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('controlled_furrowing');out=root/'runs/endpoint_mesh_functionals';obs={n:json.loads((out/f'mesh_{n}_observations.json').read_text())['variants'] for n in (12,16,24)};results=json.loads((root/'evidence/endpoint_mesh_functionals_results.json').read_text())
mid=lambda v:(v['lower']['display']+v['upper']['display'])/2
value=lambda n,w,a,k:mid(obs[n][w+'__'+a]['four_endpoints'][k])
plt.style.use('dark_background');fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained');meshes=list(obs)
for ax,key,label in ((axs[0,0],'apical_downward','Apical component'),(axs[0,1],'basal_upward','Basal component')):
 for w,a,text,c,style in [('current_disk','spatial_max','Spatial maximum','#e9b767','-'),('current_radial_line','regional_mean','Line regional mean','#77cbbb','--'),('current_disk','regional_mean','Disk regional mean','#b5a0e5',':')]:
  ax.plot(meshes,[value(n,w,a,key) for n in meshes],style+'o',color=c,label=text)
 ax.set(title=label,xlabel='Radial cells (vertical cells = half)',ylabel='Signed shortening (micrometres)',xticks=meshes);ax.legend(fontsize=8)
key='maximum_Eulerian_height_shortening';means=[value(n,'current_disk','spatial_max',key) for n in meshes];widths=np.array([[means[i]-obs[n]['current_disk__spatial_max']['four_endpoints'][key]['lower']['display'] for i,n in enumerate(meshes)],[obs[n]['current_disk__spatial_max']['four_endpoints'][key]['upper']['display']-means[i] for i,n in enumerate(meshes)]])
axs[1,0].errorbar(meshes,means,yerr=widths,color='#e9b767',marker='o',capsize=5,label='All four static variants');axs[1,0].set(title='Maximum height shortening',xlabel='Radial cells (vertical cells = half)',ylabel='Micrometres',xticks=meshes);axs[1,0].legend(fontsize=8)
for w,text,c in [('current_radial_line','Line mean','#77cbbb'),('current_disk','Disk mean','#b5a0e5')]:axs[1,1].plot(meshes,[value(n,w,'spatial_max','regional_mean_Eulerian_height_shortening') for n in meshes],'o-',color=c,label=text)
axs[1,1].set(title='Mean height: convention remains visible',xlabel='Radial cells (vertical cells = half)',ylabel='Micrometres',xticks=meshes);axs[1,1].legend(fontsize=8)
fig.suptitle('Saved high-Q gain10 states; fixed current ROI, no fit and no PDE error estimate',fontsize=12);image=out/'endpoint_mesh_values.png';fig.savefig(image,dpi=170,bbox_inches='tight');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(11,4.6),layout='constrained');colors=['#77cbbb','#e9b767'];labels=['12 to 16','16 to 24'];x=np.arange(2)
full12=json.loads((root/'evidence/joint_refinement_comparison.json').read_text())['comparison']['top_Eulerian_vertical']['supremum_absolute_discrete_height_difference_um'];full16=json.loads((root/'evidence/sparse_quadratic_24_comparison.json').read_text())['top']['supremum_absolute_discrete_height_difference_um']
for offset,arr,text,c in [(-.15,[results['comparisons'][i]['surfaces']['top']['supremum_absolute_discrete_height_difference_um'] for i in (0,1)],'Inside fixed current ROI','#77cbbb'),(.15,[full12,full16],'Full common radial domain','#e9b767')]:
 centers=[mid(v) for v in arr];err=np.array([[centers[i]-v['lower']['display'] for i,v in enumerate(arr)],[v['upper']['display']-centers[i] for i,v in enumerate(arr)]])
 axs[0].bar(x+offset,centers,width=.28,color=c,label=text);axs[0].errorbar(x+offset,centers,yerr=err,fmt='none',ecolor='white',capsize=4)
axs[0].set(xticks=x,xticklabels=labels,ylabel='Apical discrete-curve sup (micrometres)',title='ROI sensitivity differs from full-shape sensitivity');axs[0].legend(fontsize=8)
c=results['comparisons'][1];vals=[mid(c['variants']['current_disk__spatial_max'][key]['signed_fine_minus_coarse_um']) for key in ('apical_downward','basal_upward','maximum_Eulerian_height_shortening')]+[mid(c['variants'][w+'__spatial_max']['regional_mean_Eulerian_height_shortening']['signed_fine_minus_coarse_um']) for w in ('current_radial_line','current_disk')]
axs[1].bar(np.arange(5),vals,color=['#e9b767','#77cbbb','#c995ba','#b5a0e5','#8ca8de']);axs[1].axhline(0,color='white',lw=.7);axs[1].set(xticks=np.arange(5),xticklabels=['Apical\nmax','Basal\nmax','Height\nmax','Height\nline mean','Height\ndisk mean'],ylabel='24 minus16 endpoint (micrometres)',title='Small total change can hide component changes')
fig.suptitle('Same constitutive scenario; numerical mesh/branch sensitivity only',fontsize=12);image2=out/'endpoint_mesh_sensitivity.png';fig.savefig(image2,dpi=170,bbox_inches='tight');plt.close(fig)
receipt={'scope':'Original figures of saved static discrete fields; no biological calibration or PDE error bars.','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'results_sha256':hashlib.sha256((root/'evidence/endpoint_mesh_functionals_results.json').read_bytes()).hexdigest(),'images':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (image,image2)}};(root/'evidence/endpoint_mesh_functionals_visuals.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(json.dumps(receipt,indent=2))
