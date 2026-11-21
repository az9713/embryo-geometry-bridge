from pathlib import Path
import json,hashlib,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('controlled_furrowing');receipt=json.loads((root/'evidence/sparse_quadratic_24_comparison.json').read_text())
plt.style.use('dark_background');fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for j,(surface,label) in enumerate([('top','Apical'),('bottom','Basal')]):
 q=receipt[surface]['sampled_eulerian_curves'];r=[float(x['radius_um']) for x in q]
 axs[0,j].plot(r,[float(x['coarse_z_um']) for x in q],'--',color='#dbae66',label='16 x 8, high Q');axs[0,j].plot(r,[float(x['fine_z_um']) for x in q],color='#7bc7b4',label='24 x 12, high Q');axs[0,j].set(xlabel='Current physical radius (micrometres)',ylabel='Height (micrometres)',title=label+' full Eulerian graph');axs[0,j].legend(fontsize=8)
 axs[1,j].plot(r,[(float(x['difference_lower_um'])+float(x['difference_upper_um']))/2 for x in q],color='#b3a0df');axs[1,j].axhline(0,color='grey',lw=.5)
 bound=receipt[surface]['supremum_absolute_discrete_height_difference_um'];print(bound)
 # The outward-rational bounds are already computed; this is display only.
 axs[1,j].set(xlabel='Current physical radius (micrometres)',ylabel='24 minus16 height (micrometres)',title=['Apical stored-curve sup: [0.6709, 0.6803] micrometres','Basal stored-curve sup: [0.1235, 0.1331] micrometres'][j])
fig.suptitle('Same high-quadrature constitutive model and gain10; mesh/branch sensitivity, not a PDE error estimate',fontsize=10)
image=root/'runs/sparse_quadratic/24_vs16_highQ.png';fig.savefig(image,dpi=170,bbox_inches='tight');plt.close(fig)
(root/'evidence/sparse_quadratic_24_comparison_visual.json').write_text(json.dumps({'image':str(image.relative_to(root)),'sha256':hashlib.sha256(image.read_bytes()).hexdigest(),'comparison_sha256':hashlib.sha256((root/'evidence/sparse_quadratic_24_comparison.json').read_bytes()).hexdigest(),'display_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'note':'Display-only replay of saved comparison; panel bounds rounded outward. Visual inspected after title/canvas adjustment.'},indent=2),encoding='utf-8')
