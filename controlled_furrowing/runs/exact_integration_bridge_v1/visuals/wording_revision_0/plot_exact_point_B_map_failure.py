"""Receipt-backed norm comparison: old cube self-map failure at one point."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/exact_integration_bridge_v1/visuals'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=ROOT/'runs/exact_integration_bridge_v1/complete_point_B/results.json';d=json.loads(p.read_text(encoding='utf8'))
    assert d['old_map_center_test']['center_leaves_old_cube_certified']
    r=F(d['old_radius_exact']);eta=F(d['old_uniform_eta_upper_exact'])
    lo=F(d['preconditioned_force_defect_norm_lower_exact']);hi=F(d['preconditioned_force_defect_norm_upper']['upper_exact'])
    intervals=[(F(0),eta),(lo,hi),(max(F(0),lo-eta),hi+eta)]
    plt.rcParams.update({'font.size':12,'text.color':'#e6edf8','axes.labelcolor':'#e6edf8','xtick.color':'#cedcf0','ytick.color':'#cedcf0','axes.edgecolor':'#536b87','figure.facecolor':'#0b111b','axes.facecolor':'#101d2e','savefig.facecolor':'#0b111b'})
    fig,ax=plt.subplots(figsize=(12,5.9))
    labels=['Old cache map at center\ncertified upper bound','Integral/cache correction\ncomplete point norm','New exact-integral map at center\ntriangle-inequality enclosure']
    colors=['#8cd3f4','#a5bdf9','#f2bd6d']
    for y,((a,b),color) in enumerate(zip(intervals,colors)):
        a,b=float(a/r),float(b/r)
        ax.plot([a,b],[y,y],color=color,linewidth=13,solid_capstyle='butt')
        ax.plot([a,b],[y,y],'|',color=color,markersize=17)
        ax.text(b+.035,y,f'[{a:.4f}, {b:.4f}]',color=color,va='center',fontsize=11)
    ax.axvline(1,color='#ef8484',linestyle='--',linewidth=1.7)
    ax.text(1.025,-.48,'Old cube radius',color='#efaaaa',fontsize=11)
    ax.set_yticks(range(3),labels);ax.invert_yaxis();ax.set_xlim(-.025,2.13);ax.set_ylim(2.65,-.65)
    ax.set_xlabel('Displacement norm / old cube radius  (r = 1e-11 in normalized coordinates)')
    ax.grid(axis='x',alpha=.13)
    fig.suptitle('The unchanged proof map sends its center outside the old cube',fontsize=16,y=.955)
    fig.text(.5,.06,'One B-center field. All 144 bulk elements and 24 surface edges included.',ha='center',color='#c9daee')
    fig.text(.5,.02,'This rules out unchanged self-mapping; it does not rule out an equilibrium or provide a continuum error bound.',ha='center',fontsize=10.3,color='#d4c6ad')
    fig.subplots_adjust(left=.35,right=.97,top=.84,bottom=.23)
    OUT.mkdir(parents=True,exist_ok=True);image=OUT/'complete-point-old-map-failure.png'
    if image.exists():raise ValueError('write-once image exists')
    fig.savefig(image,dpi=160);plt.close(fig)
    receipt={'status':'RENDERED_PENDING_VISUAL_REVIEW','source_sha256':sha(Path(__file__)),'point_result_sha256':sha(p),'image_sha256':sha(image),'scope':'Plot of exact rational endpoint ratios, displayed as rounded floating coordinates. Root theorem and independent review remain in separate receipts.'}
    (OUT/'point_map_plot_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    print(json.dumps(receipt))
if __name__=='__main__':main()
