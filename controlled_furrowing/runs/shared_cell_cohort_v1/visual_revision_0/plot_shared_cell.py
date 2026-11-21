"""Visualize the joined68 exact uniform output intervals; no model recomputation."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
RESULT=ROOT/'evidence/shared_cell_cohort_inference.json'
EXPECTED='1f8c4dc1ad8924f8fda13e6748d5e73aac426b5771b02404803d27c1955aa677'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

if __name__=='__main__':
    if sha(RESULT)!=EXPECTED:raise ValueError('joined inference packet changed; declare a new plotting input')
    d=json.loads(RESULT.read_text());rows=d['uniform_constraints']
    ids=['optoGEF2_small_col_'+chr(c) for c in range(ord('B'),ord('R')+1)]
    keys=('maximum_shortening_um','region_averaged_shortening_um','apical_max_shortening_um','basal_max_shortening_um')
    if len(rows)!=68 or {(r['record_id'],r['observable']) for r in rows}!={(i,k) for i in ids for k in keys}:raise ValueError('actual17by4coverage required')
    plt.rcParams.update({'figure.facecolor':'#0d1521','axes.facecolor':'#132134','axes.edgecolor':'#6c819c',
        'text.color':'#eef3ff','axes.labelcolor':'#eef3ff','xtick.color':'#bdcce0','ytick.color':'#bdcce0','font.size':10})
    outputs=[];fig,axes=plt.subplots(1,4,figsize=(16,8),sharey=True)
    for key,title,ax in zip(keys,('Maximum height shortening','Disk-mean height shortening','Apical downward maximum','Basal upward maximum'),axes):
        for i,rid in enumerate(ids):
            row=next(r for r in rows if r['record_id']==rid and r['observable']==key)
            lo,hi=map(Q,row['prediction_interval_exact']);y=Q(row['observed_exact']);mid=(lo+hi)/2
            new=rid.rsplit('_',1)[-1] in ('F','K');color='#ffba70' if new else '#54d5dc'
            ax.scatter(float(y),i,marker='x',s=42,color='#fff',label='Published value' if i==0 else None,zorder=4)
            ax.errorbar(float(mid),i,xerr=[[float(mid-lo)],[float(hi-mid)]],fmt='o',ms=4.5,capsize=3,color=color,
                label='Uniform enclosure' if i==0 else 'F/K corrected proof centers' if i==4 else None,zorder=3)
            ax.plot([float(y),float(mid)],[i,i],color='#667d9a',lw=.8)
        ax.set_title(title,pad=14);ax.set_xlabel('µm');ax.grid(axis='x',alpha=.15)
    axes[0].set_yticks(range(17),[s.rsplit('_',1)[-1] for s in ids]);axes[0].invert_yaxis();axes[0].set_ylabel('Published source-column record')
    axes[0].legend(loc='lower left',fontsize=7.6,facecolor='#132134',labelcolor='#eef3ff')
    fig.suptitle('17 shared-cell stationary roots; 68 conditional static output enclosures',fontsize=16,y=.97)
    fig.text(.04,.04,'F/K are now filled using fresh root proofs and an explicit offset back to the ORIGINAL seed graphs. No physical parameter was changed.\n'
        'Bars are uniform mathematical enclosures on one tiny computational cell, not measurement error bars. Static source/ROI assumptions remain conditional.',fontsize=10,color='#c5d1e2')
    fig.subplots_adjust(left=.06,right=.985,top=.88,bottom=.16,wspace=.14)
    for suffix in ('png','svg'):
        p=OUT/('shared-cell-68-observations.'+suffix);fig.savefig(p,dpi=180,bbox_inches='tight');outputs.append(p)
    plt.close(fig)
    lo=Q(d['maximum_residual_lower_exact']);hi=Q(d['maximum_residual_upper_exact']);l=float(lo);h=float(hi)
    fig,ax=plt.subplots(figsize=(13,5));ax.set_xlim(5.07,5.23);ax.set_ylim(0,1);ax.set_yticks([])
    ax.axvspan(5.07,l,color='#ba5266',alpha=.35);ax.axvspan(h,5.23,color='#2b9693',alpha=.4)
    ax.axvline(l,color='#ffba70',lw=2);ax.axvline(h,color='#ffba70',lw=2)
    ax.scatter([5.1,5.2],[.25,.25],s=70,color=['#ff8f9f','#70e0d8'],zorder=4)
    ax.text(5.1,.34,'δ = 5.1\nselected tuple OUT\nother branches unresolved',ha='center',va='bottom',fontsize=11)
    ax.text(5.205,.34,'δ = 5.2\nC is an inner stationary-tuple witness',ha='center',va='bottom',fontsize=11)
    ax.set_xlabel('Declared common discrepancy allowance δ (µm)');ax.set_title('A validated residual enclosure, not a measured noise level',pad=18,fontsize=15)
    ax.text(.5,.89,'For every θ in C:   5.186387141964336 ≤ D(θ) ≤ 5.186387360305146 µm',transform=ax.transAxes,ha='center',fontsize=12)
    ax.text(.5,.75,'The narrow interval is unresolved by these bounds; its width is approximately 2.184 × 10⁻⁷ µm.\n'
        'Endpoints bound the residual; they are not certified attained minima or maxima.',transform=ax.transAxes,ha='center',fontsize=10,color='#d4dfed')
    fig.text(.07,.02,'All68 static constraints on one tiny five-parameter cell. No noise/stability/common-time-path/PDE/full-domain/novelty claim.',fontsize=10,color='#c5d1e2')
    fig.subplots_adjust(bottom=.19,top=.84,left=.07,right=.97)
    for suffix in ('png','svg'):
        p=OUT/('shared-cell-residual-threshold.'+suffix);fig.savefig(p,dpi=180,bbox_inches='tight');outputs.append(p)
    receipt={'status':'PLOT_GENERATED','generator_sha256':sha(__file__),'inference_sha256':sha(RESULT),
        'covered_records':17,'covered_channels':68,'F_K_recentered_labels_explicit':True,
        'residual_bounds_exact':[str(lo),str(hi)],'output_hashes':{p.name:sha(p) for p in outputs},
        'no_measurement_error_or_attained_extrema_claim':True,'visual_inspection_pending':True}
    file=OUT/'visual_receipt.json'
    if file.exists():raise ValueError('retain existing visual receipt; do not overwrite')
    file.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
