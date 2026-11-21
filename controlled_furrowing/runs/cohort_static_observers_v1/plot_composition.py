"""Plot actual saved uniform intervals; never substitute sampled errors."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent/'composition_final'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf8'))

if __name__=='__main__':
    result=load(OUT/'results.json');rows=result['aggregate']['rows']
    keys=('maximum_shortening_um','region_averaged_shortening_um','apical_max_shortening_um','basal_max_shortening_um')
    titles=('Maximum height shortening','Disk-mean height shortening','Apical downward maximum','Basal upward maximum')
    ids=[r['record_id'] for r in result['records']]
    plt.rcParams.update({'figure.facecolor':'#0d1521','axes.facecolor':'#132134','axes.edgecolor':'#6c819c',
        'text.color':'#eef3ff','axes.labelcolor':'#eef3ff','xtick.color':'#bdcce0','ytick.color':'#bdcce0','font.size':10})
    fig,axes=plt.subplots(1,4,figsize=(16,8),sharey=True)
    for key,title,ax in zip(keys,titles,axes):
        for i,rid in enumerate(ids):
            match=[r for r in rows if r['record_id']==rid and r['observable']==key]
            if not match:ax.text(.02,i,'unresolved',transform=ax.get_yaxis_transform(),color='#ffaf70');continue
            row=match[0];lo,hi=map(Q,row['prediction_interval_exact']);y=Q(row['observed_exact']);mid=(lo+hi)/2
            ax.scatter(float(y),i,marker='x',s=42,color='#ffffff',label='Published value' if i==0 else None,zorder=4)
            ax.errorbar(float(mid),i,xerr=[[float(mid-lo)],[float(hi-mid)]],fmt='o',ms=4.5,capsize=3,
                color='#54d5dc',label='Uniform static interval' if i==0 else None,zorder=3)
            ax.plot([float(y),float(mid)],[i,i],color='#667d9a',lw=.8)
        ax.set_title(title,pad=14);ax.set_xlabel('µm');ax.grid(axis='x',alpha=.15)
    axes[0].set_yticks(range(len(ids)),[s.rsplit('_',1)[-1] for s in ids]);axes[0].invert_yaxis();axes[0].set_ylabel('Published source-column record')
    axes[0].legend(loc='lower left',fontsize=8,facecolor='#132134',labelcolor='#eef3ff')
    fig.suptitle(str(len(rows))+' uniform static enclosures; '+str(len(result['aggregate']['missing_rows']))+' unresolved channels across 17 records',fontsize=16,y=.97)
    lower=result['aggregate']['maximum_residual_lower_exact'];upper=result['aggregate']['maximum_residual_upper_exact']
    scope='All68 covered' if upper is not None else 'Partial coverage: no sufficient upper bound'
    fig.text(.04,.04,scope+'. Bars enclose selected stationary roots across a tiny five-parameter cell; they are not measurement error bars.\n'
        'Static peak proxy, fixed current disk/reference and source-column pairing remain conditional. No common path, stability or PDE certificate.',fontsize=10,color='#c5d1e2')
    fig.subplots_adjust(left=.06,right=.985,top=.88,bottom=.16,wspace=.14)
    outputs=[]
    for suffix in ('png','svg'):
        p=OUT/('cohort-static-observations.'+suffix);fig.savefig(p,dpi=180,bbox_inches='tight');outputs.append(p)
    receipt={'status':'PLOT_GENERATED','generator_sha256':sha(__file__),'results_sha256':sha(OUT/'results.json'),
        'covered_rows':len(rows),'reported_lower_exact':lower,'reported_upper_exact':upper,
        'output_hashes':{p.name:sha(p) for p in outputs},'not_measurement_error_bars':True,'visual_inspection_pending':True}
    target=OUT/'visual_receipt.json'
    if target.exists():raise RuntimeError('retained visual receipt already exists')
    target.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
