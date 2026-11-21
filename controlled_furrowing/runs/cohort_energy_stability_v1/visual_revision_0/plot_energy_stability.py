from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
EXPECTED='5f590ed1bea078e2e88d17bb39f8030d144a95fb73bf96ab10b5823951883164'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

if __name__=='__main__':
    resultfile=OUT/'results.json'
    if sha(resultfile)!=EXPECTED:raise ValueError('declare a new plotting input if producer result changes')
    result=json.loads(resultfile.read_text())
    if result['status']!='CERTIFIED_ALL17_UNIFORM_LOCAL_ENERGY_MINIMA' or len(result['records'])!=17:raise ValueError('actual17 successful producer records required')
    data={};inputs={str(resultfile.relative_to(ROOT)):sha(resultfile)}
    for row in result['records']:
        p=ROOT/row['file']
        if sha(p)!=row['sha256']:raise ValueError('record output changed')
        d=json.loads(p.read_text())
        if d['status']!='CERTIFIED_UNIFORM_STRICT_LOCAL_ENERGY_MINIMUM':raise ValueError('record is not certified by producer')
        q=Q(d['symmetric_part_certificate']['q_infinity_upper_exact'])
        m=Q(d['uniform_free_coordinate_Hessian_lower_bound_exact'])
        if not 0<=q<Q('0.0000000008') or not m>Q('0.000394'):raise ValueError('safe displayed bounds do not hold')
        data[row['record_id']]=(q,m);inputs[row['file']]=sha(p)
    ids=['optoGEF2_small_col_'+chr(i) for i in range(ord('B'),ord('R')+1)]
    if set(data)!=set(ids):raise ValueError('exactB:R source cohort required')
    qvalues=[float(data[i][0]) for i in ids];mvalues=[float(data[i][1])*1e4 for i in ids]
    colors=['#ffba70' if i.rsplit('_',1)[-1] in ('F','K') else '#54d5dc' for i in ids]
    plt.rcParams.update({'figure.facecolor':'#0d1521','axes.facecolor':'#132134','axes.edgecolor':'#6c819c',
        'text.color':'#eef3ff','axes.labelcolor':'#eef3ff','xtick.color':'#bdcce0','ytick.color':'#bdcce0','font.size':10})
    fig,(left,right)=plt.subplots(1,2,figsize=(13,5.8))
    left.scatter(range(17),qvalues,c=colors,s=40,zorder=4)
    left.axhline(1,color='#ee8798',ls='--',lw=1.5,label='Acceptance boundary q = 1')
    left.axhline(8e-10,color='#bdcce0',ls=':',lw=1)
    left.set_yscale('log');left.set_ylim(1e-11,2)
    left.set_title('Exact sym(R) congruence residual upper bounds',pad=14)
    left.set_ylabel('q_i upper bound (log scale)')
    left.text(.5,.68,'Every q_i < 8 × 10⁻¹⁰ < 1',transform=left.transAxes,ha='center',fontsize=12)
    left.legend(loc='upper left',fontsize=8,facecolor='#132134',labelcolor='#eef3ff')
    right.bar(range(17),mvalues,color=colors,width=.65)
    right.axhline(3.94,color='#bdcce0',ls=':',lw=1.2)
    right.set_ylim(0,6.5);right.set_title('Uniform normalized Hessian lower bounds',pad=14)
    right.set_ylabel('m_i lower bound (×10⁻⁴, normalized coordinates)')
    right.text(.5,.93,'Every m_i > 3.94 × 10⁻⁴ > 0',transform=right.transAxes,ha='center',fontsize=12)
    for ax in (left,right):
        ax.set_xticks(range(17),[i.rsplit('_',1)[-1] for i in ids]);ax.set_xlabel('Source-column record');ax.grid(axis='y',alpha=.14)
    fig.suptitle('Producer certificates for 17 selected finite-dimensional energy minima',fontsize=15,y=.98)
    fig.text(.055,.045,'Points/bars are one-sided certificate bounds, not attained norms or eigenvalues. Orange F/K use their verified new proof centers.\\n'
        'Read with the separate independent replay receipt. No continuum, dynamical path, global-minimum or biological stability claim.',fontsize=9.5,color='#c5d1e2')
    fig.subplots_adjust(left=.07,right=.985,top=.84,bottom=.2,wspace=.25)
    outputs=[]
    for suffix in ('png','svg'):
        p=OUT/('cohort-energy-stability.'+suffix);fig.savefig(p,dpi=180,bbox_inches='tight');outputs.append(p)
    receipt={'status':'VISUAL_GENERATED_FROM_PRODUCER_RESULTS','input_hashes':inputs,'generator_sha256':sha(__file__),
        'output_hashes':{p.name:sha(p) for p in outputs},'record_count':17,
        'verified_safe_q_upper_exclusive':'0.0000000008','verified_safe_m_lower_exclusive':'0.000394',
        'independent_mathematical_review_status':'PENDING_SEPARATE_REVIEW','not_measured_eigenvalues_or_physical_rates':True,
        'visual_inspection_pending':True}
    p=OUT/'visual_receipt.json'
    if p.exists():raise ValueError('retain existing visual receipt; do not overwrite')
    p.write_text(json.dumps(receipt,indent=2)+'\\n')
