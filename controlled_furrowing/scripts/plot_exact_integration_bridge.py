"""Receipt-backed plot of local integration discrepancies, not physical errors."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/exact_integration_bridge_v1/visuals'
def read(p):return json.loads(p.read_text(encoding='utf8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def upper(v):return float(F(v['upper_exact']))
def main():
    paths=[ROOT/'evidence/exact_bulk_pilot_independent_review.json',ROOT/'runs/exact_integration_bridge_v1/bulk_pilot/results.json']
    review=read(paths[0]);assert review['status'].startswith('PASS') and all(review['checks'].values())
    results=[]
    for element in (0,6,72,66):
        p=ROOT/f'runs/exact_integration_bridge_v1/bulk_pilot/element_{element:03d}.json';paths.append(p);results.append(read(p))
    plt.rcParams.update({'font.size':12,'text.color':'#e8eef9','axes.labelcolor':'#e8eef9','axes.edgecolor':'#59708e','xtick.color':'#d0def1','ytick.color':'#d0def1','figure.facecolor':'#0b111b','axes.facecolor':'#101d2e','savefig.facecolor':'#0b111b'})
    fig,axes=plt.subplots(1,2,figsize=(13.5,6.4));x=np.arange(4)
    for ax,kind,ylabel in zip(axes,('force','H'),('Maximum local force-entry discrepancy','Maximum local Hessian-entry discrepancy')):
        total=[upper(r['max_force_defect_upper' if kind=='force' else 'max_H_entry_defect_upper']) for r in results]
        tail=[upper(r['force_max_tail_upper' if kind=='force' else 'H_max_tail_upper']) for r in results]
        ax.semilogy(x,total,'o',color='#82d6fc',markersize=9,label='Total certified absolute upper bound')
        ax.semilogy(x,tail,'D',color='#f2bd6d',markersize=7,label='Included series-remainder upper bound')
        for i,(a,b) in enumerate(zip(total,tail)):
            ax.plot([i,i],[b,a],color='#48647c',linewidth=1.3)
        ax.set_xticks(x,['Axis vertex\nT0','Near axis\nT6','Axis edge\nT72','Outer field\nT66'])
        ax.set_ylabel(ylabel+'\n(normalized energy coordinates)');ax.grid(axis='y',alpha=.15);ax.set_xlim(-.4,3.4)
        ax.legend(loc='best',fontsize=9,framealpha=.15,labelcolor='#e8eef9')
    fig.suptitle('Exactly integrated FE energy versus frozen binary quadrature cache',fontsize=18,x=.5,y=.965)
    fig.text(.5,.905,'Four bulk elements at ONE saved B-center field; order 16, 512-bit outward arithmetic',ha='center',fontsize=12,color='#b6cce8')
    fig.text(.5,.065,'Independent order-18 / 640-bit replay: all 38 force and 210 Hessian entries overlap.',ha='center',fontsize=12,color='#a2e5bc')
    fig.text(.5,.025,'Local numerical discrepancies only: surfaces, assembled fields, parameter tubes and root transfer are separate obligations.',ha='center',fontsize=10.5,color='#c5cddd')
    fig.subplots_adjust(left=.085,right=.98,bottom=.23,top=.83,wspace=.29)
    OUT.mkdir(parents=True,exist_ok=True);image=OUT/'bulk-element-integration-defects.png'
    if image.exists():raise ValueError('write-once image exists')
    fig.savefig(image,dpi=160);plt.close(fig)
    receipt={'status':'RENDERED_PENDING_VISUAL_REVIEW','source_sha256':sha(Path(__file__)),'input_hashes':{str(p.relative_to(ROOT)):sha(p) for p in paths},'image':str(image.relative_to(ROOT)),'image_sha256':sha(image),'scope':'Plot of rigorous absolute upper bounds from receipts; floating plot coordinates are illustrative. No physical measurement error or novelty claimed.'}
    (OUT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    print(json.dumps(receipt))
if __name__=='__main__':main()
