"""Plot reviewed cached derivative bounds; no spatial-integral claim."""
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'runs/transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1'
INPUT = BASE/'cached_derivative_extension_v1/independent_actual_preparation/bridge_admission.json'
EXPECTED = '269229312793ec101cc301086d634dd4de89c62a5050ccfc60f185cf761dbc7d'


def main():
    assert sha256(INPUT.read_bytes()).hexdigest() == EXPECTED
    data = json.loads(INPUT.read_text())
    bounds = [Fraction(r['reviewed_ZQ_upper_exact']) for r in data['reviewed_cases']]
    assert len(bounds) == 4 and all(0 <= z < 1 for z in bounds)
    labels = ['Candidate + | small L','Candidate + | large','Candidate - | small L','Candidate - | large']
    bg, fg = '#0b1320', '#eaf3fc'
    plt.rcParams.update({'font.size':11,'text.color':fg,'axes.labelcolor':fg,
                         'xtick.color':fg,'ytick.color':fg,'font.family':'DejaVu Sans'})
    fig, ax = plt.subplots(figsize=(10.5,5.2),facecolor=bg)
    ax.set_facecolor(bg)
    for i,z in enumerate(bounds):
        ax.barh(i,float(z),color='#77dcf4',height=.54)
        ax.barh(i,float(1-z),left=float(z),color='#243347',height=.54)
        ax.text(float(z)+.018,i,f'{float(z):.6f}',va='center',color=fg)
    ax.set_yticks(range(4),labels)
    ax.invert_yaxis()
    ax.set_xlim(0,1.03)
    ax.set_xticks([0,.25,.5,.75,1])
    ax.axvline(1,color='#ffca84',ls='--',lw=1.5)
    ax.set_xlabel('Upper bound on the cached derivative defect (dimensionless)')
    for spine in ax.spines.values():
        spine.set_color('#3b526e')
    fig.suptitle('Larger cubes retain contraction headroom',fontsize=17,x=.04,ha='left')
    fig.text(.04,.89,'Independent whole-box replay at radius 2^-34; four frozen targets',color='#a6b8ce')
    fig.text(.04,.075,'Blue: reviewed cached bound. Dark: remaining room below one.\n'
             'The spatial Hessian defect and point-force test are still pending.',color='#a6b8ce',fontsize=10)
    fig.subplots_adjust(left=.25,right=.96,top=.8,bottom=.23)
    out = BASE/'certification_checkpoint_v1'
    figure = out/'cached_derivative_headroom.png'
    with figure.open('xb') as handle:
        fig.savefig(handle,format='png',dpi=165,facecolor=bg)
    plt.close(fig)
    receipt = {'status':'PLOT_OF_REVIEWED_CACHED_BOUNDS_ONLY','input_sha256':EXPECTED,
               'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
               'figure_sha256':sha256(figure.read_bytes()).hexdigest(),
               'plotted_bounds_exact':list(map(str,bounds)),
               'floating_values_used_only_for_display':True,'new_spatial_integral_root_claimed':False}
    (out/'figure_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'figure':str(figure),'figure_sha256':receipt['figure_sha256']}))


if __name__ == '__main__':
    main()
