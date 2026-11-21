"""Figures from exact virtual-work enclosures and explicit slab constructions."""
from pathlib import Path
import importlib.util
import json
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent


def save(fig,name):
    for suffix in ['png','svg']:
        fig.savefig(HERE/(name+'.'+suffix),dpi=170,
                    metadata={'Creator':'Research Project','Date':None} if suffix=='svg' else {'Software':'Research Project'})
    plt.close(fig)


def main():
    plt.rcParams.update({'figure.facecolor':'#102035','axes.facecolor':'#102035','text.color':'#edf3fc',
                         'axes.labelcolor':'#edf3fc','axes.edgecolor':'#6c809d','xtick.color':'#d0dced',
                         'ytick.color':'#d0dced','font.size':11,'svg.fonttype':'none'})
    arb=json.loads((HERE/'virtual_work_radial_2048.json').read_text())['rows']
    replay=json.loads((HERE/'virtual_work_integer_recentered_4096.json').read_text())['rows']
    fig,ax=plt.subplots(figsize=(9,4.6),constrained_layout=True)
    for k,(rows,label,color,offset) in enumerate([(arb,'Arb: 2,048 cells','#efb96c',.12),
                                                 (replay,'Separate integer replay: 4,096 cells','#79dcd6',-.12)]):
        for j,row in enumerate(rows):
            lo,hi=row['decimal_display_only'];mid=(lo+hi)/2
            ax.errorbar(mid*1e3,j+offset,xerr=np.array([[(mid-lo)*1e3],[(hi-mid)*1e3]]),
                        fmt='o',color=color,capsize=6,label=label if j==0 else None,lw=3)
    ax.axvline(0,color='#e4d2b0',lw=1.3);ax.set_yticks([0,1],['Certified K field','Certified L field'])
    ax.set_xlabel('Reduced relaxed virtual work × 1,000');ax.set_xlim(-1.3,.07)
    ax.set_title('Both old fields admit a descent direction for the candidate relaxed energy')
    ax.grid(alpha=.13);legend=ax.legend(facecolor='#1a314b',edgecolor='#6c809d',loc='lower right')
    for t in legend.get_texts():t.set_color('#edf3fc')
    save(fig,'virtual-work-intervals')

    data=np.genfromtxt(HERE/'slab_profile.csv',delimiter=',',names=True)
    fig,axes=plt.subplots(1,2,figsize=(11.7,4.6),constrained_layout=True)
    axes[0].plot(.8*data['phase']+data['horizontal_correction'],data['vertical_correction'],color='#79dcd6',lw=3)
    axes[0].plot([0,.8],[0,0],ls=':',color='#e4c28c',label='Affine compressed face')
    axes[0].set_aspect('equal',adjustable='datalim');axes[0].set_xlabel('Current horizontal position / period')
    axes[0].set_ylabel('Current height correction / period');axes[0].set_title('One constant-area corrugation cell')
    axes[1].plot(data['phase'],np.cos(data['angle']),color='#79dcd6',lw=2.5,label='Horizontal tangent component')
    axes[1].axhline(.4,ls=':',color='#e4c28c',label='Guaranteed lower bound q/2')
    axes[1].set_xlabel('Material phase');axes[1].set_ylabel('cos(angle)');axes[1].set_ylim(.3,1.05)
    axes[1].set_title('Positive tangent permits an injective bulk extension')
    for ax in axes:
        ax.grid(alpha=.13);lg=ax.legend(facecolor='#1a314b',edgecolor='#6c809d',fontsize=9)
        for t in lg.get_texts():t.set_color('#edf3fc')
    fig.suptitle('Illustrative slab: lambda = (0.8, 1, 1), preferred area = 1; no tissue fit',fontsize=13)
    save(fig,'slab-corrugation')

    bound=json.loads((HERE/'slab_verification.json').read_text())['illustrative_bound_rows']
    b=np.array([float(__import__('fractions').Fraction(x['b_exact'])) for x in bound])
    energy=np.array([x['energy_display'] for x in bound]);position=np.array([x['position_display'] for x in bound])
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),constrained_layout=True)
    axes[0].loglog(b,energy,'o-',color='#efb96c',lw=2);axes[1].loglog(b,position,'o-',color='#79dcd6',lw=2)
    axes[0].set_ylabel('Explicit energy-error upper bound');axes[1].set_ylabel('Uniform position-error upper bound')
    for ax in axes:ax.set_xlabel('Dimensionless bending coefficient b');ax.grid(alpha=.15,which='both')
    axes[0].set_title('Constructive recovery estimate');axes[1].set_title('Same macroscopic affine limit')
    fig.suptitle('Proved upper bounds for the illustrative construction; rate 1/3 is not claimed sharp',fontsize=12)
    save(fig,'slab-recovery-bounds')
    print(json.dumps({'status':'PASS_RESULT_FIGURES','figures':3}))


if __name__=='__main__':main()
