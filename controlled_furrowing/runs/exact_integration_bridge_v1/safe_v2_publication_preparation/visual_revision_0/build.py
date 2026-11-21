"""Original source-bound figure and checks; no scientific replay or production edits."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,datetime,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):
    if p.exists():raise ValueError('preserve '+str(p))
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n')
def main():
    sources={
      'evidence/interval_polynomial_cleanup_audit.json':'fe0f5414ecf16b6ff3e804d53a6ad68d44ef6a6766824786acd4ef2645fff084',
      'evidence/exact_safe_polynomial_kernel_independent_review.json':'a2f6b6c4740e5475a818cf58320fa43cc8c8ea63c7f37afe44c1c8632bd9356d',
      'evidence/exact_bulk_independent_replay_safe_v2.json':'37515a8d8c9591bb656bc7ee06e6c9cdd94393088ce520edf1fbdb158d413af2',
      'evidence/exact_surface_safe_independent_review.json':'248ba5000b27f8e29ede2b501b00bf16e5fbdc6624f06a65a8ced339c39c37a2',
      'evidence/exact_complete_point_B_independent_review_safe_v2.json':'433125c73e1ee6856f425ab7dc1d62658d9505a7c7ff8968c5ad2d4f7363dfc2',
      'runs/exact_integration_bridge_v1/bulk_independent_replay_safe_v2/coordinate_witness.json':'98a70c8457537bda418fb94152bd5ca0536cf1731dd83a5d973d733f9c6a3c05',
      'runs/exact_integration_bridge_v1/bulk_independent_replay_safe_v2/delivery.json':'68383583a2116b9f82cff9e98746cafd6a5c878935d4e26cfd8192fcfddfa117',
      'evidence/cohort_energy_stability_independent_review.json':'d4f4c4d30d16bbbf2141cabf1cee18b47fcc35ab958a2a8a203000f0a7acaf95',
      'evidence/exact_retube_B_geometry_independent_review.json':'8114999e0dd022fb00cb8aa95040e8bf6096c636fea4978e4bb7cc40c6f40002',
      'scripts/exact_fe_target_safe.py':'889b61cfa794d432021fbd4670b781169033d4b8451554b18e7f7a569d9d3954',
      'runs/exact_integration_bridge_v1/safe_kernel_v2/polynomial_jet_safe.py':'1439e20ec127d39565467d0788173b1148b5ab47f659b3b39486b81aaf1ee373'}
    # The public-facing draft uses readable shortened hashes; this ledger binds
    # complete measured hashes.
    kernelpath='evidence/exact_safe_polynomial_kernel_independent_review.json'
    for p,h in sources.items():assert sha(ROOT/p)==h,p
    reduction=read(ROOT/'evidence/exact_complete_point_B_independent_review_safe_v2.json');c=reduction['certificate'];w=read(ROOT/'runs/exact_integration_bridge_v1/bulk_independent_replay_safe_v2/coordinate_witness.json');bulk=read(ROOT/'evidence/exact_bulk_independent_replay_safe_v2.json');surface=read(ROOT/'evidence/exact_surface_safe_independent_review.json');kernel=read(ROOT/kernelpath);old=read(ROOT/'evidence/cohort_energy_stability_independent_review.json')
    r=Q(c['old_radius_exact']);eta=Q(c['eta_used_exact']);lo=Q(c['RF_norm_lower_exact']);hi=Q(c['RF_norm_upper_exact']);maplo=Q(c['new_integral_map_center_distance_lower_exact']);shift=tuple(map(Q,w['Tintegral_u0_minus_u0_interval_exact']))
    checks={'new208_checks_literal_true':len(reduction['checks'])==208 and all(v is True for v in reduction['checks'].values()),'safe992_kernel_checks_literal_true':len(kernel['checks'])==992 and all(v is True for v in kernel['checks'].values()),'six_bulk32_checks_literal_true':len(bulk['checks'])==32 and all(v is True for v in bulk['checks'].values()),'six_surface349_checks_literal_true':len(surface['checks'])==349 and all(v is True for v in surface['checks'].values()),'all168_local_packets':len(reduction['all168_packets'])==168,'old_radius_exact_1e_minus11':r==Q(1,10**11),'positive_eta_and_exact_reverse_triangle':eta>0 and maplo==max(Q(0),lo-eta),'new_map_bound_strictly_exceeds_radius':maplo>r,'axial_coordinate_witness_outside_projection':w['free_index']==588 and w['component']=='axial'and shift[1]<-r,'witness_lower_abs_agrees_with_norm_proof':-shift[1]==maplo,'outward_readable_discrepancy_bounds':Q('1.673692093e-11')<=lo<=hi<=Q('1.673692097e-11'),'outward_readable_eta_bound':eta<=Q('5.287455044e-13'),'outward_readable_map_lower_bound':Q('1.6208175e-11')<=maplo,'old17_denseQ_energy_receipt_unchanged':old['status']=='PASS'and len(old['record_audits'])==17,'no_uniform_root_or_PDE_restoration':reduction['old_v1_claims_restored']is False,'production_files_not_written':True}
    assert all(v is True for v in checks.values()),checks
    palette={'bg':'#0d141d','panel':'#152230','ink':'#e5edf5','muted':'#a9bbcd','teal':'#65d7bd','gold':'#ffc96c','red':'#ff9398'}
    plt.rcParams.update({'font.family':'DejaVu Sans','text.color':palette['ink'],'axes.labelcolor':palette['ink'],'xtick.color':palette['muted'],'ytick.color':palette['muted'],'svg.fonttype':'none'})
    fig=plt.figure(figsize=(14,8.4),facecolor=palette['bg']);gs=fig.add_gridspec(2,2,height_ratios=[3,1.65],width_ratios=[1.5,1],left=.065,right=.96,top=.80,bottom=.12,wspace=.25,hspace=.53)
    fig.text(.065,.94,'A corrected point certificate',fontsize=25,weight='bold');fig.text(.065,.887,'The new integral map moves the old center beyond its certified cube.',fontsize=14,color=palette['muted'])
    ax=fig.add_subplot(gs[0,0]);ax.set_facecolor(palette['panel']);ax.axvspan(-1,1,color=palette['teal'],alpha=.13);ax.plot([-1,1],[.76,.76],color=palette['teal'],linewidth=12,solid_capstyle='butt');ax.text(0,.93,'Old cube projected onto coordinate588',ha='center',fontsize=11,color=palette['teal'])
    left,right=[float(x/r)for x in shift];ax.plot([left,right],[.30,.30],color=palette['gold'],linewidth=12,solid_capstyle='butt');ax.scatter([0],[.30],color=palette['ink'],s=40,zorder=5);ax.text(.02,.43,'Old center',fontsize=10);ax.text((left+right)/2,.08,'New map image interval',ha='center',fontsize=10,color=palette['gold']);ax.axvline(-1,color=palette['teal'],linestyle='--',alpha=.7);ax.axvline(1,color=palette['teal'],linestyle='--',alpha=.7);ax.axvline(0,color=palette['muted'],alpha=.2)
    ax.annotate('',xy=(-1,.54),xytext=(right,.54),arrowprops={'arrowstyle':'<->','color':palette['red'],'lw':1.6});ax.text((right-1)/2,.63,'Guaranteed gap',ha='center',fontsize=10,color=palette['red']);ax.set_xlim(-2.05,1.2);ax.set_ylim(-.02,1.13);ax.set_xticks([-2,-1,0,1]);ax.set_yticks([]);ax.set_xlabel(r'$(T_I(u_0)-u_0)_{588}/r$,  $r=10^{-11}$',labelpad=12,fontsize=11);ax.set_title('A coordinate witness; no root location is shown',loc='left',fontsize=12,pad=17)
    for spine in ax.spines.values():spine.set_visible(False)
    bar=fig.add_subplot(gs[0,1]);bar.set_facecolor(palette['panel']);bar.barh([1,0],[1,float(maplo/r)],height=.36,color=[palette['teal'],palette['gold']]);bar.set_yticks([1,0],['Old radius','New map\nlower bound'],fontsize=11);bar.set_xlim(0,2.1);bar.set_ylim(-.65,1.65);bar.set_xticks([0,1,2]);bar.axvline(1,color=palette['teal'],linestyle='--',alpha=.7);bar.text(1.04,1,'1r',va='center',fontsize=12,color=palette['teal']);bar.text(float(maplo/r)+.045,0,'≥1.6208175r',va='center',fontsize=11,color=palette['gold']);bar.set_xlabel('Normalized displacement / old radius',labelpad=12,fontsize=11);bar.set_title('Reverse-triangle bound',loc='left',fontsize=12,pad=17)
    for spine in bar.spines.values():spine.set_visible(False)
    chain=fig.add_subplot(gs[1,:]);chain.set_axis_off();cards=[('Fresh SAFE integrals','144bulk +24surface\n6 +6 independently replayed'),('Exact signed reduction','611F /6676H /611RδF\n208 checks,355bindings'),('Valid old-Q residual','Subtract max of author\nand independent η'),('Restricted conclusion','Map center exits old cube\nNo root-exclusion theorem')]
    for i,(title,body)in enumerate(cards):
        x=i*.255;patch=FancyBboxPatch((x,.24),.225,.69,boxstyle='round,pad=0.012',transform=chain.transAxes,facecolor=palette['panel'],edgecolor=palette['muted'],linewidth=.6);chain.add_patch(patch);chain.text(x+.1125,.74,title,ha='center',va='center',fontsize=11,weight='bold');chain.text(x+.1125,.47,body,ha='center',va='center',fontsize=10,color=palette['muted'],linespacing=1.5)
        if i<3:chain.annotate('',xy=(x+.247,.585),xytext=(x+.229,.585),xycoords='axes fraction',arrowprops={'arrowstyle':'->','color':palette['ink'],'lw':1.4})
    fig.text(.065,.055,'One original B field and base parameter. No uniform integration bound, new integral equilibrium, PDE or biological validation.',fontsize=10.5,color=palette['muted'])
    for suffix in ['png','svg']:
        out=HERE/f'safe_point_restoration.{suffix}';assert not out.exists();fig.savefig(out,dpi=180,facecolor=fig.get_facecolor())
    plt.close(fig)
    ledger={'status':'PASS_DRAFT_NUMERICS_AND_ORIGINAL_FIGURE','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'source_sha256':sha(Path(__file__)),'source_hashes':sources,'exact_values':{'r':str(r),'eta':str(eta),'RF_lower':str(lo),'RF_upper':str(hi),'new_map_lower':str(maplo),'map_to_radius_lower_ratio':str(maplo/r),'coordinate588_map_shift':list(map(str,shift))},'figure_hashes':{s:sha(HERE/f'safe_point_restoration.{s}')for s in ['png','svg']},'visual_values_are_approximations':True,'no_math_replay_or_new_solver_run':True,'production_files_changed':False,'parent_review_and_publication_pending':True}
    write(HERE/'receipt.json',ledger);print(json.dumps({'status':ledger['status'],'checks':len(checks),'receipt_sha256':sha(HERE/'receipt.json')}))
if __name__=='__main__':main()
