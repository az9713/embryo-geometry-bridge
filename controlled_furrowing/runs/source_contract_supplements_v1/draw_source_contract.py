"""Original audit schematic only; no source images or publication adaptations."""
from pathlib import Path
import hashlib,json,textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch

OUT=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(16,10),dpi=150);fig.patch.set_facecolor('#0d141e');ax.set_facecolor('#0d141e');ax.set(xlim=(0,16),ylim=(0,10));ax.axis('off')
white='#edf4fc';muted='#afbed0';cyan='#4ddbd1';amber='#f0b76a';blue='#172738'
def text(x,y,s,size=12,c=white,weight='normal',width=None):
    if width:s='\n'.join(textwrap.fill(line,width=width)for line in s.split('\n'))
    ax.text(x,y,s,fontsize=size,color=c,fontweight=weight,va='top',linespacing=1.5)
def card(x,y,w,h,title,body,accent=cyan):
    ax.add_patch(FancyBboxPatch((x,y-h),w,h,boxstyle='round,pad=.015,rounding_size=.13',facecolor=blue,edgecolor=accent,lw=1.3))
    text(x+.22,y-.22,title,15,accent,'bold');text(x+.22,y-.8,body,11.3,width=44)
text(.35,9.73,'From public source files to a defensible geometry test',23,white,'bold')
text(.35,9.22,'Countryman et al. (2025): source audit and observation bridge',12,muted)
card(.35,8.65,4.65,3.15,'1  Verified public package','16 official files / 111,690,608 bytes\n4 PDF + 2 XLSX + 10 AVI\n40 visible worksheets\n1,000 formula-cell elements\n922 metadata bytes remain undecoded')
card(5.65,8.65,4.65,3.15,'2  Inspected image content','903 encoded frames fully decoded\nSelected frames visually inspected\nCaptions, display clocks and scale bars\nSlices / projections / activation overlays\nPrepared views, not full raw 3D volumes')
card(10.95,8.65,4.65,3.15,'3  Needed measurement bridge','Specimen → dose → acquisition IDs\nRegistered apical / basal boundaries\nTracked response ROI + growth reference\nComponent reference / peak rules\nJustified measurement uncertainty',amber)
for x in [5.05,10.35]:
    ax.add_patch(FancyArrowPatch((x,7.15),(x+.5,7.15),arrowstyle='-|>',mutation_scale=19,color=muted,lw=1.6))
text(.45,5.12,'What the source already specifies',16,cyan,'bold')
text(.45,4.65,'32–35 µm measured cell-length selection; 40 s scan;\n2 µm z sampling and 22 s/frame for cited recordings.\nSampling pitch is not a certified contour error.\nScan-end cannot be assumed undeformed.',12,white)
text(8.45,5.12,'What the mathematics currently certifies',16,cyan,'bold')
text(8.45,4.65,'17 certified local energy minima + 68 conditional static outputs.\nOne tiny shared parameter cell; a frozen discrete energy.\nDiscrepancy remains approximately 5.186387 µm.\nThis does not identify the experimental observation operator.',12,white)
ax.plot([.45,15.55],[3.03,3.03],color='#38506b',lw=1)
text(.45,2.75,'The next useful input is a compact processed-data crosswalk',17,amber,'bold')
text(.45,2.22,'specimens  →  frames / activation times  →  registered boundaries + ROI  →  reference / summaries',13,white,'bold')
text(.45,1.67,'Request existing tables, masks and workflow settings. No new experiment is required by this request.',12,muted)
text(.45,1.12,'File integrity ≠ anatomical validity     •     numerical proof width ≠ measurement uncertainty',13,amber)
text(.45,.62,'Original schematic; no source microscopy pixels. No new fit, noise allowance, novelty claim or completed five-gate claim.',10.5,muted)
fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
for suffix in ['png','svg']:fig.savefig(OUT/f'public-source-contract.{suffix}',facecolor=fig.get_facecolor(),metadata={'Creator':'Original source-contract schematic'})
plt.close(fig)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
receipt={'status':'GENERATED_ORIGINAL_SCHEMATIC_VISUAL_INSPECTION_PENDING','source_sha256':sha(Path(__file__)),'files':{f'public-source-contract.{s}':sha(OUT/f'public-source-contract.{s}')for s in ['png','svg']},'source_image_pixels_used':False,'public_source_adaptation':False,'meaning':'Original audit/dataflow schematic. Counts are source-audit facts; it contains no empirical contour or uncertainty claim.','draft_sha256_at_generation':sha(OUT/'public-source-audit-draft.md')}
(OUT/'original_visual_generation.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
print(json.dumps(receipt))
