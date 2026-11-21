from pathlib import Path
from fractions import Fraction as Q
import sys,json,hashlib
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from run_activated_mesh_study import QuadraticGraph,compare_graphs
rows=[]
for k in range(1,13):
 a=(Q(k%3),Q(k%5-2),Q(k%7-3));b=(Q((k+1)%3),Q((k+2)%5-2),Q((k+3)%7-3));d=tuple(x-y for x,y in zip(a,b))
 def graph(c):return QuadraticGraph([[(0,c[0]),(Q(1,2),c[0]+c[1]/2+c[2]/4),(1,sum(c))]],Q(1))
 candidates=[Q(0),Q(1)]
 if d[2]!=0 and 0<-d[1]/(2*d[2])<1:candidates.append(-d[1]/(2*d[2]))
 true=max(abs(d[0]+d[1]*r+d[2]*r*r) for r in candidates)
 result=compare_graphs(graph(a),graph(b),tolerance_um=Q(1,50),maximum_cells=512);bound=result['supremum_absolute_discrete_height_difference_um']
 rows.append({'case':k,'true_exact_supremum':str(true),'lower':str(bound.lo),'upper':str(bound.hi),'contains_known_truth':bound.lo<=true<=bound.hi,'achieved_requested_width':result['requested_width_achieved']})
receipt={'status':'PASS' if all(r['contains_known_truth'] for r in rows) else 'FAIL','case_count':len(rows),'rows':rows,'source_sha256':hashlib.sha256(Path('controlled_furrowing/scripts/run_activated_mesh_study.py').read_bytes()).hexdigest(),'scope':'Independent closedform maxima of quadratic graph differences vs interval-cell algorithm; no mechanical/data inputs.'}
Path('controlled_furrowing/evidence/activated_mesh_graph_adversarial.json').write_text(json.dumps(receipt,indent=2));print(json.dumps({'status':receipt['status'],'count':len(rows),'all_widths_achieved':all(r['achieved_requested_width'] for r in rows)},indent=2))
