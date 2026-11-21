from pathlib import Path
from fractions import Fraction as Q
import sys,json,hashlib
import numpy as np
sys.path.insert(0,str(Path('controlled_furrowing/scripts').resolve()))
from robust_quadratic import build_model,Parameters
from annulus_observer import serial
root=Path('controlled_furrowing');rows=[];values={}
for mesh in ['quadratic_robust','activated_mesh']:
 for lam in [500.,5000.]:
  for gain in [1.,5.,10.]:
   p=root/f'runs/{mesh}/lambda_{lam}_target_{gain}/final.json';r=json.loads(p.read_text());m=build_model(Parameters(**r['parameters']));assert np.array_equal(m.reference,np.array(r['reference_nodes']))
   top=r['boundary_top_quadratic_controls'];bottom=r['boundary_bottom_quadratic_controls'];a=Q(r['parameters']['patch_radius']);heights=[]
   for i,(t,b) in enumerate(zip(top,bottom)):
    if Q(m.reference[m.top_edges[i,2],0])>a:continue
    q=[Q(t[k][1])-Q(b[k][1]) for k in range(3)];c0=q[0];c1=-3*q[0]+4*q[1]-q[2];c2=2*q[0]-4*q[1]+2*q[2];ts=[Q(0),Q(1)]
    if c2!=0 and 0<-c1/(2*c2)<1:ts.append(-c1/(2*c2))
    heights += [c0+c1*x+c2*x*x for x in ts]
   s=(Q(1)-min(heights))*Q('33.5');values[(mesh,lam,gain)]=s
   rows.append({'mesh':mesh,'lambda':lam,'gain':gain,'exact_stored_geometry_max_shortening_um':s,'reported_float_difference_um':abs(float(s)-r['maximum_activated_column_shortening']*33.5),'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
orders=[]
for gain in [1.,5.,10.]:
 c=values[('quadratic_robust',500.,gain)]-values[('quadratic_robust',5000.,gain)];f=values[('activated_mesh',500.,gain)]-values[('activated_mesh',5000.,gain)]
 orders.append({'gain':gain,'coarse_k500_minus_k5000_um':c,'fine_k500_minus_k5000_um':f,'strict_order_reversed':c*f<0})
checks={'all12_reported_float_maxima_match_exact_within_1e_minus10_um':all(x['reported_float_difference_um']<1e-10 for x in rows),'exact_order_reversed_at_gains5and10':all(x['strict_order_reversed'] for x in orders if x['gain'] in [5.,10.])}
receipt={'status':'PASS_DISCRETE_GEOMETRY_ONLY' if all(checks.values()) else 'FAIL','checks':checks,'rows':rows,'orders':orders,'scope':'Exact rational extrema of stored P2 material-column heights; positivegeometry andPDEerror addressedseparately. Confirmsnumericalorderinginstability, notphysicaleffect.'}
(root/'evidence/activated_mesh_exact_maxima.json').write_text(json.dumps(serial(receipt),indent=2));print(json.dumps({'status':receipt['status'],'checks':checks},indent=2))
