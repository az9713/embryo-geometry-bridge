"""Exact axis-edge basis identity only: no integration."""
from pathlib import Path
from fractions import Fraction as Q
import sys,json
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent;sys.path.insert(0,str(ROOT/'scripts'))
import exact_uniform_bulk_defect_B_safe as U
T=U.T;J=U.JET;record=U.read(ROOT/U.P.RECORD);target=T.Target(U.read(ROOT/record['point_cache_file']));rows=[]
for index in [72,74]:
 f=target.element(index);ids=list(map(int,target.a['triangles'][index]));g=T.material_basis(f['corners']);phi=g['phi'];allowed=[k for k in range(12)if k%2 or target.a['reference'][ids[k//2],0]>0];physicalB=[[]for k in range(5)];checks={}
 for k in allowed:
  p=phi[k//2];r=k%2==0
  hoop=T.mul(p,{(-1,0):Q(1)})if r else{}
  checks['physical_hoop_no_negative_exponents_'+str(k)]=all(a>=0 for a,b in hoop)
  col=[T.diff(p,0)if r else{},T.diff(p,1)if r else{},{}if r else T.diff(p,0),{}if r else T.diff(p,1),hoop]
  for b,v in zip(physicalB,col):b.append(v)
 current=[target.nodes()[i]for i in ids];v=[current[k//2][k%2]for k in allowed];q=[T.linear_combination(row,v)for row in physicalB];ch=f['chart'];x=ch['R'];y=T.sub(ch['Z'],T.constant(f['zcenter']));checks['physical_field_transforms_exactly_to_reviewed_chart']=all(T.compose(p,x,y)==old for p,old in zip(q,f['q']))
 for k in range(len(allowed)):
  actual=T.element_fields(f['corners'],[[Q(0),Q(0)]for _ in range(6)],allowed)
  checks['each_basis_transforms_'+str(k)]=all(T.compose(physicalB[j][k],x,y)==actual['B'][j][k]for j in range(5))
 for a in range(5):checks['physical_q_degree_le1_'+str(a)]=all(i+j<=1 for i,j in q[a])
 assert all(checks.values())
 jet=J.bulk_jacobian_jet(q,physicalB)
 rows.append({'element':index,'checks':checks,'physical_q_degrees':[max((a+b for a,b in p),default=0)for p in q],'physical_J_degree':max((a+b for a,b in jet['value']),default=0),'physical_J_support':len(jet['value']),'physical_B_dof_count':len(allowed)})
U.P.write(OUT/'axis_edge_basis_identity.json',{'status':'PASS_EXACT_BASIS_IDENTITIES_ONLY','rows':rows,'source_sha256':U.sha(__file__),'scope':'Algebraic candidate for future tighter integration representation; no defect moments computed and no improved integral bound claimed.'})
print(json.dumps(rows))
