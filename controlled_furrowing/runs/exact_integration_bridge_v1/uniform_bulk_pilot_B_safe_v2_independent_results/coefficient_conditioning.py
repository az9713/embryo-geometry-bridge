"""Coefficient-only conditioning probe, not an integral, no new proof packets."""
from pathlib import Path
from fractions import Fraction as Q
import sys,json
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent;sys.path.insert(0,str(ROOT/'scripts'))
import exact_uniform_bulk_defect_B_safe as U
T=U.T;J=U.JET;core=U.core
ctx.prec=512;core.A.cache_clear()
record=U.read(ROOT/U.P.RECORD);target=T.Target(U.read(ROOT/record['point_cache_file']));box=U.read(U.GEOM/'coordinate_boxes.json');Y=U.field_nodes(box)
def up(v):return max(abs(Q(str(v.lower().fmpq()))),abs(Q(str(v.upper().fmpq()))))if isinstance(v,arb)else abs(Q(v))
def rad(v):return (Q(str(v.upper().fmpq()))-Q(str(v.lower().fmpq())))/2 if isinstance(v,arb)else Q(0)
def norm(poly,f,fn):
 if f['chart']['kind']=='off_axis':
  lo=min(p[0]for p in f['corners']);hi=max(p[0]for p in f['corners']);z=max(abs(p[1]-f['zcenter'])for p in f['corners']);return sum((fn(c)*max(lo**i,hi**i)*z**j for(i,j),c in poly.items()),Q(0))
 return sum((fn(c)for c in poly.values()),Q(0))
rows=[]
for index in [0,2,6,8,72,74,66]:
 f=target.element(index,Y);j=J.bulk_jacobian_jet(f['q'],f['B']);ids=list(map(int,target.a['triangles'][index]));centers=[[Q(v)for v in box['current_centers_exact'][i]]for i in ids];radii=[[Q(v)for v in box['coordinate_radii_exact'][i]]for i in ids]
 # Form exact current-center field BEFORE interval uncertainty. Transform each
 # exact basis once, then add original independent nodal deviations once.
 allowed=[k for k in range(12)if k%2 or target.a['reference'][ids[k//2],0]>0]
 center=T.element_fields(f['corners'],centers,allowed);newq=[{k:U.A(v)for k,v in p.items()}for p in center['q']]
 for col,k in enumerate(allowed):
  r=radii[k//2][k%2];error=U.A(-r).union(U.A(r))
  for strain in range(5):newq[strain]=T.add(newq[strain],T.scale(center['B'][strain][col],error))
 newj=J.bulk_jacobian_jet(newq,f['B']);checks=[]
 # Point inclusion at center and two entire-node corners, componentwise exact.
 for sign in [-1,0,1]:
  v=[[c+sign*r for c,r in zip(cs,rs)]for cs,rs in zip(centers,radii)];pt=T.element_fields(f['corners'],v,f['free_local']);jp=J.bulk_jacobian_jet(pt['q'],pt['B'])
  checks.append(all(Q(str(newj['value'].get(k,arb(0)).lower().fmpq()))<=q<=Q(str(newj['value'].get(k,arb(0)).upper().fmpq()))for k,q in jp['value'].items()))
 olde=T.sub(j['value'],T.ONE);newe=T.sub(newj['value'],T.ONE)
 row={'element':index,'chart':f['chart']['kind'],'raw_J_coefficient_radius_majorant_exact':str(norm(olde,f,rad)),'basis_first_J_coefficient_radius_majorant_exact':str(norm(newe,f,rad)),'raw_J_minus_one_majorant_exact':str(norm(olde,f,up)),'basis_first_J_minus_one_majorant_exact':str(norm(newe,f,up)),'three_exact_nodal_witnesses_included':all(checks)}
 assert all(checks)
 rows.append(row)
U.P.write(OUT/'coefficient_conditioning.json',{'status':'PASS_COEFFICIENT_ONLY_PROBE','rows':rows,'analysis_source_sha256':U.sha(__file__),'no_integrals_cache_sums_or_new_uniform_packets':True,'scope':'Same independent nodal domain; algebraically equivalent exact-center plus exact basis times nodal errors. Does not establish improved final integral widths.'})
print(json.dumps([{'e':r['element'],'oldradius':float(Q(r['raw_J_coefficient_radius_majorant_exact'])),'newradius':float(Q(r['basis_first_J_coefficient_radius_majorant_exact'])),'ratio':float(Q(r['raw_J_coefficient_radius_majorant_exact'])/Q(r['basis_first_J_coefficient_radius_majorant_exact']))}for r in rows]))
