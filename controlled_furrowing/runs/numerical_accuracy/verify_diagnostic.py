"""Independent receipt replay plus symbolic checks; no new mechanical solves."""
from pathlib import Path
import json, hashlib
import sympy as sp
import numpy as np

out=Path(__file__).resolve().parent
root=out.parents[1]
summary=json.loads((out/'summary.json').read_text())
source=root/'scripts/axisymmetric_mechanics.py'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
checks={}
replays=[]
for s in summary['results']:
    p=out/s['artifact'];record=json.loads(p.read_text());a=record['result']
    ref=np.array(a['reference_nodes']);cur=np.array(a['current_nodes']);tri=np.array(a['triangles'])
    x=ref[tri];y=cur[tri]
    dm=np.stack((x[:,1]-x[:,0],x[:,2]-x[:,0]),axis=2)
    ds=np.stack((y[:,1]-y[:,0],y[:,2]-y[:,0]),axis=2)
    detmer=np.linalg.det(ds)/np.linalg.det(dm)
    # r=0 implies R=0: nodal ratio range over the non-axis vertices bounds
    # affine R/r everywhere inside each P1 triangle by convex weighting.
    extrema=[]
    for t,d in zip(tri,detmer):
        nz=ref[t,0]>0; ratios=cur[t[nz],0]/ref[t[nz],0]
        extrema.extend([d*ratios.min(),d*ratios.max()])
    minJ=float(min(extrema));maxJ=float(max(extrema))
    # Exact revolved volume of each straight-sided meridional triangle.
    volume=float(np.sum(np.linalg.det(ds)/2*(2*np.pi*y[:,:,0].mean(axis=1))))
    refvol=float(np.sum(np.linalg.det(dm)/2*(2*np.pi*x[:,:,0].mean(axis=1))))
    same_hash=sha(p)==s['sha256'] and record['source_sha256']==sha(source)
    stage_pass=all(v['stationarity_verified'] and v['scaled_gradient_linf']<1e-7 for v in a['continuation'])
    checks[f'bulk_{s["bulk"]}_hashes']=same_hash
    checks[f'bulk_{s["bulk"]}_stationary']=bool(a['converged'] and stage_pass and a['scaled_free_force_linf']<1e-7)
    checks[f'bulk_{s["bulk"]}_independent_J']=abs(minJ-a['minimum_J'])<1e-11 and abs(maxJ-a['maximum_J'])<1e-11 and minJ>0
    checks[f'bulk_{s["bulk"]}_independent_volume']=abs(volume-a['volume'])<1e-10 and abs(refvol-a['reference_volume'])<1e-10
    replays.append({'bulk':s['bulk'],'independent_min_J':minJ,'independent_max_J':maxJ,'independent_volume':volume,
                    'volume_absolute_difference':abs(volume-a['volume']),'file_sha256':sha(p)})
r,z,R,H=sp.symbols('r z R H',positive=True)
f=(1-r*r/(R*R))**2;g=z*z*(H-z)**2
ur=-r*f*sp.diff(g,z);uz=(2*f+r*sp.diff(f,r))*g
div=sp.simplify(sp.diff(r*ur,r)/r+sp.diff(uz,z))
checks['manufactured_axisymmetric_divergence_exact_zero']=div==0
checks['manufactured_rim_boundary_exact_zero']=sp.simplify(ur.subs(r,R))==0 and sp.simplify(uz.subs(r,R))==0
checks['manufactured_top_bottom_exact_zero']=all(sp.simplify(u.subs(z,b))==0 for u in (ur,uz) for b in (0,H))
J,k,q=sp.symbols('J k q',positive=True)
mixed=q*sp.log(J)-q*q/(2*k)
checks['mixed_elimination_exact']=sp.simplify(mixed.subs(q,k*sp.log(J))-k*sp.log(J)**2/2)==0
checks['mixed_pressure_residual_exact']=sp.simplify(sp.diff(mixed,q)-(sp.log(J)-q/k))==0
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows=summary['results'];x=[v['bulk'] for v in rows]
fig,ax=plt.subplots(1,3,figsize=(12,3.5),layout='constrained')
ax[0].plot(x,[v['minimum_J'] for v in rows],'o-',label='minimum J')
ax[0].plot(x,[v['maximum_J'] for v in rows],'s-',label='maximum J')
ax[0].axhline(1,color='0.5',ls=':');ax[0].set_ylabel('local volume ratio J');ax[0].legend(fontsize=8)
ax[1].plot(x,[100*v['relative_volume_change'] for v in rows],'o-')
ax[1].set_ylabel('global volume change (%)')
ax[2].plot(x,[v['maximum_activated_column_shortening']*33.5 for v in rows],'o-')
ax[2].set_ylabel('maximum column shortening (µm)')
for a in ax:a.set_xscale('log');a.set_xlabel('log-J penalty coefficient / µ');a.grid(alpha=.2)
fig.suptitle('Fixed activation f/µ=10, fixed 12×6 mesh: compressibility sensitivity\nDiscrete diagnostic; not a locking diagnosis or biological prediction',fontsize=11)
fig.savefig(out/'compressibility-diagnostic.png',dpi=160);plt.close(fig)
ev={'status':'PASS_EXECUTED_DIAGNOSTIC_AND_PLAN_ONLY' if all(checks.values()) else 'FAIL',
    'scope':'No mixed solver or complete locking diagnosis executed. Three fixed-load finite-compressibility diagnostics and independent geometry/hash replay.',
    'checks':checks,'results':summary['results'],'independent_replay':replays,
    'source_sha256':sha(source),'protocol_sha256':sha(out/'protocol.json'),
    'document_sha256':sha(root/'docs/NUMERICAL-ACCURACY-PATH.md'),
    'driver_sha256':sha(out/'run_diagnostic.py'),'verifier_sha256':sha(Path(__file__)),
    'visual':{'path':'runs/numerical_accuracy/compressibility-diagnostic.png','sha256':sha(out/'compressibility-diagnostic.png')},
    'recommended_next_method':'Axisymmetric P2 displacement/P1 pressure mixed saddle formulation preserving original logJ energy; verify weighted stability and nonlinear tangent.',
    'planned_not_executed':['Divergence-free manufactured forcing benchmark with P1 and mixed methods','Matched bulk5/50/500 by mesh8x4/12x6/16x8 at fixed f=.2 and separately10','Mixed finite-strain solver','Rigorous continuum-error certificate','Independent biological incompressibility measurement'],
    'literature':[{'title':'Ferrite nearly incompressible hyperelasticity tutorial','url':'https://ferrite-fem.github.io/Ferrite.jl/stable/gallery/quasi_incompressible_hyperelasticity/','used_for':'Primary software documentation: mixed Taylor-Hood finite-strain example; different volumetric law'},
                  {'title':'deal.II step-44','url':'https://dealii.org/9.7.0/doxygen/deal.II/step_44.html','used_for':'Primary software documentation: locking and three-field finite-strain alternative'},
                  {'title':'Simo Taylor Pister1985','doi':'10.1016/0045-7825(85)90033-7','used_for':'Original variational/projection volume-constraint formulation, cited by step-44; not independently line-by-line audited here'}],
    'limitations':['No continuum limit, global minimum, injectivity or biological adequacy certification','No causal attribution of discrete bulk sensitivity to locking','Global volume agreement does not exclude large local volume distortion','Zero-gap plug observer is not robust to arbitrary small surface errors']}
(root/'evidence/numerical_accuracy_plan.json').write_text(json.dumps(ev,indent=2)+'\n')
print(json.dumps({'status':ev['status'],'checks':checks,'source_sha256':ev['source_sha256'],'replay':replays},indent=2))
raise SystemExit(0 if all(checks.values()) else 1)
