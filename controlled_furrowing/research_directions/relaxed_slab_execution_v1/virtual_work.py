"""Validated relaxed virtual work in a frozen admissible P2 basis direction."""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import argparse
import importlib.util
import json
from flint import arb, ctx

HERE=Path(__file__).resolve().parent
MECH=HERE.parents[1]


def digest(p):return sha256(p.read_bytes()).hexdigest()


def inputs():
    path=MECH/'runs/nonaxisymmetric_adequacy_audit_v1/surface_stress.py'
    assert digest(path)=='f33aaf277eea10fd364a529f95d369551c99b73747f7856cee3af4318b60b3a3'
    spec=importlib.util.spec_from_file_location('_unchanged_certified_surface_inputs',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    cases,_,pins=module.inputs()
    return cases,{**pins,str(path.relative_to(MECH)).replace('\\','/'):digest(path)}


def exact(q):return arb(str(q))


def interval(lo,hi):
    assert lo<=hi
    return exact(lo).union(exact(hi))


def positive_part(x):
    if x.upper()<=0:return arb(0)
    if x.lower()>=0:return x
    return arb(0).union(x.upper())


def coefficients(values):
    a,m,b=values
    return [a,-3*a+4*m-b,2*a-4*m+2*b]


def polynomial(c,t):
    out=arb(0)
    for a in reversed(c):out=out*t+a
    return out


def evaluate(case, subdivisions, component):
    indices=list(map(int,case['a']['surface_1_edges'][10]))
    node=indices[1];dof=2*node+component
    assert dof in set(map(int,case['a']['free']))
    assert sum(node in list(map(int,e)) for side in (0,1) for e in case['a'][f'surface_{side}_edges'])==1
    ref=[case['reference'][i][0] for i in indices]
    assert ref[0]>0 and ref[1]==(ref[0]+ref[2])/2
    length=ref[2]-ref[0];free=set(map(int,case['a']['free']));radius=case['posterior']
    fields=[]
    for dim in [0,1]:
        fields.append(coefficients([interval(case['center'][2*i+dim]-radius,case['center'][2*i+dim]+radius)
                                    if 2*i+dim in free else exact(case['center'][2*i+dim]) for i in indices]))
    r,z=fields
    dr=[r[1]/exact(length),2*r[2]/exact(length)]
    dz=[z[1]/exact(length),2*z[2]/exact(length)]
    K=exact(Q(case['effective_material']['apical_area_modulus']))
    T=exact(Q(case['effective_material']['apical_surface_tension']))
    total=arb(0);crossings=0
    for i in range(subdivisions):
        t=interval(Q(i,subdivisions),Q(i+1,subdivisions))
        R=exact(ref[0])+exact(length)*t
        rv=polynomial(r,t);rr=polynomial(dr,t);zr=polynomial(dz,t)
        ell=(rr*rr+zr*zr).sqrt()
        assert ell.lower()>0 and rv.lower()>0
        J=rv/R*ell
        compression=K-T-K*J
        if compression.lower()<0<compression.upper():crossings+=1
        n=4*t*(1-t);dn=(4-8*t)/exact(length)
        # R DJ = v_r ell + r(r_R v'_r + z_R v'_z)/ell.
        tangent=n*ell+rv*rr*dn/ell if component==0 else rv*zr*dn/ell
        total+=positive_part(compression)*tangent*exact(length/Q(subdivisions))
    lo=str(total.lower().fmpq());hi=str(total.upper().fmpq())
    return {'case':case['case_id'],'face':'apical','ordered_edge':10,'node':node,
            'dof':dof,'component':'radial' if component==0 else 'vertical',
            'direction':'unit P2 midpoint basis; zero at all other nodes',
            'reduced_virtual_work_exact':[lo,hi],
            'decimal_display_only':[float(total.lower()),float(total.upper())],
            'strict_nonzero':bool(total.lower()>0 or total.upper()<0),
            'subdivisions':subdivisions,'possible_compression_crossing_cells':crossings,
            'scope':'S/(2 pi) DE_rel[v] at the certified old finite-dimensional root; no continuum or biological equilibrium claim.'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--subdivisions',type=int,default=2048)
    p.add_argument('--component',type=int,choices=[0,1],default=0)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    assert args.subdivisions>=16 and not args.output.exists()
    ctx.prec=160
    cases,pins=inputs()
    rows=[evaluate(c,args.subdivisions,args.component) for c in cases]
    result={'status':'PASS_STRICT_NONZERO_RELAXED_VIRTUAL_WORK' if all(x['strict_nonzero'] for x in rows) else 'UNRESOLVED_SIGN',
            'source_sha256':digest(Path(__file__)),'input_sha256':pins,'precision_bits':160,
            'method':'Arb interval range quadrature with the full certified nodal box; positive-part law includes transitions.',
            'rows':rows,'new_equilibrium_solved':False,'novelty_claim':False}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps({'status':result['status'],'rows':[{k:r[k] for k in ['case','component','strict_nonzero','decimal_display_only','subdivisions']} for r in rows]}))


if __name__=='__main__':main()
