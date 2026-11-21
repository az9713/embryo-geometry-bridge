"""Cold and guarded scalar total-surface-stress screen at the certified K/L roots."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
import json
from math import comb
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXACT = ROOT/"runs/transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1"
KL = EXACT/"negative_KL_active_energy_v1"
KL_SOURCE = "072ed83a30565be34bd5fe7dde559526e089a070df0ad95f0e61611db9e9e76d"
KL_PROTOCOL = "06dd35434358ad7dcb993349281ed6bf8a84fd0a34f7ca2d4932ec0792842361"
WINDOWS = tuple((Q(i,16),Q(i+1,16)) for i in range(1,15))
REVIEW = "PASS_INDEPENDENT_COLD_TRANSVERSE_SURFACE_STRESS_SCREEN"


def require(ok, why):
    if not ok:
        raise ValueError(why)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def save(path, data):
    with path.open("x",encoding="utf8",newline="\n") as f:
        json.dump(data,f,indent=2)
        f.write("\n")
    return digest(path)


def load(path, expected):
    b=path.read_bytes()
    require(sha256(b).hexdigest()==expected,f"exact source or input digest {path}")
    return json.loads(b)


def const(c):
    n=Q(c)
    return n,n


def add(a,b):
    return a[0]+b[0],a[1]+b[1]


def times(c,a):
    x,y=c*a[0],c*a[1]
    return min(x,y),max(x,y)


def mul(a,b):
    x=[v*w for v in a for w in b]
    return min(x),max(x)


def polynomial_add(p,q):
    out=[const(0)]*max(len(p),len(q))
    for v in (p,q):
        for i,a in enumerate(v):
            out[i]=add(out[i],a)
    return out


def polynomial_mul(p,q):
    out=[const(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):
            out[i+j]=add(out[i+j],mul(a,b))
    return out


def polynomial_scale(p,c):
    return [times(c,a) for a in p]


def trace_coefs(a,m,b):
    return [a,add(add(times(-3,a),times(4,m)),times(-1,b)),add(add(times(2,a),times(-4,m)),times(2,b))]


def bernstein(p,lo,hi):
    n=len(p)-1
    power=[]
    for i in range(n+1):
        z=const(0)
        for k in range(i,n+1):
            z=add(z,times(comb(k,i)*lo**(k-i)*(hi-lo)**i,p[k]))
        power.append(z)
    result=[]
    for k in range(n+1):
        z=const(0)
        for i in range(k+1):
            z=add(z,times(Q(comb(k,i),comb(n,i)),power[i]))
        result.append(z)
    return min(z[0] for z in result),max(z[1] for z in result)


def edge(ref, y):
    R0,R1=ref[0],ref[2]
    length=R1-R0
    require(R0>=0 and length>0,"ordered positive reference edge")
    r,z=(trace_coefs(*vals) for vals in zip(*y,strict=True))
    dr,dz=([times(Q(i+1)/length,a) for i,a in enumerate(v[1:])] for v in (r,z))
    # Both Cartesian meridional tangent components and radius enter J_s^2.
    tangent=polynomial_add(polynomial_mul(dr,dr),polynomial_mul(dz,dz))
    numerator=polynomial_mul(polynomial_mul(r,r),tangent)
    R=[const(R0),const(length)]
    return R, r, dr, numerator


def screen_edge(ref,y,K,T,windows):
    require(K>=0 and T>=0,"declared nonnegative surface material")
    R,r,dr,numer=edge(ref,y)
    R2=polynomial_mul(R,R)
    alpha=1-T/K if K else Q(0)
    defect=None if alpha<=0 else polynomial_add(numer,polynomial_scale(R2,-alpha**2))
    rows=[]
    for lo,hi in windows:
        Rlo=R[0][0]+R[1][0]*lo
        Rhi=R[0][0]+R[1][0]*hi
        rlo,_=bernstein(r,lo,hi)
        slope,_=bernstein(dr,lo,hi)
        require(Rlo>0,"behavior strictly away from reference axis")
        require(rlo>0 and slope>0,"strictly positive radial trace and meridional slope")
        z=None if defect is None else bernstein(defect,lo,hi)
        negative=z is not None and z[1]<0
        nonnegative=alpha<=0 or z[0]>=0
        # D/R^2 <= D_hi / Rhi^2 < 0. Since 0<J_s<alpha,
        # alpha-J_s >= -D_hi/(2 alpha Rhi^2), without a sqrt.
        tau_upper=K*z[1]/(2*alpha*Rhi**2) if negative else None
        c=(-tau_upper)*slope**2/alpha if negative else None
        rows.append({"window_t_exact":[str(lo),str(hi)],"ref_radial_interval_exact":[str(Rlo),str(Rhi)],
                     "deformed_r_lower_exact":str(rlo),"deformed_r_R_lower_exact":str(slope),
                     "D_R2_times_Js2_minus_alpha2_exact":None if z is None else list(map(str,z)),
                     "tau_negative_upper_exact":None if tau_upper is None else str(tau_upper),
                     "minus_tau_r_R_squared_over_Js_lower_exact":None if c is None else str(c),
                     "standing":"STRICT_NEGATIVE_PATCH" if negative else "NONNEGATIVE_ON_WINDOW" if nonnegative else "UNRESOLVED"})
    return rows


def fixtures():
    ref=[Q(1),Q(3,2),Q(2)]
    windows=((Q(1,4),Q(3,4)),)
    def scaled(f,xi=Q(0)):
        return [[(u*f-xi,u*f+xi),const(0)] for u in ref]
    # r=R/2, r_R=1/2, J_s=1/4, K=4, T=1; tau=-2.
    n,=screen_edge(ref,scaled(Q(1,2)),Q(4),Q(1),windows)
    assert n["standing"]=="STRICT_NEGATIVE_PATCH"
    assert Q(-2)<=Q(n["tau_negative_upper_exact"])<0
    assert 0<Q(n["minus_tau_r_R_squared_over_Js_lower_exact"])<=2
    noisy,=screen_edge(ref,scaled(Q(1,2),Q(1,10000)),Q(4),Q(1),windows)
    assert noisy["standing"]=="STRICT_NEGATIVE_PATCH"
    at,=screen_edge(ref,scaled(Q(1)),Q(4),Q(0),windows)
    assert at["standing"]=="NONNEGATIVE_ON_WINDOW"
    positive,=screen_edge(ref,scaled(Q(1,2)),Q(4),Q(4),windows)
    assert positive["standing"]=="NONNEGATIVE_ON_WINDOW"
    # A line touching 0 in one Bernstein coefficient cannot certify strict negativity.
    assert bernstein([const(-1),const(1)],Q(0),Q(1))==(Q(-1),Q(0))
    return {"hand_exact_compressed_flat_trace_with_negative_surface_stress":True,
            "interval_perturbed_trace_still_negative":True,
            "exact_zero_boundary_and_large_tension_do_not_trigger_negative":True,
            "strict_bernstein_endpoint_guard":True}




def derivative_fixtures():
    # Inert symbolic constitutive checks, with no actual model arrays.
    import sympy as s
    t,a,b,c,d,h,u,v,w=s.symbols("t a b c d h u v w", real=True)
    F=s.Matrix([[a,0,b],[0,h,0],[c,0,d]])
    G=s.Matrix([[0,0,0],[0,0,0],[u,v,w]])
    assert s.expand((F+t*G).det())==s.expand(h*(a*d-b*c+t*(a*w-b*u)))
    assert s.expand((F+t*G)[:,2].dot((F+t*G)[:,2]))==b*b+d*d+2*d*t*w+t*t*w*w
    p=(F+t*G)[:,0]
    q=(F+t*G)[:,1]
    normal=p.cross(q)
    squared=s.expand(normal.dot(normal))
    assert squared==s.expand(h*h*(a*a+(c+t*u)**2)+a*a*t*t*v*v)
    v2=s.expand(s.diff(squared,v,2))
    assert v2==2*a*a*t*t
    return {"axial_rank_one_exact_det_and_active_ignore_azimuthal_gradient":True,
            "natural_disk_area_squared_has_only_second_order_azimuthal_term":True}


def inputs():
    require(digest(KL/"active_energy.py")==KL_SOURCE,"same validated K/L proof input function")
    spec=spec_from_file_location("_existing_exact_KL_proof_identity_only",KL/"active_energy.py")
    m=module_from_spec(spec)
    spec.loader.exec_module(m)
    retained=load(KL/"protocol.json",KL_PROTOCOL)
    cases,arrays,pins=m.inputs()
    require(retained==m.protocol(cases,arrays,pins),"exact same K/L material/root/posted coordinate authority")
    pins.update({str(p.relative_to(ROOT).as_posix()):d for p,d in
                ((KL/"active_energy.py",KL_SOURCE),(KL/"protocol.json",KL_PROTOCOL))})
    geo={
        m.K:(EXACT/"remaining32_geometry_v1/independent_actual_preparation/actual_replay"/(m.K+".json"),
             "a991e1e3431c09886361a1be0d6617f20d551423d5e110567fc2bbcc63a6d952","PASS_INDEPENDENT_REMAINING32_WHOLE_BOX_GEOMETRY_ONLY"),
        m.L:(EXACT/"independent_actual/replay"/(m.L+"_box_1.json"),
             "723e71ba60138282a934a76579bd50b02e11888edd7547d4a24cd2f4d8759bf9","PASS_FRESH_SINGLETON_WHOLE_GEOMETRY_ONLY"),
    }
    out=[]
    for a in cases:
        path,h,status=geo[a["case_id"]]
        g=load(path,h)
        require(g["status"]==status and g["pilot_id"]==a["case_id"] and g["box_sha256"]==a["box_sha256"]
                and Q(g["radius_exact"])==m.RADIUS and g["all_144_body_and_24_trace_domains"] is True
                and Q(g["strong_monotonicity_lower_exact"])>0
                and Q(g["body_J_bounds_exact"][0])>Q(1,10**6) and Q(g["min_trace_hoop_exact"])>0,
                "same exact whole-box analytic, geometry and monotonicity premises")
        pins[path.relative_to(ROOT).as_posix()]=h
        # Authenticate edge membership and orientation without evaluating surface stress.
        height=Q(a["material"]["height"])
        nodes=a["reference"]
        for side,z in ((0,Q(0)),(1,height)):
            edges=a["a"][f"surface_{side}_edges"].tolist()
            require(len(edges)==12 and all(x[-1]==y[0] for x,y in zip(edges[:-1],edges[1:],strict=True))
                    and nodes[edges[0][0]][0]==0 and nodes[edges[-1][-1]][0]==Q(a["material"]["outer_radius"])
                    and all(nodes[i][1]==z for e in edges for i in e)
                    and all(nodes[i][0]<nodes[k][0]<nodes[j][0] for i,k,j in edges),
                    "same fixed ordered flat boundary")
        out.append({**m.summarized(a),"geometry_proof_sha256":h,
                    "strong_monotonicity_lower_exact":g["strong_monotonicity_lower_exact"],
                    "body_J_bounds_exact":g["body_J_bounds_exact"]})
    return cases,out,pins


def protocol(cases,pins):
    return {"status":"FROZEN_COLD_NEGATIVE_KL_TRANSVERSE_SURFACE_STRESS_SCREEN",
            "source_sha256":digest(Path(__file__)), "derivation_sha256":digest(HERE/"READBACK.md"),
            "premise_sha256":pins,"two_exact_spatial_roots":cases,
            "physical_extension":"Same objective 3D bulk and material disk area energy with full wall clamp, no added bending or couplings.",
            "fixed_windows_on_each_ordered_edge":[list(map(str,w)) for w in WINDOWS],
            "cases":2,"faces_each":2,"edges_each":12,"total_windows":672,
            "criterion":"Strict negative tau on one annular subinterval at the exact root plus matched geometry implies energy descent in the natural 3D extension.",
            "actual_stress_evaluations":0,"new_3D_solve_or_AD_or_integration":0}


def prepare():
    checks={**fixtures(),**derivative_fixtures()}
    _,cases,pins=inputs()
    checks.update(same_exact_spatial_proof_boxes_material_and_surfaces=True,
                  analytic_geometry_and_global_admissibility_premises_matched=True,
                  did_not_evaluate_actual_surface_stress_or_modes=True)
    ph=save(HERE/"protocol.json",protocol(cases,pins))
    fh=save(HERE/"preflight.json",{"status":"PASS","source_sha256":digest(Path(__file__)),
           "protocol_sha256":ph,"checks":checks,"actual_stress_evaluations":0,
           "actual_AD_FE_jobs_or_spatial_integrals":0})
    print(json.dumps({"source_sha256":digest(Path(__file__)),"protocol_sha256":ph,"preflight_sha256":fh,
                      "actual_stress_evaluations":0,"current_bound_files":len(pins)}),flush=True)


def evaluate(case):
    free=set(map(int,case["a"]["free"]))
    e=case["posterior"]
    rows=[]
    for face,i in (("basal",0),("apical",1)):
        K=Q(case["effective_material"][face+"_area_modulus"])
        T=Q(case["effective_material"][face+"_surface_tension"])
        for no,indices in enumerate(case["a"][f"surface_{i}_edges"]):
            ids=list(map(int,indices))
            ref=[case["reference"][i][0] for i in ids]
            y=[]
            for index in ids:
                y.append([times(1,(case["center"][j]-e,case["center"][j]+e)) if j in free
                          else const(case["center"][j]) for j in (2*index,2*index+1)])
            for z in screen_edge(ref,y,K,T,WINDOWS):
                rows.append({"face":face,"edge":no,**z})
    return rows


def run(args):
    source=digest(Path(__file__))
    p=load(HERE/"protocol.json",args.protocol_sha256)
    f=load(HERE/"preflight.json",args.preflight_sha256)
    g=load(HERE/"independent_prelaunch_review.json",args.review_sha256)
    phase={"source_sha256":source,"protocol_sha256":args.protocol_sha256,"preflight_sha256":args.preflight_sha256}
    require(g["status"]==REVIEW and all(g[k]==v for k,v in phase.items())
            and type(g["checks"]) is dict and g["checks"] and all(x is True for x in g["checks"].values())
            and f["source_sha256"]==source and f["status"]=="PASS" and f["protocol_sha256"]==args.protocol_sha256
            and type(f["actual_stress_evaluations"]) is int and f["actual_stress_evaluations"]==0,
            "exact reviewed cold scalar phase")
    cases,auth,pins=inputs()
    require(p==protocol(auth,pins),"same cold proof authority and window set")
    out=HERE/"actual"
    out.mkdir(exist_ok=False)
    stamp={**phase,"review_sha256":args.review_sha256,"pid":os.getpid()}
    save(out/"started.json",{"status":"RUNNING",**stamp})
    term={"status":"FAILED","exit_code":1,**stamp}
    try:
        r=[{**a,"rows":evaluate(z)} for z,a in zip(cases,auth,strict=True)]
        require(sum(len(x["rows"]) for x in r)==672,"all preregistered windows")
        results=save(out/"results.json",{"status":"CANDIDATE_TRANSVERSE_SURFACE_STRESS_AWAITS_INDEPENDENT_READBACK",
            **stamp,"input_sha256":pins,"cases":r,"evaluated_windows":672,
            "negative_windows_per_case":{x["case_id"]:sum(v["standing"]=="STRICT_NEGATIVE_PATCH" for v in x["rows"]) for x in r},
            "continuum_equilibrium_or_physical_wavelength_claimed":False,"new_3D_solve_or_AD_or_spatial_integration":0})
        term.update(status="FINISHED",exit_code=0,results_sha256=results)
        print(json.dumps({"results_sha256":results}),flush=True)
    finally:
        save(out/"terminal.json",term)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    sp=p.add_subparsers(dest="mode",required=True)
    sp.add_parser("prepare")
    a=sp.add_parser("run")
    for name in ("protocol","preflight","review"):
        a.add_argument("--"+name+"-sha256",required=True)
    args=p.parse_args()
    prepare() if args.mode=="prepare" else run(args)

