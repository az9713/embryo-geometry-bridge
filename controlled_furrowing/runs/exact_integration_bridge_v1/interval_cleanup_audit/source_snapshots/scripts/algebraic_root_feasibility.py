"""Independent Arb attempt for ONE frozen binary-cache finite-quadrature root.

No Newton iteration/equilibrium optimization is performed.  The input vector is
an existing archived field.  Binary cache constants, X+u (without rounded nodal
addition), and a binary approximate inverse define the mathematical target.
All acceptance inequalities use Arb balls; floating comparisons are diagnostics.
"""
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
from dataclasses import asdict
import argparse, base64, datetime, hashlib, io, json, sys, time, traceback
import importlib.metadata
from types import SimpleNamespace
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu
from flint import arb, arb_mat, fmpq, ctx

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "evidence/algebraic_root_certification_feasibility.json"
STATE = ROOT / "runs/activated_quadrature/lambda_5000.0_strength_1.0/reequilibrated.json"
ARCHIVE = ROOT / "runs/material_tangent_inference/strength_1/tangents.npz"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(value):
    OUT.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


@lru_cache(maxsize=200000)
def A(value):
    """Enclose an exact binary input, never its shortest decimal approximation."""
    if isinstance(value, Fraction):
        q = value
    elif isinstance(value, (float, np.floating)):
        if not np.isfinite(value):
            raise ValueError("finite binary input required")
        q = Fraction.from_float(float(value))
    else:
        q = Fraction(value)
    return arb(fmpq(q.numerator, q.denominator))


def upper(value):
    return abs(value).abs_upper()


def pack(value):
    if not value.is_finite():
        raise ValueError("nonfinite ball")
    return {"ball": str(value), "lower_exact": str(value.lower().fmpq()),
            "upper_exact": str(value.upper().fmpq())}


def bulk_local(q, mu, kappa, drive):
    """Gradient/Hessian in (a,b,c,d,h), including active b,d quadratic terms."""
    a, b, c, d, h = q
    det = a*d-b*c
    if not det > 0 or not h > 0 or not det*h > A(Fraction(1, 1000000)):
        raise ValueError("ball leaves the positive-J analytic/no-penalty region")
    logj = det.log() + h.log()
    first = [d, -c, -b, a]
    v = [x/det for x in first] + [1/h]
    Q = [[arb(0) for _ in range(5)] for _ in range(5)]
    for i in range(4):
        for j in range(4):
            determinant_second = 1 if (i,j) in ((0,3),(3,0)) else -1 if (i,j) in ((1,2),(2,1)) else 0
            Q[i][j] = determinant_second/det-first[i]*first[j]/det**2
    Q[4][4] = -1/h**2
    factor = kappa*logj-mu
    diagonal = [mu,mu+drive,mu,mu+drive,mu]
    gradient = [diagonal[i]*q[i]+factor*v[i] for i in range(5)]
    H = [[kappa*v[i]*v[j]+factor*Q[i][j]+(diagonal[i] if i==j else 0)
          for j in range(5)] for i in range(5)]
    return gradient,H,det*h


def surface_local(q, K, T, scale):
    """Gradient/Hessian for K/2*(J_s-1)^2+T*J_s; J_s=scale*R*|t|."""
    R,x,z = q
    norm2 = x*x+z*z
    if not norm2 > 0 or not R > 0:
        raise ValueError("ball includes degenerate tangent or nonpositive radius")
    length = norm2.sqrt()
    js = scale*R*length
    first = [scale*length,scale*R*x/length,scale*R*z/length]
    second = [[arb(0),scale*x/length,scale*z/length],
              [scale*x/length,scale*R*z*z/length**3,-scale*R*x*z/length**3],
              [scale*z/length,-scale*R*x*z/length**3,scale*R*x*x/length**3]]
    factor = K*(js-1)+T
    return [factor*v for v in first],[[K*first[i]*first[j]+factor*second[i][j]
            for j in range(3)] for i in range(3)],js


def synthetic_checks():
    """Exact symbolic differentiation plus ball-domain and contraction checks."""
    import sympy as s
    a,b,c,d,h,mu,kappa,drive = s.symbols("a b c d h mu kappa drive", real=True)
    q=[a,b,c,d,h];D=a*d-b*c;L=s.log(D*h)
    energy=mu/2*(sum(t*t for t in q)-3)-mu*L+kappa/2*L**2+drive/2*(b*b+d*d)
    first=s.Matrix([d,-c,-b,a]);v=list(first/D)+[1/h]
    Q=s.zeros(5)
    for i in range(4):
        for j in range(4):Q[i,j]=s.diff(D,q[i],q[j])/D-first[i]*first[j]/D**2
    Q[4,4]=-1/h**2
    diag=[mu,mu+drive,mu,mu+drive,mu]
    checks={}
    for i in range(5):
        checks[f"bulk_gradient_symbolic_{i}"]=s.simplify(s.diff(energy,q[i])-(diag[i]*q[i]+(kappa*L-mu)*v[i]))==0
        for j in range(5):
            manual=diag[i] if i==j else 0
            manual+=kappa*v[i]*v[j]+(kappa*L-mu)*Q[i,j]
            checks[f"bulk_H_symbolic_{i}_{j}"]=s.simplify(s.diff(energy,q[i],q[j])-manual)==0
    R,x,z,K,T,scale=s.symbols("R x z K T scale",positive=True)
    u=[R,x,z];length=s.sqrt(x*x+z*z);js=scale*R*length
    ef=K/2*(js-1)**2+T*js
    point=[A(2),A(3),A(4)]
    g,H,_=surface_local(point,A(3),A(7),A(Fraction(1,2)))
    values={R:2,x:3,z:4,K:3,T:7,scale:s.Rational(1,2)}
    for i in range(3):
        checks[f"surface_gradient_rational_{i}"]=g[i].contains(A(Fraction(str(s.diff(ef,u[i]).subs(values)))))
        for j in range(3):
            checks[f"surface_H_rational_{i}_{j}"]=H[i][j].contains(A(Fraction(str(s.diff(ef,u[i],u[j]).subs(values)))))
    try:bulk_local([A(0),A(0),A(0),A(1),A(1)],A(1),A(2),A(3));checks["bad_J_rejected"]=False
    except ValueError:checks["bad_J_rejected"]=True
    try:surface_local([A(1),A(0),A(0)],A(1),A(1),A(1));checks["zero_tangent_rejected"]=False
    except ValueError:checks["zero_tangent_rejected"]=True
    r=A(Fraction(1,4));eta=A(Fraction(1,12));Z=A(Fraction(1,6))
    checks["scalar_sqrt2_Banach_ball"]=bool(Z<1 and eta+Z*r<r)
    checks["scalar_too_small_ball_not_accepted"]=not bool(eta+Z*A(Fraction(1,1000000))<A(Fraction(1,1000000)))
    checks["zero_preconditioner_not_accepted"]=not bool(A(1)<1)
    if not all(checks.values()):raise AssertionError(checks)
    return {"status":"PASS","checks":checks,"scope":"Exact local derivative identities and synthetic validation logic; not an actual root certificate"}


def freeze_cache(model,u,R):
    arrays={"reference":model.reference,"free":model.free,"u":u,"triangles":model.triangles,
            "shape":model.shape.numpy(),"shape_grad":model.shape_grad.numpy(),"qr":model.qr.numpy(),
            "volume_weights":model.integration.numpy(),"active":model.active.numpy(),
            "energy_scale":np.array([model.energy_scale]),"preconditioner":R}
    for k,surface in enumerate(model.surfaces):
        for name,item in zip(("edges","N","dN","r0","l0","weights"),surface[:6]):
            arrays[f"surface_{k}_{name}"]=item.numpy()
    buff=io.BytesIO();np.savez_compressed(buff,**arrays);data=buff.getvalue()
    return {"encoding":"base64 np.savez_compressed; each stored float is an exact binary constant",
            "sha256":hashlib.sha256(data).hexdigest(),"bytes":len(data),"data":base64.b64encode(data).decode(),
            "array_shapes":{k:list(np.asarray(v).shape)for k,v in arrays.items()}}


def local_terms(model,u,strength):
    """Yield exact-linear strain maps, including fixed nodes, at every Q point."""
    if not A(model.energy_scale)>0:
        raise ValueError("positive binary normalization required")
    index={int(k):i for i,k in enumerate(model.free)}
    coords=[A(float(x)) for x in model.reference.ravel()]
    for i,k in enumerate(model.free):coords[int(k)]+=A(float(u[i]))
    grads=model.shape_grad.numpy();N=model.shape.numpy();qr=model.qr.numpy();weights=model.integration.numpy();active=model.active.numpy()
    for t,nodes in enumerate(model.triangles):
        global_dofs=[2*int(node)+c for node in nodes for c in (0,1)]
        selected=[(j,index[k]) for j,k in enumerate(global_dofs)if k in index]
        for iq in range(len(N)):
            if not A(float(qr[t,iq]))>0:
                raise ValueError("positive reference quadrature radius required")
            B=[[arb(0)for _ in global_dofs]for _ in range(5)]
            for j in range(6):
                B[0][2*j]=A(float(grads[t,iq,j,0]));B[1][2*j]=A(float(grads[t,iq,j,1]))
                B[2][2*j+1]=A(float(grads[t,iq,j,0]));B[3][2*j+1]=A(float(grads[t,iq,j,1]))
                B[4][2*j]=A(float(N[iq,j]))/A(float(qr[t,iq]))
            q=[sum((v*coords[k]for v,k in zip(row,global_dofs)),arb(0))for row in B]
            Bs=[[row[j]for j,k in selected]for row in B]
            yield {"kind":"bulk","q":q,"B":Bs,"indices":[k for j,k in selected],
                   "weight":A(float(weights[t,iq]))/A(model.energy_scale),
                   "parameters":(A(model.parameters.mu),A(model.parameters.bulk),A(float(active[t,0 if active.shape[1]==1 else iq]))*A(strength))}
    for surface in model.surfaces:
        edges,Ns,dNs,r0,l0,ws,K,T=surface
        edges,Ns,dNs,r0,l0,ws=(x.numpy()for x in (edges,Ns,dNs,r0,l0,ws))
        for e,nodes in enumerate(edges):
            dofs=[2*int(node)+c for node in nodes for c in(0,1)];selected=[(j,index[k])for j,k in enumerate(dofs)if k in index]
            for iq in range(len(Ns)):
                if not A(float(r0[e,iq]))>0 or not A(float(l0[e]))>0:
                    raise ValueError("positive surface reference radius/length required")
                B=[[arb(0)for _ in dofs]for _ in range(3)]
                for j in range(3):
                    B[0][2*j]=A(float(Ns[iq,j]));B[1][2*j]=A(float(dNs[iq,j]));B[2][2*j+1]=A(float(dNs[iq,j]))
                q=[sum((v*coords[k]for v,k in zip(row,dofs)),arb(0))for row in B]
                Bs=[[row[j]for j,k in selected]for row in B]
                yield {"kind":"surface","q":q,"B":Bs,"indices":[k for j,k in selected],
                       "weight":A(float(ws[e,iq]))/A(model.energy_scale),
                       "parameters":(A(K),A(T),1/(A(float(r0[e,iq]))*A(float(l0[e]))))}


def assemble_point(terms,n):
    H=arb_mat(n,n);g=[arb(0)for _ in range(n)];cache=[];minJ=None
    for count,term in enumerate(terms):
        if not term["weight"] > 0:
            raise ValueError("strictly positive normalized quadrature weight required")
        fun=bulk_local if term["kind"]=="bulk"else surface_local
        grad,hessian,J=fun(term["q"],*term["parameters"])
        if term["kind"]=="bulk":minJ=J.lower()if minJ is None else min(minJ,J.lower())
        B=term["B"];k=len(B);d=len(term["indices"])
        if not d:continue
        BM=arb_mat(k,d,[v for row in B for v in row]);HM=arb_mat(k,k,[v for row in hessian for v in row]);GM=arb_mat(k,1,grad)
        localg=BM.transpose()*GM;localH=BM.transpose()*HM*BM;w=term["weight"]
        for a,i in enumerate(term["indices"]):
            g[i]+=w*localg[a,0]
            for b,j in enumerate(term["indices"]):H[i,j]+=w*localH[a,b]
        term["Hcenter"]=hessian
        term["absB"]=[[upper(v)for v in row]for row in B]
        term["rowsumB"]=[sum(row,arb(0)).upper()for row in term["absB"]]
        cache.append(term)
        if count%1536==0:print(json.dumps({"phase":"point_assembly","terms":count}),flush=True)
    return arb_mat(n,1,g),H,cache,minJ


def restore_cache(receipt):
    """Rebuild only frozen tensors, without generating quadrature or using a solver."""
    import torch
    packet=receipt["frozen_binary_cache"]
    raw=base64.b64decode(packet["data"])
    if hashlib.sha256(raw).hexdigest()!=packet["sha256"]:
        raise ValueError("frozen cache checksum mismatch")
    with np.load(io.BytesIO(raw))as z:
        arrays={key:z[key].copy()for key in z.files}
    parameters=SimpleNamespace(**receipt["parameters"])
    model=SimpleNamespace(reference=arrays["reference"],free=arrays["free"],triangles=arrays["triangles"],
        shape=torch.from_numpy(arrays["shape"]),shape_grad=torch.from_numpy(arrays["shape_grad"]),
        qr=torch.from_numpy(arrays["qr"]),integration=torch.from_numpy(arrays["volume_weights"]),
        active=torch.from_numpy(arrays["active"]),energy_scale=float(arrays["energy_scale"][0]),parameters=parameters)
    model.surfaces=[]
    for k,K,T in ((0,parameters.basal_area_modulus,parameters.basal_surface_tension),
                  (1,parameters.apical_area_modulus,parameters.apical_surface_tension)):
        model.surfaces.append(tuple(torch.from_numpy(arrays[f"surface_{k}_{name}"])
            for name in("edges","N","dN","r0","l0","weights"))+(K,T))
    return model,arrays["u"],arrays["preconditioner"]


def replay():
    """Read-only complete certificate replay from the self-contained cache packet."""
    receipt=json.loads(OUT.read_text());ctx.prec=receipt["precision_bits"];A.cache_clear()
    model,u,R=restore_cache(receipt);strength=float(Fraction(receipt["activation_strength_binary_exact"]))
    n=len(u);g,H,cache,minJ=assemble_point(local_terms(model,u,strength),n)
    RM=arb_mat(n,n,[A(float(v))for v in R.ravel()]);absRM=arb_mat(n,n,[A(abs(float(v)))for v in R.ravel()])
    defect=arb_mat(n,n)
    for i in range(n):defect[i,i]=1
    defect-=RM*H
    Z0=max(sum((upper(defect[i,j])for j in range(n)),arb(0)).upper()for i in range(n))
    correction=RM*g;eta=max(upper(correction[i,0])for i in range(n));results=[]
    for row in receipt["radius_attempts"]:
        r=A(Fraction(row["radius_requested_exact"]));vr,J,domain=variation_rows(cache,n,r);v=absRM*vr
        Z=Z0+max(upper(v[i,0])for i in range(n));good=bool(Z<1 and eta+Z*r<r)
        results.append({"radius_exact":row["radius_requested_exact"],"certified":good,"Z":pack(Z),"eta_plus_Zr":pack(eta+Z*r),"analytic_domain":domain})
    if not any(x["certified"]for x in results):raise RuntimeError("cached certificate replay did not prove a root")
    print(json.dumps({"status":"READ_ONLY_FROZEN_CACHE_REPLAY_PASS","eta":pack(eta),"Z0":pack(Z0),
        "terms":len(cache),"radius_results":results,"cache_sha256":receipt["frozen_binary_cache"]["sha256"]}),flush=True)


def variation_rows(cache,n,r):
    """Bound row sums of |H(u)-H(u*)| using low-dimensional strain balls."""
    rows=[arb(0)for _ in range(n)];minJ=None;margins={};counts={"bulk":0,"surface":0}
    for term in cache:
        q=[v+arb(0,(r*s).upper())for v,s in zip(term["q"],term["rowsumB"])]
        fun=bulk_local if term["kind"]=="bulk"else surface_local
        _,Hb,J=fun(q,*term["parameters"])
        counts[term["kind"]]+=1
        if term["kind"]=="bulk":
            minJ=J.lower()if minJ is None else min(minJ,J.lower())
            domain={"bulk_det_min":q[0]*q[3]-q[1]*q[2],"bulk_hoop_min":q[4]}
        else:domain={"surface_radius_min":q[0],"surface_tangent_squared_min":q[1]*q[1]+q[2]*q[2]}
        domain["normalized_weight_min"]=term["weight"]
        for name,value in domain.items():
            lo=value.lower();margins[name]=lo if name not in margins else min(margins[name],lo)
        k=len(q);change=[[upper(Hb[a][b]-term["Hcenter"][a][b])for b in range(k)]for a in range(k)]
        projected=[sum((change[a][b]*term["rowsumB"][b]for b in range(k)),arb(0))for a in range(k)]
        w=term["weight"]
        for j,i in enumerate(term["indices"]):rows[i]+=w*sum((term["absB"][a][j]*projected[a]for a in range(k)),arb(0))
    return arb_mat(n,1,[v.upper()for v in rows]),minJ,{"term_counts":counts,"strict_lower_margins":{k:pack(v)for k,v in margins.items()}}


def attempt(precision,radii):
    import torch
    from sparse_quadratic import AxisymmetricModel,Parameters
    torch.set_num_threads(1);ctx.prec=precision;A.cache_clear()
    state=json.loads(STATE.read_text());model=AxisymmetricModel(Parameters(**state["parameters"]));u=np.array(state["free_displacement"])
    with np.load(ARCHIVE)as z:
        assert np.array_equal(z["free_vector"],u)and np.array_equal(z["free_indices"],model.free)
        H0=csr_matrix((z["H_data"],z["H_indices"],z["H_indptr"]),shape=tuple(z["H_shape"]));R=splu(H0.tocsc()).solve(np.eye(len(u)))
    if not np.all(np.isfinite(R)):raise ValueError("nonfinite proposed preconditioner")
    n=len(u);started=time.monotonic()
    receipt={"status":"RUNNING_POINT_ASSEMBLY","created_utc":stamp(),"precision_bits":precision,
             "input_state":STATE.relative_to(ROOT).as_posix(),"input_state_sha256":sha(STATE),
             "matrix_archive":ARCHIVE.relative_to(ROOT).as_posix(),"matrix_archive_sha256":sha(ARCHIVE),
             "script_sha256":sha(__file__),"source_hashes":{k:sha(Path(__file__).with_name(k))for k in
              ("sparse_quadratic.py","quadratic_quadrature_study.py","robust_quadratic.py","axisymmetric_quadratic.py","axisymmetric_mixed.py","axisymmetric_prestress.py","axisymmetric_mechanics.py")},
             "backend":{"python_flint":importlib.metadata.version("python-flint"),"numpy":np.__version__},
             "target":"Exact real smooth potential with stored binary64 quadrature/shape/reference/constants; exact X+u. Not the rounded Torch evaluation and not ideal algebraic-Gauss nodes/PDE.",
             "no_new_equilibrium_optimization":True,"free_dof":n,"integration_policy":asdict(model.integration_policy),
             "parameters":state["parameters"],"activation_strength_binary_exact":str(Fraction.from_float(float(state["activation_strength"]))),
             "synthetic_checks":synthetic_checks(),"frozen_binary_cache":freeze_cache(model,u,R),"radius_attempts":[]}
    if OUT.exists():
        old=json.loads(OUT.read_text());history=list(old.get("attempt_history",[]))
        previous=old.get("previous_attempt_summary")
        keep=("status","created_utc","completed_utc","script_sha256","error","traceback","precision_bits",
              "preconditioned_residual_eta","center_preconditioned_H_defect_Z0","radius_attempts","elapsed_seconds")
        if previous:history.append({k:previous[k]for k in keep if k in previous})
        history.append({k:old[k]for k in keep if k in old});receipt["attempt_history"]=history
        receipt["initial_engineering_failure_note"]="First point-assembly run hit IndexError from active-mask shape144x1 broadcasting; fixed direct index to column0 for broadcastmask. That run produced no certificate. Larger-radius failures and later success are retained separately."
    save(receipt)
    try:
        g,H,cache,minJ=assemble_point(local_terms(model,u,state["activation_strength"]),n)
        receipt["status"]="RUNNING_PRECONDITIONED_BOUND";receipt["terms"]=len(cache);receipt["point_min_quadrature_J"]=pack(minJ);save(receipt)
        RM=arb_mat(n,n,[A(float(v))for v in R.ravel()]);absRM=arb_mat(n,n,[A(abs(float(v)))for v in R.ravel()])
        center=arb_mat(n,n)
        for i in range(n):center[i,i]=1
        center-=RM*H
        Z0=max(sum((upper(center[i,j])for j in range(n)),arb(0)).upper()for i in range(n))
        correction=RM*g;eta=max(upper(correction[i,0])for i in range(n));gmax=max(upper(g[i,0])for i in range(n))
        receipt["preconditioned_residual_eta"]=pack(eta);receipt["center_preconditioned_H_defect_Z0"]=pack(Z0);receipt["point_gradient_linf"]=pack(gmax)
        _,floatg=model.objective(u,state["activation_strength"]);receipt["floating_gradient_linf_diagnostic"]=float(np.max(abs(floatg)))
        receipt["max_gradient_Arb_mid_vs_Torch_diagnostic"]=float(max(abs(float(g[i,0].mid())-floatg[i])for i in range(n)))
        receipt["status"]="RUNNING_RADIUS_TESTS";save(receipt)
        for radius in radii:
            r=A(Fraction(radius));row={"radius_requested_exact":str(Fraction(radius))}
            try:
                vr,minbox,domain=variation_rows(cache,n,r);prop=absRM*vr;variation=max(upper(prop[i,0])for i in range(n));Z=Z0+variation
                lhs=eta+Z*r;good=bool(Z<1 and lhs<r)
                row.update({"variation_preconditioned_bound":pack(variation),"Z_total":pack(Z),"eta_plus_Zr":pack(lhs),"radius_ball":pack(r),
                            "analytic_positive_J_box_min":pack(minbox),"analytic_domain":domain,"strict_selfmap_and_contraction":good})
                if good:row["root_distance_upper"]=pack((eta/(1-Z)).upper())
            except ValueError as error:row.update({"strict_selfmap_and_contraction":False,"failure":str(error)})
            receipt["radius_attempts"].append(row);save(receipt);print(json.dumps({"phase":"radius","r":radius,"certified":row["strict_selfmap_and_contraction"]}),flush=True)
            if row["strict_selfmap_and_contraction"]:break
        success=any(v["strict_selfmap_and_contraction"]for v in receipt["radius_attempts"])
        receipt["status"]="CERTIFIED_LOCAL_ROOT_BINARY_QUADRATURE"if success else"UNRESOLVED_CONTRACTION_NOT_ESTABLISHED"
        receipt["completed_utc"]=stamp();receipt["elapsed_seconds"]=time.monotonic()-started
        receipt["scope_limits"]=["local unique algebraic root in declared ball only","no energetic stability/global uniqueness/branch theorem","no quadrature integration or continuum discretization bound","no source temporal operator or biological validation","no uniform statement across parameters or all17samples"]
        for name,h in receipt["source_hashes"].items():assert sha(Path(__file__).with_name(name))==h
        assert sha(STATE)==receipt["input_state_sha256"]and sha(ARCHIVE)==receipt["matrix_archive_sha256"]and sha(__file__)==receipt["script_sha256"]
        save(receipt);print(json.dumps({"status":receipt["status"],"seconds":receipt["elapsed_seconds"],"eta":str(eta),"Z0":str(Z0)}),flush=True)
    except BaseException as error:
        receipt["status"]="UNRESOLVED_EXECUTION_FAILURE";receipt["error"]=str(error);receipt["traceback"]=traceback.format_exc();receipt["completed_utc"]=stamp();save(receipt);raise


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("stage",choices=["test","attempt","replay"]);parser.add_argument("--precision",type=int,default=160);parser.add_argument("--radii",nargs="+",default=["1e-11","1e-10"]);args=parser.parse_args()
    ctx.prec=args.precision
    if args.stage=="test":
        value=synthetic_checks();print(json.dumps({"status":value["status"],"checks":len(value["checks"])}))
    elif args.stage=="attempt":attempt(args.precision,args.radii)
    else:replay()
