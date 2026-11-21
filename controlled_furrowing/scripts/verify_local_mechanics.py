"""Independent geometric and symbolic checks; no biological inference."""
from pathlib import Path
import hashlib
import json
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def main():
    w, d, q, T, K = sp.symbols("w d q T K", positive=True)
    C = 3 * sp.sqrt(3) / 2
    V = C * q
    h = V / (C * (w*w + d*d/3))
    aa, ab = C*(w+d)**2, C*(w-d)**2
    al = 6*w*sp.sqrt(h*h+3*d*d)
    energy = T*(aa+ab) + K/(2*C)*(ab-C)**2 + sp.sqrt(3)*T/q*al
    hessian = sp.hessian(energy, (w,d)).subs({w:1,d:0}).applyfunc(sp.simplify)
    effective = sp.factor(hessian[0,0] - hessian[0,1]**2/hessian[1,1])
    det = sp.factor(hessian.det())
    ratio = sp.factor(-hessian[1,0]/hessian[1,1])
    expected_ratio = 3*K/(T+3*K+9*T/q**2)
    ratio_derivative = sp.factor(sp.diff(ratio,K))
    r = sp.symbols("r", positive=True)
    inverse = r*(1+9/q**2)/(3*(1-r))
    checks = {
        "unactivated_stationarity": all(sp.simplify(sp.diff(energy,v).subs({w:1,d:0})) == 0 for v in (w,d)),
        "asymmetry_ratio": sp.simplify(ratio-expected_ratio) == 0,
        "local_inverse": sp.simplify(expected_ratio.subs(K,T*inverse)-r) == 0,
        "height_zero_tilt_derivative": sp.diff(h,d).subs({w:1,d:0}) == 0,
        "height_width_derivative": sp.simplify(sp.diff(h,w).subs({w:1,d:0})+2*q) == 0,
    }
    # Independent mesh integration, rather than checking a formula against itself.
    rng = np.random.default_rng(20261119)
    geometry_errors = []
    for _ in range(100):
        width = rng.uniform(.3,3)
        delta = rng.uniform(-.8,.8)*width
        height = rng.uniform(.1,12)
        angles = np.arange(6)*np.pi/3
        xy = np.column_stack((np.cos(angles),np.sin(angles)))
        bottom = np.column_stack(((width-delta)*xy,np.zeros(6)))
        top = np.column_stack(((width+delta)*xy,np.full(6,height)))
        vertices = np.vstack([bottom,top])
        faces = []
        for j in range(1,5):
            faces.extend([(0,j+1,j),(6,6+j,6+j+1)])
        for j in range(6):
            k=(j+1)%6
            faces.extend([(j,k,6+k),(j,6+k,6+j)])
        areas=[]
        volume=0.
        for f in faces:
            a,b,c=vertices[list(f)]
            areas.append(np.linalg.norm(np.cross(b-a,c-a))/2)
            volume += np.dot(a,np.cross(b,c))/6
        lateral=sum(areas[8:])
        expected_v=float(C)*height*(width**2+delta**2/3)
        expected_l=6*width*np.sqrt(height**2+3*delta**2)
        geometry_errors.append(max(abs(volume/expected_v-1), abs(lateral/expected_l-1)))
    checks["independent_mesh_geometry_100_cases"] = max(geometry_errors) < 1e-12
    # Sylvester: Hww>0 and determinant has strictly positive polynomial terms.
    det_poly = sp.Poly(sp.expand(det*q**2),K,T,q)
    checks["positive_hessian_determinant_coefficients"] = all(c>0 for c in det_poly.coeffs())
    checks["positive_ratio_derivative"] = bool(ratio_derivative.is_positive)
    checks = {key:bool(value) for key,value in checks.items()}
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "scope": "Finite-height local cell geometry and local inference only; no tissue embedding or biological parameter identification",
        "checks": checks, "hessian": str(hessian), "determinant": str(det),
        "effective_width_stiffness": str(effective),
        "shape_to_width_response_ratio": str(ratio),
        "ratio_derivative_with_respect_to_K": str(ratio_derivative),
        "inverse_K_over_T": str(inverse),
        "maximum_relative_mesh_error": max(geometry_errors),
        "mesh_cases":100, "seed":20261119,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "novelty": "No claim: elementary mechanics and inverse-function calculation",
    }
    path=ROOT/"evidence/local_mechanics_independent.json"
    path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    plot()
    print(json.dumps(result))
    if not all(checks.values()):
        raise SystemExit(1)


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("dark_background")
    kb=np.linspace(0,10,300)
    fig,axes=plt.subplots(1,2,figsize=(11.8,4.8),constrained_layout=True)
    for aspect in [3,5,10]:
        ratio=3*kb/(1+3*kb+9/aspect**2)
        axes[0].plot(kb,ratio,label=f"H / cell edge = {aspect}")
    axes[0].set(xlabel="Basal modulus / surface tension",ylabel="Tilt response / width response",title="Same height slope can hide different shape responses")
    axes[0].legend(fontsize=8)
    r=np.linspace(0,.95,300)
    for aspect in [3,5,10]:
        k=r*(1+9/aspect**2)/(3*(1-r))
        axes[1].plot(r,k,label=f"H / cell edge = {aspect}")
    axes[1].set(xlabel="Measured tilt-to-width response ratio",ylabel="Inferred basal modulus / tension",title="A second observable resolves the local ambiguity")
    fig.suptitle("Synthetic finite-height local model: no embryo fit, no global curvature prediction",fontsize=12)
    path=ROOT/"docs/assets/local_identifiability.png"
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=170)
    plt.close(fig)


if __name__=="__main__":
    main()
