"""Conditional axisymmetric apical-slice observation, not an invented label fit.

The published classifier requires a central tissue plug separated from surrounding
tissue by perivitelline space. Missing acquisition slice height is retained as
observation uncertainty; a possible plug is not a positive label prediction.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np


def observe(top_r, top_z, bottom_r, bottom_z, uniform_surface_error=0.0, slice_height=None):
    tr,tz,br,bz = [np.asarray(a,dtype=float) for a in (top_r,top_z,bottom_r,bottom_z)]
    if any(a.ndim != 1 or len(a)<2 or not np.isfinite(a).all() for a in (tr,tz,br,bz)):
        raise ValueError("Finite one-dimensional surface arrays required")
    if len(tr)!=len(tz) or len(br)!=len(bz) or np.any(np.diff(tr)<=0) or np.any(np.diff(br)<=0):
        raise ValueError("Single-valued monotone-radius surfaces required")
    if tr[0]!=0 or br[0]!=0 or not np.isclose(tr[-1],br[-1],rtol=0,atol=1e-12):
        raise ValueError("Both surfaces must share axis and outer radius")
    if not np.isfinite(uniform_surface_error) or uniform_surface_error<0:
        raise ValueError("A finite nonnegative error bound is required")
    union=np.unique(np.r_[tr,br])
    if np.any(np.interp(union,tr,tz)<np.interp(union,br,bz)):
        raise ValueError("Crossed apical/basal surfaces")
    upper=float(min(tz[0],tz[-1]))
    lower=float(max(bz[0],bz[-1],np.min(tz)))
    gap=upper-lower
    error=2*float(uniform_surface_error)
    output={"scope":"Ideal continuous axisymmetric horizontal-slice observation; not specimen acquisition reconstruction",
            "possible_apical_void_slice_interval": [lower,upper] if gap>0 else None,
            "signed_interval_gap":gap,
            "conditional_gap_enclosure":[gap-error,gap+error],
            "uniform_surface_error_assumption":float(uniform_surface_error),
            "unknown_slice_label_set":[0] if gap+error<=0 else [0,1],
            "positive_existence_robust_under_error":gap-error>0,
            "point_label":None,
            "warning":"A nonempty existence interval does not predict an observed positive label without acquisition slice height and resolution. The supplied surface error is an assumption, not an internally proved FE error estimate."}
    if slice_height is not None:
        z=float(slice_height)
        if not np.isfinite(z):
            raise ValueError("Finite slice height required")
        output["point_label"] = int(lower<z<upper)
        output["point_label_scope"] = "Strict ideal observation at specified height for nominal piecewise-linear surfaces only"
    if gap>0:
        z=(lower+upper)/2
        output["witness_slice_height"]=z
        output["witness_empty_apical_annuli"]=empty_annuli(tr,tz,z)
    return output


def empty_annuli(r,z,slice_height):
    """Exact intersections of linear segments with the open set a(R)<z_slice."""
    segments=[]
    for x0,x1,y0,y1 in zip(r[:-1],r[1:],z[:-1],z[1:]):
        below0,below1=y0<slice_height,y1<slice_height
        if below0 and below1:
            segments.append([float(x0),float(x1)])
        elif below0 != below1:
            crossing=x0+(slice_height-y0)*(x1-x0)/(y1-y0)
            segments.append([float(x0 if below0 else crossing),float(crossing if below0 else x1)])
    merged=[]
    for left,right in segments:
        if merged and np.isclose(merged[-1][1],left,rtol=0,atol=1e-13):
            merged[-1][1]=right
        else:
            merged.append([left,right])
    return merged


def verify():
    flat=observe([0,1,2],[1,1,1],[0,1,2],[0,0,0])
    groove=observe([0,1,2],[1,.8,1],[0,1,2],[0,0,0])
    pit=observe([0,1,2],[.8,.9,1],[0,1,2],[0,0,0])
    high_base=observe([0,1,2],[1,.8,1],[0,1,2],[.9,0,.9])
    rng=np.random.default_rng(20261119)
    max_error_excess=0.
    for _ in range(1000):
        top=1+rng.uniform(-.2,.2,9)
        base=rng.uniform(-.2,.2,9)
        radial=np.arange(9,dtype=float)
        epsilon=.03
        changed_top=top+rng.uniform(-epsilon,epsilon,9)
        changed_base=base+rng.uniform(-epsilon,epsilon,9)
        before=observe(radial,top,radial,base)["signed_interval_gap"]
        after=observe(radial,changed_top,radial,changed_base)["signed_interval_gap"]
        max_error_excess=max(max_error_excess,abs(after-before)-2*epsilon)
    checks={"flat_no_existence":flat["possible_apical_void_slice_interval"] is None,
            "groove_interval_exact":np.allclose(groove["possible_apical_void_slice_interval"],[.8,1],atol=1e-15),
            "missing_slice_retains_ambiguity":groove["unknown_slice_label_set"]==[0,1],
            "pit_cannot_create_plug":pit["unknown_slice_label_set"]==[0],
            "basal_occupancy_retained":np.allclose(high_base["possible_apical_void_slice_interval"],[.9,1]),
            "witness_annulus":np.allclose(groove["witness_empty_apical_annuli"],[[.5,1.5]]),
            "specified_slice":observe([0,1,2],[1,.8,1],[0,1,2],[0,0,0],slice_height=.9)["point_label"]==1,
            "error_lipschitz_1000_perturbations":max_error_excess<=1e-14}
    result={"status":"PASS" if all(checks.values()) else "FAIL", "checks":checks,
            "scope":"Algorithm checks plus symbolic proof in SLICE-OBSERVER.md; finite-element error is NOT certified here",
            "random_seed":20261119,"max_positive_error_excess":max_error_excess,
            "source_sha256":__import__('hashlib').sha256(Path(__file__).read_bytes()).hexdigest()}
    root=Path(__file__).resolve().parents[1]
    (root/"evidence/spatial_observer_checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result))
    if not all(checks.values()):
        raise SystemExit(1)


if __name__=="__main__":
    verify()
