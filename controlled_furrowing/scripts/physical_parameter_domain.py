"""Explicit physical units and exact force-scale quotient for static mechanics.

This module does not fit biological data, impose a finite prior box, or claim
that an identified set has been computed. Lengths are micrometres, pressure
coefficients Pa, and surface coefficients Pa*micrometre. The solver is run
in one FIXED length unit; changing specimen height does not change material.
"""
from dataclasses import dataclass, asdict, replace
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import math
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
NAMES = ('mu', 'kappa', 'Ka_over_L0', 'Kb_over_L0',
         'Ta_over_L0', 'Tb_over_L0', 'gain')


@dataclass(frozen=True)
class PhysicalMaterial:
    mu: float
    kappa: float
    apical_area_modulus: float
    basal_area_modulus: float
    apical_tension: float
    basal_tension: float
    gain: float

    def __post_init__(self):
        values = tuple(asdict(self).values())
        if not all(math.isfinite(x) for x in values):
            raise ValueError('all physical coefficients must be finite')
        if self.mu <= 0 or self.kappa <= 0 or min(values[2:]) < 0:
            raise ValueError('positive shear/compression and nonnegative other coefficients required')

    def scaled(self, factor):
        if not math.isfinite(factor) or factor <= 0:
            raise ValueError('positive finite common force-scale factor required')
        scaled = {k: v*factor for k, v in asdict(self).items()}
        if any(v > 0 and scaled[k] == 0 for k, v in asdict(self).items()):
            raise ValueError('common coefficient scaling underflows; use exact arithmetic')
        return PhysicalMaterial(**scaled)

    def pressure_vector(self, length_unit_um=33.5):
        if not math.isfinite(length_unit_um) or length_unit_um <= 0:
            raise ValueError('positive finite fixed length unit required')
        result = np.array([self.mu, self.kappa,
                           self.apical_area_modulus/length_unit_um,
                           self.basal_area_modulus/length_unit_um,
                           self.apical_tension/length_unit_um,
                           self.basal_tension/length_unit_um, self.gain])
        if not np.isfinite(result).all():
            raise ValueError('coefficient conversion exceeds floating-point range')
        if any(value > 0 and converted == 0 for value, converted in zip(asdict(self).values(), result)):
            raise ValueError('coefficient conversion underflows; use exact arithmetic')
        return result

    def quotient(self, length_unit_um=33.5):
        """Floating approximation only; reject loss of any positive coefficient."""
        v = self.pressure_vector(length_unit_um)
        w = v/np.max(v)
        q = w/w.sum()
        if np.any((v > 0) & (q == 0)):
            raise ValueError('floating quotient underflows; use exact_quotient')
        return q


def exact_quotient(pressure_coefficients):
    """Rational scale quotient. String/decimal inputs retain their exact values."""
    c = tuple(Fraction(str(x)) for x in pressure_coefficients)
    if len(c) != 7 or c[0] <= 0 or c[1] <= 0 or min(c[2:]) < 0:
        raise ValueError('seven coefficients, positive mu/kappa and other entries nonnegative required')
    total = sum(c)
    return tuple(x/total for x in c)


def solver_parameters(material, *, height_um=33.5, outer_radius_um=100.5,
                      patch_radius_um=9.35512037,
                      activated_z_min_um=8.5, activated_z_max_um=28.5,
                      length_unit_um=33.5, nr=6, nz=3):
    """Return full-P2 parameters and gain in mu units, with fixed length L0.

    Activated z coordinates are heights above the basal reference surface.
    No material assumption is made by this geometry conversion.
    """
    from axisymmetric_quadratic import Parameters
    geom = (height_um, outer_radius_um, patch_radius_um,
            activated_z_min_um, activated_z_max_um, length_unit_um)
    if not all(math.isfinite(x) for x in geom):
        raise ValueError('finite geometry required')
    if not (height_um > 0 and outer_radius_um > patch_radius_um > 0 and
            0 < activated_z_min_um < activated_z_max_um < height_um and length_unit_um > 0):
        raise ValueError('positive geometry and strictly interior activation required')
    pressure = material.pressure_vector(length_unit_um)
    c = pressure/material.mu
    if not np.isfinite(c).all():
        raise ValueError('solver ratios exceed floating-point range')
    if np.any((pressure > 0) & (c == 0)):
        raise ValueError('solver coefficient ratios underflow; use exact arithmetic')
    p = Parameters(nr=nr, nz=nz, height=height_um/length_unit_um,
                   outer_radius=outer_radius_um/length_unit_um,
                   patch_radius=patch_radius_um/length_unit_um,
                   activation_z_min=activated_z_min_um/height_um,
                   activation_z_max=activated_z_max_um/height_um,
                   mu=1., bulk=float(c[1]), apical_area_modulus=float(c[2]),
                   basal_area_modulus=float(c[3]),
                   apical_surface_tension=float(c[4]), basal_surface_tension=float(c[5]),
                   quadrature=4, surface_quadrature=8)
    return p, float(c[6])


def source_hashes():
    names = ('physical_parameter_domain.py', 'quadratic_quadrature_study.py',
             'robust_quadratic.py', 'axisymmetric_quadratic.py',
             'axisymmetric_mixed.py', 'axisymmetric_prestress.py', 'axisymmetric_mechanics.py')
    return {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in names}


def verify():
    import torch
    from quadratic_quadrature_study import AxisymmetricModel
    checks = {}
    pressure = ('2', '10000', '.4', '4', '.2', '4', '12')
    q = exact_quotient(pressure)
    checks['exact_simplex_sum'] = sum(q) == 1
    checks['exact_positive_mu_kappa'] = q[0] > 0 and q[1] > 0
    checks['exact_scale_invariance'] = all(
        exact_quotient([str(Fraction(x)*scale) for x in pressure]) == q
        for scale in (Fraction(1, 101), Fraction(13, 2), Fraction(10**20)))
    checks['exact_ratio_inverse'] = all(q[i]/q[0] == Fraction(pressure[i])/Fraction(pressure[0]) for i in range(7))
    material = PhysicalMaterial(2., 10000., 13.4, 134., 6.7, 134., 12.)
    p, gain = solver_parameters(material)
    checks['physical_to_solver_expected'] = np.allclose(
        [p.mu, p.bulk, p.apical_area_modulus, p.basal_area_modulus,
         p.apical_surface_tension, p.basal_surface_tension, gain],
        [1., 5000., .2, 2., .1, 2., 6.], rtol=0, atol=2e-15)
    p_height, _ = solver_parameters(material, height_um=35., activated_z_max_um=30.)
    checks['height_changes_geometry'] = p_height.height != p.height
    checks['height_does_not_change_fixed_unit_surface_material'] = (
        p_height.apical_area_modulus == p.apical_area_modulus and
        p_height.basal_surface_tension == p.basal_surface_tension)
    checks['height_normalized_ratio_is_not_constant_material'] = (
        p_height.basal_surface_tension/p_height.height != p.basal_surface_tension/p.height)
    for scale in (.03125, 17., 1e6):
        pp, gg = solver_parameters(material.scaled(scale))
        checks[f'mapped_shape_scale_invariance_{scale}'] = all(
            math.isclose(asdict(p)[k], v, rel_tol=1e-14, abs_tol=1e-14)
            for k, v in asdict(pp).items()) and math.isclose(gain, gg, rel_tol=1e-14)
    rng = np.random.default_rng(42119)
    model = AxisymmetricModel(p)
    u = rng.normal(0, 1e-4, size=len(model.free))
    e, g = model.objective(u, gain*.37)
    parts = model.energy_components(model.coordinates(torch.tensor(u, dtype=torch.float64)), gain*.37)
    checks['physical_domain_trial_penalty_inactive'] = float(parts['trial_inversion_penalty']) == 0.
    rows = []
    # Deliberately bypass mu=1 mapping to test actual energy scaling in the
    # unchanged high-quadrature physical kernel, not only conversion arithmetic.
    for scale in (.125, 8.):
        ps = replace(p, mu=p.mu*scale, bulk=p.bulk*scale,
                     apical_area_modulus=p.apical_area_modulus*scale,
                     basal_area_modulus=p.basal_area_modulus*scale,
                     apical_surface_tension=p.apical_surface_tension*scale,
                     basal_surface_tension=p.basal_surface_tension*scale)
        other = AxisymmetricModel(ps)
        es, gs = other.objective(u, gain*.37*scale)
        pp = other.energy_components(other.coordinates(torch.tensor(u, dtype=torch.float64)), gain*.37*scale)
        raw_error = abs(float(pp['total'])-scale*float(parts['total']))
        grad_error = float(np.max(abs(gs-g)))
        checks[f'raw_energy_scales_{scale}'] = raw_error < 1e-10
        checks[f'normalized_objective_invariant_{scale}'] = abs(es-e) < 1e-12
        checks[f'normalized_gradient_invariant_{scale}'] = grad_error < 1e-11
        rows.append({'scale':scale, 'raw_energy_abs_error':raw_error,
                     'normalized_energy_abs_error':abs(es-e), 'gradient_linf_error':grad_error})
    bad = [lambda: material.scaled(0), lambda: material.scaled(float('inf')),
           lambda: material.pressure_vector(0), lambda: exact_quotient([1,0,0,0,0,0,0]),
           lambda: exact_quotient([1,1,-1,0,0,0,0]),
           lambda: solver_parameters(material, activated_z_max_um=40),
           lambda: PhysicalMaterial(1, 1, 0, 0, 0, -1, 0),
           lambda: PhysicalMaterial(1e-300, 1e300, 0, 0, 0, 0, 0).quotient(),
           lambda: PhysicalMaterial(1, 1, 1e-300, 0, 0, 0, 0).pressure_vector(1e300),
           lambda: PhysicalMaterial(1, 1, 1e-300, 0, 0, 0, 0).scaled(1e-100),
           lambda: solver_parameters(PhysicalMaterial(1e300, 1e300, 1e-300, 0, 0, 0, 0))]
    for i, fn in enumerate(bad):
        try:
            fn()
            checks[f'invalid_input_{i}'] = False
        except ValueError:
            checks[f'invalid_input_{i}'] = True
    receipt = {'status':'PASS' if all(checks.values()) else 'FAIL', 'checks':checks,
               'exact_example_quotient':dict(zip(NAMES,map(str,q))), 'kernel_scaling_rows':rows,
               'source_hashes':source_hashes(),
               'scope':'Force-scale quotient and physical-unit mapping only; no biological inference, uncertainty law, identified-set propagation or novelty claim.'}
    path = ROOT/'evidence/physical_parameter_domain_verification.json'
    path.write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'status':receipt['status'], 'checks':len(checks),
                      'failures':[k for k,v in checks.items() if not v], 'receipt':str(path)}), flush=True)
    if not all(checks.values()):
        raise RuntimeError('physical parameter mapping verification failed')


if __name__ == '__main__':
    verify()
