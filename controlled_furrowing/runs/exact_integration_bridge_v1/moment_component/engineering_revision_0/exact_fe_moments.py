"""Exact reference moments for the coherent finite-element integration bridge.

Rational Laurent moments on positive-radius triangles are rational + rational
linear combinations of logarithms of positive rationals. Axis cells must use
their already regularized Duffy chart. No Gaussian exactness is assumed.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from flint import arb, fmpq, ctx


def rational(value):
    if isinstance(value, bool):
        raise ValueError('Boolean is not a numerical input')
    return F.from_float(value) if isinstance(value, float) else F(value)


def A(value):
    q = rational(value)
    return arb(fmpq(q.numerator, q.denominator))


def powers(a, b=0):
    if any(not isinstance(x, int) or isinstance(x, bool) for x in (a, b)) or b < 0:
        raise ValueError('Integer radial and nonnegative integer axial exponents required')


@dataclass(frozen=True)
class LogMoment:
    constant: F = F(0)
    # Sorted tuples preserve a canonical, serializable symbolic expression.
    logarithms: tuple = ()

    def __add__(self, other):
        if not isinstance(other, LogMoment):
            other = LogMoment(rational(other))
        logs = dict(self.logarithms)
        for base, coefficient in other.logarithms:
            logs[base] = logs.get(base, F(0))+coefficient
        return LogMoment(self.constant+other.constant, tuple(sorted((x, c) for x, c in logs.items() if c)))

    def __mul__(self, scalar):
        c = rational(scalar)
        return LogMoment(self.constant*c, tuple((x, v*c) for x, v in self.logarithms if v*c))

    __rmul__ = __mul__

    def ball(self, precision_bits=256):
        old = ctx.prec
        try:
            ctx.prec = precision_bits
            result = A(self.constant)
            for base, coefficient in self.logarithms:
                if base <= 0:
                    raise ValueError('Positive logarithm argument required')
                result += A(coefficient)*A(base).log()
            return result
        finally:
            ctx.prec = old

    def as_dict(self):
        return {'rational': str(self.constant), 'log_terms': [
            {'positive_rational_argument': str(x), 'rational_coefficient': str(c)} for x, c in self.logarithms]}


def interval_power_moment(lower, upper, exponent):
    """Integral of R**exponent over a finite positive interval, exactly."""
    lo, hi = rational(lower), rational(upper)
    if not isinstance(exponent, int) or isinstance(exponent, bool):
        raise ValueError('Integer exponent required')
    if hi <= lo or lo < 0 or (lo == 0 and exponent < 0):
        raise ValueError('Nonintegrable/unordered power interval; regularize the axis first')
    if exponent == -1:
        return LogMoment(F(0), ((hi/lo, F(1)),))
    return LogMoment((hi**(exponent+1)-lo**(exponent+1))/F(exponent+1))


def canonical_triangle(vertices):
    if len(vertices) != 3 or any(len(p) != 2 for p in vertices):
        raise ValueError('Three meridian vertices required')
    p = tuple(tuple(rational(v) for v in x) for x in vertices)
    det = (p[1][0]-p[0][0])*(p[2][1]-p[0][1])-(p[2][0]-p[0][0])*(p[1][1]-p[0][1])
    if det == 0 or min(x[0] for x in p) < 0:
        raise ValueError('Nondegenerate triangle in R>=0 required')
    return p


@lru_cache(maxsize=250000)
def _triangle_moment(vertices, radial_exponent, axial_exponent, axial_shift):
    p = canonical_triangle(vertices)
    if radial_exponent < 0 and min(x[0] for x in p) == 0:
        raise ValueError('Negative powers on an axis triangle require an explicit regularized chart')
    cuts = sorted({x[0] for x in p})
    result = LogMoment()
    for lo, hi in zip(cuts, cuts[1:]):
        mid = (lo+hi)/2
        lines = []
        for v, w in zip(p, p[1:]+p[:1]):
            if min(v[0], w[0]) < mid < max(v[0], w[0]):
                slope = (w[1]-v[1])/(w[0]-v[0])
                intercept = v[1]-slope*v[0]-axial_shift
                lines.append((slope*mid+intercept, slope, intercept))
        if len(lines) != 2:
            raise ValueError('Two transverse triangle edges required per radial strip')
        lines.sort()
        degree = axial_exponent+1
        coefficients = [F(0)]*(degree+1)
        for sign, (_, slope, intercept) in zip((-1, 1), lines):
            for k in range(degree+1):
                coefficients[k] += sign*F(comb(degree, k), degree)*slope**k*intercept**(degree-k)
        for k, coefficient in enumerate(coefficients):
            if coefficient:
                result += interval_power_moment(lo, hi, radial_exponent+k)*coefficient
    return result


def triangle_moment(vertices, radial_exponent, axial_exponent, axial_shift=0):
    """Unweighted integral of R**a*(Z-shift)**b over an exact affine triangle."""
    if not isinstance(radial_exponent, int) or isinstance(radial_exponent, bool) or not isinstance(axial_exponent, int) or isinstance(axial_exponent, bool) or axial_exponent < 0:
        raise ValueError('Integer radial and nonnegative integer axial exponents required')
    return _triangle_moment(canonical_triangle(vertices), radial_exponent, axial_exponent, rational(axial_shift))


def normalized_bulk_moment(vertices, a, b, sigma, axial_shift=0, precision_bits=256):
    """(2*pi/sigma) integral_T R * R**a*(Z-shift)**b dR dZ."""
    powers(a, b)
    if rational(sigma) <= 0:
        raise ValueError('Positive fixed normalization required')
    old = ctx.prec
    try:
        ctx.prec = precision_bits
        return 2*arb.pi()/A(sigma)*triangle_moment(vertices, a+1, b, axial_shift).ball(precision_bits)
    finally:
        ctx.prec = old


def axis_chart_moment(kind, a, b, reference_radius, affine_determinant):
    """Integral of R*s**a*t**b times the exact material Duffy area Jacobian.

    R=s*R1. The chart Jacobian is s for one-axis-vertex cells and
    (1-s) for whole-axis-edge cells. Constitutive factors must be regularized.
    """
    if any(not isinstance(x, int) or isinstance(x, bool) or x < 0 for x in (a, b)):
        raise ValueError('Regularized chart polynomial exponents must be nonnegative integers')
    r, d = rational(reference_radius), abs(rational(affine_determinant))
    if r <= 0 or d == 0:
        raise ValueError('Positive chart geometry required')
    if kind == 'axis_vertex':
        return r*d/F((a+3)*(b+1))
    if kind == 'axis_edge':
        return r*d/F((a+2)*(a+3)*(b+1))
    raise ValueError('Explicit axis_vertex or axis_edge chart required')


def surface_moment(radius_lower, radius_upper, a):
    """Unnormalized integral R * R**a dR on a horizontal reference edge."""
    powers(a)
    return interval_power_moment(radius_lower, radius_upper, a+1)


def auxiliary_moment(points_weights, a, b=0, axial_shift=0):
    """Exact moment of ANY declared binary-rational positive auxiliary rule."""
    powers(a, b)
    total = F(0)
    for point, weight in points_weights:
        if len(point) != 2:
            raise ValueError('Two-coordinate moment point required')
        x, y, w = rational(point[0]), rational(point[1]), rational(weight)
        if w <= 0 or (a < 0 and x <= 0):
            raise ValueError('Positive rule weights and nonsingular sampled monomial required')
        total += w*x**a*(y-rational(axial_shift))**b
    return total


def bulk_moment_defect(vertices, points_weights, a, b, sigma, axial_shift=0, precision_bits=256):
    # Subtract the ACTUAL nonzero rule moment before absolute values. A nominal
    # Gaussian degree cannot replace this computation for binary cache values.
    old = ctx.prec
    try:
        ctx.prec = precision_bits
        return normalized_bulk_moment(vertices, a, b, sigma, axial_shift, precision_bits)-A(auxiliary_moment(points_weights, a, b, axial_shift))
    finally:
        ctx.prec = old


def harmonic(n):
    if not isinstance(n, int) or isinstance(n, bool) or n < 0:
        raise ValueError('Nonnegative integer harmonic index required')
    return sum((F(1, k) for k in range(1, n+1)), F(0))


def volumetric_derivative_coefficients(order, mu, kappa):
    """Power coefficients of phi'(1+x), phi''(1+x), through degree order."""
    if not isinstance(order, int) or isinstance(order, bool) or order < 0:
        raise ValueError('Nonnegative integer order required')
    mu, kappa = rational(mu), rational(kappa)
    first = [(-1)**(n+1)*(mu+kappa*harmonic(n)) for n in range(order+1)]
    second = [(-1)**n*(n+1)*(mu+kappa*harmonic(n+1)) for n in range(order+1)]
    return first, second


def volumetric_derivative_tail_bounds(order, rho, mu, kappa):
    """Rigorous scalar |phi'-P_N| and |phi''-Q_N| for |J-1|<=rho<1."""
    if not isinstance(order, int) or isinstance(order, bool) or order < 0:
        raise ValueError('Nonnegative integer order required')
    r, m, k = rational(rho), abs(rational(mu)), abs(rational(kappa))
    if not 0 <= r < 1:
        raise ValueError('0<=rho<1 required')
    gap = 1-r
    # H_(N+1+j) <= H_(N+1)+j/(N+2).
    first = r**(order+1)*((m+k*harmonic(order+1))/gap+k*r/F(order+2)/gap**2)
    # (n+1)*(mu+kappa*H_(n+1)), n=N+1+j.
    c, slope, n = m+k*harmonic(order+2), k/F(order+3), order+2
    second = r**(order+1)*(n*c/gap+(c+n*slope)*r/gap**2+slope*r*(1+r)/gap**3)
    return first, second
