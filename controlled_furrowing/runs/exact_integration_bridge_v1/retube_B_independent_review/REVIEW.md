# Independent candidate-translation review

**PASS23checks.** Frozen preparation source e118f5f8c6e9e919847bfd596b439eaf99fbe945a224da8771a8bffed47e8131, protocol b71a90b45bae25a792f602d22915791601141dcf90000fd09afefa403a09f394, receipt5620374b95deaf5362c0e4818804226e217a9ffbb4d6522f3eebd3d83a974ef0. Independent receipt SHA25671a96ef4acd3525eccf252058fd245889b6449f3dea47b7e99bca2c282db4981. No field gradient/Hessian, integration, geometry or root calculation was performed.

## Exact arithmetic

The complete point discrepancy enclosure is D=R(FI−FQ)(u0,theta0). For each coordinate, the preparation requests a binary64 midpoint decrement, rounds the subtraction from the old binary64 center, and then defines the **actually applied** decrement d=u0−u1 as an exact rational. The reviewer checks all611midpoint conversions and611subtractions using exact predecessor/successor rounding cells and ties-to-even parity, rather than merely repeating the same floating subtraction. All650authoritative exact current coefficients and serialized-node rounding gaps are also checked.

The resulting quantities are approximately

    ||d||infinity =1.673692434644991e−11,
    sup_D ||D−d||infinity ≤3.452768216230464e−18,
    r_old+||d||infinity =2.673692434644991e−11.

Exact fractions are retained. The requested decrement and applied decrement are not assumed identical. The old radius remains a provisional candidate radius; the larger number is a **common comparison-domain radius**, not a certified root radius.

Exact fingerprints verify that C, the base parameters, R, S, source dose, material model, normalization and every reference/auxiliary cache array are unchanged. Only the predictor intercept becomes u1. All611old predictor hulls match. Independent translation of the old actual coordinate hulls reproduces396escaped and254contained coordinates out of650with respect to the old body-proof boxes. This disproves automatic box containment; it does not prove inversion or any physical failure. Fresh geometry/domain arguments remain necessary.

## Conditional residual lemma

Let p0(theta)be the old predictor and p1(theta)=p0(theta)−d. Define delta(u,theta)=R(FI−FQ)(u,theta), delta0=delta(p0(theta0),theta0), with delta0 in the available interval vector D. Assume the old theorem bounds ||RFQ(p0(theta),theta)|| by eta_old. On a freshly proved common domain, require

    ||I−R HQ|| ≤ZQ_fresh,
    ||R(HI−HQ)|| ≤b_fresh,
    ||delta(p0(theta),theta)−delta0|| ≤V.

The first two bounds must cover every translation segment p0(theta)−t*d, and eventually the candidate cubes if a contraction theorem is sought. Decompose

    RFI(p1(theta),theta)
      = RFQ(p0(theta),theta)
        +[RFQ(p1,theta)−RFQ(p0,theta)+d]
        +[delta(p1,theta)−delta(p0,theta)]
        +[delta(p0,theta)−delta0]+[delta0−d].

Integrating the respective Hessians along the translation segment bounds the second and third terms by ZQ_fresh*||d|| and b_fresh*||d||. The other terms give

    eta_new ≤eta_old+sup_D||D−d||
              +(ZQ_fresh+b_fresh)||d||+V.

The parameter/predictor term V includes direct parameter dependence as well as movement of p0(theta). It cannot be supplied by a point Hessian bound alone or omitted because the cell is small. The lemma gives only a conditional residual estimate; contraction additionally needs a fresh uniform derivative bound below1and a self-map inequality on a fully admissible candidate domain.

The preparation correctly leaves ZQ_fresh, b_fresh, V and eta_new as null, with domain/root assertions false. The available old ZQis explicitly labeled non-fresh. The common radius covers the translated cube by the triangle inequality; it supplies no new analytic or geometric certificate.

A separate four-check exact-arithmetic followup verifies the lemma's eta/shift payload against the independent earlier theorem and tests a nonlinear scalar example. It introduces no model evaluation. This is an audited preparation for the next proof, not that proof or a completed physical inference.
