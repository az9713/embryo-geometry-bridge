# Independent review of physical units and the force-scale quotient

**Result:** 41 independent checks pass for the corrected `physical_parameter_domain.py`, source SHA256

```text
8469fd4532a708247100cf8c19b2833e6af1566cc18d353d0ed99fe6d00724d4
```

The unit conversion, common-force-scale symmetry, exact quotient and noncompact-domain qualifications in `PHYSICAL-PARAMETER-DOMAIN.md` are correct in their stated scope. The review found a floating underflow issue; the author corrected it, and four independent extreme-input tests now verify the guards. The complete evidence, before/after hashes, field certificate and replay program are retained in `evidence/physical_parameter_domain_independent_review.json`.

This review performs no equilibrium optimization, biological data read, parameter fit or identified-set calculation. It does not repeat only the author's verification examples.

## 1. Independent dimensional conversion

The review uses a different synthetic physical material:

\[
 \mu=3.2\ {\rm Pa},\quad\kappa=960\ {\rm Pa},\quad
 (K_a,K_b,T_a,T_b)=(37.4,136,10.2,51)\ {\rm Pa\,\mu m},
 \quad g=17.6\ {\rm Pa/input}.
\]

The fixed unit is \(L_0=17\,\mu\mathrm m\), specimen height \(H=29\,\mu\mathrm m\), outer radius 68 µm, patch radius 7 µm, and active basal-height coordinates 4–24 µm. These values are test inputs, not measured biology.

Since \(dV=L_0^3d\hat V\), \(dA=L_0^2d\hat A\), and deformation gradients and area ratios are invariant to a common coordinate-unit conversion, dividing by \(\mu L_0^3\) produces

\[
 \hat\mu=1,\quad\hat\kappa=\kappa/\mu,
 \quad \hat K_s=K_s/(\mu L_0),\quad
 \hat T_s=T_s/(\mu L_0),\quad\hat g=g/\mu.
\]

The independently expected solver values are

\[
 (1,300,11/16,5/2,3/16,15/16,11/2),
\]

in the order shear, compression, apical/basal area stiffness, apical/basal tension, gain. The implemented conversion agrees. Its dimensionless height is 29/17. Because activation inputs are basal-height coordinates, multiplying each reported activation fraction by \(\hat H\) recovers 4/17 and 24/17 exactly up to floating rounding.

Changing height to 31 µm with the same fixed \(L_0\) changes geometry and activation fractions, while all material coefficients remain unchanged. If instead one uses \(H\) itself as the length unit, surface coefficients must change as \(1/H\). Holding \(K_s/(\mu H)\) fixed across changing height would prescribe a changing physical surface material. That distinction is correctly preserved by the helper.

The unit Pa·µm³ corresponds to \(10^{-18}\) J. Surface energy is consistent: \((\mathrm{Pa}\,\mu\mathrm m)(\mu\mathrm m)^2=\mathrm{Pa}\,\mu\mathrm m^3\). No absolute biological work or force is inferred by stating those units.

## 2. A second, independent length-unit covariance test

Keeping the same physical specimen and material, the review doubles the arbitrary coordinate unit from 17 to 34 µm. The corresponding reference and displacement coordinates halve; dimensionless surface coefficients halve. On the same physical test deformation:

\[
 \mathcal E_{34}=\mathcal E_{17}/8,
 \quad \mu(34)^3\mathcal E_{34}=\mu(17)^3\mathcal E_{17}.
\]

The implementation satisfies both raw dimensionless-energy scaling and reconstructed physical-energy invariance. The volume-normalized objective remains unchanged.

A useful numerical consequence is that the gradient with respect to dimensionless displacement is **coordinate-unit dependent**. If \(\hat u'=\hat u/2\), then \(\nabla_{\hat u'}\Phi=2\nabla_{\hat u}\Phi\). The independent test verifies that relation within 2.23e−16. Stationary zeros are invariant, but a finite force-residual tolerance must be interpreted with the declared fixed \(L_0\). This is not an error in the mapper: the actual study holds \(L_0=33.5\,\mu\mathrm m\) fixed. It is a reason to retain the length-unit convention in future solver receipts.

Changing coordinate units is separate from changing the physical activation protocol. For example, holding activation depths below the apical surface fixed while varying height requires the caller to update the supplied basal-height coordinates accordingly.

## 3. Independent verification on the actual high-quadrature kernel

A different random nonzero P2 field is used on the synthetic geometry above, with actual activation \(\hat g\times0.61\). Exact correlated Bernstein certificates prove the entire stored discrete field has \(J>0.97746\), so the fixed numerical inversion penalty is inactive everywhere, not merely absent at a few quadrature points.

The reviewer bypasses the normalization-to-unit-shear helper and multiplies every actual kernel coefficient, including activation gain, by factors **3.7** and **0.23**. These non-power-of-two factors differ from the author's examples. In the full pointwise physical potential,

\[
 E_{\alpha c}=\alpha E_c,
 \quad DE_{\alpha c}=\alpha DE_c,
 \quad D^2E_{\alpha c}=\alpha D^2E_c.
\]

The kernel's objective and derivatives normalized by \(\mu V_0\) must therefore remain invariant. Observed discrepancies are:

| Scale | Raw energy scaling error | Normalized objective error | Normalized gradient max error | Normalized Hessian-vector max error |
|---:|---:|---:|---:|---:|
| 3.7 | 5.69e−14 | 1.12e−16 | 1.12e−16 | 7.11e−15 |
| 0.23 | 3.56e−15 | 2.23e−16 | 7.03e−17 | 5.69e−14 |

The identity concerns admissible physical energy. The numerical extension contains a fixed inversion-penalty coefficient. If that term is active, scaling only the seven physical coefficients gives \(\alpha E_{phys}+P_{fixed}\), which generally differs from \(\alpha(E_{phys}+P_{fixed})\). The documented exclusion is essential and correct.

The symmetry also requires no separately fixed mechanical forcing. Independently known traction, pressure or body load can supply an absolute force reference. A fixed mobility in a dynamical model does not preserve relaxation times when energy is multiplied by \(\alpha\). The proposition correctly concerns quasistatic shapes and a consistently selected branch, not a unique equilibrium or universal dynamics.

## 4. Exact quotient and the noncompact admitted domain

The independent physical example gives the pressure-unit coefficient vector

\[
 c=(16/5,960,11/5,8,3/5,3,88/5),
\]

whose exact quotient is

\[
 w=\frac1{4973}(16,4800,11,40,3,15,88).
\]

Exact arithmetic verifies its sum, all inverse ratios \(w_j/w_\mu=c_j/c_\mu\), and invariance under three new rational scaling factors, including an extreme factor \(10^{80}\).

The normalized coordinate vector represents common-positive-scale equivalence classes. It does not make the admitted domain compact. For example,

\[
 c_n=(1,n,0,0,0,0,0),\qquad
 w_n=(1/(n+1),n/(n+1),0,\ldots,0)
\]

has positive admitted \(\mu,\kappa\) for every finite \(n\), but approaches the excluded \(w_\mu=0\) face and has unbounded ratio \(\kappa/\mu=n\). Exact tests include \(n=10^{100}\). The closed simplex is compact; the admitted subset with \(w_\mu,w_\kappa>0\) is not. Neither existence of minimizers, attained prediction extrema nor uniform solver accuracy follows from compact closure alone.

The exact quotient treats its supplied rational/decimal coefficient values exactly. It does not remove uncertainty in those values or certify a prior floating-point dimensional conversion. The floating convenience method is explicitly an approximation.

## 5. Issue found, author fix, and independent regression checks

The initially reviewed source, SHA256 `30291e23d0b82f08fc6e6e2bc42db541725f1c56db5728fd07bb6c259d281123`, used floating max-rescaling followed by normalization without checking loss of positive coordinates. With admissible finite coefficients \(\mu=5\times10^{-290}\), \(\kappa=7\times10^{290}\), the reproduced original floating formula returns \(w_\mu=0\), falsely landing on an excluded face.

The author corrected four paths: floating quotient, physical-pressure conversion, common coefficient scaling, and solver-ratio conversion. Each now rejects loss of any positive input coefficient to floating zero. The independent review checks four different extreme-input examples, including loss of a positive gain under scaling and loss of a gain under division by very large shear modulus. All raise the intended `ValueError`. The exact rational helper retains strictly positive \(\mu\) for the same extreme quotient.

This fix does not make arbitrary extreme ratios accurately representable in binary64. It prevents a silent domain change and directs such cases to exact arithmetic. No other correctness blocker was found in the reviewed scope.

## 6. Limits of the conclusion

The module correctly identifies a standard structural nonidentifiability and supplies consistent units. It does not compute the conditional compatible set \(\Theta_\delta\), its larger-region prediction envelope, discrepancy allowances, a likelihood or a confidence region. The documentary caveats about excluded simplex faces, branch selection, input errors and specimen-specific nuisance variables are appropriate. No novelty claim follows from this standard scale symmetry, its quotient or the numerical implementation checks.
