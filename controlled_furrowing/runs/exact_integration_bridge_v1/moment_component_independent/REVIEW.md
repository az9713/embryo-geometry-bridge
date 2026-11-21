# Independent moment-component review

**PASS on frozen source52bd1e4b54f4ae03ca66592bf4ab295fef30713e5fe171e6fc3e955acf2b9925.** The review is of a moment/series component, not an exact-integration bridge for any611-coordinate field. No equilibrium, Hessian assembly or field-force integration was run.

## Mathematical review

For a nondegenerate affine triangle in r≥0, sorting distinct radial vertex coordinates partitions the interior into strips with two affine axial boundaries. Integrating (z−shift)^b first, expanding those boundary powers, and integrating each remaining integer power of r gives exactly the module's rational-plus-log expression. Vertical edges are strip boundaries and have zero area. Permuting vertex orientation leaves this unsigned area integral unchanged. Negative radial powers on axis triangles are deliberately refused; the helper does not claim divergent Laurent terms become integrable individually.

Independent polynomial checks use another algorithm: expand r^a(z−shift)^b in barycentric monomials and apply the Dirichlet moment det(B)∏alpha_i!/(|alpha|+2)!. Thus the600cases do not merely repeat the module's radial slicing. They cover100integrands on four rational triangles, each in all six vertex permutations. Two separate analytic Laurent identities are checked: on vertices(1,0),(2,0),(1,1), the r^-1integral equals2log2−1 and the r^-2integral equals1−log2.

The axis-vertex chart has r=R1s and area measure |det(B)|s dsdt, giving weighted monomial moment R1|det(B)|/[(a+3)(b+1)]. The axis-edge chart instead has area factor1−s, giving R1|det(B)|/[(a+2)(a+3)(b+1)]. These formulas require the stated constantR1structured chart. A general unequal-radius vertex chart needs an additional rbar(t) factor and is outside this helper. Both negative affine orientation and positive geometry are covered in48exact independent measure tests.

For phi(J)=−mu logJ+kappa(logJ)^2/2, an independent formal convolution of log(1+x) and the geometric inverse reproduces both source coefficient arrays. Differentiating the first coefficient series independently gives the second. The first tail follows from H_(N+1+j)≤H_(N+1)+j/(N+2). The second follows by multiplying H_(N+2+j)≤H_(N+2)+j/(N+3) by N+2+j and summing1,j,j²against rho^j. Absolute mu/kappa permit the module's generic signed coefficients while maintaining an absolute majorant. The physical positive-coefficient case attains the sum-of-absolute-tail maximum at x=−rho. No differentiation of an error inequality is used.

The separately proposed surface theory is also checked: the magnitudes of coefficients of(1+x)^−p satisfy c_(n+1)/c_n=(n+p)/(n+1). For p=1/2the ratio is below1, so the remaining geometric ratio is bounded by rho. For p=3/2the ratios decrease; after the first omitted coefficient use rho(N+1+p)/(N+2)<1. These24closed-form tests verify the theory, not an exported surface-tail implementation: the reviewed module does not provide that function.

## Numerical accuracy and actual geometry evidence

The full independent harness passes **1,028checks**:

* 600exact barycentric polynomial/permutation cases;
* 15Laurent/log moments against independent rigorous Arb complex quadrature;
* 48exact axis-chart weighted measures;
* 18formal-convolution coefficient comparisons and144closed-form volumetric-tail checks;
* 24surface-theory tail checks;
* 144actual element weighted-constant moment defects, plus coverage checks;
* 20adversarial input checks, direct constructor normalization, known surface/normalized moments, precision restoration and frozen source checks.

The complex integrator's callback contains only integer-power rational meromorphic expressions. There is no unchecked log/sqrt branch callback; poles automatically return nonfinite balls to analytic probes. Three explicit triangular cross-sections are supplied independently, including one with two radial strips and shifted axial powers. All returned intervals overlap the symbolic moment enclosures and have widths below2^-150. The largest quadrature real-interval width in these15cases is approximately6.36×10^-61. The analytic derivation establishes the formulas; overlap is independent numerical corroboration, not a universal correctness proof from tests.

For every one of the144actual reference elements, an independent simplex moment computes the exact affine r integral, multiplies by rigorously enclosed2pi and the declared frozen normalization, then subtracts the exact rational sum of cached volume weights. These agree with the component and are **all certifiably nonzero**. Their absolute magnitudes range approximately3.60×10^-22to1.23×10^-17per element. This directly demonstrates why nominal Gaussian exactness cannot zero even the weighted constant defect. It does not bound force or Hessian defects: derivative factors, cancellation and accumulated preconditioning still have to be evaluated.

## Guard findings and limits

Before this freeze, review identified exponent validation bypasses in wrappers and an unchecked direct LogMoment constructor. The parent fixed them before this independent execution. The frozen module now rationalizes constructor data, merges duplicate bases, rejects nonpositive log bases and Boolean numerical data, and validates exponents before incrementing them. All20independent negative tests pass. The source's axis-chart documentation now states the constant-radius geometry restriction.

The initial883-case synthetic pass is retained under `synthetic_pass_revision_0`; the current version adds the actual144element comparison and full-coverage check. Original scientific packets were not changed. Remaining work includes interval Laurent coefficient construction over parameter/root tubes, explicit derivative-specific integration and cache-map consistency, preconditioned error assembly, successful root persistence, geometry/observer retransfers and independent replay. No all17integration bridge, PDE theorem, biological/source alignment or novel scientific result follows from this component review.
