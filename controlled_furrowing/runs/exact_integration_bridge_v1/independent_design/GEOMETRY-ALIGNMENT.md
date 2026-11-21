# Canonical exactly integrated field: existing geometry alignment

**Yes.** Existing all17whole-body J bounds refer to the exact affine P1corner reference and the exact P2current field y_h=Σφ_j(X_binary,j+u_j), with stored binary nodal X retained. This is a code/receipt compatibility check, not a new root or integration-error proof.

In `root_body_geometry.py`, `triangle_box` receives only the three corner reference coordinates, constructs the exact affine inverse, and forms P2derivatives from all six current nodal intervals. It multiplies the quadratic meridional determinant by the quadratic radial polynomial. The denominator in the resulting degree-four Bernstein comparison is the **affine reference radius**, degree-elevated as Σgamma_i r_i/4. The quadratic interpolant of stored reference midpoints is not used as that denominator.

`cohort_static_observers.nodal_box` constructs exact binary X+u, covers its gap from serialized saved current nodes, then adds the exact affine-predictor hull and full corrector radius. The F/K adapter `recentered_static_observers.translated_box` uses new cache u1, covers its gap from the original saved seed nodes, and retains the same predictor hull plus normal radius. The original-to-new center translation is included in the geometric hull, not falsely equated with the root's normal radius.

The independent read-only script `check_geometry_alignment.py` verifies the hash-linked cache/root/observer packets and **11,050coordinate inclusion inequalities across17records**, using exact rational arithmetic and no geometry/Hessian assembly. The existing outer boxes therefore contain every canonical current coefficient for every already-certified parameter/root family. Stored midpoint offsets and fixed binary Dirichlet values are retained; replacing them would require a separate bridge.

The largest all17bound is rho=max|J−1|≈0.0417090846537514, attained by the conservative L lower bound. Exact comparison proves rho<21/500=0.042. This applies for positive affine reference radius. The Bernstein numerator/denominator comparison proves polynomial inequalities lower*r≤det(F_meridional)*R≤upper*r, using correlated coefficients and exact zero-axis coefficients. It does not assume R/r is globally a polynomial or has a unique limit at an axis-only vertex. The6axis-edge and6axis-vertex-only triangles require the distinct cancellations/Duffy coordinates already described in DESIGN-REVIEW.

Consequently the scalar log series is uniformly legitimate. For x=J−1 and P_N(x)=Σ(k=1…N)(−1)^(k+1)x^k/k,

    |log(1+x)−P_N(x)| ≤ rho^(N+1)/[(N+1)(1−rho)],
    |d/dJ(logJ−P_N)| ≤ rho^N/(1−rho),
    |d²/dJ²(logJ−P_N)| ≤ N rho^(N−1)/(1−rho)+rho^N/(1−rho)^2.

The derivative bounds follow by differentiating the exact first-derivative remainder (−x)^N/(1+x), rather than differentiating an inequality. `geometry_alignment.json` retains exact rational bounds for N=8,10,12,16,20using outward rho=21/500. Turning them into FE gradient/Hessian or quadrature defects still requires chain-rule factors, log-squared terms, weighted integrals, axis regularization, parameter dependence and frozen cache-map consistency. None of those defects is certified by this narrow alignment check.
