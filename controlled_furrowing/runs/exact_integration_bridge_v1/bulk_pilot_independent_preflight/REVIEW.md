# Independent prelaunch review: four-element bulk pilot

**PASS21checks.** Reviewed source91a4f9a41b69c7b9e9ae5d5e0021087a7dfa4cf9e80625273bb4945d0ac3d0b3 and protocol_v2 hash746b260f5fbe7d049958681533165784e4847a15ca7d302ae9cdf5712cda95f0. This review runs only synthetic glue tests and frozen-input checks. It does not call the driver's `run_element` or perform actual field integration.

The proposed computation is record B at its exact saved center, four elements[0,6,72,66], constitutive series order16and512-bit Arb arithmetic. It excludes surfaces, assembled global defects, parameter/root-tube uncertainty and root transfer.

## Mathematical checks

The target gradient separates the polynomial shear/active part diag(1,1+drive,1,1+drive,1)q from phi′(J)J_u. Its Hessian is the corresponding constant-map quadratic term plus phi″(J)J_uJ_v+phi′(J)J_uv. The polynomial jet and moment-dual contraction include these terms, retaining signed moment defects before absolute values.

For an omitted series remainder bounded pointwise by B1/B2, the functional integral−Qstar has absolute remainder bounded by the **sum** of its positive integral mass and positive discrete mass. The gradient multiplier is sup|J_u|; the Hessian multiplier is B2sup|J_u|sup|J_v|+B1sup|J_uv|. The driver applies these weights once. Off-axis Laurent coefficient majorants use the positive radial interval and the containing shifted-axial rectangle. Axis chart exponents are regular polynomials on the unit square. Exact-rational radial/axial fixtures verify these majorants.

`Functional` uses the same zero/off-axis vertex ordering as the target charts. Its normalized material moment includes the affine determinant, Duffy factor and axisymmetric radial weight exactly once. Its finite rule uses the original normalized binary weights at explicitly defined coherent rational auxiliary nodes. The new512-bit power cache encloses those exact rational nodes and weights before integer powers and sums; it does not substitute regenerated Gaussian nodes or zero small defects.

The finite Qstar−cache correction is valid **at this point**: both constitutive evaluations use identical exact binary X+u, all fixed nodal coefficients remain included, and the same two-index local free projection is used in both maps. Direct outward-rounded subtraction at512bits is justified for this point calculation; it does not purport to solve the much harder full-tube dependency problem. A synthetic two-free-coordinate scalar energy is differentiated independently with Sympy. Its two gradient components and four Hessian entries overlap the driver's finite correction.

The existing canonical P2-on-affine-reference geometry box contains the exact saved center; coherent auxiliary nodes are strictly interior. Hence the reviewed rho<21/500 bound applies both under the material integral and at Qstar points. Old cache analyticity is checked by the local kernel and is a separate condition. The driver does not assume that cache J and canonical J are identical.

## Provenance and nonblocking limits

The preflight verifies every bound input in both protocol collections and detects conflicting normalized paths. The run script itself does not require the independent review receipt, so the parent must bind this exact source/protocol/receipt when launching. Actual output must remain unreviewed until the numerical enclosures are independently assessed.

One diagnostic field deserves careful interpretation: `unique_moment_defects` reads the classwide LRU cache size and is cumulative across Functional instances, rather than necessarily an element-local count. This affects a performance statistic, not a mathematical enclosure.

Two independent-review engineering stops are retained: a source-hash guard stopped during the producer's prefreeze performance revision, and an overstrict test initially demanded one rounded Arb enclosure contain another. The latter was corrected to test the exact expected Fraction against the production interval endpoints. Neither stop exposed a producer mathematical error. No production source, original proof or model file was edited by the reviewer.
