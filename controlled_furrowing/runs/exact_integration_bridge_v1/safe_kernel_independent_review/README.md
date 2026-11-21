# Independent safe interval-polynomial kernel review

The immutable review receipt is `../../../evidence/exact_safe_polynomial_kernel_independent_review.json`. It is bound to the safe target SHA889b61cfa794d432021fbd4670b781169033d4b8451554b18e7f7a569d9d3954 and safe helper SHA1439e20ec127d39565467d0788173b1148b5ab47f659b3b39486b81aaf1ee373.

The review passed992 named checks. It performs no spatial integral, equilibrium solve, root proof or biological validation. It does not restore the withdrawn old integration claims.

## Mathematical policy

For a polynomial coefficient enclosure I, deleting a term is valid when I is EXACTLY the singleton zero. A nonzero-width Arb ball overlapping zero can contain nonzero coefficients and must remain. Arb `!=0` means provably nonzero, so it is not the logical complement of exact-zero detection.

The safe kernel uses finite-Arb `is_zero()` exclusively for deletion and for radial-node omission. The reviewer reproduced the unsafe negative control, tested zero-crossing constants and seven families of uncertain polynomial operations, and checked exact coefficient/value witness inclusion. Nonfinite coefficients are rejected by cleanup. Boolean powers are rejected.

Successful finite interval operations preserve coefficient/function enclosures on their nonsingular domains. This universal justification is interval arithmetic plus exact-zero-only deletion, not a claim established by vertex sampling alone.

Singular inverse powers or Laurent evaluations can return nonfinite Arb markers rather than raising. That is explicitly outside the successful finite/nonsingular domain: consuming code must verify finiteness and its analytic domain. The review tests that these singular cases do not produce finite certified values.

## Independent oracles

The review's exact polynomial oracle accepts rational coefficients ONLY and drops only exact rational zero. It does not reuse either target's sparse cleanup for its automatic differentiation.

- All144 actual B point-element dictionaries and all24 point-surface dictionaries agree EXACTLY with legacy rational algebra, coefficient for coefficient. Exact nodes, all auxiliary coordinates/weights, all17 parameter/cache records and the physical target definition agree.
- A separate rational polynomial second-order automatic-differentiation implementation checks every helper value/gradient/Hessian coefficient on four bulk elements and six surface edges. Uncertain-coefficient dual contractions are checked against exact direct signed-function evaluations for every local entry.
- Three synthetic chart classes exercise off-axis, axis-edge and axis-vertex fields. Additional axis tests keep non-axis radial coordinates whose uncertainty crosses zero, reject uncertainty on a constrained radial axis node, and verify canonical axis endpoint identity.
- Four actual B common-box elements0,2,6,8 are tested against a direct barycentric P2 gradient oracle and rational scalar second-order AD. Across synthetic and actual cases,10124 scalar strain/J/derivative witness values are included. The aggregate field tests retain85 nonzero-width zero-crossing coefficients.

Field samples supplement source/interval semantics; they do not establish root membership or whole-body physical admissibility. The coordinate boxes are numerical input domains, not experimentally measured uncertainty.

## Evidence

`review.py` is the retained independent verifier. `results.json` is the frozen local copy of its result. The evidence receipt binds both source families, safe preparation protocol/preflight, numerical inputs, test counts and scope. All input hashes were rechecked at completion. The kernel/helper remained unchanged during the audit.

The next step is separately reviewed v2 integral replay using this safe policy. Existing finite-quadrature roots and pure interval body-geometry results were not recalculated by this review; no exact-integral root or continuum/PDE claim follows from it.
