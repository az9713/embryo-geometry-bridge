# Root readback of the independent component replay preparation

This is a source/math readback, not an actual packet replay and not permission to
skip the pending producer/plan/process gates.

Inspected runner `replay.py` SHA256
`e504a67f913968888805c0385c6b3fa0501d6b3e7f5345362b453eb1dde50074`
and local arithmetic `replay_math.py` SHA256
`e7bb66066d54e131fd009d2bf48e94f3f1fb7bc1955023f66bc96cbe3c6d78e9`.
The preparation protocol/preflight are respectively
`2ec903e8088596b72d3f59e57ce6700247218b2138fbe9e78edeb1cb046cb9a1`
and `effd86b721146884ff4561e6bd24eabe3f32a260d768fd09218ef0432eb5442b`.

The reviewed implementation:

- Admits the root-observed supervisor and both real child exits, full ordered
  144+24 manifests, frozen source and both independent gates. It authenticates
  all planned and produced component bytes before parsing their payloads.
- Reopens the explicit singleton through the reviewed factory. That factory's
  authority calls bootstrap at 640 bits before any target reconstruction. It
  binds the same exact gain times dose, binary preconditioner, point anchor,
  entire radius-2^-34 box and independently reviewed cached derivative bound.
- Reconstructs canonical P2 and binary cache maps at the same full local nodal
  vector, with normalized weights applied exactly once. At the point it directly
  subtracts projected strain gradients, rather than relying on the producer's
  mean-value force correction. The nonlinear volume derivative uses the positive
  same-Y determinant path; the surface derivative uses the separate AD formula.
- Computes the point polynomial moment defect and its positive two-mass tail.
  It separately computes the entire-box Hessian moment defect and signed finite
  cache correction. The whole-box off-axis path uses the reviewed scaled radial
  chart; axis routes retain their physical regularization. Surface derivatives
  retain their separate squared-Jacobian/binomial expansion and error terms.
- Requires producer containment of the independently computed signed intervals,
  including the separately reported point cache correction. Mere overlap cannot
  pass. An enclosure that is wider because of arithmetic grouping remains HOLD;
  it is not rounded inward or replaced with an overlap test.
- Writes independent local intervals and their provenance before any eventual
  global reduction. Its successful receipt is a component-enclosure result, not
  a stationary point, minimum, continuum or biological-fit theorem.

The source shares reviewed monomial moments, polynomial jets, scaled-coordinate
transforms and unweighted third-derivative cache bounds. This is independent
composition and map/force reconstruction, not an independent re-proof of every
shared primitive. The current synthetic preflight is narrow. Actual route and
complete packet verification remain necessary, and may fail.

No blocking algebraic issue was identified in this readback. Actual replay must
still await a complete successful producer run, its root-observed true OS exit,
and the externally hashed admission described in README.md. A separate signed
reduction and independent review will be required after a passing replay.
