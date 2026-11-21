# Independent terminal analysis: safe uniform seven-element pilot

The producer completed the predeclared seven elements. This analysis passes946 read-only checks, reconstructing every stored derivative tail and the exact parameter/domain/source bindings. It does not repeat the nonlinear element integrals or the64-site cache sums. The prelaunch source proof remains a separate612-check component. Earlier unsafe integration packets remain withdrawn; these statements concern only SAFE V2.

![Local proof enclosures](local_bounds.png)

The yellow point is a conservative local maximum magnitude over the whole declared nodal/parameter box. The red cross is the largest component interval radius, not a statistical error or an attained physical variation. Blue is the independently retained corrected POINT certificate. Green bounds only the omitted constitutive-series remainder. Numeric values below are approximate; exact rational endpoints in results.json are authoritative.

|Element/chart|Uniform force norm upper|Force max radius|Point force norm upper|Uniform Hessian-entry norm upper|H max radius|Point H-entry norm upper|
|---|---:|---:|---:|---:|---:|---:|
|0/axis_vertex|1.911024e-15|5.518888e-22|1.911024e-15|1.416018e-13|9.412718e-19|1.416015e-13|
|2/axis_vertex|2.470538e-15|8.344498e-22|2.470537e-15|2.159149e-13|1.904240e-18|2.159144e-13|
|6/off_axis|2.136149e-13|9.495031e-18|2.136116e-13|3.388090e-09|1.130492e-13|3.388075e-09|
|8/off_axis|1.043893e-13|1.814847e-17|1.043711e-13|4.348787e-09|1.909907e-13|4.348761e-09|
|72/axis_edge|1.344228e-09|1.344228e-09|1.266405e-16|7.449034e-06|7.449034e-06|1.026381e-13|
|74/axis_edge|3.018187e-09|3.018186e-09|1.252350e-15|1.248463e-05|1.248463e-05|1.517503e-13|
|66/off_axis|4.943475e-13|1.819612e-18|4.943463e-13|4.201898e-12|7.638933e-15|4.196846e-12|

All seven force/Hessian packets overlap every corresponding point enclosure and have the identical local-to-global free map. Point overlap is a necessary consistency check, not a proof of uniformity. The previously reviewed interval algebra, same-Y cache correction, fresh whole-body geometry and two-mass tail bounds supply the uniform argument. All650 exact base coordinates lie inside the fresh common box, and the base gain lies in the same cell.

Elements72 and74 are the broad cases: their uniform intervals cannot resolve even the sign of any local force/H norm contribution. Their force radii exceed their point norms by about10.6million and2.41million, respectively; H-entry radii exceed the point norms by about72.6million and82.3million. The other five cases remain far more informative. The largest force remainder is below4.249e-23 and H-entry remainder below2.825e-19. Therefore the two broad boxes are not explained by series truncation. Higher series order cannot remove the displayed non-tail uncertainty.

The packets do not retain separate uniform polynomial and uniform cache contributions. Consequently this audit cannot attribute their width precisely to moment contraction, nodal coefficient dependence, or cache mean-value wrapping. The results are conservative enclosures; broadness does not demonstrate large true model discrepancy. Nothing here establishes whether the full assembled, preconditioned Hessian bound is small enough for contraction. Completing all144 bulk elements plus surfaces and the appropriate preconditioned norm is the next decision check. No611-dimensional, global root, PDE, stability, or biological fit certificate follows from these local maxima.

## A justified tightening candidate, without changing the model

The worst elements are axis-edge elements, not off-axis elements. A coefficient-only probe shows that forming the exact center and transformed basis first reduces their Jacobian coefficient-radius majorants only by factors1.35 and1.34. That is not evidence for a cure of their million-fold width growth. Physical radial shifting helps some off-axis coefficient majorants but does not address these worst cases.

There is a sharper exact structural identity. For an axis-edge P2 triangle, the radial field vanishes at all three nodes on the axis edge. The physical quadratic radial polynomial is therefore exactly divisible by the material radius. Its hoop strain is affine in physical radius and shifted height. All five strain fields can consequently be represented as physical polynomials of total degree at most1; their Jacobian has degree3/support10. The existing Duffy expansion has Jacobian degree6/support16. The new basis-only probe verifies exact coefficient equality of each basis and actual center field after transformation into the reviewed chart,48 exact checks. This supports a future physical-basis moment contraction that preserves these low-degree identities before interval expansion. It does not prove a tighter integral enclosure: no moment or integral was recomputed. Axis-vertex elements generally do not share this divisibility and must retain their regularized charts.

If the assembled bounds are inadequate, further principled options are to retain shared five-parameter correlations via a validated affine/Taylor representation, bound the centered defect variation directly by its derivative along the shared parameter path, or subdivide the exact same parameter domain with a complete cover. More floating precision addresses arithmetic rounding; it does not restore lost interval dependence. Any modified representation needs a new frozen implementation and review. No physical domain, prior, material parameter, or source tolerance is changed by these numerical proposals.

## Receipts

- results.json:946 checks, exact local statistics and independently reconstructed tails.
- coefficient_conditioning.json: same-domain coefficient-only experiment; no improved integral claim.
- axis_edge_basis_identity.json:48 exact structural identities.
- visual_generation.json and native_inspection.json: input/output hashes and observed readable PNG.
- terminal producer results SHA256:2f21379bf3e86e2d1bdb6663ca686c4585a709fbf8a2fba9e9be266f16aa7a0c.
