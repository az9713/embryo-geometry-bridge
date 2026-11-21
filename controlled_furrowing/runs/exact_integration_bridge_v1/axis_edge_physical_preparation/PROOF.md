# Same-target physical polynomial representation at a full axis edge

**Preparation only.** All six full-axis-edge elements, indices 72–77, passed 3,314 algebraic, enclosure and guard checks. No moment integral, cache job, improved error bound or root test has run.

Source: `../../../scripts/exact_axis_edge_physical.py`. API: `physical_element(target, index, nodes=None)` or `physical_fields(corners, current, free_local)`. The explicit chart kind is `axis_edge_physical`; callers must not silently send it through a positive-radius inverse-R norm routine.

## Exact equivalence

The material triangle is the same affine triangle with exact binary corner coordinates. Its coherent P2 basis is unchanged. Every current coefficient remains the original exact binary reference coefficient plus its displacement, or a supplied current nodal value in the existing common box. In particular, rounded reference midpoint offsets are not repaired.

On a full axis edge, the radial component at the two axis vertices and the axis midpoint is fixed to exactly zero. The restriction of the quadratic radial field to that edge is therefore identically zero. Hence

\[
r_h(R,\zeta)=R h_h(R,\zeta),\qquad h_h\in\mathbb P_1.
\]

The implementation divides each exact rational non-axis radial basis function by R **before** multiplying it by any interval coordinate. It requires every axis radial coefficient to be exact zero and every corresponding radial degree of freedom to be constrained. An interval that merely overlaps zero fails that requirement. No approximate coefficient cancellation or removal of zero-crossing coefficient balls supplies the divisibility proof.

Thus all five strain variables

\[
q=(r_R,r_\zeta,z_R,z_\zeta,r/R)
\]

and their selected free-coordinate basis columns are ordinary affine polynomials in physical material coordinates \((R,\zeta)\). Consequently

\[
J=(r_Rz_\zeta-r_\zeta z_R)(r/R)\in\mathbb P_3.
\]

The existing Duffy transformation is substituted into every physical polynomial using exact rational operations. For all six saved point fields, every q component and free basis column agrees coefficient by coefficient with the existing safe Duffy representation. The same identities were checked for J, its free gradient, and every upper-triangular free Hessian entry. The physical J degree is three and the Duffy J degree is six in every selected element. This is a coordinate/basis representation change, not a constitutive or geometric change.

## Uniform common-box enclosure

The same physical construction accepts interval current coefficients from the unchanged common nodal box and C. Safe arithmetic retains every nonexact zero-crossing coefficient. Each affine q coefficient is a linear functional of the full six current nodal coordinates, including fixed coordinates. Its exact rational extrema over the original local box were independently formed by signed endpoint sums and checked inside the resulting Arb coefficient enclosure.

Additionally, the exact center and 32 deterministic vertex samples per element were evaluated as rational point fields and checked against the interval coefficient enclosures. These samples support implementation verification; the uniform enclosure follows from the full linear interval construction, not sampling alone. No root membership in the explicit box is assumed.

## Integration contract for a future separately reviewed adapter

Physical integration uses the original triangle with measure

\[
\frac{2\pi}{\sigma}R\,dR\,d\zeta.
\]

All polynomial exponents are nonnegative, including those in J and its free derivatives. Ordinary physical-triangle moments are therefore integrable at R=0. A future adapter may use `normalized_bulk_moment` for these nonnegative exponents; that function already counts the radial measure R once. It must not add a second radial factor, retain a Duffy Jacobian, use negative-power off-axis moment assumptions, or call an inverse-R norm bound at zero.

This construction explicitly rejects elements with only one axis vertex. Such cells do not have the same polynomial divisibility guarantee and may retain direction-dependent hoop limits. The new chart does not reinterpret those cells.

## Retained evidence and limits

The protocol binds 156 files and the unchanged common-box/C identifiers. Guards cover nonfinite coefficients, nonexact zero-crossing constants, exact divisibility, invalid indices, free axis radial coordinates and axis-vertex rejection. Only representation tests are exposed by the command-line tool (`test` and `prepare`); no integration command is implemented.

Whether this representation materially tightens uniform force/Hessian defect bounds remains untested. Any six-element replacement requires a separately reviewed integration adapter and comparison against the full assembled bound, with all other elements and the physical target retained.
