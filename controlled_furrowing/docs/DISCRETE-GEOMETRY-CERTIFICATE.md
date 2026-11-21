# Exact geometry verification for the saved finite-element bodies

Positive Jacobians at quadrature points alone do not exclude interpenetration. This check addresses the **actual saved discrete states**, independently of the optimizer. It is standard piecewise-linear topology and rational arithmetic, not a novel theorem or a continuum-error estimate.

## Exact local Jacobian bounds

Each meridional triangle has a piecewise-affine map \((r,z)\mapsto(R,Z)\). Let its positive reference/current oriented doubled areas be \(d_0,d_1\), giving constant meridional determinant \(d_1/d_0\). Revolving the triangle adds the hoop stretch \(R/r\), so

\[
J=\frac{d_1}{d_0}\frac Rr.
\]

With axis nodes mapped to the axis and all other radii positive,

\[
\frac Rr=
\frac{\sum_i\lambda_i r_i(R_i/r_i)}{\sum_i\lambda_i r_i}
\]

is a convex combination of the non-axis vertex ratios. Thus its extrema, including axis limits, lie between those ratios. Multiplying by the positive constant determinant gives exact element Jacobian bounds. This argument is specific to the saved P1 maps; it is not valid unmodified for quadratic elements.

## Global noninterpenetration of these discrete bodies

The checker verifies that the triangulation is an oriented disk: every interior edge has opposing incidences, every vertex link is a path or cycle as appropriate, and there is one boundary cycle. It checks strict positive orientation for every reference/current triangle and exact simplicity of both boundary polygons, including nonadjacent intersections and adjacent backtracking.

For a point off the image edges, the topological degree equals the winding number of the simple oriented boundary: one inside and zero outside. Every triangle preimage contributes +1 because its determinant is positive. Therefore image interiors cannot overlap or leave an uncovered interior hole. Continuity and the disk's vertex/edge neighborhoods extend injectivity to the edges and vertices. The finite piecewise-linear meridional map is consequently a homeomorphism onto its image. Preserving the axis and positive radius elsewhere makes its axisymmetric revolution injective as a three-dimensional body.

## Exact volume ratio

For a meridional triangle, integration of the affine radial coordinate gives

\[
V_t=\frac\pi3\,d_t\sum_{i=1}^{3}r_i.
\]

The factor \(\pi/3\) cancels from the current/reference volume ratio, making that ratio rational for the stored node coordinates. The checker converts every loaded binary64 coordinate into its exact rational value, checks all geometric predicates without rounding tolerances, and independently compares rational Jacobian/volume results with the solver diagnostics.

The certificate covers 51 directly calibrated and 93 larger-region states. It proves their discrete geometry is valid, even when their volume changes are too large to justify a biological incompressibility approximation. It does not prove a continuum solution is close to those states, certify optimizer stationarity, reconstruct an embryo, or identify material parameters. [Implementation](../scripts/certify_discrete_geometry.py); [exact receipt](../evidence/discrete_geometry_certificate.json).
