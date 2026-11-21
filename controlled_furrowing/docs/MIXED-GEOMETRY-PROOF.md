# Exact geometry certificates for the quadratic mixed solver

The verifier `scripts/mixed_geometry_certificate.py` uses exact rational arithmetic on supplied serialized node coordinates. Its claims concern the **discrete P2 map**, not the unknown continuum solution. It does not establish equilibrium, physical material adequacy, incompressibility, experimental prediction, or novelty.

## 1. Reference and current elements

Let \(x_0,x_1,x_2\in\mathbb Q^2\) define a positively oriented affine reference meridional triangle. On \(\widehat T=\{(\xi,\eta):\xi,\eta\ge0,\xi+\eta\le1\}\), set

\[
 x=x_0+(x_1-x_0)\xi+(x_2-x_0)\eta,
 \quad D_m=[x_1-x_0,x_2-x_0],\quad \det D_m>0.
\]

Current coordinates \(y=(R,Z)\) are quadratic Lagrange interpolants of six nodes ordered
\([v_0,v_1,v_2,m_{01},m_{12},m_{20}]\). Each input Python float is converted with `Fraction.from_float`, which represents that binary64 value **exactly**. Fractions and rational strings are also accepted. This contains no uncertainty allowance for the mathematical solution before serialization.

For a scalar component with six nodal values \((a,b,c,d,e,f)\), its power coefficients are

\[
 y=a+(-3a-b+4d)\xi+(-3a-c+4f)\eta
 +(2a+2b-4d)\xi^2+(4a-4d+4e-4f)\xi\eta
 +(2a+2c-4f)\eta^2.
\]

Thus \(D_{\xi,\eta}y\) is affine and

\[
 d(\xi,\eta)=\det F_m=
 \frac{\det D_{\xi,\eta}y}{\det D_m}
\]

is quadratic. The code forms this polynomial by exact differentiation, multiplication and subtraction; it does not infer the determinant from quadrature samples.

## 2. Whole-element meridional determinant bounds

Use the degree-two triangular Bernstein basis

\[
 (B_0,B_1,B_2,B_{01},B_{12},B_{20})
 =(\lambda_0^2,\lambda_1^2,\lambda_2^2,
   2\lambda_0\lambda_1,2\lambda_1\lambda_2,2\lambda_2\lambda_0).
\]

Every basis function is nonnegative and their sum is one. Consequently, for any quadratic \(p=\sum b_i B_i\),

\[
 \min_i b_i\le p\le\max_i b_i \quad\hbox{on the entire closed triangle}.
\]

The vertex coefficients equal the vertex values. The edge coefficient corresponding to endpoints \(i,j\) is
\(2p(m_{ij})-[p(v_i)+p(v_j)]/2\). These identities determine all six coefficients exactly. A strictly positive minimum determinant coefficient is therefore a sufficient certificate of \(\det F_m>0\) everywhere. A nonpositive coefficient is **not** proof of an inverted element: the coefficient bound may be conservative.

For example, on a translated unit reference triangle \(r=1+\xi,z=\eta\), take

\[
 R=1+\xi+0.9\eta^2,\qquad Z=\eta+0.9\xi^2.
\]

Its determinant is \(1-3.24\xi\eta\ge1-0.81=0.19>0\), while one degree-two Bernstein coefficient is \(1-1.62=-0.62\). The verifier correctly reports `UNCERTIFIED`, not physical inversion. Adaptive exact subdivision could tighten this result; it has not been implemented in this version.

This convex-hull method is established, including its use for curved finite-element Jacobians; see Johnen, Remacle and Geuzaine, [Geometrical validity of curvilinear finite elements](https://www.sciencedirect.com/science/article/pii/S0021999112005219), *Journal of Computational Physics* 233 (2013). Our work supplies an exact arithmetic implementation for the present discrete model, rather than a new validity principle.

## 3. Hoop stretch, including the axis

The axisymmetric three-dimensional Jacobian is

\[
 J=d\,h,\qquad h=R/r.
\]

Write both numerator and denominator in the same degree-two Bernstein basis. Since reference \(r\) is affine, its elevated coefficients are

\[
 (r_0,r_1,r_2,(r_0+r_1)/2,(r_1+r_2)/2,(r_2+r_0)/2).
\]

They are nonnegative for a triangle in \(r\ge0\). Let these be \(c_i\), and let \(a_i\) be the exact current-radius Bernstein coefficients. Require \(a_i=0\) whenever \(c_i=0\). If

\[
 0<\ell_h=\min_{c_i>0}\frac{a_i}{c_i},\qquad
 u_h=\max_{c_i>0}\frac{a_i}{c_i},
\]

then \(\ell_h c_i\le a_i\le u_h c_i\) for every index. Multiplying by \(B_i\ge0\) and summing proves

\[
 \ell_h r\le R\le u_h r.
\]

For every \(r>0\), division yields \(\ell_h\le R/r\le u_h\). The inequalities extend to axis limits. On a triangle with an entire edge on the axis, the polynomial \(R\) vanishes identically along that edge, so \(R\) is divisible by its linear reference coordinate \(r\), and the axis limit exists as an affine polynomial. For a triangle touching the axis at only one vertex, the ratio bound remains valid along every approach with \(r>0\); no direction-independent pointwise axis derivative is asserted there.

With certified \(0<\ell_d\le d\le u_d\), multiplication gives the rigorous sufficient enclosure

\[
 0<\ell_d\ell_h\le J\le u_du_h.
\]

This is a whole-element enclosure for the stored discrete map. It does not equate the separate extrema of \(d\) and \(h\), and it does not claim those product endpoints are attained.

## 4. Conservative boundary and global embedding certificate

Positive local determinants alone do not rule out distant self-overlap. `certify_mesh` therefore separately verifies:

1. The reference vertex triangulation is a topological disk: consistent oriented shared edges, Euler characteristic one, one boundary cycle and path/cycle vertex links. Its P1 reference boundary is simple.
2. Each interior P2 edge uses the same global midpoint index from both adjacent elements; distinct edges use distinct midpoint indices. Together with shared endpoint nodes, this makes the P2 trace continuous.
3. The reference boundary is a flat axis rectangle. The current axis remains \(R=0\), and the rim has one fixed positive radius.
4. Current top and bottom quadratic edges are strictly increasing in \(R\); axis and rim edges are strictly increasing in \(Z\), using reference-increasing orientation. For a quadratic edge with Bézier controls \(C_0,C_1,C_2\), the derivative is \(2[(1-t)(C_1-C_0)+t(C_2-C_1)]\). Positive relevant coordinate differences certify strict monotonicity on the whole edge.
5. Top and bottom span axis to rim. A sufficient separation condition is
\(\min\{Z\hbox{ controls on top}\}>\max\{Z\hbox{ controls on bottom}\}\).
This implies pointwise top-above-bottom on their common radial domain. It is deliberately stronger than necessary: overlapping global height hulls yield `UNCERTIFIED` even if the actual graphs are disjoint.

These conditions make the entire current boundary a Jordan curve, with top/bottom separated graphs joining the vertical axis and rim only at their endpoints. Together with strictly positive whole-element meridional determinants, they imply a global embedding of the meridional disk. Here is the degree argument, rather than an inference from a drawing:

Away from images of element edges, the map has isolated regular preimages and every preimage contributes \(+1\) to the Brouwer degree. The Jordan boundary has degree one inside and zero outside. Hence every regular interior image point has exactly one preimage and no exterior regular point has a preimage. Across a shared edge, continuity and positive orientation put the adjacent element interiors on opposite local sides of the shared regular trace. At an interior vertex, the oriented element sectors form a positive local fan: the local degree is a positive integer. A multiple fan or a second preimage would force degree at least two at nearby regular values, contradicting the global degree-one count. An interior point mapping onto the boundary would similarly create preimages just outside the Jordan domain. Therefore the piecewise map is locally open, globally one-to-one, and maps the compact reference disk homeomorphically onto the enclosed current disk.

Finally, \(R/r>0\) preserves positive radial interior and the distinguished axis. Revolving an injective meridional map preserves the azimuth of each positive-radius material point and gives an injective three-dimensional axisymmetric body, after the usual identification of all angular coordinates on the axis. This does not exclude contact with an external vitelline membrane omitted from the model.

The certificate's reference argument is affine P1 geometry, with P2 displacement. It does not certify an arbitrary curved reference element. Its global boundary routine is specialized to a flat rectangular reference meridian; a curved or nonvertical-rim body must receive its own appropriate boundary check.

## 5. Executed tests and reproducibility

Run:

```powershell
python -X utf8 controlled_furrowing/scripts/mixed_geometry_certificate.py
```

The retained receipt is `evidence/mixed_geometry_synthetic.json`. The suite passes **48 checks**: exact affine meridional and 3D volume-preserving examples; forty small nonlinear maps of the form
\(R=r(1+ar+bz), Z=z+cr^2+dz^2\); independent analytic determinant/hoop/J evaluation at a rational lattice; a strictly positive but conservatively uncertified map; an inverted-map rejection; a conforming two-triangle disk embedding; nonconforming-midpoint rejection; and retention of an arbitrarily small quadratic surface extremum. Analytic lattice comparisons check the implementation, while the Bernstein proof—not lattice sampling—makes the whole-element certificate rigorous.

Mixed-solver output certification is a separate replay. The saved synthetic `runs/mixed_checks/activated.json` has now passed whole-element certificates on all **36 quadratic triangles**, the exact Jordan-boundary check, and consequently the discrete embedding check; see `evidence/mixed_geometry_activated.json`. The receipt binds the precise input-file hash and verifier-source hash. Reproduce with:

```powershell
python -X utf8 controlled_furrowing/scripts/mixed_geometry_certificate.py --input controlled_furrowing/runs/mixed_checks/activated.json
```

The replay also independently computes exact attained extrema of each quadratic apical/basal \(Z(t)\) and of aligned material-column heights. For \(q(t)=c_0+c_1t+c_2t^2\), compare endpoints and \(t_*=-c_1/(2c_2)\) when \(c_2\ne0\) and \(0<t_*<1\). No arbitrary coefficient cutoff is used; even very small extrema are retained. Material-column extrema are not minima of Eulerian thickness, since deformed apical and basal radii can differ.

No additional solver output is certified merely because this synthetic suite passes. Full bounds and exact rational coefficients must be retained for every saved state that receives a certificate. Failed sufficient bounds must be reported, not converted into a physically inverted classification or silently replaced by quadrature checks.
