# What an annulus can constrain without knowing the slice height

**Status: a conditional geometric observer is implemented and synthetically verified; a quantitative Movie 3 test is not performed.** The [movie data audit](MOVIE-GEOMETRY-FEASIBILITY.md) found no defensible anatomical inner/outer contour intervals. A dark fluorescence region is not automatically an empty tissue region, and inventing a pixel allowance would not supply the missing measurement error model. This component does not replace the 31-embryo larger-region cohort or calibrate the mechanics on an observed ring.

The earlier unknown-plane binary observer could often return both possible labels. Quantitative ring dimensions can carry more information: the same unknown horizontal plane must intersect both ring boundaries. The equality of the two boundary heights eliminates that nuisance height without reconstructing a three-dimensional object. This is elementary level-set geometry and standard interval exclusion, not a novel mathematical or biological concept.

![Synthetic annulus compatibility and exact pair curve](../evidence/annulus_observer_visual.svg)

The example is entirely synthetic. An annulus exists in the model, yet a particular inner/outer-radius pair can still be impossible for every slice height. This is the additional information lost by reducing an image to a binary ring/no-ring label.

## Primary-source audit before implementation

Contour-constrained reconstruction has a long history. [Meyers, Skinner and Sloan (1991)](https://graphicsinterface.org/wp-content/uploads/gi1991-32.pdf), printed pp246–247, define contours as surface/plane intersections and discuss the underdetermined correspondence problem. This establishes relevant ambiguity rather than a unique reconstruction from one image.

[CrossSDF (Walker et al., 2024), Section 3.1](https://arxiv.org/html/2412.04120v1#S3.SS1), takes cutting planes and segmented contours as input to reconstruct an implicit three-dimensional surface. Its plane geometry is supplied; our observer instead keeps a horizontal height unknown and tests a necessary relation. We do not import its neural reconstruction machinery or claim that sparse-image uncertainty is solved.

[Syha et al. (2013), Section 2.4](https://journals.iucr.org/j/issues/2013/04/00/ks5358/ks5358.pdf), is a directly relevant validation example outside biology: it estimates unknown cutting-plane orientation before comparing two-dimensional sections with a reconstructed material volume. Its Discussion also explains why near-tangent interfaces amplify contour shifts. [Klauschen et al. (2009)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3474339/) reconstruct cell/tissue surfaces from segmented three-dimensional microscopy and discuss lower out-of-plane resolution. Neither paper establishes that our composite movie supplies an anatomical surface contour.

Search strings, access depth, exact locators and failed accesses are retained in `evidence/annulus_observer_literature.json`. Three additional author/publisher PDFs had a size-limit, timeout or HTTP403 failure; their theorems are not silently treated as fully reviewed. This is a bounded audit supporting a standard observer, not a systematic novelty clearance.

## Conditional observation model

Let an axisymmetric body have continuous basal/apical graphs \(b(r)<a(r)\), for \(0\le r\le R\). A horizontal plane of unknown height \(\sigma\) sees ideal tissue occupancy

\[
 S_\sigma=\{r:b(r)<\sigma<a(r)\}.
\]

Assume the observed region corresponds to a central occupied disk, one empty annulus \((r_i,r_o)\), and occupied farfield. For an **apical-void** annulus the plane lies at or above the apical surface in the gap. Continuity and positive thickness imply

\[
 a(r_i)=a(r_o)=\sigma,\qquad 0\le r_i<r_o\le R.
\]

At its inner/outer boundaries, the plane exits and re-enters through the apical graph. This equality is not valid merely because an image contains a dark ring. A basal exit would instead constrain \(b\). A tilted plane has angle-dependent height, and axial fluorescence blur need not preserve solid occupancy. These are different observation models.

No time is supplied by the geometry. An equilibrium matching some slice is not a prediction of when a transient movie ring appears or recovers.

## Necessary interval test with a justified surface-error assumption

Suppose physical-radius intervals \(I_i,I_o\) enclose the anatomical boundary radii. Let the model graph \(\widehat a\) satisfy an assumed vertical error bound at fixed Eulerian radius,

\[
 \sup_r|a(r)-\widehat a(r)|\le\delta.
\]

Define its height ranges \(A_i=\widehat a(I_i)\), \(A_o=\widehat a(I_o)\). Any compatible common slice must satisfy

\[
 \sigma\in\bigl(A_i+[-\delta,\delta]\bigr)
           \cap\bigl(A_o+[-\delta,\delta]\bigr).
\]

If this intersection is empty, the data/model pair is excluded **under the stated horizontal-plane, apical-void, radius-error and surface-error assumptions**. A nonempty intersection means only that this necessary test did not exclude it. It does not prove the predicted section is an annulus or fits the image.

The bound \(\delta\) must be justified separately. It is not inferred from force stationarity, a mesh difference, fluorescence threshold sensitivity or an arbitrary pixel allowance. A Hausdorff surface bound is not automatically a vertical graph bound without suitable transversality/slope control. Similarly an unknown plane tilt or point-spread function cannot be silently replaced by a chosen numerical tolerance.

If the central/farfield tissue occupancy is independently justified at witness radii \(r_k\), an optional basal graph with the same uniform error allowance gives further necessary strict bounds:

\[
 \widehat b(r_k)-\delta<\sigma<\widehat a(r_k)+\delta.
\]

The strict endpoints matter: a slice exactly on an apical surface is not interior occupancy in this idealization. The implementation preserves open/closed interval endpoints with exact rational arithmetic. It never interprets an invalid or uncovered radial domain as a measured biological exclusion.

## Exact P1 pair curves, including topology

For a continuous piecewise-linear apical graph, height ranges are computed exactly from interval endpoints and interior profile vertices. A stronger zero-surface-error algorithm enumerates the full one-annulus family. On a descending inner edge and an ascending outer edge, their inverse radii are affine in slice height:

\[
 r_i(\sigma)=\alpha_i\sigma+\beta_i,\qquad
 r_o(\sigma)=\alpha_o\sigma+\beta_o.
\]

Each pair is restricted to the overlapping edge-height ranges and any supplied radius intervals. It also enforces:

- every profile vertex before the inner root lies strictly above \(\sigma\), so the central region has no additional apical hole;
- every vertex between the roots lies at or below \(\sigma\), so the proposed annulus contains no occupied island;
- every vertex after the outer root lies strictly above \(\sigma\), so the farfield has no additional apical hole;
- \(r_o(\sigma)-r_i(\sigma)>0\);
- if a basal graph is supplied, its entire maximum is strictly below \(\sigma\).

These are finite affine inequalities in \(\sigma\), with exact interval endpoints. For a P1 graph, extrema on each interval occur at vertices or interval endpoints, so the checks hold over whole segments. Positive graph thickness is verified on the union of apical/basal breakpoints.

Why is the basal condition necessary? Inside an apical-void annulus, \(b<a\le\sigma\); outside, occupied tissue requires \(b<\sigma\). Thus the full section has \(\max b<\sigma\). Conversely the listed inequalities and \(\max b<\sigma\) produce exactly the specified central disk/one-void-annulus/outer-ring section, up to the boundary circles. Without a supplied basal graph the result is labeled apical-only and does not certify actual tissue occupancy.

The enumeration covers all such P1 sections. Each genuine inner crossing has a descending edge and each outer crossing an ascending edge. Flat portions at the slice height are retained through their neighboring crossing edges; a flat valley at a critical height is not dropped. Duplicate descriptions at shared-edge heights are allowed and explicitly reported. Two separate void annuli are not mislabeled as one.

For each branch the central and enclosed-outer disk areas are

\[
 A_i(\sigma)=\pi(\alpha_i\sigma+\beta_i)^2,\qquad
 A_o(\sigma)=\pi(\alpha_o\sigma+\beta_o)^2,
\]

and the empty annulus area is their difference. The implementation returns the exact rational quadratic coefficients before the factor \(\pi\). Surface-error intervals use the conservative necessary-range test; the exact pair-curve claim applies only to the supplied zero-error P1 geometry.

## P2 graphs: rational radial-root enclosures

The P2 mechanical model supplies boundary segment controls in the order start/mid/end. Each segment has exact quadratic polynomials \(R(t),Z(t)\), \(0\le t\le1\), when the stored binary coordinates are interpreted as rationals. The observer first verifies radial monotonicity using the endpoint derivatives of the quadratic radial map.

For each physical radius interval, exact rational bisection encloses the corresponding parameter roots. Expanding to the outer parameter brackets can only enlarge the admissible set. The height range on that enlarged parameter interval is obtained from its endpoints and the exact critical point \(t_*=-Z_1/(2Z_2)\), if present. No small-curvature cutoff is used.

This gives a rigorous enclosure of the height range of the **supplied discrete quadratic geometry**, without floating root inversion or treating nodes as a polyline. The two enlarged height ranges feed the same necessary common-level test. It is deliberately conservative and need not reconstruct the full P2 pair curve. It does not enclose the unknown biological surface unless a separate model-error bound supplies that link.

## Unknown image center and equal-area radii

Areas are invariant under in-plane translation. For an ideal axisymmetric horizontal section, the central component has area \(\pi r_i^2\), and the area enclosed by the outer gap boundary is \(\pi r_o^2\). Anatomical area intervals therefore imply equal-area radius intervals. If only annulus area is measured, it supplies one difference constraint, not two independent radii. Additional independent central/enclosed-outer areas or a correlated two-area feasible set would be needed for a sharper test.

Observed anisotropic shapes need not be circles. Replacing their full contours by area-equivalent radii is a scalar relaxation: it ignores anisotropy and spatial arrangement, making it easier for an axisymmetric model to pass. It does not establish circularity, the correct center, anatomical segmentation or horizontal-plane alignment. Shape information may separately reject an axisymmetric model even when its scalar areas pass.

The conversion itself is rigorously enclosed. Machin's identity \(\pi=16\arctan(1/5)-4\arctan(1/239)\), with alternating-series remainder bounds, supplies rational \(\pi_-<\pi<\pi_+\). For area \(A\in[A_-,A_+]\),

\[
 \sqrt{A_-/\pi_+}\le r_{\rm eq}\le\sqrt{A_+/\pi_-}.
\]

Integer square roots supply outward rational square-root enclosures. The API returns both the area and \(\pi\) intervals. Exact arithmetic certifies the conversion, not whether the source measurement interval is justified.

## Implementation, tests and actual data limit

`scripts/annulus_observer.py` contains the standalone standard-library implementation. Inputs may be decimal strings or exact stored binary floats; both are converted to `Fraction`. Rational outputs retain numerator/denominator and a display float. Invalid nonmonotone quadratic graphs or uncovered radius intervals raise an error rather than silently becoming scientific exclusions.

```python
from annulus_observer import p1_necessary_test, p1_pair_branches
top = [(0, 2), (1, 1), (2, 0), (3, 1), (4, 2)]
bottom = [(0, -1), (4, -1)]
test = p1_necessary_test(top, ('0.4', '0.6'), ('2.9', '3.1'))
family = p1_pair_branches(top, bottom)
```

The example's exact family is \(r_i=2-\sigma\), \(r_o=2+\sigma\), \(0<\sigma<2\). The synthetic radius box above is excluded although the model can produce annuli. Adding an assumed surface error of 0.2 removes that exclusion; this demonstrates dependence on the error assumption, not an empirically justified rescue.

Run `python -B controlled_furrowing/scripts/annulus_observer.py`. The retained `evidence/annulus_observer_verification.json` records 62 synthetic checks, source hash and examples: compatible/incompatible pairs, larger-error relaxation, strict occupancy endpoints, exact P1 branches, flat valleys, multiple annuli, P2 nonlinear radial inversion, 39 rational-root checks, translation, area/π/square-root bounds and invalid-geometry rejection. The original synthetic SVG and its hash are in `evidence/annulus_observer_visual.svg` and `annulus_observer_visual.json`.

The [source movie audit](MOVIE-GEOMETRY-FEASIBILITY.md) remains decisive for empirical use. Its visible moat does not supply certified anatomical boundary areas/radii: dim labeled membranes cross apparently dark regions, the boundary is not a uniquely resolved continuous contour, and raw z/point-spread/segmentation error information is absent. Consequently **no real-data annulus exclusion or successful annulus prediction is claimed here**. The verified observer is ready for appropriately specified future data; the biological measurement and original prediction gates remain unresolved.
