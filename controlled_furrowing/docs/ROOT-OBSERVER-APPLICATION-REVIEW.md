# Independent review of the root-to-observer application

**PASS for this one finite-quadrature field.** The independent review reconstructs the exact nodal error and its transfer to the four static endpoints. All 46 arithmetic, provenance and guard checks pass. This closes an algebraic-error link; it does not close a continuum or biological-validation gate.

The reviewed field is the original-material, 12×6, dimensionless activation-strength-one state in `runs/activated_quadrature/lambda_5000.0_strength_1.0/reequilibrated.json`. It is **not** one of the seventeen new-material cohort states. The [application receipt](../evidence/algebraic_root_observer_error.json) and [independent review receipt](../evidence/root_observer_application_independent_review.json) retain the exact fractions, source hashes, extracted surface edges and self-contained independent replay program. See the [explanation of the result](ALGEBRAIC-ROOT-OBSERVER-ERROR.md) for the evidence chain and scientific limits.

## What was checked independently

The review did not call the application `main()` or write its output. It executed only the source's read-only guard prefix, then reconstructed the calculation from the frozen cache and archived field. Eleven in-memory negative probes reject missing or failed reviews, nonboolean root checks, incorrect certificate/source bindings, an uncertified producer result, or stale transfer dependencies.

The root review binds the actual producer certificate and source. The application binds that review, the conditional transfer review, their current dependencies and the exact input state. The decoded cache checksum passes. Its reference coordinates and displacement vector agree bit for bit with the archived field. The 611 unique free coordinates independently match the prescribed axis-radial and fixed-rim constraints.

The producer and independent interval calculations certify the same target, center and uniqueness cube. Their root-distance upper bounds are approximately \(2.61989954\times10^{-12}\) and \(2.62544698\times10^{-12}\). The application takes the larger bound, denoted \(r_*\). These are coordinate **infinity-norm** bounds in dimensionless free displacement coordinates; no Euclidean norm conversion or factor \(\sqrt{611}\) is appropriate.

## Exact coordinate and unit conversion

The certified center uses exact binary reference-plus-displacement sums. The observer uses the rounded archived `current_nodes`. Recomputing both with rational arithmetic gives their maximum difference

\[
 g_{\rm add}=\frac{8131}{36893488147419103232}.
\]

With the declared physical scale \(L_0=33.5\) µm, every physical nodal-coordinate error relative to the archived observer controls is bounded by

\[
 \varepsilon=L_0\bigl(r_*+g_{\rm add}\bigr).
\]

The addition gap is included before the unit conversion. The surface endpoint radial coordinates correspond to global degrees of freedom \(12,180,0,168\). Each is fixed by the actual free-index map, and the exact center agrees there with the archived coordinate. Thus the unknown root and the archived surfaces share their global radial endpoints.

## Same-radius graph transfer

The review extracts the top and bottom P2 edges directly from the reference coordinates and six-node triangle connectivity. It independently computes the exact minimum of each edge's affine radial derivative and the maximum absolute affine vertical derivative. The strict radial margin exceeds \(8\varepsilon\) on every edge.

For quadratic nodal interpolation, coordinate error is at most \(e=5\varepsilon/4\), and radial-derivative error is at most \(8\varepsilon\). Consequently the perturbed traces remain strictly radial-monotone. The independently reconstructed global graph Lipschitz upper bounds are

\[
 L_a=\frac{240694177404370}{1248097794944047},\qquad
 L_b=\frac{1898428326274565}{80363328358385792}.
\]

They apply across edge junctions, not merely within the edge containing a selected radius. At any common current physical radius, the vertical graph errors satisfy

\[
 E_a=\frac54\varepsilon(1+L_a),\qquad
 E_b=\frac54\varepsilon(1+L_b).
\]

Both surfaces cover the same fixed positive-width observation interval

\[
 I=\left[0,\frac{337053732898705073}{36028797018963968}\right]\ {m µm}.
\]

Fixed references, the same interval and fixed normalized nonnegative weights imply component maximum or mean errors at most \(E_a,E_b\). A same-radius height maximum or mean has error at most \(E_a+E_b\). No uncertainty in the ROI, reference level or observation convention is included by this theorem.

| Endpoint | Outward-rounded upper bound, µm |
|---|---:|
| Apical component maximum or regional mean | \(1.311535134918811\times10^{-10}\) |
| Basal component maximum or regional mean | \(1.1254717312172188\times10^{-10}\) |
| Same-radius height-shortening maximum | \(2.43700686613603\times10^{-10}\) |
| Same-radius height-shortening regional mean | \(2.43700686613603\times10^{-10}\) |

All four exact rational bounds agree with the application receipt. The displayed binary floats are checked to lie above their exact bounds. The outward-rounding helper also passes independent tests at zero, ordinary inexact fractions and an underflow-scale positive rational. The initial nearest-rounding issue and its correction remain archived under `runs/root_observer_application/pre_outward_display/`.

## Scope of the verdict

This is a composition of a reviewed finite-dimensional root certificate and a reviewed exact trace-transfer theorem. The tiny error concerns proximity to one frozen binary-quadrature stationary point. It does not certify the integration rule's error, the continuum PDE, energetic stability, a unique physical branch, the full-body geometry of the unknown root, all seventeen material-cohort states, an identified material domain, source measurement semantics, or cross-size prediction. There is no claim of a new mathematical concept. The original five-gate objective remains incomplete.

The review is an independent computational check within the project, not external expert peer review.
