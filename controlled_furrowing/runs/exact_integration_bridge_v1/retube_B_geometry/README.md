# Fresh B common-domain geometry: preparation only

The new wrapper is `../../../scripts/exact_recentered_geometry_B.py`. Preparation freezes a **new** common-domain nodal box; it does not run the actual body or boundary certificate. Old geometry certificates were not enlarged or relabeled.

## Domain and exact containment

Use the original B seed/cache, unchanged predictor S, unchanged exact parameter cell C/base parameters, and the actual binary-applied decrement d from the independently reviewed center preparation. Let

\[
p_0(\theta)=u_0+S(\theta-\theta_0),\qquad
p_1(\theta)=p_0(\theta)-d,
\]

\[
r_{\rm common}=r_{\rm old}+\|d\|_\infty.
\]

The radius is exactly188144378943289/7036874417766400000000000, approximately2.6736924346449912e-11. It is a new numerical proof-domain radius, not a physical parameter, measured uncertainty or established root distance.

The unchanged `cohort_static_observers.nodal_box` is called with the ORIGINAL seed/cache/S/C/base and this new explicit radius. Its exact center-to-serialization gaps and predictor hulls give650 coordinate intervals. Fixed coordinates retain zero perturbation, including radial axis and boundary endpoints.

`coordinate_boxes.json` retains every exact interval and2600 rational containment checks covering:

1. the original predictor plus old-radius cubes;
2. the translated predictor plus old-radius candidate cubes;
3. all paths p0(theta)-t*d for0<=t<=1;
4. those same translation paths swept by an old-radius cube.

For each free coordinate, the swept interval has lower endpoint p0+hull-low+min(0,-d)-r-old and corresponding upper endpoint p0+hull-high+max(0,-d)+r-old. The common radius contains it because abs(d-coordinate)+r-old<=r-common. Source gap/symmetric-hull padding can only enlarge the factory box further. These are exact inequalities, not sampled path tests or reuse of the failed old boxes.

The independent candidate audit had already found396 old-box escapes. This wrapper does not assume those old boxes cover anything new.

## Frozen inputs and preparation checks

`protocol.json` binds the retube independent review, original sources/arrays, original B state, candidate arithmetic and vector files, recursive local runtime dependencies, the common observation contract, and all old numerical evidence required by those packets. R, S, C, base parameters, geometry and source dose are checked unchanged.

`boundary_preparation.json` contains original baseline P2 controls, boundary node IDs, verified fixed radial endpoint IDs, the existing physical length67/2 micrometres and the SAME current-radius ROI. Boundary-specific nodal error budgets are computed from the new common coordinate radii. This only prepares inputs; `curve_bound` has not been called on the actual boundary.

`synthetic_and_containment_tests.json` records24 true checks, including signed/asymmetric scalar hull fixtures, review-schema guards, all650 boxes/all2600 inclusions, and an AST guard requiring the actual `body.certify_box` call to appear only in the gated run function.

## Gated future body and boundary run

A NEW `prelaunch_independent_review.json` must have status PASS, nonempty literal-True checks and exact current source/protocol/preflight hashes. Only then may the run function:

1. Call the unchanged `root_body_geometry.certify_box` on the new full nodal box, preserving canonical P1 material corners and exact P2 current coordinates/fixed binary offsets.
2. Check body embedding/axis/strong monotonicity and the stronger analytic condition J-lower>1e-6. Compute a NEW rho from the returned global J interval; never copy the old0.042 value.
3. Run `curve_bound` separately on apical and basal baseline controls with the fresh boundary-specific error budgets and proven fixed radial endpoints.
4. Certify original baseline Eulerian thickness using `positive_thickness`, then subtract both fresh graph-error bounds to obtain uniform separation on the common ROI.

Successful body monotonicity with lower constant m>0 and a fixed radial axis implies r/R>=m and dr/dR>=m on boundary traces. Consequently Js>=m squared and G>=m to the fourth power. The run records this positive-bound corollary only when its body preconditions hold.

This corollary does **not** provide a uniform upperG or surface binomial expansion rho. A NEW interval-coefficient trace wrapper is required for those bounds; the existing `positive_G_bounds` is point-only. Canonical body geometry also does not establish the separate cached/convex-strain-line analytic domain required by the uniform cache-consistency primitive.

## Scope

The future body certificate would concern every canonical P2 field in the supplied coordinate box, independent of any root existence claim. Membership of a new exact-integral root requires the missing fresh uniform Hessian/parameter-discrepancy/selfmap proof. No source observation alignment, experimental uncertainty, PDE accuracy or biological trajectory follows from geometry alone.

No actual full-body/actual boundary certificate, H assembly, integration, Newton solve, root proof, fitted parameter change, production-document edit or author contact occurred during preparation.
