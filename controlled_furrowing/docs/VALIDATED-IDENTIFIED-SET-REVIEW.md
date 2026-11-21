# Independent review of the identified-set inference layer

**PASS: 55 independent analytic/adversarial checks and 22 actual-provider/cover checks.** The actual-M cell is isolated inside the broader finite seed by nine exact cuts. Its selected root family is OUT at the declared 5.1 µm tolerance; nine other finite leaves remain unresolved. This is a software/mathematical-composition review, not a new mechanical forward proof or completion of the physical identified set.

The [review receipt](../evidence/validated_identified_set_independent_review.json) binds the engine, its residual-calculation dependency, review source and every independent fixture. The [review script](../scripts/review_validated_identified_set.py) writes only its own fixtures under `runs/review/identified_set/`. It does not run the producer's fixture generator or modify its outputs, and performs no physical-model optimization.

## Exact interval decisions

For observation \(y\), allowed discrepancy \(\delta\ge0\), and a uniform prediction enclosure \([a,b]\), the sharp absolute-residual range is

\[
 d_{\min}=\max(a-y,y-b,0),\qquad
 d_{\max}=\max(|a-y|,|b-y|).
\]

The rule is OUT if \(d_{\min}>\delta\), IN if \(d_{\max}\le\delta\), and unresolved otherwise. In particular, touching the allowed closed band is not OUT; interval overlap alone is not IN.

At the cohort level, IN requires all four declared endpoints for every declared source record. For this actual cohort that is all 68 constraints. One reviewed uniform contradiction suffices for OUT on its selected root family even when the other 67 constraints have no prediction. Their absence remains visible. The independent fixtures verify both cases, including exact threshold boundaries and a missing 68th channel.

An IN package supplies, for every common parameter vector in the cell, independently certified stationary roots for each source record. This yields a compatible **stationary root tuple**. It does not prove that the roots lie on one stable, dynamical or quasistatic path. An OUT result on one such family does not exclude other stationary branches. The engine preserves these distinctions in separate status fields and rejects unsupported common-path package claims.

## Uniform proof provenance

The engine checks the requested cell is contained in the certified outer cell, the universal root quantifier, exact source/dose association, branch domain, model/observation contract hash, and independent root and admissible-observer audit bindings. A point-only root certificate supplies no uniform bounds. Duplicate record entries invalidate the entire record rather than silently selecting one entry.

The independent suite exercises corrupted certificate hashes, empty/failed/nonboolean review checks, wrong cell and contract hashes, missing universal quantifiers, wrong doses, diagnostic inputs, unsupported branch coverage, observer branch mismatches, unknown channels and repeated records. It also independently enumerates all vertices of a two-parameter affine model and checks the prediction/residual arithmetic for all 68 channels. Record and package permutations preserve the decisions.

These audit files are a **trust boundary**: their hashes bind specific mathematical attestations. Arbitrary self-authored PASS JSON is not a mathematical proof or an authenticated reviewer identity. The actual-M audits have separate numerical and proof reviews; synthetic audit fixtures are explicitly synthetic.

## Conditional static identities

For signed static graphs on the same state, ROI, references and normalized weights, mean height shortening cannot exceed its spatial maximum. With spatial-maximum components,

\[
 \overline S\le\max S\le\max A+\max B.
\]

With regional-mean components, instead,

\[
 \overline S=\overline A+\overline B\le\max S.
\]

Measurement discrepancy sums give the corresponding necessary interval constraints. A contradiction can exclude every parameter and every branch **under that declared static operator premise**. It is not a deduction that the experimental source uses that operator. Separately temporally maximized means do not satisfy the same-state mean identity in general.

The generic library receives an operator enum and an opaque contract hash separately; it does not parse the latter to establish their semantic consistency. The caller must bind the selected enum to the reviewed model/observation contract. A mean-operator contradiction must not be applied to a spatial-maximum model. This precondition is part of this review, rather than a claim that the API automatically validates scientific semantics. The actual-M contract is explicitly the static spatial-maximum variant; its current one-channel application uses no universal static-operator shortcut.

## Cover correctness and retained historical failures

Each subdivision uses exact rational closed intervals and proves that its two children have the declared parent as their union. Shared boundary faces are retained. Tree traversal rejects gaps, invalid parent links, orphan nodes and malformed structure. Resolved leaves require a complete bound inference context and proof payload; replay recomputes their decisions against the current hashed artifacts.

The three historical cover gaps are independently retested: a proofless resolved leaf is rejected, a root-node mutation cannot silently shrink the declared root, and an unbounded complement cannot be forged EMPTY. The suite also tests live mutation of root, ambient-domain and uncovered-physics declarations, forged leaf status, mixed discrepancy contexts and standalone exported-ledger replay.

Standalone replay proves coverage of the **declared snapshot**. Binding that declaration to an earlier requested root, observations or protocol requires its externally retained hash. An internally consistent change of the entire declaration describes a different task; self-consistency is not historical authentication.

Zero coefficient faces and unbounded tails remain outside a positive finite seed unless explicitly covered. Resolving all leaves of a finite seed does not establish the full physical identified set. Every exported verification keeps global branch/path certification and full physical completion false.

## Actual provider and larger-seed cover

The actual provider's two audit adapters are exact copies of the frozen independent root and observer audit sections, with transparent derivation metadata. The certificate, canonical cell and model/observation hashes are checked. Its source-M dose and single apical endpoint are carried into all 68 original constraints; the other 67 stay unresolved.

Every declared tolerance-family decision is cold-replayed. Independently computed exact residual cutoffs reproduce the strict OUT transition: a tolerance below the lower residual bound excludes this selected family; touching that bound is unresolved, and a tolerance large enough to satisfy this one row still cannot produce full IN. The current actual application supplies no static-operator enum, so the universal mean-identity shortcut is disabled. Its reviewed contract is explicitly the static spatial-maximum-component variant. No mean-operator contradiction is used to exclude the spatial-maximum model.

The broader cover has 19 tree nodes, nine exact cuts and ten finite leaves. Independent rational parent/child comparisons establish the closed-union induction. The sum of the exact five-dimensional leaf volumes also equals the seed volume. Target leaf `0100101010` is exactly the certified parameter cell. Only that leaf receives the reviewed uniform interval and selected-family OUT decision. The other nine leaves have all 68 constraints unresolved; the tiny certificate is not extrapolated to any larger leaf.

Standalone replay recomputes the complete exported cover. Additional mutations of this **actual** export—proofless target OUT, a shrunken root node, and a forged EMPTY tail—are rejected. Every leaf retains the same inference context and leaves other stationary branches unresolved. The ambient domain includes nonnegative zero-coefficient faces and unbounded tails; its complement of the finite seed remains explicitly UNCOVERED.

The exact coverage statement therefore is: one selected source-M local stationary family is excluded on one tiny parameter cell, inside an incompletely resolved finite seed. It cannot supply full-cohort IN, exclude other branches, resolve the ambient remainder, or complete the physical identified set.

The current review binds the engine **and** `joint_calibration_contract.py`, whose endpoint list and exact residual helper the engine imports. Future reuse must preserve that dependency closure and the caller's semantic contract, rather than trusting a source filename or an old PASS label.
