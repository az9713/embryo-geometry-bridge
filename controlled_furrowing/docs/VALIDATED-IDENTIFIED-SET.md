# Validated parameter cells: what an identified-set calculation must establish

**The inference layer now distinguishes certified parameter-cell statements from checked point predictions.** It retains every one of the 17×4 source constraints, carries unresolved cells and unbounded parameter tails, and separates a selected-root-branch exclusion from exclusion of all possible roots. No actual full physical identified set has been computed.

The [implementation](../scripts/validated_identified_set.py), [25-check synthetic receipt](../evidence/validated_identified_set_synthetic.json) and [actual 68-entry audit](../evidence/validated_identified_set_actual.json) are separate from all frozen mechanical solvers. No new equilibrium was solved.

## 1. The conditional inference object

Fix measured doses \(I_i\), an explicitly stated model/geometry/observation contract \(\gamma\), and a discrepancy vector \(\delta_k\ge0\). A model compatible set would contain parameter values satisfying every observation band,

\[
 |\mathcal O^\gamma_k(u_i(\theta))-y_{ik}|\le\delta_k,
 \qquad i=1,\ldots,17,\quad k=1,\ldots,4.
\]

The discrepancy vector is an input, not an estimated noise law or confidence level. Different \(\delta\) and different source-plausible observation contracts define different conditional sets.

If the mechanics has more than one root, the branch quantifier matters. This layer uses a tuple of explicitly certified local root domains, one per source record. It does not silently prove that these roots form one global continuation path, exhaust all roots or represent experimental temporal maxima. Packages and both reviews must explicitly use the supported branch law `independent_stationary_root_per_record`. An IN result is labelled `IN_STATIONARY_ROOT_TUPLE_WITNESS`: it is not physical inner membership for a different law that requires one common stable/quasistatic path. That stronger law needs its own proof adapter; no energetic stability is inferred.

## 2. What a usable forward cell certificate contains

For a closed rational parameter box \(C\), a point solution and a point-observer interval are insufficient. A usable root proof must establish, for every \(\theta\in C\), a root in its declared local branch domain. A usable observer proof must then enclose its source-matched endpoint uniformly,

\[
 \forall\theta\in C,\qquad
 \mathcal O^\gamma_k(u_i(\theta))\in[L_{ik}(C),U_{ik}(C)].
\]

The implementation requires hash-bound root and observer artifacts and independent reviews. It checks exact parameter-box containment, the source record and dose, the gain–dose relation, the branch-domain identity, and the model/observation contract. Each independent review must explicitly attest the uniform quantifier; the observer review must also attest admissibility of the full forward family.

The expected coordinate names for the current conditional slice are apical/basal area moduli, apical/basal tensions and common gain. Shear normalization, compression coefficient 5000 and geometry are fixed by the contract; they do not become inferred physical constants. A box within this slice is not the full [physical force-ratio and geometry domain](PHYSICAL-PARAMETER-DOMAIN.md).

Artifact status text and checksums do not create a proof. The referenced independent reviews are trusted proof-checker inputs. This module checks their identities and dataflow and applies exact interval logic; it does not replace the uniform mechanical/geometry/observer proof engines.

The representative strength-one pilot has declared dose one and no covered calibration record IDs. It cannot be relabelled as source column M. Nor can the existing one-point 611-DOF certificate be promoted to a uniform parameter-cell or 17-record forward proof.

## 3. Exact interval decisions

For one row, set \(y=y_{ik}\), \(\delta=\delta_k\) and \([L,U]=[L_{ik}(C),U_{ik}(C)]\). Define

\[
 d^-=\max(L-y,y-U,0),\qquad
 d^+=\max(|L-y|,|U-y|).
\]

Then every certified output has residual in \([d^-,d^+]\). The row decision is:

| Exact condition | Row status | What it proves |
|---|---|---|
| \(d^->\delta\) | OUT | Every parameter in this cell violates this row on the specified root family |
| \(d^+\le\delta\) | IN | Every parameter in this cell satisfies this row on that family |
| Otherwise, or no uniform proof | UNRESOLVED | The current evidence does not decide the row |

Every source row is listed even when no forward proof exists.

An **IN cell** requires IN for all 68 rows and uniform root existence/admissibility. One rigorously uniform contradictory row suffices for an **OUT selected-branch cell**, while the other 67 rows may remain unresolved. An OUT result for one local branch does not exclude another admissible root. Conversely, a fully compatible certified branch gives an existence witness; it does not prove every other branch fits.

The result explicitly records both the selected-branch decision and the existential parameter decision when all possible branches are allowed. An independent theorem valid for every admissible branch can provide a stronger exclusion, but it must be named rather than inferred from a failed chosen branch.

## 4. Exact covers and the unbounded remainder

Every split is rational. Splitting coordinate interval \([a,b]\) at \(a<m<b\) produces two closed children,

\[
 [a,m]\cup[m,b]=[a,b].
\]

The other coordinates are unchanged. Induction on the stored binary tree proves that the current leaves cover the original finite box exactly. Shared boundary points are intentional: they are not gaps, and a closed cell containing both compatible and incompatible points remains unresolved.

The cover ledger retains every parent, split, child and unresolved leaf. It rejects gaps, orphan/cyclic nodes, changed proof files, forged leaf decisions and mixed data/discrepancy/operator contexts. Accepted leaf decisions are replayed from their referenced proof artifacts. Resolved leaves without complete proof/context payloads are rejected, the root node is compared with the immutable declared root box, and the outside remainder is rederived from the declared domain. Exported snapshots are deep copies and have a standalone verifier; a mutable status flag cannot erase unbounded tails.

A finite root box leaves the exact symbolic remainder

\[
 \mathcal D\setminus C_{\rm root}.
\]

If the ambient domain \(\mathcal D\) has unbounded tails, those tails stay **UNCOVERED**. A numerical support limit is never converted into a physical prior. The physical quotient permits zero surface/gain coefficients; a positive/log-coordinate local seed does not remove those boundary faces from the ambient remainder. Other geometry, compression-ratio and branch possibilities are retained separately. Even a resolved finite synthetic cover does not mark the full physical identified set complete or prove a coherent global branch path.

## 5. A legitimate necessary exclusion from the actual source data

There is one useful parameter-independent check, provided its observation assumption is made explicit.

For signed apical and basal motions \(a(R),b(R)\) at the **same static state, current-radius ROI, reference and weighting**,

\[
 s(R)=a(R)+b(R),\qquad
 \langle s\rangle=\langle a\rangle+\langle b\rangle.
\]

If the recorded component endpoints are interpreted as these regional means, compatibility requires

\[
 |y_{\rm mean}-y_{\rm apical}-y_{\rm basal}|
 \le\delta_{\rm mean}+\delta_{\rm apical}+\delta_{\rm basal}.
\]

This is the triangle inequality applied to an exact linear identity. It holds independently of mechanical parameters and root branch under this declared static operator. It is **not** an equality for separate temporal maxima.

Source column G has

\[
 |y_{\rm mean}-y_{\rm apical}-y_{\rm basal}|
 =10.42530168\ \mu{\rm m}.
\]

Consequently a common discrepancy allowance must satisfy

\[
 \delta\ge 3.47510056\ \mu{\rm m}
\]

for that static regional-mean interpretation to be compatible. This is a necessary floor, not an estimated error or a sharp mechanical fit optimum. A value below the floor excludes that conditional operator model for every parameter and branch, not biology itself. The source audit explicitly retains time/aggregation/reference uncertainty.

For the spatial-maximum-component interpretation, the universal inequalities

\[
 \langle s\rangle\le\sup s\le\sup a+\sup b
\]

supply no positive discrepancy floor for these source entries. That does not establish a good mechanical fit or identify the correct source operator.

The actual audit reports a sample tolerance family \(\delta=0,1,\ldots,6\) µm, always retaining all 68 entries. These are illustrative conditional inputs, not physiological priors or a measured uncertainty law. That foundation audit deliberately supplies **zero** actual uniform mechanical cohort envelopes to test rejection of point-only evidence. It is retained as a regression receipt, not a current-coverage aggregate. A subsequent independently reviewed actual-column-M provider now supplies one uniform channel; the separate [actual cell application](ACTUAL-M-CELL-EXCLUSION.md) records its selected-family exclusion and all remaining gaps. No numerical cell decision is fabricated from the previous point certificates.

## 6. Retained verification and interface

Twenty-five synthetic checks include an exactly solved affine forward family, a resolved finite IN cover, mixed boundary cells, selected-branch OUT, one-channel exclusion with missing channels retained, unbounded tails, changed hashes/contracts, duplicate records, mixed discrepancy contexts and forged-status/proof replay rejection. These fixtures are labelled synthetic; they are not biological evidence.

The actual-data receipt passes six checks: it verifies all 17 records and 68 rows, rejects promotion of the real point-root certificate into a uniform cell, and retains zero-coefficient faces, unbounded tails and excluded-physics remainders. Its two seed-box children retain all 68 unknown rows at the separately declared common discrepancy 5.1 µm. The separate actual source-column-M package enters only after its **actual measured dose**, new-material branch, full root family and observer bounds have been uniformly certified and independently reviewed. A contradiction at declared \(\delta=5.1\) µm would concern that one selected-branch cell; it would neither estimate noise nor finish \(\Theta_\delta\).

Run the pure inference verification with:

    python -X utf8 controlled_furrowing/scripts/validated_identified_set.py test
    python -X utf8 controlled_furrowing/scripts/validated_identified_set.py actual

The full physical compatible set, tail exclusions across all admitted physics, source-matched temporal observations and uncertainty-propagated larger-region predictions remain incomplete.

The first inference-layer review found three mutable-ledger gaps: proofless resolved leaves, a changed root-node box, and a forged empty-tail flag. The prior source and receipts are preserved under `runs/validated_identified_set_v1/pre_review_revision_0/`. The corrected standalone verifier and all three negative reproductions now pass; no actual uniform mechanical provider had been accepted before the correction.
