# Interval polynomial cleanup: critical failure audit

The old `exact_fe_target.clean` keeps a coefficient only when `v != 0`.
For the installed Arb binding, a ball enclosing zero and nonzero values has
both `v == 0` and `v != 0` false. These comparisons assert certified relations;
they are not logical complements for uncertain intervals. Consequently the
old cleanup discards an uncertain coefficient as though it were exactly zero.

For example, `arb(0,1)` contains the admissible exact coefficient 1/2.
`T.constant(arb(0,1))` returns the empty polynomial. Evaluating that output gives
zero, which does not enclose 1/2. Addition and multiplication reproduce the
same enclosure failure. This is an arithmetic soundness failure, not an API
style issue, a sampling limitation, or harmless loss of interval dependency.

The bounded audit also observes the defect in the actual new uniform element0
field/jet construction: ten uncertain coefficients disappear across nine
cleanup calls. The first omitted coefficient encloses
`[-96439867/35184372088832, 96439867/35184372088832]`, approximately ±2.75e−6.
No spatial integral or auxiliary64-site sum was run in this audit.

At512bits the complete degree16 constitutive-series construction on actual
POINT elements0 and6 did **not** encounter an uncertain coefficient discarded
by cleanup. This does not verify the full point computation. It is deliberately
reported as a negative observation, rather than claiming demonstrated corruption
of every historical point result. All full point integration certificates are
nevertheless withdrawn from reliance until a safe full replay or an equally
complete no-loss trace establishes valid local enclosures. The independent
point implementations reuse the same polynomial arithmetic and cannot rescue
the old certificates merely by agreeing with them.

In particular, the old complete-point center-escape claim is not currently
certified. Its exact signed reduction of the stored168 interval packets is
still reproducible arithmetic, but its inference requires those packets to be
valid integral enclosures. Neither hashes nor overlap between two implementations
sharing the same unsound primitive supplies that missing premise.

The independent bulk-cache AD reviewer has a separate `if v != 0` cleanup and
therefore needs correction and replay as well. The bulk-cache producer's dense
analytic third tensor and three-term mean-value identity do not perform sparse
Arb coefficient cleanup. The current surface-cache jet uses `not(v == 0)` and
passes the zero-crossing retention counterexample; archived earlier revisions
used the unsafe predicate. This is a scoped finding, not a blanket claim that
all surface-cache mathematics has now been independently reverified.

Unaffected by this particular failure are the exact rational moment formulas,
the rational coefficient/tail calculations, dense binary-quadrature root and
energy-stability proofs, and the independent tuple-interval whole-body geometry
calculations. The retube vector arithmetic remains arithmetic about the chosen
stored vectors. Its description as an integral Newton correction inherits the
withdrawn point defect and must be withheld. Fresh geometry around that chosen
arbitrary center remains a geometry theorem; it never needed equilibrium
membership. Conditional transfer and retube lemmas remain mathematical lemmas,
with no current validated integration-error premises.

The safe replacement must discard only a coefficient certified exactly zero
(`is_zero()` for Arb), and reject nonfinite balls. Fraction cleanup remains
ordinary exact arithmetic. Required tests include zero-crossing and exact-zero
balls through all sparse operations and jets, not just fixed positive input
fixtures. The receipt lists the minimal replay sequence and the hash-bound
affected/unaffected inventory. Historical sources and receipts were preserved.

Execution: twelve bounded reproducer checks and eight final scope checks pass.
The first audit incorrectly required actual POINT drops already at degree4;
that assertion failed. Its source remains in `engineering_revision_0`.
The corrected audit measures this fact without assuming the outcome.
The seven-element uniform pilot has no independent launch approval.
