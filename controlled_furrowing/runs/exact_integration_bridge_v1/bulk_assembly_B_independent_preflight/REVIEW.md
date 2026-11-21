# All144bulk point assembly: independent prelaunch review

**PASS154checks.** Frozen assembly source41fc1735566f72ec42e249aa00fbc61f9e26b9ef55ea96f3e790bb82adac83ff; protocol8b14c671ca877d96256225014dcc42fef97df6ea6f7cabc3af860b857398cb26; preflight9ffb66e73f52770825cc288163d42ae7491687f915a430edeaab630ef32a919f. The launch receipt is `bulk_assembly_B/prelaunch_independent_review.json`, SHA2567e6cd01a7faed124008d5d562c66674bd1f2cf4fe426e013013bc8f47fea0aed. No actual144force/Hessian integrals were executed by this reviewer.

## Independent selector recomputation

For all144elements, the reviewer rebuilds the canonical fields and constructs the Jacobian jet from signed cubic monomials, using the independent implementation from the four-element replay. It does not call the wrapper's selector or local-rho function.

Off axis, multiply J−1by the exact radial monomial needed to remove its negative powers. Expand r=rc+h*x by explicit binomial coefficients. On the containing rectangle |x|≤1, |zeta|≤zmax, bound the numerator by the sum of absolute scaled coefficients and divide by rlo^d. This is an exact rational upper bound, even if very conservative. On axis charts, explicitly binomial-expand around(1/2,1/2)and bound each shifted coordinate by1/2. Taking the minimum with the separately reviewed global21/500bound is valid because both bound the same exact point field on the complete element.

The positive continuous mass is independently computed from the simplex affine-r moment det(B)Σr_vertices/6, with the retained normalization and mathematical2pi. The discrete mass is the exact normalized binary weight sum. Reviewed scalar series-tail bounds and independent Ju/Juv coefficient majorants give conservative maxima over all local force and Hessian entries. The reviewer checks every allowed order4…24and confirms the **smallest** accepted degree and exact rho, coefficient bound, free map and tail budgets for every element. All144declared rational force/Hessian tail upper bounds were reproduced exactly. Selected degrees are5through16.

The targets1e−22and1e−17apply only to omitted constitutive-series tails. They do not bound the total integral-minus-cache discrepancy, coefficient arithmetic width or finite cache correction. The local kernel retains those other contributions separately.

## Assembly and norm logic

Local free indices need not be globally sorted. Each local upper-triangular Hessian entry maps to one sorted global pair; symmetric row sums then count the off-diagonal entry once in each of its two rows. This is the correct assembly for a symmetric energy Hessian. Local gradient contributions scatter add to shared global free coordinates before multiplying by the exact binary611×611preconditioner.

Independent toy blocks use unsorted local-to-global lists, shared nodes, signed force/Hessian entries and nontrivial off-diagonal values. The sparse result agrees entrywise with a separately assembled dense double loop. The signed RF vector and row norms agree, and the conservative inequality ||RH||infinity≤||R||infinity||H||infinity holds. Outward union of serialized exact interval endpoints is also checked. Taking upper endpoints before maxima avoids ambiguous comparisons between overlapping Arb balls.

The wrapper's post-run result must remain an **unreviewed point-bulk assembly** until independently assessed. The all144sum includes no surface discrepancies or parameter/root-tube uncertainty. It establishes neither equilibrium persistence nor an exact-integration stability theorem by itself.

## Retained review correction

An initial independent assertion contained an incorrectly hand-entered expected toy RF vector. Every actual144selector check already passed. The failed checker and receipt are preserved under `toy_expectation_revision_0`; after correcting that expected vector, the complete154-check run terminated with exit0. No producer source or mathematical result changed.

The reviewed local pilot, exact moment/tail component, target polynomial algebra and independent cubic-jet helper are reused with explicit hashes. This prelaunch review approves the frozen bounded computation; it is not an independent reimplementation of every previously reviewed local kernel or a completed original five-gate project.
