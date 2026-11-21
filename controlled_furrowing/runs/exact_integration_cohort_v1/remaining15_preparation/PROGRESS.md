# Remaining fifteen integral cases: execution checkpoint

The original five-gate objective remains active and incomplete. B and L have completed exact-integral root proofs; this directory concerns the other fifteen cases. A complete integral manifest is not a root certificate.

## Live queues — reuse the same handles

- Lane 1: exec session **9167**, supervisor PID **31904**, launch directory `../remaining15_integrals_preparation/launches/891e67eafef520601cf9fc500d00de82802f71b9a9b57b27a38011cd06e73fbd/`. Record order C,E,G,I,K,N,P,R; each bulk then surface.
- Lane 2: exec session **82124**, supervisor PID **23956**, launch directory `../remaining15_integrals_preparation/launches/a327e059d47fcfa01e6c1c1cc0ca891c1bbe47963e520941dcd0f3ccbf13f015/`. Record order D,F,H,J,M,O,Q; each bulk then surface.

Root subsequently observed C,D,E,F,G,H finish their bulk and surface parts, with both sessions still live. Next bulk workers were I/PID24192 and J/PID32592. Only C/D have completed the independent direct reduction described below; the newer completed parts remain unreviewed. Treat worker IDs as timestamped observations; poll the same handles and inspect actual worker and OS terminals for current state. Do not restart a queue because a tool observation expires. Each frozen queue starts one child at a time, for two operator workers total.

## C and D: complete integration, direct proof still inconclusive

Both records have 144 bulk and 24 surface packets, true worker exits and matching supervisor OS exits. The input-pin observation is `../remaining15_components_observed/CD_terminal_observation.json`, SHA `95697af788d6f746f236f634c06a103c0c5ecc0e00428b52e80817e6e804e94e`.

Generic reducer source `8e3177434299fbffeee2fa5161dd4025fc0e9cd0b7f46be85fceaeccae204ed2`; independent prelaunch `66c9a028bc45c204d5c889c4c698fef28511c63ab9df6bddf0cd02e8452d425a`. Actual reductions and independent endpoint replays are terminal. They retain 611 signed force/preconditioned-force rows and 6676 upper Hessian entries per record.

| Record | Reduction result SHA256 | Independent actual review SHA256 | Direct self-map upper bound | Original radius |
|---|---|---|---:|---:|
| C | `a38575bb4c1cb54d6ef295561828444a0fcdf4d51492d9bfbc40d5052f196f49` | `b66b821666c74ee68d2c0105df3087c4bec9e3a408d92f14cfa88231c9fc7d08` | about 7.108e-6 | 1e-11 |
| D | `6824d229c3f43527fcc50e6185d2567e57fc62c9a84a7bbc51e4ff1f243b99d5` | `03017ddf88803327d2bf24cb4d22c8b4957ae4969d999ed062c66924d1e5b411` | about 2.747e-5 | 1e-11 |

Results are under `../remaining15_reduction_preparation/RECORD/results.json`; admissions are sibling `RECORD_input_admission.json` files. Every signed preconditioned defect interval contains zero, so these data give no positive defect-norm lower bound and do not exclude roots. The independent width attribution `../remaining15_reduction_preparation/independent_actual/CD_width_results.json`, SHA `faf52e086282f86ad2632b6f932d9a438442ef86916aaf890188d8f764e386f6`, assigns more than 99.99999999 percent of the widest interval to off-axis contributions. This is enclosure-width attribution, not physical error.

## Next proof layer

The demonstrated L remedy is the next numerical step: evaluate the integral-minus-cache force at each record's own exact predictor anchor and base gain, use that record's reviewed uniform Hessian for the nodal segment, and retain the exact affine gain derivative. Preserve signs and shared parameter offsets before preconditioning/norms. If separate cached and defect norms remain inconclusive, reconstruct the signed full cached predictor residual and combine its anchor and slopes before taking the norm. Full-box derivative and predictor-only residual scopes stay distinct.

The record-neutral centered helper is frozen as preparation only: `scripts/exact_centered_cohort_bulk.py` SHA `524c22266a7e7e8138c06c8e4b444f44bac887d8404296d67d58a3541726d02c`, protocol `dbe461e20104cf60cee9de78010d2c37d086663c420e3fb47869ce5ad2ea339d`, preflight `9c4ea9ed77c22fda2a2b36ada0bc7ac3c1101f100f8d07d68a9a3bbb30af790c` under `../remaining15_centered_preparation/`. It binds 144 C and 144 D element proofs; the other 13 contexts are static routing checks, not reads of unfinished operator payloads. Module preparation review passed, `../remaining15_centered_preparation/independent_review.json` SHA `ffcebc86827a4e12cdbec8e0a7fd2b5be84e884eb83f88078d9e61356e7b3769`. Actual caller/supervisor approval remains separate. No centered C/D force computation or new root proof has yet been asserted here.

The future actual caller must set flint precision to 640 and clear `core.A` **before** `bind_reviewed_record`, use its unchanged authenticated Domain, and derive each actual local node-ball fingerprint from that Domain against the frozen box/route proof. Comparing only two stored metadata hashes is insufficient. The state segment must stay inside the reviewed Hessian domain; an arbitrary enlarged or mutated Domain is not covered. The C/D preparatory factory Domains were authenticated. Fresh domain/derivative work is mandatory for any translated cube escaping the original proof domains.

The scientific requirements remain separate: physical compatible-set coverage, source-matched reference/ROI/timing/uncertainty, continuous loading paths, continuum and nonaxisymmetric accuracy, informative larger-region source scoring and originality. No author contact or GitHub publication occurred.

## Subsequent verified progress — 2026-11-20

E, F, G and H now have terminal direct reductions and independent exact endpoint replays, not only completed integration packets. Frozen external observation `../remaining15_components_observed/EFGH_terminal_observation.json` SHA `acfae3d618422088740364a543b3d9cb2da674d7be0a2cecac914f671b817b9c` precedes their payload admission. All twelve admission/reduction/replay child OS exits were zero; the retained controller readback is `../remaining15_components_observed/EFGH_reducer_os_receipts/batch_readback.json`, SHA `b94b3ac90ff505073564cdc574a91fda7c5700ba9aa608955aef693f72383843`.

All four direct self-map tests remain **UNRESOLVED**. Their contraction bounds are below one, but force-enclosure widths prevent the self-map criterion. All 611 signed preconditioned defect intervals contain zero in each case. Successful arithmetic replay proves neither root existence nor nonexistence. Independent review hashes: E `36ae6c5709d33d41fb99fa637444aaadeffa396501b668248a644cc8283c1153`; F `7001bc8b2fe8d100096e8116ca7f644658688cbe5d3c709d0368791cfeb5a755`; G `80e53b127cbb6a7327dcbd92aff7e43c212a481a657e28bd0054293facacfa93`; H `bf195907b8ea34b953971a6ad5d38b5f3c8ba2ac35412e0318c0b03583febacc`.

The actual C/D centered-force caller now has independent PRELAUNCH approval `../remaining15_centered_actual_preparation/prelaunch_independent_review.json`, SHA `31ed5e5a15e73242414a44ed3b9412679d01d0e590eedbf4d94973dcbb556ca2`. It reconstructs 1,300 exact coordinates and all 288 local fingerprints at 640 bits before component calls. Domain growth, foreign records, wrong precision, failed children and PID mismatches are rejected. This is preparation, not a new centered-force result. Root has not launched those jobs while both original integral lanes occupy the two-worker cap.

Root's later live poll of the SAME session9167 observed I surface followed by K bulk/PID31204; session82124 remained live. Inspect current receipts rather than treating these timestamped observations as a current completion inventory.

Scientific discovery work is separate from these proof jobs. The bounded next-action memo `../transfer_ambiguity_next_v1/README.md` and receipt (`7e75e8a3f558393b7865fb0ba73d4c38fc62647e9ec1815fa5efb8b5814ab64a`) identify a full-numerical-rank but weak five-coordinate small-data tangent, and a one-large-field tangent plus signed local affine optimization as the next outcome-free diagnostic. New source `scripts/transfer_ambiguity_linear.py` is prepared under `runs/transfer_ambiguity_linear_v1/`; consult its actual terminal and independent review before claiming any calculation. Its exact affine LP certificates, if obtained, do not certify nonlinear compatible witnesses or experimental uncertainty.
