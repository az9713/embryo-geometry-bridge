# Point-integration delivery approval withdrawn

The final-delivery scientific approval is **WITHDRAWN / UNVERIFIED pending corrected arithmetic and replay**. No independent approval of the historical 422-binding delivery has been issued. The original five-gate project remains incomplete.

## Reproduced defect

The frozen `scripts/exact_fe_target.py`, SHA-256 `4fd89ddf23c59050a92f18f7b4b6e69ef8e5d08e57ba0c1949cf0b61b2c869a4`, line 34 contains:

```python
def clean(p): return {tuple(k): v for k, v in p.items() if v != 0}
```

With the installed python-flint, `x = arb(0, 1)` is not exact zero: `x.is_zero()` is false. Nevertheless, both `x != 0` and `x == 0` are false because the ball overlaps zero. The actual frozen `clean({(0, 0): x})` returns `{}`. Replacing a nonexact zero-crossing coefficient ball with exact zero need not enclose the original polynomial. This is an arithmetic-enclosure defect, irrespective of whether the omitted radii happen to be numerically small.

The executable counterexample and five passing reproduction checks are retained in `record_invalidation.py` and `../../../evidence/exact_integration_final_delivery_independent_review.json`. Receipt SHA-256: `78b8975570ccde7a94065f9fc0fced7961742702e2503d93a01eaa4db5262691`. A preliminary recorder import-path error occurred before receipt creation; its diagnostic is preserved separately in `initial_root_resolution_failure.txt`.

## Affected approval and common dependency

The point-integration polynomial jets call this shared arithmetic. Our six-edge independent surface replay also imports it and invokes its `constant`, `add`, `mul`, and `scale` operations. Higher order, higher precision, or overlapping producer/reviewer intervals do not remove a shared unsound cleanup operation.

Consequently, the point B force/Hessian integration-discrepancy enclosures and their dependent negative direct-self-map numerical certificate require corrected, versioned replay. The historical complete-point review (`bc22a23e7af96d0d1a99836401112378600c14be61d5cba957e6f15a12d699c5`), documentation review (`30acb4991fbe695f661f68275e88af017ca12b328a96e538e2b254375565a21c`), and final-delivery receipt (`63236576ffc4e518a4b80448978b3dd562a5078f560c2efbeb5c1e034decfc50`) remain provenance, not current evidence that this numerical claim is certified. The documentation's logical reduction remains conditional on valid numerical premises; those premises are now unverified.

## Scope and recovery requirement

This finding does not establish that every reported number is false, that the integral model lacks a root, or that unrelated results fail. Exact-Fraction body geometry/moment identities and old cached-Q dense root certificates are not invalidated solely by this polynomial cleanup defect; each claim must be traced to its actual arithmetic dependencies.

Recovery requires retaining nonexact zero-crossing coefficients, dropping only proven exact zeros, testing this distinction directly, and replaying every affected polynomial enclosure and dependent assembly/certificate under a versioned source closure. No corrected numerical conclusion is asserted here. This review performed no new integrations, equilibria, publication, browser action, or Space edit.
