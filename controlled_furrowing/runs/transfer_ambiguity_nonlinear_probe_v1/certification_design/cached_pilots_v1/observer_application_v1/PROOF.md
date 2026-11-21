# Cached-quadrature roots and static observer error

This is a candidate application of the independently reviewed root and geometry certificates. It awaits a separate readback.

The two large cached-quadrature roots have a positive signed current-radius contrast gap of about **0.1547551567 µm**. The exact application gives the outward interval **[0.15475515644140, 0.15475515688952] µm** for the positive-direction root minus the negative-direction root.
The two small L roots satisfy all four frozen primary comparisons at each root. The regional mean of total shortening exceeds the frozen allowance under radial-line weighting for both roots.

## Containment and coordinate error

The four certificates enclose unique local zeros of their respective exact-real frozen binary-cache quadrature potentials at the exact prescribed rational strength. The proof uses the same binary cache, preconditioner, free coordinate order, center and material tuple as the matching target. The exact center is the rational sum of reference coordinates and injected stored free displacement. All 650 coordinate hulls agree with the independently reviewed geometry box at radius \(2^{-36}\). The contraction posterior is strictly smaller than that radius for each target, so the geometry proof applies to each cached root.

The cached root error is \(e=\eta/(1-z)\). The exact maximum difference \(\xi\) between the rational center and the serialized curve nodes accounts for binary addition. The physical nodal error is \(\epsilon=L_0(e+\xi)\), with \(L_0=33.5\,\mu\mathrm m\). The radial coordinates at the axis and rim are fixed by the actual free map. Each base graph derivative has a strict positive margin over the reviewed \(8\epsilon\) bound.

| Root | Posterior error upper, dimensionless (approx.) | Binary addition gap, dimensionless (approx.) | Physical nodal error upper, µm (approx.) |
|---|---:|---:|---:|
| + small L | 1.2307616e-13 | 2.2204460e-16 | 4.1304899e-12 |
| + large | 2.1696956e-12 | 2.2204460e-16 | 7.2692242e-11 |
| − small L | 5.3509555e-12 | 2.2182776e-16 | 1.7926444e-10 |
| − large | 4.9103741e-14 | 2.1857516e-16 | 1.6522976e-12 |

## Large-patch contrast at matched optical dose

For \(Q=z_a(0)-z_a(a)\) and \(a=24.94698765\,\mu\mathrm m\), the fixed axis value is a P2 control node and contributes at most \(\epsilon\). At the same registered shoulder radius, the reviewed graph theorem gives \(E_a=\tfrac54(1+M_a)\epsilon\), where \(M_a\) bounds the slope of the stored apical graph. Each saved interval was inflated by \(\epsilon+E_a=(L_0+C_a)(e+\xi)\), and exact endpoint subtraction then gives the gap above.

| Root | Q error upper, µm (approx.) | Outward root Q interval, µm |
|---|---:|---:|
| + large | 2.1904786e-10 | [-4.13600326186264, -4.13600326142454] |
| − large | 5.0099958e-12 | [-4.29075841831406, -4.29075841830403] |

## Small L at the frozen allowance

The frozen conditional discrepancy allowance is about **5.38983567480 µm**. The two total-shortening comparisons each receive \(E_a+E_b\), while the apical and basal comparisons receive \(E_a\) and \(E_b\), respectively. A spatial maximum changes by no more than the sup norm of the graph difference. A normalized mean with common nonnegative weights has the same property. These bounds apply to the declared static operators on the same fixed radius region and unchanged reference. A comparison is IN when the maximum absolute residual is at most the allowance. It is OUT when the minimum absolute residual is strictly larger.

| Observer | + small L | − small L |
|---|---|---|
| Disk weighting, component maxima (primary) | 4 IN, 0 OUT | 4 IN, 0 OUT |
| Disk weighting, component means | 4 IN, 0 OUT | 4 IN, 0 OUT |
| Radial-line weighting, component maxima | 3 IN, 1 OUT | 3 IN, 1 OUT |
| Radial-line weighting, component means | 3 IN, 1 OUT | 3 IN, 1 OUT |

In each of the four OUT comparisons, the regional mean of total shortening is the failing observable. Its outward absolute residual interval is [6.07609207743, 6.07609207748] µm at the + small root and [6.19199948190, 6.19199948339] µm at the − small root. Changing component aggregation does not change that total mean.

The result covers only these two small states and two large states. The remaining 32 small endpoint states have no root-observer enclosure in this application. The root error accounts for the cached-quadrature algebraic error. Exact spatial-integral FE roots, stability, a continuous loading branch, PDE and source error remain open.

The exact [results](results.json) retain all 32 rows, endpoint arithmetic, proof hashes and both geometry and cache scopes. The implementation reuses the reviewed nodal-to-graph bound and saved-curve control extraction. No observer was re-integrated and no mechanical computation was run.
