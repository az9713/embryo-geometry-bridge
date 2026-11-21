# First-wave exact-spatial observer transfer

This scalar application covers four small fields in the frozen order \(+01,+13,-01,-16\). It keeps each original archived field, the same current-radius ROI and flat reference, and the four independently replayed static observers. It uses the existing conditional discrepancy threshold, approximately 5.389835675 µm. The application has no large-field observable, refit, new force evaluation, AD, or integration.

## Premise and bound

The native independent endpoint replay and the root-observed process exit supply the only exact-spatial root authority. The application matches each root's box, model contract, cached point, preconditioner, exact drive, verification center, and original archived field. It independently checks the strict self-map and computes \(e=\eta/(1-z)\). Here \(e\) bounds the root's distance from the current verification center. The exact difference \(\sigma\) between this center and the original free vector remains in both recentered cases. The original binary addition gap \(\xi\) measures the difference between exact \(X+u\) and the serialized coordinates. The nodal bound in physical units is

\[
\epsilon = (33.5\,\mu\mathrm m) (\sigma+e+\xi).
\]

The calculation rederives \(\sigma\), \(\xi\), the two P2 trace derivative bounds, and fixed radial coverage from the matching archive and cube. For a stored face with minimum parametric radial derivative \(m_R\) and graph slope cap \(M\), require \(m_R-8\epsilon>0\). The sharp P2 value interpolation norm bound is \(5/4\) and its derivative counterpart is \(8\). The vertical Eulerian error of that face is bounded by \(5(1+M)\epsilon/4\), including the change in parameter to match the same current radius. The two face errors add for total shortening. Each separate component uses its corresponding face error. Taking the maximum or integrating against fixed normalized nonnegative weights preserves the resulting uniform error bound. [Reviewed nodal transfer](../../../../../../docs/NODAL-OBSERVER-ERROR-TRANSFER.md)

Each of the four archived fields supplies 16 intervals, from four static conventions with four channels apiece. After enlargement, an observation comparison is IN if the upper endpoint of its absolute residual is at most the threshold, OUT if its lower endpoint exceeds the threshold, and UNRESOLVED otherwise. Cases retain the literal archived measurements and the original frozen thresholds. The 16 primary and 48 alternative comparisons cover these four small-field roots only. The common threshold is a model-discrepancy condition. Experimental noise, source operator equivalence, continuous loading, continuum accuracy, and full calibration remain outside this application.

## Evidence and execution

[Cold protocol](protocol.json) · [Preflight](preflight.json) · [Source](apply.py) · [Observed independent replay exit](../first_wave_reduction_v1/independent_actual_preparation/actual_replay/os_exit_observation.json) · [Independent replay batch](../first_wave_reduction_v1/independent_actual_preparation/actual_replay/receipt.json)

The preflight verifies the native and archived case links and exercises rational scalar examples. Actual execution requires a separately pinned independent source review over the frozen source, protocol, preflight, observed native replay and the four exact case IDs. The application writes its result and real process metadata under a fresh `actual/` directory. The comparisons then await an independent scalar readback.
