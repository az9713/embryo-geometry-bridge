# Cold K/L surface-stress screen

The screen uses the already certified negative small_16 K and small_17 L exact-spatial roots. The K center is the certified recentered vector. The L center is the original singleton proof center. The exact proof posteriors enclose each root in its own free-coordinate box. Fixed axis and rim coordinates remain exact. Input authentication reuses the frozen active-energy input helper only, and binds each case to its separate independent whole-body and trace geometry certificate. It never evaluates K or L active energy.

At each ordered edge use its exact affine material radius \(R(t)\) and its P2 deformed coordinates. For \(K_s>0\) and \(T_s<K_s\), set \(\alpha=1-T_s/K_s>0\) and enclose the degree-six polynomial
\[
 D(t)=r(t)^2\big(r_R(t)^2+z_R(t)^2\big)-\alpha^2 R(t)^2.
\]
The screen uses exact rational interval coefficients and Bernstein hulls on the fourteen fixed interior windows of each of twelve edges on each face. There are 672 windows in total. If \(\overline D<0\), then \(J_s<\alpha\) and \(\tau_s \le K_s\overline D/(2\alpha R_{\rm hi}^2)<0\). The strict positive radial slope is checked separately. A nonnegative lower bound for \(D\) certifies nonnegative \(\tau_s\) on the window. A crossing enclosure remains unresolved. Where \(T_s\ge K_s\) or \(K_s=0\), nonnegative \(\tau_s\) follows directly under the admitted nonnegative coefficients. For a strictly negative window the receipt also retains a positive lower bound on \(-\tau_s r_R^2/J_s\) for the thin-layer argument.

The cold preflight retains six algebra and synthetic interval checks and the authenticated 26-file proof/input map. No actual surface stress was computed. A subsequent actual result is an interval sign diagnostic pending separate arithmetic review and a genuine process-exit observation. The derivation explains the conditional implication for admissible high-angular-frequency perturbations of the explicit natural 3D law. A negative window would concern this idealization at the revolved FE field. It would not pick a physical wavelength, locate a stable alternative shape, or certify a continuum equilibrium. A screen with no negative witness leaves the transverse stability question open.

The independent cold review file is `independent_prelaunch_review.json` with status `PASS_INDEPENDENT_COLD_TRANSVERSE_SURFACE_STRESS_SCREEN`, exact `source_sha256`, `protocol_sha256`, `preflight_sha256` and a nonempty dictionary of literal-true checks. After review and root authorization, run:

```text
python -B surface_stress.py run --protocol-sha256 <sha> --preflight-sha256 <sha> --review-sha256 <sha>
```

[Second variation and scope](READBACK.md) · [Frozen inputs and windows](protocol.json) · [Cold preflight](preflight.json)
