# Independent material and initial-geometry inputs

**Conclusion: these sources justify an apical/basal asymmetry hypothesis and several geometric scales, but do not independently determine a finite, exhaustive parameter family for the Countryman experiment.** A declared exploratory family is possible. Calling it a measured interval, an identified set or a confidence region would overstate the evidence.

The analyst already knows the published results. The numerical constraints in this audit come from the independent Cheikh2023 study, Countryman's experimental protocol, and explicitly pre-activation Fig2a—not from tuning parameters to larger-region shortening or bending labels. No held-out response is fitted here.

## Primary sources and exact locators

* Cheikh et al., [eLife12:e85569 (2023)](https://doi.org/10.7554/eLife.85569), [PMC10584372](https://pmc.ncbi.nlm.nih.gov/articles/PMC10584372/), the actual Countryman reference56. Retained full primary XML: `literature/material_constraints_cheikh2023.xml`; publisher material is CC0. Inspect Results `s2-5` paragraphs1–2, `s2-8` paragraphs3–5, Fig4b, Fig7c″, Fig2—figure supplement1, Discussion `s3` paragraph10, and Appendix `s11` paragraph7 / `s13` paragraph1. The eLife webpage returned a client challenge and EuropePMC XML returned503; official NCBI S3 objects succeeded. Exact access paths/hashes are retained.
* Countryman et al., [Nature Communications16:7665 (2025)](https://doi.org/10.1038/s41467-025-62483-6), Methods `Sec14` (imaging/activation), `Sec15` (segmentation), `Sec16` (reference-corrected shortening), and Fig2a (explicit pre-activation image). The source article/figure license is CC BY-NC-ND; the downloaded Fig2 is a local inspection source, not an original project figure.

No wing-disc basement-membrane constants are used as embryonic basal-actomyosin parameters.

## What Cheikh actually estimates

Cheikh's cells have apical elastic sheets, elastic lateral structures and open basal surfaces with stiff basal edges. Their interiors/yolk are viscous fluid. This differs from an isotropic neo-Hookean bulk plus two area-elastic membranes. The experimental pulling examples have15–20µm lateral membrane lengths at mid-cellularization (Fig2—figure supplement1), whereas Countryman activates cells32–35µm long at gastrulation onset. Same species does not remove this developmental and constitutive mismatch.

The source gives an **order-of-magnitude apical2D modulus** of approximately0.001–0.01nN/µm, or10⁻⁶–10⁻⁵N/m, not a3D Young modulus in pascals and not a statistical interval. Its fitted spring examples are:

| Fig7c″ case | Apical k (nN/µm) | Lateral k | Basal k | Basal/apical spring ratio | Predicted maximum applied force (nN) |
|---|---:|---:|---:|---:|---:|
| Sim1 | 0.007 | 0.007 | 0.221 | 31.57 | 0.20 |
| Sim2 | 0.007 | 0.007 | 0.442 | 63.14 | 1.33 |
| Sim3 | 0.007 | 0.042 | 2.275 | 325 | 2.81 |
| Sim4 | 0.007 | 0.126 | 6.825 | 975 | 6.93 |

The experimental maximum-force range displayed in that comparison is5–11nN. These are **co-varied fitted discrete spring parameters**, not four measured physiological states or a measured basal/apical modulus interval. Results `s2-8` paragraph5 and Discussion paragraph10 explicitly leave lateral-versus-basal elasticity unresolved. The illustrative ratio975 is not a transferable continuum stiffness ratio, nor a validated lower or upper bound for Countryman.

The source's sheet relation is

\[
E_{2D}=\frac{2k}{\sqrt3},\qquad \nu=\frac13.
\]

That relation applies to its triangulated elastic sheet. Basal edges have a different network topology; assigning the same formula to a floppy basal network and dividing by the apical result is unjustified. A large spring constant and a large homogenized tangent modulus are different statements when geometry/topology changes.

Even for a conventional isotropic2D membrane, conversion to an area modulus requires constitutive assumptions:

\[
K_{\rm area}=\frac{E_{2D}}{2(1-\nu)},\qquad
\mu_{\rm surface}=\frac{E_{2D}}{2(1+\nu)}.
\]

Atν=1/3 these are3E₂D/4 and3E₂D/8. An area-only membrane omits this independent shear response. These formulas explain a possible model mapping; they do not prove that the proposed embryo membrane has those properties.

For a bulk/membrane model, the relevant ratios include

\[
\alpha=\frac{K_a}{\mu h},\qquad
\beta=\frac{K_b}{\mu h},\qquad
\pi_{a,b}=\frac{T_{a,b}}{\mu h}.
\]

Cheikh does not supply a Countryman-compatible elastic bulkμ, independent basal area modulus, or membrane prestress. Its cytoplasmic viscosity is not an elastic shear modulus. Consequently even an assumed transfer of the apical surface scale does not establish numerical bounds onα, β orπ.

## Prestress and basal support

No numerical pre-activation cortical tension or baseline pressure sufficient for the new constitutive model was identified in these sources. Cantilever loading forces are applied perturbation forces, not a measurement of resting prestress. Similarly, an ablation recoil reduction under RNAi does not isolate a stiffness ratio without dissipation, geometry and prestress information.

Zero prestress can be a declared reference-state null; it must not be called a measured physiological value. A stiff basal network also does not imply a vertically clamped basal surface: stiffness restricts deformation under force, whereas a displacement clamp prescribes position. That boundary condition needs an independent mechanical justification.

## Height, activation depth and cell diameter

Countryman directly specifies cell height32–35µm, activation depth5–25µm below the apical surface,1µm activation-stack spacing, and2µm imaging-stack spacing. These are usable protocol constraints. The exact experiment diameters are retained in the data audit, separately from nominal18/50µm labels.

Neither inspected primary method gives a per-specimen cell-width interval for the Countryman calibration cohort. Fig2a's pre-activation image visibly contains cellular outlines, and the18µm activation occupies a small group of cells, but the compressed illustrative image does not provide a validated cell-area distribution or a recorded segmentation uncertainty. No numerical cell-diameter prior is manufactured from an informal visual count. A biological cell diameter is also not the diameter of a cantilever tip or ferrofluid droplet.

Therefore this audit cannot establish the homogenization conditionℓcell/a≪1 for the small activated patch. A finite-thickness continuum calculation remains an explicitly conditional tissue model; a cell-resolved realization would require further initial cell-geometry information.

## Curvature: an independent scale calculation, not specimen bounds

Cheikh Fig4b uses full long/short model axes511/184µm. Assuming a prolate spheroid with semi-axesA=255.5µm andB=92µm, the equatorial principal radii are

\[
R_{AP}=\frac{A^2}{B}=709.568\ \mu\mathrm m,
\qquad R_{DV}=B=92\ \mu\mathrm m.
\]

These are derived from that paper's ideal model geometry, not measured principal curvatures in Countryman's particular pre-activation specimens. A perivitelline offset, stage differences and coverslip flattening further limit direct transfer.

For a circular section of radiusR and half-widtha, the sag isR−√(R²−a²). At the reported large-region half-width24.94699µm, the corresponding osculating-circle sag is approximately0.439µm alongAP and3.447µm alongDV. At small-region half-width9.35512µm, the values are0.0617µm and0.4769µm. The exact prolate-ellipseAP section gives nearly the same local scale (0.4396µm at the larger half-width).

This calculation shows why a flat initial condition can be a reasonable **test assumption** yet remain consequential for a near-apical optical slice: the independently estimatedDV geometric scale is comparable to or larger than the2µm Movie3 slice-depth descriptor. It does not prove that the actual experiment had that sag, or that2µm is a classification threshold. The size/curvature mechanism itself is already discussed in the source research; this arithmetic is not a novelty claim.

An axisymmetric patch cannot simultaneously represent unequalAP/DV principal curvatures. Flat and spherical-cap cases may be useful conditional sensitivities, but their numerical endpoints are not measured confidence limits on the actual anisotropic embryo.

## How to use this audit in a finite model family

1. Use the directly specified height/depth/diameter values, retaining nominal-versus-reported distinctions.
2. Label transferred apical stiffness scales and ideal embryo curvature as external, stage-mismatched assumptions. Retain units and constitutive conversions.
3. Treat basal area stiffness, basal/apical ordering in the chosen continuum law, prestress and outer support as hypotheses to vary or calibrate with small-region data. Do not use the975 spring ratio as a measured continuum parameter.
4. Include appropriate weak/strong-support or flat/curved limiting tests. Do not call a convenient finite grid exhaustive when no source supplies a finite upper bound or a justified prior distribution.
5. If fitting the small responses leaves different larger-shape outcomes possible, report those conditional predictions. Do not choose the branch agreeing with published larger-region labels.

The strongest justified numerical family here is an **explicit exploratory sensitivity ensemble**, with some protocol inputs fixed. An independently bounded material family requires further measurements or a demonstrated small-data identifiability result. This limitation does not forbid building and testing a coupled hypothesis; it limits what its result establishes about biology.

Source hashes, exact locators, stage/unit caveats and the reproducible geometric arithmetic are retained in [the material-input audit](../evidence/material_input_audit.json) and `literature/material_constraints_derived_numbers.json`. No held-out outcome fit or external contact occurred.
