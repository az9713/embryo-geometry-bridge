# Project documentation requirements

## Explain the mathematical description of shape

The user requires every future scientific document about tissue/embryonic shape to either explain the shape representation itself or prominently refer to the canonical guide:

- Source: `docs/SHAPE-DESCRIPTION.md`
- Rendered document: `docs/SHAPE-DESCRIPTION.html`

Place the explanation or link before substantive shape-dependent arguments. Use the guide's reference/material versus current/spatial coordinate distinction. If another notation is used, give its mapping explicitly.

Every new model/result must state, directly or by a clear section reference:

1. What object describes shape: a body deformation, a parametric surface, a height graph, a mesh, or an observable.
2. The reference domain, coordinates, admissibility conditions and symmetry restrictions.
3. Which geometric quantities and observation/averaging conventions are actually computed.
4. Whether the result concerns stored fields, a finite-dimensional equilibrium, the continuum, or biological data.
5. Any additional regularity needed for curvature/bending, and any coarse-graining or wrinkle-scale assumptions.

An introductory reference may be concise, for example: "Shape convention: see How we describe tissue shape mathematically; this document uses its axisymmetric deformation map and a material-ROI mean height contrast."

Keep the guide consistent with implementation and identify changes from the existing axisymmetric model. Do not turn this explanatory requirement into a new approval gate. Compile equations for HTML and check the rendered page at desktop and mobile widths. Preserve the project's privacy boundaries and the instruction not to contact authors.
