# Independent rendering-delta addendum

The rendering regression review passed 38 checks across six local fixture types at desktop1440 and mobile390 widths in Microsoft Edge157.0.4341.60 with KaTeX0.16.11. The original builder/verifier and original180-check review remain preserved byte-for-byte. Current scientific proofs/data were not changed or recomputed.

## Exact defect and fix

The fixture table contains `\(|s|\le0.001\)`. Previously Markdown tokenized the literal `|` characters inside the generated `data-tex` attribute as table separators. The old rendered table retained three nominal cells but displayed a broken span fragment in its condition cell. Crucially, the old browser probe reported math0/rendered0 and no error: its count-equality test passed despite the missing equation.

The current converter HTML-entity-escapes attribute pipes before Markdown table parsing. Browser parsing decodes the entities back to the exact original TeX `|s|\le0.001`. Both tested widths produce one math node, three correct table cells, one KaTeX render and the exact TeX annotation. The screenshots were visually inspected; the mobile table fits the viewport.

- [Current desktop screenshot](new_table_1440.png)
- [Old broken desktop screenshot](old_table_1440.png)
- [Current mobile screenshot](new_table_390.png)
- [Fixture source](table_fixture.md)
- [Retained review code](review_rendering_guard.py)

## Raw prose guard

The new verifier's exact browser probe rejects raw `\rho` and `\frac` in prose. The old probe would accept that same zero-math-node page. Explicit code examples remain literal and exempt; script/style/math blocks are also excluded. Standalone inline and display equations retain their exact TeX and render successfully.

This is a backslash-letter-control-word heuristic, not universal detection of every possible raw mathematical expression. A deliberately retained limit fixture demonstrates that plain algebra or control-symbol-only delimiters can evade it. Correct equation-node conversion and independent content roundtrip tests remain necessary.

## Binding delta and scope

The old180-check receipt is unchanged. Exactly one of its20 publication-source hashes is intentionally superseded: `controlled_furrowing/scripts/build_report.py`. This addendum/receipt records its old and new hashes. The new verifier hash is separately bound. All71 scientific combination inputs and all18 physical/observer kernel source hashes still match their retained records.

The body of the builder is unchanged in AST except its `convert` function. Exact old/new source snapshots and unified diffs are retained. No production source file, scientific Markdown, figure or certificate was edited by this reviewer.

These fixtures do not constitute the final102-page/204-viewport dossier run, live deployment verification or new scientific proof. Those require separate final build/browser receipts after the three authored Markdown formatting repairs settle.
