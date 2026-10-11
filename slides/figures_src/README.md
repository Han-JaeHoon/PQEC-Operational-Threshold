# Protocol circuit figures

Single-source TikZ drawings of the three purification protocols, in one
visual language: heavy black line art, Computer Modern for symbols,
Helvetica bold for annotation text, transparent background.

```
circuitstyle.tex        shared macros (wire, gate, control, swap cross,
                        meter, trash can, register brace)
fig_4q_protocol.tex     4Q  postselected comparator     -- full
fig_5q_protocol.tex     5Q  shared-ancilla SWAP test    -- full
fig_6q_protocol.tex     6Q  distributed LOCC-VD         -- full
fig_{4,5,6}q_compact.tex   the same circuits, trimmed for a three-up slide
```

Build one figure:

```
pdflatex fig_4q_protocol.tex
pdftocairo -png -r 600 -transp -singlefile fig_4q_protocol.pdf out
```

Outputs are exported to `../figures/protocols/` as `protocol_*.pdf`
(vector, for Keynote) and `protocol_*.png` (600 dpi, transparent).

The three `*_compact` figures share an identical bounding-box width, so
they line up when scaled to the same size on a slide.

## What each figure asserts

| figure | content |
|---|---|
| 4Q | the five CNOTs of `CNOT_SEQUENCE`; the lower pair `(A_2,B_2)` is measured and postselected on `(0,0)`; `rho^2/Tr rho^2` is left physically on `(A_1,B_1)` |
| 5Q | one shared ancilla, two controlled-SWAPs; the ancilla is read in `Z`, the retained pair `(A_1,B_1)` in the intended observable `O`, the second copy discarded |
| 6Q | wire order `a, A_1, B_1, A_2, B_2, b`; Alice's CSWAP on `{a,A_1,A_2}`, Bob's on `{b,B_1,B_2}`; the dashed `S^dagger` marks the second measurement setting -- both the `XX` and the `YY` run are needed for `Delta(O) = C_XX(O) - C_YY(O)` |
