# Slides — Purification of Noisy Bell States

Beamer deck summarising the project: physical vs. virtual purification, the
three protocols (4Q / 5Q / 6Q), their noisy circuit decompositions, the exact
Bell-diagonal maps, the repeated dynamics, and the operational thresholds.

```
slides/
  pqec_thresholds.tex          the deck (12 content slides + title + outline)
  beamerthemeqitwhite.sty      minimal white-background / black-text theme
  figures/                     vector figures copied from the result folders
  pqec_thresholds.pdf          the built deck
```

## Build

```
latexmk -pdf pqec_thresholds.tex
```

Needs `texlive-latex-recommended`, `texlive-latex-extra`,
`texlive-fonts-recommended` (beamer, booktabs, helvet, microtype).

## Figure provenance

Every figure is the vector PDF produced by the study that owns it; none is
redrawn or re-exported by hand.

| file in `figures/` | source |
|---|---|
| `4q_circuit.pdf`               | `pqec-bell-distillation/results/figures/circuit_5cnot.pdf` |
| `5q_circuit.pdf`               | `.../results/figures/locc_vd/locc_vd_5q_baseline.pdf` |
| `6q_circuit.pdf`               | `.../results/figures/locc_vd/locc_vd_6q_circuit_by_party.pdf` |
| `type3_16cnot.pdf`             | `figures/step3_16cnot_circuit.pdf` (this repository) |
| `step3_threshold.pdf`          | `figures/step3_threshold_vs_epsilon.pdf` (this repository) |
| `bell_population_dynamics.pdf` | `.../results/report_figures/fig01_bell_population_dynamics.pdf` |
| `attractor_type4.pdf`          | `.../results/report_figures/fig02_common_attractor_type4.pdf` |
| `type5_metastability.pdf`      | `.../results/report_figures/fig04_type5_metastability.pdf` |
| `qth_curves.pdf`               | `.../results/qth/qth_curves.pdf` |
| `qth_hierarchy.pdf`            | `.../results/qth/qth_hierarchy.pdf` |
| `qth_binding_round.pdf`        | `.../results/qth/qth_binding_round.pdf` |
| `type4_14cnot.png`             | `figures/step4_resynthesized_14cnot_circuit.png` (this repository) |
| `type5_14cnot.png`             | `figures/step5_learned_14cnot_circuit_two_rows.png` (this repository) |

`step3_threshold.pdf`, `qth_binding_round.pdf`, `attractor_type4.pdf`,
`bell_population_dynamics.pdf`, `type4_14cnot.png` and `type5_14cnot.png` are
not used by the current slide set; they are kept here as ready drop-ins for
backup slides.

## Numbers quoted in the deck

| slide | source |
|---|---|
| 3 — ideal purification | `pqec-bell-distillation/README.md` §§9–11 |
| 4 — self-mitigating Fredkin noise | `RESEARCH_SUMMARY.md` §2, Step 2 |
| 5 — Type 3/4/5 decompositions | `RESEARCH_SUMMARY.md` §2, Steps 3–5; `docs/type4_7cnot_factorization.md` |
| 6 — 5Q ≡ 6Q | `docs/locc_vd_analytic_derivation.md`; `docs/5q_6q_type{3,4}_noisy_equivalence_numerical.md` |
| 7 — exact maps | `RESEARCH_SUMMARY.md` §§3.3–3.4; `pqec-bell-distillation/README.md` §13b |
| 8 — repeated dynamics | `RESEARCH_SUMMARY.md` §§3.3–3.5, §4.2 |
| 9 — resources | `RESEARCH_SUMMARY.md` §2 (ladder table), §4.1 |
| 10–11 — thresholds | `docs/qth_full_circuit_numerical.md`; `docs/4q_qrel_circuit_numerical.md`; `results/qth/` |
| 12 — open problems | `docs/qth_full_circuit_numerical.md` §16; `RESEARCH_SUMMARY.md` §7 |
