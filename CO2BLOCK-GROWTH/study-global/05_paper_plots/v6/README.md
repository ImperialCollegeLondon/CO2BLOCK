# 05_paper_plots / v6 — full PNAS-style submission bundle

PNAS bundle for **"Geophysical and growth limits to the global scale-up of geological CO$_2$
storage"** (Gao, Kivi Rahimzade, Krevor). Same 5-figure **content as v5/v4-main** and the same
numbers (identical sha256-seeded common-random-number bootstrap, whole-period-feasibility
accounting, common country sets), but every figure is pushed to a **fuller, more informative
PNAS style**: the constraints each figure is judged against are now drawn *on the figure*, not
just implied. v6 supersedes v5 for the figures; v5 is kept as a fallback.

```
v6/
├── manuscript/  main.tex · si.tex · si_table_failed_cases.tex
├── code/        make_paper_figures.py
├── output/      figures/ fig{1..5}_*.{png,pdf} · tables/ table1 + fig_key_numbers.csv + fig4 bootstrap
└── si_figures/  figS1–figS7 .{png,pdf} + si_table_failed_cases.tex
```

## Regenerate

```bash
cd 05_paper_plots/v6/code
conda run -n pygis python make_paper_figures.py     # 5 main figures + Table 1 + key-number tables
```

SI figures S1–S7 come from `05_paper_plots/v4-main/code/make_si_figures.py`. Upstream reads:
`01_growth_model/output/v8_2026-06-01`, `02_co2block_screening/output`, `03_feasible_growth/output`
(nothing upstream is modified).

## What v6 adds over v5 (the "full PNAS-style" pass)

| Figure | PNAS-style enhancement |
|---|---|
| **Fig 1b** growth-form fork | shaded "both peak post-2100" quadrant + direct country labels on the four largest-demand jurisdictions |
| **Fig 2b** concentration | Lorenz **drop-lines** to the axes at the top-10 / top-20 / top-50 cut points |
| **Fig 3d** commitment collapse | faint **50% / 25% / 10%-of-goal reference guides** so the median-21% line has visual context |
| **Fig 4a** attribution matrix | per-cell **deliverable-ceiling line** (jurisdiction-summed 150-yr CO2BLOCK resource) — cumulative demand rising above it is what makes a scenario bind; delivery then flattens below it |
| **Fig 5** per-country grid | each panel is a **feasibility phase diagram**: 2-D KDE density clouds + horizontal deliverable ceiling (infeasible zone shaded) + vertical growth cap + goal-C → reconstructed-C arrows |

The shared design conventions stay: no baked-in "Figure N |" banners, no on-figure footnotes
(title + notes live in the LaTeX captions), short keyword panel titles, the unified Paul-Tol
scenario colour key, and a ≥ 6 pt legibility floor on primary labels (`FS_MIN` in the code).

## Main-text display set

| Item | File | Content |
|---|---|---|
| Fig 1 | `fig1_demand`     | a stress-pathway trajectories · b growth-form fork · c 2050 cumulative commitment |
| Fig 2 | `fig2_supply`     | a regional deliverable rate + per-basin spread · b concentration (Lorenz) |
| Fig 3 | `fig3_test`       | a runway · b shortfall severity · c capacity-used heatmaps · d commitment collapse · e retained share |
| Fig 4 | `fig4_outcome`    | a full 8-scenario demand-vs-delivered matrix · b targeted-vs-feasible ladder · c driver timing |
| Fig 5 | `fig5_commitment` | per-jurisdiction feasibility phase diagrams (KDE density + ceiling + growth cap) |
| Table 1 | `tables/table1_scenario_outcomes.{csv,md}` | scenario × (2050 rate L/G, feasibility reduction L/G, binding cases, earliest binding year, jurisdictions) |

## Headline numbers (verified against `output/tables/fig_key_numbers.csv`)

- 123 basins; 6,415 Gt deliverable over 150 yr; max aggregate rate ≈ 42.8 Gt yr⁻¹.
- 0 of 150 cases bind before 2050; earliest binding 2055; 39 of 150 bind over 2030–2179 (23 L, 16 G).
- Median shortfall in binding cases 59% (L) / 62% (G); feasible commitment retains a median 21% of design C.
- "Earliest binding year" monotone in ambition (2098 → 2069 → 2061 → 2055/2058), never before 2050.

## Semantics (identical to v4-main / v5)

Fig 4's "whole-period-feasible" keeps sampled demand futures whose lifetime commitment C is
geologically deliverable (2030–2179), every 2050 target held fixed; the one cell with no admissible
future (Canada × IPCC High × Gompertz, min C 49.9 > 38.8 Gt at g ≤ 20% yr⁻¹) is scenario-infeasible
and its target removed. This is an **admissibility** statement; 2050 delivery never binds before 2055
(Fig 3). Logistic-vs-Gompertz comparisons use common country sets (IPCC High n=6, IPCC Low n=7); the
collapse (Fig 3d/e) uses ascending-order reconstructions (descending differs ≤ 12%). Fig 3 =
**delivery**, Fig 4 = **admissibility** — never conflate the two.
