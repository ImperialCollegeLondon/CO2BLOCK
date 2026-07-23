# Quest CCS Reservoir Properties — Adversarial Review

**Review date**: 2026-07-23
**Document reviewed**: `CANADA-CCS-ALBERTA-Quest-props.md`
**Reviewer**: automated review

---

## 1. Source Tier Classification

The document lists 10 sources but does NOT explicitly tier them. Reclassification:

| Source | Tier | Rationale |
|--------|------|-----------|
| Gen-4 Report, Shell (2011) — Open Alberta | **T1-Reg** | Formal regulatory filing, Shell internal reservoir model, most detailed single source |
| D-65 Application, Shell (2011) | **T1-Reg** | Regulatory filing containing actual well-test and core data |
| AER Closure Plan (2023) | **T1-Reg** | Regulatory document, authoritative for SLA definition |
| AGS PRS 2024-001, Herbers et al. | **T1-Gov** | Alberta Geological Survey publication, independent government source |
| Rock et al., GHGT-13 (2017), Energy Procedia | **T1-PJ** | Peer-reviewed journal, Shell authors with direct data access |
| Lahvis et al., Energy Procedia (2017) | **T1-PJ** | Peer-reviewed journal, independent authors |
| Rock et al., IJGGG (2022) | **T1-PJ** | Peer-reviewed journal, BUT same author group as GHGT-13 — not independent |
| Shell & ESG, GHGT-15 (2021) | **T2-CP** | Conference presentation (less rigorous review than journal); Shell authors |
| Rock & O'Brien, AAPG #80577 (2017) | **T2-CP** | Conference presentation; same Rock author |
| ESG Solutions webpage | **T3-Web** | Industry website, cites GHGT-15, no original data |
| carbonstorage.io | **T3-Comp** | Compilation website, original sources unverifiable |

**Issue**: Sources 1, 2, 4, 6, 10 (Rock et al. × 3, Shell × 3, carbonstorage.io → GHGT-15) are **not truly independent** — most trace back to the same Shell internal dataset. Cross-validation looks stronger in the table than it is.

---

## 2. Reservoir Thickness — Review

### Reported values
- 40 m (GHGT-15), 40–48 m / mean 44 m (Gen-4), ~45 m (IJGGG), 40–80 m (AGS), 38 m net sand (AAPG), 45.7 m (carbonstorage.io)

### FINDINGS

**GREEN: Good convergence** — 5 independent-ish sources cluster at 40–45 m for the Quest site. The AGS value (40–80 m) is correctly identified as regional.

**YELLOW: Two concerns:**

1. **What does "40 m" mean?** GHGT-15 states *"bottom 30m perforated for injection."* The actual **net pay** for injectivity is 30 m, not 40 m. For storage capacity calculations, using 30 m (net-to-gross adjusted) may be more appropriate. The document recommends 40–45 m without explicitly addressing whether this is gross or net.

2. **38 m (net sand) vs 40 m (gross)** from different sources. The N/G ratio of 0.9 *almost* reconciles this (40 × 0.9 = 36 m), but this consistency check is not made explicit.

3. **carbonstorage.io value (45.7 m = 150 ft)** — this is presented alongside primary values but is a tertiary source. The units suggest it was originally reported in feet, implying a pre-metric Shell report. This should be flagged as "derived from an unverifiable compilation."

**Recommended fix**: Make net-vs-gross explicit. State: "Gross thickness 40–45 m; injection interval (perforated) 30 m; net sand ~38 m (N/G 0.9 consistent)."

---

## 3. Porosity — Review

### Reported values
- ~17% (GHGT-15, GHGT-13, IJGGG, AAPG), 14% mid-case (Gen-4), 16% (D-65 Well 8-19), 0.11–0.19 range (D-65 all wells), 16% (carbonstorage.io)

### FINDINGS

**GREEN: Strong convergence** — 17% for the high-quality injection interval is well supported, with 14–16% as conservative lower bounds.

**YELLOW: Two concerns:**

1. **The 17% consensus is not truly independent.** GHGT-15, GHGT-13, and IJGGG all have overlapping Shell authors. The AAPG presentation shares the same Rock author. Only the D-65 regulatory data and the Gen-4 report are distinct data sources. In practice, this is **2 independent data streams** (Shell model vs regulatory filings), not 5.

2. **14% (Gen-4 mid-case) vs 17% is a 3 p.u. (21% relative) spread** — notable for storage capacity calculations where each percentage point matters. The document explains this as full-zone vs injection-interval, which is fair, but the Gen-4 value being Shell's own model output should perhaps carry more weight than it's given.

3. **Measurement method not specified**: Are these log-derived porosities, core-plug measurements, or effective porosity from NMR? Each has different uncertainty (±1–2 p.u. for logs, ±0.5 p.u. for core).

---

## 4. Permeability — Review

### Reported values
- ~1000 mD (GHGT-15, IJGGG), darcy-range (GHGT-13), 150 mD (D-65 Well 8-19), 20–500 mD (D-65 injection wells), 150 mD (carbonstorage.io)

### FINDINGS

**RED: This is the weakest property in the dataset.**

1. **1000 mD vs 150 mD: 7× discrepancy.** The document attributes this to "formation-scale" vs "well-scale" measurements, which is a valid distinction, but:

2. **The 1000 mD value traces to a single source chain:** GHGT-15 (conference presentation) → cited by IJGGG, paraphrased by GHGT-13. There is no independent peer-reviewed measurement confirming 1000 mD at Quest. The **regulatory D-65 data shows 20–500 mD** across multiple wells. This is demonstrably the harder, more conservative data.

3. **The "formation-scale" justification would benefit from support**: In the Gen-4 report, what is the effective permeability in the reservoir model? If the model used 1000 mD to match historical data, that's supporting evidence. If it used lower values, that's contradictory.

4. **Recommendation**: The document should present the D-65 values (20–500 mD) as the primary data and flag the ~1000 mD as the optimistic end-member from a single conference source. Currently it does the reverse.

---

## 5. Area — Review

### Reported values
- ~3,667 km² / ~3,700 km² (AER, Lahvis)
- 666 km² (carbonstorage.io)

### FINDINGS

**GREEN: Well-defined and properly qualified.**

Strengths:
- SLA is clearly defined as the regulatory boundary.
- Calculation from townships to km² is shown stepwise.
- The conflicting carbonstorage.io value (666 km²) is flagged as unexplained and the AER value is correctly given precedence.

**YELLOW: Missing area concepts.**

The document does not clearly define:
- **SLA** (legal storage boundary, ~3,700 km²) — ✅ addressed
- **AOR / Area of Review** (10 km radial from injection wells, ~314 km² per well) — mentioned but not quantified
- **AOR** can also refer to the modeling area under the EPA/state UIC program definition, which may differ
- **CO₂ plume area** — NOT discussed. After ~7 Mt injection, the actual plume footprint is far smaller than the SLA (on the order of 10–100 km² depending on sweep). This is a critical omission for any capacity or containment analysis.
- **Gen-4 "Area of Interest"** (180 × 190 km = ~34,200 km²) — mentioned but never defined. This is vastly larger than the SLA and should not be confused with it.

**Recommended fix**: Add a row distinguishing SLA / AOR / plume area / model area, with approximate values for each.

---

## 6. Discrepancy Flagging — Review

| Discrepancy | Flagged? | Adequately? |
|-------------|----------|-------------|
| 17% (GHGT-15/GHGT-13/IJGGG) vs 14% (Gen-4) porosity | Yes | Yes — explained as interval vs full-zone |
| 1000 mD vs 20–500 mD permeability | Partially | **No** — characterized as scale difference but doesn't acknowledge that the 1000 mD source is weaker |
| 3,667 km² vs 666 km² area | Yes | Yes — correctly defers to AER |
| 40 m (GHGT-15) vs 40–80 m (AGS) thickness | Yes | Yes — explained as site-specific vs regional |
| 38 m net sand vs 40 m gross thickness | **No** | Not reconciled with N/G ratio explicitly |
| carbonstorage.io values vs all others | Yes | Yes — flagged but origin unexplained |
| 45.7 m thickness from carbonstorage.io — originally in feet, possible pre-metric report | **No** | Should note the unit conversion anomaly |

---

## 7. Cross-Validation Assessment

| Property | Truly independent sources | Verdict |
|----------|--------------------------|---------|
| Thickness | 3 (GHGT-15, Gen-4, AGS) | **Good** — Gen-4 and GHGT-15 are different Shell teams; AGS is fully independent |
| Porosity | 2 (Gen-4/Shell regulatory stream, D-65 regulatory stream) | **Adequate** — but the 5-way convergence is misleading due to shared authorship |
| Permeability | 1 (D-65 regulatory data) | **Weak** — the 1000 mD value is a single-source chain. Needs a third-party measurement or Gen-4 model confirmation. |
| Area | 2 (AER, Lahvis) | **Adequate** — AER is legally authoritative, Lahvis is independent journal confirmation |

---

## 8. Additional Issues

### 8.1 Missing Measurement Methodology
No source explicitly states whether porosity/permeability values are from:
- Core plugs (which probe ~1 in³)
- Well logs (vertical resolution ~0.3–2 ft)
- Well tests / pressure transient analysis (probe 10s–100s of meters)
- The distinction matters enormously for permeability (core: 150 mD; formation-scale: 1000 mD is plausible but unconfirmed at Quest)

### 8.2 Missing CO₂ Properties
For CCS storage calculations, CO₂ density at reservoir conditions (~21 MPa, 60 °C ≈ ~750 kg/m³) and CO₂ viscosity (~0.05 cP) are needed alongside the rock properties for any capacity or injectivity analysis.

### 8.3 Storage Complex vs BCS Only
The document notes "Storage complex thickness ~350 m (BCS through Upper Lotsberg)" but all reservoir property values are for the BCS only. The overlying Lotsberg Formation and other units are not characterized, yet they are part of the storage complex.

### 8.4 No Mention of Residual Trapping or Dissolution
For a CCS reservoir characterization, mention of residual gas saturation (Sgr) and dissolution trapping would be expected for CO₂ storage capacity calculations.

### 8.5 Kv/Kh Ratio
Anisotropy ratio (0.01–0.1) is important but rests on a single source (GHGT-15). No cross-validation.

---

## 9. Summary of Recommendations

1. **Permeability**: Reverse the weighting — present D-65 (20–500 mD) as primary, GHGT-15 (~1000 mD) as optimistic end-member with a caveat about source quality.
2. **Thickness**: State explicitly that the perforated interval is 30 m (not 40 m) and distinguish gross thickness from net pay.
3. **Porosity**: Add measurement method (log vs core vs NMR) for each source.
4. **Area**: Add a definitions table: SLA / AOR / plume / Gen-4 model area with values for each.
5. **Author independence**: Annotate which sources share authors to avoid overcounting cross-validation.
6. **carbonstorage.io**: Delete or relegate to a footnote — it's unverifiable and its values conflict without explanation.
7. **Net-to-gross**: Reconcile 38 m net sand (AAPG) with 40 m gross (GHGT-15) explicitly via the 0.9 N/G ratio.
8. **Residual trapping and dissolution**: Add Sgr and CO₂ density for completeness if this document feeds into capacity calculations.

---

## Verdict

| Criterion | Grade | Explanation |
|-----------|-------|-------------|
| Values reasonable | **B+** | All within expected ranges for BCS; some precision issues (net vs gross) |
| Well-supported | **B** | Good for porosity and thickness; weak for permeability |
| Source tiers clear | **C** | Not explicitly tiered; "peer-reviewed" label lumps journals and conferences |
| Cross-validation robust | **C+** | Appears stronger than it is due to overlapping authorship |
| Area clearly defined | **B** | SLA well-defined; plume area missing entirely |
| Discrepancies flagged | **B−** | Major area conflict flagged; permeability discrepancy under-flagged |

Document is **solid but needs refinement on permeability weighting, net-vs-gross thickness, and source independence disclosure**.
