# Adversarial Review: CANADA-CCS-WILLISTON-Midale-props.md

## 1. Reservoir Property Reasonableness

### 1.1 Thickness
- **Marly ~4 m net, Vuggy ~8 m net** — supported by Cenovus 2018 verbatim quote. Reasonable for a thin carbonate platform.
- **Total ~12 m net** — consistent with sum of component zones.
- **Evaporite seal ~5–20 m** — cited to "ScienceDirect" with no specific paper, URL, or entry in the Sources table. This is a **citation gap**. The evaporite thickness is not traceable.

### 1.2 Porosity
- **Marly 15–37% (avg 26%)** — Cenovus 2018; corroborated by Glemser (2008) at 16–38%. The two ranges are essentially identical, which is good.
- **However**: 26–38% porosity at ~1,450 m depth in a dolostone is at the **very high end** of what is physically plausible for a buried carbonate. Typical carbonate porosity at this depth is 10–25%. The file provides no explanation or caveat for why these values are so high (e.g., secondary dissolution, early dolomitization preserving porosity). This needs a comment.
- **Vuggy Shoal avg 10%; Vuggy Intershoal avg 10%** — both zones having the same average porosity despite being different facies is suspect. The Cenovus 2018 quote confirms this, but it may reflect a simplified presentation slide rather than measured reality. A caveat would help.
- **Cenovus internal inconsistency**: The line-39 verbatim quote gives "Marly 15–37% (26% avg); Vuggy Intershoal 2–15% (10% avg); Vuggy Shoal 5–20% (10% avg)". However, the line-136 verbatim quote (also Cenovus 2018) says "Porosity (%): Marly 25, Vuggy 12." The averages differ: **Marly 26% vs 25%; Vuggy 10% vs 12%**. These are two different slides/figures from the same presentation and should be reconciled or explained.

### 1.3 Permeability
- **Marly 1–100 mD, avg 6–10 mD** — supported by both Cenovus (6 mD stressed) and Glemser (1–150+ mD). The note on stress-dependent permeability (line 120) correctly explains the discrepancy between ambient and stressed measurements. Good.
- **Vuggy Shoal 1–500 mD, avg 20 mD** — 500 mD is **extremely high** for a microcrystalline limestone at 1,450 m. Even with vuggy porosity, permeabilities >100 mD at this depth are unusual. The file should flag this upper bound as potentially from a high-permeability streak or fracture, not representative of matrix permeability.
- **Vuggy (general) 0.3–500 mD, avg 6 mD stressed** — the 6 mD stressed average is consistent with the 20 mD ambient average for Shoal and 1 mD for Intershoal, assuming stress correction. The table structure is reasonable.

### 1.4 Other Parameters
- **Depth 1,450–1,500 m** — two sources agree. OK.
- **Temperature 63°C** — from Grokipedia only. Slightly high for this depth (typical gradient would give ~55–60°C). No cross-validation.
- **Pressure 12.5–18 MPa** — from Grokipedia only. At 1,450 m, hydrostatic is ~14.5 MPa, so the range covers slight sub- to over-pressure. Reasonable but lacks secondary source.
- **Incremental oil 2,300–5,800 bbl/d** — a 2.5× spread with no explanation of what drives the range. Is this seasonal? Operational? Source-specific? Unclear.
- **EOR efficiency ~3 bbl/tonne CO₂** — at the low end of the typical 3–7 bbl/tonne range. Reasonable for a mature flood.

## 2. Source Tiers

### Strengths
- Sources table correctly categorizes source types (industry, academic, government, journalism).
- Operator data (Cenovus, Cardinal, Whitecap) is appropriately weighted.
- Multiple operator sources cross-checked.

### Issues
- **Grokipedia is used as the sole source for temperature, pressure, and injection rates.** Grokipedia appears to be a wiki aggregator of unclear provenance. It is **not a primary or peer-reviewed source** and should not be the sole citation for critical reservoir parameters. These values need confirmation from primary literature or operator data.
- **"ScienceDirect" is cited for evaporite seal thickness** (line 25) but does not appear in the Sources table, and no specific paper/URL is provided. ScienceDirect is a publisher platform, not a source. This is a citation gap.
- **"AAPG" is cited for trap type** (line 80) with no specific paper title, author, or year. AAPG is a publisher, not a citation. Needs a specific reference.
- **LLNL report (2010)** appears in Sources but is never cited in any property table. If it is included for background only, that should be stated. Otherwise, data from it should be extracted.
- **Pipeline Online** (industry journalism) is the sole source for the critical CO₂ contract status update (contract end mid-2026). Acceptable for breaking news, but the document should flag this as unverifiable journalism.
- **OSTR NRCAN** forecast of 40+ Mt (line 88) cites a UUID URL but no specific document title. The URL is functional, but a proper citation would help trust.

## 3. Cross-Validation

### Well-validated
- Marly porosity: Cenovus 2018 vs Glemser 2008 — consistent.
- Marly permeability: Cenovus 2018 vs Glemser 2008 — consistent (accounting for stress vs ambient).
- Net pay: Cenovus 2018 quotes cross-reference with table values.

### Poorly validated (single-source data)
- **Temperature** — only Grokipedia.
- **Pressure** — only Grokipedia.
- **Midale injection rate** — MIT says 1,200 t/day; Grokipedia says 1,500 t/day. This is a **25% discrepancy** with no reconciliation. The table averages them (1,200–1,500 t/day) but does not acknowledge the conflict.
- **Total CO₂ stored (44 Mt)** — only Whitecap's corporate statement. Reasonable for an operator figure, but no independent audit or peer review cited.
- **Evaporite thickness** — unverifiable (ScienceDirect citation gap).

## 4. Area Definition

### Critical Gap: Midale Unit Area Not Stated
- The table gives **Weyburn Unit area** (180–210 km² from Wikipedia; 22,000 ha / 85 sq mi from Cenovus 2018) but **nowhere states the Midale Unit's area**.
- Midale OOIP is 730 MMbbls (~52% of Weyburn's 1.4 billion bbls), suggesting the Midale Unit is roughly half the size of Weyburn or has a different net-pay/porosity distribution. Without area, this cannot be evaluated.
- For any storage resource calculation, area is a critical variable. Its absence is a **blocking gap** for quantitative use of this document.

### Additional Note
- The two Weyburn area values (210 km² from Wikipedia vs 220 km² / 22,000 ha from Cenovus) are internally consistent within ~5%, which is acceptable.

## 5. Midale-Specific vs Weyburn-Borrowed Data

### Handling
- The Notes on Data Quality section (lines 110–122) addresses this issue directly and honestly. This is a **strength** of the document.
- The Sources table explicitly labels which sources cover Midale vs Weyburn vs both.

### Problems
- **In practice, nearly all detailed reservoir characterization data comes from the Weyburn Unit** (Cenovus 2018). The Cenovus presentation covers the Weyburn Unit's Midale Beds, but these are assumed (not demonstrated) to apply to the Midale field 10 km away. While the formation is regionally continuous, local heterogeneity in porosity, permeability, and thickness is expected.
- **Cardinal Energy's website provides OOIP and production data for Midale but zero reservoir properties** (no porosity, permeability, or net pay). If Cardinal has internal petrophysical data, it is not reflected here.
- **The Cenovus 2018 line-136 quote says "Midale reservoir — Average Net Pay (m): Marly 4, Vuggy 8"** within a presentation titled "Weyburn Unit Technical Presentation." It is ambiguous whether "Midale reservoir" refers to the Midale *formation* within the Weyburn Unit or the Midale *field*. This should be clarified — it changes whether these values are directly applicable to the Midale Unit or only to Weyburn.
- **Storage figures are explicitly combined** (noted in line 116), and Midale's share (~15–20%) is estimated. This is transparent and appropriate.

## Summary of Findings

| Severity | Issue |
|----------|-------|
| **Critical** | Midale Unit area is not defined anywhere in the document |
| **Critical** | Evaporite seal thickness cited to "ScienceDirect" with no specific reference |
| **Major** | Grokipedia is sole source for temperature, pressure, and injection rates |
| **Major** | Nearly all reservoir property data borrowed from Weyburn; applicability to Midale assumed, not tested |
| **Major** | Cenovus 2018 has internal inconsistency on porosity averages (Marly 26% vs 25%; Vuggy 10% vs 12%) |
| **Minor** | Marly porosity range (16–38%) at 1,450 m depth is unusually high — needs caveat |
| **Minor** | Vuggy permeability upper bound of 500 mD is extreme for matrix limestone at depth |
| **Minor** | Vuggy Shoal and Intershoal have identical average porosity (10%) — physically suspect |
| **Minor** | Midale injection rate has 25% discrepancy between MIT (1,200 t/day) and Grokipedia (1,500 t/day) with no reconciliation |
| **Minor** | Incremental oil range (2,300–5,800 bbl/d) is very wide with no explanation |
| **Minor** | AAPG cited for trap type with no specific paper reference |
| **Minor** | LLNL (2010) in source list but never used in property tables |

## Recommendations

1. **Define the Midale Unit area.** This is the single most important gap. Without it, storage resource estimates are not reproducible.
2. **Replace or supplement Grokipedia citations** with primary sources (operator publications, IEA GHG reports, peer-reviewed papers).
3. **Resolve the ScienceDirect citation** by finding the specific paper and URL for evaporite thickness.
4. **Resolve the Cenovus 2018 internal inconsistency** between the two sets of porosity averages from the same presentation.
5. **Add a caveat on Marly porosity** explaining why values exceed typical carbonate porosity at this depth.
6. **Reconcile the Midale injection rate discrepancy** between MIT (1,200 t/day) and Grokipedia (1,500 t/day).
7. **Add specific AAPG reference** for the trap type statement.
8. **Clarify whether Cenovus 2018 "Midale reservoir" applies to the Midale field or only the Weyburn Unit.**
9. **Reduce the incremental oil range** or explain the source of variability.
