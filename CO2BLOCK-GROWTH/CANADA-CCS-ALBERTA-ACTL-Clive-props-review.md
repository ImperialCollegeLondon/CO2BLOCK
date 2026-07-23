# Adversarial Review: ACTL/Clive CCS — Reservoir Properties

**Reviewer**: adversarial-reviewer
**Date**: 2026-07-23
**Document**: CANADA-CCS-ALBERTA-ACTL-Clive-props.md

## Summary

**Assessment: CONDITIONAL — FAIL in current state without major revisions.**

The document compiles a useful set of sources but exhibits a critical flaw: **all quantitative property values (porosity, permeability, net-pay thickness) trace to a single data chain — Enhance Energy's internal geo-model, cited indirectly through a 2025 peer-reviewed chapter (S7) or through ACTL operator reports (S1–S5).** There is no independent, third-party, publicly auditable measurement data. The area definition is ambiguous across four different spatial concepts that are conflated in the recommended values. Key engineering parameters needed for defensible CO2 storage resource estimation (water saturation, in-situ CO2 density, relative permeability) are absent.

---

## Property-by-Property Review

### Thickness

**Value**: 8.75 m (net pay), 35 m (reservoir interval), 250 m (total Leduc Fm)

**Assessment**: Reasonable — but the data chain is narrow.

**Issues**:
1. **Single operator source**: The 8.75 m net pay value is corroborated only by S2 and S3, both ACTL Summary Reports representing the same operator's own internal calculations. These are not independent sources; S3 (2017) and S2 (2018) likely derive from the same geo-model.
2. **8.75 m vs 35 m ambiguity**: The document recommends 8.75 m for injection and 35 m for "storage volumetric calculations," but this conflates two different physical concepts. The 8.75 m is the **perforated/net-pay interval** — the artificially completed zone open to injection. The 35 m is the **cored reservoir interval** with better-quality rock. A CO2 plume migrating buoyantly will access whatever permeable rock it encounters, which is somewhere between these values. The document should clarify which thickness applies to which calculation type (well injectivity vs storage resource).
3. **No core-based net-pay cutoff**: The 8.75 m figure is reported without the porosity or permeability cutoffs used to define it. Without knowing whether net pay was defined at φ > 6 %, k > 1 mD, or some other threshold, the number is not auditable.

**Cross-validation**: Weak. S2 and S3 are not independent. S6 provides the vertical subdivision (Massive Sandy Zone, Layered Muddy Zone, cemented baffle) but does not independently compute net pay thickness.

**Recommendation**: State the net-pay cutoff criteria (porosity and permeability thresholds used). Confirm whether the 8.75 m is arithmetic mean or thickness-weighted mean. Clarify that 8.75 m is the completion interval; the storage resource interval may differ.

---

### Porosity

**Value**: 6.1 % (average, Enhance geo-model), with range 3–9 % (locally up to 12 %)

**Assessment**: Problematic — single-source dependency and unsupported upper bound.

**Issues**:
1. **Single quantitative source**: The 6.1 % value comes exclusively from S7 (Tavallali et al. 2025, Wiley), which cites Enhance Energy's internal model [23]. The data chain is: Enhance internal model → cited in peer-reviewed chapter → reported here. There are **zero independently measured core-plug porosity values** reported or cited.
2. **No arithmetic mean calculation shown**: The 6.1 % is stated as the "average porosity from Enhance geo-model," but the type of average (arithmetic, thickness-weighted, volume-weighted) is not specified. For heterogeneous carbonate reservoirs, these can differ significantly.
3. **Unsupported 12 % upper bound**: The stated range of "locally up to 12 % in vuggy intervals" has no supporting citation. The Cooking Lake maximum is 9 % (S6), but no source supports the 12 % figure. This appears speculative.
4. **Crossplots undigitized**: S1–S5 contain porosity–permeability crossplots as images only. These cannot be used for independent quantitative verification.

**Cross-validation**: **FAIL** — only one independent quantitative source (S7, itself derivative of Enhance's model). The 4 % Bashaw Platform average (S6) is for a different scale (platform aquifer, not reservoir zone) and provides context but not validation.

**Recommendation**: Identify the average type (arithmetic vs volume-weighted). Source the 12 % upper bound or remove it. Add core-plug porosity data if available in any ACTL report appendices. Flag the single-source dependency prominently.

---

### Permeability

**Value**: 30 mD mode; range 1–100 mD; Cooking Lake avg 0.24 mD

**Assessment**: CRITICAL ISSUE — mode is not appropriate for storage resource estimation.

**Issues**:
1. **Only the mode is reported, not the mean**: This is the most significant technical issue in the document. For CO2 storage resource estimation, the appropriate permeability statistic depends on the calculation:
   - **Injectivity calculations** → need arithmetic mean or kh (flow capacity)
   - **Plume migration modeling** → need full distribution or geometric mean
   - **Permeability-thickness (kh) product** → arithmetic mean × thickness
   The mode (30 mD) is the *most frequent value* — it says nothing about the flow capacity of the system. If the permeability distribution is log-normal (typical for carbonates), the arithmetic mean may be 2–10× higher than the mode, and the geometric mean lower.
2. **Same single-source dependency as porosity**: The 30 mD mode also comes from S7 citing Enhance's internal model. No DST, well test, pressure transient analysis, or core plug permeability data is cited.
3. **No kv/kh data**: The document qualitatively notes anisotropy in the Layered Muddy Zone but provides no quantitative kv/kh ratio. For CO2 plume migration, vertical permeability is critical for assessing buoyant flow.
4. **"~3× water injectivity" claim**: The qualitative statement about CO2 injectivity (S2, S5) is not quantitatively supported by any published well test analysis.

**Cross-validation**: **FAIL** — no independent quantitative value exists. Only S7 provides a number, and it's only the mode.

**Recommendation**: 
- Add (or search for) published kh, DST, or well-test permeability data from any ACTL report or AER well-file data.
- State explicitly that no arithmetic or geometric mean is available and that the mode is insufficient for most engineering calculations.
- If mode is the only available statistic, note that storage resource calculations should use a log-normal distribution conservatively centered on 30 mD.
- Add kv/kh data or note it as a data gap.

---

### Area

**Value(s)**: 544,512 ha (lease), ~16,000–20,000 ha (Clive Unit inferred), ~40 km² (pinnacle reef pool estimate)

**Assessment**: **CRITICAL ISSUE** — conflates four fundamentally different spatial concepts.

**Issues**:
1. **SLA vs AOR vs plume vs pore space not distinguished**: The document presents three numerical values without defining which are:
   - **SLA (Storage Lease Area)**: 544,512 ha — formal Crown lease for the Origins hub
   - **AOR (Area of Review)**: Not defined or estimated
   - **CO2 Plume footprint**: Not defined or estimated
   - **Pore-space tenure relevant to Clive**: Unclear
2. **544,512 ha is the entire Origins Hub lease, not Clive-specific**: Using 544,512 ha in any volumetric calculation for the Clive Leduc D-3A pool would massively overstate storage capacity by 2–3 orders of magnitude (the Leduc D-3A pool likely covers ~40 km²; the lease area is ~5,445 km²).
3. **~40 km² pool area is unsourced**: The "Baseline Leduc D-3A pool area ~40 km²" is described as "order-of-magnitude for a single pinnacle-reef pool" but without a specific citation for the Clive field. Not every Leduc reef is the same size.
4. **Clive Unit area (16,000–20,000 ha) is inferred from well count**, not from a published source. Inferring area from 168 wellbores assumes uniform well spacing, which is not stated.

**Cross-validation**: **FAIL** — each area value serves a different definition. No single value is cross-validated by an independent source.

**Recommendation**:
- **Must define SLA, AOR, and plume footprint separately.**
- Add a dedicated section clarifying which area applies to which calculation:
  - P10/P50/P90 plume footprint for CO2 storage resource
  - SLA for regulatory/tenure context
  - Pool area for hydrocarbon-column-based volumetric comparisons
- Cite the specific pool area from an AER pool database or published geological map of the Clive Leduc D-3A pool.
- Remove or clearly caveat the 16,000–20,000 ha "inferred" value.

---

## Source Quality

### Source Tier Analysis

| Tier | Definition | Sources in this doc |
|------|-----------|---------------------|
| **Tier 1** | Direct measurements (core, logs, DST, well tests) | **None** |
| **Tier 2** | Published operator interpretations (geo-model, petrophysics) | S1, S2, S3, S4, S5, S6 (all ACTL, all same operator) |
| **Tier 3** | Peer-reviewed synthesis citing operator data | S7 (Wiley 2025, Tavallali et al.) |
| **Tier 4** | Grey literature, industry reports, news | S8, S9, S10, S11 |

### Key Finding
The document has **no Tier 1 sources**. The quantitative chain is:
Enhance internal geo-model (Tier 2) → cited in peer-reviewed chapter (Tier 3, S7) → cited in this document.

This is an inherently weak chain for a CO2 storage database intended for resource estimation. Any error or bias in Enhance's model propagates directly into the reported values without possibility of external audit.

### Missing Sources
- **AER well data**: Individual well logs, core analyses, DST results from the 168 Clive wells are publicly accessible through the AER data portal. None are cited.
- **Published petrophysical evaluations**: No independent petrophysical study of the Leduc D-3A at Clive is referenced.
- **PCOR Partnership publications**: The PCOR Partnership has published regional characterization work on the Bashaw Platform that may contain cross-validation data.
- **NRCan / GSC storage assessments**: Geological Survey of Canada carbon storage assessments for the Alberta Basin may have independent property estimates for Leduc reefs.

---

## Critical Issues

1. **AREA AMBIGUITY (Must Fix)**: The document cannot recommend 544,512 ha for storage work without defining SLA vs AOR vs plume vs pore-space footprint. Including four different area definitions without clear mapping to calculation types will produce erroneous storage resource estimates. This is the single most consequential ambiguity.

2. **PERMEABILITY MODE vs MEAN (Must Fix)**: Using mode permeability (30 mD) for any engineering calculation is inappropriate. The document must either (a) find and report the arithmetic mean from the geo-model, or (b) add a prominent caveat that mode is insufficient and a log-normal distribution assumption is needed, or (c) find DST/well-test kh values from ACTL reports.

3. **NO INDEPENDENT CROSS-VALIDATION (Must Flag)**: Every quantitative property (6.1 % φ, 30 mD k, 8.75 m net pay) traces back to Enhance Energy's internal geo-model. There is no independent audit trail. The database should clearly label this as "operator-cited values, not independently verified."

4. **WATER SATURATION NOT REPORTED (Must Fix)**: Sw is absent from the summary table (line 142–157). A CO2 storage volumetric calculation requires (1 − Sw) — without it the table is incomplete. This is a mandatory input for any storage resource equation.

5. **CO2 DENSITY NOT CALCULATED**: Pressure (~13,000 kPag) and temperature (~69 °C) are given but CO2 density at these conditions is not. This is trivially computable (~650–750 kg/m³ at these conditions) and should be stated explicitly for the user.

---

## Recommendations

### Must Fix (before revise stage)
1. **Define area concepts**: Add separate rows for SLA (544,512 ha), pool area (cite AER pool map), and estimated plume footprint. Remove the 16,000–20,000 ha inferred value or caveat it heavily.
2. **Flag the absence of arithmetic mean permeability** and add guidance for log-normal distribution assumption.
3. **Report water saturation** (Sw) from available sources — the 18.5 m oil column and 11 m gas cap imply initial saturations that should be documented.
4. **State net-pay cutoff criteria** for the 8.75 m thickness.
5. **Calculate and report CO2 density** at in-situ P/T conditions.

### Should Fix
6. **Source the 12 % vuggy porosity upper bound** or remove it.
7. **Note that S2 and S3 are not independent** — both are ACTL Summary Reports likely based on the same model run.
8. **Add a source quality flag** to each recommended value (e.g., "operator-cited, not independently verified").

### Nice to Have
9. **Digitize the ACTL Report crossplot images** (Figures 3.4.5–3.4.10 in S1) to extract the φ–k dataset for independent statistical analysis.
10. **Search AER well data** for any publicly available core analysis or DST from Clive wells.
11. **Cross-reference with PCOR Partnership Atlas** values for Bashaw Platform Leduc properties.
12. **Add kv/kh ratio** or note as data gap.
13. **Specify the averaging method** for 6.1 % porosity (arithmetic vs volume-weighted vs thickness-weighted).

---

## Minor Issues

- **S9 (CarbonStorage.io) returned 404**: The URL `https://carbonstorage.io/projects/origins` is not currently accessible. Verify the source is still online or replace with a cached/archived copy.
- **S10 (Key Facts Energy)**: This is a paywalled industry news site. The specific article about Chinook Petroleum / Origins lease area should be cited with a persistent link or archived copy (e.g., Wayback Machine).
- **S11 (Geoconvention 2025, Menger)**: Cited for P/T values but the presentation abstract is not linked. Add the DOI or direct URL if available.
- **Oil column vs CO2 column**: The document references the hydrocarbon column (11 m gas + 18.5 m oil = ~30 m) as context for thickness, but CO2 column height in a storage scenario will differ due to different capillary entry pressures and relative permeability. This should be noted.
- **No mention of the D-2 (Nisku) pool**: The summary table (lines 142–157) only lists Leduc D-3A. If the Nisku is a secondary injection target (6.4 m avg thick, S2), its properties should have a parallel table or at least a cross-reference.
- **"~3× water injectivity" for CO2**: This ratio depends on relative permeability, viscosity ratio, and saturation history. Without supporting well-test data, this is an unsupported claim. Either cite the source report page number or remove it.

---

*Review generated using document text and web-verified source status. No access was available to the ACTL Detailed Report PDF (403 on direct fetch) or to the CarbonStorage.io Origins page (404). DOI 10.1002/9781394356294.ch5 confirmed as a valid Wiley publication by Tavallali et al. (2025).*
