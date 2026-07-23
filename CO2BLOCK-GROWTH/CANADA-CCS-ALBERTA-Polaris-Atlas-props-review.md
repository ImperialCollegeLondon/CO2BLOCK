# Adversarial Review: Shell Polaris + Atlas Hub CCS — Reservoir Properties

## Summary

The file is a competent compilation of publicly available data with honest caveats about the lack of Polaris-specific subsurface information. However, the adversarial review identifies **one critical error** (Edmonton-region permeability is overestimated by ~5-10x), **three significant gaps** (no storage efficiency analysis, no geomechanical context, no discussion of Quest's injectivity issues as a risk analog), and **several source-quality nuances** that weaken the confidence of database values. The overall framework and sourcing are sound; the numerical inputs need adjustment for the Edmonton-area context.

---

## Property-by-Property Review

### 1. Reservoir Thickness

| Claim | Assessment |
|-------|------------|
| Regional 17-115 m | **Well-supported.** AGS Atlas WCSB Ch. 8 is the authoritative source. |
| CCS target 40-80 m | **Well-supported.** AGS PRS 2024-001, 28 cores. |
| Quest ~45 m | **Well-supported.** Confirmed by Halite Damage Study (SciDirect 2022) and Quest D65 Application. |
| Polaris "not disclosed" | **Correct.** Honest caveat. |

**Missing:** No estimate for Polaris from structural mapping of the BCS in the Alberta Industrial Heartland area (Twp 52-55, Rge 21-24 W4M). Regional isopach maps exist (AGS, Geoscience BC) that could constrain the Polaris lease estimate to ~40-50 m. Not a file error — a missed opportunity.

### 2. Porosity

| Claim | Assessment |
|-------|------------|
| Regional avg ~10% (Weides 2014) | **Correct but sensitive.** Based on 6 wells only. The same study notes a "zone of high porosity and permeability" in the northern well — omitting this single outlier would shift the average notably. |
| Quest ~17% | **Well-supported.** Halite Damage Study, confirmed. Quest area is anomalously good. |
| Range 8-24% | **Well-supported.** Quest D65 Application. |
| Geothermal 14.8% | **Partially misleading.** See Critical Issues #2. |
| Facies-dependent 7-15% | **Adequately supported.** Consistent with AGS and Menger 2024 descriptions. |
| 4% minimum viability | **Unsupported.** No source citation for this exact threshold. It may be a generalized CCS screening rule of thumb, not a Polaris-specific criterion. |

### 3. Permeability — THE CRITICAL ISSUE

| Claim | Assessment |
|-------|------------|
| Regional avg ~10 mD (Weides 2014) | **Partially correct but incomplete.** Weides et al. tested 6 wells across central Alberta. Importantly, Hofmann et al. (2013, Stanford Geothermal Workshop) separately evaluated wells *around Edmonton* specifically and found **average regional permeability of 1.6 mD**. The Edmonton-area dataset (~1.6 mD) is more relevant to Polaris (~45 km east of Edmonton) than the broader Weides average (~10 mD, which includes a high-permeability well in the northern study area). |
| Quest ~1,000 mD | **Correct for Quest.** Confirmed by Halite Damage Study. |
| Range 1 mD to >1 D | **Well-supported.** Quest D65 Application. |
| Geothermal: 1.9 mD horiz, 0.37 mD vert | **Mischaracterized in context.** See Critical Issues #2. |
| Regional aquifer range <10 mD to several thousand mD | **Well-supported.** Geoscience BC. |

**Bottom line on permeability:** The Polaris/Atlas lease (~45 km east of Edmonton) likely has permeability closer to **1-2 mD** (Edmonton-region average from Hofmann 2013) than to the ~10 mD cited as "regional average" from Weides 2014, which covers a larger area including a high-perm northern zone. The Quest ~1,000 mD should be treated as an upper-end outlier, not representative. The file is **500-1,000% too optimistic** on permeability for the Polaris location.

### 4. Area / Lease Data

| Claim | Assessment |
|-------|------------|
| Atlas ~45 km E of Edmonton | **Correct.** Natural Gas World; Alberta Government. |
| 22 km pipeline, 2 injection wells | **Correct.** Shell FID announcement. |
| First Alberta CSA (July 2024) | **Correct.** Alberta Government / AER. |
| ~250,000 km² BCS extent | **Correct.** Canadian Discovery Ltd. SPE SRMS assessment. |

**Missing:** The actual lease area (townships, surface area in km²). The CSA area is publicly available from Alberta Energy. Without this, storage capacity per km² cannot be cross-checked against the 300 Mt claim.

### 5. Other Key Parameters

| Claim | Assessment |
|-------|------------|
| Target: BCS | **Correct.** |
| Depth ~1,800-2,100 m | **Quest-derived analog.** Reasonable for Polaris. |
| Pressure ~20 MPa | **Quest-derived.** Hydrostatic at 2 km = ~20 MPa is reasonable. |
| Temp 65-120°C | **Correct as regional range.** But for Polaris (~2 km, gradient 35.6°C/km): expected ~76°C. File could state a site-specific estimate. |
| Triple seal system | **Quest-derived.** Reasonable analog. The Lower Lotsberg Salt is regionally extensive. |
| 300 Mt lifetime | **Severely undersupported.** See Critical Issues #3. |
| Phase 1: 650 kt/yr, 2 wells | **Correct.** Shell FID. |
| Phase 2: 7-10 Mtpa | **Speculative.** Cited to 2021 proposal; no FID. File does note this is contingent. |

---

## Source Quality

### Tier Classification

| Tier | Sources | Assessment |
|------|---------|------------|
| **Tier 1 — Peer-Reviewed** | 10 (Weides 2014, CJES), 16 (Halite Damage, 2022), 17 (Injectivity, 2020), 19 (PMC EGS) | Solid. All extant and correctly linked. |
| **Tier 2 — Gov't/Regulatory** | 9 (AGS PRS 2024-001), 13 (AGS Atlas), 14 (Quest D65), 15 (Quest Gen-4), 8 (CER), 18 (CER) | Excellent. D65 and Gen-4 are primary regulatory filings — highest reliability for technical data. |
| **Tier 3 — Industry News** | 1 (Shell FID), 2, 3 (Alberta Major Projects), 4 (CBC), 5 (NatGasWorld), 6, 7, 20, 21 | Adequate for announcements/project status. Not reliable for technical reservoir data. |

### Specific Source Issues

- **Source 19 (PMC EGS):** The file describes this as an "Enhanced Geothermal in BCS" study and lists its values in the main property tables. However, this is a **modeling study** — the permeability and porosity values are taken from literature (Chong et al., which draws from Weides et al.) and applied uniformly across a 1 km × 1 km model. They are **not** from a specific well or site. The file does not disclose this provenance in the main tables — only the source URL is given. The values (1.9 mD, 14.8%, 0.37 mD) are reasonable Edmonton-region averages but the source needs annotation: "parameters from literature, not site-specific measurements."

- **Source 10 (Weides 2014):** The 6-well sample size is small. The file should note this explicitly. The "zone of high porosity and permeability" noted in the northern well may correspond to the Quest-area trend, but the paper does not identify it as such.

- **Source 14 (Quest D65):** This is a Shell application — the data is for the Quest area, not Polaris. Correctly used as reference.

### Cross-Validation Summary

| Parameter | Cross-Validation | Verdict |
|-----------|-----------------|---------|
| Thickness | Multiple sources converge 40-80 m | **Good** |
| Porosity | Wide range, location-dependent, internally consistent with heterogeneity narrative | **Adequate** |
| Permeability | Weides ~10 mD vs Hofmann ~1.6 mD (Edmonton) — CONFLICT not resolved or even mentioned | **INADEQUATE — the single biggest weakness** |
| Pressure | Hydrostatic from depth, multiple sources agree | **Good** |
| Temp | Gradient from ~2,000 values, consistent | **Good** |
| 300 Mt | No independent calculation provided | **INADEQUATE** |

---

## Critical Issues

### CRITICAL 1: Edmonton-Region Permeability is Overestimated

**Severity: HIGH**

The file cites Weides et al. 2014 for a "regional average" of ~10 mD. However:

- Hofmann et al. (2013) — a companion study using the same core data plus additional wells from the Edmonton area — found **average regional permeability of 1.6 mD** for wells around Edmonton.
- The Polaris/Atlas lease is ~45 km east of Edmonton — squarely in the area covered by Hofmann's Edmonton-region analysis.
- Weides et al. explicitly note a "zone of high porosity and permeability" in the northern part of their study area, which inflates their average. Removing that single well would bring their average much closer to Hofmann's 1.6 mD.

This means the true permeability at Polaris is likely **50-100% lower** than the ~10 mD stated in the file. For injectivity (2 wells, 650 ktpa), this is a material difference. Quest succeeds because it has ~1,000 mD — ~500-600x higher than what Polaris can likely expect.

**Recommendation:** Add Hofmann et al. 2013 as a source and revise the Edmonton-area permeability estimate to **1-2 mD** range, with the ~10 mD from Weides presented as the central Alberta basin-wide average (which includes the high-perm northern zone).

### CRITICAL 2: Source 19 (PMC EGS) is a Modeling Study, Not a Field Study

**Severity: MEDIUM-HIGH**

The file presents Source 19 values (1.9 mD, 14.8% porosity, 0.37 mD vertical) as if from a site-specific geothermal investigation. The paper (PMC article 10835245, Zarei et al. 2024) is a **numerical simulation study** that uses properties from literature. The authors state: "Chong et al. used constant properties with 1.9 mD of horizontal permeability... Here, in the BCSU, the horizontal and vertical permeabilities and porosity are uniformly randomly distributed with same mean values as that of Chong et al." The location is a generic 1 km × 1 km model — not a specific well location.

While the values coincidentally align well with Hofmann's Edmonton-region average (1.6 mD), the source provenance needs correction.

### CRITICAL 3: 300 Mt Lifetime Storage is Unsubstantiated

**Severity: MEDIUM-HIGH**

The ~300 Mt claim cites Alberta Major Projects (Source 2), which is a government tracking website, not a technical assessment. The file should attempt a first-order verification:

- BCS thickness at Polaris: ~45 m (estimated)
- Porosity: ~10% (conservative estimate)
- CO2 density at ~20 MPa, ~75°C: ~700 kg/m³
- Storage efficiency factor (saline aquifer, no hydrodynamic trapping credit): typically 2-4%
- Area required for 300 Mt: ~350-700 km² plume area

This requires clarification. Is 300 Mt the **P50 resource**, the **theoretical capacity**, or the **commercial storage contract volume**? The file notes no storage efficiency, no areal extent of the lease, and no pore volume calculation.

**Recommendation:** Add a back-of-envelope capacity calculation with stated assumptions. Distinguish between theoretical, contingent, and commercial capacity using SPE SRMS or similar framework.

### CRITICAL 4: No Discussion of Quest Injectivity Degradation as Risk Analog

**Severity: MEDIUM**

The Halite Damage Study (Source 16) — which the file does cite for porosity/permeability numbers — is actually about **halite-induced injectivity damage at Quest**. Key findings the file does not mention:

- Quest's injectivity declined significantly due to halite precipitation, requiring fresh-water washes
- Skin values reached 147; damage zone permeability reduced to ~10 mD (a 99% reduction from native ~1,000 mD)
- If Polaris has native permeability of 1-2 mD, halite damage could reduce permeability to near-zero levels

This is critical context for any database entry about the Polaris project.

---

## Recommendations

### Immediate Corrections (Database Impact)

1. **Revise Polaris permeability** to **1-2 mD** (Edmonton-region average from Hofmann 2013). Keep the Weides 2014 ~10 mD as the **central Alberta basin-wide average** with annotation that it includes a high-perm northern zone.

2. **Annotate Source 19** as a modeling study using literature values, not in-situ measurements.

3. **Add Hofmann et al. 2013 (Stanford Geothermal Workshop)** as Source 22 for the Edmonton-region permeability data. URL: https://pangea.stanford.edu/ERE/pdf/IGAstandard/SGW/2013/Hofmann.pdf

### Strengthen (Add New Content)

4. **Add storage efficiency calculation** for the 300 Mt claim. Show the implied pore volume and efficiency factor. Distinguish between P10/P50/P90.

5. **Add Quest injectivity degradation** as an operational risk factor (halite damage reference from Source 16's main finding, not just its property values).

6. **Add site-specific temperature estimate** for Polaris: ~76°C at 2 km using 35.6°C/km gradient.

7. **Add fracture gradient / maximum injection pressure** estimate. The triple seal is described but without maximum allowable injection pressure (typically 70-80% of fracture gradient). Without this, injectivity modeling is incomplete.

### Refine (Improve Existing Content)

8. **Break the "Other Parameters" table** into logical groups: (a) Reservoir geometry, (b) Fluid/thermal, (c) Seals, (d) Operations, (e) Capacity.

9. **Add areal extent of the Atlas lease** in km² or townships (publicly available from Alberta Energy's CSA database).

10. **Add CO2 density** at reservoir conditions (~700 kg/m³ at 20 MPa, 75°C) for completeness.

---

## Minor Issues

1. **Source 12 URL (Geoscience BC / CDL):** `cdl.canadiandiscovery.com` — this is a commercial subscription site. The SPE SRMS assessment may not be publicly accessible. Verify the URL still resolves or note access restrictions.

2. **Source 7 (2021 proposal):** Pipeline distance changed from 12 km (proposal) to 22 km (FID). The file notes this in the data quality notes but doesn't flag that the Phase 2 capacity (10 Mtpa) in the same source row is from the same outdated proposal.

3. **Source 18 (CER 2022):** FID was June 2024 — the 2022 snapshot is pre-FID. Correctly used but aging.

4. **Source 15 (Quest Gen-4):** The URL is a permanent AER Open Alberta link but includes a `resource/` path that may change on file system reorganization. Add the permanent Open Alberta record ID (612332b1-475c-48c2-a071-3784f428dc9c) for robustness.

5. **Minor formatting issue:** Inline code block on line 42 (`Source [9]`) and others use `>` Markdown blockquotes, which is fine. However, values in tables are inconsistently bolded — some key numbers are bolded (`**~10%**`) while others are not. Either bold all key values or none, for consistency.

6. **Line 66** "4% minimum for viable CCS" — no source is cited for this specific threshold in the table (the Geoscience BC source row is generic). Either add a specific citation or flag it as a rule of thumb.

7. **Lines 93-94** (regional aquifer range) — the table rows for regional range and Source 12 Geoscience BC assessment are redundant with the row above (Range 1 mD to >1 D). Consider consolidating.

8. **Line 139** "Clean quartz arenite (>90% quartz)" — is well-supported by Weides 2014 petrography and AGS. No issue.

---

## Conclusion

The file is well-structured and honest about data limitations. The single most impactful correction is the **permeability estimate for the Polaris location**: Hofmann et al. (2013) indicates ~1-2 mD for the Edmonton region rather than the ~10 mD from Weides et al. (which covers a larger, heterogeneous area). This 5-10x overestimate has downstream implications for injectivity assessment and the viability of the 2-well, 650 ktpa injection plan without stimulation.

The 300 Mt capacity claim and the characterization of Source 19 (modeling study as field data) are secondary but important issues requiring revision.

After these corrections, the file will be a robust, defensible compilation for the CO2BLOCK database.
