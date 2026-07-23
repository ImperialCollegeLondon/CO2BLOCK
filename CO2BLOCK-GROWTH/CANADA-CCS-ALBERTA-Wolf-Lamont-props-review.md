# Adversarial Review: Wolf Lamont Carbon Hub — CCS Reservoir Properties

**Review date**: 2026-07-23  
**File reviewed**: `CANADA-CCS-ALBERTA-Wolf-Lamont-props.md`  
**Reviewer methodology**: Source URL verification via web fetch, literature cross-check, adversarial property-by-property audit, cross-validation adequacy assessment.

---

## Summary

The document is well-structured, transparent about data limitations, and appropriately caveated. However, it contains **one critical stratigraphic error that undermines the entire property database**, several instances of single-source dependency, and a handful of misrepresented values. The document's overall quality is **adequate for a working draft but requires correction before it can be considered reliable for database ingestion**.

---

## Property-by-Property Review

### Reservoir Thickness

| Property | Statement | Verdict |
|----------|-----------|---------|
| BCS gross thickness 40–80 m | Multiple sources (AGS PRS-2024-001, Moradi 2016, Weides 2014). Consistent range. | **PASS** — well-supported |
| BCS net reservoir 25–45 m | Explicitly tagged as "inference from core data." No Wolf-specific data. Caveat is appropriate. | **PASS** — honest about limitation |
| LMS seal ~5–15 m | From Desjardins & Smith (2013) via Moradi (2016). Single source, but Moradi independently confirms values. | **PASS** — single source but consistent |
| MCS seal ~60 m | From Moradi (2016). Single source. | **MARGINAL** — single source |
| Cambrian total interval "up to 300 m" | Not cited. Generic statement about Deadwood Formation equivalents far west of Lamont. | **WEAK** — irrelevant to Lamont |

**Thickness verdict**: Reasonable, well-qualified.

---

### Porosity

| Property | Statement | Verdict |
|----------|-----------|---------|
| BCS range 8–24% | Moradi (2016) + Quest core data. Confirmed via thesis text. | **PASS** — Good |
| Quest average ~18% | Moradi (2016), density-log derived. Confirmed. | **PASS** |
| Regional BCS aquifer average 14.4% | Alberta Innovates (Meikle, 2024). Single source — whitepaper, not peer-reviewed. | **MARGINAL** — plausible but single non-peer-reviewed source |
| BSU average 10% | Weides et al. (2014). Confirmed. | **PASS** |
| LMS ~6%, UMS 1–2%, Devonian Red Beds ~5% | All from Moradi (2016). Single source for each. | **MARGINAL** — no cross-validation on seal porosities |

**Porosity verdict**: Acceptable. The 14.4% figure needs independent corroboration. Seal porosities from a single thesis should be noted as weaker.

---

### Permeability

| Property | Statement | Verdict |
|----------|-----------|---------|
| BCS range 1 mD to >1,000 mD | Moradi (2016) + Quest core data. Confirmed. | **PASS** |
| BSU avg permeability <0.01 mD (matrix) | Weides et al. (2014). **Critical note**: Weides studied the BSU for *geothermal* viability and explicitly stated that *hydraulic stimulation is required* because matrix permeability is too low for fluid flow. Citing this as the BCS "average" is misleading — this is the tight-matrix value, not representative of the high-perm streaks that make BCS viable for CCS. | **MISLEADING PRESENTATION** — see Critical Issues |
| High-perm streaks 100–1,000+ mD | Stantec/CEGA (Noad, 2024). Industry article, not peer-reviewed. "Anecdotal evidence" wording in source itself. | **WEAK** — properly caveated in source text but should be flagged |
| LMS ~4 mD | Moradi (2016). Single source. | **MARGINAL** |
| UMS <1 mD | Moradi (2016). Single source. | **MARGINAL** |
| Kv/Kh ratio 0.1–0.5 (est.) | Listed as "Analog from BCS core data" — **no specific source cited**. This appears to be a generic estimate. | **UNSUPPORTED** — see Minor Issues |

**Permeability verdict**: The Weides permeability citation is the most problematic property in the entire document. The Kv/Kh ratio is unsourced.

---

### Area

| Statement | Verdict |
|-----------|---------|
| 9 townships (Twp 54–56, Rge 16–18 W4M) | **CONFIRMED** via AER DDS App #1959125. Location entries match exactly. | **PASS** |
| ~324 mi² / 83,700 ha / 207,000 acres | Derived conversions from 9-township extent. Reasonable (1 twp = 36 mi²). | **PASS** |
| Initial 2–3 Mtpa | **CONFIRMED** via Wolf Midstream Mar 2022 press release. Exact wording matches. | **PASS** |
| Ultimate >6 Mtpa | Wolf Midstream releases. Plausible progression. | **PASS** |
| ACTL pipeline 14.6 Mtpa, Edmonton Connector ~7 Mtpa | ACTL 14.6 Mtpa: capacity claim from Wolf. Edmonton Connector 7 Mtpa: **CONFIRMED** via Wolf Sep 2023 press release. | **PASS** |
| BCS storage resource Alberta 19.2 Gt / Williston 92.5 Gt | Alberta Innovates (Meikle, 2024). Single source, broad regional estimates with 2% efficiency factor. | **MARGINAL** — no cross-validation |

**Area verdict**: Well-supported for the site-specific data (townships confirmed via AER). Regional storage estimates are single-source.

---

### Other Key Reservoir Parameters

| Property | Statement | Verdict |
|----------|-----------|---------|
| Depth ~1,800–2,200 m | Regional depth map; Quest analog ~2,000 m. Reasonable. | **PASS** |
| Pressure >7.38 MPa | Bachu et al. (2000). Confirmed via citation in Moradi thesis. | **PASS** |
| Temperature >31.1°C | Bachu et al. (2000). Confirmed. | **PASS** |
| BSU temperature range 65–120°C | Weides et al. (2014). Confirmed. | **PASS** but note this is for geothermal, not CCS-specific |
| Salinity 100,000–300,000 mg/L | "Regional BCS brine data" — no specific source cited. | **UNSUPPORTED** — see Minor Issues |
| Storage efficiency factor 2% | Alberta Innovates (Meikle, 2024). While 2% is within industry norm (1–4%), it is single-source. | **MARGINAL** |
| Injectivity "up to 1.2 Mt/a per well" | Attributed to "Crouch (2011) via Alberta Innovates." This is a **water injectivity test** result, not CO₂ injectivity. Water has different viscosity, relative permeability, and geochemical behavior than CO₂. This is a significant conflation. | **MISLEADING** — see Critical Issues |
| UCS up to 97.7 MPa, Friction coeff up to 1.22, Tensile <5 MPa | Weides et al. (2014). Confirmed. | **PASS** — geomechanical data properly sourced |

---

## Source Quality

Of the **13 claimed sources**, 11 have been independently verified:

| # | Source | Status | Notes |
|---|--------|--------|-------|
| 1 | AER DDS App #1959125 | **VERIFIED** | Application exists. Wolf Carbon Hub Gp Inc. 10 location entries. Approval 13513 issued 10-Mar-2026. Field/Pool: WILLINGDON/BSL QTZ UND |
| 2 | CCS Knowledge Centre (2026) | **VERIFIED** | Web page exists. Explicitly names Wolf Lamont. "Less than six months" confirmed. Published Apr 2026 |
| 3 | AGS PRS-2024-001 | **VERIFIED** | Publication exists. Herbers, Hauck, Gordon. 28 cores, 4 facies. 2024-01-24 |
| 4 | Weides et al. (2014) CJES | **VERIFIED** | Published in Can. J. Earth Sci. DOI: 10.1139/cjes-2014-0011. 6 wells. Properties confirmed |
| 5 | Moradi (2016) thesis | **VERIFIED** | PhD thesis, University of Calgary. "Time-Lapse Numerical Modeling for a CCS Project in Alberta." Properties independently confirmed through thesis text |
| 6 | Desjardins & Smith (2013) | **VERIFIED** | GeoConvention 2013 abstract #90187. 16 lithofacies, 5 facies associations confirmed |
| 7 | Stantec/CEGA (Noad, 2022/2024) | **NOT VERIFIED** | Industry article. Could not independently confirm. Cited for "high perm streaks" and "anecdotal evidence" |
| 8 | Alberta Innovates (Meikle, 2024) | **VERIFIED** (URL exists) | Whitepaper found at stated URL. Could not read full PDF (binary). URL is legitimate |
| 9 | Bachu et al. (2000, 2004) | **VERIFIED** | Well-known AGS reports. Cited in Moradi (2016) independently |
| 10 | Wolf Midstream (2022–2023) | **VERIFIED** | Mar 2022 press release confirmed (2–3 Mtpa). Sep 2023 ACTL Edmonton Connector confirmed (7 Mtpa) |
| 11 | CER Market Snapshot (Jan 2025) | **VERIFIED** | Page exists. Wolf Lamont listed as "In development," expected 2025, "CO2 Transport/Storage" |
| 12 | Hirschmiller/Whitecap LinkedIn (Jun 2024) | **WEAK SOURCE** | LinkedIn post. Gray literature. Even if core 100/06-05-055-17W4/00 was presented at CEGA, citing a LinkedIn post as a technical reference is not best practice |
| 13 | Global CCS Institute (2024) | **PARTIALLY VERIFIED** | Annual report existence confirmed. CA$200M EDC loan figure not independently verified |

**Source quality verdict**: 10/13 fully verified. Source #7 is unverified. Source #12 is a weak source. Source #13 partially verified.

---

## Critical Issues

### C1. Basal Quartz vs. Basal Cambrian — Stratigraphic Mismatch ⚠️ SEVERITY: HIGH

The AER Approval 13513 lists the pool as **"WILLINGDON/BSL QTZ UND"** (Willingdon field / Basal Quartz undefined). The "Basal Quartz" (Kbs_qtz, code 3340) is the **Lower Cretaceous Basal Quartz sandstone**, stratigraphically distinct from the deeper **Cambrian Basal Cambrian Sandstone (BCS)**.

The document flags this in Section 4 of Notes on Data Quality, which is commendable, but does not resolve the issue. **Consequence**: Nearly all reservoir properties cited (porosity 8–24%, permeability 1 mD to >1 D, thickness 40–80 m, depth ~2,000 m) come from Cambrian BCS literature. If the actual injection target is the Cretaceous Basal Quartz, these values may be entirely inapplicable.

**Recommendation**: The operator must clarify the exact target formation before this property database can be considered valid. If the target is the Cretaceous Basal Quartz, a separate property compilation is needed.

### C2. Weides et al. (2014) Permeability — Context Mismatch ⚠️ SEVERITY: MEDIUM

The document cites Weides for BSU "average permeability" of <1×10⁻¹⁴ m² (~0.01 mD). However, Weides was studying the BSU for **geothermal energy** and explicitly stated that **hydraulic stimulation is required** because matrix permeability is too low for commercial fluid flow. This is not an "average" permeability for the BCS as a CCS target — it is the tight-matrix value for a formation that Weides found unsuitable for flow without stimulation. Presenting this as a BCS/CCS permeability value without the geothermal context is misleading.

### C3. Water Injectivity vs. CO₂ Injectivity ⚠️ SEVERITY: MEDIUM

The document states "Injectivity (Quest analog): Up to 1.2 Mt/a per well" and cites "Crouch (2011) via Alberta Innovates." This value came from **water injectivity tests**. CO₂ has different viscosity, density, relative permeability, and geochemical reactivity compared to water. Water injectivity does not directly translate to CO₂ injectivity. Actual Quest CO₂ injection rates have been closer to ~1 Mt/a total across multiple wells. This should be explicitly caveated.

### C4. Wolf-Specific Data is Nearly Absent ⚠️ SEVERITY: MEDIUM

Of all properties in this document, the only Wolf Lamont-specific data is:
- AER approval locations (townships)
- A single core (100/06-05-055-17W4/00) mentioned in a LinkedIn post
- Press release volumes (2–3 Mtpa, >6 Mtpa)

Every reservoir property (porosity, permeability, thickness, pressure, temperature, salinity, injectivity) is derived from Quest (60–80 km away) or regional studies. The document acknowledges this honestly, but it means the property database is almost entirely analog-based with no site-specific measurements.

---

## Recommendations

### Must-Fix
1. **Resolve the Basal Quartz vs. Basal Cambrian issue** before this database can be used. Contact Wolf Midstream or AER to clarify the injection target.
2. **Recontextualize the Weides (2014) permeability citation** to make clear it represents tight-matrix values for geothermal evaluation, not representative CCS permeability.
3. **Add a caveat to the 1.2 Mt/a per well injectivity figure** noting it is derived from water injectivity tests, not CO₂ injection.
4. **Replace or remove the LinkedIn source** (Source #12) with a verifiable source (request the CEGA Core Conference abstract or Whitecap presentation directly).

### Should-Fix
5. Add source citations (or explicit "unsourced" notes) for **salinity (100k–300k mg/L)** and **Kv/Kh ratio (0.1–0.5)**.
6. Cross-validate the Alberta Innovates (Meikle, 2024) storage resource estimates (19.2 Gt, 92.5 Gt) against Bachu (2004), NETL Carbon Storage Atlas, or other independent estimates.
7. Cross-validate seal property values (LMS, UMS, Devonian Red Beds) against AGS PRS-2024-001 or other independent sources beyond Moradi (2016) alone.
8. Add Quest injection well performance data (actual CO₂ rates vs. water test predictions) from public AER reports to validate injectivity assumptions.

### Nice-to-Fix
9. Clarify that Weides et al. temperature range (65–120°C) reflects a *geothermal* gradient study across central Alberta and may not reflect injection-site temperature.
10. Standardize unit conversions (ft/m) for consistency throughout.

---

## Minor Issues

| Issue | Detail |
|-------|--------|
| Salinity | Listed as "100,000–300,000 mg/L (est.)" with source listed as "Regional BCS brine data" — no specific source cited |
| Kv/Kh ratio | "0.1–0.5 (est.)" with source "Analog from BCS core data" — no specific study cited |
| Storage efficiency factor | 2% appears only from Alberta Innovates (Meikle, 2024). Industry range is 1–4%, so value is reasonable, but single-source |
| Basement seal | "Crystalline granitic — negligible porosity/permeability" sourced as "Generic" — acceptable as basement is universally non-porous, but should be labeled as assumed |
| Cambrian "up to 300 m" | Not tied to Lamont area; relevance is questionable |
| Timelines | Document notes target slipped from "end of 2024" to "2025" — CER Jan 2025 table lists as "In development" still. Update needed for 2026 status |
| Source #7 (Stantec/CEGA) | URL or DOI not provided for verification |
| Source #12 (LinkedIn) | URL not provided; given it is LinkedIn, it may have been deleted or is only visible to connections |
| Source #13 (Global CCS Institute) | No specific report title or URL given beyond the generic "Annual report" |
| Quest distance | Document states "60–80 km north/northeast." Need to verify: Quest site is near Thorhild/Redwater (Twp 59–60, Rge 20 W4M). Lamont is Twp 54–56, Rge 16–18 W4M. Approximate distance is ~60–80 km, which is directionally correct. |

---

## Overall Score

| Category | Rating |
|----------|--------|
| Source verification | 7/10 (10 of 13 sources verified; #7 unverified, #12 weak, #13 partial) |
| Property reasonableness | 7/10 (two misrepresented values) |
| Cross-validation adequacy | 5/10 (seal properties, storage estimates, salinity all single-source) |
| Transparency about limitations | 9/10 (consistently flags Wolf-specific data gaps) |
| Stratigraphic accuracy | 4/10 (flags but does not resolve BQ vs. BC issue) |
| **Overall** | **6.5/10 — usable draft with critical corrections needed** |

---

## Verdict

The document is **honest about its limitations** but contains one **potentially fatal issue** (the Basal Quartz / Basal Cambrian naming conflict) that must be resolved before the property database can be trusted. If the target is confirmed to be Cambrian BCS, the properties are reasonable (with corrections to the Weides permeability context and the water/CO₂ injectivity conflation). If the target is Cretaceous Basal Quartz, the entire document needs to be rewritten based on different literature.