# Review: CANADA-CCS-ALBERTA-NWR-Sturgeon-props.md

## Summary

Solid capture-facility profile for the NWR Sturgeon Refinery — the anchor CO₂ supplier for the ACTL system. The file correctly establishes this as a capture-only site and cross-references the Clive storage properties. The capture process detail (gasification → Rectisol® → compression → pipeline) is well-documented and the key participants are clearly identified. However, source specificity is weak (named organizations but few URLs or document titles), actual capture volumes are not reported (only design capacity), and the capture-to-storage chain has a monitoring/accounting gap.

---

## 1. Capture-Only Classification — Storage Clarification

**Verdict: Good, with naming concern.**

The file correctly states on lines 5, 13, and 66–68 that storage properties belong to the ACTL/Clive site. Line 13 is explicit: "This is a CO₂ capture facility, not an injection site." The cross-reference to `CANADA-CCS-ALBERTA-ACTL-Clive-props.md` appears in three separate places (lines 13, 67, 68). The Clive file I reviewed contains 161 lines of well-sourced reservoir data (Leduc D-3A: 6.1 % porosity, 30 mD mode permeability, 8.75 m net pay, 18.8 Mt capacity).

**Naming issue:** The file suffix `-props.md` conventionally implies reservoir/storage properties across the CO2BLOCK database. This file contains no storage properties — only capture parameters. Consider whether a naming convention like `-capture-props.md` or `-facility.md` would better distinguish capture profiles from storage-property profiles.

---

## 2. Capture Volumes and Rates — Source Quality

**Verdict: Design capacities are cited; actual capture rates are missing.**

| Parameter | Reported Value | Source Cited | Assessment |
|-----------|---------------|--------------|------------|
| Phase 1 design capacity | 1.2–1.3 Mtpa | NWR (2012), PCL, Wolf Midstream (2020) | Multiple consistent sources; reasonable. But no direct URLs. |
| Recent reported capacity | 1.6 Mtpa | Incorrys (2025) | **Weak.** Incorrys is a consulting firm; no document title or URL. This could be speculative. |
| ACTL pipeline capacity | 14.6 Mtpa | Wolf Midstream, CCS Knowledge Centre | Consistent with Clive file. Adequate. |
| Cumulative stored (ACTL) | >4 Mt by July 2023 | Enhance Energy, Canadian Energy Centre | Acceptable. |
| Government funding | C$558M total | MIT sequestration, Alberta Government | "MIT sequestration" is vague — which dataset? |

**Critical gap:** The file reports *design capacity* but never actual CO₂ captured per year. For an operational facility since 2018 (CO₂ recovery) / 2020 (full ops), actual annual capture volumes should be available. The ">4 Mt cumulative" figure on line 37 implies an average of ~1 Mt/yr across the ACTL system, but this includes CO₂ from other sources (e.g., the Fertilizer plant also feeding ACTL).

The 1.6 Mtpa claim (Incorrys 2025) conflicts with the 1.2–1.3 Mtpa design figure. The Notes section acknowledges this but only offers speculation ("may reflect operational optimization or Phase 2 expansion"). This discrepancy needs a stronger source or an explicit caveat that 1.6 Mtpa is unconfirmed.

---

## 3. Missing Context & Weak Evidence

1. **No URLs or direct citations.** Every source is named by organization only — no document titles, report numbers, or URLs. For a database meant to be reproducible, this is the single biggest weakness. The Clive file (by contrast) uses numbered source references (S1–S11) and provides some access notes.

2. **No actual capture volume data.** Design capacity ≠ actual throughput. The file should report annual captured CO₂ volumes (e.g., from NWR annual reports, or Alberta carbon pricing data).

3. **Emissions offset claim is unsourced.** Line 76 claims "~70% reduction compared to traditional bitumen upgrading" with no source at all.

4. **Lowest-carbon diesel claim is unsourced.** Line 79: "lowest 'wells-to-wheels' CO₂ transportation fuel based on heavy feedstock" — no source, no comparison benchmarks.

5. **"Incorrys (2025)" is insufficient.** No document accessible. This source needs a DOI, URL, or at minimum a report title.

6. **Design-for-expansion context:** Line 36 mentions 2 additional phases, but there is no indication whether these are funded, permitted, or merely theoretical. This needs qualification.

7. **No CO₂ purity/quality specifications.** Line 27 says "high-purity, dry" but gives no numbers (e.g., >99% CO₂? <10 ppm H₂S? <50 ppmv H₂O?). For a capture facility profile feeding a pipeline with multiple sources, composition matters.

8. **No capture rate (% of total emissions).** The Rectisol® unit captures CO₂ from the syngas stream, but what % of total refinery emissions does this represent? The 70% reduction claim (line 76) hints at this but doesn't specify.

9. **No discussion of energy penalty.** Pre-combustion capture via Rectisol® has significant steam/power requirements. This is relevant for net-CO₂ accounting, especially given the facility uses gasification of heavy bottoms (which themselves have upstream emissions).

---

## 4. Capture-to-Storage Chain Documentation

**Verdict: Well-structured process flow, but lacks MMV/accounting linkage.**

Strengths:
- Lines 44–49 provide a clear 6-step chain from gasification through injection.
- Lines 66–69 explicitly verify the chain with entity names.
- The cross-reference to the Clive properties file is repeated in multiple relevant sections.

Gaps:
- **No mass balance or attribution.** The >4 Mt figure (line 37) covers the entire ACTL system, not just NWR's contribution. How do we know what fraction of that came from Sturgeon?
- **No MMV discussion.** The Clive file references an MMV Plan (S6) and a geomechanical study (S11), but the Sturgeon file never mentions how the injected CO₂ is monitored or verified as permanently stored — critical for a CO₂ storage database.
- **No discussion of EOR vs. dedicated storage ratio.** CO₂ used for EOR is partially produced back with the oil and reinjected. The net storage fraction matters. The file treats "CO₂-EOR" and "permanent storage" as synonymous (line 15, 49, 67), which oversimplifies.

---

## Critical Issues

1. **No source URLs or document identifiers anywhere in the file.** Necessary for reproducibility.
2. **Actual capture volumes not reported — only design capacity.** For a facility operational since 2018/2020, this is a significant omission.
3. **1.6 Mtpa Incorrys claim is weakly supported** and conflicts with the design figure without resolution.

---

## Recommendations

1. Add URLs for every source (NWR project page, PCL case study, CCS Knowledge Centre ACTL report, Incorrys report if accessible).
2. Source and report actual annual CO₂ capture volumes (try NWR environmental reports, Alberta carbon pricing submissions, or the CCS Knowledge Centre case study).
3. Either strengthen or explicitly caveat the Incorrys (2025) 1.6 Mtpa claim.
4. Add a section or note on mass balance: how much of the ACTL system total comes from NWR vs. the adjacent fertilizer plant. Clarify the storage attribution.
5. Add a brief note on MMV — reference the Clive MMV Plan and explain that monitoring covers the stored CO₂ regardless of source.
6. Source the 70% reduction and lowest-carbon diesel claims, or remove them.
7. Clarify whether the Phase 2/3 expansion is funded, planned, or hypothetical.
8. Consider renaming to something other than `-props.md` since this file contains no reservoir properties (suggestions: `-capture.md`, `-facility.md`, or keep `-props.md` but add a note in the filename convention docs).

---

## Minor Issues

- Line 37 writes ">4 Mt" — throughout other files the convention is "~4 Mt" or "≈4 Mt". Be consistent.
- Line 38 vehicle equivalency is a useful communication device but has no technical value for the database. Consider moving to a Notes section.
- Line 35 ranges (three × 50,000 bbl/d) use inconsistent units — bbl/d vs. Mtpa — without explaining the conversion.
- Multiple date formats: lines 14, 16, 34, 35 mix m/yyyy, Month yyyy, and Month Year. Standardize.
- Line 17 says "Basin: Alberta Basin / Western Canadian Sedimentary Basin (WCSB)" — WCSB contains the Alberta Basin; they are not synonyms. WCSB is the technically correct basin name.
