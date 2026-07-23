# Review: Alberta Basin — Missing from CO2BLOCK Global.xlsx

**Review Date**: 2026-07-23
**Document**: CANADA-CCS-ALBERTA-missing.md (v1.0)
**Source verified**: CANADA-CCS-ALBERTA-revised.md (v2.0, 343 lines, 8 projects)
**Database verified**: `study-global/02_co2block_screening/input/basin_data/Global.xlsx` (203 basins, 27 columns)
**Source dataset verified**: `study-global/02_co2block_screening/input/basin_data/evenick2021.csv` (Alberta Basin UBI 13 present)

---

## Overall Assessment

**Quality**: HIGH — well-structured, methodologically sound, all claims independently verifiable.
**Severity**: CRITICAL — the omission blocks screening of Canada's largest CCS concentration.
**Confidence**: VERY HIGH — every factual claim confirmed against the source document and the database.

---

## Verdict

**Confirmed: CRITICAL GAP.** The Alberta Basin (Evenick UBI 13) is absent from Global.xlsx. All 8 Alberta CCS projects cannot be screened. This is the highest-priority missing-basin finding across all Canadian basin audits.

---

## Independent Verification Results

### Claim 1: Alberta Basin (UBI 13) is absent from Global.xlsx

| Search term | This review result |
|---|---|
| `alberta` | **0 matches** — confirmed |
| `western canada` | **0 matches** — confirmed |
| `western canadian` | **0 matches** — confirmed |
| `wcsb` | **0 matches** — confirmed |
| `williston` | **0 matches** — confirmed |
| `saskatchewan` | **0 matches** — confirmed |
| `ubi` | **0 matches** — confirmed |
| `canada` | 8 matches (only in Countries column of 5 basins, plus 3 majority-country hits) — none in Basin Name |
| `sask`, `manitoba`, `ontario`, `quebec`, `cordiller`, `rocky`, `foreland`, `interior` | **0 matches each** — confirmed |

**Extended searches (beyond the original document) also returned zero:**

- `hudson` (Hudson Bay / Foxe basins): 0
- `labrador`, `newfoundland`: 0
- `nunavut`, `northwest territ`: 0

**Verdict**: Claim is correct. The Alberta Basin is entirely absent from Global.xlsx at any column.

### Claim 2: Alberta Basin exists in Evenick 2021 source dataset

Verified: `evenick2021.csv` contains:
- **Basin Name**: Alberta
- **UBI**: 13
- **Countries**: Canada, USA
- **Setting**: Onshore
- **Type**: Foreland
- **Poly Area**: 1,602,675 km²

The source dataset has it. Global.xlsx dropped it. This is a **data pipeline loss**, not a source data issue.

### Claim 3: The 5 Canada-affiliated basins that DO exist in Global.xlsx

| Basin | Row in Global.xlsx | Countries | Verified? |
|---|---|---|---|
| Appalachian | row 8 | USA, Canada | ✅ |
| Beaufort-Mackenzie | row 19 | Canada, USA | ✅ |
| Grand Banks | row 65 | Canada, St. Pierre and Miquelon | ✅ |
| North Slope | row 120 | USA, Russia, Canada | ✅ |
| Scotian | row 164 | Canada, St. Pierre and Miquelon | ✅ |

**Verdict**: Correct. These are the only Canada-affiliated basins in the database. None overlap with the Alberta Basin / WCSB.

### Claim 4: 8 Alberta CCS Projects cannot be screened without this basin

Confirmed from `CANADA-CCS-ALBERTA-revised.md` — all 8 projects have spatial verification against Evenick UBI 13 polygon. Every project table entry marks basin verification as ✅ WCSB / Alberta Basin.

### Claim 5: While the source document is robust, the omission is in the database layer

**Important nuance**: The source document (v2.0) is excellent — the *earlier* review (`CANADA-CCS-ALBERTA-review.md`) flagged Wolf Lamont temporal anomaly as Issue #1, and the revised document added a second search round (2026-07-23) to resolve it. That fix is reflected in the current missing document's reference to the revised source.

---

## Issues

### 1. Williston Basin cross-reference is correct but could be misleading

The "Additional Notable Omissions" section lists the Williston Basin (lines 86–87), noting it hosts Boundary Dam CCS. The finding is correct. However, the Williston Basin has its own dedicated missing-document (`CANADA-CCS-WILLISTON-missing.md`). A cross-reference to that file would help downstream users connect the two documents.

**Recommendation**: Add `See also: CANADA-CCS-WILLISTON-missing.md` to the Williston row in section "Additional Notable Omissions."

### 2. The "Additional Notable Omissions" section shortlists only 3 basins

The document only lists Williston, Sverdrup, and Maritimes. However, the remaining-basins audit (`CANADA-CCS-REMAINING_BASINS_ZERO.md`) identifies **34 Canadian basins** absent from Global.xlsx. While most have zero CCS projects and are correctly outside the scope of a CCS-focused missing-basin report, a brief note would clarify that the 3 mentioned are not exhaustive of all missing Canadian basins.

**Recommendation**: Add a sentence: *"For a complete inventory of all 34 Canadian basins absent from Global.xlsx (most with zero CCS activity), see CANADA-CCS-REMAINING_BASINS_ZERO.md."*

### 3. Minor: Sverdrup Basin claim (line 87)

Line 87 states Sverdrup has "potential future storage, studied by NRCan/GSC." This is true but the studies are generic Arctic sedimentary assessments, not CCS-specific feasibility studies. The dedicated Sverdrup CCS audit (`CANADA-SVERSDRUP-CCS.md`) is more precise here.

**Recommendation**: Replace `studied by NRCan/GSC` with `mentioned in NRCan/GSC Arctic assessments (no dedicated CCS study)` for precision, and cross-reference the Sverdrup audit file.

### 4. Method section could include Evenick source dataset comparison

The Method (lines 92–97) describes Python-based searching of Global.xlsx. It does not mention that the missing document also checked `evenick2021.csv` (the source) to confirm the Alberta Basin exists in Evenick's global classification. Adding this would strengthen the case that the omission is a Global.xlsx-specific issue, not a basin-classification error.

**Recommendation**: Add: *"Cross-referenced against evenick2021.csv (Evenick 2021 global basin dataset) — Alberta Basin UBI 13 confirmed present in source classification."*

---

## Strengths

1. **Search exhaustive** — 7+ search terms (including edge cases) with zero matches. Confirmed independently in this review.

2. **Impact quantification** — the 8-project table with regulatory citations, status, formations, and cumulative storage metrics is directly actionable for database maintainers.

3. **Method transparency** — source document, target database, Python tooling all specified. Reproducible.

4. **Acknowledged additional omissions** — Williston, Sverdrup, Maritimes basins listed in section 5, creating a prioritization signal.

---

## Correctness Check

| Fact in document | Independent verification |
|---|---|
| Alberta Basin (UBI 13) | ✅ Present in evenick2021.csv (Basin Name: "Alberta", UBI: 13, Countries: Canada, USA) |
| Not in Global.xlsx | ✅ Zero matches across 7+ search terms, 203 basins, 27 columns |
| 5 Canada basins present | ✅ Confirmed — Appalachian, Beaufort-Mackenzie, Grand Banks, North Slope, Scotian |
| 8 projects affected | ✅ All verified in source document with spatial coordinates against UBI 13 polygon |
| Method uses openpyxl | ✅ Reproducible |
| "Western Canadian Sedimentary Basin" = Alberta Basin | ✅ Consistent with Evenick classification (UBI 13, Onshore, Foreland, Canada/USA) |

**No factual errors found.**

---

## Downstream Impact

If the Alberta Basin is added to Global.xlsx, the following screening capabilities are unlocked:

- **Quest CCS** — operational saline aquifer storage (Basal Cambrian Sands)
- **ACTL/Clive** — operational CO₂-EOR in carbonate reef (Leduc)
- **NWR Sturgeon** — integrated capture-transport-storage chain
- **Shell Polaris/Atlas** — under construction (~2028), purpose-built saline hub
- **Wolf Lamont** — AER-approved, open-access hub
- **Enbridge Wabamun** — sequestration agreement, deep basal clastics
- **Pathways Alliance** — largest planned CCS network globally (40 Mtpa)
- **Origins** — AER-approved, saline aquifer, pore-space precedent

Combined: **>50 Mtpa planned capacity** — a material fraction of Canada's 2030 emissions reduction targets.

---

## Recommendations

| Priority | Action |
|---|---|
| 🔴 Critical | Add Alberta Basin (UBI 13) to Global.xlsx with Evenick 2021 parameters (Foreland, Onshore, Canada/USA, 1,602,675 km²) |
| 🔴 Critical | Add Williston Basin (UBI 76) to Global.xlsx (separately documented in CANADA-CCS-WILLISTON-missing.md) |
| 🟡 Medium | Add cross-references to Williston and remaining-basins documents |
| 🟡 Medium | Add Evenick source-dataset cross-check to Method section |
| 🟢 Low | Refine Sverdrup claim for precision |
| 🟢 Low | Note remaining 34 basins are covered in separate file |

---

## Conclusion

**Fully validated. No errors found. CRITICAL severity.** The Alberta Basin omission is the single most consequential data gap in Global.xlsx for Canadian CCS screening. The source document is correct, thorough, and actionable.

---

*End of review*
