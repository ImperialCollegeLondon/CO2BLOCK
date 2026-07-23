# Review: Appalachian Basin — Present in CO2BLOCK Global.xlsx

**Review Date**: 2026-07-23
**Document**: CANADA-CCS-APPALACHIAN-missing.md (v1.0, 59 lines)
**Source verified**: CANADA-CCS-APPALACHIAN-revised.md (zero CCS projects; Tier 3–4 assessments only)
**Database verified**: `study-global/02_co2block_screening/input/basin_data/Global.xlsx`

---

## Overall Assessment

**Quality**: HIGH — well-structured, thorough subregional breakdown.
**Severity**: NONE — basin is present.
**Confidence**: VERY HIGH — confirmed.

---

## Verdict

**Confirmed: NO GAP.** The Appalachian Basin (Evenick UBI 39) is correctly present in Global.xlsx at row 8.

---

## Independent Verification

| Attribute | Global.xlsx Value | Verified? |
|---|---|---|
| **Basin Name** | Appalachian | ✅ row 8 |
| **Countries** | USA, Canada | ✅ |
| **Present in Database?** | Yes | ✅ |

### Search for edge-case variants

| Term | Matches | Notes |
|---|---|---|
| `appalachian` | **1** (row 8) | Correct entry |
| `appalachia` | **0** | Evenick name is "Appalachian" — correct |

---

## Issues

### 1. Subregional coverage table (lines 28–35) is excellent

The document correctly notes that St. Lawrence Lowlands, Gaspé, Anticosti, NB, and NS all fall within UBI 39. This level of subregional detail is the best among the "no gap" documents.

**No issue** — this is a strength worth noting.

### 2. "What IS Missing" could omit or downplay US Appalachian context

Line 48 lists "US Appalachian CCS (CAB-CS, Ohio)" as a missing detail. This is cross-border US CCS context, not Canadian. While relevant for context, it risks scope creep in a Canadian basin-audit document.

**Recommendation**: Move "US Appalachian CCS" to a footnote or remove — it belongs in a US basin audit, not a Canadian missing-basin document.

### 3. Cross-reference to St. Lawrence Lowlands could be more prominent

Line 30 cross-references `CANADA-CCS-ST_LAWRENCE_LOWLANDS-missing.md` in a footnote-style comment. Given the St. Lawrence Lowlands is a special case that generates frequent questions, a more prominent callout would help.

**Recommendation**: Add a bolded line: *"See also: CANADA-CCS-ST_LAWRENCE_LOWLANDS-missing.md for Quebec-specific CCS context (subregion of Appalachian)."*

### 4. No search methodology table

Unlike the Alberta and Williston missing documents, this document does not show exhaustive search terms. While the basin is present, adding the search confirmation would make the negative finding more rigorous.

**Recommendation**: Add a brief search-results table showing that "St. Lawrence Lowlands" and "quebec" return zero matches while "appalachian" returns 1 match.

---

## Strengths

1. **Best subregional breakdown** among all "no gap" documents — tables Canadian portions explicitly.
2. **CCS-readiness context** — correctly notes zero projects and no regulatory regimes in QC/NB/NS.
3. **Cross-references to sister documents** — St. Lawrence Lowlands file linked.
4. **Honest about project-level gaps** — INRS assessments, Clean Prosperity estimate, emitter data all documented.

---

## Correctness Check

| Fact | Independent verification |
|---|---|
| Appalachian present in Global.xlsx | ✅ Row 8: "Appalachian", Countries: "USA, Canada" |
| Zero CCS projects | ✅ Consistent with source document |
| St. Lawrence Lowlands within UBI 39 | ✅ Consistent with Evenick classification |
| UBI 39 | ✅ Consistent with evenick2021.csv |

**No factual errors found.**

---

## Recommendations

| Priority | Action |
|---|---|
| 🟢 Low | Move "US Appalachian CCS" to footnote or remove |
| 🟢 Low | Make St. Lawrence Lowlands cross-reference more prominent |
| 🟢 Low | Add brief search-results methodology table |
| 🟢 Low | None critical — document is fit for purpose |

---

## Conclusion

**Fully validated. No errors. No gap.** Best subregional documentation among "present" documents. The Appalachian Basin is correctly in Global.xlsx.

---

*End of review*
