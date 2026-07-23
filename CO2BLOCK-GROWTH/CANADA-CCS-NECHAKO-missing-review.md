# Review: Nechako Basin — Missing from CO2BLOCK Global.xlsx

**Review Date**: 2026-07-23
**Document**: CANADA-CCS-NECHAKO-missing.md (v1.0, 65 lines)
**Source verified**: CANADA-CCS-NECHAKO-revised.md (zero CCS projects; Geoscience BC 2024-08 assessment)
**Database verified**: `study-global/02_co2block_screening/input/basin_data/Global.xlsx`

---

## Overall Assessment

**Quality**: HIGH — well-structured, clearly distinguishes basin gap from CCS-readiness assessment.
**Severity**: MODERATE — basin is missing, but zero CCS projects exist.
**Confidence**: VERY HIGH — confirmed.

---

## Verdict

**Confirmed: MISSING.** The Nechako Basin (Evenick UBI 437) is absent from Global.xlsx. Impact is moderate (zero CCS projects), but the basin contains unique screening data: a definitive "not recommended" CCS assessment and an alternative carbon mineralization pathway (CO2Lock).

---

## Independent Verification

| Search term | This review result |
|---|---|
| `nechako` | **0 matches** across all columns — confirmed |

### Evenick 2021 source dataset

`evenick2021.csv` contains:
- **Basin Name**: Nechako
- **UBI**: 437
- **Countries**: Canada
- **Setting**: Onshore
- **Type**: Backarc - Marginal Sea
- **Poly Area**: 127,206 km²

Source dataset has it. Global.xlsx dropped it.

---

## Issues

### 1. "No CCS projects" is the key mitigant — could be more prominent

The verdict (line 11) correctly states "zero CCS projects," but the word "CRITICAL" is not used (correctly — the impact is lower than Alberta/Williston). However, the verdict paragraph is long and the "zero CCS projects" point could be lost. Consider a concise summary line at the top:

*"**Practical impact**: Zero CCS projects. Basin missing, but screening loss is limited to a negative CCS-readiness flag and an alternative mineralization pathway reference."*

**Recommendation**: Add a standalone "Practical impact" line.

### 2. CO2Lock carbon mineralization (lines 55–57) is unique but over-weighted

The carbon mineralization pathway (CO2Lock SAM, ultramafic peridotite) is interesting but is not a CCS project and exists only at the exploration/assessment stage. The document gives it ~5 lines in a 65-line file (~8%). This is acceptable as context but should not be interpreted as a project gap.

**No change needed** — the document clearly labels this as "not a CCS project."

### 3. Geoscience BC "not recommended" outcome is the most valuable missing data

The fact that Geoscience BC formally recommended against CCS in the Nechako Basin is unique among Canadian basin assessments. This is the most important data point that CO2BLOCK loses by not having the basin in Global.xlsx. It could serve as a negative-case benchmark for screening backarc/marginal sea basins globally.

**Recommendation**: Add a line quantifying the assessment: *"Geoscience BC 2024-08: definitive 'not recommended' for CCS — the only Canadian basin with such a formal negative assessment."* (Note: this is already stated implicitly but could be more explicit.)

### 4. Additional search terms could be tested

The document only tests `nechako`. Testing related terms would strengthen the negative finding:

| Additional search | This review result |
|---|---|
| `nechako` | **0 matches** — confirmed |
| `interior plateau` | **0 matches** — confirmed |
| `ubi 437` / `ubi437` | **0 matches** — confirmed |

**Recommendation**: Add these edge-case terms to the search results.

---

## Strengths

1. **Honest about zero CCS projects** — does not overstate severity.
2. **Unique value proposition documented** — the "not recommended" assessment and carbon mineralization alternative are correctly called out.
3. **Geological detail provided** — Triassic–Cretaceous reservoirs, zeolite cementation, 12 exploratory wells with petrophysical data.
4. **Globally relevant** — correctly identifies that this basin could serve as a negative-case benchmark for screening backarc basins.
5. **Cross-reference to Highway 16 corridor emissions** — provides source–sink context.

---

## Correctness Check

| Fact | Independent verification |
|---|---|
| Nechako Basin (UBI 437) | ✅ Present in evenick2021.csv (Basin: "Nechako", UBI: 437, Countries: Canada) |
| Not in Global.xlsx | ✅ Zero matches for "nechako" |
| Zero CCS projects | ✅ Consistent with source document |
| Geoscience BC 2024-08 "not recommended" | ✅ Consistent with published assessment |
| Area 127,206 km² | ✅ Consistent with evenick2021.csv (127,206 km²) |

**No factual errors found.**

---

## Downstream Impact

Adding Nechako Basin to Global.xlsx would primarily capture:
- **Negative screening data**: a formal "not recommended" outcome for CCS in backarc marginal sea basins
- **Carbon mineralization reference**: in-situ ultramafic peridotite storage pathway (CO2Lock SAM)
- **Geological constraints**: max 10% porosity, 3 mD permeability — a low-quality reservoir benchmark
- **No CCS projects lost** (zero exist)

---

## Recommendations

| Priority | Action |
|---|---|
| 🟡 Medium | Add edge-case search terms (interior plateau, ubi 437) |
| 🟡 Medium | Add "Practical impact" summary line |
| 🟢 Low | Explicitly quantify Geoscience BC negative assessment |
| 🟢 Low | Consider adding Nechako data completeness is acceptable |

---

## Conclusion

**Fully validated. No errors. Basin missing, but zero CCS projects mitigate practical impact.** The document correctly identifies the Nechako Basin's unique value as a negative-case benchmark and carbon mineralization reference.

---

*End of review*
