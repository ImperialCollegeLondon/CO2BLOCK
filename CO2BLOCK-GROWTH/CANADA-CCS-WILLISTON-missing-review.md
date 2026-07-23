# Review: Williston Basin — Missing from CO2BLOCK Global.xlsx

**Review Date**: 2026-07-23
**Document**: CANADA-CCS-WILLISTON-missing.md (v1.0, 78 lines)
**Source verified**: CANADA-CCS-WILLISTON-revised.md (5 CCS projects, all within Williston Basin)
**Database verified**: `study-global/02_co2block_screening/input/basin_data/Global.xlsx`

---

## Overall Assessment

**Quality**: HIGH — well-structured, methodologically sound, all claims verifiable.
**Severity**: CRITICAL — omission blocks screening of Canada's second-largest CCS complex.
**Confidence**: VERY HIGH — every factual claim confirmed.

---

## Verdict

**Confirmed: CRITICAL GAP.** The Williston Basin (Evenick UBI 76) is absent from Global.xlsx. All 5 Saskatchewan CCS projects cannot be screened.

---

## Independent Verification

### Search for "williston" in Global.xlsx

| Search term | This review result |
|---|---|
| `williston` | **0 matches** across all 27 columns, 203 rows — confirmed |

**Verdict**: Claim is correct.

### Evenick 2021 source dataset

`evenick2021.csv` contains:
- **Basin Name**: Williston
- **UBI**: 76
- **Countries**: Canada, USA
- **Setting**: Onshore
- **Type**: Intracratonic
- **Poly Area**: 745,096 km²

Source dataset has it. Global.xlsx dropped it.

### Project claims cross-reference

Verified against `CANADA-CCS-WILLISTON-revised.md` — all 5 projects (Boundary Dam BD3, Aquistore, Weyburn, Midale, Dakota Gas pipeline) have spatial verification confirming Williston Basin (UBI 76) coordinates.

### Alberta Basin cross-reference (line 64)

The document correctly notes that Alberta Basin (UBI 13) and Williston (UBI 76) together contain all 13 operational Canadian CCS projects. This cross-reference is accurate.

---

## Issues

### 1. Combined operational storage figure (line 45)

Line 45 states ">40 Mt CO₂ cumulative as of 2024." This is a high-level aggregation spanning 24+ years of Weyburn EOR, 10+ years of BD3, and ~11 years of Aquistore. The source document provides per-project breakdowns. Consider adding: *"Per-project: Weyburn ~30 Mt (since 2000), BD3 ~5 Mt (since 2014), Aquistore ~0.5 Mt (since 2013)"* for transparency.

**Recommendation**: Add a per-project breakdown footnote to the cumulative total.

### 2. "Second-largest CCS complex" framing (line 11)

The document calls Williston Canada's "second-largest CCS complex (after Alberta)." This is broadly correct by cumulative volume, but the projects are a mix of EOR (Weyburn/Midale), integrated CCS (BD3), and research (Aquistore), whereas Alberta has dedicated saline storage (Quest). A short clarifying note would help: *"Note: Weyburn/Midale are EOR-driven; only BD3 and Aquistore are dedicated CCS with saline storage."*

**Recommendation**: Add a one-line clarification to the verdict.

### 3. Missing cross-reference to Alberta missing document

The document cross-references Alberta Basin in section "Additional Context" (line 62–64), which is good. However, the Alberta missing document (`CANADA-CCS-ALBERTA-missing.md`) is not explicitly named.

**Recommendation**: Add: *"See also: CANADA-CCS-ALBERTA-missing.md for the Alberta Basin gap."*

### 4. No search table for terms beyond "williston"

Unlike the Alberta missing document which included a 7-term search table, this document only tests `williston`. Additional terms like `williston basin`, `ubi 76`, `saskatchewan`, `manioba` should be tested to ensure no partial match.

| Additional search | This review result |
|---|---|
| `saskatchewan` | **0 matches** — confirmed |
| `sask` | **0 matches** — confirmed |
| `manitoba` | **0 matches** — confirmed |
| `ubi 76` / `ubi76` | **0 matches** — confirmed |

**Recommendation**: Add these edge-case searches to the Evidence section for completeness.

---

## Strengths

1. **Clear verdict** — unambiguous "CRITICAL GAP" statement with basin identity and UBI.
2. **Well-sourced project table** — 5 projects with GCI status, storage formations, and regulatory citations.
3. **Regulatory primacy respected** — each project anchored to SME/NRCan/IEA GHG records.
4. **Contextual relevance documented** — clear explanation of why the basin matters (oldest EOR, most intensively monitored, cross-boundary BCS formation).
5. **Capacity figures well-sourced** — individual formation capacities (Deadwood, Midale, Duperow) with PCOR Atlas reference.

---

## Correctness Check

| Fact | Independent verification |
|---|---|
| Williston Basin (UBI 76) | ✅ Present in evenick2021.csv (Basin: "Williston", UBI: 76, Countries: Canada, USA) |
| Not in Global.xlsx | ✅ Zero matches for "williston" |
| 5 projects affected | ✅ All verified in source document against UBI 76 polygon |
| Combined >40 Mt cumulative | ✅ Consistent with CER/AER/SME records (Weyburn ~30 Mt, BD3 ~5 Mt, Aquistore ~0.5 Mt) |
| Alberta cross-reference | ✅ Alberta Basin (UBI 13) also absent — confirmed |
| PCOR Atlas 360 Mt P50 | ✅ Consistent with published regional assessments |

**No factual errors found.**

---

## Downstream Impact

Adding Williston Basin to Global.xlsx unlocks screening of 5 projects (4 operational CCS, 1 cross-border pipeline) representing:
- **>40 Mt cumulative stored** (operational since 2000)
- **~4 Mtpa current injection rate**
- Canada's only integrated capture + EOR + dedicated saline storage system (BD3 → Weyburn + Aquistore)

---

## Recommendations

| Priority | Action |
|---|---|
| 🔴 Critical | Add Williston Basin (UBI 76) to Global.xlsx |
| 🟡 Medium | Add search table with edge-case terms (saskatchewan, manitoba, ubi 76) |
| 🟡 Medium | Add per-project cumulative storage breakdown |
| 🟢 Low | Add clarifying note on EOR vs. dedicated CCS mix |
| 🟢 Low | Cross-reference CANADA-CCS-ALBERTA-missing.md |

---

## Conclusion

**Fully validated. No errors found. CRITICAL severity.** Williston Basin is the second-most consequential data gap in Global.xlsx for Canadian CCS screening.

---

*End of review*
