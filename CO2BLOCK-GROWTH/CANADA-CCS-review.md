# Adversarial Review: CANADA-CCS.csv

**Score: 62/100**

## Overall Verdict: Needs Rework (major structural and completeness issues)

---

## 1. Structural Integrity

**Check** | **Result** | **Notes**
---|---|---
12 data rows? | ✅ PASS | Header + 12 data rows = 13 lines
26 columns? | ✅ PASS | Header has 26 columns, all data rows match
No ragged CSV? | ✅ PASS | All 13 rows have exactly 26 fields
Quoting correct? | ✅ PASS | Fields with commas (e.g., depth values, basin_evidence) are correctly quoted

**Score: 25/25** — No structural issues.

---

## 2. Data Completeness & Consistency

### CRITICAL: All 8 `*_source_url` Columns Are Entirely Empty

| Column | Non-empty | Empty | N/A |
|--------|-----------|-------|-----|
| `depth_source_url` | 0 of 12 | 12 | 0 |
| `porosity_source_url` | 0 of 12 | 12 | 0 |
| `permeability_source_url` | 0 of 12 | 12 | 0 |
| `net_pay_source_url` | 0 of 12 | 12 | 0 |
| `temperature_source_url` | 0 of 12 | 12 | 0 |
| `pressure_source_url` | 0 of 12 | 12 | 0 |
| `storage_capacity_source_url` | 0 of 12 | 12 | 0 |
| `area_source_url` | 0 of 12 | 12 | 0 |

**Every single `*_source_url` cell in every row is empty.** These 8 columns serve no current purpose. Either:
- URLs are embedded in the value column via `value [url]` format (which is the case), making the `*_source_url` columns redundant; or
- The columns were intended for a different formatting scheme and the data extraction failed to populate them.

**Recommendation:** Decide on a single URL-storage convention. Either remove the `*_source_url` columns entirely (since URLs are embedded in values), or move the URLs from the value fields into these columns and remove the `[url]` suffix from values.

### Properties With Values but Missing Source URLs

These cases have a property value with **no source URL at all** — neither inline `[url]` nor in the corresponding `*_source_url` column:

| Row | Project | Column | Value |
|-----|---------|--------|-------|
| 6 | Wolf Lamont Carbon Hub | `pressure_mpa` | `>7.38 MPa (supercritical threshold confirmed); specific formation pressure not published` |
| 7 | Wabamun Carbon Hub | `pressure_mpa` | `~23-33 MPa (hydrostatic estimate at 2,300-3,300 m); may be overpressured` |
| 7 | Wabamun Carbon Hub | `area_km2` | `~1,500 km² (CSA North Area, BSU); full tenure ~2,971 km²` |

The Wolf Lamont pressure case is partially defensible (the text itself says "specific formation pressure not published"), but should still have a source for the supercritical threshold claim or the hydrostatic calculation method. The Wabamun cases are genuine gaps — pressure and area estimates clearly have sources that should be cited.

### Capture-Only Facilities

✅ **Well handled.** NWR Sturgeon Refinery (line 4) and Boundary Dam BD3 (line 10) are clearly identified as capture-only facilities with `N/A (capture facility)` values and proper project_type descriptions.

**Score: 10/25** — Major deduction for universally empty source_url columns + 3 missing source URLs.

---

## 3. Basin Classification Consistency

### Overall Pattern

| Basin | Project Count | Classification |
|-------|--------------|----------------|
| Alberta Basin (WCSB) | 8 | Foreland basin |
| Williston Basin | 4 | Intracratonic basin |

### Assessment

✅ **All 8 Alberta Basin** projects classified as "Foreland basin" — this is appropriate for the Western Canada Sedimentary Basin, which is a classic foreland basin formed by Laramide thrust loading.

✅ **All 4 Williston Basin** projects classified as "Intracratonic basin" — appropriate for this cratonic sag basin. The Williston Basin overlies the cratonic platform of North America and is correctly distinguished from the foreland Alberta Basin.

### Formatting Inconsistencies

| Row | Project | `basin_classification` | Issue |
|-----|---------|----------------------|-------|
| 4 | NWR Sturgeon Refinery CO₂ Capture | `Foreland basin (capture facility — storage via ACTL pipeline to Clive hub)` | Classification field should not contain operational metadata. The parenthetical note belongs in `project_type` or a separate note. |
| 10 | Boundary Dam BD3 | `Intracratonic basin (sag basin)` | Only project with parenthetical note on the classification term. Adds descriptive "sag basin" that none of the other Williston entries have. |

The NWR `basin_classification` field mixes classification with project description. The `basin_classification` column should be purely the basin tectonic classification (`Foreland basin`). The note about it being a capture facility is already present in `project_type` (`CO₂ capture-only`) and in `target_formation` (`N/A (capture facility)`).

Similarly, BD3's "sag basin" note — while geologically accurate — breaks consistency with Aquistore, Weyburn, and Midale which simply say "Intracratonic basin".

**Recommendation:** Strip parenthetical notes from `basin_classification`. Keep it purely taxonomic.

**Score: 18/20** — Classification choices are correct; minor inconsistencies in formatting.

---

## 4. URL Quality

### Spot-Check Results (6 URLs tested)

| URL | Status | Notes |
|-----|--------|-------|
| `https://ags.aer.ca/publications/all-publications/prs-2024-001` | ✅ Resolves | AGS publication page (Herbers et al. 2024) |
| `https://open.alberta.ca/publications/quest-carbon-capture-and-storage-project-d-65-application` | ❌ 403 Forbidden | Access-restricted or requires human interaction |
| `https://cdnsciencepub.com/doi/10.1139/cjes-2012-0137` | ❌ 403 Forbidden | Paywalled journal — may require subscription |
| `https://doi.org/10.1016/j.egypro.2017.03.1630` | ✅ Resolves | Elsevier paper (Jiang et al. 2017, Energy Procedia) |
| `https://open.alberta.ca/publications/` | ❌ 403 + Directory | **Generic URL** — not a specific document; this is Open Alberta's root publications page |
| `https://dds.aer.ca/iar_query/viewApplication.aspx?appnumber=1959125` | ✅ Resolves | AER Integrated Application Registry — Wolf Carbon Hub |

### Additional URLs Checked

| URL | Status | Notes |
|-----|--------|-------|
| `https://cleanprosperity.ca` | ✅ Resolves | Clean Prosperity homepage |
| `https://geoconvention.com/2025/` | ❌ 404 Not Found | Conference page no longer exists |
| `https://www.catf.us/resource/carbon-capture-storage-what-can-learn-from-project-track-record/` | ✅ Resolves | CATF report page |
| `https://rock.geosociety.org/net/gsatoday/archive/14/7/pdf/gt0407.pdf` | ✅ Resolves | PDF download |

### Critically Broken URLs

1. **Shell Polaris — `storage_capacity` uses** `https://open.alberta.ca/publications/`  
   This is Open Alberta's **generic publications directory**, not a specific document. It appears as the source URL for the 300 Mt claim. This same URL appears in the `area_km2` field and as the last entry in `primary_references` for Polaris. **This is a placeholder/broken URL.**

2. **`https://geoconvention.com/2025/`** (used by ACTL/Clive CCS for temperature, pressure)  
   Returns 404. This conference URL likely changes annually. Should be updated to a specific paper/presentation link or archived version.

### 403 Errors

Several Open Alberta URLs returned 403 Forbidden. Open Alberta may block automated crawlers. These URLs may still work for human users with a browser, but the 403 should be verified manually. If they genuinely require special access, consider suggesting alternatives or noting access restrictions.

**Score: 8/15** — One critically broken placeholder URL, one 404, and several 403s needing investigation.

---

## 5. Formatting Issues

### Value-Label Convention

✅ The `value [url]` format is used consistently throughout.

### Encoding

✅ No HTML entities (`&amp;`, `&lt;`, `&gt;`) detected.
✅ No smart/curly quotes detected.

### Primary References

✅ All `primary_references` values are properly semicolon-separated URLs.
✅ All extracted URLs begin with `http` or `https`.

### Column Truncation in Views

⚠️ Several `primary_references` fields exceed 500 characters (Quest: 506, ACTL: 488, NWR: 582, Pathways: 639, Weyburn: 678, Midale: 547). These were truncated in some CSV readers. This is not a file error but could cause confusion in spreadsheet applications. **Recommendation:** Consider whether primary_references should be its own separate reference table to avoid truncation.

**Score: 15/15** — No formatting issues in the data itself.

---

## 6. Comparison Against Source Files

Three source files were compared against the CSV:

### Quest (ALBERTA-Quest-props-revised.md)
- **Depth**: CSV says `1,800-2,132 m MD`; source says `2,105–2,132 m MD at injection wells` (GHGT-15, D-65). The `1,800 m` low end does not appear in the source file — may be from D-65 cross-referencing but is not clearly attributed.
- **All other values**: Porosity (17%/14%/8-21%), permeability (20-500 mD), temperature (60°C), pressure (21 MPa), storage (>7 Mt), area (3,667 km²) — all consistent with the source file.

### Weyburn (WILLISTON-Weyburn-props-revised.md)
- **Depth**: 1,450-1,500 m — consistent.
- **Porosity**: Marly 16-38% avg 26%, Vuggy 8-21% avg 12% — consistent.
- **Permeability**: Marly 1-50 mD, Vuggy 10-300 mD — consistent.
- **Temperature**: 60-63°C — consistent.
- **Pressure**: 12.5-18 MPa — consistent.

### Origins (ALBERTA-Origins-props-revised.md)
- **Depth**: 2,000 m (6,561 ft) — consistent.
- **Porosity**: 6% well avg — consistent.
- **Permeability**: <10 to >100 mD — consistent.
- **Temperature**: 65°C — consistent.
- **Pressure**: 15.3 MPa — consistent.

**Score: 10/15** — Quest depth discrepancy needs resolution; otherwise good fidelity from source files.

---

## Summary of All Issues

| # | Severity | Category | Issue | Row(s) | Recommendation |
|---|----------|----------|-------|--------|----------------|
| 1 | **CRITICAL** | Completeness | All 8 `*_source_url` columns completely empty across all 12 rows | All | Either remove these columns (URLs are embedded in values) or move URLs into them |
| 2 | **HIGH** | Completeness | Wabamun `pressure_mpa` has value but no source URL (neither inline nor in source column) | 7 | Add source URL for the pressure estimate |
| 3 | **HIGH** | Completeness | Wabamun `area_km2` has value but no source URL | 7 | Add source URL for the area estimate |
| 4 | **MEDIUM** | Completeness | Wolf Lamont `pressure_mpa` has no source URL | 6 | Add source for supercritical threshold or state why not available |
| 5 | **HIGH** | URL Quality | Shell Polaris uses `https://open.alberta.ca/publications/` — generic directory, not a specific document | 5 | Replace with specific document URL or mark as placeholder |
| 6 | **HIGH** | URL Quality | `https://open.alberta.ca/publications/` also appears in Polaris `primary_references` | 5 | Same fix |
| 7 | **MEDIUM** | URL Quality | `https://geoconvention.com/2025/` returns 404 | 3 | Update to a specific paper link or archived version |
| 8 | **MEDIUM** | Formatting | NWR `basin_classification` has functional note appended | 4 | Strip parenthetical; keep classification purely taxonomic |
| 9 | **LOW** | Formatting | BD3 `basin_classification` has parenthetical "sag basin" inconsistent with peers | 10 | Standardize format |
| 10 | **LOW** | Fidelity | Quest CSV depth low end (1,800 m) not found in source file | 2 | Clarify source of 1,800 m figure |
| 11 | **LOW** | URL Quality | Multiple Open Alberta URLs return 403 (not necessarily broken, but need manual verification) | Various | Verify manually; note access restrictions if applicable |

---

## Score Breakdown

| Category | Max | Score | Notes |
|----------|-----|-------|-------|
| Structural Integrity | 25 | 25 | No issues |
| Data Completeness & Consistency | 25 | 10 | All source_url columns empty + 3 missing per-property URLs |
| Basin Classification | 20 | 18 | Minor formatting inconsistencies |
| URL Quality | 15 | 8 | One critical placeholder URL, one 404, multiple 403s |
| Formatting | 15 | 15 | Clean formatting throughout |
| Total | 100 | 62 | |
