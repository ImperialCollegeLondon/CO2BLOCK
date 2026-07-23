# Origins Project (Enhance Energy) — CCS Properties: Adversarial Review

## 1. Source Quality & Tier Assessment

| # | Source | Designated Tier | Corrected Tier | Issue |
|---|--------|-----------------|----------------|-------|
| 1 | CarbonStorage.io | De facto primary | **Secondary aggregator** | Compiles AER well-file data, but no independent editorial QC. Data may contain transcription errors (see §5). OK for project-well values, but not auditable to original well-file. |
| 2 | Enhance Energy — announcement | Primary | Primary | Valid press release. No issues. |
| 3 | Enhance Energy — approval release | Primary | Primary | Valid press release. No issues. |
| 4 | AER Integrated App Registry | Primary | **Link broken (404)** | URL `https://www1.aer.ca/regulatory-application-integrated-application-registry/1956215` returns 404. AER may have moved or removed this page. Cannot verify Approval 13463 linkage. |
| 5 | AER Regulatory Appeal (PDF) | Primary | Primary | PDF is real but was not human-readable via fetch. Content cannot be independently confirmed from the link alone. |
| 6 | Alberta Major Projects | Government | Government | Reliable. No issues. |
| 7 | NRCan | Government | Government | Reliable. No issues. |
| 8 | **Grokipedia** | Implied tertiary (AGS synthesis) | **LLM-generated wiki** | "Fact-checked by Grok" is not peer review. Grokipedia is an AI-generated encyclopedia with no human editorial board. The AGS data it paraphrases is real, but users should cite the original AGS publications directly (linked in Grokipedia's references section). **Do not use Grokipedia as a primary or secondary source.** |
| 9 | AGS / SE Copernicus | Secondary | Secondary | Legitimate academic/GSCOPE literature. Proper secondary source. |
| 10 | Quest Closure Plan | Contextual | Contextual | Valid AER document for stratigraphic context. |

**Key finding**: Source 4 (AER Application) is dead. Source 8 (Grokipedia) should be replaced with direct AGS citations.

---

## 2. Reservoir Thickness

| Claim | Check | Verdict |
|-------|-------|---------|
| Gross thickness 450 ft (~137 m) at well 4-36-39-25 | Source confirms 450 ft gross thickness. 450 ft × 0.3048 = 137.16 m → "~137 m" is correct. | **PASS** |
| Leduc regional reefal buildup up to 275 m | Consistent with published AGS literature. | **PASS** (but cite AGS directly, not Grokipedia) |
| Lacombe area typical gross interval ~180–300 m | Consistent with reefal buildup range. | **PASS** |
| Net pay not published | Correctly noted. | **PASS** |

---

## 3. Porosity

| Claim | Check | Verdict |
|-------|-------|---------|
| Well average porosity 6 % | From CarbonStorage.io. Below regional avg (~8 %) but within range (5-20 %). Reasonable for a specific well location. | **PASS** |
| Regional Leduc range 5–20 % | Well-established range in AGS literature. | **PASS** (but cite AGS directly) |
| Regional avg ~8 % | Consistent with GSCOPE / SE Copernicus study. | **PASS** |
| Vuggy/moldic dolomitized porosity | Characteristic of Leduc dolomitized reef margins. | **PASS** |

**Cross-validation**: Well value (6 %) < regional avg (8 %) — plausible for a non-reef-crest location. No discrepancy.

---

## 4. Permeability

| Claim | Check | Verdict |
|-------|-------|---------|
| Well range: <10 mD to >100 mD | Original source notation: `<10-+100`. The doc's interpretation is reasonable. The dash-plus notation is unusual. Original likely intends "from <10 to >100 mD" with a typo. | **PASS (with notation)** |
| Regional range: 1–1000 mD | Standard AGS range. | **PASS** |
| Regional avg: ~1–1000 mD | Not really an "average" — it's the range restated. | **Minor: column mislabel** — value is a range, not an average. |
| Heterogeneity noted | Correctly flagged. | **PASS** |

**Cross-validation**: Well range (<10->100 mD) fits within regional range (1-1000 mD). Consistent.

---

## 5. Area — CRITICAL REVIEW SECTION

| Claim | Check | Verdict |
|-------|-------|---------|
| CO2 plume radius 4 km (max 5.6 km) | From AER Appeal 1959099. Reasonable for 28 Mt over 17.5 yr in a carbonate reef. | **PASS** — best public estimate. |
| Plume area ~50 km² (from 4 km radius) | A = πr² = π × 16 = 50.27 km². Rounding to ~50 km² is appropriate. **But**: max 5.6 km gives A = π × 31.36 = 98.5 km², which is not calculated or mentioned. | **Minor omission**: max-area scenario not presented. |
| Lease area 2,491,617 acres | 2,491,617 acres = ~10,082 km². This is ~200× the nominal plume area. Document correctly flags this as potentially total Enhance portfolio. | **PASS** (with caveat in notes). |
| Well location | 00/04-36-039-25W4/0 — valid Alberta well identifier format. LSD 4, Sec 36, Tp 39, Rg 25, W4M. | **PASS** |
| Pool: LACOMBE / TD UND | "TD UND" = pool Type D, Undefined. Consistent with a new sequestration target not previously defined as a producing pool. | **PASS** |

---

## 6. Other Reservoir Parameters — Cross-Validation

### 6.1 Depth, Temperature, Pressure — Internal Consistency

| Parameter | Value | Expected | Check |
|-----------|-------|----------|-------|
| Depth TVD | 6,561 ft (2,000 m) | — | Baseline |
| Datum depth | 6,791 ft (2,070 m) | ~2,000 m + reasonable offset | ✓ ~70 m below TVD, plausible |
| Temperature | 149°F (65°C) | At 2,000 m: 65°C → gradient ~31°C/km (Alberta basin typical: 25-35°C/km) | **PASS** |
| Initial pressure | 2,219 psia (15,300 kPa) | Freshwater hydrostatic at 2,000 m: ~2,841 psi (19,600 kPa). Actual = 2,219 psi → gradient ~0.338 psi/ft. This is significantly **sub-hydrostatic**, consistent with a "pressure-depleted saline aquifer." | **PASS** — consistent with stated depletion state. |
| Max allowable pressure | 20,000 kPa at 2,070 m | ~0.427 psi/ft (just below hydrostatic). Sensible regulatory limit. | **PASS** |

**Conversion check**: 2,219 psia × 6.89476 = 15,299 kPa → "~15,300 kPa" ✓

### 6.2 CO2 Phase

15.3 MPa, 65°C. Critical point of CO2: 7.38 MPa, 31.1°C. Well above both thresholds. Supercritical ✓

### 6.3 Injection Rate & Capacity

| Claim | Check | Verdict |
|-------|-------|---------|
| Phase 1: 1.6 Mt/yr | From CarbonStorage.io. | See §7 below. |
| Full-scale: up to 20 Mt/yr | From Enhance press releases. | **PASS** |
| Phase 1 capacity: 28 Mt | 1.6 Mt/yr × 17.5 yr = 28 Mt. **Exactly matches.** | **PASS** — internal consistency confirmed. |
| Full-scale: "several hundred million tonnes" | Consistent with 20 Mt/yr × ~20+ yr = 400 Mt+. | **PASS** |

### 6.4 Salinity

200,000 ppm TDS. Typical Leduc Fm deep brine range: 50,000–250,000 ppm. At the high end but plausible for a deep isolated reef. No reason to reject.

### 6.5 Scale-Up Check

Phase 1: 1 well × 1.6 Mt/yr → Full-scale: 20 Mt/yr requires ~12.5× the injectivity. This would need ~12-13 additional wells or significant injectivity improvement. Not flagged in the document but worth noting for context.

---

## 7. CRITICAL ISSUE: CarbonStorage.io Unit Label Error

The CarbonStorage.io page lists two conflicting unit labels for the same well:

- In the **Data/Well** section: `"Injection rate (Mt/year): 1,600,000"` — If taken literally, this is **1.6 million megatonnes per year (1.6 trillion tonnes)**, which is physically absurd.
- In the **Overview** section: `"Injection rate 1,600,000 t/yr"` — This uses **tonnes/year** (t/yr), which equals 1.6 Mt/yr.

The document correctly interprets this as **1.6 Mt/yr** (the reasonable value), but **does not flag the source's own unit inconsistency**. This is a data-entry error in CarbonStorage.io. Downstream users copying data from that site could easily misinterpret "1,600,000 Mt/year" at face value.

**Recommendation**: Add a note flagging this CarbonStorage.io unit label discrepancy explicitly. The value should be cited as "1,600,000 t/yr (1.6 Mt/yr)" with a note that the source labels are inconsistent.

---

## 8. Other Discrepancies & Concerns

### 8.1 Grokipedia as a Credible Source

The document uses Grokipedia for regional porosity (5-20 %), permeability (1-1000 mD), and thickness (up to 275 m). These values are **correct** and match established AGS literature. However, the source chain is:
> AGS peer-reviewed publication → Grokipedia (LLM-generated summary, "fact-checked by Grok") → This document

This introduces unnecessary risk of LLM hallucination or simplification error. Recommend replacing all Grokipedia citations with direct AGS references (many are already listed in Grokipedia's own references section).

### 8.2 Broken AER Link

Source 4 URL returns 404. The application registry may have moved, been taken down, or the URL has changed. This is the **primary regulatory reference** for Approval 13463 and the well/pool data. Without it, the approval cannot be independently verified from the document.

### 8.3 AER Appeal PDF Availability

The AER appeal PDF (Source 5) was fetched but returned raw binary data. While the PDF exists, its specific content (plume radius 4 km, max 5.6 km) cannot be verified without local download.

### 8.4 Permeability "Average" Colum

Line 48: `"Leduc average permeability (regional)"` = `"~1–1000 mD"`. This is a **range**, not an average. The original GSCOPE source gives a range of 10⁻¹² to 10⁻¹⁵ m², which is also a range. The column header "average" is misleading. Should read "regional range" or be re-labeled.

### 8.5 Max Plume Area Not Calculated

The document calculates plume area from the nominal 4 km radius (~50 km²) but does not calculate the max-area scenario from the 5.6 km max radius (~98.5 km²). This omission understates potential areal extent.

### 8.6 No Cross-Validation Against Analogues

The document references Quest CCS and ACTL/Clive for context but does not perform any **quantitative cross-validation** of injectivity, storage efficiency, or plume extent against these analogues. For example:
- Quest CCS injects ~1 Mt/yr into the Basal Cambrian Sands. Origins targets Leduc Fm.
- ACTL/Clive injects into Leduc D-3A pool (EOR). Origins targets Leduc TD UND pool.
A brief comparison table would strengthen confidence in the values.

---

## 9. Summary of Findings

| Severity | Issue | Location | Recommendation |
|----------|-------|----------|----------------|
| **HIGH** | Grokipedia used as a cited source | Sources §; Porosity, Permeability, Thickness rows | Replace with direct AGS citations |
| **HIGH** | AER Application URL returns 404 | Source 4 | Find current URL or remove; note the gap |
| **MEDIUM** | CarbonStorage.io unit label error not flagged | Line 79 | Add explicit note about the "Mt/year" vs "t/yr" discrepancy on CarbonStorage |
| **LOW** | "Average permeability" is actually a range | Line 48 | Re-label column as "regional range" |
| **LOW** | Max plume area (5.6 km) not calculated | Line 58 | Add max-case area: ~99 km² |
| **LOW** | No quantitative cross-validation with analogues | §Other Key Parameters | Add brief analogue comparison table |
| **INFO** | Phase 1 well data is wireline-derived, not from core/DST | Notes | Already noted — good. |
| **INFO** | Pressure-depth gradient is sub-hydrostatic (consistent with depletion) | Lines 75-77 | Correctly characterized. |
| **INFO** | Injection rate scale-up of ~12.5× for full development not discussed | Lines 79-80 | Consider noting for completeness. |

---

## 10. Overall Verdict

The document is **well-structured and mostly accurate**. Reservoir property values are internally consistent, cross-validated where possible, and the data quality notes are appropriately cautious. The three most actionable issues are:

1. **Replace Grokipedia citations** with direct AGS references — the data is correct, but the source is not credible.
2. **Fix or remove the dead AER URL** — a primary regulatory link should be functional.
3. **Flag the CarbonStorage.io unit label discrepancy** — to prevent downstream misinterpretation of "1,600,000 Mt/year."

No values need to be changed; only source hygiene and documentation quality improvements are required.
