# Nechako Basin — Missing from CO2BLOCK Global.xlsx

**Date**: 2026-07-23  
**Source document**: `CANADA-CCS-NECHAKO-revised.md` (zero CCS projects; Geoscience BC assessment concludes "not technically recommended")  
**Comparison target**: `study-global/02_co2block_screening/input/basin_data/Global.xlsx` (203 basins in database)

---

## Verdict

**The Nechako Basin** (Evenick UBI 437, Backarc–Marginal Sea, Onshore) is **absent** from `Global.xlsx`. While the Nechako Basin has zero CCS projects and the Geoscience BC Phase 1 assessment (2024) recommends against CCS, the basin's absence from the database means the CO2BLOCK screening cannot represent:

- The Geoscience BC assessment outcome (a rare definitive "not recommended" for CCS)
- The basin's geological characteristics (low-quality Triassic–Cretaceous reservoirs, zeolite cementation)
- The alternative carbon mineralization pathway (CO2Lock SAM property, ultramafic peridotite)

---

## Evidence

### 1. Target basin per source document

`CANADA-CCS-NECHAKO-revised.md` states (line 3):

> **Basin**: Nechako Basin (Evenick UBI 437) — *Backarc - Marginal Sea, Onshore*  
> **Location**: Central British Columbia (Prince George to Houston, Williams Lake to Quesnel)  
> **Polygon Area**: 127,206 km² (Evenick 2021)

### 2. Global.xlsx exhaustive search

A case-insensitive search for `nechako` across all columns of `Global.xlsx` returned **zero matches**.

---

## Impact

### Missing assessment data

| Attribute | Value from Audit | Relevance to CO2BLOCK |
|-----------|-----------------|----------------------|
| Basin type | Backarc–Marginal Sea | Unrepresented in Global.xlsx for this UBI |
| Area | 127,206 km² | 0 km² in database |
| CCS readiness | Not recommended (definitive negative) | Could set a "CCS no-go" flag for backarc basins |
| Storage potential | 197 Mt high-case (high uncertainty) | Uncaptured |
| Carbon mineralization alternative | Ultramafic peridotite (CO2Lock SAM) | Novel storage pathway — not represented in any Global.xlsx basin |
| 12 exploratory wells | All sub-reservoir grade (max 10% φ, 3 mD) | Geological constraint data absent |
| Highway 16 corridor emissions | ~4.6 Mt/yr (>60% biogenic) | Source–sink mismatch unrepresented |

### Basin unique in assessment set

The Nechako Basin is the **only Canadian basin** where a formal government-commissioned geoscience assessment (Geoscience BC 2024-08) definitively recommends against pursuing CCS. If included in Global.xlsx, this basin could serve as a negative-case benchmark for screening backarc/marginal sea basins globally.

---

## Note on CO2Lock carbon mineralization

While not a CCS project, the CO2Lock SAM property on the northeast flank of the Nechako Basin targets in-situ carbon mineralization in ultramafic rocks (peridotite, serpentinite). This is the **only carbon mineralization project** in Canada targeting engineered in-situ mineral storage. If this pathway advances, it would represent a storage class not currently captured by any basin entry in Global.xlsx.

---

## Method

- Source document: `CANADA-CCS-NECHAKO-revised.md` (v2.0, 2026-07-23) — Tier 2–4 sources: Geoscience BC Report 2024-08, BCER, Foresight Canada, CO2Lock
- Global.xlsx: `study-global/02_co2block_screening/input/basin_data/Global.xlsx`, sheet `Database`, 203 basins, 27 columns
- Cross-check via Python 3.13 with openpyxl confirmed absence
