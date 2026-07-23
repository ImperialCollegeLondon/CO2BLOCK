# Alberta Basin — Missing from CO2BLOCK Global.xlsx

**Date**: 2026-07-23  
**Source document**: `CANADA-CCS-ALBERTA-revised.md` (8 Alberta CCS projects, all within Alberta Basin / WCSB)  
**Comparison target**: `study-global/02_co2block_screening/input/basin_data/Global.xlsx` (203 basins in database)

---

## Verdict

**CRITICAL GAP**: The **Alberta Basin** (Evenick UBI 13, also called the *Western Canadian Sedimentary Basin / WCSB*) is **absent** from `Global.xlsx`. This is the geological host basin for **all 8 CCS projects** documented in `CANADA-CCS-ALBERTA-revised.md`.

---

## Evidence

### 1. Target basin per source document

`CANADA-CCS-ALBERTA-revised.md` states explicitly (line 11–12):

> **Basin**: Western Canadian Sedimentary Basin (WCSB) — Alberta Basin sub-province  
> **Evenick Basin Name**: Alberta Basin (UBI 13) — *Foreland basin, Canada/USA, Onshore*

And in the spatial verification table (lines 175–185), **all 8 projects** are confirmed within this basin.

### 2. Global.xlsx exhaustive search

The database contains **203 basins** across 27 columns. A case-insensitive search across all columns for the following terms returned **zero matches**:

| Search term | Matches in Global.xlsx |
|---|---|
| `alberta` | **0** |
| `western canada` | **0** |
| `western canadian` | **0** |
| `wcsb` | **0** |
| `williston` | **0** |
| `saskatchewan` | **0** |
| `ubi 13` / `ubi13` | **0** |

### 3. Canada-affiliated basins that DO exist in Global.xlsx

| Basin Name | Countries | Notes |
|---|---|---|
| Appalachian | USA, Canada | Covers Eastern Canada (Maritimes, Quebec) |
| Beaufort-Mackenzie | Canada, USA | Arctic offshore/nearshore |
| Grand Banks | Canada, St. Pierre and Miquelon | Offshore Newfoundland |
| North Slope | USA, Russia, Canada | Arctic Alaska extends geologically into Canada |
| Scotian | Canada, St. Pierre and Miquelon | Offshore Nova Scotia |

These 5 basins span the country, but none of them overlap with the Alberta Basin / WCSB, which covers the entire Alberta–Saskatchewan–NE BC–SW Manitoba region.

---

## Impact: 8 Alberta CCS Projects Unrepresented

Without the Alberta Basin in `Global.xlsx`, **the following 8 CCS projects (all at various stages of maturity) cannot be screened in CO2BLOCK:**

| # | Project | Status | Storage Formation | Regulatory Citation |
|---|---------|--------|-------------------|---------------------|
| 1 | **Quest CCS** (Shell) | Operational (2015+) | Basal Cambrian Sands | ERCB Decision 2012 ABERCB 008 |
| 2 | **ACTL / Clive CO₂-EOR** (Enhance Energy) | Operational (2020+) | Leduc Fm (D-3A pool) | AER Approval 12832 (Dec 2018) |
| 3 | **NWR Sturgeon Refinery Capture** | Operational (2018+) | Via ACTL → Clive | AER D-65 EOR App (Dec 2017) |
| 4 | **Shell Polaris + Atlas Hub** | Under Construction (FID Jun 2024) | Basal Cambrian Sands | Alberta Major Projects #4490 |
| 5 | **Wolf Lamont Carbon Hub** | AER Approved (Mar 2026) | Basal Quartz/Cambrian | AER Approval 13513 (App #1959125) |
| 6 | **Enbridge Wabamun Hub** | Sequestration Agreement (Oct 2025) | Basal Sandstone Unit | Alberta Hub Selection Mar 2022 |
| 7 | **Pathways Alliance CCS Hub** | Regulatory App Filed (Q1 2024) | Basal Cambrian (Cold Lake) | IAAC #89090 |
| 8 | **Origins Project** (Enhance Energy) | AER Approved (Jul 2025) | Leduc Fm (Lacombe) | AER Approval 13463 (App #1956215) |

**Combined operational storage**: >9.8 Mt cumulative as of 2024 (AER ST98 2025)  
**Planned capacity**: >50 Mtpa by 2030

### Why this matters for CO2BLOCK screening

The Alberta Basin is:
- Canada's **only onshore sedimentary basin** with active commercial CCS
- The target for Canada's largest CCS investments (~C$16.5–30B Pathways project)
- The injection formation (Basal Cambrian Sands) is the same regional aquifer targeted by multiple projects, creating interference risks that CO2BLOCK is designed to model
- Home to diverse storage settings: deep saline aquifers (BCS), carbonate reefs (Leduc), EOR transitions (Clive), and deep basal clastics (BSU)

---

## Additional Notable Omissions

| Missing Basin | Relevance |
|---|---|
| **Williston Basin** | Hosts **Boundary Dam CCS** (Saskatchewan, operational 2014) — Canada's other major CCS project. Boundary Dam is a Tier 1 CCS project often paired with Alberta projects in national CCS accounting. |
| **Sverdrup Basin** | Canadian Arctic Islands — potential future storage, studied by NRCan/GSC |
| **Maritimes Basin** | Eastern Canada offshore — some exploration interest for CO₂ storage |

---

## Method

- Source document: `CANADA-CCS-ALBERTA-revised.md` (v2.0, 2026-07-23) — 8 CCS projects traced to primary AER/IAAC regulatory records
- Global.xlsx: `study-global/02_co2block_screening/input/basin_data/Global.xlsx`, sheet `Database`, 203 basins, 27 columns
- Python 3.13 with openpyxl used for exact and fuzzy matching across all rows and columns
- No variation of "Alberta", "WCSB", "Western Canadian Sedimentary Basin", "Williston", or "Saskatchewan" found in any cell
