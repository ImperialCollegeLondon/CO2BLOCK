# Williston Basin — Missing from CO2BLOCK Global.xlsx

**Date**: 2026-07-23  
**Source document**: `CANADA-CCS-WILLISTON-revised.md` (5 CCS projects/spine infrastructure in Saskatchewan, all within Williston Basin)  
**Comparison target**: `study-global/02_co2block_screening/input/basin_data/Global.xlsx` (203 basins in database)

---

## Verdict

**CRITICAL GAP**: The **Williston Basin** (Evenick UBI 76, Intracratonic, Onshore) is **absent** from `Global.xlsx`. This basin hosts Canada's **second-largest CCS complex** (after Alberta) with an operational history dating to 2000.

> **Note on project mix**: Williston projects span both CO₂-EOR (Weyburn/Midale) and dedicated CCS with saline storage (BD3, Aquistore), whereas Alberta has dedicated saline storage (Quest). Only BD3 and Aquistore are dedicated CCS with saline storage; Weyburn and Midale are EOR-driven.

---

## Evidence

### 1. Target basin per source document

`CANADA-CCS-WILLISTON-revised.md` states (line 3):

> **Basin**: Williston Basin (Evenick UBI 76) — *Intracratonic, Onshore*  
> **Canadian Extent**: Southern Saskatchewan, southwestern Manitoba, southeastern Alberta  
> **Polygon Area**: 745,096 km² (Evenick 2021)

And in the spatial verification table (lines 127–133), **all 5 projects** are confirmed within this basin.

### 2. Global.xlsx exhaustive search

| Search term | Matches in Global.xlsx |
|---|---|
| `williston` | **0** |
| `williston basin` | **0** |
| `saskatchewan` | **0** |
| `sask` | **0** |
| `manitoba` | **0** |
| `ubi 76` / `ubi76` | **0** |

All terms returned **zero matches**, confirming absence.

---

## Impact: 5 CCS Projects Unrepresented

Without the Williston Basin in `Global.xlsx`, **the following 5 projects cannot be screened in CO2BLOCK:**

| # | Project | Status (GCI) | Storage Formation(s) | Key Regulatory Citation |
|---|---------|--------------|---------------------|-------------------------|
| 1 | **Boundary Dam Integrated CCS (BD3)** | Operational (2014+) | Weyburn-Midale Beds (Mississippian) — CO₂-EOR; Deadwood/Black Island — saline | Sask gov't approval Apr 2011; SME lease; NRCan CEF $240M |
| 2 | **Aquistore (PTRC)** | Operational Research/Demo (2013+) | Deadwood Fm / Black Island Member (Cambro-Ordovician saline aquifer) | SME research permit; NRCan ecoETI funding |
| 3 | **Weyburn CO₂-EOR** | Operational (2000+) | Mississippian Midale Beds (depleted oil field) | IEA GHG Weyburn Project (2000–2012); PTRC White Cap (2012+) |
| 4 | **Midale CO₂-EOR** | Operational (2005+) | Mississippian Midale Beds (depleted oil field) | SME approval; Whitecap Resources ops |
| 5 | **Dakota Gas → Weyburn CO₂ Pipeline** | Operational (2000+) | Transport only (320 km cross-border ND→SK) | NEB/CER pipeline regulation |

**Combined operational storage**: >40 Mt CO₂ cumulative as of 2024¹  
**Current injection rate**: ~4 Mtpa (Weyburn/Midale EOR + BD3/Aquistore saline)

> ¹ Per-project breakdown: Weyburn ~30 Mt (since 2000), BD3 ~5 Mt (since 2014), Aquistore ~0.5 Mt (since 2013). Midale EOR volumes are included in Weyburn-area accounting.

### Why this matters for CO2BLOCK screening

The Williston Basin is:
- Canada's **second most important CCS basin** with 4 operational projects + 1 cross-border pipeline
- Home to the **oldest continuously injected CO₂-EOR storage** in Canada (Weyburn since 2000)
- One of the most intensively monitored CCS sites globally (IEA GHG Weyburn Project, 12+ years of MMV)
- Target of the **Basal Cambrian (Deadwood) aquifer** — the same deep saline formation targeted by Alberta CCS projects, creating potential interference dynamics across provincial boundaries
- The site of Canada's only **integrated carbon capture + EOR + dedicated saline storage system** (BD3 linking capture to both Weyburn EOR and Aquistore saline)
- The only basin where **produced water reinjection infrastructure** has been successfully converted/adapted for CO₂ storage monitoring

---

## Additional Context

### Alberta Basin cross-reference

The Williston Basin is the second "Western Canada Sedimentary Basin" component. The Alberta Basin (UBI 13, also absent from Global.xlsx) covers AB/NEBC, while Williston (UBI 76) covers SK/MB. Together they contain all 13 operational Canadian CCS projects.

> **See also**: `CANADA-CCS-ALBERTA-missing.md` for the Alberta Basin gap documentation.

### Clean Prosperity / Transition Accelerator assessments

- Williston Basin prospective storage: ~3,890 Mt (Transition Accelerator 2023; note: this figure is WCSB-wide, not Williston-only)
- Individual formation capacities well-characterized: Deadwood (Basal Cambrian), Winnipeg-Black Island, Midale Beds, Duperow, Birdbear
- P50 technical storage capacity (Williston SK/MB portion): ~360 Mt (PCOR Partnership Atlas 2015/2017)

---

## Method

- Source document: `CANADA-CCS-WILLISTON-revised.md` (v2.0, 2026-07-23) — 5 projects traced to primary SME/NRCan/IEA GHG records
- Global.xlsx: `study-global/02_co2block_screening/input/basin_data/Global.xlsx`, sheet `Database`, 203 basins, 27 columns
- Cross-check via Python 3.13 with openpyxl confirmed absence
