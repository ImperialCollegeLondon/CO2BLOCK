# Appalachian Basin — Present in CO2BLOCK Global.xlsx

**Date**: 2026-07-23  
**Source document**: `CANADA-CCS-APPALACHIAN-revised.md` (zero CCS projects; Tier 3–4 assessments only)  
**Comparison target**: `study-global/02_co2block_screening/input/basin_data/Global.xlsx` (203 basins in database)

---

## Verdict

**NO BASIN-LEVEL GAP**: The **Appalachian Basin** (Evenick UBI 39, *Foreland, On-Offshore*) **IS present** in `Global.xlsx`. No missing basin entry to report at the basin level.

---

## Confirmed Details

| Attribute | Global.xlsx Value |
|-----------|------------------|
| **Basin Name** | Appalachian |
| **Evenick UBI** | 39 |
| **Countries** | USA, Canada |
| **Present in Database?** | ✅ Yes |

### Search verification

| Search term | Matches in Global.xlsx | Notes |
|---|---|---|
| `appalachian` | **1** (row 8) | Correct entry |
| `appalachia` | **0** | Evenick name is "Appalachian" — correct |
| `st. lawrence` | **0** | Subregion, not a separate Evenick basin |
| `quebec` | **0** | Not a basin name in Evenick taxonomy |

### Canadian portion coverage

The Canadian Appalachian Basin includes the following subregions, all of which fall within Evenick UBI 39. The basin-level entry exists, though subregional detail varies:

| Subregion | In UBI 39? | CCS Activity | Coverage Status |
|-----------|-----------|--------------|-----------------|
| St. Lawrence Lowlands (QC) | ✅ Fringe | 0 projects; assessments | Not separately listed — **see `CANADA-CCS-ST_LAWRENCE_LOWLANDS-missing.md` for Quebec-specific CCS context (subregion of Appalachian)** |
| Gaspé Peninsula (QC) | ✅ | 0 projects | Covered by basin entry |
| Anticosti Island (QC) | ✅ | 0 projects | Covered by basin entry |
| New Brunswick (onshore) | ✅ | 0 projects | Covered by basin entry |
| Nova Scotia (onshore) | ✅ | 0 projects (CCS1 well = non-prospective) | Covered by basin entry |

---

## What IS Missing (Project-Level)

Despite the basin being present, the following **project-level details** from the audit are NOT captured in Global.xlsx (which is a basin classification table, not a project registry):

| Missing Detail | Relevance |
|---------------|-----------|
| **INRS St. Lawrence Lowlands assessments** (2010–2014) | Pre-feasibility studies for Potsdam/Beekmantown/Trenton formations |
| **Clean Prosperity 2024 estimate** (2,800–3,200 Mt prospective) | Quebec's only published storage capacity estimate |
| **Saint John, NB emitter complex** (~3–5 Mt/yr) | Irving Refinery + Canaport LNG — no CCS study |
| **Bécancour industrial zone** (~2.3 Mt/yr) | Valero, Methanex, Yara — no CCS project |
| **Quebec/NB/NS regulatory status** | No CCS Acts, no pore space tenure in any province |

> **Footnote — US context**: US Appalachian CCS projects (e.g., CAB-CS, Ohio) exist but are cross-border US CCS context, not Canadian. They are outside the scope of this Canadian basin-audit document and should be tracked in a dedicated US basin audit.

None of these are basin-level gaps — they are project/policy details beyond the scope of the basin classification table.

---

## Method

- Source document: `CANADA-CCS-APPALACHIAN-revised.md` (v2.0, 2026-07-23)
- Global.xlsx: `study-global/02_co2block_screening/input/basin_data/Global.xlsx`, sheet `Database`, 203 basins, 27 columns
- Cross-reference: See also `CANADA-CCS-ST_LAWRENCE_LOWLANDS-missing.md` for St. Lawrence Lowlands subregion detail
