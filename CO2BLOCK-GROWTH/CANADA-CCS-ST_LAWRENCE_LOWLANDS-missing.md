# St. Lawrence Lowlands — Not a Separate Basin in CO2BLOCK Global.xlsx

**Date**: 2026-07-23  
**Source document**: `CANADA-CCS-ST_LAWRENCE_LOWLANDS-revised.md` (zero CCS projects; INRS pre-feasibility studies; Clean Prosperity 2024 estimate)  
**Comparison target**: `study-global/02_co2block_screening/input/basin_data/Global.xlsx` (203 basins)  
**Evenick (2021) context**: St. Lawrence Lowlands is **not listed as a separate basin**; it falls within the **Appalachian Basin (UBI 39)** fringe.

---

## Verdict

**SPECIAL CASE — NOT A STANDALONE MISSING BASIN**: The St. Lawrence Lowlands is not an independent basin in Evenick (2021). It is a geological subregion of the Appalachian Basin (UBI 39, *Foreland, On-Offshore*). **Appalachian IS present** in Global.xlsx (Countries: USA, Canada).

However, St. Lawrence Lowlands is consistently treated as a **distinct CCS assessment region** in Canadian literature:
- **INRS Chair in Geological CO₂ Sequestration** (2010–2014) — dedicated Lowlands assessment
- **Clean Prosperity (2024)** — separately identifies "St. Lawrence Lowlands" with 2,800–3,200 Mt prospective storage
- **Bédard et al. (2013)** — "CO₂ Geological Storage in the Province of Québec — St. Lawrence Lowlands"
- **Malo et al. (2010)** — "Potential CO₂ Geological Storage Sites in Québec St. Lawrence Lowlands"
- **MERN (Québec)** — treats Lowlands as a separate regulatory/geological unit

---

## Global.xlsx Search Results

| Search term | Matches |
|---|---|
| `st. lawrence` | **0** |
| `st lawrence` | **0** |
| `lowlands` | **0** |
| `quebec` | **0** |
| `ubi 39` | **0** (Appalachian present but not indexed by UBI) |
| `appalachian` | **1** (✅ present as "Appalachian", UBI 39, Countries: USA, Canada) |

---

## Impact

| Feature | Represented? |
|---------|-------------|
| Appalachian Basin (UBI 39) | ✅ Yes — in Global.xlsx |
| St. Lawrence Lowlands as separate unit | ❌ No — not an Evenick basin |
| Quebec CCS assessments (INRS, Clean Prosperity) | ❌ No project-level coverage |
| Saint-Flavien depleted field / Cairnside formations | ❌ Not captured in database |
| Quebec regulatory status (no CCS Act) | ❌ Not captured |

**Screening implication**: A user screening "St. Lawrence Lowlands" in Global.xlsx will find only "Appalachian" and will miss the Quebec-specific CCS context (2.8–3.2 Gt prospective, zero projects, no regulatory regime).

---

## Method

- Source document: `CANADA-CCS-ST_LAWRENCE_LOWLANDS-revised.md` (v2.0, 2026-07-23)
- Global.xlsx: `study-global/02_co2block_screening/input/basin_data/Global.xlsx`, sheet `Database`, 203 basins, 27 columns
- Evenick (2021): St. Lawrence Lowlands is a subregion of Appalachian UBI 39, not a standalone basin
