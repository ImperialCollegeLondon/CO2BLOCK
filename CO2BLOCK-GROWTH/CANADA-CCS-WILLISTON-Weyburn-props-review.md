# Adversarial Review: CANADA CCS WILLISTON — Weyburn-Midale props

**File reviewed:** `CANADA-CCS-WILLISTON-Weyburn-props.md` (183 lines)
**Date:** 2026-07-23
**Scope:** Reservoir properties, source tiers, cross-validation, area definition

---

## 1. Source-Tier Assessment

| Tier | Count | Examples |
|------|-------|---------|
| **Tier 1 — Peer-reviewed** | ~5 | White et al., GSA Today (2004); Smith et al., IJGGC (2013); Whittaker et al., Energy Procedia (2011); IJGGC Special Issue Vol. 16 (2013) |
| **Tier 2 — Govt/institutional reports** | ~6 | Sask. Geol. Survey (Burrowes 2001, Whittaker 2005); IEAGHG Poster (2007); BGS (Riding 2005); USask MSc thesis (2008); PTRC |
| **Tier 3 — Secondary/aggregators** | ~4 | MIT CCS Database, Wikipedia, AAPG Education, **Grokipedia** |
| **Untraceable URL** | 1 | OSTI link is a generic servlets/purl (stable but not a proper publication identifier) — borderline OK |

**Verdict: Good mix, heavily weighted toward Tier 1–2. One weak source.**

**ISSUE:** Line 147–148 cites Grokipedia (a wiki aggregation site) for reservoir temperature and pressure. These values (60–63 °C, 12.5–18 MPa) are well-established in White et al. 2004 and Whittaker 2005 — cite those directly instead. Grokipedia has no editorial review.

---

## 2. Cross-Validation Matrix

| Property | Sources | Cross-validated? |
|----------|---------|-----------------|
| Marly thickness (avg ~6 m) | White 2004, Whittaker 2005, AAPG | ✅ 3+ independent |
| Vuggy thickness (avg 15–17 m) | White 2004, Whittaker 2005, AAPG | ✅ 3+ independent |
| Marly porosity (16–38 %, avg ~26 %) | White 2004, USask MSc, OSTI, RITE 2006 | ✅ 4 independent |
| Vuggy porosity (8–21 %, avg ~12 %) | White 2004, USask MSc, OSTI, RITE 2006 | ✅ 4 independent |
| Marly perm (1–150 mD; typical 1–50 mD) | White 2004, RITE 2006, BGS 2005 | ✅ 3 independent |
| Vuggy perm (0.3–500 mD; typical 10–300 mD) | White 2004, RITE 2006, BGS 2005 | ✅ 3 independent |
| Field area (Weyburn ~180 km²) | Burrowes 2001, IEAGHG 2007, Wikipedia 210 km² | ⚠️ Discrepancy |
| Cumulative stored | IEAGHG 2007, MIT, PTRC | ⚠️ See §3 |
| Caprock integrity / 0 % escape | Whittaker 2005 (model-based) | ❌ Single source, model-dependent |

**Verdict: Strong cross-validation for the core petrophysical properties (porosity, permeability, thickness). Weaker for cumulative mass balance and model-based projections.**

---

## 3. Specific Issues

### ISSUE 1 — Fracture gradient vs. injection pressure (line 134)

> _Injection pressure: Maintained below fracture gradient (~10–11 MPa)_

**Problem:** At 1,450 m depth, a typical minimum fracture gradient would be ~15–22 MPa (10–15 kPa/m). The value labeled as fracture gradient (~10–11 MPa) is actually closer to the *injection bottomhole pressure* or the *wellhead pressure*. Injection must be below the fracture gradient (~20+ MPa), not at 10–11 MPa. This appears to confuse injection pressure with the fracture gradient itself.

**Fix:** Re-word to: "Injection bottomhole pressure maintained below fracture gradient (fracture gradient ~X MPa; injection typically at ~10–11 MPa)." Better yet, let the fracture gradient be calculated or cited from a formation integrity test.

---

### ISSUE 2 — Area discrepancy: Weyburn field size (lines 104 vs 119–120)

| Source | Area |
|--------|------|
| Burrowes (2001), Sask. Geol. Survey | ~180 km² (70 sq mi) |
| IEAGHG Poster (2007) | ~180 km² (70 sq mi) |
| Wikipedia | ~210 km² (52,000 acres) |

**Problem:** 52,000 acres = 210 km², not 180 km². A ~17 % discrepancy. The file reports 180 km² but cites Wikipedia in a quote block that says 210 km² — creating internal inconsistency.

**Possible explanation:** The Weyburn unit operated by Cenovus may be smaller than the full Weyburn field structural closure. The file should clarify whether 180 km² is the *unit boundary* (approved EOR area) vs the *structural field outline*. Recommend adding a note.

**Fix:** Add a reconciling note: "Weyburn field structural closure ~210 km²; the Cenovus-operated Weyburn unit approved for CO₂-EOR covers ~180 km²."

---

### ISSUE 3 — Cumulative stored quantity semantics (line 15)

> _>7 Mt by 2007 → ~18 Mt by 2010 → >25 Mt by 2013 → >30 Mt cumulative → ~35+ Mt by 2023_

**Problem:** The file uses "stored" (net retention after accounting for produced CO₂ that is separated and recycled). The distinction between *cumulative injected* and *cumulative stored* (net retained) is implicit, not explicit. A reader could conflate the two.

- Injection rate: ~3 Mt/yr → ~3.6 Mt/yr by 2023
- If gross injected from 2000–2023 at 3.5 Mt/yr avg = ~80 Mt gross
- Stored is ~35 Mt → retention fraction ~44 %

That retention fraction is plausible for a maturing EOR flood but should be stated explicitly.

**Fix:** Add an explicit column or parenthetical: "Cumulative *net* stored (purchased CO₂ minus produced CO₂)."

---

### ISSUE 4 — CO₂ recycle volume ambiguity (line 141)

> _CO₂ recycle: Produced CO₂ separated and re-injected (e.g., ~0.71×10⁶ m³/d at Weyburn)_

**Problem:** The units don't specify whether this is at standard or reservoir conditions. At STP, 0.71×10⁶ m³/d ≈ 1,400 t/d (assuming ρ = 1.87 kg/m³). At reservoir conditions (supercritical, ~650 kg/m³), it would be ~460,000 t/d — absurd for this project. This must be at standard conditions, but that's not stated.

**Fix:** Add "(at standard temperature and pressure)" to the volumetric flow rate.

---

### ISSUE 5 — Boundary Dam supplemental supply (line 10)

> _post-2012 supplemental supply from SaskPower Boundary Dam CCS Facility_

**Problem:** Boundary Dam Unit 3 CCS came online in October 2014. Most of its CO₂ is committed to the nearby Aquistore research project and enhanced oil recovery at other fields. The volume actually delivered to Weyburn/Midale has been minor. This claim needs qualification or a specific source showing volumes.

**Fix:** Either delete "post-2012 supplemental supply from SaskPower Boundary Dam CCS Facility" or add: "(volumes delivered to Weyburn are reported as minor relative to the Dakota Gasification supply; see [specific source])."

---

### ISSUE 6 — Model-based projections presented without uncertainty (lines 143, 153)

> _5,000-year RA shows 0 % escape above caprock_

**Problem:** This is the result of a single numerical simulation (Whittaker 2005). While it was a detailed study, calling it a "risk assessment" and presenting a number like "0 % escape" implies greater certainty than is warranted. The file's own data quality notes (line 176) acknowledge this, but the main table does not.

**Fix:** In the main table, add a footnoted caveat: "Model-based; subject to assumptions about fault-seal integrity and long-term geochemical reactions."

---

### ISSUE 7 — Porosity/permeability from BGS report (line 86)

> _Riding, British Geological Survey (2005): "Porosity and permeability are both highly variable, these parameters average 7% and 1 mD respectively"_

**Problem:** This BGS quote appears to describe different reservoir characteristics — possibly pre-Midale or a different formation within the Williston Basin — not the Midale Vuggy/Marly units. The 7 % porosity / 1 mD permeability average conflicts with every other source in the file (which report Marly φ ~26 %, Vuggy φ ~12 %). This BGS number should not be cited as supporting the Permeability section unless it is confirmed to refer to the same formation.

---

## 4. Overall Assessment

| Criterion | Grade | Notes |
|-----------|-------|-------|
| Reasonable values | B+ | Core petrophysics are well-sourced; fracture gradient confusion is a notable error |
| Source quality | B | Mostly Tier 1–2; Grokipedia and an ambiguous BGS quote weaken |
| Cross-validation | B+ | Good for φ/k/h; weak for cumulative mass, caprock integrity |
| Area clarity | C | 180 vs 210 km² discrepancy unresolved |
| Data quality notes | A | Owns limitations transparently (commercial sensitivity, pre-2012 bias) |
| Internal consistency | B | Self-contradictory on area; otherwise consistent |

**Summary:** This is a well-documented entry overall, with strong cross-validation on the core petrophysical properties. The main liabilities are: (1) the fracture gradient confusion; (2) the unresolved area discrepancy between 180 and 210 km²; (3) presenting model-based projections as hard numbers without inline caveat; and (4) the ambiguous BGS citation. Fix these and the file is solid for the CO2BLOCK database.
