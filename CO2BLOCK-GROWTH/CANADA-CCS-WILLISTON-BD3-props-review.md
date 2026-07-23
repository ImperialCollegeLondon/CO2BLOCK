# Adversarial Review: CANADA-CCS-WILLISTON-BD3-props.md

**Reviewer**: automated adversarial review  
**Date**: 2026-07-23  
**Document confidence (author's self-assessment)**: HIGH (Weyburn), MODERATE-HIGH (Aquistore)

---

## 1. Reservoir Property Reasonableness

### Weyburn-Midale Beds

| Parameter | Claimed Value | Review | Verdict |
|-----------|---------------|--------|---------|
| Marly thickness | 4–6 m (net 4 m) | Consistent across Whittaker (2005), Burrowes (2006), Jenkins (2018). Thin bed — typical for Marly. | Pass |
| Vuggy thickness | 8–15 m (net 8 m) | Consistent across sources. 0–20 m gross range is realistic for shoal/inter-shoal. | Pass |
| Marly porosity (log) | 25–26% | Supported by El-Sayed (1993) and Jenkins (2018). High but plausible for leached dolomite. | Pass |
| Marly porosity (core/stressed) | ~17% | IEA GHG Phase I (2008). Lower than log due to stress sensitivity — the gap is explained. | Pass |
| Vuggy porosity | 10–12% (avg), shoal up to 20% | Jenkins (2018), Njiekak (2013), GSA Today (2004) all agree. | Pass |
| Marly perm | 6–10 mD matrix, fracture-enhanced | Well-constrained across 5+ sources. Wide range (1→100 mD) reflects fracture contribution. | Pass |
| Vuggy perm | 6–50 mD matrix, shoal up to 500 mD | Consistent. The 500 mD upper bound comes from vuggy shoal — realistic for leached carbonate. | Pass |
| Area | 180–220 km² | Two sources say 180 km², one says 220 km². Minor spread; ~200 km² is a fair central value. | Pass |
| Depth | ~1,450 m (1,400–1,562 m) | Consistent with Williston Basin Mississippian structure. | Pass |
| Pore pressure | Initial ~14 MPa; current ~15 MPa (12.5–18) | ~9.7 kPa/m gradient — normal. Range near injectors (up to 18 MPa) reflects flood-induced overpressure. | Pass |
| Temperature | 60–63°C | ~40–43°C/km — reasonable for the basin. | Pass |
| Formation water TDS | 35–110 g/L | Wide range is realistic — mixing of original formation water with injected water. | Pass |
| **CO₂ density: ~700–800 kg/m³ at 15 MPa, 60°C** | Calculated | **ISSUE**: At 15 MPa and 60°C, CO₂ density is ~600–700 kg/m³, not 700–800. The 800 kg/m³ value is more appropriate for Aquistore conditions (39 MPa, 110°C). The range should be revised to **~600–750 kg/m³** for Weyburn conditions. | **Minor — Flag** |
| CO₂ injection rate 350 mmcf/d = ~18,500 t/d | Jenkins (2018) | 350 MMcf/d × 0.0283 m³/Mcf × 1.87 kg/m³ (std) ≈ 18,500 t/d **at surface conditions**. The entry annotates this as "at reservoir conditions" which is incorrect — reservoir volume would be ~350 MMcf/d ÷ ~150 (compression ratio) ≈ 2.3 MMcf/d. This is a units/conversion labeling error, not a property error. | **Minor — Flag** |
| Cumulative CO₂ stored | ~30 Mt | Plausible over 25+ years of EOR at ~1 Mtpa+ injection. | Pass |

### Aquistore (Deadwood / Black Island)

| Parameter | Claimed Value | Review | Verdict |
|-----------|---------------|--------|---------|
| Deadwood thickness | ~250 m gross | Zambrano-Narvaez (2015). Consistent with regional Cambrian isopachs. | Pass |
| Black Island Member | ~40 m | Same source. Consistent. | Pass |
| Perforated interval | 200 m (4 zones over 3,173–3,366 m) | Talman (2025), Jiang (2018). Well-constrained from completion data. | Pass |
| **Deadwood porosity** | **6% (core) / 7% (log)** | Low but typical for deeply buried Cambrian quartz arenite (~3.2 km). Core vs log discrepancy is small and explained. | Pass |
| **Black Island porosity** | **~21% (regional)** | Bedard et al., MDPI (2024). **ISSUE**: 21% porosity at ~3,100 m depth is unusually high for sandstone. Typical deep saline aquifers at 3+ km have 5–15% porosity. The source describes this as "regional" — it may represent shallower or updip portions of the Winnipeg Formation. No site-specific measurement is cited. This value should not be used for site-specific storage calculations without confirmation from the Aquistore well logs. | **Significant — Flag** |
| Deadwood perm (core) | 5 mD | Zambrano-Narvaez (2015). Reasonable for tight Cambrian sandstone. | Pass |
| Deadwood perm (field-scale) | 20–150 mD | White (2017), Rangriz Shokri (2019). Scale effect is plausible and acknowledged. | Pass |
| Site area (model/seismic) | 30–34 km² | Jiang (2018), Roach (2017). Well-defined survey footprint. | Pass |
| Depth | 3,130–3,350 m | Consistent across all sources. | Pass |
| Initial pressure | 34.2–35 MPa | ~10.7–11.0 kPa/m gradient — slightly overpressured. Plausible for deep Cambrian. | Pass |
| Temperature | 110–119°C | **ISSUE ACKNOWLEDGED** in Notes (§6). 9°C range may reflect measurement method (DTS vs gauge) or injection cooling. | Pass (flagged) |
| Fracture pressure | ~48 MPa (14.9 kPa/m) | Mini-frac based. Reasonable for depth. | Pass |
| Formation water TDS | ~330 g/L (halite-saturated) | Consistent with deep Williston Basin brines. | Pass |
| **Caprock identification** | **Icebox shale + Prairie Evaporite** | **The Prairie Evaporite is listed as a "secondary seal" but it is ~600 m shallower than the reservoir (~2,515 m vs ~3,200 m). This is not a "secondary seal" for the storage complex — it is a separate regional seal above. Terminology should be clarified: "secondary" could be misinterpreted as directly overlying the reservoir.** | **Minor — Flag** |

---

## 2. Source Tier Assessment

| Tier | Definition | Count in Doc | Examples | Adequate? |
|------|------------|--------------|----------|-----------|
| Tier 1 | Peer-reviewed journal | ~10+ | GSA Today, IJGGG (×6), SPE J (×2), MDPI | Yes |
| Tier 2 | Peer-reviewed conference, SPE paper, book chapter | ~6 | ACG (2015), SPE (1993), SSRN (2025), Wiley (2025) | Yes |
| Tier 3 | Industry report, operator blog, presentation | ~4 | Jenkins (2018), SaskPower blog, IEA GHG Phase I report, Whittaker (2005) | Yes |

**Issues:**

- **Bedard et al., MDPI EngProc (2024)** is labeled "Peer-reviewed journal" but *Engineering Proceedings* is a **proceedings series**, not a peer-reviewed journal. This inflates the apparent source tier. Should be reclassified as Tier 2 (conference proceedings). This is where the questionable 21% Black Island porosity originates.
- **IEA GHG Weyburn Final Report (2008)** is labeled "International research report" — this is a Tier 2 source at best (project report, not peer-reviewed in the journal sense), though it carries high authority for this specific project.
- **Jenkins, CO₂ Conference (2018)** is "Industry presentation" — Tier 3. However, it provides many of the most specific values (stressed perm, net pay, Vuggy subzones). These values are cross-checked against higher-tier sources, so the risk is contained.
- **Good practice**: The document properly cross-references lower-tier claims against higher-tier sources for most key properties.

---

## 3. Cross-Validation Status

| Property | Sources | Consistent? | Notes |
|----------|---------|-------------|-------|
| Marly thickness | Whittaker, Burrowes, Jenkins (×3) | Yes | Narrow range (4–6 m gross, 4 m net) |
| Vuggy thickness | Whittaker, Burrowes, Jenkins (×3) | Yes | 8–15 m gross, 8 m net |
| Marly porosity | Jenkins, El-Sayed, IEA GHG, Njiekak (×4) | With caveats | Log (25–26%) vs core (17%) discrepancy well-explained |
| Vuggy porosity | Jenkins, Njiekak, GSA Today (×3) | Yes | 10–12% avg with facies control |
| Marly perm | Jenkins, El-Sayed, GSA Today, Njiekak (×4) | Yes | Consistent range 1–100 mD, avg 6–10 mD |
| Vuggy perm | Jenkins, GSA Today, Njiekak (×3) | Yes | Wide but consistent range |
| Weyburn area | El-Sayed, Njiekak, Jenkins (×3) | Yes | 180 vs 220 km² — minor, central ~200 km² |
| Weyburn depth | GSA Today, Wegelin (×2) | Yes | 1,400–1,562 m |
| Aquistore porosity (Deadwood) | Zambrano-Narvaez, Roach (×2) | Yes | 6–7% — tight agreement |
| **Aquistore Black Island porosity** | **Bedard (×1 only)** | **No cross-check** | Single source — needs independent verification |
| Aquistore perm (core) | Zambrano-Narvaez (×1) | Weak | Only one core measurement cited; field-scale calcs from White and Rangriz Shokri provide cross-checks |
| Aquistore temperature | Zambrano-Narvaez, Talman, Roach (×3) | Partial | 110–119°C range acknowledged; no resolution |
| Aquistore area | Jiang, Roach (×2) | Yes | 30 vs 34 km² — consistent |

**Weak spots**: Black Island porosity (single source), Aquistore core permeability (single source), Aquistore temperature (unresolved spread).

---

## 4. Area Definition

### Weyburn-Midale: **ADEQUATE**
- Productive pool area defined by oil-water contact and field boundary.
- Two independent sources agree on ~180 km²; one says 220 km² (possibly including non-productive portions).
- ~200 km² is a reasonable working value.
- The area clearly refers to the **oil pool** extent, not the CO₂ plume or SLA.

### Aquistore: **WEAK**
- No formal storage lease area (SLA) number is presented — the document acknowledges this.
- The 30–34 km² values refer to the **3D seismic survey footprint and model domain**, not a regulatory/contractual area.
- The "regional >100,000 km²" figure is basin-scale and not useful for site-specific calculations.
- **For storage capacity calculations**, the lack of a defined SLA is a gap. If P50/P10 storage resource is needed, the undefined area makes the calculation unconstrained.
- **Recommendation**: Add a note on whether pore-space tenure is defined by the Saskatchewan *Captured Carbon Storage Act*, and if so, what area that covers.

---

## Summary of Issues

| Severity | Issue | Location | Recommendation |
|----------|-------|----------|----------------|
| **Significant** | Black Island Member porosity (21%) from single source (Bedard 2024) without site-specific verification; MDPI EngProc is a proceedings, not a journal | §2B | Downgrade source tier; add caveat that 21% represents regional (potentially shallower) porosity, not site-specific. Seek Aquistore well-log porosity for the Black Island. |
| Minor | CO₂ density range at Weyburn overstated (700–800 → 600–750 kg/m³) | §1E (Table) | Revise to 600–750 kg/m³ for 15 MPa, 60°C |
| Minor | CO₂ injection rate mislabeled as "at reservoir conditions" (should be surface/standard conditions) | §1E (Table) | Remove "at reservoir conditions" label |
| Minor | Prairie Evaporite called "secondary seal" — it's 600 m above the reservoir, not a directly overlying secondary barrier | §2E (Table) | Reword to "regional barrier / secondary containment" |
| Minor | 200 m perforated interval vs 290 m gross reservoir (250+40) inconsistency not fully explained in header | §2A | Clarify that 200 m = net perforated interval, not gross reservoir thickness |
| Minor | Marly porosity dual values (25% log vs 17% core) — no guidance on which to use for different calculation types (volumetrics vs. injectivity) | §1B | Add: "Use 25% for pore-volume-based storage calculations; use 17% for stress-sensitive core-calibrated modeling" |
| Observation | Aquistore area lacks a defined SLA — storage resource calculations are unconstrained by area | §2D | Add reference to Saskatchewan regulatory framework for pore-space tenure |

---

## Overall Assessment

The document is thorough, well-sourced, and the author has done good work flagging their own uncertainties. The Weyburn-Midale section is robust — 25+ years of research provide strong cross-validation for all key properties. The Aquistore section is necessarily weaker due to the single-well dataset, and the document acknowledges this.

**The most actionable finding**: The Black Island Member porosity of 21% at 3.1 km depth from a single proceedings paper is not credible for site-specific use without confirmation from Aquistore well logs. This value could significantly overstate storage potential if used in calculations.

**Other findings** are minor labeling/terminology issues that do not undermine the document's overall utility.

---

*Review methodology: adversarial cross-check of every property value against its source, internal consistency check, source tier verification, area definition clarity assessment.*
