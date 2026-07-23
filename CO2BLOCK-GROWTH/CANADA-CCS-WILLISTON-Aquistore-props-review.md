# Aquistore CCS Project — Reservoir Properties: Adversarial Review

**File reviewed**: `CANADA-CCS-WILLISTON-Aquistore-props.md`
**Review completed**: July 2026

---

## 1. Reservoir Property Reasonableness

### 1.1 Thickness — Adequate

The thickness values are internally consistent and well cross-validated.

| Sub-claim | Check | Verdict |
|-----------|-------|---------|
| Total interval 200–250 m, mean 219 m | Self-consistent; Deadwood (~145 m) + Black Island (~44 m) + Icebox (~25 m) = ~214 m | ✅ |
| Net-to-gross 51% → ~112 m net pay | Matches "injection zone >100 m" | ✅ |
| Perforated total 88 m | Reasonable fraction of net pay | ✅ |

No issues.

### 1.2 Porosity — Two Concerns

**Concern A: 21% porosity for Black Island Member quartz arenite (Hersi & Iqbal 2025)**

This value is cited from an *Engineering Proceedings* conference paper (MDPI Eng. Proc., variable peer review). The Black Island Member at the Aquistore injection site lies at ~3,130 m depth. At that depth, 21% intergranular porosity in a Cambro-Ordovician sandstone is **anomalously high**. For context:
- The Aquistore-specific Deadwood core plugs average 6% porosity
- Log-based Deadwood porosity averages 7–11%
- Even the high-porosity quartz arenite in comparable deep basins (e.g., Cambrian of Algeria, Ordovician of North Africa) rarely exceeds 15–18% at >3,000 m

The Hersi & Iqbal study covers the Black Island Member across a 450×500 km region. The 21% value likely comes from **shallower, less-buried portions** of the basin (near the northern subcrop edge), not from the Aquistore well location specifically. This is a **critical spatial context issue** — the value is probably correct regionally but misleading if attributed to the Aquistore injection site.

**Recommendation**: Add a note specifying that the 21% value is a regional maximum from outcrop/subcrop areas; the Black Island Member at the Aquistore well depth likely has porosity closer to 7–12%.

**Concern B: Wide porosity range without resolution**
The table lists porosities spanning 6% (core) to 21% (regional quartz arenite) — a 3.5× range. This is wider than ideal. The existing "Data Quality" notes do acknowledge discrepancies, but the table presents all values with equal weight.

**Recommendation**: Add ranking or flags (site-specific vs. regional) directly in the table.

### 1.3 Permeability — Adequate but Significant Spread

The 5 mD (core) vs. ~70 mD (effective from kh/h) vs. 300 mD (regional) spread is acknowledged in the notes. This is good. The injectivity index and transmissibility data provide operationally grounded cross-checks.

One minor concern: the calculated 70 mD uses h ~130 m — but what defines this thickness? If it is total interval (~200 m) × NTG (51%) = ~102 m, then k ≈ 88 mD. If it is the "injection zone thickness" (>100 m), similar. The 70 mD figure is plausible but the derivation should state what h is used.

### 1.4 Pressure, Temperature, CO₂ Density — One Error

**The CO₂ density calculation uses inconsistent temperature.**

The table reports:
- Reservoir temperature: 119–120°C (lines 137–138)
- CO₂ density: ~800 kg/m³ "calculated for reservoir conditions of P=39 MPa and T=110°C" (line 141)

The density calculation uses **110°C**, which is ~10°C below the stated reservoir temperature. At 39 MPa:
- At 110°C: ρCO₂ ≈ 790–830 kg/m³ (depending on impurities)
- At 120°C: ρCO₂ ≈ 700–760 kg/m³

This is a **~50–80 kg/m³ overestimate (≈7–11%)**. If applied to storage capacity calculations, it would produce a proportional overestimate.

**Recommendation**: Recalculate density at the correct reservoir temperature (119–120°C), or clarify that 110°C was a different measurement point (e.g., cooler part of the reservoir). If 110°C is correct (e.g., near-wellbore cooling from injection), note that explicitly.

### 1.5 Storage Capacity — Reasonable

The P10–P90 range of 8.4–27.1 Mt over 34 km² is broadly consistent with a quick volumetric check:
- V = 34 km² × 200 m × 0.06 (ϕ) × 0.51 (NTG) = 2.08×10⁸ m³ pore volume
- Pore volume × ρCO₂ (800 kg/m³) × E (0.05–0.15) = 8–25 Mt

The gigatonne-scale regional estimate is qualitatively plausible for a 450×500 km aquifer. The 974 km³ pore space for the Black Island Member (Hersi & Iqbal, 2025) is a **regional** number — it represents the entire formation's pore volume across SE Saskatchewan, not anything specific to Aquistore. This is not flagged in the document.

**Recommendation**: Add a note clarifying that the 974 km³ and "gigatonne-scale" numbers are regional resource estimates, not site-specific capacity for the Aquistore injection area.

### 1.6 Seal Thickness — One Flag

The Icebox Member (primary seal) thickness is given as "~15 m (some sources say 25–30 m)". A 2× range on the primary seal is significant and should be resolved if possible. Fifteen metres of shale is likely adequate caprock at this depth with low overpressure, but if a reviewer is strict, the uncertainty should be bounded.

---

## 2. Source Tier Assessment

| Source | Tier | Notes |
|--------|------|-------|
| Rostron et al. (2014) — Energy Procedia | Tier 1b | Conference proceedings, but peer-reviewed |
| White et al. (2016) — IJGG C | Tier 1 | Strong |
| White (2018) — IJGG C | Tier 1 | Strong |
| Roach et al. (2017) — Geophysics | Tier 1 | Strong |
| Haghi & Chalaturnyk (2024) — IJGG C | Tier 1 | Strong |
| Vahidinia et al. (2024) — SPE | Tier 1 | Peer-reviewed conference |
| Rangriz Shokri & Chalaturnyk (2025) — Energies | Tier 1 | MDPI, peer-reviewed |
| Peck et al. (2014) — DOE/EERC | Tier 2 | Government report, peer-adjacent |
| Zambrano-Narvaez et al. (2015) — ACG | Tier 2 | Conference paper |
| Menger (2021) — GeoConvention | Tier 2 | Conference abstract only |
| Chevrot (2025) — Thesis | Tier 2 | Primary data, but not peer-reviewed |
| Hersi & Iqbal (2025) — Eng. Proc. | Tier 2b | MDPI proceedings, limited peer review |
| Movahedzadeh et al. (2025) — SSRN | Tier 3 | Preprint server, not peer-reviewed |
| IEAGHG (2024) — Storage Catalogue | Tier 2 | Industry compilation |
| Dalkhaa (2017) — (incomplete citation) | Tier 2? | No DOI given; appears to be a report/thesis |

**Issues**:
- **Dalkhaa (2017)** is the primary source for plume diameter and breakthrough timing but has **no DOI, no URL, and no publisher** in the reference list. If this is an internal PTRC report not publicly accessible, it is effectively grey literature. **This should be resolved** — either add a full citation or note that it is a PTRC internal document.
- **Hersi & Iqbal (2025)** in Eng. Proc. is the sole source for the 21% porosity, 8.7% graywacke porosity, 78.3 mD permeability, and 974 km³ pore space. These values are central to the Black Island Member description but come from the weakest source tier that provides unique values.
- **Movahedzadeh et al. (2025) SSRN preprint** is cited but does not appear to be used for any specific table values. If not used, consider removing from the reference list to avoid the appearance of reliance.

**Overall sourcing**: Good. The core reservoir properties (thickness, porosity, permeability, pressure, temperature) are backed by Tier 1 publications. The higher and more unusual values come from lower-tier sources.

---

## 3. Cross-Validation

### 3.1 Where independent sources converge ✅

| Parameter | Source A | Source B | Match? |
|-----------|----------|----------|--------|
| Total thickness | White (2018): ~200 m | Menger (2021): 219 m mean | Good |
| Reservoir top depth | Roach (2017): ~3,130 m | White (2016): basement ~3,400 m | Consistent |
| Core plug porosity | Zambrano-Narvaez (2015): 6% | Roach (2017): 7% log | Good |
| CO₂ plume in Perf 2 | Dalkhaa (2017): 540 m (2 yr) | Roach (2017): ~200 m anomaly | Roughly consistent (different metrics) |
| Initial pressure | SPE 2024: ~35 MPa | Roach (2017): 39 MPa during inj | Consistent with expected ∆P |

### 3.2 Where independent sources diverge ⚠️

| Parameter | Source A | Source B | Gap | Resolution in notes? |
|-----------|----------|----------|-----|---------------------|
| Core k (5 mD) vs. effective k (~70 mD) | Zambrano-Narvaez | kh/h derived | 14× | ✅ Acknowledged |
| Site porosity (6–7%) vs. regional (10–15%) | Multiple site sources | Houseworth et al. | 2× | ✅ Acknowledged |
| Deadwood porosity (6%) vs. Black Island quartz arenite (21%) | Core plugs | Hersi & Iqbal | 3.5× | ❌ Not acknowledged |
| Temperature for CO₂ density: 120°C vs. 110°C | Zambrano-Narvaez (120°C) | Roach (110°C used) | 10°C | ❌ Not acknowledged |
| Icebox thickness: 15 m vs. 25–30 m | Roach (2017) | Other sources | 2× | Mentioned but unresolved |

**Key gap**: The 21% vs. 6–7% porosity divergence between the Black Island Member quartz arenite facies and the Deadwood Formation core plugs is not acknowledged in the "Data Quality Notes" section, unlike the other divergences which are explicitly discussed.

---

## 4. Area Definition — Needs Improvement

The document has an "Area" section but the definition of the **storage capacity reference area** is ambiguous:

- **3D seismic survey area**: 30 km²
- **Fine-scale model area**: ~34 km² (Peck et al., 2014)
- **Storage capacity (8.4–27.1 Mt)**: "34 km²" — same as the model area

**Questions**:

1. **Why is the model area 34 km² if the seismic survey is 30 km²?** The 4 km² difference (~13%) is unexplained. Extrapolation? Inclusion of buffer cells? This should be stated.

2. **What defines the 34 km² boundary?** Is it:
   - The area enclosed by a structural closure?
   - The simulation model boundary?
   - The area of reservoir with sufficient porosity/permeability?
   - An arbitrary modeler's cutoff?

   For a saline aquifer without a structural trap, the area used for capacity calculation is typically either the plume footprint at some future time, the "swept" area, or a regulatory boundary. None of these is stated.

3. **The 34 km² area yields 8.4–27.1 Mt storage.** If the entire Deadwood regional extent (450×500 km) were used as the area, the capacity would be orders of magnitude larger. The document correctly labels the regional numbers as gigatonne-scale, but **the factor linking 34 km² to the plume or project area is never explicitly justified**.

4. **Plume footprint after 50 yr**: ~1.5 km² (1.4 km diameter → ~1.5 km²). This is ~5% of the 34 km² model area, meaning the vast majority of the "storage capacity" is based on theoretical sweep, not demonstrated.

**Recommendation**: Add a sentence defining what the 34 km² area represents (e.g., "the areal extent of the reservoir simulation model, extending approximately 2 km beyond the seismic survey boundary to avoid boundary effects on the pressure response"). Also add a separate column or row for plume footprint area vs. model area vs. regional extent.

---

## 5. Minor Issues

| Issue | Location | Suggestion |
|-------|----------|------------|
| Dalkhaa (2017) incomplete citation | Used for plume diameter (lines 110–111), breakthrough (lines 185–187) | Add full reference or state "internal PTRC report" |
| "Verbatim: PTRC 2016 Annual Report" | Line 55, line 123 | Annual report has no DOI — state "PTRC unpublished report" |
| Roach (2015) cited for log porosity | Line 64 | Not in reference list; only Roach (2017) is listed |
| Perf 3 labelled "plugged" | Line 176 | Correct but spinner allocated ~0% — was it plugged before or after the spinner survey? Order matters. |
| 974 km³ said "subterraneous 974 km³ pore space" | Line 149 | This is both an awkward phrasing and the scale (formation-wide) should be clarified |
| BCS comparative note | Lines 209–211 | Useful but the thickness comparison (~200 m vs. ~40 m) conflates total interval at Aquistore with net pay at Quest — BCS at Quest has NTG >90% so net pay ~36 m vs. Aquistore net ~112 m; the comparison would be clearer as net-to-net |

---

## Summary

| Category | Rating | Action required |
|----------|--------|----------------|
| Reservoir thickness | ✅ Pass | None |
| Porosity | ⚠️ Minor | Flag 21% as regional, not site-specific; add spatial context for Hersi & Iqbal values |
| Permeability | ✅ Pass | Consider clarifying h used for kh→k conversion |
| Temperature/Pressure | ⚠️ Minor | Reconcile 110°C vs. 120°C in CO₂ density calculation |
| Storage capacity | ⚠️ Minor | Define what 34 km² area represents; distinguish site capacity from regional resource |
| Seal | ✅ Pass | Consider tightening Icebox thickness range |
| Area definition | ❌ Needs work | Ambiguous what 34 km² represents; add plume-scale area |
| Source citation completeness | ⚠️ Minor | Find full citation for Dalkhaa (2017); verify Roach (2015) |
| Cross-validation gaps acknowledged | ⚠️ Minor | Missing acknowledgment of Hersi & Iqbal vs. core porosity discrepancy; missing note on T inconsistency |

**Overall**: High-quality compilation. The core Tier-1 sourced properties (thickness, Deadwood core/log porosity, pressure, temperature) are well-supported. The weak points cluster around (a) the CO₂ density temperature mismatch, (b) undefined storage capacity area, and (c) the Hersi & Iqbal (2025) values needing clear spatial context.
