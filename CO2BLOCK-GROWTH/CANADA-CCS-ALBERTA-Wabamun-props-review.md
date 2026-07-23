# Adversarial Review: CANADA-CCS-ALBERTA-Wabamun-props.md

**Reviewer**: Automated adversarial review  
**Date**: 2026-07-23  
**Scope**: Reservoir property reasonableness, source tiers, area definition, discrepancies

---

## 1. Reservoir Property Review

### 1.1 Thickness (BSU — primary target)

| Issue | Severity | Finding |
|-------|----------|---------|
| Range width | **Medium** | Sources span ~33 m (Hofmann 2013) to 40–80 m (AGS PRS / IEAGHG 2024-TR03). The file converges on 30–50 m for Wabamun, then recommends ~40 m. This is reasonable as a midpoint, but the 40–80 m AGS range is cited without verifying whether that range is for the BSU specifically or the broader BCS across all of Alberta. If the latter, the upper bound (80 m) may not apply to the Wabamun area. |
| Hofmann 2013 vs. Quest analog | **Low** | The 33 m (Hofmann) and 35–46 m (Quest) are from different formations/areas at different depths. Both are presented with appropriate caveats. |
| No Wabamun-specific data | **High** | Two stratigraphic test wells drilled (2023–2024) have unpublished results. All thickness values are regional or analog. The file is transparent about this. |

**Verdict: Fair. Recommended ~40 m is a justifiable midpoint of the available regional data.**

### 1.2 Porosity (BSU)

| Issue | Severity | Finding |
|-------|----------|---------|
| Weides (2014) ~10% vs. Alberta Innovates ~14.4% | **Low** | The file correctly notes that Alberta Innovates covers a broader area including shallower, higher-porosity portions. The 10% value from Weides (6 wells, comparable depths) is more applicable to Wabamun. |
| Quest analog (16–17%) inapplicability | **Low** | The file correctly argues that Quest data at ~2,000 m depth is not directly applicable to Wabamun at 2,300–3,300 m due to compaction. The "5–10% porosity loss" estimate from an additional ~500–1,300 m of burial is qualitatively correct but not quantitatively sourced. No reference is given for this specific compaction trend. |
| Recommended range 8–12% | **Low–Medium** | The central value of 10% is well-supported by Weides (2014). However, the lower bound of 8% is not directly sourced — it is inferred from the depth correction argument. This is reasonable but should be noted as an inference. The upper bound of 12% is also an interpolation, not a measured value. |

**Verdict: Good. The 10% central value is the best available estimate. The 8–12% range is a reasonable inference but not directly measured.**

### 1.3 Permeability (BSU) — **Highest uncertainty parameter**

| Issue | Severity | Finding |
|-------|----------|---------|
| Plug-to-regional scale gap | **High** | The gap between core-plug permeability (<0.01 mD, Weides 2014, 96 samples) and regional-scale permeability (1.6 mD, Hofmann 2013) is a factor of **>160×**. This is an enormous range for a reservoir parameter and is poorly explained mechanistically. The file states the scale-up "accounts for fractures and larger-scale heterogeneity" but does not quantify the fracture contribution or cite evidence for a connected fracture network in the BSU at Wabamun depths. |
| Recommended "1 mD as representative pre-stimulation value" | **Medium–High** | The recommendation of 1 mD (midpoint of 0.01–1.6 mD range) is not well-justified. Why 1 mD rather than 0.1 mD or 0.01 mD? If the BSU has no significant natural fractures, the pre-stimulation permeability would be closer to 0.01 mD. If it does, it could approach 1.6 mD. The file does not resolve this ambiguity. Recommend presenting both end-members with explicit scenarios rather than a single "representative" value. |
| Quest analog permeability (20–1,000 mD) | **Low** | Correctly flagged as inapplicable and overly optimistic for Wabamun. |
| Stimulation requirement | **Low** | Well-supported by multiple sources (Hofmann 2013, Weides 2014) that hydraulic stimulation is required. The note about the BSU being "well-confined" is supported. |

**Verdict: Marginal. The >160× range between measurement scales is insufficiently explained. The recommended 1 mD value lacks clear justification. Recommend scenario-based presentation (low/mid/high).**

### 1.4 Geomechanical Properties

| Issue | Severity | Finding |
|-------|----------|---------|
| UCS up to 97.7 MPa | **Medium** | This is very high for a sandstone. Typical sandstone UCS ranges 30–80 MPa. 97.7 MPa is at the upper bound of plausible for a well-cemented quartz arenite. The file does not note that this is unusually high. Recommend caveat that this may represent the strongest samples, not the average. |
| Cohesion up to 69.8 MPa | **High** | This is **extremely high** for a sandstone. Typical sandstone cohesion is 5–40 MPa. 69.8 MPa is more characteristic of a low-porosity carbonate or crystalline rock than a sandstone. This value may be suspect, or it may reflect a specific well-cemented sample. The file flags the friction coefficient as "remarkably high" but does NOT flag the cohesion similarly. This is an inconsistency in the file's own quality notes. |
| Friction coefficient up to 1.22 | **Medium** | Corresponds to a friction angle of ~50.6°, which is unusually high for sandstone (typical 25–45°, μ = 0.47–1.0). The file correctly notes this is "remarkably high" but should also note these values are likely from select clean quartz arenite samples and may not represent the formation average. |
| All geomechanical values from Weides (2014) | **Medium** | These are from a single study focused on geothermal exploration, not CCS. The number of samples and whether outliers were excluded is not reported. |

**Verdict: Marginal. The geomechanical properties (especially cohesion at 69.8 MPa and friction coefficient at 1.22) are extreme values that warrant stronger caveats than currently present. Recommend flagging cohesion similarly to the friction coefficient.**

### 1.5 Other Parameters

| Parameter | Finding |
|-----------|---------|
| Depth (BSU) 2,300–3,300 m | **Good** — consistent across independent sources. |
| Depth (Nisku) 1,550–2,200 m | **Good** — single source but consistent with regional geology. |
| Temperature 80–115 °C | **Fair** — derived from regional gradient (35.6 °C/km). Reasonable but not Wabamun-specific. |
| Reservoir pressure 23–33 MPa (hydrostatic) | **Fair** — hydrostatic assumption is standard in absence of DST data, but overpressure is common in deep Alberta formations. A brief caveat about possible overpressure would strengthen the file. |
| CO₂ density 600–750 kg/m³ | **Good** — consistent with expected P/T conditions. |
| Seal properties from Quest analog | **Fair** — same regional seals should be present, but thickness and quality at Wabamun may differ. Appropriately qualified. |
| Net-to-gross 0.7–0.9 | **Medium** — unsupported wide range. Cited only as "Weides et al. (2014) — BSU has minor shale intercalations. Not quantified for Wabamun specifically." The file is transparent but 0.7–0.9 is a wide bracket. Recommend a central estimate with justification. |

---

## 2. Source Tiers

The file uses a clear three-tier system:

| Tier | Description | Sources | Assessment |
|------|-------------|---------|------------|
| 1 (Project-specific) | Lease area, target formations, depth, project status, CO₂ volumes | carbonstorage.io, Government of Alberta, Enbridge press releases, CBC, NRCan | **Good** — 18 total sources; well-referenced with verbatim quotes |
| 2 (Regional) | BSU thickness, porosity, permeability from published peer-reviewed studies in central Alberta | Weides et al. (2014), Hofmann et al. (2013), Alberta Innovates | **Good** — peer-reviewed, properly cited |
| 3 (Analog) | Quest CCS, WASP study, Project Pioneer | Shell D-65/ERCB, IEAGHG 2024-TR03, Keith & Lavoie (2009) | **Good** — correctly flagged as analog with limitations noted |

**Strengths:**
- Every value is traceable to a specific source with a verbatim quote
- Source types are clearly labeled (peer-reviewed journal, government registry, industry database, news article)
- Quest analog limitations are thoroughly documented in a dedicated comparison table
- A data quality summary table at the end provides a quick-reference quality assessment

**Weaknesses:**
- Source #13 (keyfactsenergy.com) and #14 (Chinook Petroleum blog) are low-tier industry sources (websites/blogs). Their overlap with source #1 (carbonstorage.io) is noted but not independently verified.
- Source #18 (Transition Accelerator Report 4.6) is a research report, not peer-reviewed. Its value for reservoir properties is limited (cites PCOR Atlas data).

---

## 3. Area Definition

| Issue | Finding |
|-------|---------|
| Multiple converging sources | **Good** — carbonstorage.io (~220,000 ha), keyfactsenergy.com (211,072 ha / 23 townships), and Parkland County location all converge on ~211–220k ha. |
| Discrepancy: 220,000 ha vs. 734,192 acres | **High** — These are NOT equivalent. 734,192 acres = ~297,100 ha, which is **35% larger** than 220,000 ha. The file speculates this "may represent the full sequestration lease area (post-October 2025 agreement)" but provides no evidence for this claim. This discrepancy needs resolution: either 734,192 acres is a different metric (e.g., pore space volume expressed as equivalent surface area) or it includes additional tenure. |
| DLS reference | **Good** — Well locations (15-32-050-03 W5M, 02-18-055-24 W4M, 08-19-055-22 W4M) provide concrete spatial constraints. The range of Townships 49–53, Ranges 1–6 W5M is a reasonable estimate. |
| WASP study area (5,000 km²) | **Low** — Correctly noted as a precursor study area, not the current lease. |

**Verdict: Fair. The primary values converge, but the 734,192 acre discrepancy is unresolved.** 

---

## 4. Discrepancies Log

| # | Discrepancy | Parameters | Impact | Resolution in File |
|---|-------------|------------|--------|-------------------|
| 1 | **Area**: 220,000 ha vs. 734,192 acres (~297,100 ha) | Lease area | **High** → Affects areal extent for capacity calculations | Speculative (evaluation vs. sequestration agreement); no evidence cited |
| 2 | **Permeability scale gap**: <0.01 mD (core plug) vs. 1.6 mD (regional) | BSU permeability | **High** → Factor of 160× uncertainty affects injectivity modeling | Noted but not mechanistically explained |
| 3 | **Porosity**: 10% (Weides) vs. 14.4% (Alberta Innovates) | BSU porosity | **Medium** → 44% relative difference | Explained (different study areas/scope) |
| 4 | **Thickness**: 33 m (Hofmann) vs. 40–80 m (AGS/IEAGHG) | BSU thickness | **Medium** → Factor of 2× at upper bound | Partially explained (regional vs. local); the 40–80 m range origin needs verification |
| 5 | **Quest permeability**: 20–500 mD (Shell) vs. 33–1,000 mD (IEAGHG) | BCS (analog) | **Low** — not directly applied | Both presented but correctly flagged as inapplicable to Wabamun |
| 6 | **Geomechanical extremes**: Cohesion 69.8 MPa (unusual for sandstone), μ = 1.22 (unusual for sandstone) | UCS, cohesion, friction coefficient | **Medium** → May mislead geomechanical modeling if taken as representative | Friction coefficient flagged as "remarkably high"; cohesion NOT flagged |
| 7 | **Temperature upper bound**: 120 °C (regional) vs. 115 °C (expected) | Formation temperature | **Low** → Minor 5 °C difference at upper end | The file uses 80–115 °C for Wabamun, narrower than 65–120 °C regional range. Acceptable. |
| 8 | **Depth-dependent porosity loss**: "5–10% loss" cited but not quantified | Porosity depth correction | **Medium** → No supporting reference for the specific reduction factor | The file states "an additional ~5–10% porosity loss is expected" without citing a compaction trend equation (e.g., Athy's law). This is a qualitative estimate. |

---

## 5. Summary and Recommendations

### What works well

1. **Source traceability**: Excellent — every value links to a specific source with verbatim quotes. Best practice for CO2BLOCK.
2. **Honest uncertainty communication**: The file consistently labels regional vs. Wabamun-specific data.
3. **Quest analog limitations**: Thoroughly documented with a comparison table. Correctly advises against direct application.
4. **Data quality summary table**: Provides quick reference for all key parameters.

### What needs attention

| Priority | Issue | Recommended Action |
|----------|-------|-------------------|
| **High** | Area discrepancy (220,000 ha vs. 734,192 acres / ~297,100 ha) | Resolve: Check the carbonstorage.io source directly. If 734,192 is correct, the evaluation and sequestration areas differ by ~35%. If the file's speculation is wrong, correct the values. |
| **High** | Permeability recommended value (1 mD) is poorly justified | Replace with scenario-based presentation: low (0.01 mD — matrix), mid (0.1–1 mD — with natural fractures), high (1–10 mD — stimulated). Or add explicit justification for the 1 mD choice. |
| **High** | Cohesion 69.8 MPa is an extreme value not flagged | Add caveat that this is an unusually high value typical of low-porosity carbonates/crystalline rocks, not typical sandstone, and may represent select samples. |
| **Medium** | "5–10% porosity loss" with depth lacks a reference | Add a citation for depth-porosity trends in Cambrian sandstones (e.g., Athy's law parameters from a published Alberta Basin compaction study). |
| **Medium** | Overpressure not considered for reservoir pressure estimate | Add a caveat that deep Alberta formations can be overpressured; hydrostatic assumption is a default but may be non-conservative. |
| **Medium** | Net-to-gross 0.7–0.9 is unsupported | Narrow the range or cite specific well-log data. If no data exist, state this explicitly. |
| **Low** | Hofmann (2013) temperature gradient from a geothermal study may not be ideal for CCS | Cross-check against AER Bottom-Hole Temperature database for Wabamun-area wells. Include AER BHT reference if available. |

### Overall Assessment

**The file is well-researched, transparent about uncertainties, and properly attributes all values. It is usable for CO2BLOCK with the following caveats:**

1. **Permeability and area are the two most critical unresolved issues.** The permeability uncertainty (>160× range) and the area discrepancy (35% difference) materially affect any storage capacity or injectivity calculation.
2. **Geomechanical properties from Weides (2014) should be applied with caution.** The reported maximum values (cohesion 69.8 MPa, μ = 1.22) are extreme for sandstone and may not be representative.
3. **All thickness, porosity, and permeability values remain regional estimates** until Enbridge publishes test well results. This is appropriately stated in the file.
4. **Recommended for CO2BLOCK use** with the three high-priority corrections above.
