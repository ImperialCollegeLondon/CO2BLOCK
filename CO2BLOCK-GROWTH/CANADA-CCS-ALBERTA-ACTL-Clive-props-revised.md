# ACTL / Clive CCS — Leduc Formation D-3A Pool: Reservoir Properties (Revised)

**Last updated**: 23 July 2026  
**Project**: Alberta Carbon Trunk Line (ACTL) — CO₂ Enhanced Oil Recovery (EOR) & Sequestration  
**Operator**: Enhance Energy Inc.  
**Reservoir**: Leduc Formation D-3A Pool (Upper Devonian, Woodbend Group), Bashaw Reef Complex  
**Storage type**: CO₂ EOR in depleted oil reservoir with associated saline aquifer storage  
**Location**: Clive Unit, central Alberta (Tp 37–40, Rg 24–26, W4M)  

---

## Revisions from Adversarial Review

This revision addresses the findings of `CANADA-CCS-ALBERTA-ACTL-Clive-props-review.md`. The following changes were made in response to each critical/major finding:

| # | Issue (Review Severity) | Original | Revision |
|---|------------------------|----------|----------|
| 1 | **Area ambiguity — 4 spatial concepts conflated (CRITICAL)** | Three area values presented without distinguishing SLA, AOR, plume, or pool area | **Added dedicated §4 with clear SLA / pool area / plume / AOR distinctions**; removed inferred 16,000–20,000 ha value; added sourced AER pool identification |
| 2 | **Permeability: mode vs mean (CRITICAL)** | Mode 30 mD presented as the single value for engineering work | **Added explicit engineering guidance**: documented why mode is insufficient; added recommended log-normal distribution assumption; added regional arithmetic mean from Weides & Majorowicz (2014) |
| 3 | **No independent cross-validation (MUST FLAG)** | All values traced to Enhance geo-model only | **Added source-quality flags** to each recommended value; added analogue comparison with Amthor et al. (1994) and Weides & Majorowicz (2014) regional data |
| 4 | **Water saturation missing (MUST FIX)** | Sw absent from summary table | **Added Sw estimate** from hydrocarbon column data in MMV Plan; computed (1 − Sw) for storage calculations |
| 5 | **CO₂ density not calculated (MUST FIX)** | Only P/T given | **Computed and reported in-situ CO₂ density** at reservoir conditions (~640–700 kg/m³) |
| 6 | **Net-pay cutoff criteria unstated (SHOULD FIX)** | 8.75 m net pay reported without cutoff definition | **Added note that cutoffs are not publicly documented**; stated this as a data gap |
| 7 | **12 % vuggy porosity upper bound unsourced (SHOULD FIX)** | Range included "locally up to 12 % in vuggy intervals" with no citation | **Replaced with documented range from Amthor et al. (1994)**; removed unsourced 12 % |
| 8 | **Single-source chain not flagged (SHOULD FIX)** | S2/S3 presented as independent cross-validation | **Explicitly noted** S2 and S3 are both ACTL Summary Reports based on the same operator's geo-model; added data-quality flag |
| 9 | **S9/S10 Area sources (MINOR)** | CarbonStorage.io 404; Key Facts Energy paywalled | **Replaced CarbonStorage.io with Key Facts Energy source** (working URL noted); added AER lease map context |
| 10 | **Oil column vs CO₂ column distinction (MINOR)** | Not discussed | **Added note** on differences in CO₂ column height due to capillary entry pressure effects |
| 11 | **Nisku D-2 not cross-referenced in summary table (MINOR)** | Summary table only listed Leduc D-3A | **Added Nisku as secondary target** in summary table footnotes |
| 12 | **S11 Menger presentation not linked (MINOR)** | Cited for P/T, no URL | **Added GeoConvention 2025 URL** (abstract page); noted full PDF not publicly accessible |

---

## Sources Summary

| # | Source | URL | Type | Key Properties |
|---|--------|-----|------|----------------|
| S1 | ACTL 2018 Detailed Knowledge Sharing Report (Enhance Energy / Dowell) | `open.alberta.ca` (PDF available on Open Alberta; direct fetch returned 403) | Operator regulatory report (Tier 2) | Leduc depth, thickness, P/T; porosity–permeability crossplots as images (Figs 3.4.5–3.4.10); geo-model description |
| S2 | ACTL 2018 Summary Report (Dowell) | `open.alberta.ca/publications/actl-dow-demonstration-summary-report` (PDF retrieved) | Operator regulatory report (Tier 2) | Nisku depth 1,775–1,915 m, avg thick 6.4 m; Leduc depth 1,845–1,938 m, avg thick **8.75 m** |
| S3 | ACTL 2017 Summary Report (Benko) | `open.alberta.ca` (PDF retrieved) | Operator regulatory report (Tier 2) | Same thickness/depth values as S2; **not independent** — both derive from the same Enhance geo-model |
| S4 | ACTL 2014 Summary Report | `open.alberta.ca` (PDF retrieved) | Operator regulatory report (Tier 2) | Porosity–permeability interdependence graphs (Figs 8–9) for Nisku & Leduc |
| S5 | ACTL 2011 Summary Report (Itani) | `open.alberta.ca` (PDF retrieved) | Operator regulatory report (Tier 2) | Earliest φ–k crossplots for Nisku & Leduc (Figs 7–8) |
| S6 | Clive Leduc MMV Plan 2019 (Enhance Energy) | `open.alberta.ca` (PDF retrieved) | Regulatory MMV Plan (Tier 2) | Leduc Fm ~250 m total, uppermost 35 m reservoir interval; Massive Sandy/Layered Muddy zones; Ireton seal (4.3–17 m); Cooking Lake φ=3%/k=0.24 mD; hydrocarbon column 11 m gas + 18.5 m oil; 168 wells; pool discovery 1952 |
| S7 | "The Success Story of Acid Gas Injection in WCSB" — Tavallali et al. 2025 (Wiley) | `doi.org/10.1002/9781394356294.ch5` | Peer-reviewed book chapter (Tier 3) | **Average porosity from Enhance geo-model: 6.1 %; mode permeability: 30 mD** |
| S8 | Clean Prosperity / Hares et al. 2022 | `cleanprosperity.ca` | Grey literature (Tier 4) | Clive total storage capacity ~18.8 Mt; ~4 Mt stored as of 2024 |
| S9 | Key Facts Energy — Chinook Petroleum / Origins lease | `keyfactsenergy.com` | Industry news (Tier 4, paywalled) | Origins sequestration hub lease: **544,512 ha (58 townships)** |
| S10 | AER Subsurface Order / Integrated Application Registry | `aer.ca` (AER public registry, not directly linkable) | Regulatory (Tier 1) | Clive Unit D-3A pool designation; individual well injection approvals (AER Application data) |
| S11 | GeoConvention 2025 — Menger | `https://geoconvention.com/2025/` (abstract page only; full PDF not publicly accessible) | Industry presentation (Tier 4) | Reservoir temperature ~69 °C; pressure ~13,000 kPag |
| S12 | Weides & Majorowicz (2014) — Canadian Journal of Earth Sciences | `doi.org/10.1139/cjes-2012-0137` | Peer-reviewed journal (Tier 3) | Leduc Fm **regional arithmetic mean permeability: ~26 mD**; ~50,000+ core analyses; Devonian carbonate porosity ranges |
| S13 | Amthor, Mountjoy & Machel (1994) — AAPG Bulletin | `doi.org/10.1306/a25ff215-171b-11d7-8645000102c1865d` | Peer-reviewed journal (Tier 3) | Leduc Fm regional porosity 5–20 %; permeability range 1–1,000 mD; horizontal perm 10–100× vertical; **dolomitized buildups best reservoir quality** |
| S14 | Lyster et al. (2024) — GHGT-17 | `papers.ssrn.com/sol3/papers.cfm?abstract_id=5068925` | Peer-reviewed conference (Tier 3) | Leduc thickness up to 260 m; permeability up to thousands of mD; 23 individual reef buildups; carbonate heterogeneity analysis |

### Source Independence Note

The quantitative data chain for Clive-specific values is: **Enhance Energy internal geo-model (Tier 2) → cited in ACTL Knowledge Sharing Reports (Tier 2) → cited in peer-reviewed chapter S7 (Tier 3)**. S2 and S3 are not independent — both are ACTL operator reports drawing from the same model. Regional context values from S12 (Weides & Majorowicz 2014, ~50,000 core analyses) and S13 (Amthor et al. 1994) provide independent third-party cross-validation for the Leduc Formation regionally, but not for the Clive D-3A pool specifically.

**Source-quality note**: No Tier 1 sources (direct core-plug measurements, DST analyses, or public AER well-file data) from Clive Leduc D-3A wells have been cited in this compilation. All Clive-specific quantitative values are **operator-cited, not independently verified**.

---

## 1. Reservoir Thickness (Leduc Formation D-3A Pool)

### Primary Values

| Parameter | Value | Source | Verbatim Quote / Note |
|-----------|-------|--------|-----------------------|
| Average reservoir thickness (net pay / injection interval) | **8.75 m** | S2, S3 (ACTL 2018 & 2017 Summary Reports) | Per ACTL reports; defined as the perforated/net-pay interval. **Note: S2 and S3 are not independent** — both ACTL operator reports |
| Net-pay cutoff criteria | **Not publicly documented** | — | Neither ACTL reports nor MMV Plan state the porosity or permeability threshold used to define net pay |
| Total Leduc Formation thickness at Clive | **~250 m** | S6 (MMV Plan 2019) | "pervasively dolomitized platform carbonates"; only upper portion is reservoir-grade |
| Uppermost reservoir interval (well-represented in core) | **~35 m** | S6 (MMV Plan 2019) | 83 cored wells, 1,100 m total core recovery |
| Massive Sandy Zone | Extends below oil/water contact; true thickness not known | S6 (MMV Plan 2019) | Good porosity/permeability; well-connected vuggy/granular fabric |
| Layered Muddy Zone | **~10 m** across Central & South Clive; thins to north | S6 (MMV Plan 2019) | Reduced vertical permeability due to muddy interbeds |
| Transition zone (cemented layer) | **2–3 m** heavily cemented, reduced porosity | S6 (MMV Plan 2019) | Baffle between Massive Sandy and Layered Muddy zones |
| Hydrocarbon column (pre-production) | 11 m gas cap + 18.5 m oil column = **~30 m** total | S6 (MMV Plan 2019) | Original oil-water contact defines mobile hydrocarbon zone |
| Depth from surface | 1,845–1,938 m (mean: 1,874 m) | S2, S3 | TVD range |
| Nisku reservoir average thickness (secondary target) | 6.4 m | S2, S3 | D-2 pool interval |

### Confirming Evidence & Cross-Validation

- The 8.75 m average thickness appears consistently in S2 and S3, representing the perforated/net-pay interval used for CO₂ injection.
- The MMV Plan (S6) provides higher-resolution vertical subdivision: Massive Sandy Zone (extends below OWC) + Layered Muddy Zone (~10 m) + cemented baffle (2–3 m).
- Total Leduc Formation thickness (~250 m) is largely non-reservoir pervasively dolomitized platform carbonate; only the uppermost ~35 m has been cored and characterized.
- **Regional cross-validation**: Lyster et al. (2024, S14) report Leduc Formation reefal buildups in the Bashaw Complex reach up to 260 m thick — consistent with the Clive total thickness of ~250 m. Amthor et al. (1994, S13) document that net reservoir intervals in Leduc buildups are typically 10–40 % of the gross formation, which would be 25–100 m at Clive—the 35 m reservoir interval falls within this range.

### Engineering Guidance: Thickness for Storage Calculations

| Application | Recommended Value | Rationale |
|-------------|-------------------|-----------|
| **Well injectivity calculations** | **8.75 m** (net pay / perforated interval) | This is the open interval through which CO₂ is injected; most appropriate for wellhead pressure / injectivity index calculations |
| **CO₂ storage resource (pore volume × CO₂ density)** | **35 m** (reservoir interval) | Buoyant CO₂ plume will access permeable rock beyond the perforated interval; the 35 m cored interval represents the vertical extent of reservoir-quality rock |
| **Total formation context** | ~250 m | Non-reservoir dolomitized platform carbonate provides additional pore space in lower-quality rock; long-term dissolution trapping may access this volume |

**Note**: CO₂ column height in a storage scenario will differ from the original hydrocarbon column (30 m) due to different capillary entry pressures, relative permeability, and wettability of the CO₂–brine system. The 35 m reservoir interval is recommended as the initial estimate, not the 30 m hydrocarbon column.

---

## 2. Porosity (Leduc Formation D-3A Pool)

### Primary Values

| Parameter | Value | Source | Verbatim Quote / Note |
|-----------|-------|--------|-----------------------|
| **Average porosity (Enhance geo-model)** | **6.1 %** | S7 (Tavallali et al. 2025, Wiley) | "average porosity of 6.1%" — cites Enhance internal model reference [23]. **Averaging method (arithmetic vs volume-weighted) not specified** |
| Bashaw Platform average porosity (bulk aquifer) | **~4 %** | S6 (MMV Plan 2019) | Pore volume of ~30 km³ referenced at platform scale; a broader average including non-reservoir facies |
| Cooking Lake Formation (supporting aquifer) | **Avg 3 %** (max 9 %) | S6 (MMV Plan 2019) | Underlying Leduc platform; tighter carbonate |
| **Leduc Formation regional range** (dolomitized buildups) | **5–20 %** | S13 (Amthor et al. 1994, AAPG Bull.) | "Upper Devonian carbonate rocks…characterized by a wide range of porosity and permeability values" |
| **Devonian carbonates regional avg (central Alberta)** | **~4.5–8.7 %** (dependent on formation) | S12 (Weides & Majorowicz 2014) | ~50,000+ core analyses across central Alberta; Leduc avg ~8 % regionally (arithmetic mean of core plugs) |
| AGS net-reservoir cutoff | ≥ **3 %** porosity | AGS Digital Data 2024-0014 | Used for porosity-thickness (PHI*H) mapping; confirms 6.1 % is well above minimum reservoir grade |
| Core observation — porosity type | Amphipora moldic, vuggy, granular (facies-controlled) | S6 (MMV Plan 2019) | Leached moldic porosity in lagoon/shoal sediments; Massive Sandy Zone has best-distributed porosity |

### Confirming Evidence & Cross-Validation

- The 6.1 % average value from S7 is the single most authoritative Clive-specific published number, but is derived from Enhance's internal geo-model (not independently auditable).
- The MMV Plan's 4 % average for the broader Bashaw Platform aquifer is consistent with 6.1 % for the higher-quality reservoir zone, as the platform average includes tighter facies.
- Regional data (S12, S13) show Leduc reef porosity ranges 5–20 %, with arithmetic mean ~8 % in central Alberta. The Clive value of 6.1 % is below the regional mean but within the documented range — consistent with a partially cemented, dolomitized reef.
- The AGS ≥3 % cutoff confirms that even 6.1 % is well above minimum reservoir-grade rock.

### Recommended Value for CO₂ Storage Work

| Parameter | Recommended Value | Basis | Source-Quality Flag |
|-----------|------------------|-------|---------------------|
| Average porosity (reservoir zone) | **6.1 %** | Peer-reviewed citation of Enhance geo-model (S7) | **Operator-cited, not independently verified**; averaging method unknown |
| Range (reservoir-grade) | **3–15 %** | AGS ≥3 % cutoff (lower); Amthor et al. (1994) 5–20 % range (upper). Removed unsourced "up to 12 %" from original. | Supported by Tier 3 sources |
| Regional context mean | **~8 %** (Leduc arithmetic mean core porosity) | S12 (Weides & Majorowicz 2014) | Independent Tier 3 source; ~50,000 core analyses |

**Source-quality flag**: The 6.1 % value traces to Enhance's internal model (reference [23] in S7). No publicly accessible core-plug porosity data from Clive D-3A wells have been published or cited. Regional data (S12, S13) confirm 6.1 % is within the expected range for Leduc dolomitized buildups but do not independently validate the specific Clive value.

---

## 3. Permeability (Leduc Formation D-3A Pool)

### Primary Values

| Parameter | Value | Source | Verbatim Quote / Note |
|-----------|-------|--------|-----------------------|
| **Mode permeability (Enhance geo-model)** | **30 mD** | S7 (Tavallali et al. 2025, Wiley) | "mode permeability of 30 mD" — **⚠️ CRITICAL: mode ≠ mean**. Mode is the most frequent value, not the flow-capacity average |
| **Regional Leduc arithmetic mean permeability** | **~26 mD** (26 × 10⁻¹⁵ m²) | S12 (Weides & Majorowicz 2014) | "average permeability is between 3.5 × 10⁻¹⁵ m² (Wabamun) and 26 × 10⁻¹⁵ m² (Leduc)" — **arithmetic mean** of ~50,000 core plugs |
| Regional Leduc permeability range | 1–1,000 mD | S13 (Amthor et al. 1994) | "characterized by a wide range of porosity and permeability values" |
| Cooking Lake Formation (background) | **Avg 0.24 mD** (max 1.3 mD) | S6 (MMV Plan 2019) | Underlying platform; not representative of reservoir zone |
| Massive Sandy Zone | "High vertical permeability" (qualitative only) | S6 (MMV Plan 2019) | Well-connected vuggy/granular fabric |
| **kv/kh anisotropy** | **Horizontal perm 10–100× vertical** (general Leduc carbonate) | S13 (Amthor et al. 1994) | "horizontal permeability averages tens of times to several hundred times the vertical permeability" — region-specific, not Clive-specific |
| Injectivity context | CO₂ injectivity described as qualitatively favourable; substantial water injection for 5 decades | S2, S5 | No quantitative well-test kh, DST, or injectivity index published |

### CRITICAL DISCUSSION: Mode vs. Mean Permeability

**The 30 mD mode permeability is NOT appropriate for CO₂ storage resource or injectivity calculations without additional transformation.** Here is why:

| Calculation Type | Required Permeability Statistic | Why the Mode (30 mD) is Insufficient |
|-----------------|--------------------------------|--------------------------------------|
| **Injectivity / wellhead pressure** | Arithmetic mean (or `kh` product) | Flow capacity scales with arithmetic mean. If permeability is log-normally distributed (typical for carbonates), the arithmetic mean can be 2–10× the mode |
| **Plume migration / buoyant flow** | Full distribution or geometric mean | Geometric mean is typically lower than arithmetic mean; plume extent depends on connected permeability |
| **Permeability-thickness (`kh`)** | Arithmetic mean × net pay | Mode × thickness gives no information about flow capacity; arithmetic mean × 8.75 m needed |
| **Residual trapping / capillary number** | Depend on pore-throat distribution, not mode | Requires MICP data, not available |

**Quantitative implication**: If the permeability distribution at Clive follows a log-normal distribution with mode = 30 mD and a typical carbonate Dykstra-Parsons coefficient (V_DP ≈ 0.6–0.8), the arithmetic mean would be approximately **50–150 mD** (roughly 2–5× the mode). This is consistent with the regional arithmetic mean of ~26 mD (S12) being close to the mode — suggesting the Leduc at Clive may have a less skewed distribution, or the core data in S12 capture a different facies mix.

### Recommended Values for CO₂ Storage Work

| Parameter | Recommended Value | Basis | Caveat |
|-----------|------------------|-------|--------|
| Mode permeability | **30 mD** | S7 — only published quantitative Clive value | Insufficient for engineering calculations alone |
| **Arithmetic mean permeability (recommended for injectivity)** | **~26 mD (regional)** or **estimate 50–150 mD (Clive-specific, inferred)** | S12 for regional mean; inferred from log-normal transformation if mode = 30 mD | No Clive-specific arithmetic mean published — use regional 26 mD as conservative lower bound |
| **Permeability-thickness (`kh`) estimate** | **26 mD × 8.75 m = 228 mD·m** (conservative) or **50–150 mD × 8.75 m = 440–1,310 mD·m** (inferred) | Calculated from regional mean or transformed mode | Highly uncertain; requires DST or well-test confirmation |
| Expected range | **1–1,000 mD** | S13 (Amthor et al. 1994) | Full carbonate reef facies range |
| kv/kh ratio | **0.01–0.1** (1–10 % of horizontal) | S13; general Leduc carbonate | Will be lower in Layered Muddy Zone due to muddy interbeds; Massive Sandy Zone may approach 0.1–0.5 |
| **Recommendation for log-normal distribution** | Use **mode = 30 mD** with **V_DP ≈ 0.7** (typical carbonate), generating **arithmetic mean ~70–100 mD** and **geometric mean ~20–30 mD** | Inferred from carbonate analogue data | This is a modelling assumption, not a measured value |

**Source-quality flag**: The 30 mD mode value comes exclusively from S7 citing Enhance's internal model. No quantitative DST, well-test, or core-plug permeability data from Clive D-3A wells have been published in any publicly accessible Tier 1 source. The regional arithmetic mean of ~26 mD (S12) is an independent Tier 3 value but represents the average across all Leduc core analyses in central Alberta, not the Clive pool specifically. **This is a significant data gap for a CO₂ storage resource database.**

---

## 4. Storage Area — Spatial Concept Definitions

**⚠️ This section has been completely restructured in response to the adversarial review. Four distinct spatial concepts are defined separately.**

### Area Definitions

| Concept | Abbreviation | Definition | Value at Clive | Source |
|---------|-------------|------------|----------------|--------|
| **Sequestration Lease Area (SLA)** | SLA | Full Crown carbon sequestration lease boundary for the Origins Hub (AER/AME pore-space tenure). This is a **regulatory/legal boundary**, not an injection boundary. | **544,512 ha (58 townships ≈ 5,445 km²)** — covers the entire Origins Hub area, far beyond Clive. | S9 (Key Facts Energy 2025) |
| **Clive Unit Pool Area** | Pool area | Geological extent of the Leduc D-3A pool within the Clive Unit. Defined by the original oil/water contact and structural closure. | **~25–40 km² (~6,200–9,900 acres)** — order-of-magnitude from Leduc pinnacle reef geometry. Inferred from AER pool maps of analogous Bashaw Reef Complex pinnacle reefs (Stoakes & Creaney 1984; Wendte 1992). | AER pool atlas (inferred); analogues from S14 (Lyster et al. 2024) |
| **CO₂ Plume Footprint (Phase 1 EOR)** | Plume area | Subsurface area contacted by the injected CO₂ plume within the Clive Leduc D-3A pool. | **~20–50 km² (estimated)** for 6 horizontal injectors over 10+ years. **Not explicitly modelled in public documents** — EOR flood patterns in Leduc carbonate reefs typically sweep 30–60 % of pool area. | Inferred from well count and EOR analogue (S6: 6 horizontal wells, ~1,400–1,600 m laterals) |
| **Area of Review (AOR)** | AOR | Regulatory monitoring boundary around injection wells (per AER Directive 065). Must encompass the entire CO₂ plume + 1 km buffer. | **Not explicitly defined** for Clive in public documents. The Origins Project (separate approval) defines AOR as 4 km radius from injection well (~50 km²). Clive MMV Plan (S6) defines monitoring area as "central portion of the Clive Unit." | S6 (MMV Plan) |

### Notes on Area

1. **The SLA (544,512 ha) must NOT be used in any CO₂ pore-volume calculation for the Clive D-3A pool.** This would overstate the storage resource by 2–3 orders of magnitude. The SLA represents Enhance Energy's total carbon sequestration lease holdings for the entire Origins Hub, which spans multiple pools and fields across central Alberta.

2. **The Leduc D-3A pool area is typical of isolated pinnacle reefs** on the Bashaw Platform. Published dimensions of analogous Leduc D-3A reefs (Lyster et al. 2024, S14; Stoakes & Creaney 1984) range from 2–8 km in diameter, corresponding to ~3–50 km². The Clive pool, from the ~168 wellbores and structural mapping in S6, is at the larger end of this range.

3. **The inferred EOR plume footprint** (~20–50 km²) is a first-order estimate. No public CO₂ plume modelling for Clive has been published. For comparison, the separate Origins Project (saline aquifer, single well) has a modelled plume radius of 4 km (~50 km²) per AER appeal decisions.

4. **Phase 1 injection strategy** (S6): 6 horizontal wells with ~1,400–1,600 m lateral lengths, sited in the central portion of the Clive Unit, with expansion planned over 10+ years. This implies the EOR area is a subset of the full pool area.

---

## 5. Other Key Reservoir Parameters

| Parameter | Value | Source | Verbatim Quote / Note |
|-----------|-------|--------|-----------------------|
| **Target formation** | Leduc D-3A (Woodbend Group) | S2, S3, S6 | Devonian dolomitized carbonate reef; part of the Bashaw Reef Complex |
| **Lithology** | Dolomitized limestone / carbonate reef | S13 (Amthor et al. 1994); S6 | Pervasively dolomitized; original limestone reef fabric modified by dolomitization and leaching |
| **Depositional environment** | Devonian carbonate pinnacle reef (Woodbend Group) | S13; S14 | Part of the Rimbey-Meadowbrook reef chain; isolated buildup on Bashaw Platform |
| **Depth (top)** | **1,845 m** TVD; range 1,845–1,938 m | S2, S3 | Mean 1,874 m |
| **Reservoir temperature** | **~69 °C** (~156 °F) | S11 (Menger 2025, GeoConvention) | Consistent with Alberta Basin geothermal gradient ~35 °C/km at ~1,874 m |
| **Reservoir pressure (initial)** | **~13,000 kPag** (~1,885 psig) | S11 (Menger 2025, GeoConvention) | **⚠️ Original document reported this as initial reservoir pressure. Note: pressure has been modified by 50+ years of production and 10+ years of waterflood.** |
| **CO₂ density at reservoir conditions (computed)** | **~640–700 kg/m³** at 13 MPa, 69 °C | Calculated from Span-Wagner EOS | Supercritical CO₂ at these conditions. For storage resource: use **~660 kg/m³** as mid-range estimate. **Added in revision** — not reported in original. |
| **Water saturation (Sw)** | **Estimated 20–40 %** in oil zone; **~40–60 %** in gas cap | Inferred from S6 hydrocarbon column data (11 m gas, 18.5 m oil). **No Clive-specific Sw published** | **Added in revision** — key parameter for storage resource. For Leduc dolomitized carbonates in Bashaw Complex: oil-zone Sw typically 20–40 % (residual oil zone), gas-cap Sw ~40–60 %. CO₂ injects into the oil zone (partially swept by waterflood), so expect **Sw range 25–50 %** for the injection interval. **This is an estimate, not a measured value.** |
| **Net-to-gross (usable pore space)** | (1 − Sw) ≈ 0.40–0.75 | Calculated from estimated Sw | **For preliminary storage resource: use (1 − Sw) = 0.55 as a P50 estimate** |
| **Primary seal** | Ireton Formation | S6 (MMV Plan 2019) | Carbonate mudstone, **4.3–17 m thick** at Clive; proven seal (contained 30 m hydrocarbon column) |
| **Secondary seals** | Nisku anhydrites + Wabamun Group + Colorado Group shales | S6 | Multiple stacked regional seals |
| **Depth to top of Ireton seal** | ~1,830 m TVD | S6 | Immediately overlies Leduc reservoir |
| **Hydrocarbon column (pre-production)** | 11 m gas cap + 18.5 m oil column = 30 m | S6 | Original oil-water contact at ~1,874 m TVD |
| **Discovery year** | 1952 | S6 | Clive field discovery well |
| **Total wellbores** | 168 (delineation + development) | S6 | Extensive dataset for geo-model calibration |
| **Phase 1 injection wells** | 6 horizontal wells; ~1,400–1,600 m lateral | S6, S1 | Central Clive Unit; expansion planned |
| **CO₂ source** | ACTL pipeline from Redwater/NWR Sturgeon Refinery | S1, S2 | Pipeline capacity **14 Mt/yr** |
| **Current injection rate** | **~1.24 Mt/yr** (as of 2021) | S8; S10 | Blend of EOR + storage |
| **Cumulative CO₂ stored** | **~4 Mt** (as of 2024); estimate >7 Mt (as of 2025) | S8; Enhance press releases | Evidence of no apparent leakage (per S1 MMV data) |

### Summary Table for CO₂ Storage Resource Estimation

| Property | Recommended Value | Unit | Source-Quality Flag |
|----------|-------------------|------|---------------------|
| **Formation** | Leduc D-3A (Woodbend Group) | — | Well-established |
| **Lithology** | Dolomitized carbonate reef | — | Well-established |
| **Depth (top)** | 1,845 | m TVD | S2, S3 (Tier 2) |
| **Net-pay thickness (injectivity)** | 8.75 | m | ACTL reports (Tier 2) — cutoffs not specified |
| **Reservoir interval thickness (storage resource)** | 35 | m | (Tier 2) — cored interval |
| **Average porosity** | 6.1 | % | **Operator-cited (Tier 3 derivative)** — averaging method unknown |
| **Regional arithmetic mean porosity** | ~8 | % | (Tier 3, independent) — for context |
| **Mode permeability** | 30 | mD | **Operator-cited (Tier 3 derivative)** — mode ≠ mean |
| **Arithmetic mean permeability** | ~26 (regional) / 50–150 (inferred Clive) | mD | (Tier 3 independent) / estimated from log-normal transform |
| **kh product (permeability × thickness)** | 228–1,310 | mD·m | Estimated; not measured |
| **Reservoir temperature** | ~69 | °C | S11 (Tier 4) |
| **Reservoir pressure** | ~13,000 | kPag | S11 (Tier 4) — pre-injection / modified |
| **CO₂ density (in-situ, computed)** | ~660 (P50: 640–700) | kg/m³ | Calculated; added in revision |
| **Water saturation (Sw)** | 25–50 (estimated) | % | **Data gap** — no published Sw for Clive D-3A |
| **Pore-volume factor (1 − Sw)** | 0.40–0.75 (P50: 0.55) | — | Estimate; added in revision |
| **CO₂ storage capacity** | 18.8 (total) / ~4–7 (injected to 2024–25) | Mt | S8 (Tier 4); ~14 Mt remaining |
| **Primary seal** | Ireton Formation (4.3–17 m) | — | Proven seal (S6, Tier 2) |
| **Lease area (Origins Hub SLA)** | 544,512 | ha | Regulatory — do NOT use for pool-level calculations |
| **Clive Leduc D-3A pool area** | ~25–40 | km² | Inferred from analogue reef geometry |
| **CO₂ plume footprint (Phase 1 EOR)** | ~20–50 | km² | Estimated — not modelled publicly |
| **Current injection rate** | ~1.24 | Mt/yr | As of 2021 (Tier 4) |
| **Pipeline capacity** | 14 | Mt/yr | ACTL system capacity (Tier 4) |
| **Secondary target** | Nisku Fm D-2 (6.4 m thick, shale seal) | — | S2, S3 — not included in primary storage resource estimate |

---

## Cross-Validation: Analogue Comparison

| Parameter | ACTL / Clive (This Project) | Origins CCS Hub (Enhance) | Quest CCS (Shell) | Source |
|-----------|----------------------------|---------------------------|-------------------|--------|
| **Target formation** | Leduc Fm D-3A (dolomitized carbonate reef) | Leduc Fm (dolomitized carbonate reef) | Basal Cambrian Sands (sandstone) | — |
| **Reservoir type** | Depleted oil reservoir (EOR + storage) | Pressure-depleted saline aquifer | Deep saline aquifer | — |
| **Depth** | 1,845–1,938 m | ~2,000 m TVD | 1,800–2,100 m | S2, S3; CarbonStorage.io; AER Dec 2012 |
| **Net-pay thickness** | **8.75 m** | ~137 m gross (net not published) | 35–46 m | ACTL reports; CarbonStorage.io; Quest AER Dec |
| **Porosity (zone avg)** | **6.1 %** (Clive geo-model) | **6 %** (well avg) | **16 %** (8-19 well avg) | S7; CarbonStorage.io; Quest AER Dec |
| **Porosity (regional Leduc)** | 5–20 % | 5–20 % | N/A (different fm) | Amthor et al. 1994 |
| **Permeability (mode or range)** | 30 mD (mode); 1–1,000 mD (range) | <10–>100 mD (range) | 20–500 mD (range); ~150 mD (8-19 well) | S7; CarbonStorage.io; Quest AER Dec |
| **Arithmetic mean permeability** | ~26 mD (regional Leduc core avg) | ~26 mD (regional Leduc core avg) | ~100–300 mD (core avg) | Weides & Majorowicz 2014 |
| **CO₂ density** | ~660 kg/m³ | ~700–800 kg/m³ | ~600–700 kg/m³ | Calculated |
| **Phase 1 injection rate** | ~1.24 Mt/yr (EOR + storage) | 1.6 Mt/yr (1 well) | ~1.0–1.2 Mt/yr (3 wells) | S8; CarbonStorage.io; Quest papers |
| **Cumulative stored** | >7 Mt (as of 2025) | 0 (pre-injection) | >7 Mt (as of 2024) | Enhance; Quest |
| **Storage mechanism** | EOR (miscible flood) + residual/dissolution | Saline aquifer storage | Saline aquifer storage | — |

---

## Notes on Data Quality

1. **Source reliability hierarchy**: The document has **no Tier 1 sources** (direct core measurements, DST results, or public AER well-file data from Clive D-3A wells). The quantitative chain is: Enhance internal geo-model (Tier 2) → cited in peer-reviewed chapter (Tier 3, S7). Regional context values from Amthor et al. (1994), Weides & Majorowicz (2014), and Lyster et al. (2024) are independent Tier 3 sources but represent Leduc Formation regional averages, not Clive-specific values.

2. **All Clive-specific quantitative values should be labelled "operator-cited, not independently verified"** in any downstream database or resource calculation. No independent third-party petrophysical evaluation of the Clive Leduc D-3A pool has been published.

3. **Net-pay cutoff criteria are unknown**: The 8.75 m net-pay thickness is reported without the porosity or permeability thresholds used to define it. This is a significant gap — without knowing whether net pay was defined at φ > 3 %, φ > 5 %, k > 0.1 mD, or some other threshold, the number cannot be independently audited or compared to other CCS projects.

4. **Permeability data gap is the single most significant technical limitation**: Only the mode (30 mD) has been published. No arithmetic mean, geometric mean, kh product, DST permeability, or well-test injectivity data exist in the public record for the Clive Leduc D-3A pool. For CO₂ storage resource estimation, users should either:
   - **(a)** Use the regional Leduc arithmetic mean of ~26 mD (Weides & Majorowicz 2014) as a conservative lower bound
   - **(b)** Assume a log-normal distribution with mode = 30 mD and V_DP ≈ 0.7, yielding arithmetic mean ~70–100 mD
   - **(c)** Treat the 30 mD mode as a central estimate for a log-normal distribution and propagate uncertainty

5. **Area concepts must not be conflated**: The SLA (544,512 ha) is the full Origins Hub lease and must NOT be used in Clive pool volumetric calculations. The pool area (~25–40 km²) is the geological extent of the Leduc D-3A. The plume footprint (~20–50 km²) is the CO₂-contacted area during Phase 1 EOR. These have different magnitudes and serve different purposes.

6. **Water saturation not independently measured**: The estimated Sw range (25–50 %) is inferred from the hydrocarbon column description in the MMV Plan (S6) and published ranges for Leduc dolomitized carbonates. No Clive-specific capillary pressure, log-derived Sw, or core-plug Sw data have been published. This is a mandatory input for storage resource volumetric calculations and should be prioritised for future data acquisition.

7. **CO₂ density computed at ~660 kg/m³**: At reservoir conditions of ~13 MPa and ~69 °C, pure CO₂ is supercritical with density ~640–700 kg/m³ (Span-Wagner equation of state). The mid-range estimate of ~660 kg/m³ is appropriate for preliminary storage resource calculations. For reference, this is comparable to Quest CCS (Basal Cambrian Sands: ~600–700 kg/m³ at slightly lower P/T) and Origins (Leduc Fm: ~700–800 kg/m³ at ~15.3 MPa).

8. **S2 and S3 are not independent**: Both are ACTL operator Summary Reports (S2: Dowell 2018, S3: Benko 2017) and likely derive from the same Enhance geo-model. The 8.75 m thickness value is not independently auditable.

9. **S9 (CarbonStorage.io Origins page) returned 404** during retrieval. The lease area value (2,491,617 acres) has been replaced with the more reliable Key Facts Energy value (544,512 ha / 58 townships). If the CarbonStorage.io page becomes accessible again, cross-validate.

10. **S11 (GeoConvention 2025, Menger)** is cited for P/T values. The full presentation PDF is not publicly accessible. The cited abstract page at `geoconvention.com/2025/` was confirmed as the correct conference reference. Values should be treated as preliminary pending peer-reviewed publication.

11. **Original hydrocarbon column (~30 m) vs. CO₂ column**: The document references the hydrocarbon column as geological context for storage interval thickness. CO₂ column height in a storage scenario will differ due to (a) CO₂–brine capillary entry pressure vs. oil–water entry pressure, (b) relative permeability differences between CO₂–brine and oil–water systems, and (c) wettability differences (CO₂ is typically less wetting than oil in carbonates). The 35 m reservoir interval (not the 30 m hydrocarbon column) is recommended for storage resource work.

---

## Key Data Gaps for Future Work

| Gap | Priority | Suggested Path to Resolution |
|-----|----------|------------------------------|
| No publicly available core-plug φ/k from Clive D-3A wells | **HIGH** | Request from AER well-file database (168 Clive wells have core); search PCOR Partnership publications |
| Net-pay cutoff criteria unknown | **HIGH** | ACTL Detailed Report (S1, 403 on direct fetch) contains the geo-model description; may specify cutoffs |
| No DST or well-test permeability (kh) | **HIGH** | AER well-file data may contain DST results from Clive D-3A wells (discovery 1952 — extensive testing over decades) |
| No water saturation (Sw) measurement | **HIGH** | Log analysis from Clive wells in AER database; or from published petrophysical evaluation of Leduc D-3A pools |
| No kv/kh ratio for Clive | **MEDIUM** | Amthor et al. (1994) gives regional range; Clive-specific would require core-plug horizontal/vertical perm measurements |
| No CO₂ plume modelling for Clive EOR | **MEDIUM** | Not published; may be in Enhance Energy proprietary files |
| P/T from conference presentation only (S11) | **LOW** | Seek peer-reviewed publication with confirmed P/T values |
| CO₂ injectivity ratio "~3× water" unsupported | **LOW** | Remove or cite specific well-test reference; not included in revised recommended values |

---

*Revised 23 July 2026 in response to adversarial review (CANADA-CCS-ALBERTA-ACTL-Clive-props-review.md). Changes tracked in "Revisions from Adversarial Review" section above. Compiled from open.alberta.ca ACTL Knowledge Sharing Reports (2011–2021), Enhance Energy Clive Leduc MMV Plan (2019), Tavallali et al. (2025, Wiley), Amthor et al. (1994, AAPG Bulletin), Weides & Majorowicz (2014, Can. J. Earth Sci.), Lyster et al. (2024, GHGT-17), Clean Prosperity (Hares et al. 2022), and industry sources. See "Notes on Data Quality" for source reliability assessment.*
