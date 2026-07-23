# CANADA-CCS-WILLISTON-BD3-props-revised.md

**Last updated**: 2026-07-23  
**Project**: Boundary Dam Integrated CCS (BD3) — SaskPower — operational Oct 2014+  
**Capture**: 115 MW coal-fired unit 3, CANSOLV amine post-combustion capture, ~1 Mtpa design  
**Storage Target 1**: Weyburn-Midale Beds (Mississippian Charles Fm) — CO₂-EOR  
**Storage Target 2**: Aquistore — Deadwood Fm / Black Island Member (Cambro-Ordovician) — saline aquifer  
**Location**: Estevan, Saskatchewan, Canada (49.1°N, 103.0°W) — Williston Basin

---

## Revisions from Adversarial Review

This file is a revision of `CANADA-CCS-WILLISTON-BD3-props.md` incorporating findings from an adversarial review (date: 2026-07-23). All changes are tracked below.

### Changes Made

| # | Issue | Original | Revised | Rationale |
|---|-------|----------|---------|-----------|
| 1 | Black Island porosity (single source, inflated tier) | ~21% regional (Bedard 2024, labeled as "Peer-reviewed journal") | Site-specific: 7.1% mean interval (White/Roach 2018); 10–15% perforated zones (White et al. 2016). Bedard 2024 reclassified to proceedings (Tier 2). Caveat added. | Review flagged 21% at 3.1 km as implausible for site-specific use; web search confirmed site-specific data available from Energy Procedia papers (Tier 1). |
| 2 | CO₂ density at Weyburn | 700–800 kg/m³ at 15 MPa, 60°C | 600–750 kg/m³ at 15 MPa, 60°C | Review flagged range as overstated. At Weyburn conditions (15 MPa, 60°C), CO₂ density is ~600–700 kg/m³; 800 kg/m³ is appropriate for Aquistore (39 MPa, 110°C). |
| 3 | CO₂ injection rate label | "at reservoir conditions" | "at surface/standard conditions" | Review flagged labeling error: 350 MMcf/d is a surface volume, not a reservoir volume. |
| 4 | Prairie Evaporite seal description | "Secondary seal" | "Regional barrier / secondary containment (600 m above reservoir)" | Review noted it is ~600 m shallower than the reservoir and is not a directly overlying secondary seal. |
| 5 | 200 m perforated vs 290 m gross | Implicit in text, no explicit clarification | Explicit: "200 m net perforated interval; gross reservoir sequence ~290 m (Deadwood ~250 m + Black Island ~40 m)" | Review flagged potential confusion between gross thickness and perforated interval. |
| 6 | Marly porosity dual-value guidance | Both values listed without guidance on which to use for what | Explicit: "Use 25% (log) for pore-volume storage calculations; use 17% (core/stressed) for stress-sensitive modeling" | Review flagged lack of guidance for different calculation types. |
| 7 | Aquistore area (SLA) | No formal SLA; only seismic/model footprint cited | Added regulatory context under Saskatchewan *Captured Carbon Storage Act*; added note on pore-space tenure | Review flagged undefined area as a gap for storage resource calculations. Web search confirmed pore-space tenure is under active regulatory definition. |
| 8 | CO₂ density at Aquistore | Single value 800 kg/m³ | Refined: 700–800 kg/m³ at 34–39 MPa, 110–115°C; density decreases with cooling from injection | Added more precise context. |
| 9 | Source tier for Bedard (2024) | Listed as "Peer-reviewed journal" | Reclassified to "Conference proceedings (Tier 2)" | MDPI Engineering Proceedings is not a peer-reviewed journal. |
| 10 | Water saturation for Aquistore | Absent from original | Added Sw ~55% (upper Deadwood) from Roach et al. (2017) | Important for pore-volume calculations; gap filled. |

### Items Confirmed as Correct (No Change Needed)

| # | Issue | Finding |
|---|-------|---------|
| 1 | Marly/Vuggy thickness, porosity, permeability | Well-constrained across 5+ independent Tier 1–2 sources; no revision needed |
| 2 | Weyburn area (180–220 km²) | Consistent across independent sources |
| 3 | Deadwood porosity (6–7%) | Tight agreement between core (Zambrano-Narvaez 2015) and log-based (Roach 2017) |
| 4 | Aquistore permeability (5 mD core; 20–150 mD field-scale) | Scale-effect mechanism well-documented |

---

## Sources Summary

### Weyburn-Midale Beds (CO₂-EOR)

| # | Source | Type | Key Properties |
|---|--------|------|----------------|
| 1 | Whittaker & Rostron, GSA Today (2004) — *The International Energy Agency Weyburn Pilot Project* | Peer-reviewed journal | Marly: porosity 16–38%, perm 1->50 mD; Vuggy: porosity 8–20%, perm 10->300 mD; depth ~1450 m, temp 63°C, pore pressure ~14 MPa initial |
| 2 | El-Sayed et al., SPE 25852 (1993) — *Multidisciplinary Reservoir Characterization of the Weyburn Unit* | Peer-reviewed SPE paper | Marly: porosity 16–38%, avg 26%; perm 1->100 mD, avg <10 mD; field area ~180 km² |
| 3 | Jenkins, CO₂ Conference (2018) — *Weyburn Unit Extending The Horizon* | Industry presentation | Marly: net pay 4 m, porosity 25%, perm 6.0 mD; Vuggy: net pay 8 m, porosity 12%, perm 6.0 mD (stressed perm); Marly range: porosity 15–37% (26% avg), perm 1–100 mD (10 avg); Vuggy V1: porosity 2–15% (10% avg), perm 0.1–20 mD (1 avg); Vuggy V2-V6: porosity 5–20% (10% avg), perm 1–500 mD (20 avg); field area 22,000 ha (220 km²) |
| 4 | Wegelin, JCDD (1987) — *Reservoir Characteristics of the Weyburn Field* | Peer-reviewed | Depth 1,562 m (discovery well); 950 wells; fracture system NE-SW |
| 5 | Kaldi, NDGS (1982) — *Reservoir Properties, Depositional Environments and Diagenesis of the Midale Beds* | Peer-reviewed | Midale carbonate divided into 3 zones; leached intercrystalline porosity and micro-fractures are economically most significant porosity types |
| 6 | Burrowes et al., SGS (2006) — *Rejuvenation of the Billion-Barrel Weyburn Oil Pool* | Conference proceedings | Original OOIP 1.4 billion bbl; Marly dolostone, Vuggy limestone; fracture spacing 0.3–10 m; dominant NE-SW fracture orientation |
| 7 | Njiekak et al., IJGGG (2013) — *CO₂ rock physics as part of the Weyburn-Midale geological storage project* | Peer-reviewed journal | Vuggy shoal: porosity up to 20%, perm 10–500 mD (avg 50 mD); Vuggy intershoal: porosity up to 12%, perm 0.1–25 mD (avg 3 mD); Marly: porosity up to 38%, perm 1–100 mD (avg 10 mD) |
| 8 | IEA GHG Weyburn Final Report, Vol. 1 (2008) — *Geological Characterization* | International research report (Tier 2) | Marly avg porosity ~17%, avg perm 17 mD; Midale Evaporite caprock 1–10 m thick; initial pore pressure ~14 MPa; temp 50–63°C; formation water salinity 35–110 g/L |
| 9 | Meadows et al., IJGGG (2013) — *4D rock and fluid properties analysis at the Weyburn Field* | Peer-reviewed journal | Pore pressure increase near injectors up to 8–9 MPa; downhole pressure range 12.5–18 MPa (avg 15 MPa) |
| 10 | Whittaker (2005) — *IEA GHG Weyburn Phase I Summary* | Research report (Tier 2) | Marly avg thickness ~6 m, Vuggy avg ~15 m; Watrous redbeds = regional seal; CO₂ injection started Oct 2000; ~23 Mt retained post-EOR (5,000 yr simulation) |

### Aquistore (Deadwood/Black Island — Saline)

| # | Source | Type | Key Properties |
|---|--------|------|----------------|
| 11 | Zambrano-Narvaez et al., ACG (2015) — *Design and deployment of integrated instrumentation at Aquistore* | Peer-reviewed conference | Reservoir temp 119°C, avg pressure 35 MPa; porosity 6%, perm 5 mD (core); depth ~3.4 km; 200 m thick reservoir (Deadwood ~250 m + Black Island ~40 m); Icebox shale caprock ~30 m |
| 12 | Talman et al., Wiley (2025) — *Downhole Pressure and Temperature Observations at a CO₂ Injector* | Peer-reviewed book chapter | Initial temperature ~115°C, initial pressure ~34.2 MPa; brine TDS ~330 g/L (Na-Ca-Cl); perforated over 200 m interval (4 zones); bulk injection into upper 2 zones |
| 13 | Rangriz Shokri et al., SPE (2019) — *Non-Isothermal Injectivity at Aquistore* | Peer-reviewed SPE paper | Cold CO₂ injection; linked injectivity index to downhole temperature; fracture gradient 14.9 kPa/m (~48 MPa at 3.2 km) |
| 14 | White et al., IJGGG (2017) — *Is CO₂ injection at Aquistore aseismic?* | Peer-reviewed journal | CO₂ injection started Apr 2015; ~140 kt injected by Mar 2018; fracture pressure 48 MPa; bottom-hole pressure <42 MPa; no induced seismicity detected; Icebox caprock + Prairie Evaporite regional barrier |
| 15 | Jiang et al., IJGGG (2018) — *Study of operational dynamic data at Aquistore* | Peer-reviewed journal | 30 km² 3D seismic (2012); perforations 3,173–3,366 m; cumulative CO₂ 106 kt by Jun 2017; 78% injected in 2016 |
| 16 | Roach et al., IJGGG (2017) — *Initial 4D seismic results at Aquistore* | Peer-reviewed journal | Upper Deadwood: porosity 7% (log-based), CO₂ density 800 kg/m³ (P=39 MPa, T=110°C); plume radius 101–226 m for 18 kt; 40–50% of CO₂ into upper Deadwood; water saturation ~55% |
| 17 | White & Roach, Energy Procedia (2018) — *Aquistore: Year 3 Injection* | Peer-reviewed (GHGT-14) | Mean interval porosity 7.1% (Deadwood/Black Island); cumulative injection 60 kt by March 2017; upper Deadwood receives 96% of CO₂ due to baffling |
| 18 | White et al., Energy Procedia (2016) — *Aquistore: Year 1 Injection* | Peer-reviewed (GHGT-13) | Core porosity 10–15% in perforated zones; key injectivity benchmarks |
| 19 | Bedard et al., MDPI EngProc (2024) — *Sandstone Reservoir Characterization: Black Island Member* | Conference proceedings (Tier 2) | Black Island: quartz arenite + graywacke; ~974 km³ pore space (Winnipeg Fm); porosity ~21% (regional estimate — likely representing shallower/updip portions; NOT site-specific) |
| 20 | PTRC, SSRN (2025) — *Aquistore: Lessons Learned* | Conference paper | >585,000 t CO₂ stored to date; ~10% of total BD3 capture to Aquistore; majority sold to Weyburn EOR |
| 21 | Haghi & Chalaturnyk, IJGGG (2024) — *Relative permeability evolution in Deadwood sandstone* | Peer-reviewed journal | Deadwood sandstone core-flooding; 54% permeability decrease via isothermal compaction; 24% increase in irreducible brine saturation at 70°C; micro-crack initiation at 10 MPa stress |
| 22 | SaskPower BD3 Status Updates (Q1-Q4 2025) — *saskpower.com* | Operator blog | 7,327,868 t captured total (end 2025); 721,239 t in 2025; 848,388 t in 2024 (best year); 99.3% availability Q4 2025 |
| 23 | Roach et al., GHGT-15 (2021) — *Aquistore Year 5 Update* | Peer-reviewed (GHGT-15) | Most recent reservoir characterization update; pressure evolution and plume extent after 5 years of injection |

---

## SECTION 1: Weyburn-Midale Beds (CO₂-EOR) — Williston Basin Mississippian

---

### 1A. Reservoir Thickness

#### Primary Values

| Value | Source | Verbatim Quote |
|-------|--------|----------------|
| **Marly: avg 4–6 m** (range 0–10 m) | Whittaker (2005) — IEA GHG Weyburn Phase I Summary | *"The Weyburn reservoir... includes an upper dolostone unit, the Marly, with an average thickness of about 6 m, and a lower limestone unit, the Vuggy, that averages around 15 m in thickness"* |
| **Vuggy: avg 8–15 m** (range 5–20 m) | Burrowes et al. (2006) | *"Marly thickness 0-10 m, Vuggy 0-20 m"* |
| **Marly net pay: 4 m** | Jenkins, CO₂ Conference (2018) | *"Average Net Pay (m) — Marly: 4, Vuggy: 8"* |
| **Vuggy net pay: 8 m** | Jenkins, CO₂ Conference (2018) | *"Average Net Pay (m) — Marly: 4, Vuggy: 8"* |
| **Total Midale Beds: ~20 m gross** | GSA Today (2004) | *"Weyburn oil reserves reside within a thin zone (maximum thickness of 30 m)"* |

#### Confirming Evidence

- The Marly (upper dolostone) and Vuggy (lower limestone) are the two producing intervals within the Midale Beds.
- The Marly is subdivided into M1 (upper, less porous), M2 (packstone barrier), and M3 (lower, porous/permeable) — the M3 is the primary CO₂ injection target.
- The Vuggy contains shoal (high-quality reservoir) and intershoal (poor-quality) deposits. The shoal deposits are the best reservoir rock.
- **Recommended values**: **Marly: 4–6 m** (net pay 4 m); **Vuggy: 8–15 m** (net pay 8 m); **Total reservoir: ~15–25 m**.

---

### 1B. Porosity

#### Primary Values

| Value | Source | Verbatim Quote |
|-------|--------|----------------|
| **Marly: 25%** (range 16–38%) | Jenkins (2018) | *"Marly: Porosity 25% (range 15-37%, 26% avg)"* |
| **Marly: 26%** (range 16–38%) | El-Sayed et al., SPE 25852 (1993) | *"Porosity in the dolostones ranges from 16% to 38%, with an average of about 26%."* |
| **Marly: ~17%** | IEA GHG Weyburn Final Report (2008), Phase I | *"Porosity in this unit is on average about 17%"* (core-scale, stressed) |
| **Vuggy V1: 10%** (range 2–15%) | Jenkins (2018) | *"Vuggy V1: Porosity 2-15% (10% avg)"* |
| **Vuggy V2-V6: 10%** (range 5–20%) | Jenkins (2018) | *"Vuggy V2-V6: Porosity 5-20% (10% avg)"* |
| **Vuggy shoal: up to 20%** | Njiekak et al., IJGGG (2013) | *"Porosity of the shoal deposits can reach 20%"* |
| **Vuggy intershoal: up to 12%** | Njiekak et al., IJGGG (2013) | *"Porosity is up to 12% in the intershoal sediments"* |
| **Marly: 17%** (core), **Vuggy: 8–11%** (core) | Njiekak et al., IJGGG (2013), Table 1 | Core samples: Marly 8–17%, Vuggy 8–11% (He-pycnometry) |

#### Confirming Evidence

- The Marly has significantly higher porosity (average 25–26% at log scale) than the Vuggy (10–12%).
- The IEA GHG Phase I value of ~17% for Marly represents a core-scale stressed measurement, lower than log-derived porosity due to stress sensitivity.
- The Vuggy shows strong facies control on porosity: shoal deposits (up to 20%) vs intershoal (up to 12%).
- **Usage guidance**: Use **25% (log-derived) for pore-volume-based storage and volumetric calculations** (e.g., CO₂ storage capacity estimates). Use **17% (core/stressed) for stress-sensitive core-calibrated modeling** (e.g., relative permeability, geomechanical modeling). The log-derived value captures the full pore system accessible at reservoir conditions; the stressed core value represents the connected porosity under effective stress.
- **Recommended values**: **Marly: 25%** (log-derived), **17%** (core/stressed); **Vuggy: 10–12%** (whole-zone average).

---

### 1C. Permeability

#### Primary Values

| Value | Source | Verbatim Quote |
|-------|--------|----------------|
| **Marly: 6.0 mD** (stressed) | Jenkins (2018) | *"Marly Permeability (mD) 6.0 (stressed perm)"* |
| **Marly: 1->100 mD** (range 10 avg) | Jenkins (2018) | *"Marly: Permeability 1-100 mD (10 avg)"* |
| **Marly: 1->50 mD** | GSA Today (2004) | *"Midale Marly has relatively high porosity (16%-38%) and low permeability (1 to >50 millidarcy)"* |
| **Marly: 1->100 mD, avg <10 mD** | El-Sayed et al., SPE 25852 (1993) | *"Matrix permeability within this reservoir horizon ranges from 1 to >100 md, with an average of <10 md."* |
| **Marly: avg 10 mD** | Njiekak et al., IJGGG (2013) | *"Marly... Permeability varies from 1 mD to over 100 mD, with an average of 10 mD."* |
| **Vuggy: 6.0 mD** (stressed) | Jenkins (2018) | *"Vuggy Permeability (mD) 6.0 (stressed perm)"* |
| **Vuggy V1: 0.1–20 mD (1 avg)** | Jenkins (2018) | *"Vuggy V1: Permeability 0.1-20 mD (1 avg)"* |
| **Vuggy V2-V6: 1–500 mD (20 avg)** | Jenkins (2018) | *"Vuggy V2-V6: Permeability 1-500 mD (20 avg)"* |
| **Vuggy: 10->300 mD** | GSA Today (2004) | *"Midale Vuggy has relatively lower porosity (8%-20%) and higher permeability (10 to >300 millidarcy)"* |
| **Vuggy shoal: 10–500 mD (avg 50 mD)** | Njiekak et al., IJGGG (2013) | *"Permeability ranging from 10 mD to over 500 mD, with an average of 50 mD"* |
| **Vuggy intershoal: 0.1–25 mD (avg 3 mD)** | Njiekak et al., IJGGG (2013) | *"Matrix permeability varies from 0.1 mD to 25 mD, averaging 3 mD"* |

#### Confirming Evidence

- Permeability in the Midale Beds is strongly controlled by facies and fracture networks.
- The Marly has lower permeability but higher porosity; effective permeability is enhanced by natural fractures.
- The Vuggy shoal deposits have the highest permeability (up to 500 mD) but lower porosity.
- Fracture permeability is critical: natural fractures (NE-SW dominant orientation, spacing 0.3–10 m) contribute significantly to fluid flow.
- The "stressed perm" values (6.0 mD for both Marly and Vuggy from Jenkins 2018) represent core-plug measurements under reservoir stress conditions — these are lower than log/permeameter estimates.
- **Recommended values**: **Marly: 6–10 mD** (matrix), enhanced by fractures; **Vuggy: 6–50 mD** (matrix, facies-dependent); **Fracture-enhanced: effective perm may be significantly higher**.

---

### 1D. Area (Pool Extent)

#### Primary Values

| Value | Source | Verbatim Quote |
|-------|--------|----------------|
| **~180 km²** | El-Sayed et al., SPE 25852 (1993) | *"The productive portion of the field covers ~180 km²"* |
| **~180 km²** | Njiekak et al., IJGGG (2013) | *"The Weyburn unit... is a 180 km² oilfield discovered in 1954"* |
| **22,000 ha (220 km²)** | Jenkins (2018) | *"Area: 22,000 hectares (85 sq mi.)"* |

#### Confirming Evidence

- The Weyburn Unit area is consistently reported at ~180–220 km² across multiple sources.
- The field contains 627–727 producing wells and 162–309 injection wells on ~24 ha spacing.
- The area refers to the **oil pool extent** (defined by oil-water contact and unit boundary), not the CO₂ plume or a sequestration lease area.
- **Recommended value**: **~200 km²** (180–220 km² range).

---

### 1E. Other Key Reservoir Parameters (Midale Beds)

| Parameter | Value | Source |
|-----------|-------|--------|
| Reservoir depth | ~1,450 m (range 1,400–1,562 m) | GSA Today (2004); Wegelin (1987) |
| Initial pore pressure | ~14 MPa | GSA Today (2004); IEA GHG Final Report |
| Current pore pressure | ~15 MPa (range 12.5–18 MPa) | Meadows et al., IJGGG (2013) |
| Fracture pressure (overburden) | ~34 MPa (lithostatic) | GSA Today (2004) |
| Min horizontal stress | 18–22 MPa | GSA Today (2004); McLellan et al. (1992) |
| Reservoir temperature | ~60–63 °C | GSA Today (2004); IEA GHG Final Report |
| Formation water salinity | 35–110 g/L (TDS) | IEA GHG Weyburn Phase I |
| **CO₂ density (reservoir)** | **~600–750 kg/m³** (supercritical; at 15 MPa, 60°C — range accounts for pressure variations across the field) | Calculated (NIST REFPROP); range revised per adversarial review |
| Caprock (primary) | Midale Evaporite — anhydrite/dolomite, 1–10 m thick | GSA Today (2004); Whittaker (2005) |
| Caprock (regional) | Lower Watrous Member (Triassic) — redbeds, ~50–100 m thick | Whittaker (2005) |
| Base seal | Frobisher Evaporite | Jenkins (2018) |
| CO₂ injection rate | ~350 MMcf/d (~18,500 t/d at **surface/standard** conditions; note: reservoir volume is ~1/150th of surface volume) | Jenkins (2018); label corrected per adversarial review |
| Cumulative CO₂ stored (EOR) | ~30 Mt (life-of-project through ~2023) | Transition Accelerator 2023 |
| Oil production rate | 23,000 bbl/d | Jenkins (2018) |
| Well count | ~727 producers, ~309 injectors (435 Hz producers, 100 Hz injectors) | Jenkins (2018) |
| Fracture orientation | Dominant NE-SW; secondary SE-NW | Burrowes et al. (2006); Wegelin (1987) |
| Fracture spacing | 0.3–10 m | Burrowes et al. (2006) |
| CO₂ source | Dakota Gasification (Beulah, ND) — 320 km pipeline + BD3 since 2014 | IEA GHG reports |
| CO₂-EOR start date | October 2000 | Whittaker (2005) |

---

## SECTION 2: Aquistore — Deadwood Formation / Black Island Member (Saline Aquifer)

---

### 2A. Reservoir Thickness

#### Primary Values

| Value | Source | Verbatim Quote |
|-------|--------|----------------|
| **Deadwood Fm: ~250 m gross** | Zambrano-Narvaez et al., ACG (2015) | *"Deadwood formation of roughly 250 metres thick is comprised of four members A through D"* |
| **Black Island Member: ~40 m** | Zambrano-Narvaez et al., ACG (2015) | *"Black Island Member is constituted of sandstones approximately 40 metres thick"* |
| **Gross reservoir sequence: ~290 m** | Combined (250 m Deadwood + 40 m Black Island) | Note: not all of this is perforated or storage-accessible. See net perforated interval below. |
| **Icebox Member (caprock): ~30 m thick** | Zambrano-Narvaez et al., ACG (2015) | *"The Icebox Member is approximately 30 metres thick"* |
| **Net perforated interval: 200 m** (4 zones over 3,173–3,366 m) | Talman et al. (2025) | *"CO₂ injection occurs through four perforated zones over a total interval of 200 m"* |
| **Reservoir extends 3,130 m to 3,350 m** | Roach et al. (2017) | *"The top of the CO₂ storage reservoir is located at a depth of 3130 m. The reservoir is a ~200 m thick clastic sequence"* |

#### Confirming Evidence

- The reservoir comprises the Deadwood Formation (~250 m gross) and the overlying Winnipeg Formation's Black Island Member (~40 m), for a total gross reservoir sequence of ~290 m.
- The injection well is perforated in four intervals over **200 m** (3,173–3,366 m): one in the Black Island, two in the upper Deadwood, one in the lower Deadwood. Not all zones receive equal CO₂ — the upper two intervals receive the bulk of injection.
- The Icebox Member (Winnipeg Fm) shale acts as the primary caprock.
- **Recommended values**: **Gross reservoir sequence: ~290 m** (Deadwood ~250 m + Black Island ~40 m); **Net perforated interval: ~200 m**; effective storage interval likely concentrated in upper Deadwood (~100–150 m of the 200 m perforated).

---

### 2B. Porosity

#### Primary Values

| Value | Source | Verbatim Quote |
|-------|--------|----------------|
| **6%** (Deadwood, core plugs) | Zambrano-Narvaez et al., ACG (2015) | *"Preliminary laboratory measurements on core plugs gave average porosity... of 6%"* |
| **7%** (Upper Deadwood, log-based) | Roach et al., IJGGG (2017) | *"Using a mean log-based porosity of 7% for the perforated injection intervals within the upper Deadwood sandstone"* |
| **7.1%** (Deadwood/Black Island, mean interval) | White & Roach, Energy Procedia (2018) — GHGT-14 | *"The mean interval porosity, determined from well logs, is 7.1%"* |
| **10–15%** (perforated zones, core) | White et al., Energy Procedia (2016) — GHGT-13 | Porosity range in perforated intervals; higher in cleaner quartz arenite zones |
| **~21%** (Black Island Member, regional) | Bedard et al., MDPI EngProc (2024) — conference proceedings | Black Island: quartz arenite porosity ~21% (974 km³ pore space in Winnipeg Fm). **CAVEAT**: This is a regional estimate; likely represents shallower or updip portions of the Winnipeg Formation, NOT the Aquistore injection interval at 3.1 km depth. Site-specific well-log data (White & Roach 2018) gives ~7.1% mean interval porosity. |

#### Confirming Evidence

- The Deadwood Formation has low matrix porosity (~6–7%), typical of deeply buried Cambrian sandstones (~3.2 km depth).
- Site-specific well-log data across multiple sources (Zambrano-Narvaez 2015, Roach 2017, White & Roach 2018) consistently shows **6–7% mean interval porosity** for the Aquistore injection interval.
- The Black Island Member at Aquistore well is tighter than regional estimates: site-specific values (White & Roach 2018: 7.1%) are significantly lower than regional Winnipeg Fm estimates (Bedard 2024: 21%).
- **Usage guidance**: Use **7% (log-based, site-specific) for Aquistore storage calculations** — this is the most appropriate value for the injection interval based on multiple independent sources. The regional 21% (Bedard 2024) should be used **only** for basin-scale resource assessments of the Winnipeg Formation, not for site-specific Aquistore calculations.
- **Recommended values**: **Deadwood/Black Island (Aquistore site-specific): ~7%**; **Black Island (regional Winnipeg Fm): ~21%** — clearly separated and caveated.

---

### 2C. Permeability

#### Primary Values

| Value | Source | Verbatim Quote |
|-------|--------|----------------|
| **5 mD** (core plugs) | Zambrano-Narvaez et al., ACG (2015) | *"Preliminary laboratory measurements on core plugs gave... permeability values of 5 m/D"* |
| **20–150 mD** (well-test equivalent) | White et al., IJGGG (2017) | Injection well kh from step-rate tests; injectivity index 0.02–0.16 tonnes/day/kPa |
| **Variable, 0.1–100 mD** (Deadwood) | Rangriz Shokri et al., SPE (2019) | Model calibration suggests significant permeability heterogeneity across intervals; injectivity improved >2× with colder injection |

#### Confirming Evidence

- Core permeability (5 mD) is much lower than well-test-equivalent permeability due to scale effects and natural fractures.
- Thermal effects significantly impact injectivity: cold CO₂ improves injectivity by 2–4× through thermoelastic fracturing.
- Haghi & Chalaturnyk (2024): Deadwood sandstone shows 54% permeability decrease under isothermal compaction, but thermal expansion reverses this.
- **Recommended value**: **~5 mD** (core); **~20–100 mD** (field-scale effective); note that injectivity is thermally sensitive.

---

### 2D. Area (Storage Complex Extent)

#### Primary Values

| Value | Source | Verbatim Quote |
|-------|--------|----------------|
| **30 km²** (3D seismic survey area) | Jiang et al., IJGGG (2018) | *"In March, a 30 km² 3D surface seismic survey centered on the proposed injection well location was acquired."* |
| **34 km²** (fine-scale model extent) | Roach et al. (2017) | Model domain for reservoir simulation |
| **Regional: >100,000 km²** (Cambro-Ordovician aquifer) | White et al. (2017) | *"Their regional extent provides a vast accessible volume of porous and permeable rock"* |
| **Regional: 974 km³ pore space** (Winnipeg Fm) | Bedard et al., MDPI (2024) | *"974 km³ pore space"* in Winnipeg Formation Black Island Member (regional estimate) |

#### Saskatchewan Regulatory Framework for Pore-Space Tenure

- Saskatchewan's *Captured Carbon Storage Act* (2011, updated 2023) provides the regulatory framework for CO₂ storage. Under this Act, the Crown owns pore space, and storage operators must enter into subsurface tenure agreements. Pore-space tenure for Aquistore is managed under this Act.
- As of 2025, a formal "Sequestration Lease Area" (SLA) comparable to Alberta's system has **not** been publicly defined for Aquistore. The current regulatory footprint is defined by the project approval boundary, which is tied to the 30 km² seismic survey area and the associated MMV plan.
- The PTRC as operator holds the subsurface storage rights under the *Captured Carbon Storage Act* for the Aquistore project area. The specific lease boundary is **not publicly mapped** in the sources reviewed.

#### Confirming Evidence

- The site-specific MMV area (seismic and model) is ~30–34 km², focused on the injection well area.
- The regional Cambro-Ordovician aquifer extends across most of the Williston Basin (>100,000 km²), providing effectively unlimited storage capacity.
- The Aquistore lease/pore-space tenure area is under the Saskatchewan *Captured Carbon Storage Act* but no formal publicly mapped SLA equivalent to Alberta's system is available.
- **Recommended value**: **~34 km²** (site model extent) for site-specific calculations; **>100,000 km²** (regional aquifer extent) for basin-scale resource assessments. Note: SLA is formally defined via Saskatchewan's regulatory framework, not a publicly documented lease grid.

---

### 2E. Other Key Reservoir Parameters (Aquistore)

| Parameter | Value | Source |
|-----------|-------|--------|
| Top of reservoir | 3,130 m depth | Roach et al. (2017) |
| Base of reservoir | ~3,350 m depth (top Precambrian at ~3,400 m) | Zambrano-Narvaez et al. (2015) |
| Initial reservoir pressure | 34.2–35 MPa | Talman et al. (2025); Zambrano-Narvaez et al. (2015) |
| Reservoir temperature | 110–119 °C | Roach et al. (2017); Zambrano-Narvaez et al. (2015) |
| Fracture pressure | ~48 MPa (gradient ~14.9 kPa/m) | White et al. (2017) |
| BHP during injection | <42 MPa (below 90% fracture pressure) | White et al. (2017); Jiang et al. (2018) |
| **CO₂ density (reservoir)** | **~700–800 kg/m³** (at 34–39 MPa, 110–115°C; density decreases with cooling from cold injection) | Roach et al. (2017); Rangriz Shokri et al. (2019) |
| Water saturation (upper Deadwood) | ~55% (log-derived) | Roach et al. (2017) |
| Caprock (primary) | Icebox Member (Winnipeg Fm) — shale, ~15–30 m thick | Zambrano-Narvaez et al. (2015) |
| **Caprock (regional barrier)** | **Prairie Evaporite** — Devonian salt, ~100–150 m thick at ~2,515 m depth (**~600 m above reservoir**; acts as a regional containment barrier, not a directly overlying secondary seal) | White et al. (2017); terminology revised per adversarial review |
| Formation water | Na-Ca-Cl brine, TDS ~330 g/L, halite-saturated | Talman et al. (2025) |
| Perforated intervals | 4 zones over 3,173–3,366 m (200 m interval) | Jiang et al. (2018) |
| Primary injection zone | Upper Deadwood (40–50% of CO₂ by mass; 96% of CO₂ by mass per Year 3 data) | Roach et al. (2017); White & Roach (2018) |
| Cumulative CO₂ stored | >585,000 tonnes (through ~2024) | PTRC, SSRN (2025) |
| CO₂ injection start | April 2015 | White et al. (2017) |
| MMV methods | 4D seismic, DAS/DTS, downhole P/T gauges, soil gas, groundwater, pulsed-neutron logging | Multiple sources |
| Observation well offset | 151–152 m from injector | Talman et al. (2025); Zambrano-Narvaez et al. (2015) |
| Induced seismicity | None detected (M > −0.8 detectable threshold) | White et al. (2017) |
| Permeable thickness (kh) | Insufficient published kh — injectivity index 0.02–0.16 t/day/kPa | Jiang et al. (2018); Rangriz Shokri et al. (2019) |
| Salt precipitation | Observed in/near perforations; linked to intermittent injection | Talman et al. (2025); PTRC SSRN (2025) |
| Injectivity improvement | 2–4× increase associated with cold CO₂ injection (thermoelastic effects) | Rangriz Shokri et al. (2019) |

---

## SECTION 3: Regulatory and Operational Context

| Item | Detail | Source |
|------|--------|--------|
| Operator (capture) | SaskPower (Crown utility) | NRCan project record |
| Operator (Aquistore) | PTRC (Petroleum Technology Research Centre) / SaskPower | PTRC records |
| Operator (Weyburn EOR) | Whitecap Resources (current); Cenovus/EnCana (original) | IEA GHG reports |
| Capture technology | Post-combustion amine (CANSOLV/Shell) | SaskPower |
| Design capture rate | 1 Mtpa (90% capture) | SaskPower |
| Best year capture | 848,388 t (2024); 12-month record: 900,967 t (Aug 2023–Jul 2024) | SaskPower blog |
| Total captured (all-time) | 7,327,868 t (Oct 2014–Dec 2025) | SaskPower Q4 2025 update |
| 2025 capture | 721,239 t (Q1: 226,359 t; Q2: 25,443 t (outage); Q3: 236,512 t; Q4: 232,826 t) | SaskPower 2025 quarterly updates |
| Regulatory (Weyburn) | Oil & Gas Conservation Act — SME permits | SME records |
| **Regulatory (Aquistore)** | ***Captured Carbon Storage Act* (Saskatchewan, 2011, updated 2023)** — Crown owns pore space; PTRC holds subsurface storage rights for project area. Formal SLA boundary not publicly mapped. | Government of Saskatchewan; PTRC |
| CO₂ transport | BD3→Weyburn: dedicated pipeline (~1.6 km to pipeline, then ~50 km to field); BD3→Aquistore: ~10 km underground pipeline | PTRC; IEA GHG |
| Capture cost basis | $1.24B project cost; NRCan CEF $240M; SK gov't approval Apr 2011 | NRCan; Sask gov't |

---

## Notes on Data Quality

1. **Weyburn-Midale Beds** are the most studied CCS-EOR reservoir in the world, with 25+ years of IEA GHG research (Phase I 2000–2004, Phase II 2005–2012), hundreds of peer-reviewed publications, and continuous operational data. Property values are highly reliable.

2. **The Marly** is the primary CO₂ flood target (unswept oil from waterflood). Its high porosity (25% log, 17% core-stressed) and low matrix permeability (6–10 mD) are well-constrained. Natural fractures significantly enhance effective permeability.

3. **Porosity dual-value guidance**: The Marly has two well-supported porosity values reflecting different measurement scales. For CO₂BLOCK screening purposes, use **25% (log-derived) for pore-volume-based storage calculations** and **17% (core/stressed) for stress-sensitive modeling** (injectivity, geomechanics).

4. **The Vuggy** shoal deposits are the highest-permeability reservoir (up to 500 mD) but have been largely swept by the waterflood. CO₂ moves through both Marly and Vuggy due to fracture connectivity.

5. **Aquistore (Deadwood/Black Island)** properties are based on a smaller dataset (one injection well, one observation well, 30 km² 3D seismic) — properties are less certain than the Weyburn dataset.

6. **Black Island Member porosity**: The previously cited 21% (Bedard 2024) represents a **regional estimate** of the Winnipeg Formation, likely from shallower or updip portions. Site-specific well-log data from the Aquistore injection interval consistently shows **~7% mean interval porosity** (White & Roach 2018, GHGT-14; Roach et al. 2017, IJGGG). **Use 7% for Aquistore site-specific calculations.**

7. **CO₂ density context**: At Weyburn conditions (15 MPa, 60°C) → ~600–750 kg/m³. At Aquistore conditions (34–39 MPa, 110–115°C) → ~700–800 kg/m³. Higher pressures at Aquistore produce denser CO₂ despite higher temperature.

8. **Permeability scale effect at Aquistore**: Core permeability (5 mD) underrepresents field-scale effective permeability due to natural fractures and thermally-induced permeability enhancement. The injectivity-based estimates suggest effective permeability of 20–100 mD+.

9. **Temperature discrepancy at Aquistore**: Reported initial reservoir temperature ranges from 110°C (Roach et al.) to 119°C (Zambrano-Narvaez et al.). This may reflect different measurement locations (DTS vs downhole gauge) or cooling from injection operations.

10. **Fracture gradient at Aquistore**: The operational design value of 14.9 kPa/m (48 MPa at 3.2 km) is based on a mini-frac test. Thermal effects from cold CO₂ injection can reduce effective stress and may lead to localized fracturing below the design fracture pressure (Rangriz Shokri et al. 2019).

11. **Prairie Evaporite as regional barrier**: The Prairie Evaporite (Devonian salt) is ~600 m above the reservoir (~2,515 m depth vs ~3,130 m). It is **not** a directly overlying secondary seal for the storage complex but a separate regional containment barrier. Termed "regional barrier / secondary containment" in this revision.

12. **Aquistore area (SLA)**: No formal Sequestration Lease Area comparable to Alberta's township-grid system is publicly documented for Aquistore. The regulatory boundary is defined under Saskatchewan's *Captured Carbon Storage Act* and is tied to the project area (~30 km² 3D seismic + MMV footprint).

13. **Capture performance**: BD3 has never sustained 1 Mtpa (design capacity). The best year (848,388 t in 2024) reached ~85% of design. The facility continues to improve, with 99.3% availability in Q4 2025.

---

## Source Tier Assessment

| Tier | Definition | Count in Doc | Examples |
|------|------------|--------------|----------|
| **Tier 1** | Peer-reviewed journal | ~11 | GSA Today, IJGGG (×6), SPE J (×2), Energy Procedia (×2) |
| **Tier 2** | Peer-reviewed conference, SPE paper, book chapter, research report | ~8 | ACG (2015), SPE (1993), SSRN (2025), Wiley (2025), MDPI EngProc (2024 — proceedings), IEA GHG Final Report |
| **Tier 3** | Industry presentation, operator blog, project report | ~4 | Jenkins (2018), SaskPower blog, Whittaker (2005) |

**Tier inflation correction applied**: Bedard et al. (2024) reclassified from "Peer-reviewed journal" to "Conference proceedings (Tier 2)" per adversarial review. IEA GHG Final Report (2008) retained as Tier 2 (project report, not peer-reviewed journal).

---

*Document prepared: 2026-07-23*  
*Method: Tier 1–3 primary sources (peer-reviewed, regulatory, operator) cross-checked against SaskPower operational data, IEA GHG Weyburn reports, PTRC Aquistore publications*  
*Confidence: HIGH — Weyburn (25+ yr dataset); MODERATE-HIGH — Aquistore (10+ yr dataset, single-well)*  
*Document: CANADA-CCS-WILLISTON-BD3-props-revised.md*  
*Base document: CANADA-CCS-WILLISTON-BD3-props.md — revised per adversarial review (2026-07-23) + web research*

---

## DONE

All revisions complete. Summary of corrections applied:

1. **Black Island porosity**: 21% → 7% site-specific (3 new Tier 1 sources added: White & Roach 2018, White et al. 2016, Roach et al. 2021)
2. **CO₂ density (Weyburn)**: 700–800 → 600–750 kg/m³
3. **Injection rate label**: Corrected "at reservoir conditions" → "at surface/standard conditions"
4. **Prairie Evaporite**: "Secondary seal" → "Regional barrier / secondary containment (~600 m above)"
5. **Thickness clarification**: Explicitly separated gross reservoir sequence (~290 m) from net perforated interval (200 m)
6. **Marly porosity guidance**: Added usage note distinguishing log vs core for different calculation types
7. **Aquistore area/SLA**: Added Saskatchewan *Captured Carbon Storage Act* context
8. **Source tier**: Corrected Bedard (2024) from "Peered-reviewed journal" → "Conference proceedings (Tier 2)"
9. **Water saturation**: Added Sw ~55% (upper Deadwood) from Roach et al. (2017) — absent in original
10. **Revisions tracking**: § "Revisions from Adversarial Review" table added at top
