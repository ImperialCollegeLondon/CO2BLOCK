#!/usr/bin/env python3
"""Compile CANADA-CCS.csv from extracted project data."""

import csv

# Define CSV columns
FIELDS = [
    "site_name", "basin_name", "basin_classification", "basin_evidence", "basin_evidence_url",
    "project_type", "target_formation", "status", "operator",
    "depth_m", "depth_source_url",
    "porosity_pct", "porosity_source_url",
    "permeability_md", "permeability_source_url",
    "net_pay_m", "net_pay_source_url",
    "temperature_c", "temperature_source_url",
    "pressure_mpa", "pressure_source_url",
    "storage_capacity", "storage_capacity_source_url",
    "area_km2", "area_source_url",
    "primary_references"
]

def s(val, url):
    """Format a property value with its source URL."""
    if url:
        return f"{val} [{url}]"
    return val

def link_list(urls_str):
    """Normalize semicolon-separated URLs into a readable list."""
    if not urls_str:
        return ""
    parts = [u.strip() for u in urls_str.split(";") if u.strip()]
    return "; ".join(parts)

rows = []

# ============ 1. QUEST CCS ============
rows.append({
    "site_name": "Quest CCS (Shell)",
    "basin_name": "Alberta Basin (Western Canada Sedimentary Basin)",
    "basin_classification": "Foreland basin",
    "basin_evidence": "Western Canada Sedimentary Basin — classic foreland basin formed by thrust loading during the Laramide orogeny; BCS is a Cambrian-Ordovician cratonic sandstone at the base of the sedimentary succession.",
    "basin_evidence_url": "https://ags.aer.ca/publications/all-publications/prs-2024-001",
    "project_type": "CCS — deep saline aquifer storage",
    "target_formation": "Basal Cambrian Sands (BCS), Cambrian-Ordovician saline aquifer",
    "status": "Operational (injection since Nov 2015)",
    "operator": "Shell Canada Products",
    "depth_m": s("~2,000 m (range 1,800-2,132 m MD)", "https://open.alberta.ca/publications/quest-carbon-capture-and-storage-project-d-65-application"),
    "porosity_pct": s("17% (high-quality interval); 14% full-formation avg; core range 8-21%", "https://doi.org/10.1016/j.egypro.2017.03.1630"),
    "permeability_md": s("20-500 mD (D-65 well-test range); ~1,000 mD (optimistic end-member, single conference source)", "https://open.alberta.ca/publications/quest-carbon-capture-and-storage-project-d-65-application"),
    "net_pay_m": s("30 m (perforated injection interval); gross 40-45 m; net sand ~38 m", "https://open.alberta.ca/publications/quest-carbon-capture-and-storage-project-gen-4-reservoir-modeling-report"),
    "temperature_c": s("~60°C", "https://open.alberta.ca/publications/quest-carbon-capture-and-storage-project-d-65-application"),
    "pressure_mpa": s("~21 MPa (initial hydrostatic)", "https://open.alberta.ca/publications/quest-carbon-capture-and-storage-project-d-65-application"),
    "storage_capacity": s(">7 Mt cumulative (as of 2023); target 25 Mt over 25 years", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "area_km2": s("SLA ~3,667 km² (906,240 acres); AoR 461 km²; CO₂ plume ~29.3 km²", "https://static.aer.ca/prd/documents/decisions/2012/2012-ABERCB-008.pdf"),
    "primary_references": link_list(
        "https://open.alberta.ca/publications/quest-carbon-capture-and-storage-project-d-65-application;"
        "https://open.alberta.ca/publications/quest-carbon-capture-and-storage-project-gen-4-reservoir-modeling-report;"
        "https://ags.aer.ca/publications/all-publications/prs-2024-001;"
        "https://doi.org/10.1016/j.egypro.2017.03.1630;"
        "https://doi.org/10.1016/j.ijggc.2022.103718;"
        "https://static.aer.ca/prd/documents/decisions/2012/2012-ABERCB-008.pdf;"
        "https://www.searchanddiscovery.com/documents/2017/80577rock/ndx_rock.pdf"
    )
})

# ============ 2. ACTL / Clive CCS ============
rows.append({
    "site_name": "ACTL / Clive CCS (Enhance Energy)",
    "basin_name": "Alberta Basin (Western Canada Sedimentary Basin)",
    "basin_classification": "Foreland basin",
    "basin_evidence": "Located on the Bashaw Reef Complex of the Rimbey-Meadowbrook reef chain, within the foreland-basin fill of the WCSB. The Leduc Formation is a Late Devonian carbonate platform developed on the Peace River Arch architecture.",
    "basin_evidence_url": "https://cdnsciencepub.com/doi/10.1139/cjes-2012-0137",
    "project_type": "CO₂-EOR + associated saline aquifer storage",
    "target_formation": "Leduc Formation D-3A Pool (Upper Devonian, Woodbend Group), Bashaw Reef Complex",
    "status": "Operational (ACTL pipeline 2020; Clive CO₂-EOR ongoing)",
    "operator": "Enhance Energy Inc.",
    "depth_m": s("1,845-1,938 m TVD (mean 1,874 m)", "https://open.alberta.ca/publications/actl-dow-demonstration-summary-report"),
    "porosity_pct": s("6.1% (Enhance geo-model avg); regional Leduc range 5-20%; arithmetic mean ~8%", "https://doi.org/10.1306/a25ff215-171b-11d7-8645000102c1865d"),
    "permeability_md": s("30 mD (mode); regional arithmetic mean ~26 mD; range 1-1,000 mD", "https://doi.org/10.1139/cjes-2012-0137"),
    "net_pay_m": s("8.75 m (perforated/net-pay); reservoir interval 35 m for storage", "https://open.alberta.ca/publications/actl-dow-demonstration-summary-report"),
    "temperature_c": s("~69°C", "https://geoconvention.com/2025/"),
    "pressure_mpa": s("~13 MPa (~13,000 kPag initial, modified by production)", "https://geoconvention.com/2025/"),
    "storage_capacity": s("18.8 Mt total; ~4-7 Mt injected (2024-25); ~14 Mt remaining", "https://cleanprosperity.ca"),
    "area_km2": s("Pool ~25-40 km²; SLA 544,512 ha (~5,445 km² for Origins hub)", "https://keyfactsenergy.com"),
    "primary_references": link_list(
        "https://open.alberta.ca/publications/actl-dow-demonstration-summary-report;"
        "https://doi.org/10.1306/a25ff215-171b-11d7-8645000102c1865d;"
        "https://doi.org/10.1139/cjes-2012-0137;"
        "https://cleanprosperity.ca;"
        "https://geoconvention.com/2025/;"
        "https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5068925;"
        "https://open.alberta.ca/dataset/90f61413-0ef1-45a4-9e1c-6bff7c23fd7e/resource/99fae078-1316-432b-8447-acd74bd1abb0/download/energy-actl-enhance-energy-clive-leduc-field-mmv-plan-2019-07.pdf"
    )
})

# ============ 3. NWR Sturgeon Refinery ============
rows.append({
    "site_name": "NWR Sturgeon Refinery CO₂ Capture",
    "basin_name": "Alberta Basin (Western Canada Sedimentary Basin)",
    "basin_classification": "Foreland basin (capture facility — storage via ACTL pipeline to Clive hub)",
    "basin_evidence": "Sturgeon County ~45 km NE of Edmonton, Alberta, within the WCSB foreland basin. The CO₂ captured here is transported via the ACTL pipeline to the Clive Leduc D-3A pool for CO₂-EOR and storage.",
    "basin_evidence_url": "https://nwrsturgeonrefinery.com/project/carbon-capture-and-storage/",
    "project_type": "CO₂ capture-only (captured CO₂ transported via ACTL for EOR at Clive)",
    "target_formation": "N/A (capture facility — CO₂ injected into Leduc D-3A at Clive)",
    "status": "Operational (capture started 2020; refinery operational since 2017)",
    "operator": "North West Redwater Partnership (NWRP) — 50/50 JV: CNRL + APMC",
    "depth_m": "N/A (capture facility)",
    "depth_source_url": "",
    "porosity_pct": "N/A (capture facility)",
    "porosity_source_url": "",
    "permeability_md": "N/A (capture facility)",
    "permeability_source_url": "",
    "net_pay_m": "N/A (capture facility)",
    "net_pay_source_url": "",
    "temperature_c": "N/A (capture facility)",
    "temperature_source_url": "",
    "pressure_mpa": "N/A (capture facility)",
    "pressure_source_url": "",
    "storage_capacity": s("Design capture 1.2-1.3 Mtpa; cumulative >5 Mt captured (2020-2024)", 
        "https://nwrsturgeonrefinery.com/assets/uploads/2025/06/NWRP_MISSIONMATTERS_REPORT_2024_June-2_FINAL.pdf"),
    "area_km2": "N/A (capture facility)",
    "area_source_url": "",
    "primary_references": link_list(
        "https://nwrsturgeonrefinery.com/project/carbon-capture-and-storage/;"
        "https://nwrsturgeonrefinery.com/assets/uploads/2025/06/NWRP_MISSIONMATTERS_REPORT_2024_June-2_FINAL.pdf;"
        "https://www.catf.us/resource/carbon-capture-storage-what-can-learn-from-project-track-record/;"
        "https://www.apmc.ca/public/download/files/247620;"
        "https://ccsknowledge.com/technical-studies/the-alberta-carbon-trunk-line-case-study/;"
        "https://open.alberta.ca/dataset/51dd6f66-30a0-4ae9-908f-db3730172fed/resource/2fa85ef5-6612-4dd5-b710-89807049263a/download/energy-actl-knowledge-sharing-2021-summary-report.pdf"
    )
})

# ============ 4. Polaris + Atlas Hub ============
rows.append({
    "site_name": "Shell Polaris + Atlas Hub CCS",
    "basin_name": "Alberta Basin (Western Canada Sedimentary Basin)",
    "basin_classification": "Foreland basin",
    "basin_evidence": "BCS saline aquifer in the WCSB foreland basin, gentle regional SW dip (~5-10 m/km); located in the industrial heartland region of Alberta near Fort Saskatchewan.",
    "basin_evidence_url": "https://ags.aer.ca/publications/all-publications/prs-2024-001",
    "project_type": "CCS — deep saline aquifer storage",
    "target_formation": "Basal Cambrian Sands (BCS), Cambrian-Ordovician",
    "status": "Under Construction (FID June 26, 2024; operations expected end-2028)",
    "operator": "Shell Canada Products (Polaris) / Shell-ATCO EnPower JV (Atlas)",
    "depth_m": s("~2,000 m (range 1,800-2,100 m)", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "porosity_pct": s("~10% (regional avg for Edmonton area); range <1% to >25%", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "permeability_md": s("1-2 mD (Edmonton-region avg, Hofmann 2013); basin-wide avg ~10 mD (Weides 2014)", "https://pangea.stanford.edu/ERE/pdf/IGAstandard/SGW/2013/Hofmann.pdf"),
    "net_pay_m": s("~45 m (P50 estimate); regional range 17-115 m", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "temperature_c": s("~76°C (site-specific estimate); regional range 65-120°C", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "pressure_mpa": s("~20 MPa (estimated hydrostatic)", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "storage_capacity": s("300 Mt (lifetime claimed, Contingent Resource); Phase 1 contracted ~20 Mt", "https://open.alberta.ca/publications/"),
    "area_km2": s("10,312 km² (1,031,232 ha / 111 townships — CSA lease area)", "https://open.alberta.ca/publications/"),
    "primary_references": link_list(
        "https://ags.aer.ca/publications/all-publications/prs-2024-001;"
        "https://ags.aer.ca/publications/all-publications/prs-2023-001;"
        "https://pangea.stanford.edu/ERE/pdf/IGAstandard/SGW/2013/Hofmann.pdf;"
        "https://doi.org/10.1139/cjes-2012-0137;"
        "https://open.alberta.ca/publications/"
    )
})

# ============ 5. Wolf Lamont Carbon Hub ============
rows.append({
    "site_name": "Wolf Lamont Carbon Hub",
    "basin_name": "Alberta Basin (Western Canada Sedimentary Basin)",
    "basin_classification": "Foreland basin",
    "basin_evidence": "BCS in the Lloydminster Embayment of the WCSB foreland basin; Cambrian-Ordovician cratonic siliciclastic transgressive systems tract.",
    "basin_evidence_url": "https://ags.aer.ca/publications/all-publications/prs-2024-001",
    "project_type": "CCS — deep saline aquifer storage",
    "target_formation": "Basal Cambrian Sandstone (BCS), Cambrian-Ordovician",
    "status": "AER Approved 10-Mar-2026 (Approval 13513, Directive 065); pre-injection",
    "operator": "Wolf Midstream / Wolf Carbon Solutions / Whitecap Resources / FNCIP / Heart Lake First Nation",
    "depth_m": s("~1,800-2,200 m (~5,900-7,200 ft)", "https://dds.aer.ca/iar_query/viewApplication.aspx?appnumber=1959125"),
    "porosity_pct": s("~17% (formation-scale avg, Quest analog); regional variance 11-19%", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "permeability_md": s("~1,000 mD (formation-scale avg, Quest analog); range 1 mD to >1,000 mD", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "net_pay_m": s("40-50 m (gross thickness inferred); regional BCS range 40-80 m", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "temperature_c": s("65-80°C (expected at Lamont depth ~2 km)", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "pressure_mpa": s(">7.38 MPa (supercritical threshold confirmed); specific formation pressure not published", ""),
    "storage_capacity": s("19.2 Gt (Alberta BCS total, 2% efficiency factor); initial injection 2-3 Mtpa", "https://albertainnovates.ca/wp-content/uploads/2024/03/CCUS-Whitepaper-Supplement-2-Carbon-Storage.pdf"),
    "area_km2": s("~839 km² (~324 mi² / ~9 townships — AER application boundary)", "https://dds.aer.ca/iar_query/viewApplication.aspx?appnumber=1959125"),
    "primary_references": link_list(
        "https://dds.aer.ca/iar_query/viewApplication.aspx?appnumber=1959125;"
        "https://ags.aer.ca/publications/all-publications/prs-2024-001;"
        "https://ags.aer.ca/publications/all-publications/prs-2023-001;"
        "https://wolfmidstream.com/carbon/;"
        "https://albertainnovates.ca/wp-content/uploads/2024/03/CCUS-Whitepaper-Supplement-2-Carbon-Storage.pdf"
    )
})

# ============ 6. Wabamun Carbon Hub ============
rows.append({
    "site_name": "Wabamun Carbon Hub (Enbridge)",
    "basin_name": "Alberta Basin (Western Canada Sedimentary Basin)",
    "basin_classification": "Foreland basin",
    "basin_evidence": "WCSB foreland basin; targeting BSU (Cambrian-Ordovician) + Nisku Formation (Devonian) in the Wabamun Lake area, central Alberta.",
    "basin_evidence_url": "https://ags.aer.ca/publications/all-publications/prs-2024-001",
    "project_type": "CCS — deep saline aquifer storage (BSU primary; Nisku secondary)",
    "target_formation": "Basal Sandstone Unit (BSU) — Cambrian-Ordovician deep saline aquifer (primary); Nisku Formation — Devonian carbonate (secondary)",
    "status": "CSEA awarded March 2022; CSA signed October 2025; stratigraphic test wells drilled 2023-2024; pre-injection",
    "operator": "Enbridge Inc. (through Enbridge Wabamun Holdings Inc.)",
    "depth_m": s("2,300-3,300 m (top BSU); Nisku secondary 1,550-2,200 m", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "porosity_pct": s("8-12% (central value 10%, Weides 2014); Alberta Innovates regional 14.4%", "https://doi.org/10.1139/cjes-2012-0137"),
    "permeability_md": s("Low (matrix) ~0.01 mD; Mid (effective) ~0.1-1.6 mD; High (post-stim) ~1-10 mD (speculative)", "https://doi.org/10.1139/cjes-2012-0137"),
    "net_pay_m": s("~40 m (regional midpoint); BSU range 30-50 m for Wabamun area", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "temperature_c": s("~80-115°C (expected at Wabamun); geothermal gradient 35.6°C/km", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "pressure_mpa": s("~23-33 MPa (hydrostatic estimate at 2,300-3,300 m); may be overpressured", ""),
    "storage_capacity": s("Up to 4 Mt/yr total target (Heidelberg 780 kt/yr + Capital Power 3 Mt/yr discontinued May 2024); no public P50/P90 estimate", "https://www.nrcan.gc.ca/climate-change/canadas-green-future/carbon-capture-utilization-and-storage"),
    "area_km2": s("~1,500 km² (CSA North Area, BSU); full tenure ~2,971 km²", ""),
    "primary_references": link_list(
        "https://ags.aer.ca/publications/all-publications/prs-2024-001;"
        "https://doi.org/10.1139/cjes-2012-0137;"
        "https://www.nrcan.gc.ca/climate-change/canadas-green-future/carbon-capture-utilization-and-storage"
    )
})

# ============ 7. Pathways Alliance ============
rows.append({
    "site_name": "Pathways Alliance CCS Hub",
    "basin_name": "Alberta Basin (Western Canada Sedimentary Basin)",
    "basin_classification": "Foreland basin",
    "basin_evidence": "WCSB foreland basin formed by Laramide orogeny; hub located in northeastern Alberta (Cold Lake/Fort McMurray region) targeting the BCS saline aquifer.",
    "basin_evidence_url": "https://ags.aer.ca/publications/all-publications/prs-2024-001",
    "project_type": "CCS — deep saline aquifer storage (BCS primary); secondary Mannville Group",
    "target_formation": "Basal Cambrian Sandstone (BCS) — Cambrian-Ordovician; secondary: McMurray Fm / Clearwater Fm (Lower Cretaceous)",
    "status": "Pre-FID; Tripartite MOU signed July 2026 (federal-provincial-industry); regulatory application filed Q1 2024; IAAC suspended Dec 2024",
    "operator": "Canadian Natural Resources Ltd. (lead), Cenovus, ConocoPhillips, Imperial Oil, Suncor",
    "depth_m": s("1,000-1,500 m (Cold Lake area)", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "porosity_pct": s("N/A (not publicly available for Cold Lake hub; Quest analogue ~17%)", "https://doi.org/10.1016/j.ijggc.2022.103718"),
    "permeability_md": s("N/A (not publicly available; estimated 1-10 mD based on AGS diagenetic framework)", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "net_pay_m": s("40-80 m gross (BCS regional range)", "https://ags.aer.ca/publications/all-publications/prs-2024-001"),
    "temperature_c": s("N/A (site-specific not public; regional gradient applies)", ""),
    "pressure_mpa": s("N/A (site-specific not public; hydrostatic to slightly overpressured regionally)", ""),
    "storage_capacity": s("MOU Phase 1: 6 Mtpa by 2035; total target: 16 Mtpa net by 2045; full build-out: up to 40 Mtpa by 2050", "https://www.cbc.ca/news/canada/edmonton/alberta-ottawa-pathways-oilsands-carbon-capture-9.7268392"),
    "area_km2": s("N/A (not publicly delineated for hub; BCS regional extent ~250,000 km²)", "https://cdl.canadiandiscovery.com/study/regional-geological-ccs-assessment-of-the-basal-cambrian-sandstone"),
    "primary_references": link_list(
        "https://ags.aer.ca/publications/all-publications/prs-2024-001;"
        "https://ags.aer.ca/publications/all-publications/prs-2023-001;"
        "https://ags.aer.ca/publications/all-publications/prs-2022-001;"
        "https://doi.org/10.1016/j.ijggc.2022.103718;"
        "https://majorprojects.alberta.ca/details/Pathways-Alliance-Carbon-Capture-Storage-Hub-Phase-1/10695;"
        "https://iaac-aeic.gc.ca/050/evaluations/proj/89090;"
        "https://www.cbc.ca/news/canada/edmonton/alberta-ottawa-pathways-oilsands-carbon-capture-9.7268392;"
        "https://oilsandsalliance.ca/pathways-project/;"
        "https://cdl.canadiandiscovery.com/study/regional-geological-ccs-assessment-of-the-basal-cambrian-sandstone"
    )
})

# ============ 8. Origins CCS Hub ============
rows.append({
    "site_name": "Origins CCS Hub (Enhance Energy)",
    "basin_name": "Alberta Basin (Western Canada Sedimentary Basin)",
    "basin_classification": "Foreland basin",
    "basin_evidence": "Located in central Alberta within the WCSB foreland basin; targeting Leduc Formation carbonate reefs of the Rimbey-Meadowbrook reef trend (Bashaw Reef Complex).",
    "basin_evidence_url": "https://cdnsciencepub.com/doi/10.1139/cjes-2012-0137",
    "project_type": "CCS — pressure-depleted saline aquifer in dolomitized carbonate reef",
    "target_formation": "Leduc Formation (Upper Devonian, Woodbend Group), Lacombe Field",
    "status": "AER Approval No. 13463 issued 30 July 2025; pre-injection",
    "operator": "Enhance Energy Inc.",
    "depth_m": s("2,000 m (6,561 ft TVD)", "https://www.carbonstorage.io/storage/origins-project"),
    "porosity_pct": s("6% (well avg)", "https://www.carbonstorage.io/storage/origins-project"),
    "permeability_md": s("<10 to >100 mD (well range); regional Leduc avg ~26 mD", "https://www.carbonstorage.io/storage/origins-project"),
    "net_pay_m": s("137 m gross (450 ft)", "https://www.carbonstorage.io/storage/origins-project"),
    "temperature_c": s("65°C (149°F)", "https://www.carbonstorage.io/storage/origins-project"),
    "pressure_mpa": s("15.3 MPa (2,219 psia initial); max allowable 20 MPa @ 2,070 m TVD", "https://www.carbonstorage.io/storage/origins-project"),
    "storage_capacity": s("Phase 1: 28 Mt (1.6 Mt/yr × 17.5 yr); full-scale: up to 20 Mt/yr, several hundred Mt total", "https://enhanceenergy.com/origins-ccs-hub-receives-approval/"),
    "area_km2": s("Plume: ~50 km² nominal (4 km radius), ~99 km² max (5.6 km radius); lease: ~10,082 km²", "https://static.aer.ca/prd/documents/decisions/regulatory-appeal-decisions/1959099-20260414.pdf"),
    "primary_references": link_list(
        "https://www.carbonstorage.io/storage/origins-project;"
        "https://enhanceenergy.com/origins-ccs-hub-receives-approval/;"
        "https://enhanceenergy.com/enhance-announces-origins-project;"
        "https://static.aer.ca/prd/documents/decisions/regulatory-appeal-decisions/1959099-20260414.pdf;"
        "https://cdnsciencepub.com/doi/10.1139/cjes-2012-0137;"
        "https://doi.org/10.1306/a25ff215-171b-11d7-8645000102c1865d;"
        "https://majorprojects.alberta.ca/details/Origins-Carbon-Sequestration-Hub-Project/8539"
    )
})

# ============ 9. BD3 (Boundary Dam) ============
rows.append({
    "site_name": "Boundary Dam Integrated CCS (BD3)",
    "basin_name": "Williston Basin",
    "basin_classification": "Intracratonic basin (sag basin)",
    "basin_evidence": "The Williston Basin is a structural basin overlying the cratonic platform of North America, with Paleozoic carbonate and clastic reservoirs. BD3 captures CO₂ from a coal-fired power plant near Estevan, Saskatchewan.",
    "basin_evidence_url": "https://rock.geosociety.org/net/gsatoday/archive/14/7/pdf/gt0407.pdf",
    "project_type": "CO₂ capture + EOR (Weyburn-Midale) + saline aquifer injection (Aquistore)",
    "target_formation": "CO₂ sent to: Midale Beds (~1,450 m, Mississippian) and Deadwood Fm / Black Island Mb (~3,130 m, Cambro-Ordovician)",
    "status": "Operational (capture since Oct 2014; >7.3 Mt captured through Dec 2025)",
    "operator": "SaskPower (capture); Whitecap Resources / PTRC (storage)",
    "depth_m": s("N/A (capture facility; injection at Weyburn ~1,450 m and Aquistore ~3,130 m)", ""),
    "porosity_pct": s("N/A (capture facility; see Aquistore and Weyburn rows for storage properties)", ""),
    "permeability_md": s("N/A (capture facility; see Aquistore and Weyburn rows)", ""),
    "net_pay_m": s("N/A (capture facility)", ""),
    "temperature_c": s("N/A (capture facility)", ""),
    "pressure_mpa": s("N/A (capture facility)", ""),
    "storage_capacity": s("7,327,868 t cumulative CO₂ captured (Oct 2014-Dec 2025); design 1 Mtpa (de-rated from 1.6 Mtpa)", "https://ptrc.ca/aquistore"),
    "area_km2": s("N/A (capture facility)", ""),
    "primary_references": link_list(
        "https://ptrc.ca/aquistore;"
        "https://www.saskpower.com/projects/carbon-capture-and-storage"
    )
})

# ============ 10. Aquistore ============
rows.append({
    "site_name": "Aquistore CCS Research Project",
    "basin_name": "Williston Basin",
    "basin_classification": "Intracratonic basin",
    "basin_evidence": "Williston Basin — intracratonic, onshore setting; targeting the Deadwood Formation / Black Island Member (Cambro-Ordovician) at ~3.1-3.4 km depth near Estevan, SK.",
    "basin_evidence_url": "https://doi.org/10.1016/j.egypro.2014.11.320",
    "project_type": "CCS — deep saline aquifer storage (NOT EOR; research/demo)",
    "target_formation": "Deadwood Formation / Black Island Member (Cambro-Ordovician saline aquifer)",
    "status": "Operational Research/Demo (injection Apr 2015-ongoing; >585,000 t stored)",
    "operator": "Petroleum Technology Research Centre (PTRC) / SaskPower",
    "depth_m": s("3,130 m (reservoir top); gross reservoir ~290 m", "https://doi.org/10.1190/geo2016-0488.1"),
    "porosity_pct": s("6-11% site-specific; best estimate 7% (log-based, injection interval); 10-15% in perforated zones", "https://doi.org/10.1190/geo2016-0488.1"),
    "permeability_md": s("5 mD (core plug avg); ~70-88 mD (effective reservoir-scale); 20-150 mD field-scale effective", "https://doi.org/10.1190/geo2016-0488.1"),
    "net_pay_m": s("~112 m (51% NTG of 219 m mean thickness); perforated interval ~200 m", "https://geoconvention.com/wp-content/uploads/abstracts/2021/67360-closing-the-loop.pdf"),
    "temperature_c": s("119-120°C (far-field); 112.8°C at 3,173 m", "https://papers.acg.uwa.edu.au/d/1508_52_Zambrano-Narvaez/"),
    "pressure_mpa": s("34.1-35 MPa (initial reservoir); fracture ~48 MPa; BHP <42 MPa", "https://papers.acg.uwa.edu.au/d/1508_52_Zambrano-Narvaez/"),
    "storage_capacity": s("8.4-27.1 Mt (P10-P90 static, 34 km² model area); >585,000 t stored through ~2024", "https://doi.org/10.2172/1874344"),
    "area_km2": s("34 km² (fine-scale simulation model / 3D seismic survey area)", "https://doi.org/10.2172/1874434"),
    "primary_references": link_list(
        "https://doi.org/10.1016/j.egypro.2014.11.320;"
        "https://doi.org/10.1016/j.ijggc.2016.10.001;"
        "https://doi.org/10.1016/j.ijggc.2018.02.009;"
        "https://doi.org/10.1190/geo2016-0488.1;"
        "https://doi.org/10.2118/196118-ms;"
        "https://doi.org/10.2172/1874344;"
        "https://doi.org/10.2172/1874434;"
        "https://papers.acg.uwa.edu.au/d/1508_52_Zambrano-Narvaez/;"
        "https://ptrc.ca/aquistore;"
        "https://geoconvention.com/wp-content/uploads/abstracts/2021/67360-closing-the-loop.pdf"
    )
})

# ============ 11. Weyburn ============
rows.append({
    "site_name": "Weyburn CO₂-EOR Project",
    "basin_name": "Williston Basin",
    "basin_classification": "Intracratonic basin",
    "basin_evidence": "Williston Basin — intracratonic, onshore setting; the Weyburn field is located on the northern margin of the Williston Basin in southeastern Saskatchewan, producing from Mississippian Midale carbonates.",
    "basin_evidence_url": "https://rock.geosociety.org/net/gsatoday/archive/14/7/pdf/gt0407.pdf",
    "project_type": "CO₂-EOR (Enhanced Oil Recovery) with associated storage",
    "target_formation": "Mississippian Charles Formation, Midale Beds — Marly dolostone (upper) + Vuggy limestone (lower)",
    "status": "Active (injection since Sept 2000; ongoing under Cenovus)",
    "operator": "Cenovus Energy (Weyburn unit); co-funded by IEA GHG / PTRC",
    "depth_m": s("~1,450-1,500 m (mean ~1,450 m)", "https://rock.geosociety.org/net/gsatoday/archive/14/7/pdf/gt0407.pdf"),
    "porosity_pct": s("Marly 16-38% (avg ~26% log; 17% stressed core); Vuggy 8-21% (avg ~12%)", "https://rock.geosociety.org/net/gsatoday/archive/14/7/pdf/gt0407.pdf"),
    "permeability_md": s("Marly 1-50 mD (matrix, fracture-enhanced); Vuggy 10-300 mD (matrix, up to 500 mD with fractures)", "https://rock.geosociety.org/net/gsatoday/archive/14/7/pdf/gt0407.pdf"),
    "net_pay_m": s("<30 m (Marly ~6 m, Vuggy ~15-17 m; total net pay ~20-25 m)", "https://education.aapg.org/carbonsequestration/white.pdf"),
    "temperature_c": s("60-63°C", "https://rock.geosociety.org/net/gsatoday/archive/14/7/pdf/gt0407.pdf"),
    "pressure_mpa": s("~14 MPa (range 12.5-18 MPa)", "https://rock.geosociety.org/net/gsatoday/archive/14/7/pdf/gt0407.pdf"),
    "storage_capacity": s(">35 Mt net stored by 2023; projected ~40-55 Mt over project life; injection rate ~3.6 Mt/yr combined", "https://sequestration.mit.edu/tools/projects/weyburn.html"),
    "area_km2": s("~180 km² (Weyburn unit); ~210 km² (full structural closure); combined Weyburn+Midale ~284 km²", "https://fossil.energy.gov/archives/cslf/sites/default/files/documents/IEAGHGWeyburnProjectPoster0307.pdf"),
    "primary_references": link_list(
        "https://rock.geosociety.org/net/gsatoday/archive/14/7/pdf/gt0407.pdf;"
        "https://www.sciencedirect.com/science/article/pii/S1876610211008915;"
        "https://fossil.energy.gov/archives/cslf/sites/default/files/documents/IEAGHGWeyburnProjectPoster0307.pdf;"
        "https://sequestration.mit.edu/tools/projects/weyburn.html;"
        "https://pubsaskdev.blob.core.windows.net/pubsask-prod/88854/88854-Whittaker_2005vol1.pdf;"
        "https://education.aapg.org/carbonsequestration/white.pdf;"
        "https://harvest.usask.ca/items/721c70e8-061a-42c9-b9b0-0f5070a0e9d5;"
        "https://www.osti.gov/servlets/purl/1016981;"
        "https://www.sciencedirect.com/science/article/pii/S1750583612003386;"
        "https://ptrc.ca/past-projects/weyburn-midale"
    )
})

# ============ 12. Midale ============
rows.append({
    "site_name": "Midale CO₂-EOR Project",
    "basin_name": "Williston Basin",
    "basin_classification": "Intracratonic basin",
    "basin_evidence": "Williston Basin — intracratonic, onshore setting; Midale field is ~10 km east of Weyburn along the northern margin of the Williston Basin, producing from the same Midale Beds formation.",
    "basin_evidence_url": "https://rock.geosociety.org/gsatoday/archive/14/7/pdf/i1052-5173-14-7-4.pdf",
    "project_type": "CO₂-EOR and associated storage",
    "target_formation": "Mississippian Charles Formation, Midale Beds — Marly dolostone + Vuggy limestone",
    "status": "Active (injection since Sept 2005; ongoing under Cardinal Energy Ltd., 77% WI)",
    "operator": "Cardinal Energy Ltd. (acquired Apache's Midale Unit ~2021)",
    "depth_m": s("~1,350-1,500 m (reference 1,440 m)", "https://rock.geosociety.org/gsatoday/archive/14/7/pdf/i1052-5173-14-7-4.pdf"),
    "porosity_pct": s("Marly 16-38% (avg ~25-26%); Vuggy 2-20% (avg ~10-12%)", "https://rock.geosociety.org/gsatoday/archive/14/7/pdf/i1052-5173-14-7-4.pdf"),
    "permeability_md": s("Marly 1-100 mD (stressed ~6-10 mD); Vuggy <0.01-500 mD (shoal ~20 mD, intershoal ~1 mD)", "https://www.co2conference.net/wp-content/uploads/2018/12/Th8-Update-on-the-Weyburn-Project-in-Canada-Dec-6-2018.pdf"),
    "net_pay_m": s("~10-12 m total (Marly ~4 m, Vuggy ~8 m)", "https://www.co2conference.net/wp-content/uploads/2018/12/Th8-Update-on-the-Weyburn-Project-in-Canada-Dec-6-2018.pdf"),
    "temperature_c": s("63°C", "https://rock.geosociety.org/gsatoday/archive/14/7/pdf/i1052-5173-14-7-4.pdf"),
    "pressure_mpa": s("12.5-18 MPa (avg ~15 MPa; initial ~14 MPa)", "https://rock.geosociety.org/gsatoday/archive/14/7/pdf/i1052-5173-14-7-4.pdf"),
    "storage_capacity": s(">3 Mt injected by 2012; forecast 10 Mt (gross) over project life; cumulative 44+ Mt combined with Weyburn (as of 2025)", "https://natural-resources.canada.ca/energy/publications/16459"),
    "area_km2": s("~125 km² (31,000 acres unitized); ~260 km² (total field outline)", "https://doi.org/10.2118/22947-ms"),
    "primary_references": link_list(
        "https://rock.geosociety.org/gsatoday/archive/14/7/pdf/i1052-5173-14-7-4.pdf;"
        "https://doi.org/10.2118/22947-ms;"
        "https://doi.org/10.2118/15635-pa;"
        "https://doi.org/10.2118/22946-pa;"
        "https://www.osti.gov/biblio/1016981;"
        "https://natural-resources.canada.ca/energy/publications/16459;"
        "https://sequestration.mit.edu/tools/projects/weyburn.html;"
        "https://www.co2conference.net/wp-content/uploads/2012/12/1-2_Whittaker_Sask_CO2_Floods_2010.pdf;"
        "https://cardinalenergy.ca/operations/midale-sk/;"
        "https://www.wcap.ca/operations/core-areas/conventional-division"
    )
})

# Write CSV
output_path = "/Users/melizare/GitHub/CO2BLOCK/CO2BLOCK-GROWTH/CANADA-CCS.csv"
with open(output_path, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=FIELDS)
    w.writeheader()
    for r in rows:
        w.writerow(r)

print(f"CSV written: {output_path}")
print(f"Rows: {len(rows)}")
print(f"Columns: {len(FIELDS)}")
