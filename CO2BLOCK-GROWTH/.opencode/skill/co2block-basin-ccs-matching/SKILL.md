# Skill: Matching CO₂ Storage Sites to Basins via Primary/Academic Web Sources

**Purpose**: Independently verify or discover CCS project–basin correspondences using original government, regulatory, industry, and peer-reviewed sources — NOT pre-compiled databases (e.g., Evenick 2021, Global CCS Institute project lists, CO2BLOCK basin tables).

---

## When to Use This Skill
- User asks: "Find CCS projects in Basin X" and wants **traceable primary evidence**
- User needs to **validate or refute** a database-derived basin–project mapping
- Basin is **frontier/underexplored** (Arctic, deep offshore, international) where compiled databases are incomplete
- User requires **audit-ready citations** (regulatory filings, government program records, peer-reviewed storage assessments)

---

## Step-by-Step Procedure

### 1. Extract Target Basin Names Precisely
- Use the **exact basin name** from the geological dataset (e.g., Evenick `Basin Name` column)
- Note **synonyms/aliases** (e.g., "Sverdrup Basin" = "Sverdrup Basin Province" USGS; "Arctic Islands Basin" informal)
- Record **basin polygons/coordinates** if available for spatial cross-checks

### 2. Build Targeted Search Queries (One Basin Per Query)
| Query Pattern | Example |
|---------------|---------|
| `"<Basin Name>" "carbon capture" OR "carbon storage" OR "CCS" OR "CO2 storage"` | `"Sverdrup Basin" "carbon capture"` |
| `"<Basin Name>" "CO2" OR "CO₂" "storage" OR "sequestration" site:gov.ca OR site:gc.ca` | `"Sverdrup Basin" CO2 storage site:gc.ca` |
| `"<Basin Name>" "pore space" OR "sequestration lease" OR "storage complex"` | `"Sverdrup Basin" "pore space"` |
| `"<Basin Name>" "Geological Survey" "Open File" "CO2"` | `"Sverdrup Basin" "Geological Survey" "Open File" CO2` |
| `site:nrcan.gc.ca "<Basin Name>" CCS` | `site:nrcan.gc.ca "Sverdrup" CCS` |

**Rule**: One basin per search call. Do NOT combine multiple basin names in one query.

### 3. Prioritize Source Tiers (Highest → Lowest Credibility)
| Tier | Source Types | Weight |
|------|-------------|--------|
| **1. Regulatory Filings** | AER/ERCB Directive 065 applications, SK pore space leases, BCER CCS permits, CER project assessments | ★★★★★ |
| **2. Government Program Records** | NRCan Clean Energy Fund/CCS Fund awards, CCS Knowledge Centre project database, provincial CCS funding lists | ★★★★★ |
| **3. Geological Survey Assessments** | GSC Open Files (CO₂ storage atlases), USGS Professional Papers (storage resource assessments), provincial atlas publications | ★★★★☆ |
| **4. Peer-Reviewed Storage Studies** | *IJGGC*, *Marine and Petroleum Geology*, *AAPG Bulletin* — papers with volumetric CO₂ capacity estimates using SRMS/DOE methodology | ★★★★☆ |
| **5. Industry Project Announcements** | Corporate press releases, investor presentations, Pathways Alliance filings — **only if matched to regulatory Tier 1/2** | ★★★☆☆ |
| **6. Strategic Assessments / Hypotheticals** | Federal SEA/EIS documents (e.g., Beaufort SEA), academic scenario modeling — **flag as "not a project"** | ★★☆☆☆ |
| **7. Media / DAC Concepts** | News articles, startup websites (e.g., TerraFixing) — **capture only, no storage site** | ★☆☆☆☆ |

### 4. Evidence Extraction Checklist (Per Candidate Project)
For each potential project hit, extract and record:
- [ ] **Project name** (official regulatory name)
- [ ] **Proponent/operator**
- [ ] **Storage formation(s)** (stratigraphic unit, depth, lithology)
- [ ] **Basin assignment** (explicit in source, or derived from coordinates)
- [ ] **Status** (Operational / Construction / FID / FEED / Concept / Cancelled) — use **Global CCS Institute status taxonomy**
- [ ] **Capacity** (Mtpa design, total Mt)
- [ ] **Regulatory approval reference** (AER Scheme Approval #, SK Lease #, BCER Permit #, CER Order #)
- [ ] **Pore space tenure instrument** (Alberta: *Mines and Minerals Act* lease; SK: *Captured Carbon Storage Act* lease)
- [ ] **CO₂ source** (facility name, type, distance to injection)
- [ ] **Coordinates / well UWIs** (for spatial basin verification)
- [ ] **Citation** (URL, document title, date, page/section)

### 5. Cross-Validate Basin Assignment
- **Spatial check**: Plot project injection well coordinates (or facility location) against basin polygon (Evenick/GSC/USGS). Do NOT trust the source's basin label — verify geometrically.
- **Stratigraphic check**: Does the storage formation belong to the basin's stratigraphic column? (e.g., Basal Cambrian Sands = WCSB, not Sverdrup)
- **Jurisdictional check**: Provincial regulator (AER/SME/BCER) → WCSB. Federal/territorial → Arctic offshore.

### 6. Document Negative Findings Explicitly
If **no projects found** after Tier 1–4 search:
- State: "No Tier 1–4 evidence of CCS projects in Basin X"
- List search queries run and date
- Note any Tier 5–7 mentions with **"NOT A PROJECT"** flags (e.g., hypothetical scenario, DAC concept only)
- Record geological storage potential assessments (Tier 3) separately from project lists

---

## Caveats & Common Mistakes to Avoid

| Mistake | Why It's Wrong | Correct Approach |
|---------|----------------|------------------|
| **Trusting Global CCS Institute "basin" field** | GCGI assigns projects to *countries/regions*, not geological basins; often mislabels WCSB projects as "Alberta Basin" or "Western Canada" | Spatially verify each project against basin polygons |
| **Using Evenick 2021 basin table as project list** | Evenick is a *basin catalog*, not a *project database*; basin-level attributes (type, countries) ≠ project-level data | Use Evenick only for basin geometry/names; search projects separately |
| **Conflating "CCS mentioned in SEA/EIS" with "CCS project"** | Strategic assessments list CCS as *hypothetical mitigation*; no proponent, no lease, no schedule | Tag as "Scenario/Mitigation Only" — exclude from project counts |
| **Counting DAC (Direct Air Capture) as basin storage** | DAC captures CO₂ from air; storage may be 1000+ km away in different basin (e.g., TerraFixing Fermont → WCSB) | Separate **capture location** from **storage basin**; only storage basin matters for basin matching |
| **Assuming provincial CCS regulations apply to federal/Arctic areas** | AER/SME/BCER jurisdiction ends at provincial boundary; Arctic offshore = federal (CER, CIRNAC) + Inuit co-management | Check jurisdiction per project; no CCS regulatory regime exists for Sverdrup/Nunavut offshore |
| **Treating "gas field discovered" as "depleted field storage candidate"** | Discovered ≠ produced ≠ depleted. Sverdrup fields (Drake, Hecla) never produced; no well infrastructure | Require **production history + abandonment status** for depleted-field storage claims |
| **Ignoring Indigenous governance requirements** | Nunavut Agreement, NIRB screening, IIBA mandatory for any Arctic project; adds 5-10+ years | Flag as "Governance Barrier" in viability assessment |
| **Using outdated basin names** | "Arctic Islands Basin" (old) vs "Sverdrup Basin" (modern); "Beaufort Sea Basin" vs "Beaufort-Mackenzie Basin" | Use current GSC/USGS/CSPG nomenclature; note aliases in search |

---

## Output Template (Per Basin)

```markdown
# CCS Projects in [Basin Name] — Primary Source Audit

## Search Metadata
- Date: YYYY-MM-DD
- Queries: [list exact strings]
- Tiers searched: 1–4 (primary), 5–7 (supplementary)
- Spatial verification: [basin polygon source]

## Findings Summary
| Project | Status | Storage Formation | Basin Verified? | Tier | Key Citation |
|---------|--------|-------------------|-----------------|------|--------------|
| None    | —      | —                 | —               | —    | —            |

## Negative Result Statement
"No Tier 1–4 evidence of CCS projects in [Basin] as of [date]."

## Tier 5–7 Mentions (Non-Projects)
| Item | Type | Why Excluded |
|------|------|--------------|
| ...  | ...  | ...          |

## Geological Storage Potential (Tier 3 only)
| Assessment | Formation | Capacity Estimate | Methodology | Citation |
|------------|-----------|-------------------|-------------|----------|
| ...        | ...       | ...               | ...         | ...      |

## Barriers to Development
| Category | Details |
|----------|---------|
| Emission Sources | ... |
| Infrastructure | ... |
| Regulation | ... |
| Governance | ... |
| Economics | ... |
```

---

## Example Application: Sverdrup Basin (2025-07-23)
See `CANADA-SVERSDRUP-CCS.md` — zero projects found; only hypothetical SEA mention (Beaufort, not Sverdrup) and DAC concept (TerraFixing, storage unspecified).

---

## Maintenance Notes
- Re-run Tier 1–2 searches **quarterly** (regulatory filings update)
- Re-run Tier 3 searches **annually** (GSC Open Files, USGS assessments)
- Update jurisdictional framework table when federal Arctic CCS regulation enacted
- Archive all source PDFs/HTML locally (regulatory sites change URLs)