# Skill: CO2BLOCK Missing Project Reservoir Properties Searcher

## Purpose
Given a list of CCS projects missing from the CO2BLOCK `Global.xlsx` database (from a `*-missing-revised.md` file), independently research and compile the four key reservoir properties — **thickness, porosity, permeability, and reservoir/storage area** — for each project's target storage formation, using primary and academic web sources.

## When to use
- Use this skill when the task involves filling in reservoir property data for CCS projects that need to be added to CO2BLOCK's screening database.
- Use this skill when you have a `CANADA-CCS-*-missing-revised.md` file (or equivalent) listing projects by name, status, target formation, and regulatory citation.
- Do NOT use this skill for basin-scale properties (use the basin classification / Evenick matching skill instead).

## Workflow

### Phase 1: Identify the Target
1. Open the `*-missing-revised.md` file and extract the list of CCS projects.
2. For each project, note: **Project name**, **Status** (Operational/Approved/Construction), **Storage formation**, and **Regulatory citation**.
3. Select one project at a time. Start with the most mature (Operational) project — it has the most publicly documented data.
4. Confirm the storage formation name, depth range, and geological context before searching for properties.

### Phase 2: Search for Reservoir Properties
For EACH property (thickness, porosity, permeability, area), run targeted web searches using these patterns:

**Thickness:**
```
"<Project Name>" "<Formation Name>" thickness
"<Formation Name>" "average thickness" "<Project>" CCS
"<Formation Name>" "gross thickness" OR "net pay"
```

**Porosity:**
```
"<Project Name>" "<Formation Name>" porosity
"<Project Name>" "average porosity" storage reservoir
```

**Permeability:**
```
"<Project Name>" "<Formation Name>" permeability
"<Formation Name>" "darcy" OR "millidarcy" "<Project>"
"<Project>" injectivity permeability
```

**Area:**
```
"<Project Name>" sequestration lease area OR SLA
"<Project Name>" "leased acres" OR "lease boundary"
ERCB OR AER "<Project>" townships lease approval
```

### Phase 3: Source Prioritization
Rank sources in this order (1 = best):

| Tier | Source Type | Examples |
|------|-------------|----------|
| **1** | Government geological survey reports | AGS PRS, USGS bulletins, provincial surveys |
| **2** | Regulatory agency filings | AER Closure Plans, D-65 Applications, ERCB Decisions, MMV Plans |
| **3** | Peer-reviewed journal articles | Energy Procedia / GHGT proceedings, IJGGG, ScienceDirect |
| **4** | Industry conference presentations | AAPG Search & Discovery, GHGT presentations, SPE papers |
| **5** | Project operator websites / fact sheets | Shell, Chevron, Enhance Energy project pages |
| **6** | Secondary compilations | carbonstorage.io, CSLF, GCCSI, MIT sequestration pages |
| **7** | Encyclopedia / Wiki | Wikipedia (triangulation only — never sole source) |

**Key principle**: Every property value MUST be traceable to at least one Tier 1–3 source. Tier 6–7 may be used for cross-checking and discovery but never as the primary citation.

### Phase 4: Verification & Cross-Validation
1. For each property, find **at least 2 independent sources** that agree within a reasonable range. If they disagree:
   - Record ALL values with their sources
   - Note which source tier each belongs to
   - The higher-tier source should be treated as more authoritative
   - If same-tier sources disagree, flag the discrepancy explicitly
2. Capture **verbatim quotes** from each source — do not paraphrase property values.
3. Record the **exact URL** and retrieval date for every source.

### Phase 5: Resolve Area Discrepancies
Reservoir/storage area is the most commonly discrepant property because:
- **Sequestration Lease Area (SLA)** = the full legal lease boundary (often large, e.g., thousands of km²)
- **Area of Review (AOR)** = the monitoring radius around injection wells (often 10 km radius per well)
- **Pore space tenure / leased acres** = the subsurface pore space rights (may differ from surface SLA)
- **Storage complex area** = the geological extent of the formation

When different numbers appear:
1. Look for the **regulatory source** (AER, ERCB, state oil & gas board) first — this defines the legal boundary.
2. Government documents often describe the SLA in terms of **townships and sections** (Alberta) or **lease blocks** (other jurisdictions).
3. Calculate the area yourself from the township/section description if needed (1 township = 36 sections = ~93.24 km²; 1 section = 640 acres = 2.59 km²).
4. Clearly distinguish between "SLA," "AOR," "plume area," and "pore space acreage" in the output — these are NOT interchangeable.

### Phase 6: Compile Output File
Write a `PROJECT_NAME-props.md` file with this structure:

```markdown
# <Project Name> — Reservoir Properties (<Formation Name>)

**Last updated**: <date>  
**Project**: <full description>  
**Reservoir**: <formation details>  
**Location**: <geographic/legal coordinates>

---

## Sources Summary

| # | Source | Type | Key Properties |
|---|--------|------|----------------|
| 1 | ... | ... | ... |

---

## Reservoir Thickness

### Primary Values
| Value | Source | Verbatim Quote |
|-------|--------|----------------|
| ... | ... | *"..."* |

### Confirming Evidence
- Cross-validation notes
- **Recommended value**: ...

---

## Porosity

### Primary Values
...

---

## Permeability

### Primary Values
...

---

## Area

### Primary Values
...

---

## Other Key Reservoir Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| Depth, pressure, temperature, etc. | ... | ... |

---

## Notes on Data Quality
1. Source reliability assessment
2. Discrepancy explanations
3. Applicability caveats
```

## Key Source Repositories (Canada-Alberta Focus)

### Government / Regulatory (Tier 1–2)
| Repository | URL Pattern | Contents |
|------------|-------------|----------|
| **Open Alberta** | `https://open.alberta.ca/` | Government publications; search by project name + "report" |
| **AER (Alberta Energy Regulator)** | `https://static.aer.ca/prd/documents/by-topic/ccus/` | Closure Plans, MMV Plans, D-65 filings |
| **AGS (Alberta Geological Survey)** | `https://ags.aer.ca/publications/` | PRS bulletins, basin characterization, CCS case studies |
| **Natural Resources Canada** | `https://natural-resources.canada.ca/` | CCS project summaries and funding records |

### Peer-Reviewed (Tier 3)
| Repository | Key Search Strategy |
|------------|---------------------|
| **ScienceDirect** | Search for "Energy Procedia GHGT-13" + project name + formation name |
| **ResearchGate** | Search for author name + "Quest" or project name + publication title |
| **OSTI.gov** | US Department of Energy — hosts GHGT proceedings PDFs as open access |

### Industry (Tier 4–5)
| Repository | Notes |
|------------|-------|
| **AAPG Search & Discovery** | Conference presentations with formation properties |
| **CSEG RECORDER** | Canadian geophysical articles |
| **ESG Solutions** | Industry web pages often reproduce project data from operator presentations |

## Known Source-Specific Caveats

### carbonstorage.io
- Good for **discovery** of projects and general parameters
- Often has unexplained or conflicting values (e.g., "Leased acres" may not match the full SLA)
- **Always verify** every value from their data against a primary regulatory or peer-reviewed source

### AER Directives / D-65 Applications
- Contain **actual well-test data** (not just model averages) — these are point measurements
- Well-scale measurements (core plugs, DSTs) may differ significantly from formation-scale averages
- Document both values with appropriate context

### Energy Procedia / GHGT Papers
- Many papers are **open access** under Creative Commons (CC BY-NC-ND 4.0)
- Shell-authored GHGT papers (Rock, O'Brien, Bourne, Harvey) are the most reliable for Quest properties
- Property values in these papers are typically **formation-scale averages**, not point measurements

### Gen-4 Integrated Reservoir Modeling Reports (Shell)
- Filed with Open Alberta as regulatory submissions
- Extremely detailed (200+ pages) with property tables, grid statistics, and model assumptions
- Best single source for understanding **vertical heterogeneity** and **range of values** across the model area
- Look for "Table 40" style entries in appendices

## Output File Naming Convention
`<TERRITORY>-CCS-<BASIN>-<PROJECT>-props.md`

Examples:
- `CANADA-CCS-ALBERTA-Quest-props.md`
- `CANADA-CCS-ALBERTA-ACTL-Clive-props.md`
- `CANADA-CCS-WILLISTON-BoundaryDam-props.md`

## Specific Search Queries That Worked (Quest CCS Case Study)

Save these as reusable templates. Replace `<Project>` and `<Formation>` with the target names.

```text
# For government documents
open.alberta.ca <Project> report
static.aer.ca ccus <Project>
ags.aer.ca "<Formation>" CCS

# For peer-reviewed papers
"<Formation>" "<Project>" "average" porosity OR permeability
"<Project>" "1st Year" OR "first year" review injection
"<Project>" GHGT-13 OR GHGT-15 OR GHGT-14
"<Project>" Energy Procedia

# For area (Alberta-specific)
"<Project>" townships lease sequestration AER
"ERCB" "<Project>" Decision townships
"Carbon Sequestration Tenure" "<Project>"

# For area (general)
"<Project>" sequestration lease area" acres OR km2
"<Project>" "lease boundary" OR "SLA" storage
```

## Verification Checklist

Before marking a property file as complete, confirm:

- [ ] Each of the 4 properties (thickness, porosity, permeability, area) has at least one Tier 1–3 source
- [ ] At least 2 independent sources exist for each property (not different pages of the same document)
- [ ] Verbatim quotes are captured for each primary value
- [ ] Each value includes the source URL and date accessed
- [ ] Any discrepancies between sources are explicitly flagged and explained
- [ ] The area value is clearly described as SLA / AOR / pore space / other
- [ ] The file includes other key context (depth, pressure, temperature, etc.)
- [ ] Data quality notes explain source hierarchies and any reliability concerns

## Lessons from Quest CCS Research (July 2026)

1. **The GHGT conference proceedings (Energy Procedia) are the single best peer-reviewed source** for CCS project properties across Canada. Multiple Shell-authored papers (GHGT-13 2016, GHGT-15 2021) contain formation-scale average values for the Basal Cambrian Sands.

2. **The AER Closure Plan is the definitive source for storage area** — it defines the SLA in terms of townships and sections. This is more authoritative than any secondary compilation.

3. **Well-test data (D-65 Applications) consistently shows lower permeability than formation-scale averages** — this reflects the difference between a core-plug measurement (~10⁻² m scale) and a formation averaged value (~10³–10⁴ m scale). Both are valid for different purposes.

4. **AGS publications** (PRS 2024-001) give the **regional context** but may use broader ranges than site-specific project data. Use them to validate that project-specific values fall within the regional envelope.

5. **Never take "Leased acres" from carbonstorage.io at face value** — verify against the actual regulatory lease document. The Quest SLA from AER is ~3,667 km² (906,240 acres) but carbonstorage.io reports 164,639 acres. The discrepancy is unresolved but the regulatory source takes precedence.

6. **Shell's Gen-4 Integrated Reservoir Modeling Report** (2011, Open Alberta) is a 247-page document with comprehensive property tables. It gives the full range (40–48 m thickness, 14% mid-case porosity) rather than just the high-quality interval averages (40 m, 17%). Both are useful — document both.
