# Review: CCS Projects in Alberta — Primary Source Audit

**Review Date**: 2025-07-23  
**Document**: CANADA-CCS-ALBERTA.md  
**Reviewer**: Automated quality and completeness review

---

## Overall Assessment

**Quality**: HIGH — well-structured, methodologically sound, citation-rich.  
**Completeness**: HIGH for Tier 1–5 primary sources; some gaps in interconnection analysis.  
**Confidence**: HIGH — all projects traced to AER Directive 065 approvals, Alberta Hub Selection records, or IAAC filings.

---

## Strengths

1. **Tiered search methodology** (Tiers 1–5) is systematic and reproducible. Each query is documented with source and date.

2. **Regulatory primacy respected** — every project is anchored to a specific AER Directive 065 approval, ERCB decision, or IAAC filing. This is the gold standard for basin CCS auditing.

3. **Spatial verification** — coordinates check against Evenick UBI 13 polygon for every project. Stratigraphic column cross-check confirms all storage formations belong to the WCSB.

4. **Project status granularity** — uses GCI taxonomy correctly (Operational, Under Construction, In Development, Proposed) with clear differentiation between AER-approved vs. evaluation agreement vs. application stages.

5. **Barriers/Risks section** — honest per-project risk identification (capture rate shortfalls, EOR dependency, injectivity uncertainty, cost overruns, Indigenous consultation). This is a genuine analytical contribution.

6. **Negative findings** — explicitly documents what was searched and not found (no projects outside WCSB, no offshore, no DAC with AB storage). This eliminates false negatives.

7. **Geological storage potential** — Tier 3 basin-scale assessments (WASP, HARP, ASAP, Athabasca Area, NE BC Granite Wash) are well-sourced.

---

## Issues & Recommendations

### 1. Temporal Anomaly

| Issue | Detail |
|-------|--------|
| Wolf Lamont approval date | Document records "AER Approval 13513 (Mar 10, 2026)" but the search date is **2025-07-23**. This is a future date. Either the document was updated after the stated search date, or the entry needs correction. |
| **Recommendation** | Verify the Wolf Lamont approval date against the AER DDS record App #1959125. If it is indeed March 2026, update the search metadata to reflect a later search round. |

### 2. Missing Capacity Citations

| Project | Gap |
|---------|-----|
| **Origins** | "Capacity: Not publicly specified" — but the document infers "commercial-scale." The AER Approval 13463 PDF (via ABlawg.ca) likely contains permitted injection rates or pressure limits. |
| **Wabamun** | "Phase 1: ~4 Mtpa" — no specific AER document cited for this number. The Enbridge PR and Major Projects entry should be checked. |
| **Pathways** | "$16.5B capex" — not directly cited. The document references IAAC #89090 and Pathways regulatory overview, but a page/paragraph should be specified. |
| **Recommendation** | Add specific page/paragraph citations for all capacity and cost figures. Where figures come from news articles, note "as reported by" and cross-reference against regulatory filings. |

### 3. Interconnection & Competition Risks

The document treats each project independently but does not address:

- **Pore space competition**: Quest + Polaris/Atlas + Wolf Lamont all target the Basal Cambrian Sands in the same regional area (Radway to Fort Saskatchewan). There is no discussion of pressure interference between these projects. AER Directive 065 requires a "containment" analysis but cumulative pressure effects across multiple nearby sequestration schemes are not discussed.
- **ACTL trunk line capacity**: ACTL is rated at 14.6 Mtpa pipeline design but currently injects ~2–3 Mtpa. Wolf Lamont plans to connect via ACTL. What happens when multiple projects share the same trunk line and disposal formation?
- **Origins vs. ACTL/Clive**: Origins is a new saline hub in the Leduc Formation near Lacombe — the same formation as ACTL/Clive (Leduc D-3A). Whether Origins taps into a separate pool or competes for the same reservoir pressure support is not addressed.

**Recommendation**: Add an "Interference & Interconnection" subsection to the Barriers/Risks section.

### 4. Induced Seismicity

No mention of induced seismicity risk for any project. Key context:

- The WCSB has documented induced seismicity from hydraulic fracturing in the Duvernay and Montney plays (Fox Creek, etc.).
- CO₂ injection into deep carbonates (Leduc) and basal clastics (BCS) can alter pore pressure and potentially trigger slip on pre-existing faults.
- Basel monitoring is standard for Quest (microseismic array), but the document does not compare monitoring frameworks across projects.

**Recommendation**: Add an induced seismicity risk row to each project's risk table, or a cross-cutting note in the geological storage potential section.

### 5. Regulatory Gaps

The document thoroughly catalogs the existing regulatory framework (Directive 065, MMA, OGCA) but does not identify gaps:

- **Commercial-scale closure framework**: Quest is the only project with a published MMV plan. Directive 065 App P provides closure criteria, but no project has reached the closure stage. The regulatory framework for post-closure liability transfer is untested.
- **Long-term liability**: Canada lacks a federal CCS liability mechanism (unlike the US EPA Class VI primacy or UK CO₂ Storage Licensing). Alberta's MMA transfers liability to the Crown post-closure, but the financial assurance requirements for large hubs (Pathways at $16.5B) are not discussed.
- **Cross-border storage**: Alberta CCS projects are within a single province, but the document could note that the WCSB extends into BC, Saskatchewan, and Montana/North Dakota. No cross-jurisdictional transport/storage is identified (correctly), but this should be explicit.

**Recommendation**: Add a "Regulatory Gaps & Uncertainties" subsection.

### 6. Tier 5–7 Exclusions

The "Tier 5–7 Mentions" section is brief (5 items). For completeness, consider documenting:

- **Alberta Innovates / Emissions Reduction Alberta (ERA) funded studies** that did not result in projects (e.g., Carbon XPRIZE, various feasibility studies).
- **University research pilot injections** (e.g., University of Calgary / CMCR pilot projects) — these are not commercial CCS but are often confused in secondary literature.
- **Methane pyrolysis / blue hydrogen projects** that capture CO₂ for dedicated storage vs. EOR (e.g., Air Products Net-Zero H₂ complex — which is listed as a Wolf Lamont capture source but not as a standalone project).

**Recommendation**: Expand the exclusions list with a brief rationale for each, and note any that appear in secondary sources as "Alberta CCS" but do not meet Tier 1 criteria.

### 7. Minor Edits

| Line | Issue | Suggestion |
|------|-------|------------|
| 51 | ">9 Mt stored as of 2024" | CER Market Snapshot 2025 Table 1 should be explicitly cited for cumulative injection totals. |
| 64 | ">7 Mt stored by YE 2021" | Update to latest available year. The document is dated 2025-07-23 but uses 2021 data. |
| 89 | Atlas Hub coordinates | "Twp 54–56, Rge 16–18 W4M" — appears to overlap with Wolf Lamont coordinates. Confirm these are distinct project areas within the same township range. |
| 152 | "GSC Open File 8996 (Carey & Durling 2023)" | This reference is for **Atlantic Canada** per the citation, but it is used to support WASP capacity estimates. Verify applicability or correct the citation. |
| 231 | "Alberta Major Projects #10970 (ACTL Edmonton Connector)" | Not listed in the findings summary table. Either promote to the main table or note why it is excluded. |

### 8. Missing Metadata

- **Document version** — no version number or revision history. First version is fine, but add a version field for future updates.
- **Reviewer** — no reviewer attribution. Add a reviewer line.

---

## Verdict

**Fit for purpose**. This is a rigorous, Tier-1-anchored audit of Alberta CCS projects within the WCSB. The eight identified projects are correctly classified and verified. The document meets its stated objectives.

**Priority fixes** (before publication or downstream use):
1. Verify the Wolf Lamont approval date anomaly (Mar 2026 vs. Jul 2025 search date).
2. Add induced seismicity and pore-space interference analysis.
3. Address the WASP/GSC Open File 8996 citation mismatch.
4. Add specific capacity/cost figure citations.

**Nice-to-have improvements**:
- Regulatory gaps subsection.
- Expanded Tier 5–7 exclusions.
- Version history table.

---

*End of review*
