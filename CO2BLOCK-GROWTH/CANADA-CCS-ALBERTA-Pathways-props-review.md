# Adversarial Review: Pathways Alliance CCS Hub — Reservoir & Project Properties

## Summary

The file is a competent compilation of publicly available data with consistent caveats about the project's pre-FID status. However, the adversarial review identifies **three critical issues** (conflicting capacity targets that are not reconciled, missing cost escalation, and a seal description that conflates stratigraphic intervals), **four significant gaps** (no storage efficiency analysis, no source-tier classification, no cross-reference of Quest's injectivity issues as a risk analog, and no treatment of the consortium governance risk), and **several source-quality issues** that weaken the database's defensibility.

---

## 1. Project Hub / Collective Nature — Correct but Under-Analyzed

| Claim | Assessment |
|-------|------------|
| "Operators: Canadian Natural (lead), Cenovus Energy, ConocoPhillips Canada, Imperial Oil, Suncor Energy" | **Correct.** Matches Oil Sands Alliance Project page and AER filings. |
| "MEG Energy was founding member; current proponent list on AER filings = 5" | **Correct.** Good attention to discrepancy. MEG appears on earlier docs (e.g., Vision Deck Nov 2023) but not on current AER filings. |
| Consortium described as 95% of oilsands production | **Not in this file** (it is in ALBERTA-revised.md). Minor omission. |

**Missing — Governance Risk:** The collective nature has direct implications for FID timing and project risk that should be noted:
- Five (originally six) separate companies must agree on cost-sharing, liability allocation, and offtake terms. This is structurally harder than a single-operator project.
- Canadian Natural is the *applicant of record* for AER/IAAC filings, but the other members are not legally committed until individual FIDs. The consortium could fragment.
- The July 2026 tripartite MOU is described as a milestone, but it is a *non-binding political agreement* (MOU = Memorandum of Understanding), not a binding commercial contract or FID. The file should clarify this distinction.

**Recommendation:** Add a governance risk note similar to the ACTL/Clive props file: consortium FID risk, liability allocation not yet public, and the distinction between the MOU (political) and a binding commercial FID.

---

## 2. Reservoir / Storage References

### 2.1 BCS Depth at Cold Lake — Reasonable Range But Could Be Tighter

| Claim | Assessment |
|-------|------------|
| BCS depth at Cold Lake: ~1,000–2,000 m | **Broadly correct but imprecise.** The AGS PRS 2024-001 regional range is 1,000–3,600 m. Cold Lake is near the eastern subcrop edge of the BCS, where the formation pinches out against the Precambrian Shield. Depth at the storage hub (e.g., Twp 54–56, Rge 1–4 W4M) is likely **1,000–1,500 m**, not 1,000–2,000 m. The BCS at 2,000 m would be significantly farther west (toward Edmonton). The 2,000 m upper bound may be a Quest-analogue transfer. |
| BCS depth range (Alberta-wide): 1,000–3,600 m | **Correct.** AGS PRS 2024-001. |
| Regional gradient applies; site-specific not public | **Correct.** Honest caveat. |

**Recommendation:** Tighten to **~1,000–1,500 m** for the Cold Lake hub location based on the BCS subcrop map (AGS PRS 2024-001, Figure 1). Note that the 2,000 m bound is the formation's western extent, not the Cold Lake area.

### 2.2 Quest Analogue — Useful but Overextended

| Claim | Assessment |
|-------|------------|
| "Same BCS formation, >2,000 m depth, >9 Mt stored since 2015" | **Correct factually.** Quest is in the BCS. |
| Quest injectivity/permeability cited as context | **Potentially misleading without lateral-variability caveats.** Quest is ~200 km west of Cold Lake. The file itself documents extreme BCS heterogeneity (2–3% clay → 3 orders of magnitude permeability reduction). Quest's ~1,000 mD is from the high-quality western BCS fairway. Cold Lake is near the eastern pinch-out where BCS properties are expected to be significantly poorer (thinner, shallower, more cemented). The file uses Quest as an operational analogue in Section 6 but should add a **lateral variability caveat** to every Quest-derived value. |

**Recommendation:** Add a note: "Quest is ~200 km west in a different BCS depositional fairway. Properties (porosity ~17%, permeability ~1,000 mD) are not representative of the Cold Lake area and should not be used as a direct analogue."

### 2.3 Seal Description — Stratigraphically Imprecise

| Claim | Assessment |
|-------|------------|
| "Multiple overlying layers of salt formations that act as seals" | **Partially misleading.** The BCS is Cambrian (~500 Ma). The Prairie Evaporite (Middle Devonian, ~390 Ma) is a salt formation — but it is **~100+ Myr younger** and separated from the BCS by the Middle Cambrian shale, the Winnipeg Formation, the Red River Formation, and other Ordovician–Silurian–Devonian carbonates. The "overlying salt formations" do not sit *directly* on the BCS. The immediate seal for the BCS is the **Middle Cambrian shale** (regionally ~30–50 m thick). The Prairie Evaporite and other Devonian salts provide a *secondary/regional top seal* but are not the primary containment. |
| Seal thickness not publicly specified for Cold Lake hub | **Correct.** Honest caveat. |

**Recommendation:** Restructure the seal description to distinguish:
- **Primary seal:** Middle Cambrian shale (direct top seal for BCS)
- **Secondary/regional seal:** Devonian evaporite sequence (Prairie Evaporite, Lotsberg salts)

### 2.4 Secondary Targets — Mannville Group Depth Issue

| Claim | Assessment |
|-------|------------|
| Lower Grand Rapids Fm net pay >10 m at depth ~375–500 m | **Correct** per CSEG Recorder 2015. |
| "Too shallow for primary CCS" | **Correct.** CO₂ is supercritical only below ~800 m. At 375–500 m, CO₂ would be in gas phase, drastically reducing storage efficiency. |

**Recommendation:** Good. No change needed.

---

## 3. Source Tier Classification (Per CO2BLOCK Matching Skill)

The file does not classify its sources by tier. This should be added for database defensibility.

| # | Source | CO2BLOCK Tier | Confidence | Notes |
|---|--------|---------------|------------|-------|
| 1 | Alberta Major Projects Database | **Tier 2** ★★★★★ | High | Government program record — primary for project status/cost |
| 2 | IAAC Registry #89090 | **Tier 2** ★★★★★ | High | Federal regulatory registry — primary for IA timeline |
| 3 | AGS PRS 2024-001 (Herbers et al.) | **Tier 3** ★★★★☆ | High | Geological Survey — authoritative for BCS regional properties |
| 4 | AGS PRS 2023-001 (Herbers et al.) | **Tier 3** ★★★★☆ | High | Same, petrography focus |
| 5 | CDL Regional CCS Assessment (Mar 2023) | **Tier 3** ★★★★☆ | Medium | Geological Survey-quality methodology. **Proprietary; confidential until Mar 2028.** Cannot be independently verified within CO2BLOCK's audit window. |
| 6 | Pathways Alliance Vision Deck (MEG, Nov 2023) | **Tier 5** ★★★☆☆ | Medium | Industry announcement — MEG's investor deck (not Pathways' own). MEG was a founding member but appears to have left. |
| 7 | Pathways Project Overview (Aug 2023) | **Tier 5** ★★★☆☆ | Medium | Industry announcement — published before Competition Act website takedown. Primary for Phase 1 configuration. |
| 8 | Pathways CCS Fact Sheet (June 2025) | **Tier 5** ★★★☆☆ | Medium | Industry announcement — from oilsandsalliance.ca (rebranded). |
| 9 | Oil Sands Alliance Project page | **Tier 5** ★★★☆☆ | Medium | Industry self-published. |
| 10 | Global News explainer (May 2026) | **Tier 7** ★☆☆☆☆ | Low | Media. Useful for project status/context, not technical data. |
| 11 | Enerdata (Jan 2023) | **Tier 7** ★☆☆☆☆ | Low | Media/aggregator. Source of the 22 Mtpa claim — see Critical Issue #1. |
| 12 | EnergyNow (Dec 2024) | **Tier 7** ★☆☆☆☆ | Low | Media. Useful for test well and regulatory update. |
| 13 | Canadian Indigenous Investment Forum (Dec 2025) | **Tier 5** ★★★☆☆ | Medium | Industry fact sheet. Source of cost breakdown. |
| 14 | CBC News (July 2026) | **Tier 7** ★☆☆☆☆ | Low | Media. Source of MOU announcement and 16 Mtpa target. |
| 15 | The Narwhal (Oct 2024) | **Tier 7** ★☆☆☆☆ | Low | Media/investigative journalism. Useful for St. Paul evaluation agreement. |
| 16 | OurCommons.ca Pathways Brief | **Tier 5** ★★★☆☆ | Medium | Industry submission to Parliamentary Committee. |

**Assessment:** The file relies heavily on Tier 5 (industry) and Tier 7 (media) sources for technical claims. This is not the file's fault — Pathways is pre-FID and has not filed detailed technical data. However, the database should flag all technical parameter values (porosity, permeability, thickness) as having **low confidence** because no Tier 1–2 source has published site-specific Cold Lake data. The file's existing caveats are honest but should be **explicitly mapped to tier categories**.

### 3.1 Critical Source-Specific Issues

- **Source 6 (MEG Vision Deck):** This is hosted on MEG Energy's Q4 investor relations page (s203.q4cdn.com), not on a Pathways Alliance website or a regulatory filing. It is MEG's *own* investor presentation, not a consortium-sanctioned technical document. The file should note this provenance.

- **Source 7 (Project Overview PDF):** The URL (thenarwhal.ca/wp-content/uploads/2026/01/...) means the document was published/released by The Narwhal (a nonprofit investigative news outlet), not by Pathways Alliance directly. The Narwhal likely obtained it through sources. This does not invalidate the content but means it was **not publicly released by the proponent in this form**. The file should note this secondary distribution pathway.

- **Source 5 (CDL):** CDL (Canadian Discovery Ltd.) is a commercial geoscience consulting firm. The study is proprietary and has a **confidentiality period until March 2028**. The file cannot independently verify the porosity mapping, temperature mapping, or pressure mapping attributed to this source. Any database value sourced solely from CDL should be flagged as **PROPRIETARY — UNVERIFIABLE UNTIL 2028**.

- **Source 11 (Enerdata, Jan 2023):** The 22 Mtpa target originates from this single news aggregator article. No other source in the file confirms 22 Mtpa. The 2026 MOU (Source 14) specifies 16 Mtpa by 2045. The Phase 1 target is 10–12 Mtpa. This 22 Mtpa number is an outlier.

---

## 4. Cross-Validation

### 4.1 Capacity Targets — CONFLICTING AND NOT RECONCILED

| Value | Source | Context |
|-------|--------|---------|
| 10–12 Mtpa (Phase 1) | MEG Vision Deck (Nov 2023) | Consistent across most sources |
| 22 Mtpa (long-term target) | Enerdata (Jan 2023) | **Outlier** — appears in no other source |
| 16 Mtpa by 2045 | CBC (May 2026), MOU (July 2026) | **Lower than 22 Mtpa** — suggests target reduction or different scope |
| 40 Mtpa full build-out | CANADA-CCS-ALBERTA-revised.md | **Higher than both** — appears only in the revised master doc |

These numbers do not reconcile. Possible explanations:
- The 22 Mtpa (Enerdata 2023) and 16 Mtpa (MOU 2026) refer to different project phases or scopes
- The target was reduced between 2023 and 2026 (project scope reassessment)
- The 40 Mtpa figure includes third-party volumes (open-access hub expansion)
- The Enerdata figure may be a misinterpretation of the Phase 1 "initial capacity" vs full build-out

**Recommendation:** Add a cross-validation note explaining the discrepancy. The file should list **all three figures** with an explanation: "The 22 Mtpa from Enerdata (Jan 2023) may reflect an initial/pre-MOU target. The 16 Mtpa from the July 2026 MOU is the current political target. The 40 Mtpa referenced in CO2BLOCK's master Alberta assessment appears to be a theoretical full build-out including non-Pathways third-party volumes. These are not directly comparable."

### 4.2 Cost — 2025 vs 2026 Figures Not Reconciled

| Value | Source | Year |
|-------|--------|------|
| Phase 1: $16.5B | Alberta Major Projects | Pre-2026 |
| Additional tech: $7.6B | Canadian Indigenous Investment Forum | Dec 2025 |
| Total Phase 1: $24B | Canadian Indigenous Investment Forum | Dec 2025 |
| Revised to C$20–30B | CarbonCredits.com (per ALBERTA-revised.md) | Jul 2026 |

The file uses the $16.5B + $7.6B = $24B breakdown from the Dec 2025 source but does not mention the **2026 cost escalation** to C$20–30B reported by other sources. Given the file's stated date ambition ("Last updated" implied by the 2026 MOU date), this is a significant omission.

**Recommendation:** Add the cost revision and note: "Original $16.5B Phase 1 estimate (Alberta Major Projects) was reported as C$20–30B by mid-2026 due to inflation and scope changes. The $24B total from the Indigenous Investment Forum (Dec 2025) sits between these estimates."

### 4.3 No Storage Efficiency or Pore Volume Calculation

The file records storage capacity in Mtpa (flow rate) but does not:
- Estimate total pore volume at Cold Lake
- Apply a storage efficiency factor (typically 1–4% for saline aquifers)
- Distinguish between theoretical, contingent, and commercial capacity (SPE SRMS framework)

This is a gap in every CO2BLOCK props file, but it matters most for Pathways because:
- The Phase 1 target (10–12 Mtpa × decades) implies hundreds of Mt total
- Without a pore volume check, there's no way to tell if 10 Mtpa for 30 years (300 Mt) is feasible
- The BCS at Cold Lake is thinner and shallower than at Quest — storage efficiency will be lower

**Recommendation:** Add a first-order calculation:
- BCS thickness at Cold Lake: ~40 m (AGS PRS 2024-001)
- Porosity: ~10% (conservative, based on AGS diagenesis discussion)
- CO₂ density at ~1,200 m, ~45°C: ~650–750 kg/m³
- Storage efficiency: 2–4% (saline aquifer, no hydrodynamic trapping)
- Implied pore volume per km²: ~2.4 × 10⁶ m³/km² × 0.10 = 240,000 m³ pore volume/km²
- Storage per km²: 240,000 m³ × 700 kg/m³ × 0.03 efficiency = ~5,040 tonnes/km²
- At 10 Mtpa × 30 years = 300 Mt: requires ~60,000 km² of BCS — far more than a single hub lease

This back-of-envelope suggests the **10 Mtpa Phase 1 rate requires a very large injection area** or significantly better properties than the conservative estimate above. The file should acknowledge this tension.

### 4.4 No Quest Injectivity Degradation Cross-Reference

The Quest halite damage study (Rock et al. 2022, IJGGG) documents:
- Skin values reaching 147
- Near-wellbore permeability reduced to ~10 mD (from ~1,000 mD native)
- Required freshwater washes to restore injectivity

This is not mentioned in the Pathways file. If Cold Lake BCS has native permeability in the **1–10 mD range** (likely, due to diagenesis and proximity to subcrop), halite-induced damage could reduce injectivity to **near-zero levels**. This is arguably the single most important risk factor for Pathways that could affect the viability of the 10–12 Mtpa Phase 1 target.

**Recommendation:** Add a cross-reference: "Quest has experienced significant halite-induced injectivity degradation despite ~1,000 mD native permeability. For Cold Lake BCS, where native permeability is expected to be 1–3 orders of magnitude lower, halite damage could severely impact injectivity. This risk is heightened by the high-TDS formation water expected in the BCS at the Cold Lake subcrop margin."

### 4.5 Timeline Conflict: "Suspended Dec 2024" vs. "Tripartite MOU July 2026"

| Claim | Source | Assessment |
|-------|--------|------------|
| IAAC timeline suspended Dec 2024 per proponent request | IAAC Registry | **Correct.** |
| Tripartite MOU signed July 13, 2026 | CBC News July 2026 | **Correct but requires context.** |

These are not contradictory — the IAAC federal IA process was paused at Pathways' request in Dec 2024, and the July 2026 MOU is a separate political/funding agreement between the companies and governments. The MOU does not restart the IAAC process. The file should explain this distinction explicitly: the AER pipeline applications (Q1 2024) and the IAAC process (paused) are on separate tracks; the MOU addresses fiscal/policy framework, not regulatory approvals.

### 4.6 Number of Proponents — MEG Status

| Claim | Assessment |
|-------|------------|
| "MEG Energy was founding member (listed on earlier docs; current AER filings show 5 proponents)" | **Correct identification of discrepancy but should resolve it.** The MEG Vision Deck (Nov 2023) lists 6 members. The Oil Sands Alliance page (2025) lists 5. MEG's departure is widely reported. **Recommendation:** State definitively: "MEG Energy appears to have withdrawn from the Pathways Alliance between 2023 and 2025. The Oil Sands Alliance (rebranded Pathways) currently lists 5 proponents (CNRL, Cenovus, ConocoPhillips, Imperial, Suncor). |

### 4.7 Test Wells — Status

| Claim | Assessment |
|-------|------------|
| 2 test wells drilled (2024); plans for 2–3 more | **Source: EnergyNow Dec 2024.** No independent confirmation from AER well licenses or IAAC. This is a media-reported claim from a Tier 7 source. The file should note: "Test well drilling not independently confirmed via AER well license database search." |

---

## 5. Minor Issues

1. **Line 22 (22 Mtpa target):** Cites Enerdata Jan 2023. This is a Tier 7 (media) source for a technical capacity claim. The 22 Mtpa figure appears in no regulatory filing, no Pathways official document, and no other source in the file. Strongly consider removing or adding a "requires verification — single source" caveat.

2. **Line 23 (16 Mtpa by 2045 from CBC May 2026):** This is the MOU target from CBC. The file's reference to "May 2026" (CBC) and "July 13, 2026" (MOU) are slightly confusing — the CBC article is dated "July 2026" in the reference list (Source 14) but says "May 2026" in the inline table (line 23). Verify the publication date.

3. **Line 42 (Regulatory application filed Q1 2024):** This correctly references AER pipeline applications. However, the narrative note (line 8) says "Regulatory Application Filed (Q1 2024)" while line 43 says "Federal IAAC designation Requested Nov 2024." These are **different regulatory processes** (AER provincial vs IAAC federal). The file should explicitly separate these: "AER pipeline agreement applications: filed Q1 2024. IAAC federal impact assessment: requested Nov 2024, suspended Dec 2024."

4. **Line 39 ($16.5B capital cost):** The reference to "Alberta Major Projects database" is correct. However, the API/database field may say "estimated cost" which is a pre-FEED estimate. Add caveat: "Pre-FID estimate; not a sanctioned budget."

5. **Line 152 (Quest analogue):** Says "Quest analogue: Same BCS formation, >2,000 m depth, >9 Mt stored since 2015." Quest's injection depth is ~2,000 m (specifically 1,800–2,100 m). The Cold Lake BCS is at 1,000–2,000 m according to the file's own range. The statement "same formation" is correct but the **depth difference matters** for CO₂ density, injectivity, and storage efficiency. Add: "Note: Quest injects at ~2,000 m (higher CO₂ density, higher injectivity). Cold Lake BCS may be 500–1,000 m shallower."

6. **Line 178 (Total storage resource proprietary):** The CDL study is confidential until March 2028. The file should state the confidentiality expiration date explicitly: "Proprietary results — cannot be independently verified within CO2BLOCK until after Mar 2028."

7. **Missing the Rebranding Narrative:** The file mentions the June 2024 website takedown (Competition Act amendments) and the rebranding to Oil Sands Alliance (line 217). However, it does not connect this to the document itself — some sources are referenced as "Pathways Project Overview" (pre-takedown) and others as "Oil Sands Alliance" (post-takedown). This historical context means **documents published before mid-2024 may not reflect current proponent list or scope**. MEG's departure, for example, happened during this period. Add: "Due to the Competition Act-related website shutdown (June 2024) and subsequent rebranding as Oil Sands Alliance (2025), pre-2024 documents may include MEG Energy as a proponent. Post-2024 documents consistently list 5 proponents."

---

## 6. Recommendations

### Immediate Corrections (Database Impact)

1. **Add source-tier classification** to all table rows (Tier 1–7 per CO2BLOCK matching skill). Flag all technical parameters as low-confidence (pre-FID, no public site-specific data).

2. **Reconcile capacity targets** — explain the discrepancy between 22 Mtpa (Enerdata 2023), 16 Mtpa (MOU 2026), and the 40 Mtpa figure elsewhere. Do not present them as equivalent or additive.

3. **Add cost escalation** — note the revision from $16.5B to $20–30B by mid-2026.

4. **Fix seal description** — distinguish primary seal (Middle Cambrian shale) from secondary/regional seal (Devonian evaporites).

5. **Flag CDL proprietary data** with explicit confidentiality expiration (March 2028).

6. **Add Quest halite injectivity cross-reference** as critical risk context for Cold Lake BCS injectivity.

### Strengthen (Add New Content)

7. **Add first-order storage efficiency calculation** showing implied pore volume vs. the 10–12 Mtpa × decades target. Call out the tension between a large capacity target and limited publicly known pore space.

8. **Add governance risk note** — consortium FID risk, liability allocation, MOU (non-binding) vs. commercial FID distinction.

9. **Add lateral-variability caveat to every Quest-derived value** — Quest is ~200 km west in a different BCS depositional fairway.

10. **Add depth caveat to Quest analogue** — Quest injects at ~2,000 m; Cold Lake BCS may be 500–1,000 m shallower.

### Refine (Improve Existing Content)

11. **Tighten Cold Lake BCS depth range** to ~1,000–1,500 m (eastern subcrop margin) instead of ~1,000–2,000 m.

12. **Add cross-reference** for test wells — flag that EnergyNow (Tier 7) is the sole source; note absence of AER well license confirmation.

13. **Separate AER and IAAC regulatory timelines** explicitly — they are different processes on different tracks.

14. **Add CDL proprietary note** to Section 5 storage capacity table (line 178): "Confidential until March 2028."

---

## Conclusion

The file is well-structured, honestly caveated, and contains accurate data from the best available public sources given Pathways' pre-FID status. The most impactful corrections are:

1. **Capacity targets need reconciliation** — three different numbers (22, 16, 40 Mtpa) from different sources/dates are presented as if equivalent
2. **Cost escalation missing** — the C$20–30B revision (mid-2026) postdates the file's latest cost reference (Dec 2025)
3. **Seal stratigraphy is imprecise** — the evaporite seal is Devonian, not directly overlying the Cambrian BCS; the Middle Cambrian shale is the primary seal
4. **Quest analogue needs stronger caveats** — 200 km lateral separation, 500–1,000 m depth difference, and halite injectivity risk
5. **No source-tier classification** — the file should map all sources to CO2BLOCK's Tier 1–7 system for database defensibility
6. **Governance risk untreated** — the consortium structure, MOU non-binding nature, and MEG departure have implications for project viability that the database should capture
