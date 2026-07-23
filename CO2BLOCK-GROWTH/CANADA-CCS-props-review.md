# CANADA-CCS-props.md — Source Citation Audit

**Audit date:** 23 July 2026  
**Audit method:** Grep for `https?://` and domain-name patterns across the aggregate file and all 12 constituent `*-props-revised.md` files.  
**Finding:** **100% of citations in the aggregate file are name-only.** Zero hyperlinks/URLs.

---

## Results

| Pattern searched | Matches in aggregate file |
|---|---|
| `https?://` (any URL) | **0** |
| `http` (any occurrence) | **0** |
| `www\.` | **0** |
| `\.com`, `\.org`, `\.gov`, `\.ca` | **2** (both bare domain names, not hyperlinks) |
| Markdown link syntax `[text](url)` | **0** |

The two bare domain references found:
- Line 252: `keyfactsenergy.com` (cited as source text, not a clickable link)
- Line 631: `ABlawg.ca` (cited as source text, not a clickable link)

---

## Contrast with constituent files

All 12 `*-props-revised.md` files contain working hyperlinks in their source-attribution tables (DOI links, direct PDF URLs, government registry links). The aggregate file strips these out entirely, leaving only human-readable citations (author names, publication titles, organization names, report codes).

---

## Examples of name-only citations in the aggregate

| Location | Citation style | Example |
|---|---|---|
| Property tables | Author/publisher + year | `Rock et al. GHGT-13 (2017)` |
| Property tables | Report code | `Shell D-65 Application` |
| Reference Index §5 | Full bibliographic | `Rock et al. (2022) IJGGG — Quest halite damage` |

None of these are clickable or carry a URL.

---

## Verdict

**CRITICAL FINDING**: The aggregate file is a readable summary but is **not independently verifiable** from within the document itself. A reader cannot click through to any source. This is acceptable for a compiled overview that cross-references the detailed constituent files (where URLs live), but the document should make this dependency explicit — ideally with a note at the top stating "All source URLs are in the corresponding `*-props-revised.md` file listed in each section."

**Recommendation**: Consider adding a header notice: *"Source hyperlinks are maintained in the individual project property files referenced in each section. This aggregate document cites sources by name only for readability."*
