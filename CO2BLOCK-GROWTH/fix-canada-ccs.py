#!/usr/bin/env python3
"""
Fix CANADA-CCS.csv per adversarial review.

Issues addressed:
  1. Remove all _source_url columns (8 columns)
  2. Wabamun pressure: add AGS PRS 2024-001 source URL
  3. Wabamun area: add carbonstorage.io source URL
  4. Shell Polaris: replace generic open.alberta.ca/publications/ with specific URLs
  5. ACTL/Clive: replace geoconvention.com/2025/ with specific Open Alberta document
  6. Wolf Lamont pressure: add AGS PRS 2024-001 source URL
  7. NWR: strip parenthetical from basin_classification
  8. BD3: strip parenthetical from basin_classification
  9. Quest depth: add attribution clarifying 1,800 m source
"""

import csv
import sys
import io
import re

INPUT = "CANADA-CCS.csv"
OUTPUT = "CANADA-CCS-revised.csv"

# ── Row identifiers (must match CSV exactly) ──
QUEST_SITE = "Quest CCS (Shell)"
ACTL_SITE = "ACTL / Clive CCS (Enhance Energy)"
NWR_SITE = "NWR Sturgeon Refinery CO\u2082 Capture"
POLARIS_SITE = "Shell Polaris + Atlas Hub CCS"
WOLF_SITE = "Wolf Lamont Carbon Hub"
WABAMUN_SITE = "Wabamun Carbon Hub (Enbridge)"
BD3_SITE = "Boundary Dam Integrated CCS (BD3)"

# ── Expected before/after values ──
# (Format: (site, col, old_value, new_value_or_transform)
# Using None for transform means "replace old with new_value"
FIX_TABLE = {
    # Issue 9: Quest depth attribution
    (QUEST_SITE, "depth_m"): (
        "~2,000 m (range 1,800-2,132 m MD) [https://open.alberta.ca/publications/quest-carbon-capture-and-storage-project-d-65-application]",
        "~2,000 m (range 1,800-2,132 m MD; D-65 documents 2,105-2,132 m at injection wells; ~1,800 m at periphery) [https://open.alberta.ca/publications/quest-carbon-capture-and-storage-project-d-65-application]",
    ),
    # Issue 7: NWR basin_classification
    (NWR_SITE, "basin_classification"): (
        "Foreland basin (capture facility \u2014 storage via ACTL pipeline to Clive hub)",
        "Foreland basin",
    ),
    # Issue 8: BD3 basin_classification
    (BD3_SITE, "basin_classification"): (
        "Intracratonic basin (sag basin)",
        "Intracratonic basin",
    ),
    # Issue 2: Wabamun pressure — append URL
    (WABAMUN_SITE, "pressure_mpa"): (
        "~23-33 MPa (hydrostatic estimate at 2,300-3,300 m); may be overpressured",
        "~23-33 MPa (hydrostatic estimate at 2,300-3,300 m); may be overpressured [https://ags.aer.ca/publications/all-publications/prs-2024-001]",
    ),
    # Issue 3: Wabamun area — append URL
    (WABAMUN_SITE, "area_km2"): (
        "~1,500 km\u00b2 (CSA North Area, BSU); full tenure ~2,971 km\u00b2",
        "~1,500 km\u00b2 (CSA North Area, BSU); full tenure ~2,971 km\u00b2 [https://carbonstorage.io/storage/wabamun-hub]",
    ),
    # Issue 6: Wolf Lamont pressure — append URL
    (WOLF_SITE, "pressure_mpa"): (
        ">7.38 MPa (supercritical threshold confirmed); specific formation pressure not published",
        ">7.38 MPa (supercritical threshold confirmed); specific formation pressure not published [https://ags.aer.ca/publications/all-publications/prs-2024-001]",
    ),
}

# URL replacements that should scan columns for a substring and replace
URL_REPLACEMENTS = [
    # Issue 5: ACTL — replace geoconvention URL everywhere it appears
    {
        "site": ACTL_SITE,
        "columns": ["temperature_c", "pressure_mpa", "primary_references"],
        "old_url": "https://geoconvention.com/2025/",
        "new_url": "https://open.alberta.ca/dataset/81932c1f-1e33-4a81-be18-27373cec105e/resource/daae5c3c-5480-458a-96a8-ae93b2f9609b/download/ccsactlreport2015.pdf",
    },
    # Issue 4: Shell Polaris — replace generic URL with specific ones per column
    {
        "site": POLARIS_SITE,
        "columns": ["storage_capacity"],
        "old_url": "https://open.alberta.ca/publications/",
        "new_url": "https://majorprojects.alberta.ca/details/Shell-Polaris-Carbon-Capture-Project/4490",
    },
    {
        "site": POLARIS_SITE,
        "columns": ["area_km2"],
        "old_url": "https://open.alberta.ca/publications/",
        "new_url": "https://chinookpetroleum.com/carbon-sequestration-projects-canada/",
    },
    {
        "site": POLARIS_SITE,
        "columns": ["primary_references"],
        "old_url": "https://open.alberta.ca/publications/",
        "new_url": "https://majorprojects.alberta.ca/details/Atlas-Carbon-Storage-Hub/11237",
    },
]

# ═══════════════════════════════════════════════════════════════
# Read input
# ═══════════════════════════════════════════════════════════════

with open(INPUT, "r", newline="", encoding="utf-8") as f:
    reader = csv.reader(f)
    rows = list(reader)

header = rows[0]
data = rows[1:]

print(f"Read {len(data)} data rows, {len(header)} columns")

# ── Issue 1: Find and drop _source_url columns ──
src_url_indices = [i for i, col in enumerate(header) if col.endswith("_source_url")]
print(f"Issue 1: Dropping {len(src_url_indices)} _source_url columns: {[header[i] for i in src_url_indices]}")

keep_indices = [i for i in range(len(header)) if i not in src_url_indices]
print(f"  Remaining columns: {len(keep_indices)}")

# Helper functions
def old_idx(col_name):
    return header.index(col_name)

def get_field(row, col_name):
    return row[old_idx(col_name)]

def set_field(row, col_name, new_value):
    row[old_idx(col_name)] = new_value

def drop_cols(row):
    return [row[i] for i in keep_indices]

# ═══════════════════════════════════════════════════════════════
# Apply changes
# ═══════════════════════════════════════════════════════════════

changes_log = []
any_fix_mismatch = False

for idx, row in enumerate(data):
    site = get_field(row, "site_name")
    row_num = idx + 2  # 1-indexed + header

    # ── Apply exact-match transforms from FIX_TABLE ──
    for (fix_site, fix_col), (old_val, new_val) in FIX_TABLE.items():
        if site == fix_site:
            current = get_field(row, fix_col)
            if current == old_val:
                set_field(row, fix_col, new_val)
                changes_log.append(f"Row {row_num} ({site}): {fix_col} updated")
            else:
                print(f"  WARN [{fix_site} / {fix_col}]: expected value mismatch — "
                      f"expected={repr(old_val[:60])}..., got={repr(current[:60])}...")
                any_fix_mismatch = True

    # ── Apply URL substring replacements ──
    # Uses regex anchoring to avoid mid-URL prefix concatenation:
    # only replaces when old_url is followed by whitespace, ';', ']', or end-of-string
    for rule in URL_REPLACEMENTS:
        if site == rule["site"]:
            for col_name in rule["columns"]:
                current = get_field(row, col_name)
                pattern = re.escape(rule["old_url"]) + r"(?=\s|;|\]|$)"
                if re.search(pattern, current):
                    new_val = re.sub(pattern, rule["new_url"], current)
                    set_field(row, col_name, new_val)
                    changes_log.append(f"Row {row_num} ({site}): {col_name} URL replaced")
                else:
                    print(f"  Note [{site} / {col_name}]: '{rule['old_url']}' not found — skipping")

# If any exact match failed, abort early
if any_fix_mismatch:
    print("\n❌ One or more expected values did not match input. Aborting to avoid silent corruption.")
    sys.exit(1)

# ═══════════════════════════════════════════════════════════════
# Write output
# ═══════════════════════════════════════════════════════════════

new_header = drop_cols(header)
new_rows = [drop_cols(row) for row in data]

print(f"\nChanges made ({len(changes_log)}):")
for c in changes_log:
    print(f"  ✓ {c}")

print(f"\nOutput: {len(new_rows)} data rows, {len(new_header)} columns")

with open(OUTPUT, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(new_header)
    writer.writerows(new_rows)

print(f"Written to: {OUTPUT}")

# ═══════════════════════════════════════════════════════════════
# Validation
# ═══════════════════════════════════════════════════════════════

errors = []

# Basic structure
if len(new_rows) != 12:
    errors.append(f"Expected 12 data rows, got {len(new_rows)}")
if len(new_header) != 18:
    errors.append(f"Expected 18 columns, got {len(new_header)}")
for col in new_header:
    if col.endswith("_source_url"):
        errors.append(f"Leftover _source_url column: {col}")

# Round-trip integrity
with open(OUTPUT, "r", newline="", encoding="utf-8") as f:
    try:
        reader = csv.reader(f)
        verify_rows = list(reader)
        if len(verify_rows) != 13:
            errors.append(f"Round-trip: expected 13 rows, got {len(verify_rows)}")
        if len(verify_rows[0]) != 18:
            errors.append(f"Round-trip: expected 18 columns, got {len(verify_rows[0])}")
        print(f"Round-trip OK: {len(verify_rows)} rows × {len(verify_rows[0])} cols")
    except Exception as e:
        errors.append(f"CSV parse error on round-trip: {e}")

# Field-level post-condition checks
with open(OUTPUT, "r", newline="", encoding="utf-8") as f:
    reader = csv.reader(f)
    verify_header = next(reader)
    verify_data = list(reader)

def find_row(site_pattern):
    """Find first data row whose site_name contains pattern."""
    si = verify_header.index("site_name")
    for r in verify_data:
        if site_pattern in r[si]:
            return r
    return None

def col_index(name):
    return verify_header.index(name)

# Issue 5: ACTL — all three columns should have replacement URL
actl_row = find_row("ACTL")
if actl_row:
    replacement_present = True
    for col in ["temperature_c", "pressure_mpa"]:
        val = actl_row[col_index(col)]
        if "open.alberta.ca/dataset/81932c1f" not in val:
            errors.append(f"Issue 5 FAIL: ACTL {col} does not contain replacement URL")
            replacement_present = False
    refs = actl_row[col_index("primary_references")]
    if "open.alberta.ca/dataset/81932c1f" not in refs:
        errors.append(f"Issue 5 FAIL: ACTL primary_references does not contain replacement URL")
        replacement_present = False
    if replacement_present:
        print("Issue 5 OK: ACTL temperature/pressure/refs all have replacement URL")
    # Also verify no geoconvention URL remains
    for col in ["temperature_c", "pressure_mpa", "primary_references"]:
        if "geoconvention.com/2025/" in actl_row[col_index(col)]:
            errors.append(f"Issue 5 FAIL: geoconvention.com/2025/ still present in ACTL {col}")
else:
    errors.append("Issue 5 FAIL: Could not find ACTL row for field-level check")

# Issue 4: Polaris — field-level check on replacement URLs
polaris_row = find_row("Polaris")
if polaris_row:
    storage_val = polaris_row[col_index("storage_capacity")]
    area_val = polaris_row[col_index("area_km2")]
    refs_val = polaris_row[col_index("primary_references")]

    if "majorprojects.alberta.ca/details/Shell-Polaris" in storage_val:
        print("Issue 4 OK: Polaris storage_capacity has Major Projects URL")
    else:
        errors.append("Issue 4 FAIL: Polaris storage_capacity missing expected URL")

    if "chinookpetroleum.com/carbon-sequestration-projects-canada" in area_val:
        print("Issue 4 OK: Polaris area_km2 has Chinook URL")
    else:
        errors.append("Issue 4 FAIL: Polaris area_km2 missing expected URL")

    if "majorprojects.alberta.ca/details/Atlas-Carbon-Storage-Hub" in refs_val:
        print("Issue 4 OK: Polaris primary_references has Atlas Hub URL")
    else:
        errors.append("Issue 4 FAIL: Polaris primary_references missing expected URL")

    # Verify no bare open.alberta.ca/publications/ remains in Polaris row
    for col in ["storage_capacity", "area_km2", "primary_references"]:
        val = polaris_row[col_index(col)]
        # Match bare URL (without specific document path after /publications/)
        if re.search(r"https://open\.alberta\.ca/publications/(?:\s|;|\]|$)", val):
            errors.append(f"Issue 4 FAIL: bare open.alberta.ca/publications/ still in Polaris {col}")
else:
    errors.append("Issue 4 FAIL: Could not find Polaris row for field-level check")

# Issue 2: Wabamun pressure — row-level check
wabamun_row = find_row("Wabamun")
if wabamun_row:
    p_val = wabamun_row[col_index("pressure_mpa")]
    if "may be overpressured [https://ags.aer.ca/publications/all-publications/prs-2024-001]" in p_val:
        print("Issue 2 OK: Wabamun pressure_mpa has PRS 2024-001 URL")
    else:
        errors.append("Issue 2 FAIL: Wabamun pressure_mpa missing PRS 2024-001 URL")

    # Issue 3: Wabamun area — row-level check
    a_val = wabamun_row[col_index("area_km2")]
    if "carbonstorage.io/storage/wabamun-hub" in a_val:
        print("Issue 3 OK: Wabamun area_km2 has carbonstorage.io URL")
    else:
        errors.append("Issue 3 FAIL: Wabamun area_km2 missing carbonstorage.io URL")
else:
    errors.append("Issue 2/3 FAIL: Could not find Wabamun row")

# Issue 6: Wolf Lamont pressure — row-level check
wolf_row = find_row("Wolf")
if wolf_row:
    p_val = wolf_row[col_index("pressure_mpa")]
    if "prs-2024-001" in p_val and "supercritical" in p_val:
        print("Issue 6 OK: Wolf Lamont pressure_mpa has PRS 2024-001 URL")
    else:
        errors.append("Issue 6 FAIL: Wolf Lamont pressure_mpa missing PRS 2024-001 URL")
else:
    errors.append("Issue 6 FAIL: Could not find Wolf Lamont row")

# Issue 9: Quest depth — row-level check
quest_row = find_row("Quest CCS")
if quest_row:
    d_val = quest_row[col_index("depth_m")]
    if "D-65 documents 2,105-2,132 m at injection wells" in d_val:
        print("Issue 9 OK: Quest depth_m has D-65 attribution")
    else:
        errors.append("Issue 9 FAIL: Quest depth_m missing D-65 attribution")
else:
    errors.append("Issue 9 FAIL: Could not find Quest row")

# Issue 7: NWR basin_classification — row-level check
nwr_row = find_row("NWR")
if nwr_row:
    bc_val = nwr_row[col_index("basin_classification")]
    if bc_val == "Foreland basin":
        print("Issue 7 OK: NWR basin_classification = 'Foreland basin'")
    else:
        errors.append(f"Issue 7 FAIL: NWR basin_classification = '{bc_val}' (expected 'Foreland basin')")
else:
    errors.append("Issue 7 FAIL: Could not find NWR row")

# Issue 8: BD3 basin_classification — row-level check
bd3_row = find_row("BD3")
if bd3_row:
    bc_val = bd3_row[col_index("basin_classification")]
    if bc_val == "Intracratonic basin":
        print("Issue 8 OK: BD3 basin_classification = 'Intracratonic basin'")
    else:
        errors.append(f"Issue 8 FAIL: BD3 basin_classification = '{bc_val}' (expected 'Intracratonic basin')")
else:
    errors.append("Issue 8 FAIL: Could not find BD3 row")

# ═══════════════════════════════════════════════════════════════
# Final verdict
# ═══════════════════════════════════════════════════════════════

if errors:
    print("\n❌ VALIDATION ERRORS:")
    for e in errors:
        print(f"  ✗ {e}")
    sys.exit(1)
else:
    print("\n✅ All validations passed!")

# Print summary count
print(f"\nTotal changes: {len(changes_log)}")