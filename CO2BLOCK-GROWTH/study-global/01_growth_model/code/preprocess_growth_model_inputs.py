from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


TAKEOFF_YEAR = 2030
INFLECTION_AFTER_YEAR = 2050

KNOWN_COUNTRY_ALIASES = {
    "australia": "Australia",
    "canada": "Canada",
    "china": "China",
    "indonesia": "Indonesia",
    "thailand": "Thailand",
    "uk": "UK",
    "unitedkingdom": "UK",
    "united_kingdom": "UK",
    "britain": "UK",
    "us": "US",
    "usa": "US",
    "unitedstates": "US",
    "united_states": "US",
    "middleeast": "Middle East",
    "middle_east": "Middle East",
    "eu": "EU",
    "europe": "EU",
    "norway": "EU",
}


@dataclass
class ParsedCountryFile:
    source_file: str
    country: str
    country_inferred: bool
    year_col: str | None
    annual_col: str | None
    cumulative_mt_col: str | None
    cumulative_gt_col: str | None
    data_year_start: float
    data_year_end: float
    p2030_mt: float
    annual_2030_mtpy: float
    cumulative_2030_gt: float
    notes: list[str]
    structure_ambiguous: bool = False

    def partial_row(self) -> dict[str, object]:
        return {
            "country": self.country,
            "source_file": self.source_file,
            "p2030_mt": self.p2030_mt,
            "annual_2030_mtpy": self.annual_2030_mtpy,
            "cumulative_2030_gt": self.cumulative_2030_gt,
            "data_year_start": self.data_year_start,
            "data_year_end": self.data_year_end,
            "notes": "; ".join(self.notes),
        }

    def template_row(self) -> dict[str, object]:
        is_china = self.country == "China"
        notes = list(self.notes)
        if not self.country_inferred:
            notes.append("country requires manual confirmation before applying any China-specific override")
        notes.append("resource_central_gt must be supplied manually")
        return {
            "country": self.country,
            "p2030_mt": self.p2030_mt,
            "resource_central_gt": math.nan,
            "resource_min_gt": math.nan,
            "resource_max_gt": math.nan,
            "growth_cap_reference": 0.25 if is_china else 0.20,
            "growth_cap_minimum": 0.10,
            "growth_cap_maximum": 0.25 if is_china else 0.20,
            "growth_cap_growth10": 0.10,
            "takeoff_year": TAKEOFF_YEAR,
            "inflection_must_be_after_year": INFLECTION_AFTER_YEAR,
            "notes": "; ".join(notes),
        }


def normalise_token(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", text.lower())


def infer_country_from_filename(path: Path) -> tuple[str, bool]:
    stem = path.stem.lower()
    compact = normalise_token(stem)
    for alias, country in KNOWN_COUNTRY_ALIASES.items():
        alias_compact = normalise_token(alias)
        if alias_compact and alias_compact in compact:
            return country, True
    return "", False


def clean_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.astype(str).str.strip().str.replace(",", "", regex=False), errors="coerce")


def score_column(name: str, all_terms: Iterable[str], any_terms: Iterable[str] = ()) -> int:
    label = normalise_token(name)
    if not all(term in label for term in all_terms):
        return -1
    score = 10
    for term in any_terms:
        if term in label:
            score += 1
    if "unnamed" in label:
        score -= 3
    return score


def choose_column(columns: Iterable[str], all_terms: Iterable[str], any_terms: Iterable[str] = ()) -> tuple[str | None, bool]:
    scored = [(col, score_column(col, all_terms, any_terms)) for col in columns]
    scored = [(col, score) for col, score in scored if score >= 0]
    if not scored:
        return None, False
    scored.sort(key=lambda item: item[1], reverse=True)
    ambiguous = len(scored) > 1 and scored[0][1] == scored[1][1]
    return scored[0][0], ambiguous


def choose_year_column(columns: Iterable[str]) -> tuple[str | None, bool]:
    preferred = []
    fallback = []
    for col in columns:
        label = normalise_token(col)
        if label in {"year", "years"}:
            preferred.append(col)
        elif "year" in label and "storage" not in label and "rate" not in label and "cumulative" not in label:
            fallback.append(col)
    if len(preferred) == 1:
        return preferred[0], False
    if len(preferred) > 1:
        return preferred[0], True
    if len(fallback) == 1:
        return fallback[0], False
    if len(fallback) > 1:
        return fallback[0], True
    return None, False


def detect_main_block(df: pd.DataFrame) -> tuple[dict[str, str | None], bool]:
    year_col, amb_year = choose_year_column(df.columns)
    annual_col, amb_annual = choose_column(df.columns, ["storage", "rate"], ["mtyear", "mtyr", "mtyear", "mt", "year"])
    cumulative_mt_col, amb_cum_mt = choose_column(df.columns, ["cumulative", "storage"], ["mt"])
    cumulative_gt_col, amb_cum_gt = choose_column(df.columns, ["cumulative", "storage"], ["gt"])
    return {
        "year": year_col,
        "annual": annual_col,
        "cumulative_mt": cumulative_mt_col,
        "cumulative_gt": cumulative_gt_col,
    }, any([amb_year, amb_annual, amb_cum_mt, amb_cum_gt])


def parse_country_file(path: Path) -> ParsedCountryFile:
    df = pd.read_csv(path)
    mapping, ambiguous = detect_main_block(df)
    notes: list[str] = []
    country, country_inferred = infer_country_from_filename(path)
    if not country_inferred:
        notes.append("country could not be confidently inferred from filename")

    if not mapping["year"] or not mapping["annual"] or not mapping["cumulative_mt"]:
        notes.append("main time-series block not fully detected")
        if not mapping["cumulative_gt"]:
            notes.append("cumulative storage in Gt column not found")
        return ParsedCountryFile(
            source_file=str(path.relative_to(path.parents[2])),
            country=country,
            country_inferred=country_inferred,
            year_col=mapping["year"],
            annual_col=mapping["annual"],
            cumulative_mt_col=mapping["cumulative_mt"],
            cumulative_gt_col=mapping["cumulative_gt"],
            data_year_start=math.nan,
            data_year_end=math.nan,
            p2030_mt=math.nan,
            annual_2030_mtpy=math.nan,
            cumulative_2030_gt=math.nan,
            notes=notes,
            structure_ambiguous=ambiguous,
        )

    work = df.copy()
    work["_year"] = clean_numeric(work[mapping["year"]])
    work["_annual"] = clean_numeric(work[mapping["annual"]])
    work["_cumulative_mt"] = clean_numeric(work[mapping["cumulative_mt"]])
    work["_cumulative_gt"] = clean_numeric(work[mapping["cumulative_gt"]]) if mapping["cumulative_gt"] else math.nan
    work = work[work["_year"].notna()].copy()
    work = work.sort_values("_year")

    if ambiguous:
        notes.append("column detection was ambiguous; please review mapping")
    if mapping["cumulative_gt"] is None:
        notes.append("cumulative storage in Gt column not found")

    year_start = float(work["_year"].min()) if not work.empty else math.nan
    year_end = float(work["_year"].max()) if not work.empty else math.nan

    row_2030 = work[work["_year"] == TAKEOFF_YEAR]
    if row_2030.empty:
        notes.append("Year = 2030 missing from main time series")
        p2030_mt = math.nan
        annual_2030_mtpy = math.nan
        cumulative_2030_gt = math.nan
    else:
        row = row_2030.iloc[0]
        p2030_mt = float(row["_cumulative_mt"]) if pd.notna(row["_cumulative_mt"]) else math.nan
        annual_2030_mtpy = float(row["_annual"]) if pd.notna(row["_annual"]) else math.nan
        cumulative_2030_gt = float(row["_cumulative_gt"]) if pd.notna(row["_cumulative_gt"]) else math.nan
        if pd.isna(p2030_mt):
            notes.append("2030 cumulative storage (Mt) missing")
        if pd.isna(annual_2030_mtpy):
            notes.append("2030 annual storage rate (Mt/yr) missing")
        if pd.isna(cumulative_2030_gt):
            notes.append("2030 cumulative storage (Gt) missing")

    notes.append("resource_central_gt not present in raw time-series block")

    return ParsedCountryFile(
        source_file=str(path.relative_to(path.parents[2])),
        country=country,
        country_inferred=country_inferred,
        year_col=mapping["year"],
        annual_col=mapping["annual"],
        cumulative_mt_col=mapping["cumulative_mt"],
        cumulative_gt_col=mapping["cumulative_gt"],
        data_year_start=year_start,
        data_year_end=year_end,
        p2030_mt=p2030_mt,
        annual_2030_mtpy=annual_2030_mtpy,
        cumulative_2030_gt=cumulative_2030_gt,
        notes=notes,
        structure_ambiguous=ambiguous,
    )


def write_yaml(path: Path) -> None:
    text = """# Scenario configuration for first-stage top-down growth modelling
# Policy scenarios are intentionally excluded at this preprocessing stage.
# resource_upper_bound describes how each scenario should interpret resource_central_gt.

takeoff_year: 2030
inflection_must_be_after_year: 2050
growth_rate_definition: annualised growth rate of the annual storage-rate curve from take-off year to inflection year

scenarios:
  Reference:
    resource_upper_bound: central estimate
    growth_cap_default: 0.20
    growth_cap_overrides:
      China: 0.25
    takeoff_year: 2030
    inflection_year_rule: greater than 2050

  Minimum:
    resource_upper_bound: 0.1 x central estimate
    growth_cap_default: 0.10
    growth_cap_overrides: {}
    takeoff_year: 2030
    inflection_year_rule: greater than 2050

  Maximum:
    resource_upper_bound: 10 x central estimate
    growth_cap_default: 0.20
    growth_cap_overrides:
      China: 0.25
    takeoff_year: 2030
    inflection_year_rule: greater than 2050

  Growth10:
    resource_upper_bound: central estimate
    growth_cap_default: 0.10
    growth_cap_overrides: {}
    takeoff_year: 2030
    inflection_year_rule: greater than 2050
"""
    path.write_text(text)


def write_report(path: Path, parsed: list[ParsedCountryFile]) -> None:
    success = [p for p in parsed if not any("main time-series block not fully detected" in n for n in p.notes)]
    missing_2030 = [p for p in parsed if any("Year = 2030 missing" in n or "2030 " in n for n in p.notes)]
    ambiguous = [p for p in parsed if p.structure_ambiguous]
    manual_resource = list(parsed)

    def fmt(items: list[ParsedCountryFile], extra: str = "") -> str:
        if not items:
            return "- None\n"
        lines = []
        for item in items:
            country_label = item.country if item.country else "[manual country mapping required]"
            lines.append(f"- `{item.source_file}` -> `{country_label}`{extra}")
        return "\n".join(lines) + "\n"

    lines = [
        "# Parsing Report",
        "",
        "## Successfully parsed main country-level time-series block",
        fmt(success).rstrip(),
        "",
        "## Files with missing or incomplete 2030 extraction",
        fmt(missing_2030).rstrip(),
        "",
        "## Files with ambiguous structure",
        fmt(ambiguous).rstrip(),
        "",
        "## Files still needing manual `resource_central_gt` input",
        fmt(manual_resource).rstrip(),
        "",
        "## Notes",
    ]
    for item in parsed:
        note_text = "; ".join(item.notes) if item.notes else "no issues detected"
        country_label = item.country if item.country else "[manual country mapping required]"
        lines.append(f"- `{item.source_file}` -> `{country_label}`: {note_text}")
    path.write_text("\n".join(lines) + "\n")


def write_readme(path: Path) -> None:
    text = """# Growth Preprocessing Outputs

This folder contains preprocessing outputs that convert country-specific raw CSV files into the input structure needed for the first-stage top-down growth modelling.

## What was extracted from raw files

For each raw file, the parser tries to detect the main country-level time-series block using:
- `Year`
- `Storage Rate (Mt/Year)`
- `Cumulative Storage (Mt)`
- `Cumulative Storage (Gt)`

Auxiliary facility tables and `Unnamed:*` columns are ignored unless they interfere with detection.

The parser extracts:
- `p2030_mt`: cumulative storage in Mt at Year = 2030
- `annual_2030_mtpy`: annual storage rate in Mt/yr at Year = 2030
- `cumulative_2030_gt`: cumulative storage in Gt at Year = 2030
- data coverage start and end year
- notes for missing or ambiguous fields

## Files in this folder

- `country_growth_inputs_partial.csv`
  - Direct extraction from raw files.
- `country_growth_inputs_template.csv`
  - Template for first-stage growth modelling with scenario growth caps prefilled.
- `scenario_config.yaml`
  - Scenario definitions for Reference, Minimum, Maximum, and Growth10.
- `parsing_report.md`
  - Parsing status and manual follow-up items.

## What still needs manual supplementation

- `resource_central_gt`
  - This is not expected to be present in the raw time-series files.
  - It must be added later from another resource table or from manual compilation.
- `resource_min_gt` and `resource_max_gt`
  - These remain blank until `resource_central_gt` is available.
  - Once available:
    - `resource_min_gt = 0.1 * resource_central_gt`
    - `resource_max_gt = 10 * resource_central_gt`
- `country`
  - If a filename does not clearly identify the country, manual mapping is still required.

## Scenario configuration structure

`scenario_config.yaml` stores scenario-level rules only:
- resource upper-bound logic
- growth cap defaults
- China-specific override where applicable
- fixed `takeoff_year = 2030`
- fixed `inflection_year > 2050`

It does not contain country-specific resource values.
"""
    path.write_text(text)


def build_outputs(input_dir: Path, output_dir: Path) -> None:
    raw_files = sorted(input_dir.glob("*.csv"))
    parsed = [parse_country_file(path) for path in raw_files]

    partial_df = pd.DataFrame([item.partial_row() for item in parsed])
    template_df = pd.DataFrame([item.template_row() for item in parsed])

    output_dir.mkdir(parents=True, exist_ok=True)
    partial_df.to_csv(output_dir / "country_growth_inputs_partial.csv", index=False)
    template_df.to_csv(output_dir / "country_growth_inputs_template.csv", index=False)
    write_yaml(output_dir / "scenario_config.yaml")
    write_report(output_dir / "parsing_report.md", parsed)
    write_readme(output_dir / "README.md")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Preprocess raw country files for top-down growth modelling.")
    parser.add_argument(
        "--input-dir",
        default="01_growth_model/input/raw",
        help="Directory containing raw country-specific CSV files.",
    )
    parser.add_argument(
        "--output-dir",
        default="01_growth_model/output/preprocessing",
        help="Directory to write preprocessing outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_outputs(Path(args.input_dir), Path(args.output_dir))


if __name__ == "__main__":
    main()
