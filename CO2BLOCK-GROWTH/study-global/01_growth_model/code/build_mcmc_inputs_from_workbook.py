from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

BASE_YEAR = 2030
END_YEAR = 2430
WORKBOOK_DEFAULT = "01_growth_model/input/2026_global_raw_data_20260601.xlsx"

# v8 (2026-06-01) country pool. EU is sourced from the "EU without UK"
# workbook sheet but labelled "EU" in all outputs/figures.
SHEET_TO_COUNTRY = {
    "UK": "UK",
    "US": "US",
    "EU without UK": "EU",
    "China": "China",
    "Middle East": "Middle East",
    "Australia": "Australia",
    "Canada": "Canada",
    "Indonesia": "Indonesia",
    "Thailand": "Thailand",
    "Brazil": "Brazil",
}

COUNTRY_ALIASES = {
    "EU": ["EU"],
}


def parse_number(value: object) -> float:
    if pd.isna(value):
        return float("nan")
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip().replace(",", "")
    match = re.search(r"-?\d+(?:\.\d+)?", text)
    return float(match.group(0)) if match else float("nan")


def load_capacity_table(workbook: Path) -> pd.DataFrame:
    df = pd.read_excel(workbook, sheet_name="Geological Storage Capacity")
    df.columns = [str(c).strip() for c in df.columns]
    df["Country"] = df["Country"].astype(str).str.strip()
    df["Storage capacity (Gt)"] = df["Storage resource estimate (Gt)"].map(parse_number)
    return df[["Country", "Storage capacity (Gt)"]].copy()


def resolve_capacity_country(capacity_index: pd.Index, country: str) -> str:
    if country in capacity_index:
        return country
    for alias in COUNTRY_ALIASES.get(country, []):
        if alias in capacity_index:
            return alias
    raise KeyError(f"Capacity table missing country: {country}")


def extract_2030_from_sheet(workbook: Path, sheet_name: str) -> dict[str, float]:
    raw = pd.read_excel(workbook, sheet_name=sheet_name, header=None)
    summary = pd.DataFrame(
        {
            "Year": pd.to_numeric(raw.iloc[:, 0], errors="coerce"),
            "2030 Storage Rate (Mt/Year)": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
            "2030 Cumulative Storage (Mt)": pd.to_numeric(raw.iloc[:, 2], errors="coerce"),
            "2030 Cumulative Storage (Gt)": pd.to_numeric(raw.iloc[:, 3], errors="coerce"),
        }
    )
    summary = summary[
        summary["Year"].notna()
        & summary["2030 Storage Rate (Mt/Year)"].notna()
        & summary["2030 Cumulative Storage (Mt)"].notna()
    ].copy()
    row_2030 = summary[summary["Year"] == BASE_YEAR]
    if row_2030.empty:
        raise ValueError(f"{sheet_name}: missing summary row for {BASE_YEAR}")

    row = row_2030.iloc[-1]
    cumulative_gt = row["2030 Cumulative Storage (Gt)"]
    if pd.isna(cumulative_gt):
        cumulative_gt = row["2030 Cumulative Storage (Mt)"] / 1000.0

    return {
        "2030 Storage Rate (Mt/Year)": float(row["2030 Storage Rate (Mt/Year)"]),
        "2030 Cumulative Storage (Mt)": float(row["2030 Cumulative Storage (Mt)"]),
        "2030 Cumulative Storage (Gt)": float(cumulative_gt),
    }


def build_country_base_table(workbook: Path) -> pd.DataFrame:
    capacity = load_capacity_table(workbook).set_index("Country")
    rows: list[dict[str, object]] = []

    for sheet_name, country in SHEET_TO_COUNTRY.items():
        metrics = extract_2030_from_sheet(workbook, sheet_name)
        capacity_country = resolve_capacity_country(capacity.index, country)
        central_capacity = float(capacity.loc[capacity_country, "Storage capacity (Gt)"])
        reference_growth_pct = 25.0 if country == "China" else 20.0

        rows.append(
            {
                "Country": country,
                "Storage capacity central (Gt)": central_capacity,
                "2030 Cumulative Storage (Gt)": metrics["2030 Cumulative Storage (Gt)"],
                "2030 Cumulative Storage (Mt)": metrics["2030 Cumulative Storage (Mt)"],
                "2030 Storage Rate (Mt/Year)": metrics["2030 Storage Rate (Mt/Year)"],
                "Reference growth to Tn (%)": reference_growth_pct,
                "BASE_YEAR": BASE_YEAR,
                "END_YEAR": END_YEAR,
                "Source sheet": sheet_name,
            }
        )

    return pd.DataFrame(rows)


def scenario_frame(base: pd.DataFrame, scenario: str) -> pd.DataFrame:
    out = base.copy()
    if scenario == "reference":
        out["Storage capacity (Gt)"] = out["Storage capacity central (Gt)"]
        out["Growth_to_Tn_%"] = out["Reference growth to Tn (%)"]
    elif scenario == "minimum":
        out["Storage capacity (Gt)"] = out["Storage capacity central (Gt)"] * 0.1
        out["Growth_to_Tn_%"] = 10.0
    elif scenario == "maximum":
        out["Storage capacity (Gt)"] = out["Storage capacity central (Gt)"] * 10.0
        out["Growth_to_Tn_%"] = out["Reference growth to Tn (%)"]
    elif scenario == "growth10":
        out["Storage capacity (Gt)"] = out["Storage capacity central (Gt)"]
        out["Growth_to_Tn_%"] = 10.0
    else:
        raise ValueError(f"Unknown scenario: {scenario}")

    out = out[
        [
            "Country",
            "Storage capacity (Gt)",
            "2030 Cumulative Storage (Gt)",
            "2030 Cumulative Storage (Mt)",
            "2030 Storage Rate (Mt/Year)",
            "Growth_to_Tn_%",
            "BASE_YEAR",
            "END_YEAR",
        ]
    ].copy()
    out["Storage capacity (Gt)"] = out["Storage capacity (Gt)"].round(3)
    out["2030 Cumulative Storage (Gt)"] = out["2030 Cumulative Storage (Gt)"].round(5)
    out["2030 Cumulative Storage (Mt)"] = out["2030 Cumulative Storage (Mt)"].round(2)
    out["2030 Storage Rate (Mt/Year)"] = out["2030 Storage Rate (Mt/Year)"].round(2)
    out["Growth_to_Tn_%"] = out["Growth_to_Tn_%"].round(1)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Build four scenario-specific MCMC input CSVs directly from the workbook.")
    parser.add_argument("--workbook", default=WORKBOOK_DEFAULT, help="Path to source workbook.")
    parser.add_argument("--output-dir", default="01_growth_model/output/mcmc_inputs", help="Output directory.")
    args = parser.parse_args()

    workbook = Path(args.workbook)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    base = build_country_base_table(workbook)
    scenario_frame(base, "reference").to_csv(outdir / "mcmc_input_reference.csv", index=False)
    scenario_frame(base, "minimum").to_csv(outdir / "mcmc_input_minimum.csv", index=False)
    scenario_frame(base, "maximum").to_csv(outdir / "mcmc_input_maximum.csv", index=False)
    scenario_frame(base, "growth10").to_csv(outdir / "mcmc_input_growth10.csv", index=False)
    # README is hand-maintained at ../output/mcmc_inputs/README.md; do not overwrite.


if __name__ == "__main__":
    main()
