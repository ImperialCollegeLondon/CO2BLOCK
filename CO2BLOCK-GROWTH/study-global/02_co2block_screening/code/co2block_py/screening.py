from __future__ import annotations

import re

import pandas as pd


COUNTRY_ALIAS_MAP: dict[str, dict[str, object]] = {
    "US": {"mode": "country", "values": ["USA"]},
    "EU": {
        "mode": "region_excluding_country",
        "values": ["Europe"],
        "exclude": ["UK"],
    },
    "Middle East": {"mode": "region", "values": ["Middle East"]},
}


def slugify(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", value.strip()).strip("_").lower()
    return slug or "unnamed"


def normalize_country_name(name: object) -> str:
    if pd.isna(name):
        return "Unassigned"
    text = str(name).strip()
    if not text:
        return "Unassigned"
    if text == "US":
        return "USA"
    return text


def assign_screening_case(
    majority_country: object, region_group: object, scope: str = "country"
) -> str | None:
    region_text = "" if pd.isna(region_group) else str(region_group).strip()
    if scope == "region":
        return region_text or None

    normalized_country = normalize_country_name(majority_country)
    if normalized_country in {
        "Australia",
        "Brazil",
        "Canada",
        "China",
        "Indonesia",
        "Thailand",
        "UK",
        "USA",
    }:
        return "US" if normalized_country == "USA" else normalized_country
    if region_text == "Europe" and normalized_country != "UK":
        return "EU"
    if region_text == "Middle East":
        return "Middle East"
    return None


def resolve_region_nos(
    metadata_df: pd.DataFrame, case_name: str, scope: str = "country"
) -> tuple[list[int], str]:
    if scope == "region":
        subset = metadata_df[metadata_df["RegionGroup"] == case_name]
        if not subset.empty:
            return subset["Region_no"].astype(int).tolist(), f"RegionGroup == {case_name}"
        return [], "No basin match found"

    subset = metadata_df[metadata_df["AssignedCountry"] == case_name]
    if not subset.empty:
        return subset["Region_no"].astype(int).tolist(), f"AssignedCountry == {case_name}"

    alias = COUNTRY_ALIAS_MAP.get(case_name)
    if alias:
        if alias["mode"] == "country":
            subset = metadata_df[metadata_df["AssignedCountry"].isin(alias["values"])]
            if not subset.empty:
                return subset["Region_no"].astype(int).tolist(), f"AssignedCountry in {alias['values']}"
        if alias["mode"] == "region_excluding_country":
            subset = metadata_df[
                metadata_df["RegionGroup"].isin(alias["values"])
                & ~metadata_df["AssignedCountry"].isin(alias.get("exclude", []))
            ]
            if not subset.empty:
                return subset["Region_no"].astype(int).tolist(), (
                    f"RegionGroup in {alias['values']} and AssignedCountry not in {alias.get('exclude', [])}"
                )
        if alias["mode"] == "region":
            subset = metadata_df[metadata_df["RegionGroup"].isin(alias["values"])]
            if not subset.empty:
                return subset["Region_no"].astype(int).tolist(), f"RegionGroup in {alias['values']}"

    subset = metadata_df[metadata_df["RegionGroup"] == case_name]
    if not subset.empty:
        return subset["Region_no"].astype(int).tolist(), f"RegionGroup == {case_name}"

    return [], "No basin match found"
