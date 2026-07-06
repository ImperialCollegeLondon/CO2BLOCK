"""Python port of the MATLAB CO2BLOCK workflows."""

from .allocation import AllocationConfig, run_allocation_workflow
from .core import (
    CalculationConfig,
    CalculationResult,
    RegionalSummaryRow,
    build_regional_summary,
    calculate_site,
    read_site_data,
)
from .screening import assign_screening_case, normalize_country_name, resolve_region_nos, slugify

__all__ = [
    "AllocationConfig",
    "CalculationConfig",
    "CalculationResult",
    "RegionalSummaryRow",
    "build_regional_summary",
    "calculate_site",
    "assign_screening_case",
    "normalize_country_name",
    "read_site_data",
    "resolve_region_nos",
    "run_allocation_workflow",
    "slugify",
]
