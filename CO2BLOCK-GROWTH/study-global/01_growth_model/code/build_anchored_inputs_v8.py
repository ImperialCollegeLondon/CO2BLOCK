"""Build the four anchored-scenario MCMC input CSVs for v8 (2026-06-01).

The four *technical* scenario CSVs (reference/minimum/maximum/growth10) are
produced by ``build_mcmc_inputs_from_workbook.py``. The four *anchored*
scenarios below add per-country 2050 targets on top of the **reference** base
columns (central storage resource C, growth cap 20% except China 25%):

  - us1gt     : US fixed at p2050 = 1.0 Gt/yr; all others unconstrained (Type A)
  - policy    : published national storage targets
  - ipcc_low  : IPCC demand low/high (model reads p2050_lo)
  - ipcc_high : identical file to ipcc_low (model reads p2050_hi)

Policy and IPCC targets are read straight from the workbook so the CSVs stay in
sync with the source. ``ipcc_low`` and ``ipcc_high`` are written from a single
DataFrame, so they are byte-identical by construction (the v8 notebook asserts
this).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parents[2]
WORKBOOK = BASE_DIR / "01_growth_model" / "input" / "2026_global_raw_data_20260601.xlsx"
INPUT_DIR = BASE_DIR / "01_growth_model" / "output" / "mcmc_inputs"
REFERENCE_CSV = INPUT_DIR / "mcmc_input_reference.csv"

# us1gt: only the US carries a fixed target (README 6.1).
US1GT_TARGETS = {"US": 1.0}
EXPECTED_COUNTRIES = [
    "UK",
    "US",
    "EU",
    "China",
    "Middle East",
    "Australia",
    "Canada",
    "Indonesia",
    "Thailand",
    "Brazil",
]


def load_workbook_targets() -> pd.DataFrame:
    """Return per-country policy / IPCC targets from the IPCCreference sheet."""
    df = pd.read_excel(WORKBOOK, sheet_name="IPCCreference")
    df.columns = [str(c).strip() for c in df.columns]
    df["Country"] = df["Country"].astype(str).str.strip()

    def num(x):
        if isinstance(x, (int, float, np.integer, np.floating)):
            return float(x)
        return np.nan  # "None identified", "No data", blanks -> unconstrained

    out = pd.DataFrame(
        {
            "Country": df["Country"],
            "ipcc_lo": df["IPCC Projected Demands: Low [Gt/yr]"].map(num),
            "ipcc_hi": df["IPCC Projected Demands: High [Gt/yr]"].map(num),
            "policy": df["Published storage target [Gt/yr]"].map(num),
        }
    )
    return out.set_index("Country")


def write_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False)
    print(f"  wrote {path.name}  ({len(df)} rows)")


def main() -> None:
    base = pd.read_csv(REFERENCE_CSV)
    base["Country"] = base["Country"].astype(str).str.strip()
    countries = list(base["Country"])
    tgt = load_workbook_targets()

    # --- us1gt ---
    us1gt = base.copy()
    us1gt["p2050_lo"] = [US1GT_TARGETS.get(c, np.nan) for c in countries]
    us1gt["p2050_hi"] = us1gt["p2050_lo"]
    write_csv(us1gt, INPUT_DIR / "mcmc_input_us1gt.csv")

    # --- policy ---
    policy = base.copy()
    pol = [float(tgt.loc[c, "policy"]) if c in tgt.index else np.nan for c in countries]
    policy["p2050_lo"] = pol
    policy["p2050_hi"] = pol
    write_csv(policy, INPUT_DIR / "mcmc_input_policy.csv")

    # --- ipcc_low / ipcc_high (single source, identical files) ---
    ipcc = base.copy()
    ipcc["p2050_lo"] = [float(tgt.loc[c, "ipcc_lo"]) if c in tgt.index else np.nan for c in countries]
    ipcc["p2050_hi"] = [float(tgt.loc[c, "ipcc_hi"]) if c in tgt.index else np.nan for c in countries]
    write_csv(ipcc, INPUT_DIR / "mcmc_input_ipcc_low.csv")
    write_csv(ipcc, INPUT_DIR / "mcmc_input_ipcc_high.csv")

    # --- assertions ---
    lo = pd.read_csv(INPUT_DIR / "mcmc_input_ipcc_low.csv")
    hi = pd.read_csv(INPUT_DIR / "mcmc_input_ipcc_high.csv")
    pd.testing.assert_frame_equal(lo, hi, check_exact=True)
    assert countries == EXPECTED_COUNTRIES, f"Unexpected country pool: {countries}"
    assert policy.loc[policy["Country"] == "Canada", "p2050_lo"].iloc[0] == 0.06
    print("  [assert] ipcc_low == ipcc_high; pool = 10 countries; Canada policy = 0.06")
    print("Done.")


if __name__ == "__main__":
    main()
