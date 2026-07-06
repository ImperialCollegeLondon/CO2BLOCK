#!/usr/bin/env python3
"""PNAS main display items — V5 (medium-size set: 5 compact figures + Table 1).

Title: "Carbon storage can scale up for climate mitigation for decades before its
geological limits begin to bind."

PNAS asks for MEDIUM-SIZE graphical elements (column widths 8.7 / 11.4 / 17.8 cm;
modest heights). V5 keeps V4's data semantics and numbers EXACTLY (same sha256-seeded
common-random-number bootstrap, same whole-period-feasibility accounting, same common
country sets) but re-lays-out the display items:

  F1  Demand        17.8 cm x 8.4 cm   (V4: 17 cm tall -> one-row, rain dots dropped)
  F2  Supply        17.8 cm x 7.6 cm   (unchanged content, tightened)
  F3  Test          17.8 cm x 10.9 cm  (centerpiece keeps the most space)
  F4  Outcome       17.8 cm x 8.6 cm   (NEW: exemplar cell + ladder + drivers;
                                        the V4 8-scenario matrix moves to SI = V4 fig4)
  F5  Where it binds 11.4 cm x 9.9 cm  (NEW 1.5-column "commitment collapse" scatter;
                                        the V4 10-country grid moves to SI = V4 fig5)
  T1  Scenario outcomes table           (tables/table1_scenario_outcomes.{csv,md})

Fig-4 semantics (identical to V4 — do not re-read as delivery): "whole-period-feasible"
keeps sampled demand futures whose lifetime commitment C is geologically deliverable
(2030-2179), holding every 2050 target fixed; the one cell with NO admissible future
(Canada x IPCC High x Gompertz) is scenario-infeasible and its target is removed from
the feasible world (annotated). 2050 delivery itself never binds before 2055 (Fig 3).
"""
from __future__ import annotations
import hashlib, traceback
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde

REPO = Path(__file__).resolve().parents[3]
G1 = REPO / "01_growth_model" / "output" / "v8_2026-06-01"
G2 = REPO / "02_co2block_screening" / "output"
G3 = REPO / "03_feasible_growth" / "output"
FIG = REPO / "05_paper_plots" / "v6" / "output" / "figures"
TAB = REPO / "05_paper_plots" / "v6" / "output" / "tables"
SIFIG = REPO / "05_paper_plots" / "v6" / "si_figures"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)
SIFIG.mkdir(parents=True, exist_ok=True)
S2F = G2 / "step2_screening" / "ascending" / "final"
S2D = G2 / "step2_screening" / "descending" / "final"

mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.0, "axes.titlesize": 8.0, "axes.labelsize": 7.0,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "figure.dpi": 160, "savefig.dpi": 400, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})
# PNAS legibility floor: no annotation below this point size at final print width.
FS_MIN = 6.0
# Figures carry NO baked-in "Figure N | ..." banner and NO on-figure footnote block:
# the journal sets the figure number, and the title + notes live in the LaTeX caption
# (manuscript/main.tex). This removes figure<->caption duplication and keeps all text >= FS_MIN.
SCEN = ["minimum", "growth10", "ipcc_low", "policy", "us1gt", "reference", "ipcc_high", "maximum"]
SL = {"minimum": "Minimum", "growth10": "Growth 10%", "ipcc_low": "IPCC Low", "policy": "Policy",
      "us1gt": "US 1 Gt", "reference": "Reference", "ipcc_high": "IPCC High", "maximum": "Maximum"}
TOL = {"minimum": "#888888", "growth10": "#DDCC77", "ipcc_low": "#88CCEE", "policy": "#44AA99",
       "us1gt": "#117733", "reference": "#332288", "ipcc_high": "#AA4499", "maximum": "#CC6677"}
RA = {"Middle East": "ME", "Brazil": "BR", "US": "US", "EU": "EU", "Australia": "AU",
      "Indonesia": "ID", "China": "CN", "UK": "UK", "Canada": "CA", "Thailand": "TH"}
MC_L, MC_G = "#3B6FB6", "#E8833A"
PASS_C, FAIL_C = "#2E7D32", "#C62828"
N_BOOT = 10_000
KEY: list[tuple[str, object, str]] = []


def savefig(fig, stem):
    fig.savefig(FIG / f"{stem}.png", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig); print(f"  saved {stem}")


def savefig_si(fig, stem):
    fig.savefig(SIFIG / f"{stem}.png", bbox_inches="tight")
    fig.savefig(SIFIG / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig); print(f"  saved {stem}")


def _seed(*parts) -> int:
    return int.from_bytes(hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:4], "little")


def cap_by_region():
    bp = pd.read_csv(G2 / "step1_precompute" / "precompute" / "basin_period_resource_long.csv")
    return bp, bp[bp.Period_yr == 150].groupby("ScreeningCase")["Max_Capacity_Gt"].sum()


def first_binding_year() -> tuple[int, int, int]:
    sw = pd.read_csv(S2F / "case_model_shortfall_windows.csv")
    fail = sw[~sw.whole_period_pass]
    return int(fail.first_shortfall_year.min()), int((fail.first_shortfall_year <= 2050).sum()), len(sw)


# ============================================ whole-period-feasibility bootstrap (CRN; == V4)
def _lifetime_ceilings() -> pd.DataFrame:
    pit = pd.read_csv(G3 / "part3_input_table.csv")
    rec = pd.read_csv(G3 / "smooth_reconstruction_summary.csv")
    grp = ["scenario", "country", "model"]
    out = pd.DataFrame({
        "passed": pit.groupby(grp)["whole_period_pass"].all(),
        "C_scenario": pit.groupby(grp)["C_scenario_gt"].first(),
    })
    rmin = rec.groupby(grp)["C_recon_order_gt"].min().reindex(out.index)
    rmean = rec.groupby(grp)["C_recon_order_gt"].mean().reindex(out.index)
    out["C_cap_min"] = np.where(out.passed, out.C_scenario, rmin)
    out["C_cap_mean"] = np.where(out.passed, out.C_scenario, rmean)
    return out


_BOOT_CACHE: dict | None = None


def feasible_bootstrap() -> dict:
    """Identical to V4 (same seeds, same draws): see 05_paper_plots/v4/README.md."""
    global _BOOT_CACHE
    if _BOOT_CACHE is not None:
        return _BOOT_CACHE
    caps = _lifetime_ceilings()
    draws: dict[tuple, dict] = {}
    summary, removed = [], []
    for s in SCEN:
        samp = pd.read_csv(G1 / f"samples_{s}.csv", usecols=["Country", "Model", "C_sampled", "rate_2050"])
        for m in ["Logistic", "Gompertz"]:
            sm = samp[samp.Model == m]
            if sm.empty:
                continue
            tot = {"targeted": np.zeros(N_BOOT), "min": np.zeros(N_BOOT), "mean": np.zeros(N_BOOT)}
            nc = {"targeted": 0, "min": 0, "mean": 0}
            for c in sorted(sm.Country.unique()):
                pool = sm.loc[sm.Country == c, ["C_sampled", "rate_2050"]].to_numpy()
                rates = pool[:, 1]
                idx = np.random.default_rng(_seed(s, c, m)).integers(0, len(rates), N_BOOT)
                tot["targeted"] += rates[idx]; nc["targeted"] += 1
                for kind in ["min", "mean"]:
                    cap = caps["C_cap_" + kind].get((s, c, m), np.inf)
                    keep = pool[:, 0] <= cap * (1 + 1e-12)
                    if keep.all():
                        tot[kind] += rates[idx]; nc[kind] += 1
                    elif keep.any():
                        kept = rates[keep]
                        idx2 = np.random.default_rng(_seed(s, c, m, "feasible", kind)).integers(0, len(kept), N_BOOT)
                        tot[kind] += kept[idx2]; nc[kind] += 1
                    else:
                        removed.append(dict(scenario=s, model=m, country=c, ceiling_kind=kind,
                                            min_C_sampled=float(pool[:, 0].min()), C_cap_gt=float(cap),
                                            targeted_rate_2050=float(np.median(rates))))
            draws[(s, m)] = tot
            for kind in ["min", "mean"]:
                mu, mf = float(np.median(tot["targeted"])), float(np.median(tot[kind]))
                summary.append(dict(
                    scenario=s, model=m, ceiling_kind=kind, med_targeted=round(mu, 6),
                    med_feasible=round(mf, 6), reduction_pct=round((mu - mf) / mu * 100, 4) if mu else 0.0,
                    n_countries_targeted=nc["targeted"], n_countries_feasible=nc[kind],
                    removed_countries=";".join(r["country"] for r in removed
                                               if r["scenario"] == s and r["model"] == m and r["ceiling_kind"] == kind)))
    summ = pd.DataFrame(summary)
    rem = pd.DataFrame(removed)
    summ.to_csv(TAB / "fig4_feasible_bootstrap_summary.csv", index=False)
    rem.to_csv(TAB / "fig4_removed_cells.csv", index=False)
    _BOOT_CACHE = dict(draws=draws, summary=summ, removed=rem)
    return _BOOT_CACHE


# ================================================================= FIGURE 1 (17.8 x 8.4 cm)
def fig1():
    gc = pd.read_csv(G1 / "global_central.csv")
    sc = pd.read_csv(G1 / "scalars_central.csv")
    # Roomy 2-row layout (v4-main style): a trajectories + b fork on top, c full-width raincloud below.
    fig = plt.figure(figsize=(7.2, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.3], width_ratios=[1.5, 1.0],
                          hspace=0.40, wspace=0.30, left=0.075, right=0.985, top=0.95, bottom=0.155)
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1]); axC = fig.add_subplot(gs[1, :])

    # a — global stress-pathway trajectories (resource-exhaustion demand screened in Figs 3-4)
    HERO = {"minimum", "reference", "maximum"}
    for s in SCEN:
        gcs = gc[(gc.Scenario == s) & (gc.Year.between(2030, 2180))]
        if s in HERO:
            for m, ls in [("Logistic", "-"), ("Gompertz", (0, (3.5, 2)))]:
                gm = gcs[gcs.Model == m].sort_values("Year")
                axA.plot(gm.Year, gm.rate_central_global, color=TOL[s], lw=1.6, ls=ls, zorder=5)
        else:
            gm = gcs[gcs.Model == "Logistic"].sort_values("Year")
            axA.plot(gm.Year, gm.rate_central_global, color=TOL[s], lw=0.9, alpha=0.85, zorder=3)
    axA.axvline(2050, color="#666", lw=0.7, ls=(0, (2, 2)), alpha=0.6)
    axA.text(2051, 6.5e-2, "2050", fontsize=FS_MIN, color="#666", va="bottom")
    axA.set_yscale("log"); axA.set_xlim(2030, 2180); axA.set_ylim(5e-2, 1e3)
    axA.set_xticks([2050, 2100, 2150])
    axA.set_xlabel("Year", fontsize=7.0); axA.set_ylabel("Global storage rate (Gt yr$^{-1}$)", fontsize=7.0)
    axA.set_title("a   Demand trajectories", loc="left", fontweight="bold", fontsize=8.0)

    # b — country-level growth-form fork (peak years, one point per country x scenario)
    sb = sc.copy()
    sb["d2050"] = sb[["rate_2050_L", "rate_2050_G"]].mean(axis=1)
    LO, HI = 2063, 2175
    xl = sb.tp_L.clip(LO + 1, HI - 2); yl = sb.tp_G.clip(LO + 1, HI - 2)
    # "both peak post-2100" quadrant: the late-peaking pathways that load the subsurface latest
    axB.add_patch(plt.Rectangle((2100, 2100), HI - 2100, HI - 2100, facecolor="#C62828", alpha=0.05, lw=0, zorder=0))
    axB.text(2107, 2146, "both peak\npost-2100", fontsize=5.6, color="#b06", ha="left", va="top",
             style="italic", linespacing=1.05, zorder=1)
    axB.plot([LO, HI], [LO, HI], color="#bbb", lw=0.8, ls=(0, (4, 3)), zorder=1)
    axB.axvline(2100, color="#d8d8d8", lw=0.7, ls=(0, (2, 2)), zorder=1)
    axB.axhline(2100, color="#d8d8d8", lw=0.7, ls=(0, (2, 2)), zorder=1)
    axB.scatter(xl, yl, s=8 + np.sqrt(sb.d2050.clip(lower=0)) * 52, c=[TOL[s] for s in sb.Scenario],
                alpha=0.80, edgecolor="white", linewidth=0.3, zorder=4)
    # direct-label one representative (max-demand) point per major jurisdiction
    off = {"US": (6, -8), "China": (7, -6), "Canada": (-15, 5), "Middle East": (7, 1), "Brazil": (7, 0)}
    for c in ["US", "China", "Canada", "Brazil", "Middle East"]:
        sub = sb[sb.Country == c]
        if sub.empty:
            continue
        r = sub.loc[sub.d2050.idxmax()]
        t = axB.annotate(RA[c], (min(max(r.tp_L, LO + 1), HI - 2), min(max(r.tp_G, LO + 1), HI - 2)),
                         textcoords="offset points", xytext=off[c], fontsize=5.8, fontweight="bold", color="#222", zorder=6)
        t.set_path_effects([pe.withStroke(linewidth=1.6, foreground="white")])
    n_later = int((sb.tp_G > sb.tp_L).sum())
    axB.text(0.96, 0.06, f"{n_later}/{len(sb)} above diagonal:\nGompertz peaks later\n"
             f"(median +{(sb.tp_G - sb.tp_L).median():.0f} yr)", transform=axB.transAxes, fontsize=FS_MIN,
             color="#555", ha="right", va="bottom", style="italic", linespacing=1.25)
    axB.set_xlim(LO, HI); axB.set_ylim(LO, HI)
    axB.set_xticks([2080, 2120, 2160]); axB.set_yticks([2080, 2120, 2160])
    axB.set_xlabel("Logistic peak year", fontsize=7.0); axB.set_ylabel("Gompertz peak year", fontsize=7.0)
    axB.set_title("b   Growth-form fork", loc="left", fontweight="bold", fontsize=8.0)
    KEY.extend([("fig1b_pathways", len(sb), "scalars_central.csv"),
                ("fig1b_gompertz_peaks_later", n_later, "tp_G > tp_L"),
                ("fig1b_median_peak_delay_yr", round(float((sb.tp_G - sb.tp_L).median()), 1), "median tp_G - tp_L")])

    # c — 2050 cumulative commitment: MC raincloud (cloud + rain dots + box) + deterministic stress
    #     inputs, summed over the MC-covered jurisdictions
    rob = pd.read_csv(G1 / "robustness_sample_paths.csv",
                      usecols=["Scenario", "Country", "Model", "k_index", "Year", "rate"])
    mc = (rob[rob.Year.between(2030, 2050)].groupby(["Scenario", "Model", "k_index"])["rate"].sum()
          .reset_index(name="cum"))
    cov = rob.groupby(["Scenario", "Model"])["Country"].agg(set)
    sps = pd.read_csv(S2F / "case_model_screened_paths.csv")
    s50 = sps[sps.year == 2050]
    keep = [r.country in cov.get((r.scenario, r.model), set()) for r in s50.itertuples()]
    det = s50[keep].groupby(["scenario", "model"])["raw_cumulative_gt"].sum()
    rng = np.random.default_rng(1)
    for i, s in enumerate(SCEN):
        v = mc.loc[mc.Scenario == s, "cum"].to_numpy(); v = v[v > 0]
        if v.size < 5:
            continue
        col = TOL[s]; lv = np.log10(v)
        xs = np.linspace(lv.min(), lv.max(), 160)
        dens = gaussian_kde(lv)(xs); dens = dens / dens.max() * 0.36
        axC.fill_between(10 ** xs, i + 0.05, i + 0.05 + dens, color=col, alpha=0.20, lw=0, zorder=2)
        axC.scatter(v, i - 0.08 - rng.random(v.size) * 0.24, s=1.3, color=col, alpha=0.13,
                    linewidths=0, zorder=1, rasterized=True)
        p5, p25, p50, p75, p95 = np.percentile(v, [5, 25, 50, 75, 95])
        axC.plot([p5, p95], [i, i], color=col, lw=1.0, zorder=4)
        axC.plot([p25, p75], [i, i], color=col, lw=3.4, zorder=4, solid_capstyle="round")
        axC.plot([p50, p50], [i - 0.085, i + 0.085], color="white", lw=1.3, zorder=6, solid_capstyle="butt")
        for m, mk in [("Logistic", "o"), ("Gompertz", "^")]:
            dv = det.get((s, m), np.nan)
            if pd.notna(dv):
                axC.scatter([dv], [i], marker=mk, s=30, facecolor=col, edgecolor="#111", linewidth=0.8, zorder=7)
        if s == "reference":
            KEY.append(("fig1c_reference_median_cum2050_gt", round(float(p50), 1), "robustness_sample_paths.csv"))
    axC.set_yticks(range(len(SCEN))); axC.set_yticklabels([SL[s] for s in SCEN], fontsize=6.8)
    axC.set_ylim(-0.55, len(SCEN) - 0.25)
    axC.set_xscale("log"); axC.set_xlim(3, 180)
    axC.set_xlabel("Cumulative CO$_2$ stored by 2050 (GtCO$_2$)", fontsize=7.0)
    axC.set_title("c   2050 cumulative commitment", loc="left", fontweight="bold", fontsize=8.0)
    axC.grid(axis="x", color="#eee", lw=0.5)

    fig.legend(handles=[Line2D([0], [0], color="#444", ls="-", marker="o", lw=1.4, markersize=4.5, label="Logistic (line / ● input)"),
                        Line2D([0], [0], color="#444", ls=(0, (3.5, 2)), marker="^", lw=1.4, markersize=4.5, label="Gompertz (line / ▲ input)"),
                        Line2D([0], [0], color="#888", lw=0, marker="|", markersize=8, markeredgewidth=1.4, label="MC median"),
                        Line2D([0], [0], color="#888", lw=3.2, label="MC 25–75%"),
                        Line2D([0], [0], color="#888", lw=1.0, label="MC 5–95%")],
               loc="lower center", bbox_to_anchor=(0.5, 0.075), ncol=5, frameon=False,
               fontsize=FS_MIN, handlelength=1.9, columnspacing=1.2)
    fig.legend(handles=[Line2D([0], [0], color=TOL[s], lw=2.2, label=SL[s]) for s in SCEN],
               loc="lower center", bbox_to_anchor=(0.5, 0.025), ncol=8, frameon=False,
               fontsize=FS_MIN, handlelength=1.2, columnspacing=1.0)
    savefig(fig, "fig1_demand")


# ================================================================= FIGURE 2 (17.8 x 7.6 cm)
def fig2():
    bp, cap = cap_by_region()
    b150 = bp[bp.Period_yr == 150]
    rate = b150.groupby("ScreeningCase")["Region_Q_Mt_yr"].sum().sort_values()
    order = list(rate.index)
    _blues = mpl.colormaps["Blues"]
    RC2 = {r: _blues(0.34 + 0.62 * i / (len(order) - 1)) for i, r in enumerate(order)}
    fig = plt.figure(figsize=(7.0, 3.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.92, 0.80, 1.18],
                          wspace=0.13, left=0.085, right=0.975, top=0.945, bottom=0.135)
    axBar = fig.add_subplot(gs[0]); axBox = fig.add_subplot(gs[1], sharey=axBar)
    axB = fig.add_subplot(gs[2])

    y = np.arange(len(order))
    axBar.barh(y, rate[order].values, color=[RC2[r] for r in order], height=0.7, edgecolor="white", linewidth=0.4)
    axBar.set_yticks(y); axBar.set_yticklabels(order, fontsize=6.0)
    axBar.set_xscale("log"); axBar.set_xlim(150, 2e4)
    axBar.set_xlabel("Deliverable rate at 150 yr (Mt yr$^{-1}$)", fontsize=6.2)
    axBar.set_title("a   Regional supply", loc="left", fontweight="bold", fontsize=7.6)
    nbreg = b150.groupby("ScreeningCase")["Region_no"].nunique()
    for yi, r in zip(y, order):
        axBar.text(rate[r] * 1.15, yi, f"{int(nbreg[r])}", va="center", fontsize=5.5, color="#777")
    rng = np.random.default_rng(0)
    bd = [b150[b150.ScreeningCase == r]["Region_Q_Mt_yr"].values for r in order]
    bx = axBox.boxplot(bd, positions=y, vert=False, widths=0.55, patch_artist=True, showfliers=False,
                       medianprops=dict(color="#222", lw=0.9),
                       whiskerprops=dict(color="#999", lw=0.7), capprops=dict(color="#999", lw=0.7))
    for patch, r in zip(bx["boxes"], order):
        patch.set_facecolor(RC2[r]); patch.set_alpha(0.40); patch.set_edgecolor("none")
    for r, yi in zip(order, y):
        v = b150[b150.ScreeningCase == r]["Region_Q_Mt_yr"].values
        w = b150[b150.ScreeningCase == r]["Number_of_Wells"].values
        axBox.scatter(v, yi + (rng.random(len(v)) - 0.5) * 0.45, s=np.sqrt(np.clip(w, 1, None)) * 1.3,
                      color=RC2[r], alpha=0.5, edgecolor="white", linewidth=0.2, zorder=4)
    axBox.axvline(1.0, color="#C62828", lw=0.7, ls=":")
    axBox.set_xscale("log"); axBox.set_xlim(0.2, 2e4)
    axBox.set_xlabel("Basin rate (Mt yr$^{-1}$)", fontsize=6.2)
    plt.setp(axBox.get_yticklabels(), visible=False); axBox.tick_params(axis="y", length=0)
    axBar.set_yticks(y); axBar.set_yticklabels(order, fontsize=6.0)

    bdf = (b150.groupby("Region_no").agg(cap=("Max_Capacity_Gt", "sum"), reg=("ScreeningCase", "first"))
                .sort_values("cap", ascending=False).reset_index())
    nb = len(bdf); cumv = np.cumsum(bdf.cap.values) / bdf.cap.sum() * 100
    xf = (np.arange(nb) + 1) / nb * 100; cumd = np.concatenate([[0], cumv])
    axB.plot([0, 100], [0, 100], color="#bbb", lw=0.7, ls=(0, (3, 3)))
    axB.plot(np.concatenate([[0], xf]), cumd, color="#555", lw=0.9, zorder=2)
    axB.scatter(xf, cumv, c=[RC2[r] for r in bdf.reg], s=9, zorder=4, edgecolor="white", linewidth=0.2)
    for k, ty in [(10, 22), (20, 36), (50, 50)]:   # text order matches point order -> leaders stay parallel
        xk = k / nb * 100; yk = cumd[k]
        axB.plot([xk, xk], [0, yk], color="#999", lw=0.5, ls=(0, (2, 2)), zorder=1)   # Lorenz drop-lines
        axB.plot([0, xk], [yk, yk], color="#999", lw=0.5, ls=(0, (2, 2)), zorder=1)
        axB.scatter([xk], [yk], s=18, facecolor="none", edgecolor="#333", lw=0.7, zorder=5)
        axB.annotate(f"top {k}: {yk:.1f}%", xy=(xk, yk), xytext=(54, ty), fontsize=FS_MIN, color="#333",
                     ha="left", va="center", arrowprops=dict(arrowstyle="-", lw=0.4, color="#888", alpha=0.7))
        KEY.append((f"fig2b_top{k}_capacity_share_pct", round(float(yk), 1), "basin_period_resource_long.csv"))
    for i in range(4):
        axB.annotate(RA[bdf.reg[i]], (xf[i], cumv[i]), textcoords="offset points", xytext=(2.5, -6),
                     fontsize=FS_MIN, color=RC2[bdf.reg[i]], fontweight="bold")
    axB.set_xlim(0, 100); axB.set_ylim(0, 100)
    axB.set_xlabel("Basins ranked by capacity (%)", fontsize=6.4); axB.set_ylabel("Cumulative capacity (%)", fontsize=6.4)
    axB.set_title("b   Resource concentration", loc="left", fontweight="bold", fontsize=7.6)
    KEY.extend([("fig2_total_capacity_gt_150yr", round(float(b150.Max_Capacity_Gt.sum()), 0), "Period_yr==150 sum"),
                ("fig2_total_rate_mt_yr_150yr", round(float(b150.Region_Q_Mt_yr.sum()), 0), "Period_yr==150 sum"),
                ("fig2_n_basins", nb, "basin_period_resource_long.csv")])

    savefig(fig, "fig2_supply")


# ================================================================= FIGURE 3 (17.8 x 10.9 cm, centerpiece)
def fig3():
    sw = pd.read_csv(S2F / "case_model_shortfall_windows.csv")
    sw["bind"] = np.where(sw.whole_period_pass, 2180, sw.first_shortfall_year)
    sw["shortfall_pct"] = (sw.cumulative_demand_full_gt - sw.cumulative_screened_full_gt) / sw.cumulative_demand_full_gt * 100
    sp = pd.read_csv(S2F / "case_model_screened_paths.csv")
    _, cap = cap_by_region()
    sp["util"] = sp["screened_cumulative_gt"] / sp["country"].map(cap) * 100
    bind0, n_b50, n_tot = first_binding_year()          # derived, not hardcoded (V3: 2055 / "0 of 150")
    KEY.extend([("fig3_first_binding_year", bind0, "case_model_shortfall_windows.csv"),
                ("fig3_cases_binding_by_2050", n_b50, "first_shortfall_year <= 2050"),
                ("fig3_n_cases", n_tot, "ascending shortfall table")])
    rec = pd.read_csv(G3 / "smooth_reconstruction_summary.csv")          # panels d/e: the binding cases
    asc = rec[rec.ordering == "ascending"].copy()
    asc["ratio"] = asc.C_recon_order_gt / asc.C_scenario_gt
    med = float(asc.ratio.median())

    fig = plt.figure(figsize=(7.4, 7.4))
    gsA = fig.add_gridspec(1, 1, left=0.075, right=0.45, top=0.92, bottom=0.455)   # panel a (own block)
    axA = fig.add_subplot(gsA[0])
    gsR = fig.add_gridspec(2, 2, left=0.575, right=0.925, top=0.92, bottom=0.455,  # right group, with a clear gutter
                           height_ratios=[0.82, 1.12], hspace=0.62, wspace=0.20)
    axB = fig.add_subplot(gsR[0, :]); axCL = fig.add_subplot(gsR[1, 0]); axCG = fig.add_subplot(gsR[1, 1], sharey=axCL)
    gsB = fig.add_gridspec(1, 2, left=0.075, right=0.84, top=0.36, bottom=0.10,    # bottom row: d (collapse) + e (by form)
                           width_ratios=[1.55, 0.42], wspace=0.18)
    axD = fig.add_subplot(gsB[0]); axE = fig.add_subplot(gsB[1])

    # A: runway barcode, SPLIT into Logistic (top) and Gompertz (bottom) blocks so the two growth
    #     forms can be compared directly (Logistic binds more often, and usually earlier).
    axA.axvspan(2030, 2050, color="#2E7D32", alpha=0.08, lw=0, zorder=0)   # before binding: runway
    axA.axvspan(bind0, 2100, color="#C62828", alpha=0.07, lw=0, zorder=0)  # binding window (starts at the DATA's first bind)
    axA.axvspan(2100, 2180, color="#C62828", alpha=0.03, lw=0, zorder=0)   # late tail
    GAP = 9; y = 0; binfo = []
    for mod in ["Gompertz", "Logistic"]:                                   # bottom block = Gompertz, top = Logistic
        sub = sw[sw.model == mod].sort_values("bind").reset_index(drop=True)
        ys = np.arange(len(sub)) + y
        axA.add_collection(LineCollection([[(2030, yy), (br.bind, yy)] for yy, (_, br) in zip(ys, sub.iterrows())],
                                          colors=[TOL[br.scenario] for _, br in sub.iterrows()], linewidths=1.0, zorder=2))
        nb = ys[sub.whole_period_pass.values]
        axA.scatter([2179] * len(nb), nb, marker=">", s=4, color="#2E7D32", zorder=4, clip_on=False)
        binfo.append((mod, ys[0], ys[-1], int((~sub.whole_period_pass).sum())))
        y = ys[-1] + 1 + GAP
    ytop = y - GAP
    axA.axhline((binfo[0][2] + binfo[1][1]) / 2, color="#bbb", lw=0.7, zorder=3)   # divider between blocks
    axA.axvline(2050, color="#222", lw=1.0, ls=(0, (3, 2)), zorder=5)
    axA.set_xlim(2030, 2180); axA.set_ylim(-3, ytop + 1)
    axA.set_yticks([(b[1] + b[2]) / 2 for b in binfo])                     # block centres carry the model + count labels
    axA.set_yticklabels([f"{b[0]}\n{b[3]} bind" for b in binfo], fontsize=6.0, fontweight="bold")
    for tl, b in zip(axA.get_yticklabels(), binfo):
        tl.set_color(MC_G if b[0] == "Gompertz" else MC_L)
    axA.tick_params(axis="y", length=0)
    axA.text(2053, ytop * 0.30, f"{n_b50} of {n_tot} bind\nby 2050", fontsize=6.0, ha="left", va="center",
             color="#111", fontweight="bold", linespacing=0.95, zorder=6,
             bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#222", lw=0.5, alpha=0.92))
    axA.set_xlabel("Year")
    axA.set_ylabel(f"{n_tot} cases (country × scenario × form)", fontsize=6.2)
    axA.set_title("a   Feasibility runway", loc="left", fontweight="bold", fontsize=7.6)
    KEY.extend([(f"fig3_bind_total_{b[0]}", b[3], "whole_period_pass == False") for b in binfo])

    # B: shortfall severity by growth form (drop bind-year): x = % of cumulative demand unmet by 2180,
    #     one dot per binding case (coloured by scenario), median tick + IQR bar per form.
    bnd = sw[~sw.whole_period_pass]; rngb = np.random.default_rng(3)
    for i, mod in enumerate(["Gompertz", "Logistic"]):
        d = bnd[bnd.model == mod]
        p25, p50, p75 = np.percentile(d.shortfall_pct, [25, 50, 75])
        axB.plot([p25, p75], [i, i], color="#999", lw=3.0, alpha=0.4, solid_capstyle="round", zorder=2)
        axB.scatter(d.shortfall_pct, i + (rngb.random(len(d)) - 0.5) * 0.5, c=[TOL[s] for s in d.scenario],
                    s=13, edgecolor="white", linewidth=0.3, zorder=3)
        axB.plot([p50, p50], [i - 0.30, i + 0.30], color="#222", lw=1.4, zorder=4)
        axB.text(p50, i + 0.40, f"med {p50:.0f}%", fontsize=FS_MIN, ha="center", va="bottom", color="#333")
        axB.text(99, i + 0.40, f"{len(d)} bind", fontsize=FS_MIN, ha="right", va="bottom",
                 color=(MC_L if mod == "Logistic" else MC_G), fontweight="bold")
        KEY.append((f"fig3b_median_shortfall_pct_{mod}", round(float(p50), 1), "binding cases, demand unmet by 2180"))
    axB.set_yticks([0, 1]); axB.set_yticklabels(["Gompertz", "Logistic"], fontsize=6.2)
    for tl, mod in zip(axB.get_yticklabels(), ["Gompertz", "Logistic"]):
        tl.set_color(MC_G if mod == "Gompertz" else MC_L); tl.set_fontweight("bold")
    axB.set_xlim(0, 100); axB.set_ylim(-0.55, 1.78)
    axB.set_xlabel("Shortfall — cumulative demand unmet by 2180 (%)", fontsize=6.5)
    axB.set_title("b   Shortfall severity", loc="left", fontweight="bold", fontsize=7.6)
    axB.grid(axis="x", color="#eee", lw=0.5)

    # C: capacity used (max across cases, %), Logistic | Gompertz side by side: shows Logistic
    #     front-loads capacity use (high by 2100) while Gompertz defers it; both saturate by 2180.
    chk = [2050, 2100, 2179]; xl = ["2050", "2100", "2180"]
    rows = ["Thailand", "Indonesia", "Middle East", "Brazil", "EU", "Australia", "UK", "US", "Canada", "China"]
    def ugrid(model):
        return np.array([[min(sp[(sp.country == c) & (sp.year == yv) & (sp.model == model)]["util"].max(), 100)
                          for yv in chk] for c in rows])
    im = None
    for ax, model in [(axCL, "Logistic"), (axCG, "Gompertz")]:
        U = ugrid(model)
        im = ax.imshow(U, aspect="auto", cmap="OrRd", vmin=0, vmax=100, extent=[0, 3, 0, len(rows)])
        ax.set_xticks([0.5, 1.5, 2.5]); ax.set_xticklabels(xl, fontsize=5.4)
        ax.set_title(model, fontsize=6.5, fontweight="bold", color=(MC_L if model == "Logistic" else MC_G), pad=2)
        for i in range(len(rows)):
            for j in range(3):
                v = U[len(rows) - 1 - i, j]
                ax.text(j + 0.5, i + 0.5, f"{v:.0f}", ha="center", va="center", fontsize=5.2,
                        color="white" if v > 55 else "#333")
    axCL.set_yticks(np.arange(len(rows)) + 0.5); axCL.set_yticklabels(rows[::-1], fontsize=5.4)
    plt.setp(axCG.get_yticklabels(), visible=False); axCG.tick_params(axis="y", length=0)
    fig.text(0.75, axCL.get_position().y1 + 0.038, "c   Capacity used (%)",
             fontsize=7.6, fontweight="bold", ha="center")
    cpos = axCG.get_position()
    cax = fig.add_axes([0.94, cpos.y0, 0.013, cpos.height * 0.85])
    cb = fig.colorbar(im, cax=cax); cb.set_label("capacity used (%)", fontsize=FS_MIN); cb.ax.tick_params(labelsize=5.2, length=2)

    # D: commitment collapse: every binding case as goal C -> reconstructed feasible C
    LIM = (12, 7000)
    axD.plot(LIM, LIM, color="#999", lw=0.9, ls=(0, (4, 3)), zorder=1)
    axD.text(1250, 1600, "no loss (recon = goal)", fontsize=FS_MIN, color="#888", rotation=23,
             ha="center", va="bottom", rotation_mode="anchor")
    gx = np.array(LIM)
    # retained-share reference guides (recon = f x goal) give the median line visual context;
    # the 10% label anchors mid-line so it stays clear of the AU/US country labels at right
    for f, xa in [(0.5, LIM[1]), (0.25, LIM[1]), (0.1, 2600)]:
        axD.plot(gx, gx * f, color="#cfcfcf", lw=0.5, ls=(0, (2, 2)), zorder=1.4)
        axD.text(xa * 0.93, xa * 0.93 * f * 0.93, f"{int(f*100)}%", fontsize=5.2, color="#aaa",
                 ha="right", va="top", rotation=23, rotation_mode="anchor", zorder=1.4)
    axD.plot(gx, gx * med, color="#C62828", lw=1.0, alpha=0.85, zorder=2)
    axD.text(900, 900 * med * 0.90, f"median: {med*100:.0f}% of goal", fontsize=FS_MIN, color="#C62828",
             rotation=23, ha="center", va="top", rotation_mode="anchor")
    for m, mkk in [("Logistic", "o"), ("Gompertz", "^")]:
        d = asc[asc.model == m]
        axD.scatter(d.C_scenario_gt, d.C_recon_order_gt, s=30, marker=mkk,
                    c=[TOL[s] for s in d.scenario], edgecolor="white", linewidth=0.4, alpha=0.92, zorder=4)
    inf = asc[asc.g_recon_order > asc.g_cap]                 # the scenario-infeasible cell (Fig 4)
    for _, r in inf.iterrows():
        axD.scatter([r.C_scenario_gt], [r.C_recon_order_gt], s=100, facecolor="none",
                    edgecolor="#C62828", linewidth=0.9, zorder=5)
        axD.annotate("Canada · IPCC High · Gompertz:\nreconstruction needs g > cap\n(the scenario-infeasible cell, Fig. 4)",
                     xy=(r.C_scenario_gt, r.C_recon_order_gt), xytext=(15.5, 1500), fontsize=FS_MIN,
                     color="#C62828", va="center", linespacing=1.25,
                     arrowprops=dict(arrowstyle="->", lw=0.5, color="#C62828"))
    lab = asc.sort_values("C_scenario_gt", ascending=False).groupby("country").head(1)
    for _, r in lab.iterrows():
        t = axD.annotate(RA.get(r.country, r.country), (r.C_scenario_gt, r.C_recon_order_gt),
                         textcoords="offset points", xytext=(4, 4), fontsize=6.2, fontweight="bold", color="#222")
        t.set_path_effects([pe.withStroke(linewidth=1.8, foreground="white")])
    axD.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor="#888", markersize=4.5, label="Logistic"),
                        Line2D([0], [0], marker="^", color="w", markerfacecolor="#888", markersize=4.5, label="Gompertz")],
               loc="lower right", fontsize=FS_MIN, frameon=False, handletextpad=0.3, borderaxespad=0.25)
    axD.set_xscale("log"); axD.set_yscale("log"); axD.set_xlim(*LIM); axD.set_ylim(*LIM)
    axD.set_xlabel("Goal resource $C$ — scenario design (Gt)", fontsize=6.6)
    axD.set_ylabel("Reconstructed feasible $C$ (Gt)", fontsize=6.6)
    axD.set_title("d   Commitment collapse", loc="left", fontweight="bold", fontsize=7.6)
    axD.grid(color="#f3f3f3", lw=0.4); axD.tick_params(labelsize=6.0)
    KEY.extend([("fig3d_n_binding_cases", len(asc), "ascending reconstructions"),
                ("fig3d_recon_over_goal_median_pct", round(med * 100, 1), "ascending reconstructions")])

    # E: retained share of the goal resource, split by growth form
    rngE = np.random.default_rng(7)
    for i, (m, cc) in enumerate([("Logistic", MC_L), ("Gompertz", MC_G)]):
        d = asc[asc.model == m]
        axE.scatter(i + (rngE.random(len(d)) - 0.5) * 0.55, d.ratio * 100, s=16,
                    c=[TOL[s] for s in d.scenario], edgecolor="white", linewidth=0.3, zorder=3)
        mm = float(d.ratio.median()) * 100
        axE.plot([i - 0.30, i + 0.30], [mm, mm], color="#222", lw=1.3, zorder=4)
        xt, ha_ = (i + 0.36, "left") if m == "Logistic" else (i - 0.36, "right")   # beside the tick, clear of dots
        yt, va_ = (mm + 1.5, "bottom") if m == "Logistic" else (mm - 1.5, "top")   # offset off the red median dash
        axE.text(xt, yt, f"{mm:.0f}%", fontsize=FS_MIN, ha=ha_, va=va_, color="#222")
        KEY.append((f"fig3e_median_retained_pct_{m}", round(mm, 1), "ascending reconstructions"))
    axE.axhline(med * 100, color="#C62828", lw=0.9, ls=(0, (4, 2)), zorder=2)
    axE.set_xlim(-0.6, 1.6); axE.set_ylim(0, 100)
    axE.set_xticks([0, 1]); axE.set_xticklabels(["Logistic", "Gompertz"], fontsize=6.2, fontweight="bold")
    for tl, cc in zip(axE.get_xticklabels(), [MC_L, MC_G]):
        tl.set_color(cc)
    axE.set_ylabel("retained share of goal $C$ (%)", fontsize=6.4)
    axE.set_title("e   Retained share", loc="left", fontweight="bold", fontsize=7.6)
    axE.grid(axis="y", color="#f2f2f2", lw=0.4); axE.tick_params(labelsize=6.0)

    fig.legend(handles=[Line2D([0], [0], color=TOL[s], lw=2.4, label=SL[s]) for s in SCEN],
               loc="lower center", bbox_to_anchor=(0.5, 0.012), ncol=8, frameon=False,
               fontsize=FS_MIN, handlelength=1.2, columnspacing=1.0)
    savefig(fig, "fig3_test")


# ================================================================= FIGURE 4 (17.8 x 8.6 cm)
def fig4():
    bb = feasible_bootstrap()
    draws, summ = bb["draws"], bb["summary"]
    smin = summ[summ.ceiling_kind == "min"]
    piv_t = smin.pivot(index="scenario", columns="model", values="med_targeted").reindex(SCEN)
    piv_f = smin.pivot(index="scenario", columns="model", values="med_feasible").reindex(SCEN)
    red = smin.pivot(index="scenario", columns="model", values="reduction_pct").reindex(SCEN)
    sp = pd.read_csv(S2F / "case_model_screened_paths.csv")
    bind0, _, _ = first_binding_year()

    # COMMON country set per scenario for every Logistic-vs-Gompertz comparison (IPCC High:
    # Brazil is Logistic-only; IPCC Low: Canada is Logistic-only); coverage must never be
    # read as a growth-form effect.
    cset = {s: (set(sp.loc[(sp.scenario == s) & (sp.model == "Logistic"), "country"]) &
                set(sp.loc[(sp.scenario == s) & (sp.model == "Gompertz"), "country"])) for s in SCEN}
    full = {s: set(sp.loc[sp.scenario == s, "country"]) for s in SCEN}
    spc = sp[[c in cset[s] for s, c in zip(sp.scenario, sp.country)]]
    _, cap = cap_by_region()                                                # 150-yr deliverable resource per jurisdiction

    yrs = [y for y in sorted(sp.year.unique()) if y <= 2180]

    fig = plt.figure(figsize=(6.6, 8.7))
    outer = fig.add_gridspec(5, 2, height_ratios=[1, 1, 1, 1, 0.88], hspace=0.44, wspace=0.16,
                             left=0.105, right=0.95, top=0.90, bottom=0.07)

    # ---- a: ONE CELL PER SCENARIO = cumulative trajectory (left) + 2050-rate MC marginal (right) ----
    # main: deterministic stress-path cumulative Gt (Σ common countries), L blue solid / G orange
    # dashed, shaded shortfall (delivery -> demand). marginal: MC 2050-rate distribution,
    # targeted fill / whole-period-feasible step.
    for k, s in enumerate(SCEN):
        row, col = k // 2, k % 2
        cg = outer[row, col].subgridspec(1, 2, width_ratios=[2.2, 1.0], wspace=0.05)
        axM = fig.add_subplot(cg[0]); axR = fig.add_subplot(cg[1])
        for m, cc, ls in [("Logistic", MC_L, "-"), ("Gompertz", MC_G, (0, (3, 2)))]:
            sub = spc[(spc.scenario == s) & (spc.model == m)]
            if sub.empty:
                continue
            dem = sub.groupby("year")["raw_cumulative_gt"].sum().reindex(yrs).to_numpy()
            dlv = sub.groupby("year")["screened_cumulative_gt"].sum().reindex(yrs).to_numpy()
            axM.fill_between(yrs, dlv, dem, where=dem > dlv, color=cc, alpha=0.14, lw=0, zorder=2)
            axM.plot(yrs, dlv, color=cc, ls=ls, lw=1.2, zorder=3)
        dcell = float(sum(cap.get(c, 0.0) for c in cset[s]))                # geological supply ceiling for this cell
        if dcell > 0:
            axM.axhline(dcell, color="#C62828", lw=0.7, ls=(0, (3, 2)), alpha=0.85, zorder=2.5)
            if k == 0:
                axM.text(2034, dcell * 1.35, "deliverable ceiling", fontsize=5.0, color="#C62828", va="bottom")
        axM.axvline(2050, color="#999", lw=0.6, ls=(0, (2, 2)), zorder=1)
        axM.set_yscale("log"); axM.set_xlim(2030, 2180); axM.set_ylim(1, 3e4)
        axM.set_xticks([2050, 2100, 2150]); axM.set_yticks([1, 1e2, 1e4])
        axM.grid(axis="y", color="#f4f4f4", lw=0.4, zorder=0)
        axM.set_title(SL[s], fontsize=7.0, color=TOL[s], fontweight="bold", pad=2, loc="left")
        axM.tick_params(labelsize=FS_MIN)
        if len(cset[s]) < len(full[s]):
            axM.text(0.03, 0.96, f"common n={len(cset[s])}", transform=axM.transAxes,
                     fontsize=FS_MIN, color="#777", va="top")
        if row < 3:
            axM.set_xticklabels([])
        if col == 1:
            axM.set_yticklabels([])
        if s == "maximum":
            axM.annotate("geology caps\nthis path late", xy=(2105, 9e3), xytext=(2034, 2.3e3), fontsize=FS_MIN,
                         color=MC_L, ha="left", va="center", arrowprops=dict(arrowstyle="->", lw=0.5, color=MC_L))
        # marginal: 2050-rate MC distribution (rate on y), targeted = fill / feasible = step
        hi = 1e-3
        for m in ["Logistic", "Gompertz"]:
            if (s, m) in draws:
                hi = max(hi, draws[(s, m)]["targeted"].max(), draws[(s, m)]["min"].max())
        bins = np.linspace(0, hi, 22)
        for m, cc in [("Logistic", MC_L), ("Gompertz", MC_G)]:
            if (s, m) not in draws:
                continue
            axR.hist(draws[(s, m)]["targeted"], bins=bins, orientation="horizontal", color=cc, alpha=0.26, lw=0)
            axR.hist(draws[(s, m)]["min"], bins=bins, orientation="horizontal", histtype="step", color=cc, lw=0.8)
        axR.set_ylim(0, hi); axR.set_xticks([])
        axR.yaxis.tick_right(); axR.set_yticks([0, hi]); axR.set_yticklabels(["0", f"{hi:.0f}"], fontsize=5.2)
        axR.tick_params(axis="y", length=1.5, pad=1)
        for sp_ in ["top", "left"]:
            axR.spines[sp_].set_visible(False)
    fig.text(0.020, 0.57, "Cumulative CO$_2$ stored (Gt)", rotation=90, va="center", ha="center", fontsize=7.2)
    fig.text(0.50, 0.205, "Year   (right strip of each cell = 2050 rate, Gt yr$^{-1}$)", ha="center", fontsize=6.4)
    fig.text(0.105, 0.945, "a   Demand vs delivered (all scenarios)",
             fontsize=7.6, fontweight="bold")

    # ---- b: targeted vs whole-period-feasible 2050 rate across scenarios ----
    SH = {"minimum": "Min", "growth10": "G10", "ipcc_low": "IPL", "policy": "Pol",
          "us1gt": "US1", "reference": "Ref", "ipcc_high": "IPH", "maximum": "Max"}
    axB = fig.add_subplot(outer[4, 0])
    xpos = np.arange(len(SCEN))
    for m, cc in [("Logistic", MC_L), ("Gompertz", MC_G)]:
        un = piv_t[m].to_numpy(); fe = piv_f[m].to_numpy(); rd = red[m].to_numpy()
        axB.plot(xpos, un, color=cc, lw=1.3, marker="o", ms=3.0, zorder=3)
        for x, u, f, r in zip(xpos, un, fe, rd):
            if pd.notna(u) and pd.notna(f) and r > 0.5:          # CRN makes true no-ops exactly 0
                axB.plot([x, x], [u, f], color=cc, lw=0.9, alpha=0.6, zorder=2)
                axB.scatter([x], [f], color=cc, marker="v", s=15, facecolors="white", linewidths=0.8, zorder=4)
                if r >= 5:
                    if x == len(SCEN) - 1:      # last scenario: flip inward so the label is not clipped
                        axB.text(x - 0.13, f, f"−{r:.0f}%", fontsize=FS_MIN, color=cc, va="center", ha="right")
                    else:
                        axB.text(x + 0.13, f, f"−{r:.0f}%", fontsize=FS_MIN, color=cc, va="center")
    axB.set_xticks(xpos); axB.set_xticklabels([SH[s] for s in SCEN], fontsize=FS_MIN)
    axB.set_ylim(0, None); axB.set_ylabel("2050 rate (Gt yr$^{-1}$)", fontsize=6.2)
    axB.tick_params(labelsize=FS_MIN); axB.grid(axis="y", color="#f1f1f1", lw=0.4)
    axB.set_title("b   2050-rate ladder", loc="left", fontweight="bold", fontsize=7.6)
    axB.annotate("demand scenario =\nthe main driver", xy=(6.5, 5.2), xytext=(0.2, 7.1), fontsize=FS_MIN, color="#555",
                 va="top", arrowprops=dict(arrowstyle="->", lw=0.5, color="#999"))
    # the one scenario-infeasible cell, annotated as a highlighted binding case
    rem = bb["removed"]
    rem_min = rem[rem.ceiling_kind == "min"]
    for _, rr in rem_min.iterrows():
        x = SCEN.index(rr.scenario)
        f = piv_f.loc[rr.scenario, rr.model]
        axB.annotate(f"{RA.get(rr.country, rr.country)} target removed\n(no feasible pathway)",
                     xy=(x, f), xytext=(7.25, 1.55), fontsize=FS_MIN, color=MC_G, fontweight="bold",
                     ha="right", va="center", arrowprops=dict(arrowstyle="->", lw=0.5, color=MC_G))
    for _, r in smin.iterrows():
        KEY.append((f"fig4b_{r.scenario}_{r.model}_med_targeted_gt_yr", r.med_targeted, "CRN bootstrap, seed sha256"))
        KEY.append((f"fig4b_{r.scenario}_{r.model}_reduction_pct", r.reduction_pct, "whole-period-feasible vs targeted"))

    # c: WHEN each driver acts: growth-form choice vs geology vs allocation order, on cumulative
    # stored CO2 over time, ALL evaluated on the common country sets. Growth-form effect =
    # |Logistic - Gompertz| delivered; geology effect = shortfall (demand - delivered), each form a
    # separate case; order effect = |ascending - descending|. Median across scenarios (line) +
    # full range (band, geology only). The point: growth form acts from the start, geology is ~0
    # through 2050 and only emerges in the long run.
    axC = fig.add_subplot(outer[4, 1])
    g4 = spc.groupby(["scenario", "model", "year"])[["raw_cumulative_gt", "screened_cumulative_gt"]].sum()
    yy = np.array(yrs)
    spd = pd.read_csv(S2D / "case_model_screened_paths.csv")
    spdc = spd[[c in cset[s] for s, c in zip(spd.scenario, spd.country)]]
    g4d = spdc.groupby(["scenario", "model", "year"])["screened_cumulative_gt"].sum()
    GF, GEO, ORD = [], [], []
    for s in SCEN:
        if (s, "Logistic") not in g4.index.droplevel(2) or (s, "Gompertz") not in g4.index.droplevel(2):
            continue
        Ld = g4.loc[(s, "Logistic"), "screened_cumulative_gt"].reindex(yy).to_numpy()
        Gd = g4.loc[(s, "Gompertz"), "screened_cumulative_gt"].reindex(yy).to_numpy()
        Lr = g4.loc[(s, "Logistic"), "raw_cumulative_gt"].reindex(yy).to_numpy()
        Gr = g4.loc[(s, "Gompertz"), "raw_cumulative_gt"].reindex(yy).to_numpy()
        Ldd = g4d.loc[(s, "Logistic")].reindex(yy).to_numpy() if (s, "Logistic") in g4d.index.droplevel(2) else Ld
        Gdd = g4d.loc[(s, "Gompertz")].reindex(yy).to_numpy() if (s, "Gompertz") in g4d.index.droplevel(2) else Gd
        GF.append(np.abs(Ld - Gd))                  # growth-form effect = DIFFERENCE |L - G| (never a sum)
        GEO.append(Lr - Ld); GEO.append(Gr - Gd)    # geology shortfall, each form a separate case (not mixed)
        ORD.append(np.abs(Ld - Ldd)); ORD.append(np.abs(Gd - Gdd))   # allocation order = |ascending - descending|
    GF = np.vstack(GF); GEO = np.vstack(GEO); ORD = np.vstack(ORD); flo = 0.5
    drivers = [(GF,  "#6A3D9A", "growth-form choice  (|L − G|)",   1.8, "-",             False),
               (GEO, "#C62828", "geology  (shortfall)",            1.8, "-",             True),
               (ORD, "#777777", "allocation order  (asc vs desc)", 1.2, (0, (1.4, 1.4)), False)]
    for arr, col, lab, lw_, ls_, band in drivers:
        if band:
            axC.fill_between(yy, np.clip(arr.min(0), flo, None), np.clip(arr.max(0), flo, None),
                             color=col, alpha=0.07, lw=0, zorder=2)
        axC.plot(yy, np.clip(np.median(arr, 0), flo, None), color=col, lw=lw_, ls=ls_, zorder=4, label=lab)
    axC.axvline(2050, color="#888", lw=0.7, ls=(0, (2, 2)), zorder=1)
    axC.set_yscale("log"); axC.set_xlim(2030, 2180); axC.set_ylim(flo, 3e4)
    axC.set_xlabel("Year", fontsize=6.2); axC.set_ylabel("effect on stored CO$_2$ (Gt)", fontsize=6.0)
    axC.set_xticks([2050, 2100, 2150]); axC.tick_params(labelsize=FS_MIN); axC.grid(axis="y", color="#f1f1f1", lw=0.4)
    axC.set_title("c   When each driver acts", loc="left", fontweight="bold", fontsize=7.6)
    axC.legend(loc="upper left", fontsize=FS_MIN, frameon=False, handlelength=1.4, borderaxespad=0.3)
    axC.annotate(f"geology ≈ 0\nthrough 2050\n(first bind {bind0})", xy=(2050, 1.3), xytext=(2080, 6.5), fontsize=FS_MIN,
                 color="#C62828", va="center", ha="left", arrowprops=dict(arrowstyle="->", lw=0.5, color="#C62828"))
    KEY.append(("fig4c_ord_median_max_gt", round(float(np.median(ORD, 0).max()), 2), "common-set |asc-desc|, median across cases"))
    KEY.append(("fig4c_ord_case_max_gt", round(float(ORD.max()), 1), "common-set |asc-desc|, single-case max"))
    KEY.append(("fig4c_display_floor_gt", flo, "values below are clipped to the floor for log display"))

    # ---- L/G key (top-right) ----
    fig.legend(handles=[Line2D([0], [0], color=MC_L, ls="-", lw=1.8, label="Logistic"),
                        Line2D([0], [0], color=MC_G, ls=(0, (3, 2)), lw=1.8, label="Gompertz")],
               loc="upper right", bbox_to_anchor=(0.95, 0.945), ncol=1, frameon=False,
               fontsize=6.2, handlelength=2.0)
    savefig(fig, "fig4_outcome")


# ================================================================= FIGURE 5 (11.4 x 9.9 cm, 1.5-column)
def fig5():
    mc = pd.read_csv(G3 / "global_country_scatter" / "country_mc_scatter_points.csv")
    recon = pd.read_csv(G3 / "global_country_scatter" / "country_reconstruction_points.csv")
    pit = pd.read_csv(G3 / "part3_input_table.csv")
    goal = pit.groupby(["country", "scenario"])["C_scenario_gt"].first()
    passed = pit.groupby(["country", "scenario"])["whole_period_pass"].all()
    # ascending allocation order only, the figures' stated convention (V3 averaged asc & desc)
    crec = (recon[recon.ordering == "ascending"]
            .groupby(["country", "scenario", "model"])["C_recon_gt"].first())
    _, cap = cap_by_region()                                                # 150-yr deliverable resource per jurisdiction
    FREE = ["minimum", "growth10", "reference", "maximum"]                  # the free-growth MC clouds

    REFLINES = {"minimum": "Min", "reference": "Ref", "maximum": "Max"}
    order = ["US", "China", "Canada", "Australia", "UK", "Brazil", "EU", "Thailand", "Indonesia", "Middle East"]
    fig, axes = plt.subplots(2, 5, figsize=(7.4, 4.8), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.06, right=0.99, top=0.80, bottom=0.115, hspace=0.32, wspace=0.08)

    YLO = 0.02            # below the data minimum (0.04 Gt) so no MC point is cropped (V3 floor 0.3 cropped 6%)
    # MC clouds drawn as 2-D kernel-density fields (smooth density beats grainy alpha-scatter
    # at print size) in (growth rate, log10 C) space; one soft cloud per growth form.
    gx = np.linspace(0, 30, 90)
    gly = np.linspace(np.log10(YLO), np.log10(2e4), 90)
    GX, GLY = np.meshgrid(gx, gly); GY = 10 ** GLY
    DLEV = [0.12, 0.35, 0.65, 0.9]
    for ax, c in zip(axes.flat, order):
        sub = mc[(mc.Country == c) & (mc.scenario.isin(FREE))]
        for m, col, cmap in [("Logistic", MC_L, "Blues"), ("Gompertz", MC_G, "Oranges")]:
            d = sub[sub.Model == m]
            x = d.growth_rate_pct.to_numpy(); y = d.resource_required_gt.to_numpy()
            ok = np.isfinite(x) & np.isfinite(y) & (y > 0); x, y = x[ok], y[ok]
            if len(x) < 25:                                     # too few points for a stable KDE
                ax.scatter(x, y, s=1.6, alpha=0.25, color=col, linewidths=0, rasterized=True, zorder=2)
                continue
            try:
                dens = gaussian_kde(np.vstack([x, np.log10(y)]))(np.vstack([GX.ravel(), GLY.ravel()])).reshape(GX.shape)
            except Exception:
                ax.scatter(x, y, s=1.6, alpha=0.25, color=col, linewidths=0, rasterized=True, zorder=2)
                continue
            dens /= dens.max()
            ax.contourf(GX, GY, dens, levels=DLEV + [1.0], cmap=cmap, alpha=0.42, zorder=1.5)
            ax.contour(GX, GY, dens, levels=DLEV, colors=[col], linewidths=0.35, alpha=0.75, zorder=1.6)
        # geological constraints turn the panel into a feasibility phase diagram:
        #   horizontal = 150-yr deliverable ceiling (C above it is geologically infeasible);
        #   vertical   = admissible growth cap (20%, 25% for China).
        cdel = float(cap.get(c, np.nan))
        if np.isfinite(cdel):
            ax.axhspan(cdel, 2e4, color="#C62828", alpha=0.05, lw=0, zorder=0.8)
            ax.axhline(cdel, color="#333", lw=1.0, zorder=3.6)
            # label at panel LEFT (the reconstruction arrows live at x = 19-29); flip below the
            # line when a goal line sits just above it (e.g. Australia, Brazil, EU)
            gvals = [goal.get((c, s), np.nan) for s in REFLINES]
            near_above = any(np.isfinite(g) and cdel < g < cdel * 1.8 for g in gvals)
            yl_, va_ = (cdel * 0.72, "top") if near_above else (cdel * 1.25, "bottom")
            tdel = ax.text(0.7, yl_, f"deliverable {cdel:.0f} Gt", fontsize=5.6, color="#333",
                           ha="left", va=va_, fontweight="bold", zorder=3.7)
            tdel.set_path_effects([pe.withStroke(linewidth=1.6, foreground="white")])
        gcap = 25.0 if c == "China" else 20.0
        ax.axvline(gcap, color="#888", lw=0.7, ls=(0, (2, 2)), zorder=1.2)
        # min/ref/max goal lines, labelled with the NUMBER (green = passes, red = fails)
        for s, tag in REFLINES.items():
            g = goal.get((c, s), np.nan)
            if pd.isna(g):
                continue
            lc = PASS_C if passed.get((c, s), True) else FAIL_C
            ax.axhline(g, color=lc, lw=0.8, ls=(0, (5, 3)), alpha=0.85, zorder=3)
            yg, vag = g * 1.2, "bottom"
            if np.isfinite(cdel) and g * 1.2 < cdel < g * 1.6:   # ceiling line would cross the label -> drop below the goal line
                yg, vag = g * 0.83, "top"
            tg = ax.text(0.4, yg, f"{tag} {g:.0f}", fontsize=FS_MIN, color=lc, va=vag, ha="left", fontweight="bold")
            tg.set_path_effects([pe.withStroke(linewidth=1.6, foreground="white")])
        # reconstruction shown IN the scatter: goal cap -> arrow down to reconstructed C (the shortfall)
        failed = [s for s in SCEN if (c, s) in passed.index and not passed[(c, s)]]
        if failed:
            xs = np.linspace(19, 29, len(failed)) if len(failed) > 1 else np.array([24.0])
            for x, s in zip(xs, failed):
                g = goal.get((c, s), np.nan)
                if pd.isna(g):
                    continue
                ax.plot([x - 0.9, x + 0.9], [g, g], color=TOL[s], lw=1.0, zorder=5)            # goal cap
                for m, mk in [("Logistic", "o"), ("Gompertz", "^")]:
                    cr = crec.get((c, s, m), np.nan)
                    if pd.notna(cr):
                        ax.annotate("", xy=(x, cr), xytext=(x, g),
                                    arrowprops=dict(arrowstyle="-|>", lw=0.7, color=TOL[s], alpha=0.9,
                                                    shrinkA=0, shrinkB=0), zorder=4)
                        ax.scatter([x], [cr], marker=mk, s=15, color=TOL[s], edgecolor="white", linewidth=0.3, zorder=6)
        else:
            ax.text(0.5, 0.05, "passed all\nscreening", transform=ax.transAxes, fontsize=FS_MIN, color=PASS_C,
                    ha="center", va="bottom", style="italic", linespacing=0.9)
        ax.set_yscale("log"); ax.set_ylim(YLO, 2e4); ax.set_xlim(0, 30); ax.set_xticks([0, 10, 20, 30])
        npass = int(passed.loc[c].sum()) if c in passed.index.get_level_values(0) else 0
        ntot = int(passed.loc[c].size) if c in passed.index.get_level_values(0) else 0
        ax.set_title(f"{c}  ({npass}/{ntot} pass)", fontsize=6.9, fontweight="bold", pad=2)
        ax.tick_params(labelsize=5.6); ax.grid(axis="y", color="#f2f2f2", lw=0.4)
    for ax in axes[:, 0]:
        ax.set_ylabel("Required storage C (Gt)", fontsize=6.6)
    for ax in axes[1, :]:
        ax.set_xlabel("Growth rate (%)", fontsize=6.4)
    ratio = (crec / goal.reindex(crec.index.droplevel("model")).set_axis(crec.index)).dropna()
    KEY.extend([("fig5_recon_over_goal_median_pct", round(float(ratio.median()) * 100, 1), "ascending reconstructions"),
                ("fig5_cells", int(passed.size), "country x scenario deterministic cells"),
                ("fig5_cells_pass", int(passed.sum()), "all models & orderings pass")])

    mk = [Patch(facecolor=MC_L, alpha=0.5, edgecolor="none", label="Logistic density"),
          Patch(facecolor=MC_G, alpha=0.5, edgecolor="none", label="Gompertz density"),
          Line2D([0], [0], color="#333", lw=1.0, label="deliverable ceiling (150 yr)"),
          Line2D([0], [0], color="#888", lw=0.7, ls=(0, (2, 2)), label="growth cap (20%, CN 25%)"),
          Line2D([0], [0], color=PASS_C, lw=1.0, ls=(0, (5, 3)), label="goal $C$ — passes"),
          Line2D([0], [0], color=FAIL_C, lw=1.0, ls=(0, (5, 3)), label="goal $C$ — fails"),
          Line2D([0], [0], marker=r"$\downarrow$", color="#666", lw=0, markersize=6, label="→ reconstructed $C$ (● Logistic, ▲ Gompertz)")]
    sc = [Line2D([0], [0], color=TOL[s], lw=2.4, label=SL[s]) for s in SCEN]
    leg1 = fig.legend(handles=mk, loc="upper left", bbox_to_anchor=(0.06, 0.965), ncol=5, frameon=False,
                      fontsize=FS_MIN, handlelength=1.5, columnspacing=1.0)
    fig.add_artist(leg1)
    fig.legend(handles=sc, loc="upper left", bbox_to_anchor=(0.06, 0.915), ncol=8, frameon=False,
               fontsize=FS_MIN, handlelength=1.1, columnspacing=0.9)
    savefig(fig, "fig5_commitment")


# ================================================================= TABLE 1
def table1():
    """Main-text Table 1: scenario definitions x outcomes (PNAS-friendly display item)."""
    bb = feasible_bootstrap()
    smin = bb["summary"].query("ceiling_kind == 'min'").set_index(["scenario", "model"])
    sw = pd.read_csv(S2F / "case_model_shortfall_windows.csv")
    rows = []
    for s in SCEN:
        g = sw[sw.scenario == s]
        fail = g[~g.whole_period_pass]
        rows.append({
            "Scenario": SL[s],
            "Targeted 2050 rate, Logistic (Gt/yr)": round(smin.loc[(s, "Logistic"), "med_targeted"], 2),
            "Targeted 2050 rate, Gompertz (Gt/yr)": round(smin.loc[(s, "Gompertz"), "med_targeted"], 2),
            "Feasibility reduction, L (%)": round(smin.loc[(s, "Logistic"), "reduction_pct"], 1),
            "Feasibility reduction, G (%)": round(smin.loc[(s, "Gompertz"), "reduction_pct"], 1),
            "Binding cases (of n)": f"{len(fail)} of {len(g)}",
            "Earliest binding year": int(fail.first_shortfall_year.min()) if len(fail) else "—",
            "Modeled jurisdictions (L/G)": f"{smin.loc[(s,'Logistic'),'n_countries_targeted']}/{smin.loc[(s,'Gompertz'),'n_countries_targeted']}",
        })
    t = pd.DataFrame(rows)
    t.to_csv(TAB / "table1_scenario_outcomes.csv", index=False)
    with open(TAB / "table1_scenario_outcomes.md", "w") as f:
        f.write("**Table 1 | Scenario demand levels, whole-period geological feasibility, and binding onset.**\n"
                "Targeted rates are medians of the global 2050 storage-rate bootstrap (10 modeled jurisdictions; "
                "IPCC scenarios cover 7–8). 'Feasibility reduction' removes sampled demand futures whose lifetime "
                "commitment exceeds the deliverable resource (2050 targets held fixed); zeros are exact by "
                "common-random-number construction. Binding cases/years from the year-by-year delivery test "
                "(ascending allocation order). The IPCC High Gompertz reduction is one scenario-infeasible cell "
                "(Canada) whose target admits no feasible pathway.\n\n")
        cols = list(t.columns)
        f.write("| " + " | ".join(cols) + " |\n")
        f.write("|" + "|".join(["---"] * len(cols)) + "|\n")
        for _, r in t.iterrows():
            f.write("| " + " | ".join(str(r[c]) for c in cols) + " |\n")
    print("  saved table1_scenario_outcomes.{csv,md}")


if __name__ == "__main__":
    for n, fn in [("F1", fig1), ("F2", fig2), ("F3", fig3), ("F4", fig4), ("F5", fig5), ("T1", table1)]:
        print(f"[{n}]")
        try:
            fn()
        except Exception:
            print(f"  !! {n} FAILED"); traceback.print_exc()
    pd.DataFrame(KEY, columns=["name", "value", "source"]).to_csv(TAB / "fig_key_numbers.csv", index=False)
    print("  wrote tables: fig_key_numbers.csv (+ fig4 bootstrap tables, table1)")
    print("done.")
