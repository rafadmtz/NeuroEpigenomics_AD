#!/usr/bin/env python3
"""
09_te_analysis.py — TE family enrichment among differentially accessible peaks
==============================================================================
 
Tests whether some TE families are over- or underrepresented among the
differentially accessible (DA) peaks of each subclass, taking into account how
common each family is among all the tested peaks that overlap a TE.
 
Two levels are analysed:
  1. Families   : every repFamily of LINE, LTR and SINE, compared within its
                  repClass, in all subclasses.
  2. ERV1       : every repName of LTR/ERV1 (subfamilies), compared within the
                  whole subclass, in L2-3-IT and L4-IT only.
 
Definitions:
  - DA peak: padj < 0.05 and |log2FoldChange| > 0.1 (is_up if positive,
    is_down if negative). These flags are recomputed here because the ones
    of step 07 only use padj.
  - Background: all tested peaks of the subclass (and of the repClass at the
    family level) that overlap a TE. Each peak counts once per subclass and
    TE group, so a peak with two LTR families counts in both.
  - log2_oe = log2(fraction of DA peaks in the group / fraction of tested
    peaks in the group). Positive means over-represented among DA peaks,
    negative means under-represented.
  - p-value from Fisher's exact test; fdr is Benjamini-Hochberg over all the
    tests of a run.
  - Groups with fewer than MIN_DA DA peaks are dropped.
 
Inputs:
  - results/05_te_annotation/peak_te_overlaps.csv
  - results/07_build_peak_table/peak_table_{subclass}_trend.parquet   (step 07)
 
Outputs:
  - results/09_te_analysis/family_enrichment.csv
  - results/09_te_analysis/erv1_repname_enrichment.csv
  - results/figures/09_te_analysis/family_bubble_{up,down}.png
  - results/figures/09_te_analysis/erv1_repname_up.png
  The output folders are created if they do not exist.
 
Usage:
  Run from the repository root (environment: see environment.yml):
    python src/09_te_analysis.py
"""
 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from pathlib import Path
 
PROJECT_ROOT = Path(__file__).resolve().parents[1]
TE_DIR = PROJECT_ROOT / "results" / "05_te_annotation"
PEAK_DIR = PROJECT_ROOT / "results" / "07_build_peak_table"
FIG_DIR = PROJECT_ROOT / "results" / "figures" / "09_te_analysis"
OUT_DIR = PROJECT_ROOT / "results" / "09_te_analysis"
 
TES = ["LINE", "LTR", "SINE"]
SUBCLASSES = ["L2-3-IT", "L4-IT", "L5-IT", "Vip", "Sst", "Astrocyte"]
DIRECTIONS = ["is_up", "is_down"]
MIN_DA = 5
 
 
def get_peak_overlap_df(dedup_col: str = "repFamily") -> pd.DataFrame:
    """
    Build the table of peaks that overlap a TE, with their differential accessibility.
 
    Redefines `is_up` and `is_down` as padj < 0.05 and log2FoldChange > 0.1
    (or < -0.1). Peaks without a TE overlap are dropped.
 
    Parameters
    ----------
    dedup_col : str
        Column used to count each peak once per subclass: "repFamily" or
        "repName". Default "repFamily".
 
    Returns
    -------
    pd.DataFrame
        One row per peak, subclass and `dedup_col`: peak statistics
        (baseMean, log2FoldChange, padj, is_up, is_down, gc) and TE annotation
        (repClass, repFamily, repName, milliDiv).
 
    Raises
    ------
    ValueError
        If `dedup_col` is not valid or no peak overlaps a TE.
    FileNotFoundError
        If the TE overlaps or a peak table is missing.
    KeyError
        If an input table lacks a required column.
    """
    if dedup_col not in ("repFamily", "repName"):
        raise ValueError(f"dedup_col must be 'repFamily' or 'repName', got '{dedup_col}'")
 
    try:
        peak_overlap = pd.read_csv(TE_DIR / "peak_te_overlaps.csv")
        peak_overlap = peak_overlap[peak_overlap["repClass"].isin(TES)][
            ["peak_id", "repClass", "repFamily", "repName", "milliDiv"]
        ]
    except FileNotFoundError as e:
        raise FileNotFoundError(f"TE overlaps not found in {TE_DIR}") from e
    except KeyError as e:
        raise KeyError(f"Missing column in peak_te_overlaps.csv: {e}") from e
 
    peak_tables = []
    for subclass in SUBCLASSES:
        try:
            peak_tables.append(pd.read_parquet(PEAK_DIR / f"peak_table_{subclass}_trend.parquet"))
        except FileNotFoundError as e:
            raise FileNotFoundError(f"Peak table of {subclass} not found in {PEAK_DIR}") from e
    da_peaks = pd.concat(peak_tables)
 
    try:
        da_peaks = da_peaks[[
            "peak_id",
            "subclass",
            "baseMean",
            "log2FoldChange",
            "padj",
            "is_up",
            "is_down",
            "is_ns",
            "gc"
        ]]
    except KeyError as e:
        raise KeyError(f"Missing column in the peak tables: {e}") from e
 
    da_peaks["is_up"] = (da_peaks["padj"] < 0.05) & (da_peaks["log2FoldChange"] > 0.1)
    da_peaks["is_down"] = (da_peaks["padj"] < 0.05) & (da_peaks["log2FoldChange"] < -0.1)
 
    df_peak_annotation = da_peaks.merge(peak_overlap, on="peak_id", how="inner")
    if df_peak_annotation.empty:
        raise ValueError("No peak overlaps a TE. Check that peak_id has the same format in both tables.")
 
    df_peak_annotation = df_peak_annotation.drop_duplicates(
        ["peak_id", "subclass", dedup_col]
    )
    return df_peak_annotation
 
 
def compute_hypergeometric(
    df_peak_annotation: pd.DataFrame,
    te_level: str = "repFamily",
    background_class: str | None = "repClass",
) -> pd.DataFrame:
    """
    Compare how often each TE group is DA against how often it is tested.
 
    For each subclass (and `background_class`), `log2_oe` is the log2 of the
    fraction of DA peaks in the TE group over the fraction of all tested peaks
    in that group. Positive means over-represented among DA peaks. The p-value
    comes from Fisher's exact test. Groups with fewer than `MIN_DA` DA peaks
    are dropped.
 
    Parameters
    ----------
    df_peak_annotation : pd.DataFrame
        Output of get_peak_overlap_df, one row per peak, subclass and TE group.
    te_level : str
        Column with the TE groups to compare: "repFamily" or "repName".
        Default "repFamily".
    background_class : str or None
        Column that defines the background together with `subclass`
        (e.g. "repClass"). If None, the background is the whole subclass.
        Default "repClass".
 
    Returns
    -------
    pd.DataFrame
        One row per subclass, background, TE group and direction: subclass,
        `background_class`, `te_level`, direction, n_da, log2_oe, odds_ratio,
        p, fdr (Benjamini-Hochberg).
 
    Raises
    ------
    KeyError
        If a required column is missing in `df_peak_annotation`.
    ValueError
        If no TE group reaches `MIN_DA` DA peaks.
    """
    background_cols = ["subclass"] + ([background_class] if background_class else [])
 
    missing = [c for c in background_cols + [te_level] + DIRECTIONS if c not in df_peak_annotation.columns]
    if missing:
        raise KeyError(f"Missing columns in df_peak_annotation: {missing}")
 
    rows = []
    for background_key, g_background in df_peak_annotation.groupby(background_cols):
        background_key = background_key if isinstance(background_key, tuple) else (background_key,)
        n_background = len(g_background)
        for direction in DIRECTIONS:
            n_da_background = g_background[direction].sum()
            for te_name, g_te in g_background.groupby(te_level):
 
                n_da_te = g_te[direction].sum()
                if n_da_te < MIN_DA:
                    continue
                n_non_da_te = len(g_te) - n_da_te
                n_da_outside = n_da_background - n_da_te
                n_non_da_outside = (n_background - n_da_background) - n_non_da_te
 
                obs = n_da_te / n_da_background
                exp = len(g_te) / n_background
 
                odds_ratio, p_value = fisher_exact([[n_da_te, n_non_da_te],
                                                    [n_da_outside, n_non_da_outside]])
 
                rows.append(dict(
                    zip(background_cols, background_key),
                    **{te_level: te_name},
                    direction=direction,
                    n_da=n_da_te,
                    log2_oe=np.log2(obs / exp),
                    odds_ratio=odds_ratio,
                    p=p_value,
                ))
 
    if not rows:
        raise ValueError(f"No {te_level} has at least {MIN_DA} DA peaks. Lower MIN_DA or check the input.")
 
    result = pd.DataFrame(rows)
    result["fdr"] = multipletests(result["p"], method="fdr_bh")[1]
    return result
 
 
def bubble(result: pd.DataFrame, direction: str, out_path: Path) -> None:
    """
    Plot TE family enrichment as bubbles, one panel per TE class.
 
    Color is `log2_oe`, size is `n_da` and a black edge marks fdr < 0.05.
    The figure is saved to `out_path` and closed.
 
    Parameters
    ----------
    result : pd.DataFrame
        Output of compute_hypergeometric at repFamily level with repClass.
    direction : str
        "is_up" or "is_down".
    out_path : Path
        File where the figure is saved (PNG).
 
    Raises
    ------
    ValueError
        If `direction` is not valid or `result` has no rows for it.
    """
    if direction not in DIRECTIONS:
        raise ValueError(f"direction must be one of {DIRECTIONS}, got '{direction}'")
 
    r = result[result["direction"] == direction]
    if r.empty:
        raise ValueError(f"No results for direction '{direction}'")
 
    vmax = r["log2_oe"].abs().max()
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), layout="constrained")
    for ax, cl in zip(axes, TES):
        d = r[r["repClass"] == cl]
        fams = sorted(d["repFamily"].unique())
        x = d["subclass"].map({s: i for i, s in enumerate(SUBCLASSES)})
        y = d["repFamily"].map({f: i for i, f in enumerate(fams)})
        sc = ax.scatter(x, y, s=30 + 25 * np.sqrt(d["n_da"]), c=d["log2_oe"],
                        cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                        edgecolors=np.where(d["fdr"] < 0.05, "black", "none"), linewidths=1.5)
        ax.set_xticks(range(len(SUBCLASSES)))
        ax.set_xticklabels(SUBCLASSES, rotation=45, ha="right")
        ax.set_yticks(range(len(fams))); ax.set_yticklabels(fams)
        ax.set_xlim(-0.5, len(SUBCLASSES) - 0.5); ax.set_ylim(-0.5, len(fams) - 0.5)
        ax.set_title(cl)
    fig.colorbar(sc, ax=axes, label="log2(obs/exp)", shrink=0.8)
    fig.suptitle(f"TE families among {direction} peaks (padj < 0.05, |log2FC| > 0.1)")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
 
 
def plot_names(result_name: pd.DataFrame, out_path: Path, direction: str = "is_up") -> None:
    """
    Plot ERV1 subfamily enrichment (`log2_oe`) in L2-3-IT and L4-IT.
 
    Size is `n_da`, color is `log2_oe` and a black edge marks fdr < 0.05.
    The figure is saved to `out_path` and closed.
 
    Parameters
    ----------
    result_name : pd.DataFrame
        Output of compute_hypergeometric at repName level.
    out_path : Path
        File where the figure is saved (PNG).
    direction : str
        "is_up" or "is_down". Default "is_up".
 
    Raises
    ------
    ValueError
        If `direction` is not valid or `result_name` has no rows for it.
    """
    if direction not in DIRECTIONS:
        raise ValueError(f"direction must be one of {DIRECTIONS}, got '{direction}'")
 
    r = result_name[result_name["direction"] == direction]
    if r.empty:
        raise ValueError(f"No results for direction '{direction}'")
 
    names = sorted(r["repName"].unique())
    vmax = r["log2_oe"].abs().max()
    fig, axes = plt.subplots(1, 2, figsize=(10, max(4, 0.35 * len(names))), sharey=True)
    for ax, sc in zip(axes, ["L2-3-IT", "L4-IT"]):
        d = r[r["subclass"] == sc]
        ax.scatter(d["log2_oe"], d["repName"].map({n: i for i, n in enumerate(names)}),
                   s=30 + 25 * np.sqrt(d["n_da"]), c=d["log2_oe"], cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax,
                   edgecolors=np.where(d["fdr"] < 0.05, "black", "none"), linewidths=1.5)
        ax.axvline(0, color="grey", lw=0.8)
        ax.set_yticks(range(len(names))); ax.set_yticklabels(names)
        ax.set_xlabel("log2(obs/exp)"); ax.set_title(sc)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
 
 
def main() -> None:
    """Run the family-level enrichment and the ERV1 repName-level enrichment."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
 
    df_fam = get_peak_overlap_df(dedup_col="repFamily")
    result_fam = compute_hypergeometric(df_fam)
    result_fam.to_csv(OUT_DIR / "family_enrichment.csv", index=False)
    bubble(result_fam, "is_up", FIG_DIR / "family_bubble_up.png")
    bubble(result_fam, "is_down", FIG_DIR / "family_bubble_down.png")
 
    try:
        df_name = get_peak_overlap_df(dedup_col="repName")
        erv1 = df_name[
            (df_name["repFamily"] == "LTR/ERV1")
            & df_name["subclass"].isin(["L2-3-IT", "L4-IT"])
        ]
        if erv1.empty:
            raise ValueError("No LTR/ERV1 peaks found. Check the values of repFamily.")
        result_name = compute_hypergeometric(erv1, te_level="repName", background_class=None)
        result_name.to_csv(OUT_DIR / "erv1_repname_enrichment.csv", index=False)
        plot_names(result_name, FIG_DIR / "erv1_repname_up.png", "is_up")
    except ValueError as e:
        print(f"ERV1 analysis skipped: {e}")
 
 
if __name__ == "__main__":
    main()