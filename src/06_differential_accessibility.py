#!/usr/bin/env python3
"""
06_differential_accessibility.py — Differential accessibility with PyDESeq2
==========================================================================

Runs a per-subclass differential accessibility (DA) analysis on the per-donor
pseudobulk counts produced in step 03, using PyDESeq2 (negative binomial GLM,
donors as biological replicates).

Two questions are asked for every subclass:
  1. trend          : accessibility change per step of AD neuropathological
                      change (Not AD=0, Low=1, Intermediate=2, High=3), with
                      the stage as a continuous variable. Uses all donors in a
                      single coefficient.
                      Design: ~ sex + age + pmi + stage
  2. High_vs_NotAD  : the extreme contrast, with the condition as a factor and
                      "Not AD" as reference (all four levels are kept in the
                      model so every donor informs the dispersion estimates).
                      Design: ~ sex + age + pmi + condition

Age at death and PMI are z-scored to help the GLM converge.

Filters:
  - A subclass is analysed only if every condition has >= MIN_DONORS_PER_CONDITION
    donors (the donor-level filter of >= 10 nuclei was already applied in step 03).
  - A peak is tested only if it has >= MIN_COUNT counts in at least as many donors
    as the smallest condition group of that subclass.

Inputs:
  - data/processed/pseudobulk/pseudobulk_{subclass}.h5ad   (step 03)

Outputs (in results/06_differential_accessibility/):
  - da_{subclass}_trend.csv.gz           One row per tested peak.
  - da_{subclass}_High_vs_NotAD.csv.gz   One row per tested peak.
        Columns: peak_id, chr, start, end, baseMean, log2FoldChange, lfcSE,
        stat, pvalue, padj, log2FC_shrunk (apeGLM shrinkage).
  - vst_{subclass}.h5ad                  Donors x tested peaks, VST values
                                         (blind to the design), size factors in obs.
  - da_summary.csv                       One row per subclass x contrast: donors
                                         per condition, peaks tested, significant
                                         peaks up/down, status.
  - 06_da_log.txt
Figures (in results/figures/06_differential_accessibility/):
  - pca_{subclass}.png, volcano_{subclass}.png, ma_{subclass}.png,
    dispersion_{subclass}.png, n_significant_per_subclass.png

Usage:
  Run from the repository root (environment: neuroAD, see environment.yml):
    python src/06_differential_accessibility.py

  Run in the background and save output to a log:
    nohup python -u src/06_differential_accessibility.py > logs/06_differential_accessibility.log 2>&1 &
"""

import importlib.metadata
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

import matplotlib
matplotlib.use("Agg")  # headless / server-safe backend
import matplotlib.pyplot as plt

#==================================================
# Variable definitions
#==================================================
PROJECT_ROOT = Path(__file__).resolve().parents[1]
PSEUDOBULK_DIR = PROJECT_ROOT / "data" / "processed" / "pseudobulk"
OUT_DIR = PROJECT_ROOT / "results" / "06_differential_accessibility"
FIG_DIR = PROJECT_ROOT / "results" / "figures" / "06_differential_accessibility"
SUMMARY_CSV = OUT_DIR / "da_summary.csv"
LOG_PATH = OUT_DIR / "06_da_log.txt"

# Metadata columns written by step 03.
COND_COL = "Overall AD neuropathological Change"
SEX_COL = "Sex"
AGE_COL = "Age at Death"
PMI_COL = "PMI"

# Ordered stages; the first level is the reference of the factor model.
GROUPS = ["Not AD", "Low", "Intermediate", "High"]

ALPHA = 0.05                  # padj threshold for counting significant peaks
MIN_DONORS_PER_CONDITION = 3  # every condition needs at least this many donors
MIN_COUNT = 10                # peak filter: counts >= MIN_COUNT ...
                              # ... in >= (smallest condition group) donors
N_CPUS = 16
N_TOP_VAR_PCA = 5000          # most variable VST peaks used for the PCA

CONTRASTS = ["trend", "High_vs_NotAD"]
PEAK_RE = re.compile(r"^(chr[\w]+)[:\-](\d+)[:\-](\d+)$")

#==================================================
# Functions
#==================================================
def log_print(msg):
    """Prints a message and appends it to the log file."""
    print(msg)
    with open(LOG_PATH, "a") as f:
        f.write(msg + "\n")


def subclass_name(path: Path) -> str:
    """Recover the subclass label from a pseudobulk file."""
    obs = ad.read_h5ad(path, backed="r").obs
    return str(obs["subclass"].iloc[0])


def load_pseudobulk(path: Path):
    """
    Load one pseudobulk file into a counts DataFrame (donors x peaks, int) and a
    design metadata DataFrame (sex, age, pmi, condition, stage).
    """
    pb = ad.read_h5ad(path)
    donors = pb.obs_names.astype(str)
    counts = pd.DataFrame(
        np.rint(pb.X.toarray()).astype(np.int64),
        index=donors, columns=pb.var_names.astype(str),
    )

    obs = pb.obs
    # Numeric covariates may be stored as categoricals in the .h5ad.
    age = pd.to_numeric(obs[AGE_COL].astype(str)).to_numpy()
    pmi = pd.to_numeric(obs[PMI_COL].astype(str)).to_numpy()
    meta = pd.DataFrame({
        "sex": obs[SEX_COL].astype(str).to_numpy(),
        "age": (age - age.mean()) / age.std(),
        "pmi": (pmi - pmi.mean()) / pmi.std(),
        "condition": pd.Categorical(obs[COND_COL].astype(str), categories=GROUPS),
    }, index=donors)
    meta["stage"] = meta["condition"].cat.codes.astype(float)
    return counts, meta


def donors_per_condition(meta: pd.DataFrame) -> pd.Series:
    """Number of donors in each condition, in GROUPS order."""
    return meta["condition"].value_counts().reindex(GROUPS, fill_value=0)


def filter_peaks(counts: pd.DataFrame, n_min: int) -> pd.DataFrame:
    """Keep peaks with >= MIN_COUNT counts in at least n_min donors."""
    keep = (counts >= MIN_COUNT).sum(axis=0) >= n_min
    return counts.loc[:, keep]


def build_design(meta: pd.DataFrame, term: str) -> str:
    """Design formula; drops covariates without variation in this subclass."""
    covs = [c for c in ("sex", "age", "pmi") if meta[c].nunique() > 1]
    return "~" + " + ".join(covs + [term])


def run_contrast(counts, meta, contrast_name):
    """
    Fit one PyDESeq2 model and return (results DataFrame, fitted dds).
    Results keep the MLE log2FC and add the apeGLM-shrunk log2FC.
    """
    if contrast_name == "trend":
        design = build_design(meta, "stage")
    else:
        design = build_design(meta, "condition")

    dds = DeseqDataSet(counts=counts, metadata=meta, design=design,
                       n_cpus=N_CPUS, quiet=True)
    dds.deseq2()

    columns = list(dds.obsm["design_matrix"].columns)
    if contrast_name == "trend":
        contrast = np.array([1.0 if c == "stage" else 0.0 for c in columns])
        coeff = "stage"
    else:
        contrast = ["condition", "High", "Not AD"]
        coeff = "condition[T.High]"

    stats = DeseqStats(dds, contrast=contrast, alpha=ALPHA,
                       n_cpus=N_CPUS, quiet=True)
    stats.summary()
    res = stats.results_df.copy()
    stats.lfc_shrink(coeff=coeff)
    res["log2FC_shrunk"] = stats.results_df["log2FoldChange"]
    return res, dds


def format_results(res: pd.DataFrame) -> pd.DataFrame:
    """Add peak coordinates and order the columns."""
    out = res.rename_axis("peak_id").reset_index()
    coords = out["peak_id"].str.extract(PEAK_RE)
    out.insert(1, "chr", coords[0])
    out.insert(2, "start", pd.to_numeric(coords[1]))
    out.insert(3, "end", pd.to_numeric(coords[2]))
    return out[["peak_id", "chr", "start", "end", "baseMean", "log2FoldChange",
                "lfcSE", "stat", "pvalue", "padj", "log2FC_shrunk"]]


def write_vst(dds, meta: pd.DataFrame, safe: str) -> np.ndarray:
    """Blind VST of the tested peaks, saved as donors x peaks .h5ad."""
    dds.vst(use_design=False)
    vst = np.asarray(dds.layers["vst_counts"], dtype=np.float32)
    obs = meta.copy()
    obs["condition"] = obs["condition"].astype(str)
    obs["size_factors"] = dds.obs["size_factors"].to_numpy()
    ad.AnnData(X=vst, obs=obs, var=pd.DataFrame(index=dds.var_names)).write_h5ad(
        OUT_DIR / f"vst_{safe}.h5ad", compression="gzip"
    )
    return vst


#==================================================
# Figures
#==================================================
COND_COLORS = dict(zip(GROUPS, ["#4C78A8", "#F2B701", "#E45756", "#72241B"]))


def plot_pca(vst: np.ndarray, meta: pd.DataFrame, subclass: str, safe: str):
    """PCA of donors on the most variable VST peaks, by condition and by sex."""
    top = np.argsort(vst.var(axis=0))[::-1][:N_TOP_VAR_PCA]
    mat = vst[:, top] - vst[:, top].mean(axis=0)
    u, s, _ = np.linalg.svd(mat, full_matrices=False)
    pcs = u[:, :2] * s[:2]
    var_exp = 100 * s[:2] ** 2 / np.sum(s ** 2)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for axis, col in zip(axes, ["condition", "sex"]):
        labels = meta[col].astype(str).to_numpy()
        levels = GROUPS if col == "condition" else sorted(set(labels))
        for lev in levels:
            m = labels == lev
            if m.any():
                axis.scatter(pcs[m, 0], pcs[m, 1], s=28, label=lev,
                             color=COND_COLORS.get(lev), alpha=0.85)
        axis.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)")
        axis.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)")
        axis.set_title(f"{subclass}: donors by {col}")
        axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"pca_{safe}.png", dpi=150)
    plt.close(fig)


def plot_volcano(results: dict, subclass: str, safe: str):
    """Volcano (shrunk log2FC vs -log10 padj) for both contrasts."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for axis, name in zip(axes, CONTRASTS):
        res = results[name].dropna(subset=["padj"])
        y = -np.log10(res["padj"].clip(lower=1e-300))
        up = (res["padj"] < ALPHA) & (res["log2FC_shrunk"] > 0)
        down = (res["padj"] < ALPHA) & (res["log2FC_shrunk"] < 0)
        ns = ~(up | down)
        axis.scatter(res.loc[ns, "log2FC_shrunk"], y[ns], s=3, color="lightgrey")
        axis.scatter(res.loc[up, "log2FC_shrunk"], y[up], s=4, color="#C0392B",
                     label=f"up ({up.sum()})")
        axis.scatter(res.loc[down, "log2FC_shrunk"], y[down], s=4, color="#2471A3",
                     label=f"down ({down.sum()})")
        axis.axhline(-np.log10(ALPHA), color="black", lw=0.6, ls="--")
        xlabel = "log2FC per stage (shrunk)" if name == "trend" else "log2FC High vs Not AD (shrunk)"
        axis.set_xlabel(xlabel)
        axis.set_ylabel("-log10(padj)")
        axis.set_title(f"{subclass}: {name}")
        axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"volcano_{safe}.png", dpi=150)
    plt.close(fig)


def plot_ma(results: dict, subclass: str, safe: str):
    """MA plot (baseMean vs shrunk log2FC) for both contrasts."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for axis, name in zip(axes, CONTRASTS):
        res = results[name]
        sig = (res["padj"] < ALPHA).to_numpy()
        axis.scatter(res.loc[~sig, "baseMean"], res.loc[~sig, "log2FC_shrunk"],
                     s=3, color="lightgrey")
        axis.scatter(res.loc[sig, "baseMean"], res.loc[sig, "log2FC_shrunk"],
                     s=4, color="#C0392B", label=f"padj < {ALPHA} ({sig.sum()})")
        axis.axhline(0, color="black", lw=0.6)
        axis.set_xscale("log")
        axis.set_xlabel("baseMean (normalized counts)")
        axis.set_ylabel("log2FC (shrunk)")
        axis.set_title(f"{subclass}: {name}")
        axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"ma_{safe}.png", dpi=150)
    plt.close(fig)


def plot_dispersion(dds, subclass: str, safe: str):
    """Gene-wise, fitted and final dispersions against mean accessibility."""
    var = dds.var
    order = np.argsort(var["_normed_means"].to_numpy())
    fig, axis = plt.subplots(figsize=(6, 4.5))
    axis.scatter(var["_normed_means"], var["genewise_dispersions"], s=2,
                 color="lightgrey", label="gene-wise")
    axis.scatter(var["_normed_means"], var["dispersions"], s=2,
                 color="#2471A3", alpha=0.5, label="final (MAP)")
    axis.plot(var["_normed_means"].to_numpy()[order],
              var["fitted_dispersions"].to_numpy()[order],
              color="#C0392B", lw=1.2, label="fitted trend")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("Mean of normalized counts")
    axis.set_ylabel("Dispersion")
    axis.set_title(f"{subclass}: dispersion (factor model)")
    axis.legend(fontsize=8, markerscale=4)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"dispersion_{safe}.png", dpi=150)
    plt.close(fig)


def plot_n_significant(summary: pd.DataFrame):
    """Significant peaks up/down per subclass, one panel per contrast."""
    ok = summary[summary["status"] == "ok"]
    if ok.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)
    for axis, name in zip(axes, CONTRASTS):
        sub = ok[ok["contrast"] == name].sort_values("subclass")
        x = np.arange(sub.shape[0])
        axis.bar(x, sub["n_sig_up"], color="#C0392B", label="up")
        axis.bar(x, -sub["n_sig_down"], color="#2471A3", label="down")
        axis.axhline(0, color="black", lw=0.6)
        axis.set_xticks(x)
        axis.set_xticklabels(sub["subclass"], rotation=90, fontsize=7)
        axis.set_ylabel(f"Significant peaks (padj < {ALPHA})")
        axis.set_title(name)
        axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "n_significant_per_subclass.png", dpi=150)
    plt.close(fig)


#==================================================
# Per-subclass pipeline
#==================================================
def analyse_subclass(path: Path) -> list:
    """Run both contrasts for one subclass; return summary rows."""
    subclass = subclass_name(path)
    safe = path.stem.replace("pseudobulk_", "")
    counts, meta = load_pseudobulk(path)
    n_cond = donors_per_condition(meta)
    base = {"subclass": subclass,
            **{f"n_{g.replace(' ', '')}": int(n_cond[g]) for g in GROUPS},
            "n_peaks_total": counts.shape[1]}

    if (n_cond < MIN_DONORS_PER_CONDITION).any():
        reason = (f"< {MIN_DONORS_PER_CONDITION} donors in: "
                  + ", ".join(n_cond[n_cond < MIN_DONORS_PER_CONDITION].index))
        log_print(f"[skip] {subclass}: {reason}")
        return [{**base, "contrast": c, "status": "skipped", "reason": reason}
                for c in CONTRASTS]

    n_min = int(n_cond.min())
    counts = filter_peaks(counts, n_min)
    log_print(f"\n[{subclass}] donors per condition {n_cond.to_dict()} | "
              f"peaks tested: {counts.shape[1]} / {base['n_peaks_total']} "
              f"(>= {MIN_COUNT} counts in >= {n_min} donors)")

    results, rows, dds_factor = {}, [], None
    for name in CONTRASTS:
        res, dds = run_contrast(counts, meta, name)
        res = format_results(res)
        res.to_csv(OUT_DIR / f"da_{safe}_{name}.csv.gz", index=False)
        results[name] = res
        if name == "High_vs_NotAD":
            dds_factor = dds

        sig = res["padj"] < ALPHA
        n_up = int((sig & (res["log2FoldChange"] > 0)).sum())
        n_down = int((sig & (res["log2FoldChange"] < 0)).sum())
        log_print(f"  {name}: {n_up} up, {n_down} down (padj < {ALPHA})")
        rows.append({**base, "contrast": name, "status": "ok", "reason": "",
                     "n_peaks_tested": counts.shape[1],
                     "n_peaks_filtered": base["n_peaks_total"] - counts.shape[1],
                     "n_sig_up": n_up, "n_sig_down": n_down})

    vst = write_vst(dds_factor, meta, safe)
    plot_pca(vst, meta, subclass, safe)
    plot_volcano(results, subclass, safe)
    plot_ma(results, subclass, safe)
    plot_dispersion(dds_factor, subclass, safe)
    return rows


#==================================================
# Main
#==================================================
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    with open(LOG_PATH, "w") as f:
        f.write("=== DIFFERENTIAL ACCESSIBILITY (PyDESeq2) START ===\n")

    versions = {p: importlib.metadata.version(p)
                for p in ("pydeseq2", "anndata", "numpy", "pandas")}
    log_print(f"Versions: {versions}")
    log_print(f"Parameters: ALPHA={ALPHA}, MIN_DONORS_PER_CONDITION="
              f"{MIN_DONORS_PER_CONDITION}, MIN_COUNT={MIN_COUNT}, N_CPUS={N_CPUS}")

    files = sorted(PSEUDOBULK_DIR.glob("pseudobulk_*.h5ad"))
    if not files:
        raise FileNotFoundError(f"No pseudobulk files in {PSEUDOBULK_DIR}; run step 03.")
    log_print(f"Pseudobulk files found: {len(files)}")

    rows = []
    for path in files:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=UserWarning)
                rows.extend(analyse_subclass(path))
        except Exception as e:  # keep going with the other subclasses
            log_print(f"[error] {path.name}: {e!r}")
            rows.extend({"subclass": path.stem.replace("pseudobulk_", ""),
                         "contrast": c, "status": "error", "reason": repr(e)}
                        for c in CONTRASTS)
        # Write the summary after every subclass so partial runs are usable.
        pd.DataFrame(rows).to_csv(SUMMARY_CSV, index=False)

    summary = pd.DataFrame(rows)
    plot_n_significant(summary)
    log_print(f"\n-> Summary saved to: {SUMMARY_CSV.name}")
    log_print("=== DIFFERENTIAL ACCESSIBILITY FINISHED ===")


if __name__ == "__main__":
    main()
