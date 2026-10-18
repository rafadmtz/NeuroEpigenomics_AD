
"""
03_pseudobulk_by_donor.py
==========================================================================

Build per-donor pseudobulk count files (one per cell subclass) to be used as
input for differential accessibility (DA) analysis (DESeq2 / edgeR).

Rationale
---------
For a multi-donor, multi-condition design, the biological unit of replication
is the *donor*, not the individual nucleus. Testing at the nucleus level is
pseudoreplication and inflates false positives (Squair et al. 2021, Nat Commun).
The recommended strategy is therefore to aggregate (sum) the peak counts of all
nuclei belonging to the same donor within a given subclass -> one "bulk-like"
count vector per donor per subclass. The condition (AD neuropathological stage)
enters the downstream model as a donor-level design variable.

Memory strategy
---------------
The input .h5ad is ~18 GB and a single subclass can hold hundreds of thousands
of nuclei, so it is NEVER loaded fully into memory. Instead the count matrix is
read in row-chunks straight from the backed file and summed on-the-fly into a
dense accumulator of shape (groups x peaks), where a group is a
(subclass, donor) pair. Memory is bounded by:
    accumulator  ~= groups_in_pass * n_peaks * 4 bytes   (float32)
    one chunk    ~= nnz(chunk)
Subclasses are processed in batches sized so the accumulator stays under
MEM_BUDGET_GB; each batch is one streaming pass over the file.

What this script does
---------------------
  1. Reads nucleus metadata (obs), keeps the selected conditions, and assigns
     each nucleus to a (subclass, donor) group.
  2. Streams the backed count matrix in row-chunks, accumulating per-group sums.
  3. Builds a per-donor sample table (colData): condition, covariates and the
     number of aggregated nuclei; drops donors with too few nuclei.
  4. Writes one pseudobulk .h5ad per subclass (obs = donors, var = peaks,
     X = summed integer counts) plus a colData CSV, ready for DA in R.
  5. Draws QC plots of the counts (nuclei per donor x subclass, pseudobulk
     library sizes per condition, donors passing the filter per condition).

Use:
  Run from the project root:
    python src/03_pseudobulk_by_donor.py

  Run in the background and save output to a log:
    nohup python -u src/03_pseudobulk_by_donor.py > logs/03_pseudobulk_by_donor.log 2>&1 &
"""

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
import scanpy as sc
import anndata as ad
from anndata._io.specs import read_elem

import matplotlib
matplotlib.use("Agg")  # headless / server-safe backend
import matplotlib.pyplot as plt

# ==================================================
# Configuration
# ==================================================
PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_FILE = PROJECT_ROOT / "data" / "SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad"
OUTPUT_DIR = PROJECT_ROOT / "data" / "processed" / "pseudobulk"
FIG_DIR = PROJECT_ROOT / "results" / "figures" / "pseudobulk"

# Metadata (obs) column names. Adjust here if the dataset uses different labels.
COND_COL = "Overall AD neuropathological Change"
SUBCLASS_COL = "Subclass"
DONOR_COL = "Donor ID"
# Donor-level covariates to carry into colData (only those present are kept).
COVARIATES = ["Sex", "Age at Death", "PMI"]

# Conditions to keep (the four AD neuropathological-change stages).
GROUPS = ["Not AD", "Low", "Intermediate", "High"]

# Minimum nuclei for a (donor, subclass) pseudobulk sample to be retained.
MIN_NUCLEI = 10

# --- Memory controls ---
# Nuclei read per streaming step. Smaller -> less peak memory, more I/O.
CHUNK_SIZE = 20_000
# Upper bound for the dense per-pass accumulator (groups x peaks, float32).
# Subclasses are batched so this is never exceeded. ~8 GB is safe on 64 GB RAM.
MEM_BUDGET_GB = 8.0


# ==================================================
# IO helpers
# ==================================================
def load_metadata(path: Path) -> pd.DataFrame:
    """Read only the nucleus metadata (obs) from a local .h5ad file."""
    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"])
    return obs


def check_columns(obs: pd.DataFrame) -> list:
    """Validate required columns and return the covariates actually present."""
    for col in (COND_COL, SUBCLASS_COL, DONOR_COL):
        if col not in obs.columns:
            raise KeyError(
                f"Required column '{col}' not found in obs. "
                f"Available columns include: {list(obs.columns)[:20]} ..."
            )
    present = [c for c in COVARIATES if c in obs.columns]
    missing = [c for c in COVARIATES if c not in obs.columns]
    if missing:
        print(f"[warn] Covariates not found and skipped: {missing}")
    return present


# ==================================================
# Grouping
# ==================================================
def build_groups(obs: pd.DataFrame):
    """
    Assign each nucleus to a (subclass, donor) group.

    Returns:
        codes     : int array (n_obs,), group code per nucleus (-1 = excluded).
        group_meta: DataFrame indexed by group code with subclass, donor,
                    condition, covariates and n_nuclei.
    """
    keep = (
        obs[COND_COL].isin(GROUPS)
        & obs[SUBCLASS_COL].notna()
        & obs[DONOR_COL].notna()
    ).to_numpy()

    pair = (
        obs[SUBCLASS_COL].astype(str) + " ||| " + obs[DONOR_COL].astype(str)
    ).to_numpy()
    cats = pd.Categorical(pair[keep])

    codes = np.full(obs.shape[0], -1, dtype=np.int64)
    codes[keep] = cats.codes

    kept = obs.loc[keep].copy()
    kept["__code"] = cats.codes
    grp = kept.groupby("__code", observed=True)
    group_meta = pd.DataFrame({
        SUBCLASS_COL: grp[SUBCLASS_COL].first(),
        DONOR_COL: grp[DONOR_COL].first(),
        COND_COL: grp[COND_COL].first(),
        "n_nuclei": grp.size(),
    })
    for cov in [c for c in COVARIATES if c in obs.columns]:
        group_meta[cov] = grp[cov].first()
    group_meta = group_meta.sort_index()
    return codes, group_meta


# ==================================================
# Streaming aggregation
# ==================================================
def accumulate_pass(path: Path, codes: np.ndarray, batch_codes: np.ndarray,
                    n_peaks: int, chunk_size: int) -> np.ndarray:
    """
    One streaming pass over the backed matrix, summing counts for the groups in
    `batch_codes` only. Returns a dense (len(batch_codes) x n_peaks) accumulator.
    """
    n_groups = batch_codes.shape[0]
    acc = np.zeros((n_groups, n_peaks), dtype=np.float32)

    # Map global group code -> local accumulator row (-1 if not in this batch).
    code2local = np.full(codes.max() + 1 if codes.max() >= 0 else 1, -1,
                         dtype=np.int64)
    code2local[batch_codes] = np.arange(n_groups)

    adata = sc.read_h5ad(path, backed="r")
    try:
        n_obs = adata.n_obs
        X = adata.X
        for start in range(0, n_obs, chunk_size):
            stop = min(start + chunk_size, n_obs)
            cc = codes[start:stop]
            local = np.where(cc >= 0, code2local[np.clip(cc, 0, None)], -1)
            sel = np.nonzero(local >= 0)[0]
            if sel.size == 0:
                continue
            block = X[start:stop]
            block = block if sp.issparse(block) else sp.csr_matrix(block)
            block = block[sel]
            indicator = sp.csr_matrix(
                (np.ones(sel.size, dtype=np.float32),
                 (local[sel], np.arange(sel.size))),
                shape=(n_groups, sel.size),
            )
            prod = (indicator @ block).tocoo()
            np.add.at(acc, (prod.row, prod.col), prod.data)
    finally:
        adata.file.close()
    return acc


def plan_batches(group_meta: pd.DataFrame, n_peaks: int) -> list:
    """Group subclasses into batches whose accumulator fits MEM_BUDGET_GB."""
    max_groups = max(1, int(MEM_BUDGET_GB * 1e9 / (n_peaks * 4)))
    batches, current, used = [], [], 0
    for subclass, sub_meta in group_meta.groupby(SUBCLASS_COL, observed=True):
        size = sub_meta.shape[0]
        if size > max_groups:
            print(f"[warn] Subclass '{subclass}' needs {size} groups "
                  f"(> {max_groups} for {MEM_BUDGET_GB} GB); processed alone, "
                  f"may exceed the budget.")
        if current and used + size > max_groups:
            batches.append(current)
            current, used = [], 0
        current.append(subclass)
        used += size
    if current:
        batches.append(current)
    return batches


# ==================================================
# Writing outputs
# ==================================================
def write_subclass(subclass: str, acc: np.ndarray, local_codes: np.ndarray,
                   group_meta: pd.DataFrame, var: pd.DataFrame,
                   covariates: list, output_dir: Path) -> pd.DataFrame:
    """Slice the accumulator for one subclass, filter and write outputs."""
    sub_meta = group_meta.loc[local_codes].copy()
    donors = sub_meta[DONOR_COL].astype(str).to_numpy()
    mat = acc  # rows already aligned to local_codes order

    total_counts = mat.sum(axis=1)
    col_data = pd.DataFrame(index=pd.Index(donors, name=DONOR_COL))
    col_data[COND_COL] = sub_meta[COND_COL].to_numpy()
    for cov in covariates:
        col_data[cov] = sub_meta[cov].to_numpy()
    col_data["subclass"] = subclass
    col_data["n_nuclei"] = sub_meta["n_nuclei"].to_numpy()
    col_data["total_counts"] = total_counts

    keep = col_data["n_nuclei"].to_numpy() >= MIN_NUCLEI
    n_dropped = int((~keep).sum())
    if n_dropped:
        print(f"[{subclass}] dropping {n_dropped}/{len(donors)} donors "
              f"with < {MIN_NUCLEI} nuclei.")

    pb = ad.AnnData(
        X=sp.csr_matrix(np.rint(mat[keep])),
        obs=col_data.loc[keep].copy(),
        var=var.copy(),
    )
    safe = subclass.replace(" ", "-").replace("/", "-")
    output_dir.mkdir(parents=True, exist_ok=True)
    pb.write_h5ad(output_dir / f"pseudobulk_{safe}.h5ad", compression="gzip")
    col_data.loc[keep].to_csv(output_dir / f"colData_{safe}.csv",
                              index=True, index_label=DONOR_COL)
    print(f"[ok] {subclass}: {int(keep.sum())} donor pseudobulks x "
          f"{pb.n_vars} peaks.")

    summary = col_data.copy()
    summary["passed_filter"] = keep
    return summary.reset_index()


# ==================================================
# QC plots of the counts
# ==================================================
def plot_counts(summary: pd.DataFrame, fig_dir: Path) -> None:
    """Draw QC figures describing the pseudobulk counts."""
    fig_dir.mkdir(parents=True, exist_ok=True)
    kept = summary[summary["passed_filter"]].copy()
    if kept.empty:
        print("[warn] Nothing passed the filter; skipping plots.")
        return

    # 1) Nuclei per (donor x subclass) as a heatmap.
    mat = kept.pivot_table(index=DONOR_COL, columns="subclass",
                           values="n_nuclei", aggfunc="sum")
    fig, axis = plt.subplots(
        figsize=(max(6, 0.4 * mat.shape[1] + 4), max(4, 0.18 * mat.shape[0] + 2))
    )
    im = axis.imshow(np.log10(mat.to_numpy() + 1), aspect="auto", cmap="viridis")
    axis.set_xticks(range(mat.shape[1]))
    axis.set_xticklabels(mat.columns, rotation=90, fontsize=7)
    axis.set_yticks(range(mat.shape[0]))
    axis.set_yticklabels(mat.index, fontsize=5)
    axis.set_title("Nuclei per donor x subclass (log10)")
    fig.colorbar(im, ax=axis, label="log10(n_nuclei + 1)")
    fig.tight_layout()
    fig.savefig(fig_dir / "nuclei_per_donor_subclass.png", dpi=150)
    plt.close(fig)

    # 2) Pseudobulk library size per condition (log scale).
    order = [g for g in GROUPS if g in kept[COND_COL].unique()]
    fig, axis = plt.subplots(figsize=(6, 4))
    axis.boxplot([kept.loc[kept[COND_COL] == g, "total_counts"] for g in order])
    axis.set_xticks(range(1, len(order) + 1))
    axis.set_xticklabels(order)
    axis.set_yscale("log")
    axis.set_ylabel("Pseudobulk total counts (library size)")
    axis.set_xlabel("Condition")
    axis.set_title("Pseudobulk library size per condition")
    fig.tight_layout()
    fig.savefig(fig_dir / "library_size_per_condition.png", dpi=150)
    plt.close(fig)

    # 3) Donors passing the filter per condition and subclass.
    counts = (kept.groupby(["subclass", COND_COL], observed=True)[DONOR_COL]
              .nunique().unstack(COND_COL).reindex(columns=order).fillna(0))
    fig, axis = plt.subplots(figsize=(max(6, 0.5 * counts.shape[0] + 3), 4))
    bottom = np.zeros(counts.shape[0])
    for cond in order:
        axis.bar(counts.index, counts[cond], bottom=bottom, label=cond)
        bottom += counts[cond].to_numpy()
    axis.set_ylabel("Donors (pseudobulk samples)")
    axis.set_xlabel("Subclass")
    axis.set_title(f"Donors passing filter (>= {MIN_NUCLEI} nuclei) per subclass")
    axis.set_xticks(range(counts.shape[0]))
    axis.set_xticklabels(counts.index, rotation=90, fontsize=7)
    axis.legend(title="Condition", fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "donors_per_condition_subclass.png", dpi=150)
    plt.close(fig)

    print(f"[ok] QC figures written to {fig_dir}")


# ==================================================
# Main
# ==================================================
def main():
    print("Starting per-donor pseudobulk aggregation (streaming) ...")
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Input file: {INPUT_FILE}")
    print(f"Output directory: {OUTPUT_DIR}")

    obs = load_metadata(INPUT_FILE)
    covariates = check_columns(obs)
    codes, group_meta = build_groups(obs)

    # Peak metadata and matrix width (no counts loaded here).
    with h5py.File(INPUT_FILE, "r") as f:
        var = read_elem(f["var"])
    n_peaks = var.shape[0]

    n_groups = group_meta.shape[0]
    acc_gb = n_groups * n_peaks * 4 / 1e9
    print(f"{n_groups} (subclass, donor) groups x {n_peaks} peaks "
          f"(full accumulator would be {acc_gb:.1f} GB).")

    batches = plan_batches(group_meta, n_peaks)
    print(f"Processing in {len(batches)} streaming pass(es) "
          f"(budget {MEM_BUDGET_GB} GB).")

    summaries = []
    for i, subclasses in enumerate(batches, 1):
        batch_meta = group_meta[group_meta[SUBCLASS_COL].isin(subclasses)]
        batch_codes = batch_meta.index.to_numpy()
        print(f"[pass {i}/{len(batches)}] subclasses: {subclasses} "
              f"({batch_codes.size} groups)")
        acc = accumulate_pass(INPUT_FILE, codes, batch_codes, n_peaks, CHUNK_SIZE)

        for subclass in subclasses:
            sub_meta = group_meta[group_meta[SUBCLASS_COL] == subclass]
            sub_codes = sub_meta.index.to_numpy()
            rows = np.searchsorted(batch_codes, sub_codes)
            summaries.append(
                write_subclass(subclass, acc[rows], sub_codes, group_meta,
                               var, covariates, OUTPUT_DIR)
            )
        del acc  # free before the next pass

    summary = pd.concat([s for s in summaries if not s.empty], ignore_index=True)
    summary.to_csv(OUTPUT_DIR / "pseudobulk_summary.csv", index=False)
    print(f"[ok] Wrote summary: {OUTPUT_DIR / 'pseudobulk_summary.csv'}")

    plot_counts(summary, FIG_DIR)
    print("Per-donor pseudobulk aggregation completed.")


if __name__ == "__main__":
    main()
