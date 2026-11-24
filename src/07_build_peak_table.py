#!/usr/bin/env python3
"""
07_build_peak_table.py — Peak table for TE enrichment
=====================================================

Joins the differential accessibility results of step 06 with the TE annotation
of step 05 and the GC content of each peak. The output has one row per tested
peak and subclass, and is the input of step 08.

Workflow:
  1. Computes GC per peak with `bedtools nuc` (only if peaks_nuc.txt is
     missing or empty).
  2. For a list of subclasses, reads the DA results, drops peaks without padj and
     joins GC and the binary TE matrices (peaks without TE are kept as
     background).
  3. Flags each peak as up, down or not significant (padj < 0.05).
  4. Prints counts per TE class and plots them.

Inputs:
  - data/genome/hg38.fa
  - results/01_metadata_output/atac_peaks.bed
  - results/05_te_annotation/binary_matrix_peaks_te_class.csv
  - results/05_te_annotation/binary_matrix_peaks_te_family.csv
  - results/06_differential_accessibility/da_{subclass}_{CONTRAST}.csv.gz
  - results/05_te_annotation/peak_te_overlaps.csv

Outputs (in results/07_build_peak_table/):
  - peaks_nuc.txt
  - peak_table_{subclass}_{CONTRAST}.parquet
        Columns: peak_id, subclass, baseMean, log2FoldChange, log2FC_shrunk,
        padj, is_up, is_down, is_ns, gc, any_TE, TE classes and TE families.
  - figures/07_build_peak_table/te_differential_accessibility*.png
  - figures/07_build_peak_table/hist_log2FoldChange.png

Usage:
  Execute directly from the command line from the repository root:
  python src/07_build_peak_table.py
  nohup python src/07_build_peak_table.py > logs/07_build_peak_table.log 2>&1 &
"""

import shutil
import subprocess
import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DA_DIR = PROJECT_ROOT / "results" / "06_differential_accessibility"
TE_DIR = PROJECT_ROOT / "results" / "05_te_annotation"
OUT_DIR = PROJECT_ROOT / "results" / "07_build_peak_table"
FIG_DIR = PROJECT_ROOT / "results" / "figures" / "07_build_peak_table"

GENOME_FA = PROJECT_ROOT / "data" / "genome" / "hg38.fa"
PEAKS_BED = PROJECT_ROOT / "results" / "01_metadata_output" / "atac_peaks.bed"
CLASS_MATRIX = TE_DIR / "binary_matrix_peaks_te_class.csv"
FAMILY_MATRIX = TE_DIR / "binary_matrix_peaks_te_family.csv"
NUC_PATH = OUT_DIR / "peaks_nuc.txt"

CONTRAST = "trend"
ALPHA = 0.05

TEs= ["LINE", "LTR", "SINE"]

SUBCLASSES = [
    # Excitatory neurons (hypothesis from Liu et al., ZBP1 in excitatory neurons)
    "L2-3-IT",
    "L4-IT",
    "L5-IT",
    "L6-IT",
    # Inhibitory neurons
    "Vip",
    "Sst",
    # Glia
    "Astrocyte",
    "Oligodendrocyte",
    "Microglia-PVM"
]

def compute_gc() -> pl.DataFrame:
    """
    Read the GC fraction of each peak, running `bedtools nuc` if needed.

    Returns
    -------
    pl.DataFrame
        One row per peak: peak_id, gc.

    Raises
    ------
    EnvironmentError
        If bedtools is needed but is not in PATH.
    FileNotFoundError
        If the genome or the peaks BED is missing.
    RuntimeError
        If `bedtools nuc` fails.
    ValueError
        If peaks_nuc.txt does not have the bedtools nuc columns.
    """

    if not NUC_PATH.exists() or NUC_PATH.stat().st_size == 0:

        if not shutil.which("bedtools"):
            raise EnvironmentError("bedtools is not in PATH (see environment.yml)")

        for path in [GENOME_FA, PEAKS_BED]:
            if not path.exists():
                raise FileNotFoundError(f"Input of bedtools nuc not found: {path}")

        cmd = ["bedtools", "nuc", "-fi", str(GENOME_FA), "-bed", str(PEAKS_BED)]

        try:
            with open(NUC_PATH, "w") as out_file:
                subprocess.run(cmd, stdout=out_file, check=True)
        except subprocess.CalledProcessError as e:
            NUC_PATH.unlink(missing_ok=True)  # do not leave an empty file behind
            raise RuntimeError("bedtools nuc failed, check the genome and the peaks BED") from e

    try:
        df_gc = pl.read_csv(NUC_PATH, separator="\t").select(
            pl.col("4_usercol").alias("peak_id"),
            pl.col("6_pct_gc").alias("gc")
        )
    except pl.exceptions.ColumnNotFoundError as e:
        raise ValueError(f"{NUC_PATH} does not look like bedtools nuc output") from e

    return df_gc

def read_te_matrix() -> pl.DataFrame:
    """
    Read the binary TE matrices of step 05 (one row per peak).

    Returns
    -------
    pl.DataFrame
        peak_id, any_TE and one 0/1 column per TE class and family.

    Raises
    ------
    FileNotFoundError
        If a matrix is missing.
    ValueError
        If the family matrix has no any_TE column.
    """

    try:
        df_class = pl.read_csv(CLASS_MATRIX)
        df_family = pl.read_csv(FAMILY_MATRIX).drop("any_TE")
    except FileNotFoundError as e:
        raise FileNotFoundError(f"TE matrix not found, run step 05 first: {e}") from e
    except pl.exceptions.ColumnNotFoundError as e:
        raise ValueError(f"Unexpected columns in the TE matrices: {e}") from e

    return df_class.join(df_family, on="peak_id", how="left")

def read_da(subclass: str, te_matrix: pl.DataFrame, df_gc: pl.DataFrame) -> pl.DataFrame:
    """
    Build the peak table of one subclass.

    Peaks without padj (not tested by DESeq2) are dropped. Peaks without TE
    are kept with every TE column set to 0.

    Parameters
    ----------
    subclass : str
        Subclass name, as in the DA file names.
    te_matrix : pl.DataFrame
        Output of read_te_matrix.
    df_gc : pl.DataFrame
        Output of compute_gc.

    Returns
    -------
    pl.DataFrame
        One row per tested peak, with is_up, is_down and is_ns.

    Raises
    ------
    FileNotFoundError
        If the DA file of `subclass` is missing.
    ValueError
        If the file has no padj column or no tested peaks.
    """

    da_path = DA_DIR / f"da_{subclass}_{CONTRAST}.csv.gz"

    try:
        df_da = pl.read_csv(da_path)
        df_da = df_da.filter(
            pl.col("padj").is_not_null() &
            pl.col("padj").is_not_nan()
        )
    except FileNotFoundError as e:
        raise FileNotFoundError(f"DA results of {subclass} not found: {da_path}") from e
    except pl.exceptions.ColumnNotFoundError as e:
        raise ValueError(f"{da_path.name} has no padj column") from e
    except pl.exceptions.InvalidOperationError as e:
        raise ValueError(f"padj is not numeric in {da_path.name}") from e

    if df_da.is_empty():
        raise ValueError(f"No tested peaks (padj) for {subclass}")

    te_cols = [col for col in te_matrix.columns if col != "peak_id"]

    df_merged = (
        df_da
        .join(df_gc, on="peak_id", how="left")
        .join(te_matrix, on="peak_id", how="left")
        .with_columns(
            pl.col(te_cols).fill_null(0)
        )
    )

    df_merged = df_merged.with_columns(
        ((pl.col("padj") <= ALPHA) & (pl.col("log2FoldChange") > 0)).alias("is_up"),
        ((pl.col("padj") <= ALPHA) & (pl.col("log2FoldChange") < 0)).alias("is_down"),
        (pl.col("padj") > ALPHA).alias("is_ns"),
        pl.lit(subclass).alias("subclass")
    )

    first_cols = ["peak_id", "subclass", "baseMean", "log2FoldChange",
                  "log2FC_shrunk", "padj", "is_up", "is_down", "is_ns", "gc"]

    return df_merged.select(first_cols + te_cols)

def summarize_results(df_summary: pl.DataFrame):
    """
    Print the number of significant peaks per subclass and TE class.

    Parameters
    ----------
    df_summary : pl.DataFrame
        Peak table joined with peak_te_overlaps.csv (one row per peak and repClass).
    """
    
    thresholds= [0.05, 0.03, 0.01]
    
    for threshold in thresholds:
        
        df_filtered= df_summary.filter(
            pl.col("padj") <= threshold
        )

        df_counts = (
            df_filtered
            .group_by(["subclass", "repClass"])
            .len()
            .pivot(
                on="repClass",
                index="subclass",
                values="len",
            )
            .fill_null(0)
        )

        df_status = (
            df_filtered
            .group_by(["subclass", "repClass"])
            .agg(
                pl.col("is_up").sum().alias("up"),
                pl.col("is_down").sum().alias("down")
            )
        )

        print(f"Threshold: {threshold}")
        print()
        print("Counts:")
        print(df_counts)
        print()
        print("Status:")
        print(df_status)

def plot_tes(df_summary: pl.DataFrame, lfc: float =0):
    """
    Plot up and down peaks per TE class and subclass (padj < 0.05).

    Parameters
    ----------
    df_summary : pl.DataFrame
        Peak table joined with peak_te_overlaps.csv (one row per peak and repClass).
    lfc : float
        Minimum |log2FoldChange|; 0 means no cutoff. Default 0.
    """
    
    if lfc:
        df_counts = (
            df_summary
            .filter(
                (pl.col("padj") < 0.05) &
                (pl.col("log2FoldChange").abs() > lfc)
            )
            .group_by(["subclass", "repClass"])
            .agg(
                pl.col("is_up").sum().alias("up"),
                pl.col("is_down").sum().alias("down")
            )
        )

    else:
        df_counts = (
            df_summary
            .filter(
                (pl.col("padj") < 0.05)
            )
            .group_by(["subclass", "repClass"])
            .agg(
                pl.col("is_up").sum().alias("up"),
                pl.col("is_down").sum().alias("down")
            )
        )

    repclasses = sorted(df_counts["repClass"].unique().to_list())

    fig, ax = plt.subplots(figsize=(14, 6))

    width = 0.35
    gap = 1.5

    xticks = []
    xlabels = []

    for i, subclass in enumerate(SUBCLASSES):

        df_sub = df_counts.filter(
            pl.col("subclass") == subclass
        )

        counts = {
            row["repClass"]: (row["up"], row["down"])
            for row in df_sub.iter_rows(named=True)
        }

        up = [counts.get(rep, (0, 0))[0] for rep in repclasses]
        down = [counts.get(rep, (0, 0))[1] for rep in repclasses]

        x = np.arange(len(repclasses)) + i * (len(repclasses) + gap)

        ax.bar(
            x - width / 2,
            up,
            width,
            label="Up" if i == 0 else None
        )

        ax.bar(
            x + width / 2,
            -np.array(down),
            width,
            label="Down" if i == 0 else None
        )

        xticks.extend(x)
        xlabels.extend(repclasses)

        center = (x[0] + x[-1]) / 2

        ax.text(
            center, -0.26, subclass,
            transform=ax.get_xaxis_transform(),
            ha="center",
            fontweight="bold"
        )

    ax.axhline(0, color="black", linewidth=0.8)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=45, ha="right")

    ax.set_ylabel("Number of differentially accessible peaks")

    title = f"Differential accessibility by TE class (padj < 0.05, |log2FC| > {lfc})" if lfc \
        else "Differential accessibility by TE class (padj < 0.05)"

    ax.set_title(title)

    ax.legend()

    plt.subplots_adjust(bottom=0.3)
    plt.tight_layout()

    if lfc:
        filename = f"te_differential_accessibility_lfc{lfc}.png"
    else:
        filename = "te_differential_accessibility.png"

    fig.savefig(
        FIG_DIR / filename,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()

def hist_log2FoldChange(df_summary: pl.DataFrame):
    """
    Plot the log2FoldChange distribution of significant peaks (padj <= 0.05).

    Parameters
    ----------
    df_summary : pl.DataFrame
        Concatenated output of read_da.
    """

    df_filtered = df_summary.filter(
        pl.col("padj") <= 0.05
    )

    df_pandas = df_filtered.to_pandas()

    g =sns.displot(
        df_pandas,
        x="log2FoldChange",
        binwidth=0.05,
        color= "#D43636"
        )

    g.savefig(
        FIG_DIR / "hist_log2FoldChange.png",
        dpi=300,
        bbox_inches="tight"
        )

    plt.show()


def main():

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    df_gc = compute_gc()

    te_matrix = read_te_matrix()

    df_summary=[]

    for subclass in SUBCLASSES:

        try:
            df_subclass= read_da(subclass, te_matrix, df_gc)
        except (FileNotFoundError, ValueError) as e:
            print(f"Skipping {subclass}: {e}")
            continue

        df_subclass.write_parquet(
            OUT_DIR / f"peak_table_{subclass}_{CONTRAST}.parquet"
        )

        df_summary.append(df_subclass)

    if not df_summary:
        raise RuntimeError("No subclass could be processed, check results/06_differential_accessibility")

    df_summary = pl.concat(df_summary)

    try:
        peak_overlap= (
            pl.read_csv(TE_DIR / "peak_te_overlaps.csv")
            .select("peak_id", "repClass")
            .unique()
        )
    except FileNotFoundError as e:
        raise FileNotFoundError(f"TE overlaps not found, run step 05 first: {e}") from e
    except pl.exceptions.ColumnNotFoundError as e:
        raise ValueError(f"peak_te_overlaps.csv needs the columns peak_id and repClass: {e}") from e

    df_plot = df_summary.join(
        peak_overlap,
        on="peak_id",
        how="inner"
    )

    if df_plot.is_empty():
        raise ValueError("No peak overlaps a TE, check that peak_id has the same format in both tables")

    summarize_results(df_plot)
    
    lfcs= [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    
    for lfc in lfcs:
        plot_tes(df_plot, lfc= lfc)
    
    hist_log2FoldChange(df_summary)

if __name__ == "__main__":
    main()