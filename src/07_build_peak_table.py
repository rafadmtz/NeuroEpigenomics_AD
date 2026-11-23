import polars as pl
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DA_DIR = PROJECT_ROOT / "results" / "06_differential_accessibility"
OUT_DIR = PROJECT_ROOT / "results" / "07_build_peak_table"
FIG_DIR = PROJECT_ROOT / "results" / "figures" / "07_build_peak_table"

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

def read_da(subclass: str, peak_overlap: pl.DataFrame) -> pl.DataFrame:
    
    df_da= pl.read_csv(
        DA_DIR / f"da_{subclass}_trend.csv.gz"
        )
    
    df_merged = df_da.join(
        peak_overlap,
        on="peak_id",
        how="left"
    )

    df_merged = df_merged.filter(
        pl.col("repClass").is_not_null()
    )

    df_merged = df_merged.with_columns(
        (pl.col("log2FoldChange") > 0).alias("is_up"),
        (pl.col("log2FoldChange") < 0).alias("is_down"),
        pl.lit(subclass).alias("subclass")
    )
    
    return df_merged

def resume_results(df_summary: pd.DataFrame):
    
    thresholds= [0.05, 0.03, 0.01]
    
    for threshold in thresholds:
        
        df_filtered= df_summary.filter(
            pl.col("padj").is_not_nan() &
            (pl.col("padj") <= threshold)
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

def plot_tes(df_summary: pl.DataFrame):
    
    df_counts = (
        df_summary
        .filter(
            (pl.col("padj") < 0.05) &
            pl.col("repClass").is_not_null()
        )
        .group_by(["subclass", "repClass"])
        .agg(
            pl.col("is_up").fill_null(False).sum().alias("up"),
            pl.col("is_down").fill_null(False).sum().alias("down")
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
    ax.set_title("Differential accessibility by TE class (padj < 0.05)")
    ax.legend()

    plt.subplots_adjust(bottom=0.3)
    plt.tight_layout()

    plt.savefig(
        FIG_DIR / "te_differential_accessibility.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()

def hist_log2FoldChange(df_summary):
    
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
    
    PEAK_OVERLAP= pl.read_csv(
        "results/05_te_annotation/peak_te_overlaps.csv"
    )
    
    df_summary=[]
    
    for subclass in SUBCLASSES:
        
        df_subclass= read_da(subclass, PEAK_OVERLAP)
        
        df_summary.append(df_subclass)
        
    df_summary = pl.concat(df_summary)
    
    resume_results(df_summary)
    
    plot_tes(df_summary)
    
    hist_log2FoldChange(df_summary)
    
if __name__ == "__main__":
    main()