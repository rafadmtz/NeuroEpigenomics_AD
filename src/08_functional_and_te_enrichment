#!/usr/bin/env python3
"""
08_functional_and_te_enrichment.py — Functional and TE enrichment of differentially accessible peaks
==========================================================================

Takes the differential accessibility (DA) results of step 06 (trend contrast:
change in accessibility per step of AD neuropathological change) and asks two
questions for each subclass in SUBCLASSES (L2-3-IT, L4-IT, L5-IT, Vip, Sst
and Astrocyte, the populations with enough DA peaks and of interest to the
project). Every table and graph covers only these subclasses, and the
Benjamini-Hochberg correction is applied over them.

  A. Functional enrichment: are the genes near the peaks that gain (up) or lose
     (down) accessibility over-represented in Hallmark or GO Biological
     Process gene sets?
  B. Transposable elements (TEs): do peaks whose center falls inside a TE
     (step 05) change more than comparable peaks without a TE?

Why covariates are needed in B: in step 06, weakly accessible peaks gain
accessibility with stage much more often than strongly accessible ones (in
L2-3-IT, 91% of the DA peaks of the lowest baseMean quintile are up, against
63% in the highest quintile). TEs fall mostly in weak, distal peaks, so a crude
test could attribute this global effect to the TEs. Every TE test is therefore
also run adjusted for accessibility (log baseMean and its square), GC content,
peak width and distance to the nearest TSS. The same covariates are used in
the functional enrichment (see step 5).

Requires having previously run src/00, src/01, src/05 and src/06.

Workflow:
  1. Reads the protein-coding genes of the GENCODE GTF and assigns each peak
     center to its nearest TSS (gene, distance and a promoter / proximal /
     distal category).
  2. Computes the GC content of every peak with `bedtools nuc` (the output is
     reused if it already exists).
  3. Builds a per-peak TE table: class and family flags (step 05), a flag for
     young TE subfamilies and the divergence (milliDiv) of the TE.
  4. For each subclass of SUBCLASSES, keeps the tested peaks (padj not NaN) as universe and
     flags peaks as up or down (padj < ALPHA).
  5. Functional enrichment, for each direction with >= MIN_DA_PEAKS peaks:
     a gene is a hit if it has >= 1 DA peak, and its probability of being a
     hit by chance comes from a logistic model of its number of peaks and
     their covariates (in the spirit of GOseq). Genes of the same gene set
     whose TSSs are less than LOCUS_DIST apart (gene clusters such as
     olfactory receptors or histones) are merged into one locus, and the test
     counts hit loci, so one region that changes accessibility counts once.
     Run with all peaks and with promoter peaks only.
     Discarded alternatives: a peak-level hypergeometric test gave ~16% of
     p < 0.05 with random gene sets (genes with many weak, distal peaks are
     hit more often), and counting genes instead of loci made single gene
     clusters (olfactory receptors, histones) look like enriched pathways.
  6. TE enrichment, for each direction with >= MIN_DA_PEAKS peaks: Fisher's
     exact test (crude) and logistic regression adjusted for covariates.
  7. TE shift (no padj threshold): linear regression of the
     DESeq2 Wald statistic on the TE flag plus covariates.
  8. TE age: Spearman correlation between the covariate-adjusted Wald
     statistic and the milliDiv of the TE, per TE class.
  9. Benjamini-Hochberg correction and tables.
 10. Two functional enrichment dotplots (up and down peaks), made from the
     saved table of significant terms.

Inputs:
  - results/06_differential_accessibility/da_{subclass}_trend.csv.gz   (step 06,
                                                    one per subclass of SUBCLASSES)
  - results/05_te_annotation/binary_matrix_peaks_te_class.csv          (step 05)
  - results/05_te_annotation/binary_matrix_peaks_te_family.csv         (step 05)
  - results/05_te_annotation/peak_te_overlaps.csv                      (step 05)
  - results/01_metadata_output/atac_peaks.bed                          (step 01)
  - data/genome/hg38.fa                                                (step 00)
  - data/annotation/gencode.v50.basic.annotation.gtf.gz                (step 00)
  - data/genesets/{h.all, c5.go.bp}.v2026.2.Hs.symbols.gmt           (step 00)

Outputs (in results/08_enrichment/):
  - peaks_nuc.txt                    bedtools nuc output (GC per peak).
  - peak_annotation.csv.gz           One row per peak: nearest gene, distance
                                     to TSS, TSS category, GC, TE flags.
  - functional_enrichment.csv.gz     Every gene set tested (subclass x direction
                                     x scope x collection).
  - functional_enrichment_significant.csv
                                     padj < ALPHA and >= MIN_LOCI_HIT loci.
  - te_enrichment_class.csv          Fisher and adjusted logistic regression,
                                     TE classes, any_TE and young_TE.
  - te_enrichment_family.csv         Same, for TE families.
  - te_shift_ols.csv                 Shift of the Wald statistic (classes and
                                     families).
  - te_age.csv                       Spearman stat vs milliDiv per TE class.
  - enrichment_summary.csv           One row per subclass: peaks tested,
                                     up/down and which analyses were run.
  - 08_enrichment_log.txt

Every table has a `significant` column (padj < ALPHA, plus >= MIN_LOCI_HIT
hit loci for the functional enrichment), used by the graphs.

Graphs (in results/08_graphs/):
  1.- 08_functional_dotplot_up.png     Peaks that gain accessibility with stage.
  2.- 08_functional_dotplot_down.png   Peaks that lose accessibility with stage.
  Each graph has one panel per collection (Hallmark, GO_BP) and one column
  per subclass of SUBCLASSES, with the top N_TERMS_DOTPLOT significant terms
  (scope 'all') of each subclass; dot size = hit loci, colour = -log10(padj)
  (same scale in both panels). An empty column means no significant terms.

Usage:
  Execute directly from the command line from the repository root:
  python src/08_functional_and_te_enrichment.py

  Run in the background and save output to a log:
  nohup python -u src/08_functional_and_te_enrichment.py > logs/08_functional_and_te_enrichment.log 2>&1 &

  Or import as a module and call the main() function.
"""

#==================================================
# Libraries used
#==================================================

import gzip
import os
import re
import shutil
import subprocess
from pathlib import Path

# Limit the threads of the linear algebra libraries (OpenBLAS / MKL / OpenMP)
# used by numpy and scipy. By default they start one thread per core of the
# server. This must be set before numpy is imported; a value already exported
# in the shell (e.g. OMP_NUM_THREADS=2 python ...) takes precedence.
N_THREADS = 4
for _var in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ.setdefault(_var, str(N_THREADS))

import matplotlib
matplotlib.use("Agg")  # Headless environment for server

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import sparse, special, stats

#==================================================
# Variable definitions
#==================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"

# Inputs
DA_DIR = RESULTS_DIR / "06_differential_accessibility"
TE_DIR = RESULTS_DIR / "05_te_annotation"
CLASS_MATRIX_CSV = TE_DIR / "binary_matrix_peaks_te_class.csv"
FAMILY_MATRIX_CSV = TE_DIR / "binary_matrix_peaks_te_family.csv"
OVERLAPS_CSV = TE_DIR / "peak_te_overlaps.csv"
PEAKS_BED_PATH = RESULTS_DIR / "01_metadata_output" / "atac_peaks.bed"
GENOME_FA = DATA_DIR / "genome" / "hg38.fa"
GTF_PATH = DATA_DIR / "annotation" / "gencode.v50.basic.annotation.gtf.gz"
GENESETS_DIR = DATA_DIR / "genesets"
MSIGDB_VERSION = "2026.2.Hs"
COLLECTIONS = {
    "Hallmark": GENESETS_DIR / f"h.all.v{MSIGDB_VERSION}.symbols.gmt",
    "GO_BP": GENESETS_DIR / f"c5.go.bp.v{MSIGDB_VERSION}.symbols.gmt",
}

# Outputs
OUT_DIR = RESULTS_DIR / "08_enrichment"
GRAPHS_DIR = RESULTS_DIR / "08_graphs"
NUC_PATH = OUT_DIR / "peaks_nuc.txt"   # kept between runs (bedtools nuc is slow)
ANNOT_CSV = OUT_DIR / "peak_annotation.csv.gz"
FUNC_CSV = OUT_DIR / "functional_enrichment.csv.gz"
FUNC_SIG_CSV = OUT_DIR / "functional_enrichment_significant.csv"
TE_CLASS_CSV = OUT_DIR / "te_enrichment_class.csv"
TE_FAMILY_CSV = OUT_DIR / "te_enrichment_family.csv"
TE_OLS_CSV = OUT_DIR / "te_shift_ols.csv"
TE_AGE_CSV = OUT_DIR / "te_age.csv"
SUMMARY_CSV = OUT_DIR / "enrichment_summary.csv"
LOG_PATH = OUT_DIR / "08_enrichment_log.txt"

OUTPUTS = [
    ANNOT_CSV, FUNC_CSV, FUNC_SIG_CSV, TE_CLASS_CSV, TE_FAMILY_CSV,
    TE_OLS_CSV, TE_AGE_CSV, SUMMARY_CSV, LOG_PATH
]

# Subclasses analysed (tables and graphs)
SUBCLASSES = ["L2-3-IT", "L4-IT", "L5-IT", "Vip", "Sst", "Astrocyte"]

# Differential accessibility
CONTRAST = "trend"            # contrast of step 06 used as the consensus result
ALPHA = 0.05                  # padj threshold for DA peaks and for enrichment
MIN_DA_PEAKS = 50             # minimum up (or down) peaks to run threshold tests

# Peak-to-gene assignment
GENE_TYPES = ["protein_coding"]
PROMOTER_DIST = 2_000         # |distance to TSS| <= 2 kb   -> promoter
PROXIMAL_DIST = 50_000        # |distance to TSS| <= 50 kb  -> proximal, else distal
MAX_GENE_DIST = 1_000_000     # peaks farther than 1 Mb get no gene (as in GREAT)
CANONICAL_CHR = r"^chr([\d]+|[XY])$"

# Functional enrichment
MIN_SET_SIZE = 10             # gene set size, counted on genes of the universe
MAX_SET_SIZE = 500
LOCUS_DIST = 500_000          # genes of a set with TSSs < 500 kb apart (chained) = 1 locus
MIN_LOCI_HIT = 3              # a significant term needs >= 3 independent hit loci
N_GENES_REPORTED = 20         # hit genes listed per significant term

# Graphs
N_TERMS_DOTPLOT = 5           # top terms per subclass and collection in each dotplot
MAX_TERM_LABEL = 70           # characters of a term name shown on the axis
DPI = 300

# TE tests
TE_CLASSES = ["any_TE", "LINE", "SINE", "LTR", "DNA", "Retroposon"]
MIN_TE_PEAKS = 30             # minimum TE peaks in the universe to test a feature
# Young (recently active) subfamilies; L1HS is listed even if few peaks reach it.
YOUNG_TE_RE = re.compile(r"^(L1HS|L1PA[2-8]|AluY[a-z0-9]*|HERVK-int|LTR5_Hs|SVA_[A-F])$")
RIDGE = 1e-3                  # ridge penalty of the logistic models (see fit_logistic)
COVARIATES = ["log_baseMean", "log_baseMean_sq", "gc", "log_width", "proximal", "distal"]

# Columns of the output tables (an empty table is still written with them)
FUNC_COLUMNS = ["subclass", "direction", "scope", "collection", "term", "set_size",
                "n_loci", "n_genes_hit", "n_loci_hit", "expected_loci_hit",
                "fold_enrichment", "z", "pvalue", "padj", "significant", "genes"]
TE_THRESHOLD_COLUMNS = ["subclass", "direction", "feature", "n_te_peaks", "n_te_da", "n_da",
                        "pct_da_te", "pct_da_no_te", "crude_log2OR", "crude_pvalue",
                        "adj_log2OR", "adj_log2OR_ci_low", "adj_log2OR_ci_high",
                        "adj_pvalue", "crude_padj", "adj_padj", "significant"]
TE_OLS_COLUMNS = ["subclass", "level", "feature", "n_te_peaks", "n_peaks", "mean_stat_te",
                  "mean_stat_no_te", "shift", "shift_se", "pvalue", "padj", "significant"]
TE_AGE_COLUMNS = ["subclass", "te_class", "n_te_peaks", "spearman_rho", "pvalue", "padj",
                  "significant"]

#==================================================
# Functions
#==================================================

def log_print(msg):
    """
    Print a message and append it immediately to the log file.

    Parameters
    ----------
    msg : str
        Message to print and log.
    """
    print(msg)
    with open(LOG_PATH, 'a') as f:
        f.write(msg + "\n")

def check_dependencies():
    """
    Check that every input exists before starting.

    bedtools, atac_peaks.bed and hg38.fa are only required when the GC
    content has not been computed yet (peaks_nuc.txt missing or empty).

    Raises
    ------
    EnvironmentError
        If bedtools is needed and not in PATH.
    FileNotFoundError
        If an input file or the DA table of a subclass of SUBCLASSES is missing.
    """
    inputs = [CLASS_MATRIX_CSV, FAMILY_MATRIX_CSV, OVERLAPS_CSV, GTF_PATH,
              *COLLECTIONS.values()]
    if not NUC_PATH.exists() or NUC_PATH.stat().st_size == 0:
        inputs += [PEAKS_BED_PATH, GENOME_FA]
        if not shutil.which("bedtools"):
            raise EnvironmentError(
                "CRITICAL ERROR: 'bedtools' is not installed or not in PATH "
                "(needed for bedtools nuc). Install it with "
                "'conda install -c bioconda bedtools'."
            )
    missing = [str(p) for p in inputs if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing inputs:\n  " + "\n  ".join(missing))
    missing_da = [sc for sc in SUBCLASSES if not da_path(sc).exists()]
    if missing_da:
        raise FileNotFoundError(f"Missing DA results (step 06) in {DA_DIR} for: "
                                + ", ".join(missing_da))

def da_path(subclass):
    """
    Build the path of the DA table of step 06 for one subclass.

    Parameters
    ----------
    subclass : str
        Subclass name, as in the DA file names (e.g. "L2-3-IT").

    Returns
    -------
    Path
        results/06_differential_accessibility/da_{subclass}_{CONTRAST}.csv.gz
    """
    return DA_DIR / f"da_{subclass}_{CONTRAST}.csv.gz"

def benjamini_hochberg(pvalues):
    """
    Adjust p-values with the Benjamini-Hochberg (FDR) procedure.

    Parameters
    ----------
    pvalues : array-like
        Raw p-values; NaN values are ignored.

    Returns
    -------
    np.ndarray
        Adjusted p-values in the same order; NaN where the input was NaN.
    """
    p = np.asarray(pvalues, dtype=float)
    padj = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    n = ok.sum()
    if n == 0:
        return padj
    order = np.argsort(p[ok])
    ranked = p[ok][order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    padj[ok] = out
    return padj

def concat_tables(tables, columns):
    """
    Concatenate per-subclass result tables, skipping empty ones.

    Parameters
    ----------
    tables : list of pd.DataFrame
        Result tables to concatenate.
    columns : list of str
        Columns of the empty table returned when there is nothing to concatenate.

    Returns
    -------
    pd.DataFrame
        Concatenated table, or an empty table with `columns`.
    """
    tables = [t for t in tables if not t.empty]
    if not tables:
        return pd.DataFrame(columns=columns)
    return pd.concat(tables, ignore_index=True)

#--------------------------------------------------
# Peak annotation
#--------------------------------------------------

def read_tss():
    """
    Read one transcription start site (TSS) per gene from the GENCODE GTF.

    Only `gene` lines of the types in GENE_TYPES on canonical chromosomes are
    kept. The GTF is 1-based: the TSS is `start` on the + strand and `end` on
    the - strand, converted to a 0-based coordinate.

    Returns
    -------
    pd.DataFrame
        One row per gene, sorted by position: chr, tss, strand, gene.
    """
    log_print(f"\nReading gene TSSs from: {GTF_PATH.name}...")
    name_re = re.compile(r'gene_name "([^"]+)"')
    type_re = re.compile(r'gene_type "([^"]+)"')
    rows = []
    with gzip.open(GTF_PATH, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "gene":
                continue
            gene_type = type_re.search(fields[8]).group(1)
            if gene_type not in GENE_TYPES:
                continue
            # GTF is 1-based inclusive: the TSS is `start` on + and `end` on -.
            tss = int(fields[3]) - 1 if fields[6] == "+" else int(fields[4]) - 1
            rows.append((fields[0], tss, fields[6], name_re.search(fields[8]).group(1)))

    tss = pd.DataFrame(rows, columns=["chr", "tss", "strand", "gene"])
    tss = tss[tss["chr"].str.match(CANONICAL_CHR)]
    tss = tss.sort_values(["chr", "tss"]).reset_index(drop=True)
    log_print(f"Genes ({', '.join(GENE_TYPES)}) on canonical chromosomes: {tss['gene'].nunique()}")
    return tss

def assign_nearest_tss(peaks, tss):
    """
    Assign each peak center to its nearest TSS on the same chromosome.

    Parameters
    ----------
    peaks : pd.DataFrame
        Indexed by peak_id, with chr, start and end.
    tss : pd.DataFrame
        Output of read_tss.

    Returns
    -------
    pd.DataFrame
        Indexed by peak_id: gene (NA if farther than MAX_GENE_DIST), dist_tss
        (signed, negative upstream of the TSS following the gene strand) and
        tss_category (promoter / proximal / distal).
    """
    log_print("\nAssigning peaks to the nearest TSS...")
    centers = (peaks["start"] + peaks["end"]) // 2   # same center as step 05
    out = pd.DataFrame({"gene": pd.NA, "dist_tss": np.nan}, index=peaks.index)

    for chrom, idx in peaks.groupby("chr").groups.items():
        chr_tss = tss[tss["chr"] == chrom]
        if chr_tss.empty:
            continue
        pos = chr_tss["tss"].to_numpy()
        c = centers.loc[idx].to_numpy()
        right = np.clip(np.searchsorted(pos, c), 0, len(pos) - 1)
        left = np.clip(right - 1, 0, len(pos) - 1)
        nearest = np.where(np.abs(c - pos[left]) <= np.abs(c - pos[right]), left, right)
        dist = c - pos[nearest]
        dist = np.where(chr_tss["strand"].to_numpy()[nearest] == "+", dist, -dist)
        out.loc[idx, "gene"] = chr_tss["gene"].to_numpy()[nearest]
        out.loc[idx, "dist_tss"] = dist

    abs_dist = out["dist_tss"].abs()
    out["tss_category"] = np.select(
        [abs_dist <= PROMOTER_DIST, abs_dist <= PROXIMAL_DIST],
        ["promoter", "proximal"], default="distal")
    out.loc[abs_dist > MAX_GENE_DIST, "gene"] = pd.NA

    log_print("Peaks per TSS category: "
              + out["tss_category"].value_counts().to_dict().__repr__())
    log_print(f"Peaks without a gene within {MAX_GENE_DIST:,} bp: {out['gene'].isna().sum()}")
    return out

def compute_gc():
    """
    Read the GC fraction of each peak, running `bedtools nuc` if needed.

    The bedtools output (peaks_nuc.txt) is kept between runs, because
    computing it again takes a while.

    Returns
    -------
    pd.Series
        GC fraction, indexed by peak_id, named "gc".
    """
    if not NUC_PATH.exists() or NUC_PATH.stat().st_size == 0:
        log_print(f"\nRunning bedtools nuc on {PEAKS_BED_PATH.name} (this may take a while)...")
        cmd = ["bedtools", "nuc", "-fi", str(GENOME_FA), "-bed", str(PEAKS_BED_PATH)]
        try:
            with open(NUC_PATH, "w") as out_file:
                subprocess.run(cmd, stdout=out_file, check=True)
        except subprocess.CalledProcessError as e:
            NUC_PATH.unlink(missing_ok=True)
            log_print(f"ERROR executing bedtools nuc: {e}")
            raise
    else:
        log_print(f"\nReusing GC content from: {NUC_PATH.name}")

    nuc = pd.read_csv(NUC_PATH, sep="\t", usecols=["4_usercol", "6_pct_gc"])
    return nuc.set_index("4_usercol")["6_pct_gc"].rename("gc")

def build_te_table():
    """
    Build the per-peak TE table from the outputs of step 05.

    A peak center can overlap more than one TE fragment; the minimum (youngest)
    milliDiv per TE class is kept.

    Returns
    -------
    te : pd.DataFrame
        Indexed by peak_id: any_TE, one 0/1 column per TE class and family,
        young_TE and one milliDiv_{class} column per TE class.
    families : list of str
        Names of the TE family columns.
    """
    log_print("\nReading TE annotation (step 05)...")
    te_class = pd.read_csv(CLASS_MATRIX_CSV, index_col="peak_id")
    te_family = pd.read_csv(FAMILY_MATRIX_CSV, index_col="peak_id").drop(columns="any_TE")
    overlaps = pd.read_csv(OVERLAPS_CSV)

    young_peaks = overlaps.loc[overlaps["repName"].str.match(YOUNG_TE_RE), "peak_id"].unique()
    te_class["young_TE"] = te_class.index.isin(young_peaks).astype(int)

    # A peak center can overlap more than one TE fragment: keep the youngest.
    milli = (overlaps.groupby(["peak_id", "repClass"])["milliDiv"].min()
             .unstack().add_prefix("milliDiv_"))

    te = te_class.join(te_family).join(milli)
    log_print(f"Peaks with TE: {int(te['any_TE'].sum())} | young TE: {int(te['young_TE'].sum())} | "
              f"families: {te_family.shape[1]}")
    return te, list(te_family.columns)

def read_gmt(path):
    """
    Read an MSigDB GMT file.

    Parameters
    ----------
    path : Path
        GMT file (one gene set per line: name, URL, genes...).

    Returns
    -------
    dict
        {gene set name: set of gene symbols}.
    """
    gene_sets = {}
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            gene_sets[fields[0]] = set(fields[2:])
    return gene_sets

def membership_matrix(gene_sets, gene_index):
    """
    Build the gene set membership matrix.

    Parameters
    ----------
    gene_sets : dict
        Output of read_gmt.
    gene_index : np.ndarray
        Gene symbols (columns of the matrix); genes not in it are ignored.

    Returns
    -------
    scipy.sparse.csr_matrix
        (gene sets x genes) 0/1 matrix.
    """
    pos = {g: i for i, g in enumerate(gene_index)}
    rows, cols = [], []
    for i, genes in enumerate(gene_sets.values()):
        hits = [pos[g] for g in genes if g in pos]
        rows.extend([i] * len(hits))
        cols.extend(hits)
    data = np.ones(len(rows), dtype=np.float64)
    return sparse.csr_matrix((data, (rows, cols)), shape=(len(gene_sets), len(gene_index)))

def locus_matrices(gene_sets, gene_index, gene_pos):
    """
    Group the genes of every gene set into loci.

    Genes of the same set on the same chromosome whose TSSs are less than
    LOCUS_DIST apart (chaining consecutive genes) form one locus, so a gene
    cluster (e.g. olfactory receptors, histones) counts as a single locus.

    Parameters
    ----------
    gene_sets : dict
        Output of read_gmt.
    gene_index : np.ndarray
        Gene symbols (columns of L).
    gene_pos : pd.DataFrame
        chr and tss of every gene of gene_index, in the same order.

    Returns
    -------
    L : scipy.sparse.csr_matrix
        (loci x genes) 0/1 matrix.
    S : scipy.sparse.csr_matrix
        (gene sets x loci) 0/1 matrix.
    """
    pos = {g: i for i, g in enumerate(gene_index)}
    chrom = gene_pos["chr"].to_numpy()
    tss = gene_pos["tss"].to_numpy()
    l_rows, l_cols, s_rows, s_cols = [], [], [], []
    n_loci = 0
    for i, genes in enumerate(gene_sets.values()):
        idx = np.array([pos[g] for g in genes if g in pos], dtype=int)
        if idx.size == 0:
            continue
        idx = idx[np.lexsort((tss[idx], chrom[idx]))]
        new_locus = np.ones(idx.size, dtype=bool)
        new_locus[1:] = (chrom[idx][1:] != chrom[idx][:-1]) | (np.diff(tss[idx]) >= LOCUS_DIST)
        locus = n_loci + np.cumsum(new_locus) - 1
        l_rows.extend(locus)
        l_cols.extend(idx)
        s_rows.extend([i] * int(new_locus.sum()))
        s_cols.extend(np.unique(locus))
        n_loci = int(locus[-1]) + 1
    L = sparse.csr_matrix((np.ones(len(l_rows)), (l_rows, l_cols)),
                          shape=(n_loci, len(gene_index)))
    S = sparse.csr_matrix((np.ones(len(s_rows)), (s_rows, s_cols)),
                          shape=(len(gene_sets), n_loci))
    return L, S

#--------------------------------------------------
# Per-subclass peak table
#--------------------------------------------------

def load_da(path, annotation):
    """
    Build the peak table of one subclass.

    Peaks without padj (not tested by DESeq2) are dropped. Peaks are flagged
    as up or down (padj < ALPHA) and joined with the per-peak annotation.

    Parameters
    ----------
    path : Path
        DA table of step 06 (output of da_path).
    annotation : pd.DataFrame
        Per-peak annotation: nearest gene, TSS category, GC and TE columns.

    Returns
    -------
    pd.DataFrame
        One row per tested peak, indexed by peak_id, with is_up, is_down,
        the annotation and the covariates (log_baseMean, gc, log_width,
        proximal, distal).
    """
    da = pd.read_csv(path)
    da = da[da["padj"].notna()].set_index("peak_id")
    da["is_up"] = ((da["padj"] < ALPHA) & (da["log2FoldChange"] > 0)).astype(int)
    da["is_down"] = ((da["padj"] < ALPHA) & (da["log2FoldChange"] < 0)).astype(int)
    df = da.join(annotation, how="left")

    df["log_baseMean"] = np.log10(df["baseMean"])
    df["log_baseMean_sq"] = 0.0          # filled in covariate_matrix (after z-scoring)
    df["log_width"] = np.log10(df["end"] - df["start"])
    df["proximal"] = (df["tss_category"] == "proximal").astype(int)
    df["distal"] = (df["tss_category"] == "distal").astype(int)
    df["gc"] = df["gc"].fillna(df["gc"].median())
    return df

def covariate_matrix(df):
    """
    Build the design matrix of the covariates shared by the TE models.

    Continuous covariates are z-scored, and the square of log_baseMean is
    added because the effect of accessibility is not linear.

    Parameters
    ----------
    df : pd.DataFrame
        Output of load_da.

    Returns
    -------
    np.ndarray
        (peaks x (1 + len(COVARIATES))) matrix, first column = intercept.
    """
    cov = df[COVARIATES].astype(float).copy()
    for c in ["log_baseMean", "gc", "log_width"]:
        cov[c] = (cov[c] - cov[c].mean()) / cov[c].std()
    cov["log_baseMean_sq"] = cov["log_baseMean"] ** 2   # accessibility effect is not linear
    return np.column_stack([np.ones(len(df)), cov.to_numpy()])

#--------------------------------------------------
# Statistical models
#--------------------------------------------------

def fit_logistic(X, y, max_iter=100, tol=1e-8):
    """
    Fit a logistic regression by Newton-Raphson with step-halving.

    A small ridge penalty (RIDGE, not applied to the intercept) keeps the
    coefficients finite when a covariate (almost) separates hits from
    non-hits, which happens in small groups; with the sample sizes used here
    its effect on the estimates is negligible. Step-halving guarantees that
    every step improves the fit, so the algorithm cannot diverge.

    Parameters
    ----------
    X : np.ndarray
        Design matrix, first column = intercept.
    y : np.ndarray
        0/1 outcome.
    max_iter : int
        Maximum Newton iterations. Default 100.
    tol : float
        Relative change of the objective to stop. Default 1e-8.

    Returns
    -------
    beta : np.ndarray or None
        Coefficients (log odds ratios), or None if the fit fails.
    se : np.ndarray or None
        Standard errors, or None if the fit fails.
    """
    penalty = np.full(X.shape[1], RIDGE)
    penalty[0] = 0.0

    def objective(b):
        eta = X @ b
        return np.sum(np.logaddexp(0, eta) - y * eta) + 0.5 * np.sum(penalty * b ** 2)

    beta = np.zeros(X.shape[1])
    obj = objective(beta)
    try:
        for _ in range(max_iter):
            p = special.expit(X @ beta)
            hessian = X.T @ (X * (p * (1 - p))[:, None]) + np.diag(penalty)
            step = np.linalg.solve(hessian, X.T @ (y - p) - penalty * beta)
            t = 1.0
            while t > 1e-10:
                new_obj = objective(beta + t * step)
                if new_obj <= obj:
                    break
                t /= 2
            beta = beta + t * step
            converged = obj - new_obj < tol * (abs(obj) + 1)
            obj = new_obj
            if converged:
                break
        p = special.expit(X @ beta)
        hessian = X.T @ (X * (p * (1 - p))[:, None]) + np.diag(penalty)
        se = np.sqrt(np.diag(np.linalg.inv(hessian)))
    except np.linalg.LinAlgError:
        return None, None
    if not np.all(np.isfinite(beta)) or not np.all(np.isfinite(se)):
        return None, None
    return beta, se

def fit_ols(X, y):
    """
    Fit an ordinary least squares regression with robust standard errors.

    The standard errors are heteroscedasticity-robust (HC1), because the
    variance of the Wald statistic is not the same for every peak.

    Parameters
    ----------
    X : np.ndarray
        Design matrix, first column = intercept.
    y : np.ndarray
        Continuous outcome.

    Returns
    -------
    beta : np.ndarray
        Coefficients.
    se : np.ndarray
        HC1 standard errors.
    resid : np.ndarray
        Residuals (y - X @ beta).
    """
    xtx_inv = np.linalg.inv(X.T @ X)
    beta = xtx_inv @ (X.T @ y)
    resid = y - X @ beta
    n, k = X.shape
    meat = (X * resid[:, None]).T @ (X * resid[:, None])
    cov = xtx_inv @ meat @ xtx_inv * n / (n - k)
    return beta, np.sqrt(np.diag(cov)), resid

#--------------------------------------------------
# Analyses
#--------------------------------------------------

def gene_level_table(df, flag):
    """
    Collapse the peaks of one scope to one row per nearest gene.

    Parameters
    ----------
    df : pd.DataFrame
        Peaks of one scope (all or promoter) with a nearest gene.
    flag : str
        DA flag column: "is_up" or "is_down".

    Returns
    -------
    pd.DataFrame
        Indexed by gene: hit (>= 1 DA peak), n_peaks and the mean covariates
        of its peaks (log_baseMean, gc, log_width, proximal, distal).
    """
    return df.groupby("gene").agg(
        hit=(flag, "max"), n_peaks=(flag, "size"),
        log_baseMean=("log_baseMean", "mean"), gc=("gc", "mean"),
        log_width=("log_width", "mean"), proximal=("proximal", "mean"),
        distal=("distal", "mean"))

def expected_hit_probability(genes):
    """
    Estimate the probability that each gene is a hit by chance.

    Logistic regression of hit on the number of peaks of the gene and the
    covariates of its peaks, fitted on every gene of the universe. Genes with
    many weak or distal peaks are hit more often regardless of their function.
    Columns without variation (e.g. TSS category in the promoter scope) are
    dropped; if the full model fails, it falls back to the number of peaks
    only, and then to a constant rate.

    Parameters
    ----------
    genes : pd.DataFrame
        Output of gene_level_table.

    Returns
    -------
    np.ndarray
        P(hit) for each gene, in the order of `genes`.

    Raises
    ------
    RuntimeError
        If no model converges.
    """
    def z(x):
        return ((x - x.mean()) / x.std()).to_numpy()

    log_n = np.log(genes["n_peaks"])
    full = [log_n, genes["log_baseMean"], genes["gc"], genes["log_width"],
            genes["proximal"], genes["distal"]]
    full = [z(c) for c in full if c.std() > 0]
    if log_n.std() > 0 and genes["log_baseMean"].std() > 0:
        full += [z(log_n) ** 2, z(genes["log_baseMean"]) ** 2]
    simple = [z(log_n)] if log_n.std() > 0 else []

    y = genes["hit"].to_numpy(dtype=float)
    for cols in [full, simple, []]:
        X = np.column_stack([np.ones(len(genes))] + cols)
        beta, _ = fit_logistic(X, y)
        if beta is not None:
            if cols is not full:
                log_print("    WARNING: full gene-level model failed; using a simpler one")
            return special.expit(X @ beta)
    raise RuntimeError("Gene-level logistic model did not converge")

def functional_enrichment(df, subclass, collections, gene_index):
    """
    Run the covariate-adjusted, locus-level enrichment of one subclass.

    Runs for up and down peaks, with all peaks and with promoter peaks only.
    A locus is hit if any of its genes is hit, with
    P(locus hit) = 1 - prod(1 - P(gene hit)). Observed and expected hit loci
    of each gene set are compared with a one-sided normal approximation to the
    Poisson-binomial distribution (with continuity correction). Treating the
    genes of a locus as independent can only overestimate P(locus hit), so
    clustered loci make the test conservative, never liberal.

    Parameters
    ----------
    df : pd.DataFrame
        Output of load_da.
    subclass : str
        Subclass name.
    collections : dict
        {collection name: (gene_sets, M, L, S)} from read_gmt,
        membership_matrix and locus_matrices.
    gene_index : np.ndarray
        Gene symbols shared by every matrix.

    Returns
    -------
    pd.DataFrame
        One row per gene set tested, with the columns of FUNC_COLUMNS
        (padj is per subclass x direction x scope x collection).
    """
    rows = []
    for direction in ["up", "down"]:
        flag = f"is_{direction}"
        for scope in ["all", "promoter"]:
            in_scope = df["gene"].notna()
            if scope == "promoter":
                in_scope &= df["tss_category"] == "promoter"
            sub = df[in_scope]
            if sub[flag].sum() < MIN_DA_PEAKS:
                log_print(f"  functional {direction}/{scope}: skipped "
                          f"(< {MIN_DA_PEAKS} {direction} peaks)")
                continue

            genes = gene_level_table(sub, flag)
            try:
                genes["P"] = expected_hit_probability(genes)
            except RuntimeError as e:
                log_print(f"  functional {direction}/{scope}: skipped ({e})")
                continue
            genes = genes.reindex(gene_index)
            present = genes["hit"].notna().to_numpy(dtype=float)
            hit = genes["hit"].fillna(0).to_numpy(dtype=float)
            prob = genes["P"].fillna(0).to_numpy(dtype=float)
            log_miss = np.log1p(-np.clip(prob, 0, 1 - 1e-12))   # log P(gene not hit)
            log_print(f"  functional {direction}/{scope}: {int(hit.sum())} hit genes "
                      f"of {int(present.sum())} genes with tested peaks")

            for name, (gene_sets, M, L, S) in collections.items():
                locus_present = (L @ present > 0).astype(float)
                locus_hit = (L @ hit > 0).astype(float)
                locus_prob = -np.expm1(L @ log_miss)
                set_size = M @ present
                observed = S @ locus_hit
                expected = S @ locus_prob
                variance = S @ (locus_prob * (1 - locus_prob))
                keep = ((set_size >= MIN_SET_SIZE) & (set_size <= MAX_SET_SIZE)
                        & (variance > 0))
                if not keep.any():
                    continue
                zscore = (observed[keep] - 0.5 - expected[keep]) / np.sqrt(variance[keep])
                pval = stats.norm.sf(zscore)
                padj = benjamini_hochberg(pval)
                rows.append(pd.DataFrame({
                    "subclass": subclass, "direction": direction, "scope": scope,
                    "collection": name, "term": np.array(list(gene_sets))[keep],
                    "set_size": set_size[keep].astype(int),
                    "n_loci": (S @ locus_present)[keep].astype(int),
                    "n_genes_hit": (M @ hit)[keep].astype(int),
                    "n_loci_hit": observed[keep].astype(int),
                    "expected_loci_hit": expected[keep],
                    "fold_enrichment": observed[keep] / expected[keep],
                    "z": zscore, "pvalue": pval, "padj": padj,
                    "significant": (padj < ALPHA) & (observed[keep] >= MIN_LOCI_HIT),
                    "genes": "",
                    "_row": np.flatnonzero(keep),
                }))
    if not rows:
        return pd.DataFrame(columns=FUNC_COLUMNS)

    res = pd.concat(rows, ignore_index=True)
    # List the hit genes only for significant terms (sorted by number of DA peaks).
    for i in res.index[res["significant"]]:
        r = res.loc[i]
        M = collections[r["collection"]][1]
        sub = df[df["gene"].notna() & (df[f"is_{r['direction']}"] == 1)]
        if r["scope"] == "promoter":
            sub = sub[sub["tss_category"] == "promoter"]
        set_genes = gene_index[M[r["_row"]].indices]
        hits = sub.loc[sub["gene"].isin(set_genes), "gene"].value_counts()
        res.at[i, "genes"] = ";".join(hits.index[:N_GENES_REPORTED])
    return res[FUNC_COLUMNS]

def te_threshold_tests(df, subclass, features):
    """
    Test whether TE peaks are DA more or less often than peaks without TE.

    For each TE feature and direction with >= MIN_DA_PEAKS DA peaks: Fisher's
    exact test (crude) and logistic regression is_DA ~ TE + covariates
    (adjusted). Features with < MIN_TE_PEAKS peaks are skipped.

    Parameters
    ----------
    df : pd.DataFrame
        Output of load_da.
    subclass : str
        Subclass name.
    features : list of str
        TE columns to test (classes, any_TE, young_TE or families).

    Returns
    -------
    pd.DataFrame
        One row per direction x feature: counts, crude and adjusted log2 OR,
        95% CI of the adjusted log2 OR and raw p-values.
    """
    rows = []
    X_cov = covariate_matrix(df)
    for direction in ["up", "down"]:
        y = df[f"is_{direction}"].to_numpy()
        if y.sum() < MIN_DA_PEAKS:
            continue
        for feat in features:
            te = df[feat].to_numpy()
            if te.sum() < MIN_TE_PEAKS:
                continue
            a = int((te & y).sum())                   # TE and DA
            b = int(te.sum() - a)                     # TE, not DA
            c = int(y.sum() - a)                      # no TE, DA
            d = int(len(y) - a - b - c)               # no TE, not DA
            crude_or, crude_p = stats.fisher_exact([[a, b], [c, d]])

            beta, se = (None, None) if a == 0 else fit_logistic(
                np.column_stack([X_cov, te]), y.astype(float))
            if beta is None:
                adj_log_or = adj_se = adj_p = np.nan
            else:
                adj_log_or, adj_se = beta[-1], se[-1]
                adj_p = 2 * stats.norm.sf(abs(adj_log_or / adj_se))
            rows.append({
                "subclass": subclass, "direction": direction, "feature": feat,
                "n_te_peaks": int(te.sum()), "n_te_da": a, "n_da": int(y.sum()),
                "pct_da_te": 100 * a / (a + b), "pct_da_no_te": 100 * c / (c + d),
                "crude_log2OR": np.log2(crude_or) if 0 < crude_or < np.inf else np.nan,
                "crude_pvalue": crude_p,
                "adj_log2OR": adj_log_or / np.log(2),
                "adj_log2OR_ci_low": (adj_log_or - 1.96 * adj_se) / np.log(2),
                "adj_log2OR_ci_high": (adj_log_or + 1.96 * adj_se) / np.log(2),
                "adj_pvalue": adj_p,
            })
    return pd.DataFrame(rows, columns=TE_THRESHOLD_COLUMNS[:-3])

def te_shift_tests(df, subclass, features, level):
    """
    Test the shift of the Wald statistic in TE peaks (no padj threshold).

    OLS of the DESeq2 Wald statistic on the TE flag plus covariates, using
    every tested peak. A positive shift means that TE peaks gain accessibility
    with stage more than comparable peaks without TE.

    Parameters
    ----------
    df : pd.DataFrame
        Output of load_da.
    subclass : str
        Subclass name.
    features : list of str
        TE columns to test.
    level : str
        "class" or "family" (the BH correction is applied per level).

    Returns
    -------
    pd.DataFrame
        One row per feature: mean stat with and without TE, shift, its
        standard error and raw p-value.
    """
    rows = []
    X_cov = covariate_matrix(df)
    y = df["stat"].to_numpy(dtype=float)
    for feat in features:
        te = df[feat].to_numpy()
        if te.sum() < MIN_TE_PEAKS:
            continue
        beta, se, _ = fit_ols(np.column_stack([X_cov, te]), y)
        rows.append({
            "subclass": subclass, "level": level, "feature": feat,
            "n_te_peaks": int(te.sum()), "n_peaks": len(y),
            "mean_stat_te": y[te == 1].mean(), "mean_stat_no_te": y[te == 0].mean(),
            "shift": beta[-1], "shift_se": se[-1],
            "pvalue": 2 * stats.norm.sf(abs(beta[-1] / se[-1])),
        })
    return pd.DataFrame(rows, columns=TE_OLS_COLUMNS[:-2])

def te_age_tests(df, subclass):
    """
    Test whether younger TEs change differently from older ones.

    Spearman correlation between milliDiv and the Wald statistic adjusted for
    covariates (OLS residuals), within the peaks of each TE class. A negative
    rho means that younger (less diverged) TEs gain more accessibility.

    Parameters
    ----------
    df : pd.DataFrame
        Output of load_da.
    subclass : str
        Subclass name.

    Returns
    -------
    pd.DataFrame
        One row per TE class: number of TE peaks, Spearman rho and raw p-value.
    """
    rows = []
    _, _, resid = fit_ols(covariate_matrix(df), df["stat"].to_numpy(dtype=float))
    resid = pd.Series(resid, index=df.index)
    for te_class in TE_CLASSES[1:]:
        col = f"milliDiv_{te_class}"
        if col not in df.columns:
            continue
        mask = df[col].notna()
        if mask.sum() < MIN_TE_PEAKS:
            continue
        rho, p = stats.spearmanr(df.loc[mask, col], resid[mask])
        rows.append({"subclass": subclass, "te_class": te_class,
                     "n_te_peaks": int(mask.sum()), "spearman_rho": rho, "pvalue": p})
    return pd.DataFrame(rows, columns=TE_AGE_COLUMNS[:-2])

#--------------------------------------------------
# Graphs
#--------------------------------------------------

def graph_direction_dotplot(func_sig, direction):
    """
    Plot the significant terms of one direction.

    One panel per collection and one column per subclass of SUBCLASSES, with
    the top N_TERMS_DOTPLOT terms (scope "all") of each subclass. Dot size =
    hit loci, colour = -log10(padj), with the same scale in both panels. A
    collection without significant terms gets no panel.

    Parameters
    ----------
    func_sig : pd.DataFrame
        Significant functional terms (functional_enrichment_significant.csv).
    direction : str
        "up" or "down".
    """
    sub = func_sig[(func_sig["direction"] == direction) & (func_sig["scope"] == "all")
                   & func_sig["subclass"].isin(SUBCLASSES)]
    panels = []
    for collection in COLLECTIONS:
        col = sub[sub["collection"] == collection]
        if col.empty:
            log_print(f"  No significant {collection} terms for {direction} peaks.")
            continue
        top = col.sort_values("padj").groupby("subclass").head(N_TERMS_DOTPLOT)
        terms = list(top.sort_values("padj")["term"].unique())
        panels.append((collection, col[col["term"].isin(terms)], terms))
    if not panels:
        log_print(f"No significant terms for {direction} peaks; graph skipped.")
        return

    log_p = -np.log10(pd.concat([p[1] for p in panels])["padj"].clip(lower=1e-300))
    vmin, vmax = log_p.min(), log_p.max()
    if vmax - vmin < 1:   # one value (or almost): widen the scale so it stays readable
        vmin, vmax = max(0.0, vmin - 0.5), vmax + 0.5
    size = lambda n: np.clip(n, 0, 50) * 8
    x_pos = {sc: i for i, sc in enumerate(SUBCLASSES)}
    heights = [len(terms) + 1 for _, _, terms in panels]

    fig, axes = plt.subplots(len(panels), 1, squeeze=False, constrained_layout=True,
                             figsize=(9, 0.28 * sum(heights) + 2.5),
                             gridspec_kw={"height_ratios": heights})
    axes = axes[:, 0]
    for ax, (collection, data, terms) in zip(axes, panels):
        y = data["term"].map({t: i for i, t in enumerate(terms)})
        dots = ax.scatter(data["subclass"].map(x_pos), y, s=size(data["n_loci_hit"]),
                          c=-np.log10(data["padj"].clip(lower=1e-300)),
                          cmap="viridis", vmin=vmin, vmax=vmax)
        labels = [re.sub(r"^(HALLMARK|GOBP)_", "", t)[:MAX_TERM_LABEL] for t in terms]
        ax.set_yticks(range(len(terms)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xticks(range(len(SUBCLASSES)))
        ax.set_xticklabels(SUBCLASSES, rotation=45, ha="right", fontsize=8)
        ax.set_xlim(-0.5, len(SUBCLASSES) - 0.5)
        ax.set_ylim(len(terms) - 0.5, -0.5)
        ax.grid(color="lightgrey", lw=0.4)
        ax.set_axisbelow(True)
        ax.set_title(collection, loc="left", fontsize=9)

    fig.suptitle(f"Functional enrichment of DA {direction} peaks")
    cbar = fig.colorbar(dots, ax=list(axes), label="-log10(padj)", shrink=0.6)
    cbar.formatter.set_useOffset(False)   # plain tick labels, no "1e-6+3.03" offset
    cbar.update_ticks()
    handles = [axes[0].scatter([], [], s=size(n), color="grey") for n in (3, 10, 30)]
    axes[0].legend(handles, ["3", "10", "30"], title="hit loci", fontsize=7,
                   title_fontsize=7, loc="lower left", bbox_to_anchor=(1.0, 1.0),
                   frameon=False, ncol=3)

    out_path = GRAPHS_DIR / f"08_functional_dotplot_{direction}.png"
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    log_print(f"-> Saved: {out_path.name}")

def make_graphs():
    """
    Make the up and down dotplots from the saved table of significant terms.
    """
    for p in GRAPHS_DIR.glob("08_*.png"):
        p.unlink()
    log_print(f"\nGenerating graphs for: {', '.join(SUBCLASSES)}...")
    func_sig = pd.read_csv(FUNC_SIG_CSV)
    for direction in ["up", "down"]:
        graph_direction_dotplot(func_sig, direction)
    log_print(f"-> Graphs saved to: {GRAPHS_DIR.relative_to(PROJECT_ROOT)}")

def run_analyses():
    """
    Run the annotation, the statistics and save the tables.

    Steps 1-9 of the workflow: per-peak annotation, gene sets and loci,
    functional and TE tests for each subclass of SUBCLASSES, Benjamini-Hochberg
    correction and output tables in OUT_DIR.
    """
    log_print(f"Subclasses: {', '.join(SUBCLASSES)}")
    log_print(f"Parameters: CONTRAST={CONTRAST}, ALPHA={ALPHA}, MIN_DA_PEAKS={MIN_DA_PEAKS}, "
              f"MIN_TE_PEAKS={MIN_TE_PEAKS}, PROMOTER_DIST={PROMOTER_DIST}, "
              f"PROXIMAL_DIST={PROXIMAL_DIST}, set size {MIN_SET_SIZE}-{MAX_SET_SIZE}, "
              f"LOCUS_DIST={LOCUS_DIST}, MIN_LOCI_HIT={MIN_LOCI_HIT}, "
              f"BLAS threads={os.environ['OPENBLAS_NUM_THREADS']}")

    # --- Per-peak annotation (shared by all subclasses) ---
    te, families = build_te_table()
    coords = te.index.to_series().str.extract(r"^(chr[\w]+)[:\-](\d+)[:\-](\d+)$")
    peaks = pd.DataFrame({"chr": coords[0], "start": coords[1].astype(int),
                          "end": coords[2].astype(int)}, index=te.index)
    tss = read_tss()
    annotation = assign_nearest_tss(peaks, tss).join(compute_gc()).join(te)
    annotation.to_csv(ANNOT_CSV)
    log_print(f"-> Peak annotation saved to: {ANNOT_CSV.name}")

    # --- Gene set collections ---
    log_print("\nReading gene set collections (MSigDB) and grouping genes into loci...")
    gene_index = np.array(sorted(tss["gene"].unique()))
    # One position per gene (genes in both PARs keep their chrX copy)
    gene_pos = tss.drop_duplicates("gene").set_index("gene").loc[gene_index]
    collections = {}
    for name, path in COLLECTIONS.items():
        gene_sets = read_gmt(path)
        M = membership_matrix(gene_sets, gene_index)
        L, S = locus_matrices(gene_sets, gene_index, gene_pos)
        collections[name] = (gene_sets, M, L, S)
        n_genes = int((np.asarray(M.sum(axis=0)).ravel() > 0).sum())
        n_clustered = int((np.asarray(S.sum(axis=1)).ravel()
                           < np.asarray(M.sum(axis=1)).ravel()).sum())
        log_print(f"  {name}: {len(gene_sets)} sets; {n_genes} of {len(gene_index)} GTF "
                  f"genes found; {n_clustered} sets contain clustered genes")

    # --- Per-subclass analyses ---
    func_res, te_class_res, te_family_res, ols_res, age_res, summary = [], [], [], [], [], []
    for subclass in SUBCLASSES:
        df = load_da(da_path(subclass), annotation)
        n_up, n_down = int(df["is_up"].sum()), int(df["is_down"].sum())
        log_print(f"\n[{subclass}] tested peaks: {len(df)} | up: {n_up} | down: {n_down}")
        run_threshold = max(n_up, n_down) >= MIN_DA_PEAKS

        if run_threshold:
            func_res.append(functional_enrichment(df, subclass, collections, gene_index))
            te_class_res.append(te_threshold_tests(df, subclass, TE_CLASSES + ["young_TE"]))
            te_family_res.append(te_threshold_tests(df, subclass, families))
        else:
            log_print(f"  < {MIN_DA_PEAKS} DA peaks in both directions: "
                      "only threshold-free TE tests are run")
        ols_res.append(te_shift_tests(df, subclass, TE_CLASSES + ["young_TE"], "class"))
        ols_res.append(te_shift_tests(df, subclass, families, "family"))
        age_res.append(te_age_tests(df, subclass))
        summary.append({"subclass": subclass, "n_tested": len(df), "n_up": n_up,
                        "n_down": n_down, "threshold_tests": run_threshold})

    # --- Multiple testing correction and tables ---
    log_print("\nApplying Benjamini-Hochberg correction and saving tables...")
    func = concat_tables(func_res, FUNC_COLUMNS)
    func.to_csv(FUNC_CSV, index=False)
    func_sig = func[func["significant"].astype(bool)]
    func_sig.sort_values(["subclass", "direction", "scope", "collection", "padj"]) \
        .to_csv(FUNC_SIG_CSV, index=False)
    log_print(f"-> Functional enrichment: {len(func)} tests, {len(func_sig)} significant")

    for results, out in [(te_class_res, TE_CLASS_CSV), (te_family_res, TE_FAMILY_CSV)]:
        table = concat_tables(results, TE_THRESHOLD_COLUMNS)
        table["crude_padj"] = benjamini_hochberg(table["crude_pvalue"])
        table["adj_padj"] = benjamini_hochberg(table["adj_pvalue"])
        table["significant"] = table["adj_padj"] < ALPHA
        table[TE_THRESHOLD_COLUMNS].to_csv(out, index=False)
        log_print(f"-> {out.name}: {len(table)} tests, "
                  f"{int(table['significant'].sum())} significant after adjustment")

    te_ols = concat_tables(ols_res, TE_OLS_COLUMNS)
    te_ols["padj"] = np.nan
    for level in ["class", "family"]:
        m = te_ols["level"] == level
        te_ols.loc[m, "padj"] = benjamini_hochberg(te_ols.loc[m, "pvalue"])
    te_ols["significant"] = te_ols["padj"] < ALPHA
    te_ols[TE_OLS_COLUMNS].to_csv(TE_OLS_CSV, index=False)
    log_print(f"-> {TE_OLS_CSV.name}: {len(te_ols)} tests, "
              f"{int(te_ols['significant'].sum())} significant")

    te_age = concat_tables(age_res, TE_AGE_COLUMNS)
    te_age["padj"] = benjamini_hochberg(te_age["pvalue"])
    te_age["significant"] = te_age["padj"] < ALPHA
    te_age[TE_AGE_COLUMNS].to_csv(TE_AGE_CSV, index=False)
    log_print(f"-> {TE_AGE_CSV.name}: {len(te_age)} tests, "
              f"{int(te_age['significant'].sum())} significant")

    pd.DataFrame(summary).to_csv(SUMMARY_CSV, index=False)
    log_print(f"-> Summary saved to: {SUMMARY_CSV.name}")

#==================================================
# Main
#==================================================
def main():
    # Safe cleanup of specific script outputs to prevent leftover files
    # (peaks_nuc.txt is kept on purpose: bedtools nuc is slow)
    """
    Run the whole step 08: clean previous outputs, check inputs, run the
    analyses and make the graphs.
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    GRAPHS_DIR.mkdir(parents=True, exist_ok=True)
    for p in OUTPUTS:
        p.unlink(missing_ok=True)

    # Initialize fresh log file
    with open(LOG_PATH, 'w') as f:
        f.write("=== FUNCTIONAL AND TE ENRICHMENT PIPELINE START ===\n")

    try:
        check_dependencies()
        run_analyses()
        make_graphs()

    except Exception as e:
        log_print(f"\nENRICHMENT PIPELINE FAILED WITH ERROR: {str(e)}")
        raise

    log_print("\n=== FUNCTIONAL AND TE ENRICHMENT PIPELINE FINISHED ===")

if __name__ == "__main__":
    main()