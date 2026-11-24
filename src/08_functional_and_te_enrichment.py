#!/usr/bin/env python3
"""
08_functional_and_te_enrichment.py — Functional and TE enrichment of differentially accessible peaks
==========================================================================

Takes the differential accessibility (DA) results of step 06 (trend contrast: change in accessibility per step of AD neuropathological change) and asks two
questions for each subclass in SUBCLASSES (the Benjamini-Hochberg correction is applied over them):

  A. Functional enrichment: are the genes near the peaks that gain (up) or lose (down) accessibility over-represented in Hallmark or GO Biological Process gene sets?
  B. Transposable elements (TEs): do peaks whose center falls inside a TE (step 05) change more than comparable peaks without a TE?

Why covariates are needed: in step 06, weakly accessible peaks gain accessibility with stage much more often than strongly accessible ones (in
L2-3-IT, 91% of the DA peaks of the lowest baseMean quintile are up, against 63% in the highest quintile). TEs fall mostly in weak, distal peaks, so a crude
test could attribute this global effect to the TEs. Every test is therefore adjusted for accessibility (log baseMean and its square), GC content, peak
width and distance to the nearest TSS.

Methods: peaks are assigned to the nearest protein-coding TSS (GENCODE), GC comes from `bedtools nuc`, and per subclass the tested peaks (padj not NaN)
are the universe. A: locus-level, GOseq-like test (see functional_enrichment).
B: Fisher + adjusted logistic regression (DA vs TE), OLS shift of the Wald statistic (no threshold) and Spearman with TE divergence (milliDiv).

Requires having previously run src/00, src/01, src/05 and src/06.

Inputs:
  - results/06_differential_accessibility/da_{subclass}_trend.csv.gz   (step 06)
  - results/05_te_annotation/binary_matrix_peaks_te_{class,family}.csv (step 05)
  - results/05_te_annotation/peak_te_overlaps.csv                      (step 05)
  - results/01_metadata_output/atac_peaks.bed                          (step 01)
  - data/genome/hg38.fa                                                (step 00)
  - data/annotation/gencode.v50.basic.annotation.gtf.gz                (step 00)
  - data/genesets/{h.all, c5.go.bp}.v2026.2.Hs.symbols.gmt             (step 00)

Outputs (in results/08_enrichment/):
  - peaks_nuc.txt                          bedtools nuc output (kept between runs).
  - peak_annotation.csv.gz                 Nearest gene, TSS distance/category, GC, TE flags.
  - functional_enrichment.csv.gz           Every gene set tested.
  - functional_enrichment_significant.csv  padj < ALPHA and >= MIN_LOCI_HIT loci.
  - te_enrichment_{class,family}.csv       Fisher and adjusted logistic regression.
  - te_shift_ols.csv                       Shift of the Wald statistic.
  - te_age.csv                             Spearman stat vs milliDiv per TE class.
  - enrichment_summary.csv                 Peaks tested, up/down and analyses run.
  - 08_enrichment_log.txt
Every table has a `significant` column, used by the graphs.

Graphs (in results/08_graphs/): 08_functional_dotplot_{up,down}.png, one panel
per collection and one column per subclass, top N_TERMS_DOTPLOT terms (scope 'all'); dot size = hit loci, colour = -log10(padj).

Usage (from the repository root):
  python src/08_functional_and_te_enrichment.py
  nohup python -u src/08_functional_and_te_enrichment.py > logs/08_functional_and_te_enrichment.log 2>&1 &
Or import as a module and call main().
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

# Limit the BLAS / OpenMP threads of numpy and scipy (by default, one per core).
# Must be set before importing numpy; a value exported in the shell takes precedence.
N_THREADS = 4
for _v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ.setdefault(_v, str(N_THREADS))

import matplotlib
matplotlib.use("Agg")  # Headless environment for server
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, special, stats
#==================================================
# Variable definitions
#==================================================
ROOT = Path(__file__).resolve().parents[1]
DATA, RES = ROOT / "data", ROOT / "results"
DA_DIR, TE_DIR = RES / "06_differential_accessibility", RES / "05_te_annotation"
OUT_DIR, GRAPHS_DIR = RES / "08_enrichment", RES / "08_graphs"

# Inputs
CLASS_CSV = TE_DIR / "binary_matrix_peaks_te_class.csv"
FAMILY_CSV = TE_DIR / "binary_matrix_peaks_te_family.csv"
OVERLAPS_CSV = TE_DIR / "peak_te_overlaps.csv"
PEAKS_BED = RES / "01_metadata_output" / "atac_peaks.bed"
GENOME_FA = DATA / "genome" / "hg38.fa"
GTF_PATH = DATA / "annotation" / "gencode.v50.basic.annotation.gtf.gz"
MSIGDB_VERSION = "2026.2.Hs"
COLLECTIONS = {"Hallmark": DATA / "genesets" / f"h.all.v{MSIGDB_VERSION}.symbols.gmt", "GO_BP": DATA / "genesets" / f"c5.go.bp.v{MSIGDB_VERSION}.symbols.gmt",}

# Outputs
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
OUTPUTS = [ANNOT_CSV, FUNC_CSV, FUNC_SIG_CSV, TE_CLASS_CSV, TE_FAMILY_CSV,
           TE_OLS_CSV, TE_AGE_CSV, SUMMARY_CSV, LOG_PATH]

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
RIDGE = 1e-3                  # ridge penalty of the logistic models (see fit_logistic)
# Young (recently active) subfamilies; L1HS is listed even if few peaks reach it.
YOUNG_TE_RE = re.compile(r"^(L1HS|L1PA[2-8]|AluY[a-z0-9]*|HERVK-int|LTR5_Hs|SVA_[A-F])$")
COVARIATES = ["log_baseMean", "log_baseMean_sq", "gc", "log_width", "proximal", "distal"]

# Columns of the output tables (an empty table is still written with them)
FUNC_COLS = ["subclass", "direction", "scope", "collection", "term", "set_size","n_loci", "n_genes_hit", "n_loci_hit", "expected_loci_hit",
             "fold_enrichment", "z", "pvalue", "padj", "significant", "genes"]
TE_THR_COLS = ["subclass", "direction", "feature", "n_te_peaks", "n_te_da", "n_da", "pct_da_te", "pct_da_no_te", "crude_log2OR", "crude_pvalue",
               "adj_log2OR", "adj_log2OR_ci_low", "adj_log2OR_ci_high", "adj_pvalue", "crude_padj", "adj_padj", "significant"]
TE_OLS_COLS = ["subclass", "level", "feature", "n_te_peaks", "n_peaks", "mean_stat_te", "mean_stat_no_te", "shift", "shift_se", "pvalue", "padj", "significant"]
TE_AGE_COLS = ["subclass", "te_class", "n_te_peaks", "spearman_rho", "pvalue", "padj", "significant"]
#==================================================
# Functions
#==================================================
def log_print(msg):
    """Print a message and append it to the log file."""
    print(msg)
    with open(LOG_PATH, "a") as f:
        f.write(msg + "\n")

def da_path(sc):
    """DA table of step 06 for one subclass (e.g. "L2-3-IT")."""
    return DA_DIR / f"da_{sc}_{CONTRAST}.csv.gz"

def check_dependencies():
    """
    Check that every input exists before starting. bedtools, atac_peaks.bed and
    hg38.fa are only needed if peaks_nuc.txt is missing or empty.
    """
    inputs = [CLASS_CSV, FAMILY_CSV, OVERLAPS_CSV, GTF_PATH, *COLLECTIONS.values()]
    if not NUC_PATH.exists() or NUC_PATH.stat().st_size == 0:
        inputs += [PEAKS_BED, GENOME_FA]
        if not shutil.which("bedtools"):
            raise EnvironmentError("CRITICAL ERROR: 'bedtools' is not installed or not in PATH "
                                   "(needed for bedtools nuc). Install it with "
                                   "'conda install -c bioconda bedtools'.")
    missing = [str(p) for p in inputs if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing inputs:\n  " + "\n  ".join(missing))
    missing_da = [sc for sc in SUBCLASSES if not da_path(sc).exists()]
    if missing_da:
        raise FileNotFoundError(f"Missing DA results (step 06) in {DA_DIR} for: "
                                f"{', '.join(missing_da)}")

def benjamini_hochberg(pvalues):
    """Benjamini-Hochberg adjusted p-values, same order; NaN stays NaN."""
    p = np.asarray(pvalues, dtype=float)
    padj = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    n = ok.sum()
    if n:
        order = np.argsort(p[ok])
        ranked = np.minimum.accumulate((p[ok][order] * n / np.arange(1, n + 1))[::-1])[::-1]
        out = np.empty(n)
        out[order] = np.minimum(ranked, 1.0)
        padj[ok] = out
    return padj

def concat_tables(tables, cols):
    """Concatenate the non-empty tables; empty table with `cols` if none."""
    tables = [t for t in tables if not t.empty]
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame(columns=cols)

def save_table(tbl, cols, path, suffix=""):
    """Save a test table with `cols` and log how many tests are significant."""
    tbl[cols].to_csv(path, index=False)
    log_print(f"-> {path.name}: {len(tbl)} tests, {int(tbl['significant'].sum())} significant{suffix}")

#-----------------------------------------------------
# Annotation & Gene Sets

def read_tss():
    """
    One TSS per gene of GENE_TYPES on canonical chromosomes (GENCODE GTF).
    The GTF is 1-based: TSS = start on +, end on -, converted to 0-based.
    """
    log_print(f"\nReading gene TSSs from: {GTF_PATH.name}...")
    name_re, type_re = re.compile(r'gene_name "([^"]+)"'), re.compile(r'gene_type "([^"]+)"')
    rows = []
    with gzip.open(GTF_PATH, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                if fields[2] == "gene" and type_re.search(fields[8]).group(1) in GENE_TYPES:
                    tss = int(fields[3]) - 1 if fields[6] == "+" else int(fields[4]) - 1
                    rows.append((fields[0], tss, fields[6], name_re.search(fields[8]).group(1)))
    tss = pd.DataFrame(rows, columns=["chr", "tss", "strand", "gene"])
    tss = tss[tss["chr"].str.match(CANONICAL_CHR)].sort_values(["chr", "tss"]).reset_index(drop=True)
    log_print(f"Genes ({', '.join(GENE_TYPES)}) on canonical chromosomes: {tss['gene'].nunique()}")
    return tss

def assign_nearest_tss(peaks, tss):
    """
    Assign each peak center to its nearest TSS on the same chromosome:
    gene (NA if farther than MAX_GENE_DIST), dist_tss (signed, negative
    upstream following the gene strand) and tss_category.
    """
    log_print("\nAssigning peaks to the nearest TSS...")
    centers = (peaks["start"] + peaks["end"]) // 2   # same center as step 05
    out = pd.DataFrame({"gene": pd.NA, "dist_tss": np.nan}, index=peaks.index)
    for chrom, idx in peaks.groupby("chr").groups.items():
        ct = tss[tss["chr"] == chrom]
        if ct.empty:
            continue
        pos, c = ct["tss"].to_numpy(), centers.loc[idx].to_numpy()
        right = np.clip(np.searchsorted(pos, c), 0, len(pos) - 1)
        left = np.clip(right - 1, 0, len(pos) - 1)
        nearest = np.where(np.abs(c - pos[left]) <= np.abs(c - pos[right]), left, right)
        dist = c - pos[nearest]
        out.loc[idx, "gene"] = ct["gene"].to_numpy()[nearest]
        out.loc[idx, "dist_tss"] = np.where(ct["strand"].to_numpy()[nearest] == "+", dist, -dist)
    abs_d = out["dist_tss"].abs()
    out["tss_category"] = np.select([abs_d <= PROMOTER_DIST, abs_d <= PROXIMAL_DIST],
                                    ["promoter", "proximal"], default="distal")
    out.loc[abs_d > MAX_GENE_DIST, "gene"] = pd.NA
    log_print(f"Peaks per TSS category: {out['tss_category'].value_counts().to_dict()}")
    log_print(f"Peaks without a gene within {MAX_GENE_DIST:,} bp: {out['gene'].isna().sum()}")
    return out

def compute_gc():
    """GC fraction per peak (Series "gc" indexed by peak_id), running `bedtools nuc` if needed."""
    if not NUC_PATH.exists() or NUC_PATH.stat().st_size == 0:
        log_print(f"\nRunning bedtools nuc on {PEAKS_BED.name} (this may take a while)...")
        try:
            with open(NUC_PATH, "w") as f:
                subprocess.run(["bedtools", "nuc", "-fi", str(GENOME_FA), "-bed", str(PEAKS_BED)],
                               stdout=f, check=True)
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
    Per-peak TE table from step 05: any_TE, 0/1 columns per TE class and family, young_TE and milliDiv_{class}. 
    Also returns the family names.
    """
    log_print("\nReading TE annotation (step 05)...")
    te_class = pd.read_csv(CLASS_CSV, index_col="peak_id")
    te_family = pd.read_csv(FAMILY_CSV, index_col="peak_id").drop(columns="any_TE")
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

def collection_matrices(gene_sets, gene_index, gene_pos):
    """
    Sparse 0/1 matrices of one collection (genes not in gene_index are ignored):
      M (sets x genes): membership.
      L (loci x genes) and S (sets x loci): genes of the same set and chromosome
        whose TSSs are < LOCUS_DIST apart (chained) form one locus, so a gene
        cluster (olfactory receptors, histones) counts once.
    """
    pos = {g: i for i, g in enumerate(gene_index)}
    chrom, tss = gene_pos["chr"].to_numpy(), gene_pos["tss"].to_numpy()
    m_r, m_c, l_r, l_c, s_r, s_c = [], [], [], [], [], []
    n_loci = 0
    for i, genes in enumerate(gene_sets.values()):
        idx = np.array([pos[g] for g in genes if g in pos], dtype=int)
        if not idx.size:
            continue
        m_r.extend([i] * idx.size)
        m_c.extend(idx)
        idx = idx[np.lexsort((tss[idx], chrom[idx]))]
        new_loc = np.ones(idx.size, dtype=bool)
        new_loc[1:] = (chrom[idx][1:] != chrom[idx][:-1]) | (np.diff(tss[idx]) >= LOCUS_DIST)
        locus = n_loci + np.cumsum(new_loc) - 1
        l_r.extend(locus)
        l_c.extend(idx)
        s_r.extend([i] * int(new_loc.sum()))
        s_c.extend(np.unique(locus))
        n_loci = int(locus[-1]) + 1
    n_sets, n_genes = len(gene_sets), len(gene_index)
    M = sparse.csr_matrix((np.ones(len(m_r)), (m_r, m_c)), shape=(n_sets, n_genes))
    L = sparse.csr_matrix((np.ones(len(l_r)), (l_r, l_c)), shape=(n_loci, n_genes))
    S = sparse.csr_matrix((np.ones(len(s_r)), (s_r, s_c)), shape=(n_sets, n_loci))
    return M, L, S

#-----------------------------------------------------
# Data Loading & Statistical Models

def load_da(path, annotation):
    """
    Peak table of one subclass: tested peaks only (padj not NaN), is_up / is_down flags (padj < ALPHA), the per-peak annotation and the covariates.
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
    Design matrix shared by the TE models: intercept, z-scored continuous covariates and the square of log_baseMean (its effect is not linear).
    """
    cov = df[COVARIATES].astype(float).copy()
    for c in ["log_baseMean", "gc", "log_width"]:
        cov[c] = (cov[c] - cov[c].mean()) / cov[c].std()
    cov["log_baseMean_sq"] = cov["log_baseMean"] ** 2
    return np.column_stack([np.ones(len(df)), cov.to_numpy()])

def fit_logistic(X, y, max_iter=100, tol=1e-8):
    """
    Logistic regression by Newton-Raphson with step-halving (X includes the intercept). A small ridge penalty (RIDGE, not on the intercept) keeps the
    coefficients finite under (quasi-)separation in small groups; step-halving makes every step improve the fit, so it cannot diverge.
    Returns (beta, se) as log odds ratios, or (None, None) if the fit fails.
    """
    pen = np.array([0.0] + [RIDGE] * (X.shape[1] - 1))

    def obj_fn(b):
        eta = X @ b
        return np.sum(np.logaddexp(0, eta) - y * eta) + 0.5 * np.sum(pen * b ** 2)

    beta = np.zeros(X.shape[1])
    obj = obj_fn(beta)
    try:
        for _ in range(max_iter):
            p = special.expit(X @ beta)
            H = X.T @ (X * (p * (1 - p))[:, None]) + np.diag(pen)
            step = np.linalg.solve(H, X.T @ (y - p) - pen * beta)
            t = 1.0
            while t > 1e-10:
                new_obj = obj_fn(beta + t * step)
                if new_obj <= obj:
                    break
                t /= 2
            beta = beta + t * step
            converged = obj - new_obj < tol * (abs(obj) + 1)
            obj = new_obj
            if converged:
                break
        p = special.expit(X @ beta)
        H = X.T @ (X * (p * (1 - p))[:, None]) + np.diag(pen)
        se = np.sqrt(np.diag(np.linalg.inv(H)))
    except np.linalg.LinAlgError:
        return None, None
    return (beta, se) if np.all(np.isfinite(beta)) and np.all(np.isfinite(se)) else (None, None)

def fit_ols(X, y):
    """
    OLS with heteroscedasticity-robust (HC1) standard errors, because the variance of the Wald statistic differs between peaks. Returns (beta, se, resid).
    """
    xtx_inv = np.linalg.inv(X.T @ X)
    beta = xtx_inv @ (X.T @ y)
    resid = y - X @ beta
    n, k = X.shape
    meat = (X * resid[:, None]).T @ (X * resid[:, None])
    cov = xtx_inv @ meat @ xtx_inv * n / (n - k)
    return beta, np.sqrt(np.diag(cov)), resid

#--------------------------------------------------
# Enrichment & TE Analyses

def expected_hit_probability(genes):
    """
    P(gene is a hit by chance): logistic regression of hit on the number of peaks of the gene and the mean covariates of its peaks (genes with many
    weak or distal peaks are hit more often regardless of function). Columns without variation are dropped; if the model fails, falls back to the
    number of peaks only, then to a constant rate.
    """
    z = lambda x: ((x - x.mean()) / x.std()).to_numpy()
    log_n = np.log(genes["n_peaks"])
    full = [z(c) for c in [log_n, genes["log_baseMean"], genes["gc"], genes["log_width"],
                           genes["proximal"], genes["distal"]] if c.std() > 0]
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
    Covariate-adjusted, locus-level enrichment of one subclass (up / down peaks, all / promoter peaks). A locus is hit if any of its genes is hit,
    with P(locus hit) = 1 - prod(1 - P(gene hit)); observed vs expected hit loci are compared with a one-sided, continuity-corrected normal
    approximation to the Poisson-binomial. Treating the genes of a locus as independent can only overestimate P(locus hit), so the test is
    conservative. padj is per subclass x direction x scope x collection. Discarded alternatives: a peak-level hypergeometric test gave ~16% of
    p < 0.05 with random gene sets, and counting genes instead of loci made single gene clusters (olfactory receptors, histones) look like pathways.
    """
    rows = []
    for direction in ["up", "down"]:
        flag = f"is_{direction}"
        for scope in ["all", "promoter"]:
            sub = df[df["gene"].notna()]
            if scope == "promoter":
                sub = sub[sub["tss_category"] == "promoter"]
            if sub[flag].sum() < MIN_DA_PEAKS:
                log_print(f"  functional {direction}/{scope}: skipped "
                          f"(< {MIN_DA_PEAKS} {direction} peaks)")
                continue
            # One row per nearest gene: hit (>= 1 DA peak), n_peaks and mean covariates
            genes = sub.groupby("gene").agg(
                hit=(flag, "max"), n_peaks=(flag, "size"),
                **{c: (c, "mean") for c in ["log_baseMean", "gc", "log_width", "proximal", "distal"]}
            )
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
                locus_prob = -np.expm1(L @ log_miss)
                set_size = M @ present
                obs = S @ (L @ hit > 0).astype(float)
                expected = S @ locus_prob
                var = S @ (locus_prob * (1 - locus_prob))
                keep = (set_size >= MIN_SET_SIZE) & (set_size <= MAX_SET_SIZE) & (var > 0)
                if not keep.any():
                    continue
                zscore = (obs[keep] - 0.5 - expected[keep]) / np.sqrt(var[keep])
                pval = stats.norm.sf(zscore)
                padj = benjamini_hochberg(pval)
                rows.append(pd.DataFrame({
                    "subclass": subclass, "direction": direction, "scope": scope,
                    "collection": name, "term": np.array(list(gene_sets))[keep],
                    "set_size": set_size[keep].astype(int),
                    "n_loci": (S @ (L @ present > 0).astype(float))[keep].astype(int),
                    "n_genes_hit": (M @ hit)[keep].astype(int),
                    "n_loci_hit": obs[keep].astype(int),
                    "expected_loci_hit": expected[keep],
                    "fold_enrichment": obs[keep] / expected[keep],
                    "z": zscore, "pvalue": pval, "padj": padj,
                    "significant": (padj < ALPHA) & (obs[keep] >= MIN_LOCI_HIT),
                    "genes": "", "_row": np.flatnonzero(keep),
                }))
    if not rows:
        return pd.DataFrame(columns=FUNC_COLS)
    res = pd.concat(rows, ignore_index=True)
    # List the hit genes only for significant terms (sorted by number of DA peaks).
    for i in res.index[res["significant"]]:
        r = res.loc[i]
        M = collections[r["collection"]][1]
        sub = df[df["gene"].notna() & (df[f"is_{r['direction']}"] == 1)]
        if r["scope"] == "promoter":
            sub = sub[sub["tss_category"] == "promoter"]
        hits = sub.loc[sub["gene"].isin(gene_index[M[r["_row"]].indices]), "gene"].value_counts()
        res.at[i, "genes"] = ";".join(hits.index[:N_GENES_REPORTED])
    return res[FUNC_COLS]

def te_threshold_tests(df, subclass, features):
    """
    Are TE peaks DA more or less often than peaks without TE? Per direction and feature: Fisher's exact test (crude) and logistic regression
    is_DA ~ TE + covariates (adjusted log2 OR with 95% CI). Raw p-values.
    """
    rows, X_cov = [], covariate_matrix(df)
    for direction in ["up", "down"]:
        y = df[f"is_{direction}"].to_numpy()
        if y.sum() < MIN_DA_PEAKS:
            continue
        for feat in features:
            te = df[feat].to_numpy()
            if te.sum() < MIN_TE_PEAKS:
                continue
            a = int((te & y).sum())        # TE and DA
            b = int(te.sum() - a)          # TE, not DA
            c = int(y.sum() - a)           # no TE, DA
            d = int(len(y) - a - b - c)    # no TE, not DA
            crude_or, crude_p = stats.fisher_exact([[a, b], [c, d]])
            beta, se = (None, None) if a == 0 else fit_logistic(
                np.column_stack([X_cov, te]), y.astype(float))
            adj_or, adj_se = (np.nan, np.nan) if beta is None else (beta[-1], se[-1])
            rows.append({
                "subclass": subclass, "direction": direction, "feature": feat,
                "n_te_peaks": int(te.sum()), "n_te_da": a, "n_da": int(y.sum()),
                "pct_da_te": 100 * a / (a + b), "pct_da_no_te": 100 * c / (c + d),
                "crude_log2OR": np.log2(crude_or) if 0 < crude_or < np.inf else np.nan,
                "crude_pvalue": crude_p,
                "adj_log2OR": adj_or / np.log(2),
                "adj_log2OR_ci_low": (adj_or - 1.96 * adj_se) / np.log(2),
                "adj_log2OR_ci_high": (adj_or + 1.96 * adj_se) / np.log(2),
                "adj_pvalue": np.nan if beta is None else 2 * stats.norm.sf(abs(adj_or / adj_se)),
            })
    return pd.DataFrame(rows, columns=TE_THR_COLS[:-3])

def te_shift_tests(df, subclass, features, level):
    """
    Shift of the Wald statistic in TE peaks, without padj threshold: OLS of stat ~ TE + covariates. A positive shift means TE peaks gain accessibility
    with stage more than comparable peaks without TE. `level` groups the BH.
    """
    rows, X_cov, y = [], covariate_matrix(df), df["stat"].to_numpy(dtype=float)
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
    return pd.DataFrame(rows, columns=TE_OLS_COLS[:-2])

def te_age_tests(df, subclass):
    """
    Do younger TEs change differently? Spearman correlation between milliDiv and the covariate-adjusted Wald statistic (OLS residuals), per TE class.
    Negative rho = younger (less diverged) TEs gain more accessibility.
    """
    _, _, resid = fit_ols(covariate_matrix(df), df["stat"].to_numpy(dtype=float))
    resid = pd.Series(resid, index=df.index)
    rows = []
    for te_class in TE_CLASSES[1:]:
        col = f"milliDiv_{te_class}"
        if col in df.columns and (mask := df[col].notna()).sum() >= MIN_TE_PEAKS:
            rho, p = stats.spearmanr(df.loc[mask, col], resid[mask])
            rows.append({"subclass": subclass, "te_class": te_class,
                         "n_te_peaks": int(mask.sum()), "spearman_rho": rho, "pvalue": p})
    return pd.DataFrame(rows, columns=TE_AGE_COLS[:-2])

#-------------------------------------------------
# Graphs & Orchestration

def graph_direction_dotplot(func_sig, direction):
    """
    Dotplot of the significant terms of one direction: one panel per collection, one column per subclass, top N_TERMS_DOTPLOT terms (scope "all") per subclass.
    Size = hit loci, colour = -log10(padj), shared scale.
    """
    sub = func_sig[(func_sig["direction"] == direction) & (func_sig["scope"] == "all")
                   & func_sig["subclass"].isin(SUBCLASSES)]
    panels = []
    for col_name in COLLECTIONS:
        col = sub[sub["collection"] == col_name]
        if col.empty:
            log_print(f"  No significant {col_name} terms for {direction} peaks.")
            continue
        top = col.sort_values("padj").groupby("subclass").head(N_TERMS_DOTPLOT)
        terms = list(top.sort_values("padj")["term"].unique())
        panels.append((col_name, col[col["term"].isin(terms)], terms))
    if not panels:
        log_print(f"No significant terms for {direction} peaks; graph skipped.")
        return

    log_p = -np.log10(pd.concat([p[1] for p in panels])["padj"].clip(lower=1e-300))
    vmin, vmax = log_p.min(), log_p.max()
    if vmax - vmin < 1:   # one value (or almost): widen the scale so it stays readable
        vmin, vmax = max(0.0, vmin - 0.5), vmax + 0.5
    size = lambda n: np.clip(n, 0, 50) * 8
    x_pos = {sc: i for i, sc in enumerate(SUBCLASSES)}
    heights = [len(t) + 1 for _, _, t in panels]

    fig, axes = plt.subplots(len(panels), 1, squeeze=False, constrained_layout=True,
                             figsize=(9, 0.28 * sum(heights) + 2.5),
                             gridspec_kw={"height_ratios": heights})
    axes = axes[:, 0]
    for ax, (col_name, data, terms) in zip(axes, panels):
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
        ax.set_title(col_name, loc="left", fontsize=9)

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
    """Up and down dotplots from the saved table of significant terms."""
    for p in GRAPHS_DIR.glob("08_*.png"):
        p.unlink()
    log_print(f"\nGenerating graphs for: {', '.join(SUBCLASSES)}...")
    func_sig = pd.read_csv(FUNC_SIG_CSV)
    for direction in ["up", "down"]:
        graph_direction_dotplot(func_sig, direction)
    log_print(f"-> Graphs saved to: {GRAPHS_DIR.relative_to(ROOT)}")

def run_analyses():
    """Workflow steps 1-9: annotation, per-subclass tests, BH correction and tables."""
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

    # --- Gene set collections (GMT: one set per line: name, URL, genes...) ---
    log_print("\nReading gene set collections (MSigDB) and grouping genes into loci...")
    gene_index = np.array(sorted(tss["gene"].unique()))
    # One position per gene (genes in both PARs keep their chrX copy)
    gene_pos = tss.drop_duplicates("gene").set_index("gene").loc[gene_index]
    collections = {}
    for name, path in COLLECTIONS.items():
        with open(path) as f:
            gene_sets = {p[0]: set(p[2:]) for p in (line.rstrip("\n").split("\t") for line in f)}
        M, L, S = collection_matrices(gene_sets, gene_index, gene_pos)
        collections[name] = (gene_sets, M, L, S)
        n_genes = int((np.asarray(M.sum(axis=0)).ravel() > 0).sum())
        n_clust = int((np.asarray(S.sum(axis=1)).ravel() < np.asarray(M.sum(axis=1)).ravel()).sum())
        log_print(f"  {name}: {len(gene_sets)} sets; {n_genes} of {len(gene_index)} GTF "
                  f"genes found; {n_clust} sets contain clustered genes")

    # --- Per-subclass analyses ---
    func_res, te_cls_res, te_fam_res, ols_res, age_res, summary = [], [], [], [], [], []
    for sc in SUBCLASSES:
        df = load_da(da_path(sc), annotation)
        n_up, n_down = int(df["is_up"].sum()), int(df["is_down"].sum())
        log_print(f"\n[{sc}] tested peaks: {len(df)} | up: {n_up} | down: {n_down}")
        run_thr = max(n_up, n_down) >= MIN_DA_PEAKS
        if run_thr:
            func_res.append(functional_enrichment(df, sc, collections, gene_index))
            te_cls_res.append(te_threshold_tests(df, sc, TE_CLASSES + ["young_TE"]))
            te_fam_res.append(te_threshold_tests(df, sc, families))
        else:
            log_print(f"  < {MIN_DA_PEAKS} DA peaks in both directions: "
                      "only threshold-free TE tests are run")
        ols_res.append(te_shift_tests(df, sc, TE_CLASSES + ["young_TE"], "class"))
        ols_res.append(te_shift_tests(df, sc, families, "family"))
        age_res.append(te_age_tests(df, sc))
        summary.append({"subclass": sc, "n_tested": len(df), "n_up": n_up,
                        "n_down": n_down, "threshold_tests": run_thr})

    # --- Multiple testing correction and tables ---
    log_print("\nApplying Benjamini-Hochberg correction and saving tables...")
    func = concat_tables(func_res, FUNC_COLS)
    func.to_csv(FUNC_CSV, index=False)
    func_sig = func[func["significant"].astype(bool)]
    func_sig.sort_values(["subclass", "direction", "scope", "collection", "padj"]) \
        .to_csv(FUNC_SIG_CSV, index=False)
    log_print(f"-> Functional enrichment: {len(func)} tests, {len(func_sig)} significant")

    for res_list, out in [(te_cls_res, TE_CLASS_CSV), (te_fam_res, TE_FAMILY_CSV)]:
        tbl = concat_tables(res_list, TE_THR_COLS)
        tbl["crude_padj"] = benjamini_hochberg(tbl["crude_pvalue"])
        tbl["adj_padj"] = benjamini_hochberg(tbl["adj_pvalue"])
        tbl["significant"] = tbl["adj_padj"] < ALPHA
        save_table(tbl, TE_THR_COLS, out, " after adjustment")

    te_ols = concat_tables(ols_res, TE_OLS_COLS)
    te_ols["padj"] = np.nan
    for level in ["class", "family"]:
        m = te_ols["level"] == level
        te_ols.loc[m, "padj"] = benjamini_hochberg(te_ols.loc[m, "pvalue"])
    te_ols["significant"] = te_ols["padj"] < ALPHA
    save_table(te_ols, TE_OLS_COLS, TE_OLS_CSV)

    te_age = concat_tables(age_res, TE_AGE_COLS)
    te_age["padj"] = benjamini_hochberg(te_age["pvalue"])
    te_age["significant"] = te_age["padj"] < ALPHA
    save_table(te_age, TE_AGE_COLS, TE_AGE_CSV)

    pd.DataFrame(summary).to_csv(SUMMARY_CSV, index=False)
    log_print(f"-> Summary saved to: {SUMMARY_CSV.name}")
#==================================================
# Main
#==================================================
def main():
    """
    Run step 08: clean previous outputs (peaks_nuc.txt is kept on purpose, bedtools nuc is slow), check inputs, run the analyses and make the graphs.
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    GRAPHS_DIR.mkdir(parents=True, exist_ok=True)
    for p in OUTPUTS:
        p.unlink(missing_ok=True)
    with open(LOG_PATH, "w") as f:
        f.write("=== FUNCTIONAL AND TE ENRICHMENT PIPELINE START ===\n")
    try:
        check_dependencies()
        run_analyses()
        make_graphs()
    except Exception as e:
        log_print(f"\nENRICHMENT PIPELINE FAILED WITH ERROR: {e}")
        raise
    log_print("\n=== FUNCTIONAL AND TE ENRICHMENT PIPELINE FINISHED ===")

if __name__ == "__main__":
    main()