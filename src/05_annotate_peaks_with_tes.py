#!/usr/bin/env python3
"""
05_annotate_peaks_with_tes.py — Peak Annotation with Transposable Elements
==========================================================================

Takes the accessibility peak coordinates (atac_peaks.bed) generated in step 01 
and intersects them with the RepeatMasker database (rmsk.txt.gz) from step 00.

Workflow:
  1. Reads and formats rmsk.txt.gz (handling UCSC headers) into an 8-column BED.
  2. Filters by exact TE classes (LINE, SINE, LTR, DNA, Retroposon) without 
     altering ambiguous names, and restricts to canonical chromosomes.
  3. Calculates the geometric center of each peak and strictly sorts them to 
     avoid length biases and bedtools errors.
  4. Runs `bedtools intersect -sorted` mapping only the peak center.
  5. Generates a detailed overlap table (using original peak coordinates) and 
     two complete binary matrices (Peaks x Class, Peaks x Family).

Inputs:
  - data/rmsk.txt.gz
  - results/01_metadata_output/atac_peaks.bed

Outputs (in results/05_te_annotation/):
  - atac_peaks_centers.bed
  - rmsk_filtered.bed
  - raw_intersect.bed
  - peak_te_overlaps.csv
  - binary_matrix_peaks_te_class.csv
  - binary_matrix_peaks_te_family.csv
  - 05_te_annotation_log.txt
  
Usage:
  Execute directly from the command line from the repository root:
  python src/05_annotate_peaks_with_tes.py
  
  Or import as a module and call the main() function.
"""

import gzip
import shutil
import subprocess
import pandas as pd
from pathlib import Path

#==================================================
# Variable definitions
#==================================================
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"

# Inputs
RMSK_PATH = DATA_DIR / "rmsk.txt.gz"
PEAKS_BED_PATH = RESULTS_DIR / "01_metadata_output" / "atac_peaks.bed"

# Outputs
OUT_DIR = RESULTS_DIR / "05_te_annotation"
CENTERS_BED_PATH = OUT_DIR / "atac_peaks_centers.bed"
RMSK_BED_PATH = OUT_DIR / "rmsk_filtered.bed"
INTERSECT_OUT = OUT_DIR / "raw_intersect.bed"
DETAILS_CSV = OUT_DIR / "peak_te_overlaps.csv"
CLASS_MATRIX_CSV = OUT_DIR / "binary_matrix_peaks_te_class.csv"
FAMILY_MATRIX_CSV = OUT_DIR / "binary_matrix_peaks_te_family.csv"
LOG_PATH = OUT_DIR / "05_te_annotation_log.txt"

OUTPUTS = [
    CENTERS_BED_PATH, RMSK_BED_PATH, INTERSECT_OUT, DETAILS_CSV,
    CLASS_MATRIX_CSV, FAMILY_MATRIX_CSV, LOG_PATH
]

# Exact classes as defined in the project report
TE_CLASSES = ['LINE', 'SINE', 'LTR', 'DNA', 'Retroposon']

#==================================================
# Functions
#==================================================
def log_print(msg):
    """Prints a message and instantly appends it to the log file."""
    print(msg)
    with open(LOG_PATH, 'a') as f:
        f.write(msg + "\n")

def check_dependencies():
    """Ensures bedtools is available and inputs exist."""
    if not shutil.which("bedtools"):
        raise EnvironmentError(
            "CRITICAL ERROR: 'bedtools' is not installed or not in PATH. "
            "Please add the 'bioconda' channel and install it: "
            "'conda install -c bioconda bedtools', or update environment.yml."
        )
    if not PEAKS_BED_PATH.exists():
        raise FileNotFoundError(f"Peak BED not found: {PEAKS_BED_PATH}")
    if not RMSK_PATH.exists():
        raise FileNotFoundError(f"RepeatMasker file not found: {RMSK_PATH}")

def format_rmsk():
    """
    Reads UCSC rmsk dump, detecting headers dynamically, filters it, 
    and converts it to an 8-column BED keeping milliDiv for young TE divergence.
    """
    log_print(f"\nReading raw RepeatMasker annotation: {RMSK_PATH.name}...")
    
    rmsk_cols = [
        'bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns', 'genoName', 
        'genoStart', 'genoEnd', 'genoLeft', 'strand', 'repName', 'repClass', 
        'repFamily', 'repStart', 'repEnd', 'repLeft', 'id'
    ]
    
    with gzip.open(RMSK_PATH, "rt") as f:
        has_header = f.readline().startswith("#bin")
        
    df = pd.read_csv(RMSK_PATH, sep='\t', header=0 if has_header else None, 
                     names=rmsk_cols, low_memory=False)
    log_print(f"Total raw records: {len(df)}")
    
    df = df[df['genoName'].str.match(r'^chr([\d]+|[XYM])$')]
    df = df[df['repClass'].isin(TE_CLASSES)]
    log_print(f"Records after chromosome and strict class filtering: {len(df)}")
    
    # Prepend the class to the family to avoid naming collisions across classes
    df['repFamily'] = df['repClass'] + "/" + df['repFamily']
    
    bed_rmsk = df[['genoName', 'genoStart', 'genoEnd', 'repName', 
                   'milliDiv', 'strand', 'repFamily', 'repClass']]
    
    bed_rmsk = bed_rmsk.sort_values(['genoName', 'genoStart', 'genoEnd'])
    
    bed_rmsk.to_csv(RMSK_BED_PATH, sep='\t', header=False, index=False)
    log_print(f"-> Filtered TE BED saved to: {RMSK_BED_PATH.name}")

def extract_peak_centers():
    """Calculates the geometric center of each peak and sorts them strictly."""
    log_print("\nExtracting sorted peak centers (1 bp) to prevent size bias...")
    peaks = pd.read_csv(PEAKS_BED_PATH, sep="\t", header=None, names=["chr", "start", "end", "peak_id"])
    mid = (peaks["start"] + peaks["end"]) // 2
    
    # Assignment and sort to prevent out-of-order records due to overlapping varying-length peaks
    centers = peaks.assign(start=mid, end=mid + 1).sort_values(["chr", "start", "end"])
    centers.to_csv(CENTERS_BED_PATH, sep="\t", header=False, index=False)
    
    log_print(f"-> Peak centers BED saved to: {CENTERS_BED_PATH.name}")

def run_intersect():
    """Runs bedtools intersect mapping the 1bp peak center against the TEs."""
    log_print("\nRunning bedtools intersect...")
    log_print(f"  - File A (Peak Centers): {CENTERS_BED_PATH.name}")
    log_print(f"  - File B (TEs):          {RMSK_BED_PATH.name}")
    log_print("  - Criterion:             Geometric center of the peak (1 bp overlap)")
    
    cmd = [
        "bedtools", "intersect", 
        "-a", str(CENTERS_BED_PATH), 
        "-b", str(RMSK_BED_PATH),
        "-sorted",
        "-wo"
    ]
    
    try:
        with open(INTERSECT_OUT, "w") as out_file:
            subprocess.run(cmd, stdout=out_file, check=True)
        log_print(f"-> Intersection completed. Raw output: {INTERSECT_OUT.name}")
    except subprocess.CalledProcessError as e:
        log_print(f"ERROR executing bedtools: {e}")
        raise

def process_intersection():
    """Processes the intersection to create full binary matrices and detailed tables."""
    log_print("\nProcessing intersection results...")
    
    cols = [
        'chr_center', 'start_center', 'end_center', 'peak_id',
        'chr_te', 'start_te', 'end_te', 'repName', 'milliDiv', 'strand_te', 
        'repFamily', 'repClass', 'overlap_bp'
    ]
    
    try:
        df_int = pd.read_csv(INTERSECT_OUT, sep='\t', header=None, names=cols)
    except pd.errors.EmptyDataError:
        df_int = pd.DataFrame(columns=cols)

    # Load ALL peaks to ensure matrices include peaks with 0 TEs and to restore original coordinates
    all_peaks_df = pd.read_csv(PEAKS_BED_PATH, sep="\t", header=None, names=['chr_p', 'start_p', 'end_p', 'peak_id'])
    all_peaks_series = all_peaks_df['peak_id'].rename("peak_id")
    
    if df_int.empty:
        log_print("WARNING: The intersection is empty. No peaks overlap with TEs given the current criteria.")
        df_details = pd.DataFrame(columns=['peak_id', 'chr_p', 'start_p', 'end_p', 'repName', 
                                           'repFamily', 'repClass', 'milliDiv', 'strand_te'])
        df_details.to_csv(DETAILS_CSV, index=False)
    else:
        log_print(f"Peaks mapped to a TE (via their center): {df_int['peak_id'].nunique()}")
        
        # Merge to recover original start/end of the peaks, dropping the 1bp center coordinates
        df_merged = df_int.merge(all_peaks_df, on='peak_id', how='left')
        
        df_details = df_merged[['peak_id', 'chr_p', 'start_p', 'end_p', 'repName', 
                                'repFamily', 'repClass', 'milliDiv', 'strand_te']]
        df_details.to_csv(DETAILS_CSV, index=False)
        log_print(f"-> Detailed Peak-TE table saved to: {DETAILS_CSV.name}")
    
    log_print("Generating binary matrices (including non-overlapping peaks)...")
    
    if not df_int.empty:
        mat_class = (pd.crosstab(df_int["peak_id"], df_int["repClass"]) > 0).astype(int)
        mat_fam = (pd.crosstab(df_int["peak_id"], df_int["repFamily"]) > 0).astype(int)
    else:
        mat_class = pd.DataFrame(index=all_peaks_series)
        mat_fam = pd.DataFrame(index=all_peaks_series)
        
    # --- Class Matrix ---
    mat_class = mat_class.reindex(all_peaks_series, fill_value=0)
    mat_class["any_TE"] = mat_class.max(axis=1) if not df_int.empty else 0
    mat_class.to_csv(CLASS_MATRIX_CSV)
    
    # --- Family Matrix ---
    mat_fam = mat_fam.reindex(all_peaks_series, fill_value=0)
    mat_fam["any_TE"] = mat_fam.max(axis=1) if not df_int.empty else 0
    mat_fam.to_csv(FAMILY_MATRIX_CSV)
    
    log_print(f"-> Binary Class matrix saved to: {CLASS_MATRIX_CSV.name}")
    log_print(f"-> Binary Family matrix saved to: {FAMILY_MATRIX_CSV.name}")

#==================================================
# Main
#==================================================
def main():
    # Safe cleanup of specific script outputs to prevent leftover files
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for p in OUTPUTS:
        p.unlink(missing_ok=True)
    
    # Initialize fresh log file
    with open(LOG_PATH, 'w') as f:
        f.write("=== TE ANNOTATION PIPELINE START ===\n")
        
    try:
        check_dependencies()
        format_rmsk()
        extract_peak_centers()
        run_intersect()
        process_intersection()
    except Exception as e:
        log_print(f"\nTE ANNOTATION PIPELINE FAILED WITH ERROR: {str(e)}")
        raise
    
    log_print("\n=== TE ANNOTATION PIPELINE FINISHED ===")

if __name__ == "__main__":
    main()