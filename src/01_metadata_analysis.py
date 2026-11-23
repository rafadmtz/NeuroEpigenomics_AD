#!/usr/bin/env python3
"""
01_metadata_analysis.py — Metadata analysis for the NeuroEpigenomics_AD project
==========================================================================

Analyzes snATAC-seq metadata from the middle temporal gyrus (MTG) of SEA-AD,
from Alzheimer's (AD) donors and controls.

SEA-AD Source: AWS S3 public bucket "sea-ad-single-cell-profiling"
(Gabitto et al. 2024, Nat Neurosci, "Data availability" section).

Requires having previously run src/00_download_data.py

Inputs:
  1.- data/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad
      Filtered data (final nuclei only), opened in read-only mode (backed='r').
  2.- data/metadata/SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv
      Metadata for all nuclei, including those removed by filtering.

Steps:
  1.- Load the h5ad and csv and verify that the nuclei IDs match.
  2.- Merge both metadata sets. If a column is in both, the final h5ad annotation is kept.
  3.- Define the severity group (Severity_Group) from 'Overall AD neuropathological Change' and 'Neurotypical reference'.
  4.- Check if ATAC quality metrics are calculated per nucleus or per library.
  5.- Summarize clinical information by donor and verify that each donor has a single clinical profile.
  6.- Extract ATAC peak coordinates in BED format (chromosomes chr1-22, X, Y and M) to cross-reference with the rmsk table.

Outputs (in results/01_metadata_output/):
  1.- complete_nuclei_metadata.csv        Merged metadata, one row per nucleus.
  2.- extended_donor_metadata.csv         Clinical variables per donor and number of nuclei (n_nuclei), sorted by Continuous Pseudo-progression Score.
  3.- atac_peaks.bed                      Peak coordinates (chr, start, end, peak_id), sorted. Used in steps 05, 07 and 08.
  4.- 00_qc_and_merge_log.txt             Log of all the above.

Usage:
  The file can be executed directly from the command line, or
  imported as a module to call the main() function.

"""

#==================================================
# Libraries used
#==================================================

import scanpy as sc
import pandas as pd
from pathlib import Path
 
#==================================================
# Variable definitions
#==================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
OUT_DIR = PROJECT_ROOT / "results" / "01_metadata_output"

h5ad_path = DATA_DIR / "SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad"
csv_path = DATA_DIR / "metadata" / "SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv"

nuclei_out = OUT_DIR / "complete_nuclei_metadata.csv"
qc_log_path = OUT_DIR / "00_qc_and_merge_log.txt"

log_lines = []
 
#==================================================
# Functions
#==================================================

def log_print(msg):
    """
    Print a message and keep it in log_lines for the run log.

    Parameters
    ----------
    msg : str
        Message to print and log.
    """
    print(msg)
    log_lines.append(msg)
 
def load_data(h5ad_path, csv_path):
    """
    Load the final-nuclei h5ad (read-only) and the all-nuclei metadata CSV.

    Parameters
    ----------
    h5ad_path : Path
        SEAAD_MTG_ATACseq_final-nuclei h5ad file.
    csv_path : Path
        SEAAD_MTG_ATACseq_all-nuclei_metadata CSV file.

    Returns
    -------
    adata : anndata.AnnData
        h5ad object opened in backed mode (backed='r').
    obs_h5ad : pd.DataFrame
        Copy of adata.obs, one row per final nucleus.
    df_meta : pd.DataFrame
        All-nuclei metadata, indexed by nucleus ID.
    """
    log_print(f"Loading object (backed='r'): {h5ad_path.name}")
    adata = sc.read_h5ad(h5ad_path, backed="r")
    obs_h5ad = adata.obs.copy()

    log_print(f"Loading metadata CSV: {csv_path.name}")
    df_meta = pd.read_csv(csv_path, index_col=0, low_memory=False)

    return adata, obs_h5ad, df_meta
 
def compare_indices(obs_h5ad, df_meta):
    """
    Check that the nuclei IDs of the h5ad and the CSV match.

    Parameters
    ----------
    obs_h5ad : pd.DataFrame
        Nuclei metadata of the h5ad (output of load_data).
    df_meta : pd.DataFrame
        All-nuclei metadata of the CSV (output of load_data).

    Raises
    ------
    AssertionError
        If no nucleus ID is shared by both tables.
    """
    n_h5ad = len(obs_h5ad)
    n_csv = len(df_meta)
    matches = obs_h5ad.index.intersection(df_meta.index)

    log_print(f"Nuclei in h5ad: {n_h5ad}")
    log_print(f"Nuclei in CSV: {n_csv}")
    log_print(f"Exact matches: {len(matches)}")

    assert len(matches) > 0, "Error: The CSV index does not match the h5ad index."
    if len(matches) < n_h5ad:
        log_print(f"Warning: Missing {n_h5ad - len(matches)} h5ad nuclei in the CSV.")
 
def merge_metadata(obs_h5ad, df_meta):
    """
    Add the CSV columns that are missing in the h5ad metadata.

    When a column is in both tables, the final h5ad annotation is kept.

    Parameters
    ----------
    obs_h5ad : pd.DataFrame
        Nuclei metadata of the h5ad.
    df_meta : pd.DataFrame
        All-nuclei metadata of the CSV.

    Returns
    -------
    pd.DataFrame
        Merged metadata, one row per nucleus of the h5ad.
    """
    new_cols = [c for c in df_meta.columns if c not in obs_h5ad.columns]
    obs = obs_h5ad.join(df_meta[new_cols], how='left')
    return obs
 
def clean_and_define_severity(obs):
    """
    Define the severity group of every nucleus.

    Severity_Group is taken from 'Overall AD neuropathological Change';
    donors flagged as 'Neurotypical reference' become 'Not AD (Neurotypical)',
    and missing values become 'Unknown'.

    Parameters
    ----------
    obs : pd.DataFrame
        Merged metadata (output of merge_metadata).

    Returns
    -------
    pd.DataFrame
        The same table with the Severity_Group column.
    """
    log_print("\nBuilding severity groups...")
    obs['Severity_Group'] = obs['Overall AD neuropathological Change'].astype(str)
    mask_neurotypical = obs['Neurotypical reference'].isin(['Yes', 'True', True, 1, '1'])
    obs.loc[mask_neurotypical, 'Severity_Group'] = 'Not AD (Neurotypical)'
    obs['Severity_Group'] = obs['Severity_Group'].replace({'nan': 'Unknown', 'None': 'Unknown'})

    log_print("Severity Counts:")
    log_print(obs['Severity_Group'].value_counts(dropna=False).to_string())
    return obs
 
def save_nuclei_level_table(obs, output_path):
    """
    Save the nuclei-level metadata table.

    Parameters
    ----------
    obs : pd.DataFrame
        Metadata with Severity_Group, one row per nucleus.
    output_path : Path
        Output CSV (complete_nuclei_metadata.csv).
    """
    obs.to_csv(output_path)
    log_print(f"-> Nuclei-level metadata exported to: {output_path.name}")
 
def atac_seq_quality(obs):
    """
    Report whether each ATAC quality metric is per nucleus or per library.

    A metric with a single value per library ('library_prep', or 'sample_id'
    if missing) is computed at the library level; otherwise, per nucleus.

    Parameters
    ----------
    obs : pd.DataFrame
        Nuclei-level metadata.
    """
    log_print("\nAnalyzing ATAC-seq quality metrics...")
    atac_qc_cols = [
        'ATAC_Median_high_quality_fragments_per_cell', 
        'ATAC_TSS_enrichment_score', 
        'ATAC_Fraction_of_high_quality_fragments_in_cells',
        'Doublet score'
    ]
    valid_qc_cols = [c for c in atac_qc_cols if c in obs.columns]

    lib_col = 'library_prep' if 'library_prep' in obs.columns else 'sample_id'
    for col in valid_qc_cols:
        max_unique_per_lib = obs.groupby(lib_col, observed=True)[col].nunique().max()
        level_type = "LIBRARY" if max_unique_per_lib == 1 else "NUCLEUS"
        log_print(f" - {col}: Metric calculated at the {level_type} level")
 
def export_donor_metadata(obs, out_dir):
    """
    Build and save the donor-level clinical table.

    Clinical columns that are not found are skipped with a warning. Each
    donor must have a single clinical profile. The number of nuclei per donor
    (n_nuclei) is added and donors are sorted by Continuous Pseudo-progression
    Score.

    Parameters
    ----------
    obs : pd.DataFrame
        Nuclei-level metadata.
    out_dir : Path
        Output folder of extended_donor_metadata.csv.

    Raises
    ------
    ValueError
        If a donor has more than one clinical profile.
    """
    log_print("\nGenerating donor-level table...")
    desired_donor_cols = [
        'Donor ID', 'Continuous Pseudo-progression Score', 'Overall AD neuropathological Change', 
        'Braak', 'Thal', 'CERAD score', 'Sex', 'Age at Death', 'APOE Genotype',
        'Severity_Group', 'Neurotypical reference'
    ]
    
    missing = [c for c in desired_donor_cols if c not in obs.columns]
    if missing:
        log_print(f"Warning: The following clinical columns were not found and will be skipped: {missing}")

    donor_cols = [c for c in desired_donor_cols if c in obs.columns]
    donor_metadata = obs[donor_cols].drop_duplicates()

    cps_counts = donor_metadata.groupby('Donor ID').size()
    if (cps_counts > 1).any():
        raise ValueError(f"Error: Donors with multiple clinical profiles detected:\n{cps_counts[cps_counts > 1]}")

    nuc_counts = obs['Donor ID'].value_counts().rename("n_nuclei")
    donor_metadata = donor_metadata.merge(nuc_counts, left_on='Donor ID', right_index=True)
    donor_metadata = donor_metadata.sort_values('Continuous Pseudo-progression Score').reset_index(drop=True)

    donor_out = out_dir / "extended_donor_metadata.csv"
    donor_metadata.to_csv(donor_out, index=False)
    log_print(f"-> Donor metadata exported to: {donor_out.name}")
 
def extract_BED(adata, out_dir):
    """
    Extract the peak coordinates of the h5ad as a sorted BED file.

    Peak IDs are parsed as chr:start-end; peaks that cannot be parsed or that
    are not on chr1-22, X, Y or M are skipped.

    Parameters
    ----------
    adata : anndata.AnnData
        h5ad object (peaks in adata.var).
    out_dir : Path
        Output folder of atac_peaks.bed (chr, start, end, peak_id).
    """
    log_print("\nProcessing BED peaks...")
    var_df = adata.var.copy().reset_index().rename(columns={'index': 'peak_id'})

    coord_split = var_df['peak_id'].str.extract(r'(chr[\w\d]+)[:\-](\d+)[:\-](\d+)')
    var_df['chr'] = coord_split[0]
    var_df['start'] = coord_split[1]
    var_df['end'] = coord_split[2]

    failures = var_df['chr'].isna().sum()
    log_print(f"Total variables (peaks): {len(var_df)}")
    if failures > 0:
        log_print(f"Warning: {failures} peaks could not be parsed and will be skipped.")

    bed_df = var_df.dropna(subset=['chr', 'start', 'end']).copy()
    bed_df = bed_df[bed_df['chr'].str.match(r'^chr([\d]+|[XYM])$')]
    bed_df = bed_df[['chr', 'start', 'end', 'peak_id']]

    bed_df['start'] = pd.to_numeric(bed_df['start'], errors='coerce')
    bed_df['end'] = pd.to_numeric(bed_df['end'], errors='coerce')
    bed_df = bed_df.dropna(subset=['start', 'end'])
    bed_df['start'] = bed_df['start'].astype(int)
    bed_df['end'] = bed_df['end'].astype(int)

    bed_df = bed_df.sort_values(['chr', 'start', 'end'])

    var_out = out_dir / "atac_peaks.bed"
    bed_df.to_csv(var_out, sep='\t', header=False, index=False)
    log_print(f"Validated BED file ({len(bed_df)} clean peaks) saved to: {var_out.name}")
 
#==================================================
# Main
#==================================================
def main():

    """
    Run the metadata analysis: load and merge the metadata, define the
    severity groups, export the nuclei and donor tables and the peak BED,
    and save the log.
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load and safely merge data (h5ad vs CSV)
    adata, obs_h5ad, df_meta = load_data(h5ad_path, csv_path)

    # Verify integrity in ID matching
    compare_indices(obs_h5ad, df_meta)

    # Integrate data prioritizing the final h5ad annotation
    obs = merge_metadata(obs_h5ad, df_meta)

    # Clean and define severity
    obs = clean_and_define_severity(obs)

    # Save the nuclei-level table
    save_nuclei_level_table(obs, nuclei_out)

    # ATAC quality exploration 
    atac_seq_quality(obs)

    # Export donor metadata
    export_donor_metadata(obs, OUT_DIR)
    
    # Extract BED for intersection with rmsk
    extract_BED(adata, OUT_DIR)

    # Save logs
    with open(qc_log_path, 'w') as f:
        f.write("\n".join(log_lines))


if __name__ == "__main__":
    main()