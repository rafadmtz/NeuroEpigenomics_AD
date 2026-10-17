
"""
03_group_atac_by_condition_and_subclass.py
==========================================================================

Group nuclei by condition and subclass:
  1. Reads nucleus metadata, filters the selected conditions, and saves
     the metadata to a CSV file.
  2. Saves ATAC-seq peak identifiers and groups nuclei by condition
     and subclass, retaining all peaks.
  3. Writes one .h5ad file per condition-subclass combination,
     preserving the original counts and nucleus metadata.
     
Use:
  Run from the directory containing the script:
    python 03_group_atac_by_condition_and_subclass.py

  Run in the background and save output to a log:
    nohup python -u src/03_group_atac_by_condition_and_subclass.py > logs/03_group_atac_by_condition_and_subclass.log 2>&1 &

  Monitor progress:
    tail -f 03_group_atac_by_condition_and_subclass.log

  The script can also be loaded as a module using importlib and executed
  by calling main(). No command-line arguments are required; input and
  output paths are defined in the script.

"""

from pathlib import Path

import scanpy as sc
import anndata as ad
import h5py
import pandas as pd
from anndata._io.specs import read_elem
from joblib import Parallel, delayed

PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_FILE = PROJECT_ROOT / "data" / "SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad"
OUTPUT_DIR = PROJECT_ROOT / "data" / "processed"

GROUPS= ["Not AD", "Low", "Intermediate", "High"]

def load_scanpy(path: str) -> ad.AnnData:
    """
    Load scanpy.
    
    Returns:
        ad.AnnData: An AnnData object containing the ATACseq data.
    """

    adata = sc.read_h5ad(path, backed='r')
    # Return the AnnData object
    return adata


def load_metadata(path: str) -> pd.DataFrame:
    """
    Read and load only metadata  (obs).
    
    Returns:
        pd.DataFrame
    """
    with h5py.File(path) as f:
        obs = read_elem(f["obs"])
    return obs

def load_peaks(path: str) -> pd.DataFrame:
    """
    Load peak metadata from a local .h5ad file.

    Returns:
        pd.DataFrame: Peak metadata, indexed by genomic region.
    """
    with h5py.File(path, "r") as f:
        peaks = read_elem(f["var"])
    
    return peaks


def save_matrix_patient_type(row, path, output_dir) -> None:
    """
    Read and load matrix for each condition and subclass.
    Args:
        row: A row from the DataFrame containing "Overall AD neuropathological Change" and "Subclass".
        adata: The AnnData object containing the dataset.

    Returns:
        pd.DataFrame
    """
    
    adata = load_scanpy(path)
    
    try:
        unique_adata = adata[
            (
                adata.obs["Overall AD neuropathological Change"]
                == row["Overall AD neuropathological Change"]
            )
            & (adata.obs["Subclass"] == row["Subclass"]),
            :
        ].to_memory()
    finally:
        adata.file.close()                       
    
    condition = row["Overall AD neuropathological Change"].replace(" ", "_")
    subclass = row["Subclass"].replace(" ", "-").replace("/", "-")
    
    unique_adata.write_h5ad(
        output_dir / f"{condition}_{subclass}.h5ad",
        compression="gzip"
    )
    
    print(f"Saved matrix for condition {condition} and subclass {subclass}.")
    
def main():
    n_threads = 6
    print("Starting grouping by condition and subclass ...")
    
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print(f"Input file: {INPUT_FILE}")
    print(f"Output directory: {OUTPUT_DIR}")
    
    metadata = load_metadata(INPUT_FILE)

    metadata = metadata.loc[
        metadata["Overall AD neuropathological Change"].isin(GROUPS)
    ]

    metadata.to_csv(
        OUTPUT_DIR / "metadata_test.csv",
        index=True,
        index_label="nucleus_id"
    )

    load_peaks(INPUT_FILE).to_csv(
        OUTPUT_DIR / "peaks.txt",
        index=True,
        header=False
    )

    df_mtg = metadata[
        ["Overall AD neuropathological Change", "Subclass"]
    ].drop_duplicates()
    
    Parallel(n_jobs=n_threads)(
        delayed(save_matrix_patient_type)(row, INPUT_FILE, OUTPUT_DIR)
        for _, row in df_mtg.iterrows()
    )
    
    print("Grouping by condition and subclass completed.")

if __name__ == "__main__":
    main()