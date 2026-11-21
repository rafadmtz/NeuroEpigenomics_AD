"""
04_pca_umap_atac_by_subclass.py
==========================================================================

Calculate UMAP by subclass and for all subclasses combined:
  1. Reads the AnnData files for each condition and groups them by subclass.
  2. Calculates raw-count QC metrics and removes all-zero nuclei and peaks.
  3. Selects the most accessible peaks and calculates a spectral embedding
     with SnapATAC2.
  4. Calculates neighbors and UMAP on the spectral embedding.
  5. Saves the UMAP of each subclass colored by total counts and supertype,
     and a separate figure with one panel per condition.
  6. Concatenates all subclasses, repeats steps 3 and 4, and saves the
     UMAP colored by subclass and supertype, and a separate figure with
     one panel per condition.

Use:
  Run the script from the project root:

      python -u src/04_umap_atac_by_subclass.py

  Run the script in the background:

      nohup python -u src/04_umap_atac_by_subclass.py > logs/04_umap_atac_by_subclass.py.log 2>&1 &

  The script can also be imported as a module to call its functions.
"""

from pathlib import Path

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import polars as pl
import snapatac2 as snap
from scipy import sparse
import matplotlib.pyplot as plt
from importlib import import_module

module = import_module("02_group_atac_by_condition_and_subclass")
load_scanpy = module.load_scanpy

PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = PROJECT_ROOT / "data" / "processed"
OUTPUT_DIR = PROJECT_ROOT / "results" / "figures"

GROUPS= ["Not AD", "Low", "Intermediate", "High"]
CONDITION_KEY = "Overall AD neuropathological Change"
WSPACE = 0.6

if not OUTPUT_DIR.exists():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_subclass_data(
    subclass: str
    ) -> ad.AnnData:

    """
    Read and load AnnData for each condition and group by subclass
    in a single AnnData.

    Args:
        subclass: The subclass name used in the filenames, such as
            "Vip", "Sst" or "L2-3-IT".

    Returns:
        ad.AnnData: An AnnData containing the nuclei from all conditions
            and the peaks shared by the input files.
    """

    grouped_anndata=[]

    for group in GROUPS:
        group_name = group.replace(" ", "_")

        path = INPUT_DIR / f"{group_name}_{subclass}.h5ad"

        adata = load_scanpy(path)

        try:
            grouped_anndata.append(adata.to_memory())
        finally:
            adata.file.close()

    return ad.concat(grouped_anndata, axis=0, join="inner")


def prepare_data(
    adata: ad.AnnData
    ) -> ad.AnnData:

    """
    Calculate raw-count QC metrics and remove all-zero nuclei and peaks.

    Args:
        adata: An AnnData containing the original counts in adata.X,
            with nuclei as rows and peaks as columns.

    Returns:
        ad.AnnData: A filtered copy containing nuclei and peaks with
            positive total counts. Per-nucleus metrics are stored in
            obs["total_counts"] and obs["n_peaks_detected"].
    """

    adata.obs["total_counts"] = np.asarray(
        adata.X.sum(axis=1)
    ).ravel()

    adata.obs["n_peaks_detected"] = np.asarray(
        (adata.X > 0).sum(axis=1)
    ).ravel()

    keep_nuclei = adata.obs["total_counts"].to_numpy() > 0

    keep_peaks = np.asarray(
        adata.X.sum(axis=0)
    ).ravel() > 0

    return adata[keep_nuclei, keep_peaks].copy()

def run_umap(
    adata: ad.AnnData,
    rep= "X_spectral",
    n_neighbors=30,
    random_state=2
    ):

    """
    Select features, calculate the spectral embedding, neighbors and UMAP,
    and save the results in the AnnData.

    Args:
        adata: An AnnData containing the original counts in adata.X,
            after removing all-zero nuclei and peaks.
        rep: The key in adata.obsm used to find neighbors.
        n_neighbors: The neighborhood size used to build the graph.
        random_state: The seed used for the random calculations.

    Returns:
        None. Results are stored directly in adata:
            var["selected"]: Peaks used for the spectral embedding.
            obsm["X_spectral"]: Spectral embedding coordinates.
            uns["neighbors"]: Neighbor settings.
            obsp["distances"]: Neighbor distances.
            obsp["connectivities"]: Graph connections.
            obsm["X_umap"]: UMAP coordinates.
    """

    snap.pp.select_features(adata, n_features=50_000)

    snap.tl.spectral(adata)

    sc.pp.neighbors(
        adata,
        use_rep=rep, 
        n_neighbors=n_neighbors, 
        random_state=random_state
        )

    sc.tl.umap(
        adata, 
        random_state=random_state
        )
    
def plot_umap(
    adata: ad.AnnData,
    name: str,
    color: list[str]
    ):

    """
    Plot the calculated UMAP and save the figure as a PNG.

    Args:
        adata: The AnnData containing the UMAP coordinates
            and the nucleus metadata.
        name: The name used in the output filename.
        color: The obs columns used to color each panel.

    Returns:
        None. The figure is saved in OUTPUT_DIR.
    """

    fig = sc.pl.umap(
        adata,
        color=color,
        ncols=2,
        wspace=WSPACE,
        show=False,
        return_fig=True,
    )

    fig.savefig(OUTPUT_DIR / f"{name}_umap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

def plot_umap_by_condition(adata: ad.AnnData, name: str):

    """
    Plot one UMAP panel per condition. The nuclei of the condition
    are colored and the rest are shown in gray in the background.

    Args:
        adata: The AnnData containing the UMAP coordinates
            and the nucleus metadata.
        name: The name used in the output filename.

    Returns:
        None. The figure is saved in OUTPUT_DIR.
    """

    fig, axes = plt.subplots(1, len(GROUPS), figsize=(5 * len(GROUPS), 5))

    for ax, group in zip(axes, GROUPS):
        sc.pl.umap(
            adata,
            color=CONDITION_KEY,
            groups=[group],
            na_color="lightgray",
            title=group,
            legend_loc=None,
            ax=ax,
            show=False,
        )

    fig.savefig(
        OUTPUT_DIR / f"{name}_umap_by_condition.png",
        dpi=300,
        bbox_inches="tight"
        )
    plt.close(fig)

def main():

    """
    Run the UMAP workflow for the selected subclasses.

    For each subclass, load and prepare the data, calculate the UMAP
    and save the figures. Then concatenate all subclasses, calculate
    the UMAP of the combined data and save the figures.

    Returns:
        None.
    """

    subclasses=[
        "Vip",
        "Sst",
        "L2/3 IT"
    ]

    all_subclasses=[]

    for subclass in subclasses:

        subclass = subclass.replace(" ", "-").replace("/", "-")

        adata = prepare_data(load_subclass_data(subclass))

        all_subclasses.append(adata.copy())

        run_umap(
            adata=adata
            )

        plot_umap(
            adata=adata,
            name=subclass,
            color=["total_counts", "Supertype"]
            )

        plot_umap_by_condition(
            adata=adata,
            name=subclass
            )

    combined = ad.concat(
        all_subclasses, 
        axis=0, 
        join="inner"
        )

    run_umap(adata=combined)

    plot_umap(
        adata=combined,
        name="all_subclasses",
        color=["Subclass", "Supertype"]
        )

    plot_umap_by_condition(
        adata=combined,
        name="all_subclasses"
        )


if __name__ == "__main__":
    main()