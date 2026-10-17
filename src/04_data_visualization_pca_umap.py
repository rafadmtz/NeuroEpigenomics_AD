"""
04_pca_umap_atac_by_subclass.py
==========================================================================

Calculate PCA and UMAP by subclass:
  1. Reads the AnnData files for each condition and groups them by subclass.
  2. Calculates raw-count QC metrics and removes all-zero nuclei and peaks.
  3. Calculates TF-IDF and saves the transformed matrix in an AnnData layer.
  4. Calculates PCA, neighbors and UMAP for raw counts and TF-IDF values.
  5. Saves the plots colored by condition, donor, total counts and supertype.

Use:
  Run the script from the project root:

      python -u src/04_pca_umap_atac_by_subclass.py

  Run the script in the background:

      nohup python -u src/04_pca_umap_atac_by_subclass.py > logs/04_pca_umap_atac_by_subclass.log 2>&1 &

  The script can also be imported as a module to call its functions.
"""

from pathlib import Path

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import polars as pl
from scipy import sparse
import matplotlib.pyplot as plt
from importlib import import_module

module = import_module("03_group_atac_by_condition_and_subclass")
load_scanpy = module.load_scanpy

PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = PROJECT_ROOT / "data" / "processed"
OUTPUT_DIR = PROJECT_ROOT / "results" / "figures"

GROUPS= ["Not AD", "Low", "Intermediate", "High"]

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


def run_pca_umap(
    tag: str,
    adata: ad.AnnData,
    n_components: int = 50,
    n_neighbors: int = 30,
    layer=None,
    random_state: int = 2
    ):

    """
    Calculate PCA, neighbors and UMAP and save the results in the
    AnnData for each tag assigned.

    Args:
        tag: The label used to identify the results, such as "raw"
            or "tfidf".
        adata: The AnnData containing the data to analyze.
        n_components: The number of PCA components to calculate and
            use to find neighbors. Must be smaller than both the
            number of nuclei and the number of peaks.
        n_neighbors: The neighborhood size used to build the graph.
        layer: The name of the layer to use for PCA. If None,
            use adata.X.
        random_state: The seed used for the random calculations.

    Returns:
        None. Results are stored directly in adata:
            obsm[f"X_pca_{tag}"]: PCA coordinates.
            uns[f"neighbors_{tag}"]: Neighbor settings.
            obsp[f"neighbors_{tag}_distances"]: Neighbor distances.
            obsp[f"neighbors_{tag}_connectivities"]: Graph connections.
            obsm[f"X_umap_{tag}"]: UMAP coordinates.
    """

    sc.pp.pca(
        data=adata,
        n_comps=n_components,
        layer=layer,
        random_state=random_state,
        key_added=f"X_pca_{tag}",
        svd_solver="arpack",
        zero_center=True
    )

    sc.pp.neighbors(
        adata=adata,
        n_neighbors=n_neighbors,
        use_rep=f"X_pca_{tag}",
        n_pcs=n_components,
        random_state=random_state,
        key_added=f"neighbors_{tag}"
    )

    sc.tl.umap(
        adata=adata,
        neighbors_key=f"neighbors_{tag}",
        n_components=2,
        random_state=random_state,
        key_added=f"X_umap_{tag}"
    )


def compute_tfidf(
    adata: ad.AnnData,
    scale_factor: float = 10000
    ) -> ad.AnnData:

    """
    Calculate TF-IDF from the original counts and save it in an
    AnnData layer.

    TF divides each count by the total counts of its nucleus.
    IDF divides the number of nuclei by the total counts of each peak.
    The final values are calculated as log1p(scale_factor * TF * IDF).

    Args:
        adata: An AnnData containing a sparse matrix of original counts
            in adata.X, after removing all-zero nuclei and peaks.
        scale_factor: The factor applied before the log1p transformation.

    Returns:
        ad.AnnData: The same AnnData with the transformed matrix in
            layers["tfidf"] and the original counts in adata.X.
    """

    n_nuclei=adata.n_obs

    counts_per_nucleus= np.asarray(
        (adata.X).sum(axis=1)
    ).ravel()

    counts_per_peak= np.asarray(
        (adata.X).sum(axis=0)
    ).ravel()

    tf = adata.X.multiply(
        (1 / counts_per_nucleus)[:, None]
    )

    idf= n_nuclei / counts_per_peak

    tfidf = tf.multiply(idf)

    tfidf.data = np.log1p(tfidf.data * scale_factor)

    adata.layers["tfidf"] = tfidf.tocsr()

    return adata


def plot_embeddings(
    adata: ad.AnnData,
    tag: str,
    typeplot: str,
    subclass: str
    ):

    """
    Plot the calculated PCA or UMAP and save the figure as a PNG.

    Each figure contains panels colored by condition, donor,
    total counts and supertype.

    Args:
        adata: The AnnData containing the calculated coordinates
            and the nucleus metadata.
        tag: The label identifying the results, such as "raw"
            or "tfidf".
        typeplot: The representation to plot, either "pca" or "umap".
        subclass: The subclass name used in the output filename.

    Returns:
        None. The figure is saved in OUTPUT_DIR.
    """

    plot= sc.pl.embedding(
        adata=adata,
        basis=f"{typeplot}_{tag}",
        color=[
            "Overall AD neuropathological Change",
            "Donor ID",
            "total_counts",
            "Supertype"
        ],
        dimensions=(0, 1),
        ncols=2,
        show=False,
        return_fig=True
    )

    plot.savefig(
        OUTPUT_DIR / f"{subclass}_{typeplot}_{tag}.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close(plot)


def main():

    """
    Run the PCA and UMAP workflow for the selected subclasses.

    For each subclass, load and prepare the data, calculate TF-IDF,
    calculate both representations and save the figures.

    Returns:
        None.
    """

    subclasses=[
        "Vip",
        "Sst",
        "L2/3 IT"
    ]

    plots=[
        "pca",
        "umap"
    ]

    for subclass in subclasses:

        subclass = subclass.replace(" ", "-").replace("/", "-")

        anndata= compute_tfidf(prepare_data(load_subclass_data(subclass)))

        run_pca_umap(
            tag="raw",
            adata=anndata,
            layer=None
        )

        run_pca_umap(
            tag="tfidf",
            adata=anndata,
            layer="tfidf"
        )

        for plot in plots:
            plot_embeddings(
                adata=anndata,
                tag="raw",
                typeplot=plot,
                subclass=subclass
            )

            plot_embeddings(
                adata=anndata,
                tag="tfidf",
                typeplot=plot,
                subclass=subclass
            )

if __name__ == "__main__":
    main()