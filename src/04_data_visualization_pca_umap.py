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
    
def load_subclass_data (
    subclass: str
    )  -> ad.AnnData:
    
    """
    Read and load AnnData for each condition and subclass and group by subclass in a single AnnData.
    Args:
        subclass: The name of the subclass to group by.

    Returns:
        ad.AnnData
    """
    
    for group in GROUPS:
        group_name = group.replace(" ", "_")

        path = INPUT_DIR / f"{group_name}_{subclass}.h5ad"
        
        grouped_anndata=[]
        
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
    layer = None,
    random_state: int = 2
    ):
    
    """
    Calculate pca and umap and save it in the anndata for each tag assigned.
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
        n_cols=2,
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
                subclass= subclass
            )
        
            plot_embeddings(
                adata=anndata,
                tag="tfidf",
                typeplot=plot,
                subclass= subclass
            )
            
if __name__ == "__main__":
    main()