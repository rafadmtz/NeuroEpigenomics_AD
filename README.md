# Epigenomic dynamics during Alzheimer's disease progression

Exploratory project that analyzes single-nucleus ATAC-seq (snATAC-seq) data to
track chromatin accessibility changes associated with the reactivation of
transposable elements (TEs) in Alzheimer's disease (AD).

For the project rationale, research questions, methodology and results, see the
[project report](docs/project_report.md).

## Features

- Download of the SEA-AD snATAC-seq dataset (middle temporal gyrus, MTG), its
  nucleus metadata, the RepeatMasker hg38 annotation, the hg38 genome, the
  GENCODE gene annotation and MSigDB gene sets.
- Merging and validation of nucleus metadata, definition of AD
  neuropathological change (ADNC) groups and a per-donor clinical table.
- Export of ATAC peak coordinates in BED format for intersection with
  RepeatMasker.
- Splitting of the count matrix into one file per ADNC category and cell
  subclass.
- Per-donor pseudobulk aggregation by subclass, with quality-control plots.
- De novo UMAP embeddings of the L2/3 IT, Sst and Vip subclasses with SnapATAC2.
- Annotation of each peak with the transposable elements (class and family)
  that overlap its center.
- Differential accessibility per subclass with PyDESeq2: trend across ADNC
  stages and High vs Not AD contrast.
- Per-peak table joining differential accessibility, TE annotation and GC
  content.
- Functional (Hallmark, GO Biological Process) and TE enrichment of
  differentially accessible peaks, adjusted for accessibility, GC content,
  peak width and distance to the nearest TSS.
- Over- and under-representation of TE families (LINE, LTR, SINE) and LTR/ERV1
  subfamilies among differentially accessible peaks.

## Repository structure

```text
NeuroEpigenomics_AD/
├── data/            # Local data (large files are not tracked by Git)
├── docs/            # Report and citation metadata (CITATION.cff)
├── notebooks/       # Exploratory analyses
├── results/         # Generated tables and figures
├── src/             # Source code (pipeline scripts)
├── tests/           # Tests
├── environment.yml  # Conda environment definition
├── LICENSE          # License
└── README.md        # Introduction and quick-start guide
```

## Background and motivation

Recent studies show that, in tauopathy models (PS19 mice) and in brains of
patients with AD, tau aggregation destabilizes H3K9me3-marked constitutive
heterochromatin, displacing HP1α from nucleosomes and eroding the nuclear
lamina. This loss of condensation reactivates normally silenced transposable
elements (mainly LINE and LTR), producing double-stranded transcripts that adopt
a Z conformation (Z-RNA) and activate ZBP1, triggering neuronal necroptosis
(Liu et al., *Nature Neuroscience*, 2026). In parallel, cell type-specific
atlases of AD brain tissue make it possible to place these phenomena within the
cellular diversity of the brain.

Building on this, the project uses snATAC-seq data to explore whether chromatin
accessibility changes at single-nucleus resolution overlap the genomic regions
where these transposable elements are located, as an indirect way to track their
possible reactivation during disease progression.

## Requirements

All dependencies are listed in [`environment.yml`](environment.yml) and are
installed when the Conda environment is created (see [Usage](#usage)). The main
ones are:

- Python 3.11, NumPy, SciPy, pandas, polars, h5py and joblib.
- scanpy and anndata, to read and process the `.h5ad` files.
- SnapATAC2, for the UMAP embeddings (step 04).
- bedtools, for the peak–TE intersection and the GC content of the peaks
  (steps 05, 07 and 08).
- PyDESeq2, for the differential accessibility analysis (step 06).
- matplotlib and seaborn, for the plots.

The project uses only this Conda environment; no other virtual environment
(e.g. a local `.venv`) is needed.

The input `.h5ad` file is about 18 GB, so the pipeline is designed to run on a
server; the pseudobulk step is tuned for roughly 64 GB of RAM.

## Data

The snATAC-seq data come from the **SEA-AD** (Seattle Alzheimer's Disease Brain
Cell Atlas) initiative, described in Gabitto *et al.* (2024). The download script
retrieves the following files from the public AWS S3 bucket
`sea-ad-single-cell-profiling` and from UCSC:

| File | Source | Destination |
|---|---|---|
| `SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad` | SEA-AD (MTG, final nuclei) | `data/` |
| `SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv` | SEA-AD (all-nuclei metadata) | `data/metadata/` |
| `rmsk.txt.gz` | UCSC Genome Browser, RepeatMasker hg38 | `data/` |

No primary SEA-AD data are redistributed through this repository. Users are
responsible for complying with the access conditions, licenses and citation
requirements of the original datasets.

## Usage

### Environment setup

From the root of the `NeuroEpigenomics_AD/` directory, create the environment:

```bash
conda env create -f environment.yml
```

Activate it before running any script in `src/`:

```bash
conda activate neuroAD
```

If the environment already exists and `environment.yml` has changed, update it:

```bash
conda env update -n neuroAD -f environment.yml
```

Deactivate it when you are done:

```bash
conda deactivate
```

### Running the pipeline

Run the scripts from the repository root, in this order:

```bash
python src/00_download_data.py
python src/01_metadata_analysis.py
python src/01b_filtered_data_graphs.py
python src/02_group_atac_by_condition_and_subclass.py
python src/03_pseudobulk_by_donor.py
python src/04_umap_atac_by_subclass.py
python src/05_annotate_peaks_with_tes.py
python src/06_differential_accessibility.py
python src/07_build_peak_table.py
python src/08_functional_and_te_enrichment.py
python src/09_te_analysis.py
```

The longer steps can be run in the background while saving a log, for example:

```bash
nohup python -u src/03_pseudobulk_by_donor.py > logs/03_pseudobulk_by_donor.log 2>&1 &
```

Each script can also be imported as a module and run through its `main()`
function.

## Outputs

| Step | Output | Location |
|---|---|---|
| `00_download_data.py` | SEA-AD `.h5ad` and metadata, RepeatMasker, hg38 genome, GENCODE annotation, MSigDB gene sets | `data/` and its subfolders `metadata/`, `genome/`, `annotation/`, `genesets/` |
| `01_metadata_analysis.py` | Per-nucleus metadata, per-donor clinical table, ATAC peak BED file and log | `results/01_metadata_output/` |
| `01b_filtered_data_graphs.py` | Subclass and supertype composition plots | `results/01b_graphs/` |
| `02_group_atac_by_condition_and_subclass.py` | One `.h5ad` per ADNC category × subclass | `data/processed/condition_subclass_grouped/` |
| `03_pseudobulk_by_donor.py` | Per-subclass pseudobulk `.h5ad`, sample table (colData) and summary | `data/processed/pseudobulk/` |
| `03_pseudobulk_by_donor.py` | Pseudobulk QC plots | `results/figures/pseudobulk/` |
| `04_umap_atac_by_subclass.py` | UMAPs by subclass, supertype, total counts and ADNC category | `results/figures/` |
| `05_annotate_peaks_with_tes.py` | Peak–TE overlap table and binary peak × TE class / family matrices | `results/05_te_annotation/` |
| `06_differential_accessibility.py` | Per-subclass DA results (trend and High vs Not AD), VST matrices and summary | `results/06_differential_accessibility/` |
| `06_differential_accessibility.py` | PCA, volcano, MA, dispersion and per-subclass significance plots | `results/figures/06_differential_accessibility/` |
| `07_build_peak_table.py` | Per-peak tables (DA + TE + GC) for each subclass and contrast | `results/07_build_peak_table/` |
| `08_functional_and_te_enrichment.py` | Peak annotation, functional and TE enrichment tables | `results/08_enrichment/` |
| `08_functional_and_te_enrichment.py` | Functional enrichment dotplots | `results/08_graphs/` |
| `09_te_analysis.py` | TE family and ERV1 subfamily enrichment tables | `results/09_te_analysis/` |
| `09_te_analysis.py` | Family bubble plots and ERV1 subfamily plot | `results/figures/09_te_analysis/` |

## Documentation

- [Project report](docs/project_report.md)
- [Citation metadata](docs/CITATION.cff)

## Team

- Rafael Díaz Martínez
- Santiago Martínez Enciso
- Sixto Netzahualcóyotl González Morales

## Citation

If you use this software, please cite the repository (release v1.0.0):

> Díaz Martínez, R., Martínez Enciso, S. and González Morales, S. N. (2026).
> *NeuroEpigenomics_AD: Epigenomic dynamics during Alzheimer's disease
> progression*. GitHub repository.
> [https://github.com/rafadmtz/NeuroEpigenomics_AD](https://github.com/rafadmtz/NeuroEpigenomics_AD)

Citation metadata are provided in [CITATION.cff](docs/CITATION.cff); GitHub uses
this file to show the **Cite this repository** option on the repository page.

If you use the SEA-AD data, please also cite the original atlas:

> Gabitto, M. I. *et al.* Integrated multimodal cell atlas of Alzheimer's
> disease. *Nature Neuroscience* **27**, 2366–2383 (2024).
> [https://doi.org/10.1038/s41593-024-01774-5](https://doi.org/10.1038/s41593-024-01774-5)

This section will be updated with a versioned DOI and a publication citation
when they become available.

## License

The source code in this repository is distributed under the
[MIT License](LICENSE). External tools and datasets (SEA-AD, RepeatMasker/UCSC)
remain subject to their respective licenses and terms of use.
