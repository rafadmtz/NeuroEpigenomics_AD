# Epigenomic dynamics during Alzheimer's disease progression

Exploratory project that analyzes single-nucleus ATAC-seq (snATAC-seq) data to
track chromatin accessibility changes associated with the reactivation of
transposable elements (TEs) in Alzheimer's disease (AD).

For the project rationale, research questions, methodology and results, see the
[project report](docs/reporte-proyecto.md).

## Features

- Download of the SEA-AD snATAC-seq dataset (middle temporal gyrus, MTG), its
  nucleus metadata and the RepeatMasker hg38 annotation.
- Merging and validation of nucleus metadata, definition of AD
  neuropathological change (ADNC) groups and a per-donor clinical table.
- Export of ATAC peak coordinates in BED format for intersection with
  RepeatMasker.
- Splitting of the count matrix into one file per ADNC category and cell
  subclass.
- Per-donor pseudobulk aggregation by subclass, with quality-control plots.
- De novo UMAP embeddings of the L2/3 IT, Sst and Vip subclasses with SnapATAC2.

## Repository structure

```text
NeuroEpigenomics_AD/
├── data/            # Local data (large files are not tracked by Git)
├── docs/            # Report and documentation
├── notebooks/       # Exploratory analyses
├── results/         # Generated tables and figures
├── src/             # Source code (pipeline scripts)
├── tests/           # Tests
├── environment.yml  # Conda environment definition
├── CITATION.cff     # Citation metadata
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
ones are Python 3.11, scanpy, anndata, SnapATAC2, pandas, NumPy, SciPy, h5py,
joblib and matplotlib.

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

Deactivate it when you are done:

```bash
conda deactivate
```

### Running the pipeline

Run the scripts from the repository root, in this order:

```bash
python src/00_descargar_datos.py
python src/01_analisis_metadatos.py
python src/01b_graficos_datos_filtrados.py
python src/02_group_atac_by_condition_and_subclass.py
python src/03_pseudobulk_by_donor.py
python src/04_umap_atac_by_subclass.py
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
| `01_analisis_metadatos.py` | Per-nucleus metadata, per-donor clinical table, ATAC peak BED file and log | `results/01_metadata_output/` |
| `01b_graficos_datos_filtrados.py` | Subclass and supertype composition plots | `results/01b_graficas/` |
| `02_group_atac_by_condition_and_subclass.py` | One `.h5ad` per ADNC category × subclass | `data/processed/` |
| `03_pseudobulk_by_donor.py` | Per-subclass pseudobulk `.h5ad`, sample table (colData) and summary | `data/processed/pseudobulk/` |
| `03_pseudobulk_by_donor.py` | Pseudobulk QC plots | `results/figures/pseudobulk/` |
| `04_umap_atac_by_subclass.py` | UMAPs by subclass, supertype, total counts and ADNC category | `results/figures/` |

## Documentation

- [Project report](docs/reporte-proyecto.md)
- [Citation metadata](CITATION.cff)

## Team

- Rafael Díaz Martínez
- Santiago Martínez Enciso
- Sixto Nezahualcoyotl Gonzales Morales

## Citation

This software is under active development. If you use it before a formal
release or associated publication is available, please cite the repository:

> Díaz Martínez, R., Martínez Enciso, S. and Gonzales Morales, S. N. (2026).
> *NeuroEpigenomics_AD: Epigenomic dynamics during Alzheimer's disease
> progression*. GitHub repository.
> [https://github.com/rafadmtz/NeuroEpigenomics_AD](https://github.com/rafadmtz/NeuroEpigenomics_AD)

Citation metadata are provided in [CITATION.cff](CITATION.cff); GitHub uses
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