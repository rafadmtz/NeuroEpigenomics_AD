# Epigenomic dynamics of transposable elements during Alzheimer's disease progression: an snATAC-seq analysis of SEA-AD

## General information

| Item | Information |
|:--|:--|
| Member | Rafael Díaz Martínez (rafadiaz@lcg.unam.mx) |
| Member | Santiago Martínez Enciso (santiago.marenc@gmail.com) |
| Member | Sixto Netzahualcóyotl González Morales (sixtogm@lcg.unam.mx) |
| Creation date | 25/08/2026 |
| Last updated | 23/11/2026 |
| Repository | https://github.com/rafadmtz/NeuroEpigenomics_AD |
| Status | Complete: workflow implemented (scripts 00–09) and results obtained |

## Project summary

In tauopathy mouse models and in the brains of Alzheimer's disease (AD) patients, tau aggregation disrupts heterochromatin and reactivates silenced transposable elements (TEs), whose transcripts activate ZBP1 and trigger neuronal necroptosis (Liu et al., 2026). We asked whether this loss of silencing can be detected as a change in chromatin accessibility in human brain tissue. Using snATAC-seq data from the middle temporal gyrus (MTG) of the SEA-AD atlas (Gabitto et al., 2024), we built a reproducible Python workflow that aggregates counts per donor and cell subclass, tests differential accessibility (DA) across AD neuropathological stages with PyDESeq2, and intersects the peaks with RepeatMasker TE annotations. DA was concentrated in excitatory IT neurons (L2/3, L4 and L5 IT), where most changing peaks gained accessibility with disease stage. Among these peaks, LTR/ERV1 elements (especially the MER52 subfamilies) were over-represented, while LINE-1 was under-represented, so the data point to changes in specific ERV-derived regions rather than a general derepression of TEs.

## 1. Context and background

About half of the human genome derives from transposable elements (LINE, SINE and LTR retrotransposons, and DNA transposons). In neurons, most of them are kept silent in constitutive heterochromatin marked by H3K9me3 and bound by HP1α.

Liu et al. (2026) showed in the PS19 tauopathy mouse model, and in AD brains, that tau aggregates sequester H3K9me3-marked chromatin away from HP1α. Heterochromatin decondenses, LINE and LTR elements are reactivated, and their double-stranded transcripts adopt a Z conformation (Z-RNA) that activates ZBP1 and leads to necroptosis.

snATAC-seq measures which genomic regions are open in each nucleus. Since chromatin opening accompanies transcription, increased accessibility over TEs could be an indirect sign of their derepression. The SEA-AD atlas (Gabitto et al., 2024) profiled the MTG of donors across the full AD neuropathological spectrum (ADNC: Not AD, Low, Intermediate, High) with a detailed cell-type annotation, which makes it possible to look for these changes in specific cell populations.

## 2. Problem statement

Evidence that tau reactivates TEs comes mostly from animal models and from transcriptomic or bulk-tissue data. It is not known whether this leaves a footprint on chromatin accessibility in human neurons, in which cell types, or at which disease stage. SEA-AD can address this, but the dataset is large (~18 GB, hundreds of thousands of nuclei), has a multi-donor design that requires donor-level statistics, and must be combined with an external repeat annotation.

## 3. Rationale

Finding TE-associated accessibility changes in specific neuronal subclasses would test whether the mechanism described in mice has a correlate in human tissue and would point to the most vulnerable populations. A version-controlled workflow lets the analysis be repeated with the same data versions and parameters, and reused for other brain regions or SEA-AD modalities.

## 4. General objective

To characterize chromatin accessibility changes in neuronal populations along Alzheimer's disease progression and evaluate their relationship with transposable elements, as an indirect approach to their reactivation.

### Specific objectives

1. Analyze how chromatin accessibility changes across ADNC stages in MTG cell subclasses.
2. Identify which subclasses concentrate these changes and which biological processes are associated with them.
3. Evaluate whether peaks overlapping TEs (with emphasis on LINE-1 and LTR/ERV) gain accessibility as disease progresses.

## 5. Research questions

### Question 1

**Question:** Does chromatin accessibility change along ADNC stages, and in which subclasses?
**Required evidence:** Per-subclass DA results from donor-level pseudobulk counts.

### Question 2

**Question:** Are TE classes or families over-represented among peaks that change accessibility with disease stage?
**Required evidence:** Peak–TE annotation and enrichment of TE families among DA peaks.

## 6. Scope and limitations

### Included

- SEA-AD MTG snATAC-seq data (final-nuclei, version 2024-12-06), the four ADNC stages and the peaks defined by the authors.
- DA in nine subclasses (L2/3, L4, L5 and L6 IT; Vip; Sst; astrocytes, oligodendrocytes and microglia), and TE analyses in the six with enough signal.
- TE annotation from the UCSC RepeatMasker table (hg38), LINE, SINE, LTR, DNA and Retroposon classes.

### Not included

- Other brain regions or modalities (snRNA-seq, multiome), new peak calling, TE expression and experimental validation.

### Known limitations

- Accessibility is an indirect measure; it does not show that TEs are transcribed.
- Young TE copies are nearly identical, so reads on them are often lost at alignment, and the UCSC annotation (RepeatMasker 3.0.1, Repbase 2010) classifies them less accurately than newer versions.
- The Not AD group is small (about 9 donors per subclass), and only the IT subclasses have thousands of DA peaks.
- The family-level enrichment (step 09) does not adjust for peak accessibility or GC content, so it should be taken as exploratory.

## 7. Proposed solution

A workflow of numbered Python scripts in `src/`, run in order from the repository root:

| Step | Script | What it does |
|:--|:--|:--|
| 00 | download_data | Downloads SEA-AD data, RepeatMasker, hg38, GENCODE and MSigDB |
| 01 / 01b | metadata_analysis, filtered_data_graphs | Merges and validates metadata, builds a donor table and cohort plots |
| 02 | group_atac_by_condition_and_subclass | Splits the matrix by ADNC stage and subclass |
| 03 | pseudobulk_by_donor | Sums counts per donor and subclass (≥10 nuclei) |
| 04 | umap_atac_by_subclass | Spectral embedding and UMAP with SnapATAC2 |
| 05 | annotate_peaks_with_tes | Intersects peak centers with TEs (bedtools) |
| 06 | differential_accessibility | PyDESeq2 per subclass |
| 07 | build_peak_table | Joins DA, TE annotation and GC content |
| 08 | functional_and_te_enrichment | Gene-set and covariate-adjusted TE enrichment |
| 09 | te_analysis | TE family and ERV1 subfamily enrichment among DA peaks |

### 7.1 Expected result or product

An MIT-licensed repository with the scripts, the Conda environment, documentation, result tables (pseudobulk, DA, TE annotation and enrichment) and reproducible figures.

## 8. Data

### 8.1 Data sources

| Source | Identifier or version | URL | Access date | License or terms |
|:--|:--|:--|:--|:--|
| SEA-AD (Allen Institute) | SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad | https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad | 15/10/2026 | Allen Institute Terms of Use; cite Gabitto et al. (2024) |
| SEA-AD (Allen Institute) | SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv | https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq/SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv | 15/10/2026 | Allen Institute Terms of Use |
| UCSC Genome Browser | hg38 rmsk table (RepeatMasker 3.0.1, Repbase 2010-03-02) | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz | 15/10/2026 | Free for academic use |
| UCSC Genome Browser | hg38 genome FASTA | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | [dd/mm/2026] | Free for academic use |
| GENCODE | Release 50, basic annotation | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/ | [dd/mm/2026] | Open access |
| MSigDB (Broad Institute) | 2026.2.Hs, Hallmark and GO BP | https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.2.Hs | [dd/mm/2026] | MSigDB license |

### 8.2 Data characteristics

- *Homo sapiens*, post-mortem MTG. Up to 84 donors per subclass after filtering (about 9 Not AD, 12 Low, 21 Intermediate and 42 High).
- snATAC-seq: AnnData object (~18 GB) with a sparse nuclei × peaks count matrix and per-nucleus metadata (donor, cell type, neuropathology, sex, age, PMI). Coordinates in GRCh38.
- RepeatMasker: 5,683,690 repeats; 4,505,027 TE fragments remain after keeping the five TE classes on chr1–22, X and Y.

### 8.3 Data organization

`data/`, `logs/` and `results/` are generated locally and listed in `.gitignore`.

```text
NeuroEpigenomics_AD/
├── data/        # Raw downloads and processed matrices (not tracked)
├── docs/        # Report and CITATION.cff
├── notebooks/   # Exploratory analyses
├── results/     # Tables and figures per step (not tracked)
├── src/         # Scripts 00–09
├── tests/
├── environment.yml
├── LICENSE
└── README.md
```

### 8.4 Data dictionary or format

| File or dataset | Field/column | Type | Description | Values or units |
|:--|:--|:--|:--|:--|
| colData_{subclass}.csv | Overall AD neuropathological Change | categorical | ADNC stage of the donor | Not AD, Low, Intermediate, High |
| colData_{subclass}.csv | Sex, Age at Death, PMI | mixed | Covariates of the DA model | Female/Male; years; hours |
| peak_te_overlaps.csv | repClass, repFamily, repName | string | TE containing the peak center | e.g. LTR / ERV1 / MER52D |
| da_{subclass}_trend.csv.gz | log2FoldChange, padj | numeric | Change per ADNC step and adjusted p-value | log2 units |
| family_enrichment.csv | log2_oe, fdr | numeric | log2(observed/expected) of a TE family among DA peaks and its FDR | log2 units |

## 9. Methodology

### 9.1 Analysis or development stages

1. **Preprocessing (00–02).** Data download with pinned versions; metadata merging and validation (matching identifiers, one clinical profile per donor); split of the matrix by ADNC stage and subclass.
2. **Pseudobulk (03).** Counts are summed per donor and subclass, reading the matrix in chunks of 20,000 nuclei to stay within memory. Donors with fewer than 10 nuclei are removed.
3. **Exploration (04).** Spectral embedding and UMAP of L2/3 IT, Sst and Vip with SnapATAC2 (50,000 features, 30 neighbors).
4. **TE annotation (05).** A peak is assigned to a TE when its center falls inside a RepeatMasker fragment (`bedtools intersect`).
5. **Differential accessibility (06).** PyDESeq2 per subclass, with donors as replicates and the design `~ sex + age + PMI + stage`, where stage is coded 0–3 (trend). A High vs Not AD contrast was also run. Subclasses need ≥3 donors per stage; significance at padj < 0.05.
6. **Enrichment (07–09).** Functional enrichment (Hallmark, GO BP) of genes near DA peaks, adjusted for accessibility, GC, peak width and distance to the TSS. For TE families, the fraction of DA peaks (padj < 0.05, |log2FC| > 0.1) in each family was compared to its fraction among all tested TE peaks with Fisher's exact test and Benjamini–Hochberg correction.

### 9.2 Tools and technologies

| Tool | Version | Purpose |
|:--|:--|:--|
| Python | 3.11 | Workflow language |
| pandas, NumPy, SciPy, polars, statsmodels | environment.yml | Data handling and statistics |
| anndata, scanpy | environment.yml | Reading `.h5ad` files |
| SnapATAC2 | environment.yml | UMAP embeddings |
| PyDESeq2 | environment.yml | Differential accessibility |
| bedtools | environment.yml | Peak–TE intersection and GC content |
| matplotlib, seaborn | environment.yml | Figures |

### 9.3 Validation strategy

Each script checks its inputs and required columns before running. Metadata consistency is verified in step 01. Donor-level pseudobulk avoids treating nuclei as independent replicates, and results were reviewed with quality control plots (library size per stage, PCA by stage and sex, dispersion and MA plots). Multiple testing is corrected with Benjamini–Hochberg. Automated unit tests are not yet implemented.

## 10. Work plan

### 10.1 Division of responsibilities

| Member | Main responsibility | Review responsibility |
|:--|:--|:--|
| [Name] | [Responsibility] | [What or whom they will review] |

### 10.2 Risks and alternatives

| Risk | Probability | Impact | Prevention or alternative |
|:--|:--|:--|:--|
| Dataset too large for memory | Medium | High | Backed reading and chunked aggregation on a server |
| Few donors in some groups | High | Medium | Donor filters and analysis restricted to subclasses with enough replicates |
| A global accessibility shift is mistaken for a TE effect | High | High | Covariate-adjusted tests and effect-size threshold |

## 11. Results

**Cohort and quality control.** L2/3 IT is the most abundant subclass (>140,000 nuclei), followed by L4 and L5 IT. Pseudobulk library sizes were similar across stages. In the UMAPs, nuclei from the four stages overlap, so disease stage is not the main source of variation. PCA showed separation by sex in L2/3 IT, which supports including sex as a covariate.

**Differential accessibility.** With the stage trend, L2/3 IT had 22,354 DA peaks, followed by L4 IT (~20,000) and L5 IT (~12,000). Vip (4,312), astrocytes (~3,400) and Sst (2,053) had fewer, and most other subclasses had almost none. In all subclasses with signal, more peaks gained than lost accessibility. Effects per stage were small (|log2FC| mostly below 0.3). The High vs Not AD contrast gave the same ranking with fewer peaks (3,087 in L2/3 IT).

**Functional enrichment.** Few terms were significant and none was shared between subclasses. The clearest was TNFα signaling via NF-κB among peaks losing accessibility in Vip neurons.

**TE families.** Among peaks gaining accessibility, LTR/ERV1 and LINE/L2 were over-represented and LINE/L1 under-represented in L2/3, L4 and L5 IT neurons. Within ERV1, MER52D, MER52C and MER52A were the most enriched subfamilies in L2/3 IT (MER52D about 6-fold), and MER52D and MER52A were also enriched in L4 IT. Among peaks losing accessibility, SINE/MIR was over-represented and Alu under-represented in the IT subclasses.

## 12. Discussion

Accessibility changes with AD progression were concentrated in excitatory IT neurons, consistent with the vulnerability of these populations, although part of this difference is statistical power, since they are also the most abundant subclasses. The small effect sizes suggest a gradual shift rather than the strong opening of a few regions.

The TE results only partly match the mouse model. The enrichment of LTR/ERV1 among gaining peaks agrees with the reactivation of LTR elements described by Liu et al. (2026), but LINE-1, also central in that model, was depleted, and LINE/L2 and MIR are old families that are no longer active. This fits better with changes in TE-derived regulatory regions than with broad TE derepression. The loss of young copies at alignment may also hide part of the LINE-1 signal.

## 13. Conclusions

We built a reproducible workflow that runs on the full SEA-AD MTG snATAC-seq dataset and relates accessibility changes to TEs at the cell-subclass level. Accessibility changes along AD progression mainly in L2/3, L4 and L5 IT neurons, and the peaks that gain accessibility are enriched in LTR/ERV1 elements, particularly MER52 subfamilies, but not in LINE-1. These results support a selective association between AD progression and ERV-derived regions in excitatory neurons, not direct evidence of TE reactivation. Next steps are to compare these loci with those reported in the mouse model, control the family-level test for accessibility and GC, and check TE expression in the matching SEA-AD snRNA-seq data.

## 14. Availability, license and citation

**Code:** https://github.com/rafadmtz/NeuroEpigenomics_AD  
**Data:** Not redistributed; obtained with `python src/00_download_data.py` from the public sources in section 8.1.  
**Code license:** MIT (see LICENSE). Source data keep their own terms of use.  
**How to cite:** Díaz Martínez, R., Martínez Enciso, S. and González Morales, S. N. (2026). *NeuroEpigenomics_AD: Epigenomic dynamics during Alzheimer's disease progression*. GitHub repository. https://github.com/rafadmtz/NeuroEpigenomics_AD (metadata in docs/CITATION.cff).  
**Version or release:** [URL]

## 15. References

1. Gabitto, M. I. et al. (2024). Integrated multimodal cell atlas of Alzheimer's disease. *Nature Neuroscience*, 27, 2366–2383. https://doi.org/10.1038/s41593-024-01774-5
2. Liu, et al. (2026). Tau aggregates cause reactivation of transposable DNA elements, leading to Z-RNA–ZBP1-mediated neuronal death. *Nature Neuroscience*. https://doi.org/10.1038/s41593-026-02299-9
