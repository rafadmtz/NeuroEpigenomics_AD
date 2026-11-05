# Epigenomic dynamics of transposable elements during Alzheimer's disease progression: an snATAC-seq analysis of SEA-AD

## General information

| Item | Information |
|:--|:--|
| Member | [Rafael Díaz Martínez — rafadiaz@lcg.unam.mx] |
| Member | [Santiago Martínez Enciso — santiago.marenc@gmail.com] |
| Member | [Sixto Netzahualcóyotl González Morales — sixtogm@lcg.unam.mx] |
| Creation date | [25/08/2026] |
| Last updated | [20/10/2026] |
| Repository | [https://github.com/rafadmtz/NeuroEpigenomics_AD] |
| Status | [In development: preprocessing pipeline implemented (scripts 00–04); the intersection with transposable elements and the differential accessibility analysis are still pending] |

## Project summary

<!-- HELP: Explain in one paragraph the problem, its relevance, the data, the
proposed solution and the main result. It is best to update it at the end.
EXAMPLE: E. coli genomes available at NCBI will be analyzed to identify and
compare resistance genes. A reproducible Python workflow will be developed that
retrieves the data and generates comparative tables and visualizations. -->

In tauopathy models and in the brains of patients with Alzheimer's disease (AD), tau aggregation destabilizes heterochromatin and reactivates normally silenced transposable elements (TEs). These TEs produce double-stranded RNA that activates the Z-RNA–ZBP1 pathway and neuronal death (Liu et al., 2026). This project explores whether this loss of silencing is reflected in chromatin accessibility at single-nucleus resolution. To do so, it uses snATAC-seq data from the middle temporal gyrus (MTG) of the SEA-AD atlas (Gabitto et al., 2024) and the RepeatMasker repeat annotation for hg38. A reproducible Python workflow is being developed that: downloads the data; integrates and validates nucleus and donor metadata; groups nuclei by neuropathological stage (ADNC) and cell subclass; generates per-donor pseudobulk profiles for a differential accessibility analysis; and produces visualizations of cell composition and UMAP embeddings of the L2/3 IT, Sst and Vip subclasses. The next step is to intersect the accessibility peaks with TE coordinates.

## 1. Context and background

<!-- HELP: Present the information needed to understand the project. Define
biological and computational concepts, describe what is known and cite related
work, data or tools.
QUESTIONS: What is the phenomenon of interest? What does the reader need to
know? What methods or tools have been used before?
EXAMPLE: Antimicrobial resistance is a public health problem. Although
specialized databases exist, comparing several genomes requires integrating
data from different files. -->

**Transposable elements and heterochromatin.** About half of the human genome is made up of sequences derived from transposable elements: LINE, SINE and LTR retrotransposons (including endogenous retroviruses), SVA elements and DNA transposons. In somatic cells, most of them remain silenced in constitutive heterochromatin marked by H3K9me3 and bound by HP1α.

**Tau and TE reactivation.** In the PS19 tauopathy model and in the brains of patients with AD, tau aggregation destabilizes this heterochromatin: it displaces HP1α from nucleosomes and erodes the nuclear lamina. The loss of condensation reactivates normally silenced TEs, mainly LINEs and LTRs. Their double-stranded transcripts adopt a Z conformation (Z-RNA) and activate ZBP1, which triggers neuronal necroptosis (Liu et al., 2026).

**snATAC-seq.** Single-nucleus assay for transposase-accessible chromatin sequencing (snATAC-seq) measures which regions of the genome are open in each nucleus. The output is a sparse nuclei × accessibility peaks matrix. Because chromatin opening precedes or accompanies transcription, increased accessibility over TEs can serve as an indirect indicator of their derepression.

**SEA-AD atlas.** The Seattle Alzheimer's Disease Brain Cell Atlas (Gabitto et al., 2024) profiled the MTG of donors spanning the full neuropathological spectrum of AD, with hierarchical cell annotation (class, subclass, supertype). Each donor is characterized by neuropathological measures (ADNC, Braak, Thal, CERAD) and a Continuous Pseudo-progression Score (CPS) that orders them along a progression continuum. This makes it possible to place epigenomic changes within specific cell types.

**Prior tools.** Single-cell data are handled with the AnnData/scanpy ecosystem (Wolf et al., 2018). For ATAC, SnapATAC2 (Zhang et al., 2024) provides spectral dimensionality reduction suited to binary, sparse matrices. To compare conditions across donors, the recommended approach is to aggregate counts per donor (pseudobulk) and model them with DESeq2 or edgeR. This avoids the pseudoreplication of nucleus-level analyses (Squair et al., 2021). TE coordinates are obtained from the UCSC rmsk table, generated with RepeatMasker and Repbase.

## 2. Problem statement

<!-- HELP: Describe the difficulty, need or knowledge gap you want to address,
who it affects and what happens if it is not solved. Do not confuse the problem
with the tool or with the solution.
EXAMPLE: Manual identification of resistance genes across several genomes is
slow, error-prone and hard to reproduce. -->

The evidence that tau reactivates TEs comes mainly from animal models and from transcriptomic or bulk-tissue analyses. It is not known whether TE derepression leaves a detectable footprint on chromatin accessibility in humans, in which cell types it occurs, or at what point in disease progression it appears. The SEA-AD snATAC-seq data could answer this question, but their size (an ~18 GB .h5ad object with hundreds of thousands of nuclei), the multi-donor structure of the design and the need to intersect them with an external repeat annotation make a manual, statistically sound and reproducible analysis difficult.

## 3. Rationale

<!-- HELP: Explain why the project is worth doing, its biological, scientific,
technical or social relevance and who could benefit.
EXAMPLE: An automated workflow will reduce errors and allow the analysis to be
repeated with the same data, parameters and software versions. -->

From a biological standpoint, locating possible TE opening in specific neuronal subclasses and along the progression continuum would help assess whether the tau → TEs → Z-RNA–ZBP1 mechanism described in models has a correlate in human tissue. It would also help identify which populations might be most vulnerable. From a technical standpoint, an automated, version-controlled workflow allows the analysis to be repeated with the same data (dated versions of SEA-AD and rmsk), parameters and dependencies (environment.yml). It also handles the data volume on a server without loading the full matrix into memory. Groups studying the epigenomics of neurodegeneration could benefit, as could anyone wishing to reuse the workflow with other brain regions or SEA-AD modalities.

## 4. General objective

<!-- HELP: State the overall result using an infinitive verb. It must be
achievable within the semester.
EXAMPLE: To develop a reproducible Python workflow to identify and compare
resistance genes in a set of E. coli genomes. -->

To characterize chromatin accessibility changes in neuronal populations along Alzheimer's disease progression, and to evaluate their relationship with regions containing transposable elements, as an indirect approach to their possible reactivation.

## 5. Research questions

<!-- HELP: Formulate biological or computational questions that can be answered
with the available data and methods. State what evidence would be needed.
EXAMPLE: Which resistance genes appear in each genome? Evidence: a
presence/absence table obtained from the annotations. -->

### Question 1

**Question:** [Write the question.]  
**Required evidence:** [Data or results that will allow it to be answered.]

### Question 2

**Question:** [Write the question.]  
**Required evidence:** [Data or results that will allow it to be answered.]

## 6. Scope and limitations

<!-- HELP: Define the organisms, samples, data, analyses and expected result.
State what will not be addressed and any time, computing, access or quality
constraints.
EXAMPLE: At most 20 complete RefSeq genomes will be analyzed. No clinical data
will be used and no experimental validation will be performed. -->

### Included

- Human snATAC-seq data from the middle temporal gyrus (MTG) of SEA-AD, using the object already filtered by the authors (final-nuclei, version 2024-12-06).
- The four AD neuropathological change stages: Not AD, Low, Intermediate and High.
- Detailed analysis of three subclasses: L2/3 IT (excitatory), Sst and Vip (inhibitory). Grouping and pseudobulk are generated for all subclasses.
- TE annotation from the UCSC rmsk table (hg38), restricted to the main chromosomes.
- Pseudobulk aggregation by donor and subclass, and differential accessibility between stages.
- Quality control, cell composition and UMAP figures.

### Not included

- Other SEA-AD brain regions (e.g. dorsolateral prefrontal cortex) or other modalities (snRNA-seq, multiome, spatial transcriptomics).
- De novo peak calling or read realignment: the peaks and nuclei defined by SEA-AD are used.
- Quantification of TE expression or detection of polymorphic insertions per donor.
- Experimental validation.

### Known limitations

- The UCSC rmsk annotation uses RepeatMasker 3.0.1 and Repbase from 2010, and has not been recomputed. Classification of young subfamilies (L1HS, AluY, SVA) is less accurate than with Dfam.
- Recent TE copies are nearly identical to each other; short reads falling on them are often discarded during alignment. These TEs may appear closed even when they are not.
- rmsk describes the reference genome, not each donor's own insertions.
- Each rmsk row is a fragment, not necessarily a complete element.
- Accessibility is an indirect indicator of TE reactivation, not a measure of their transcription.
- With the ≥10 nuclei per donor and subclass filter (MIN_NUCLEI), some subclass–stage combinations may end up with few replicates.
- The .h5ad object is ~18 GB. The workflow requires a server (the pseudobulk step is calibrated for ~64 GB of RAM).
- There are no automated tests yet.

## 7. Proposed solution

<!-- HELP: Describe the product or strategy that could solve the problem. This
is an initial proposal and may change. Explain its components, not just the
technologies.
EXAMPLE: A modular program will receive identifiers, download files, extract
genes, store results and generate visualizations. -->

A modular workflow of numbered scripts in src/ that are run in order from the repository root. Each one can be run from the terminal or imported as a module through its main() function. Components:

1. Download (00_download_data.py): retrieves the snATAC-seq .h5ad object, the all-nuclei metadata CSV and the rmsk table, with versions pinned in the code.
2. Metadata integration and validation (01_metadata_analysis.py): merges the metadata from the .h5ad and the CSV, checks that identifiers match and that each donor has a unique clinical profile, and defines the severity groups. It also builds a per-donor table and exports the peak coordinates in BED format.
3. Cohort description (01b_filtered_data_graphs.py): plots subclass composition by severity and supertype composition by donor.
4. Matrix partitioning (02_group_atac_by_condition_and_subclass.py): writes one .h5ad per ADNC stage and subclass combination, in parallel.
5. Pseudobulk aggregation (03_pseudobulk_by_donor.py): sums counts by donor and subclass, reading the matrix in chunks under a memory limit. It generates the sample tables (colData) for the differential analysis and quality control figures.

### 7.1 Expected result or product

<!-- HELP: State the concrete deliverable: program, package, analysis workflow,
database, visualizations or other product.
EXAMPLE: Executable repository with scripts, test data, documentation,
comparative table and reproducible figures. -->

An executable, MIT-licensed repository including: the workflow scripts; the Conda environment (environment.yml); the documentation (README, this report, CITATION.cff); per-nucleus and per-donor metadata tables; per-subclass pseudobulk matrices ready for differential analysis; a table of TE-annotated peaks; differential accessibility results; and reproducible figures (cell composition, pseudobulk quality control, UMAP).

## 8. Data

### 8.1 Data sources

<!-- HELP: Include institution, database, URL, identifier, version or access
date and terms of use. Not every project will use NCBI.
EXAMPLE: NCBI RefSeq, GCF_000005845.2, accessed on dd/mm/yyyy. -->

| Source | Identifier or version | URL | Access date | License or terms |
|:--|:--|:--|:--|:--|
| [SEA-AD (Allen Institute), public AWS S3 bucket sea-ad-single-cell-profiling] | [SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad (version 2024-12-06)] | [https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad] | [15/10/2026] | [Allen Institute Terms of Use (https://alleninstitute.org/legal/terms-use/); cite Gabitto et al. (2024)] |
| [SEA-AD (Allen Institute), public AWS S3 bucket sea-ad-single-cell-profiling] | [SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv (version 2024-12-06)] | [https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq/SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv] | [15/10/2026] | [Allen Institute Terms of Use (https://alleninstitute.org/legal/terms-use/); cite Gabitto et al. (2024)] |
| [UCSC Genome Browser, hg38 rmsk table] | [RepeatMasker 3.0.1 + Repbase 2010-03-02 (last updated at UCSC: 18/10/2022)] | [https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz] | [15/10/2026] | [Free for academic use according to UCSC; cite RepeatMasker and Repbase] |

### 8.2 Data characteristics

<!-- HELP: Describe organisms, samples, variables, formats, versions, size and
other attributes needed to interpret the data.
EXAMPLE: FASTA and GFF3 files of 20 complete E. coli genomes. -->

- Organism and tissue: *Homo sapiens*, middle temporal gyrus (MTG), post-mortem tissue from SEA-AD atlas donors spanning the full spectrum of AD neuropathological change, including neurotypical reference donors.
- snATAC-seq (.h5ad, ~18 GB): AnnData object. X is a sparse nuclei × peaks count matrix; obs contains per-nucleus metadata (donor, hierarchical cell annotation, neuropathological and clinical variables, quality metrics); var contains the peaks, whose identifier encodes their coordinates. Nuclei were already filtered by the authors (final-nuclei). Coordinates are in GRCh38, with chromosome names chr1, chr2…
- Metadata (.csv): one row per nucleus, including those removed by SEA-AD quality control (all-nuclei). The workflow keeps only the nuclei present in the .h5ad.
- RepeatMasker (rmsk.txt.gz, ~152 MB compressed): 5,683,690 repeats in 17 tab-separated columns. Coordinates are 0-based, half-open and unsorted. The hgdownload version has no header. Keeping the LINE, SINE, LTR, DNA and Retroposon classes on chr1–chr22, chrX and chrY leaves 4,505,027 TE fragments.

### 8.3 Data organization

<!-- HELP: Show the planned structure. Do not upload sensitive data, tokens,
passwords or large files. Use .gitignore and document how to obtain anything
not stored in Git.
EXAMPLE: data/raw keeps the originals and data/processed the derived files. -->

data/, logs/ and results/ are listed in .gitignore: data and results are generated locally by the scripts and are not version-controlled.

```text
NeuroEpigenomics_AD/
├── data/                                   # Not version-controlled
│   ├── SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad    # 00
│   ├── rmsk.txt.gz                                       # 00
│   ├── metadata/
│   │   └── SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv  # 00
│   └── processed/
│       ├── <Condition>_<Subclass>.h5ad                   # 02
│       ├── metadata_test.csv, peaks.txt                  # 02
│       └── pseudobulk/
│           ├── pseudobulk_<Subclass>.h5ad                # 03
│           ├── colData_<Subclass>.csv                    # 03
│           └── pseudobulk_summary.csv                    # 03
├── docs/
│   ├── reporte-proyecto.md
│   └── CITATION.cff
├── logs/                                   # Not version-controlled (nohup output)
├── notebooks/
├── results/                                # Not version-controlled
│   ├── 01_metadata_output/                               # 01
│   ├── 01b_graphs/                                     # 01b
│   └── figures/                                          # 03 (pseudobulk/) and 04
├── src/                                    # Scripts 00–04
├── tests/
├── environment.yml
├── LICENSE
└── README.md
```

### 8.4 Data dictionary or format

<!-- HELP: Describe relevant fields or columns. Include small snippets when they
help explain the format, but not complete files.
EXAMPLE: In GFF3, seqid identifies the sequence; type indicates gene, CDS, etc. -->

| File or dataset | Field/column | Type | Description | Values or units |
|:--|:--|:--|:--|:--|
| [File] | [Field] | [Type] | [Description] | [Values] |

## 9. Methodology

<!-- HELP: This section will evolve. First describe the plan and then update it
with what you actually ran, including parameters and decisions. -->

### 9.1 Analysis or development stages

<!-- HELP: Describe the sequence from data retrieval to result validation. Link
each stage to a question or objective.
EXAMPLE: download, validation, transformation, analysis, visualization and
testing. -->

1. Data retrieval (00_download_data.py). Downloads with wget -c (or with urllib if wget is not installed) and skips files that already exist. General objective.
2. Metadata integration and validation (01_analisis_metadatos.py). Reads the .h5ad in backed='r' mode; checks the intersection of identifiers with the CSV; merges the metadata, giving priority to the .h5ad annotation; defines Severity_Group; classifies quality metrics as nucleus-level or library-level; builds a per-donor table (fails if a donor has more than one clinical profile); and exports atac_peaks.bed.
3. Cohort description (01b_graficos_datos_filtrados.py). Generates stacked bar plots of nuclei by subclass and severity. It also plots the proportion of supertypes per donor, ordered by CPS, for L2/3 IT, Sst and Vip: supertypes below 3% are grouped into "Other" and donors with fewer than 50 nuclei are flagged.
4. Partitioning by stage and subclass (02_group_atac_by_condition_and_subclass.py). Writes one compressed .h5ad per ADNC × subclass combination, using 6 parallel processes (joblib).
5. Pseudobulk aggregation (03_pseudobulk_by_donor.py). Sums counts by (subclass, donor), reading the matrix in chunks of 20,000 nuclei through a sparse indicator matrix. Subclasses are processed in batches so that the accumulator does not exceed 8 GB, and donors with fewer than 10 nuclei are discarded. Exports the matrices and colData, plus three quality control figures.
6. Unsupervised exploration (04_umap_atac_by_subclass.py). Per subclass and for all subclasses combined: removes nuclei and peaks with zero counts, selects 50,000 peaks (snap.pp.select_features), computes the spectral embedding (snap.tl.spectral), the 30-neighbor graph and the UMAP (seed 2). Colors by total counts, supertype, subclass and condition.

### 9.2 Tools and technologies

<!-- HELP: Record languages, libraries and programs with their versions and
purpose. Do not include credentials.
EXAMPLE: Python 3.x; Biopython to read biological formats; Seaborn for
visualization. -->

| Tool | Version | Purpose |
|:--|:--|:--|
| [Tool] | [Version] | [Use] |

### 9.3 Validation strategy

<!-- HELP: Explain how you will check the code and results: unit tests, known
data, comparison with another tool or manual review.
EXAMPLE: Five known annotations will be compared, and valid inputs, nonexistent
identifiers and incomplete files will be tested. -->

[Describe the validation here.]

## 10. Work plan

### 10.1 Division of responsibilities

<!-- HELP: Define initial responsibilities without isolating each member. Every
major contribution must be reviewed by another person through a Pull Request.
EXAMPLE: Ana develops the download module and reviews the visualization module. -->

| Member | Main responsibility | Review responsibility |
|:--|:--|:--|
| [Name] | [Responsibility] | [What or whom they will review] |

### 10.2 Risks and alternatives

<!-- HELP: Identify situations that could prevent or delay the project and
define an alternative.
EXAMPLE: The data require too much storage; alternative: reduce the number of
genomes using documented criteria. -->

| Risk | Probability | Impact | Prevention or alternative |
|:--|:--|:--|:--|
| [Risk] | Low/Medium/High | Low/Medium/High | [Action] |

## 11. Results

<!-- HELP: Present results linked to the questions and objectives. Include
tables or figures with titles, captions and source files. Describe here what was
obtained; interpret it in the Discussion.
EXAMPLE: Presence/absence table generated by src/compare_genes.py.
FIRST SESSION: leave this section as pending. -->

> Status: pending. It will be completed during development.

## 12. Discussion

<!-- HELP: Interpret the results, explain whether they answer the questions,
compare them with the background and point out limitations. Do not just repeat
the values.
EXAMPLE: The observed distribution suggests..., although interpretation is
limited by the quality of the annotations. -->

> Status: pending. It will be completed after results are obtained.

## 13. Conclusions

<!-- HELP: Summarize what was learned, which questions were answered and whether
the objectives were met. Include contributions, limitations and future work. Do
not introduce new results.
EXAMPLE: The workflow made it possible to identify..., but it will be necessary
to incorporate... -->

> Status: pending. It will be completed at the end of the project.

## 14. Availability, license and citation

<!-- HELP: State where the code is, under which license it can be reused and how
to cite it. Link this section to LICENSE, CITATION.cff, codemeta.json, the final
release and, where applicable, a DOI.
EXAMPLE: Code on GitHub under MIT; citation available in CITATION.cff. -->

**Code:** [https://github.com/rafadmtz/NeuroEpigenomics_AD]  
**Data:** [Not redistributed in the repository. They are obtained with python src/00_download_data.py from the public SEA-AD AWS S3 bucket (sea-ad-single-cell-profiling, version 2024-12-06) and from UCSC (hg38 rmsk); see section 8.1.]  
**Code license:** [MIT (© 2026 rafadmtz; see LICENSE). The SEA-AD and UCSC/RepeatMasker data retain their own terms of use.]  
**How to cite:** [docs/CITATION.cff (CFF 1.2.0; release date 19/10/2026):
Díaz Martínez, R., Martínez Enciso, S. and González Morales, S. N. (2026). NeuroEpigenomics_AD: Epigenomic dynamics during Alzheimer's disease progression. GitHub repository. https://github.com/rafadmtz/NeuroEpigenomics_AD
]  
**Version or release:** [URL]

## 15. References

<!-- HELP: Record publications, data, software and documents consulted in a
consistent format. Include a DOI, URL or persistent identifiers. Every citation
in the text must appear here.
EXAMPLE: Blattner, F. R. et al. (1997). The complete genome sequence of
Escherichia coli K-12. Science, 277(5331), 1453–1462.
https://doi.org/10.1126/science.277.5331.1453 -->

1. [Reference 1]
2. [Reference 2]

---

<!-- GUIDANCE FOR THE FIRST TWO SESSIONS:
Complete General information, the provisional Summary, sections 1 to 7,
preliminary data sources and the work plan. Methodology, Results, Discussion,
Conclusions and Availability will evolve during the semester. Replace the
bracketed prompts with the team's content. -->