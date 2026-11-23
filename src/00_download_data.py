#!/usr/bin/env python3
"""
00_download_data.py — Data and metadata download for the NeuroEpigenomics_AD project
==========================================================================

Downloads to the server:
  1. Pre-filtered snATAC-seq dataset from AD patients (middle temporal gyrus) from SEA-AD, in h5ad format,
     with filters established by the dataset owners. Used in steps 01, 02 and 03.
  2. SEA-AD metadata (middle temporal gyrus) and controls, in CSV format.
     Used in step 01.
  3. Transposable elements annotation from RepeatMasker (UCSC hg38), in
     the same coordinates (GRCh38, chromosomes chr1, chr2...) as SEA-AD.
     Used in step 05.
  4. Reference genome sequence (UCSC hg38 FASTA), in the same coordinates
     (GRCh38, chromosomes chr1, chr2...) as SEA-AD. Used to compute the GC
     content of the accessibility peaks with bedtools nuc (steps 07 and 08).
  5. Gene annotation from GENCODE (basic annotation, reference chromosomes
     only), in GTF format and in the same coordinates (GRCh38, chromosomes
     chr1, chr2...) as SEA-AD. Used to assign each peak to its nearest gene
     TSS and to classify peaks as promoter, proximal or distal (step 08).
  6. Gene set collections from MSigDB (human, gene symbols), in GMT format:
     Hallmark and GO Biological Process. Used for the functional enrichment
     of differentially accessible peaks (step 08).

SEA-AD Source: AWS S3 public bucket "sea-ad-single-cell-profiling"
(Gabitto et al. 2024, Nat Neurosci, "Data availability" section).
GENCODE Source: EBI FTP (Frankish et al. 2023, Nucleic Acids Res).
MSigDB Source: Broad Institute (Liberzon et al. 2015, Cell Systems). 

Usage:
  The file can be executed directly from the command line, or
  imported as a module to call the main() function.

"""

#==================================================
# Libraries used
#==================================================

import os
import shutil
import subprocess
import urllib.error
import urllib.request

#==================================================
# Variable definitions
#==================================================

URL_SEAAD_ATAC_SEQ= "https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad"
VERSION_SEAAD= "2024-12-06"

URL_SEAAD_METADATA= "https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq/SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv"
VERSION_SEAAD_METADATA= "2024-12-06"

URL_REPEATMASKER= "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
VERSION_RMSK= "RepeatMasker 3.0.1, RepBase 2010-03-02, UCSC hg38 rmsk table"

URL_HG38_FASTA= "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
VERSION_HG38= "UCSC hg38 (GRCh38) genome FASTA"

VERSION_GENCODE= "50"
URL_GENCODE_GTF= (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    f"release_{VERSION_GENCODE}/gencode.v{VERSION_GENCODE}.basic.annotation.gtf.gz"
)

VERSION_MSIGDB= "2026.2.Hs"
URL_MSIGDB= f"https://data.broadinstitute.org/gsea-msigdb/msigdb/release/{VERSION_MSIGDB}"
# Collections used in step 08 (gene symbols): Hallmark and GO Biological Process.
MSIGDB_COLLECTIONS= ["h.all", "c5.go.bp"]

DATA_OUTPUT_DIR = "data/"

METADATA_OUTPUT_DIR = "data/metadata/"

GENOME_OUTPUT_DIR = "data/genome/"

ANNOTATION_OUTPUT_DIR = "data/annotation/"

GENESETS_OUTPUT_DIR = "data/genesets/"

#==================================================
# Functions
#==================================================

def download(url, destination):
    """Downloads a file. Uses wget -c (resumes if interrupted) and, if not
    installed, urllib.request.urlretrieve."""
    os.makedirs(os.path.dirname(destination), exist_ok=True)

    if shutil.which("wget"):
        subprocess.run(["wget", "-c", "-q", "--show-progress", "-O", destination, url],
                       check=True)
    else:
        print("wget is not installed; using urllib (will not resume downloads)")
        urllib.request.urlretrieve(url, destination)

#==================================================
# Main
#==================================================
def main():

    print("Starting data download...")
    seaad_filename = os.path.join(DATA_OUTPUT_DIR, f"SEAAD_MTG_ATACseq_final-nuclei.{VERSION_SEAAD}.h5ad")
    seaad_metadata_filename = os.path.join(METADATA_OUTPUT_DIR, f"SEAAD_MTG_ATACseq_all-nuclei_metadata.{VERSION_SEAAD_METADATA}.csv")
    repeatmasker_filename = os.path.join(DATA_OUTPUT_DIR, f"rmsk.txt.gz")
    hg38_filename_gz = os.path.join(GENOME_OUTPUT_DIR, "hg38.fa.gz")
    hg38_filename = os.path.join(GENOME_OUTPUT_DIR, "hg38.fa")
    # The GTF is kept compressed; step 08 reads it directly from the .gz file.
    gencode_filename = os.path.join(ANNOTATION_OUTPUT_DIR, f"gencode.v{VERSION_GENCODE}.basic.annotation.gtf.gz")

    try:
        if not os.path.exists(seaad_filename):
            print(f"Downloading filtered SEA-AD ATAC-seq dataset to {seaad_filename}...")
            download(URL_SEAAD_ATAC_SEQ, seaad_filename)
        if not os.path.exists(seaad_metadata_filename):
            print(f"Downloading SEA-AD metadata to {seaad_metadata_filename}...")
            download(URL_SEAAD_METADATA, seaad_metadata_filename)
        if not os.path.exists(repeatmasker_filename):
            print(f"Downloading RepeatMasker annotation to {repeatmasker_filename}...")
            download(URL_REPEATMASKER, repeatmasker_filename)
        if not os.path.exists(hg38_filename):
            print(f"Downloading hg38 reference genome to {hg38_filename}...")
            download(URL_HG38_FASTA, hg38_filename_gz)
            subprocess.run(["gunzip", hg38_filename_gz], check=True)
        if not os.path.exists(gencode_filename):
            print(f"Downloading GENCODE v{VERSION_GENCODE} gene annotation to {gencode_filename}...")
            download(URL_GENCODE_GTF, gencode_filename)
        for collection in MSIGDB_COLLECTIONS:
            gmt_name = f"{collection}.v{VERSION_MSIGDB}.symbols.gmt"
            gmt_filename = os.path.join(GENESETS_OUTPUT_DIR, gmt_name)
            if not os.path.exists(gmt_filename):
                print(f"Downloading MSigDB {collection} gene sets to {gmt_filename}...")
                download(f"{URL_MSIGDB}/{gmt_name}", gmt_filename)
    except (subprocess.CalledProcessError, urllib.error.URLError) as e:
        print(f"Download error: {e}")
        print("Check your internet connection or the file URL.")


if __name__ == "__main__":
    main()