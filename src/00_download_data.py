#!/usr/bin/env python3
"""
00_download_data.py — Data and metadata download for the NeuroEpigenomics_AD project
==========================================================================

Downloads to the server:
  1. Pre-filtered snATAC-seq dataset from AD patients (middle temporal gyrus) from SEA-AD, in h5ad format,
     with filters established by the dataset owners.
  2. SEA-AD metadata (middle temporal gyrus) and controls, in CSV format.
  3. Transposable elements annotation from RepeatMasker (UCSC hg38), in
     the same coordinates (GRCh38, chromosomes chr1, chr2...) as SEA-AD.

SEA-AD Source: AWS S3 public bucket "sea-ad-single-cell-profiling"
(Gabitto et al. 2024, Nat Neurosci, "Data availability" section).

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

DATA_OUTPUT_DIR = "data/"

METADATA_OUTPUT_DIR = "data/metadata/"

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
    except subprocess.CalledProcessError as e:
        print(f"Download error: {e}")
        print("Check your internet connection or the file URL.")


if __name__ == "__main__":
    main()