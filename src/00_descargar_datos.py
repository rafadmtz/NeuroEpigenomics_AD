#!/usr/bin/env python3
"""
00_descarga_de_datos.py — Descarga de datos y metadatos para el proyecto NeuroEpigenomics_AD
==========================================================================

Descarga al servidor:
  1. Dataset ya filtrado de snATAC-seq de pacientes de AD(giro temporal medio) de SEA-AD, en formato h5ad,
     con los filtros establecidos por los dueños del dataset.
  2. Metadatos de SEAD-AD (giro temporal medio) y controles, en formato CSV.
  3. La anotación de elementos transponibles de RepeatMasker (UCSC hg38), en
     las mismas coordenadas (GRCh38, cromosomas chr1, chr2...) que SEA-AD.

Fuente de SEA-AD: bucket público de AWS S3 "sea-ad-single-cell-profiling"
(Gabitto et al. 2024, Nat Neurosci, sección "Data availability").

Uso: 
  El archivo se puede ejecutar directamente desde la línea de comandos, o bien
  importarse como módulo y llamar a la función main().


"""

#==================================================
#Liberias utilizadas
#==================================================

import os
import shutil
import subprocess
import urllib.request

#==================================================
#Definición de variables
#==================================================

URL_SEAAD_ATAC_SEQ= "https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad"
VERSION_SEAAD= "2024-12-06" 

URL_SEAAD_METADATA= "https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq/SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv"
VERSION_SEAAD_METADATA= "2024-12-06"

URL_REPEATMASKER= "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
VERSION_RMSK= "RepeatMasker 3.0.1, RepBase 2010-03-02, tabla rmsk de UCSC hg38"

DIR_SALIDA_DATA = "data/"

DIR_SALIDA_METADATA = "data/metadata/"
#==================================================
#Funciones
#==================================================

def descargar(url, destino):
    """Descarga un archivo. Usa wget -c (reanuda si se corta) y, si no está
    instalado, urllib.request.urlretrieve."""
    os.makedirs(os.path.dirname(destino), exist_ok=True)

    if shutil.which("wget"):
        subprocess.run(["wget", "-c", "-q", "--show-progress", "-O", destino, url],
                       check=True)
    else:
        print("wget no está instalado; usando urllib (no reanuda descargas)")
        urllib.request.urlretrieve(url, destino)

#==================================================
#Main
#==================================================
def main():

    print("Iniciando descarga de datos...")
    nombre_archivo_seaad = os.path.join(DIR_SALIDA_DATA, f"SEAAD_MTG_ATACseq_final-nuclei.{VERSION_SEAAD}.h5ad")
    nombre_archivo_seaad_metadata = os.path.join(DIR_SALIDA_METADATA, f"SEAAD_MTG_ATACseq_all-nuclei_metadata.{VERSION_SEAAD_METADATA}.csv")
    nombre_archivo_repeatmasker = os.path.join(DIR_SALIDA_DATA, f"rmsk.txt.gz")
    try:
        if not os.path.exists(nombre_archivo_seaad):
            print(f"Descargando SEA-AD ATAC-seq dataset filtrado a {nombre_archivo_seaad}...")
            descargar(URL_SEAAD_ATAC_SEQ, nombre_archivo_seaad)
        if not os.path.exists(nombre_archivo_seaad_metadata):
            print(f"Descargando metadatos de SEA-AD a {nombre_archivo_seaad_metadata}...")
            descargar(URL_SEAAD_METADATA, nombre_archivo_seaad_metadata)
        if not os.path.exists(nombre_archivo_repeatmasker):
            print(f"Descargando anotación de RepeatMasker a {nombre_archivo_repeatmasker}...")
            descargar(URL_REPEATMASKER, nombre_archivo_repeatmasker)
    except subprocess.CalledProcessError as e:
        print(f"Error al descargar: {e}")
        print("Verifica tu conexión a internet o la URL del archivo.")


if __name__ == "__main__":
    main()