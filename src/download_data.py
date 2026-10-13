#!/usr/bin/env python3
"""
download_data.py — Descarga de datos para el proyecto NeuroEpigenomics_AD
==========================================================================

Descarga al servidor:
  1. Objetos de snATAC-seq por donador del atlas SEA-AD (giro temporal medio),
     con la matriz núcleos x picos y la metadata de cada núcleo.
  2. BigWigs de accesibilidad agregada por categoría ADNC x subclase celular.
  3. La anotación de elementos transponibles de RepeatMasker (UCSC hg38), en
     las mismas coordenadas (GRCh38, cromosomas chr1, chr2...) que SEA-AD.

Fuente de SEA-AD: bucket público de AWS S3 "sea-ad-single-cell-profiling"
(Gabitto et al. 2024, Nat Neurosci, sección "Data availability").

Uso:
  1. Editar la sección CONFIGURACIÓN.
  2. Correr primero con SOLO_REVISAR = True para ver qué existe y cuánto pesa.
  3. Cambiar a SOLO_REVISAR = False y corre:  python download_data.py
     Si una descarga se interrumpe, volver a correrlo: reanuda donde se quedó.

Requisitos: pandas. Recomendado: wget (permite reanudar descargas).
"""

import os
import shutil
import subprocess
import urllib.request

import pandas as pd


# =============================================================================
# CONFIGURACIÓN
# =============================================================================

DIR_SALIDA = "data/raw"     # carpeta base (data/ está en .gitignore)
SOLO_REVISAR = False         # True: solo verifica existencia y tamaños

URL_SEAAD = "https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/ATACseq"
VERSION = "2024-12-06"      # versión de los objetos en el bucket

URL_RMSK = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"

# --- Donadores (objetos de ~20–620 MB cada uno) ------------------------------
# Lista de donadores con objeto en el bucket (consultada en oct-2026).
TODOS_LOS_DONADORES = [
    "H19.30.002", "H19.30.004", "H19.33.004",
    "H20.33.001", "H20.33.002", "H20.33.004", "H20.33.005", "H20.33.008",
    "H20.33.011", "H20.33.012", "H20.33.013", "H20.33.014", "H20.33.015",
    "H20.33.016", "H20.33.017", "H20.33.018", "H20.33.019", "H20.33.020",
    "H20.33.024", "H20.33.025", "H20.33.026", "H20.33.027", "H20.33.028",
    "H20.33.029", "H20.33.030", "H20.33.031", "H20.33.032", "H20.33.033",
    "H20.33.034", "H20.33.035", "H20.33.036", "H20.33.037", "H20.33.038",
    "H20.33.039", "H20.33.040", "H20.33.041", "H20.33.043", "H20.33.044",
    "H20.33.045", "H20.33.046",
    "H21.33.001", "H21.33.002", "H21.33.003", "H21.33.004", "H21.33.005",
    "H21.33.006", "H21.33.007", "H21.33.008", "H21.33.009", "H21.33.010",
    "H21.33.011", "H21.33.012", "H21.33.013", "H21.33.014", "H21.33.015",
    "H21.33.016", "H21.33.017", "H21.33.018", "H21.33.019", "H21.33.020",
    "H21.33.021", "H21.33.022", "H21.33.023", "H21.33.025", "H21.33.026",
    "H21.33.027", "H21.33.028", "H21.33.029", "H21.33.030", "H21.33.031",
    "H21.33.032", "H21.33.033", "H21.33.034", "H21.33.035", "H21.33.036",
    "H21.33.037", "H21.33.038", "H21.33.039", "H21.33.040", "H21.33.041",
    "H21.33.042", "H21.33.043", "H21.33.044", "H21.33.045", "H21.33.046",
    "H21.33.047",
]
DONADORES = TODOS_LOS_DONADORES   # prueba; para todo: TODOS_LOS_DONADORES

# --- BigWigs (señal agregada por ADNC x subclase) ----------------------------
# Subclases: Astro, Chandelier, Endo, L2-3IT, L4IT, L5-6NP, L5ET, L5IT, L6CT,
# L6IT, L6ITCar3, Lamp5, Lamp5_Lhx6, Micro-PVM, OPC, Oligo, Pax6, Pvalb, Sncg,
# Sst, VLMC, Vip.   (L2-3IT pesa 2.4–4.9 GB por archivo)
ADNC = ["ADNC0", "ADNC3"]
SUBCLASES = ["L2-3IT"]


# =============================================================================
# 1. Lista de archivos
# =============================================================================

def construir_lista_archivos():
    """Regresa una lista de diccionarios {tipo, url, destino} según la configuración."""
    archivos = []

    for donador in DONADORES:
        nombre = f"{donador}_SEAAD_MTG_ATACseq_final-nuclei.{VERSION}.h5ad"
        archivos.append({
            "tipo": "h5ad",
            "url": f"{URL_SEAAD}/donor_objects/{nombre}",
            "destino": os.path.join(DIR_SALIDA, "ATAC", "donadores", nombre),
        })

    for grupo in ADNC:
        for subclase in SUBCLASES:
            nombre = f"{grupo}{subclase}.bw"
            archivos.append({
                "tipo": "bigwig",
                "url": f"{URL_SEAAD}/bw/{nombre}",
                "destino": os.path.join(DIR_SALIDA, "ATAC", "bw", nombre),
            })

    archivos.append({
        "tipo": "gzip",
        "url": URL_RMSK,
        "destino": os.path.join(DIR_SALIDA, "referencia", "hg38_rmsk.txt.gz"),
    })

    return archivos


# =============================================================================
# 2. Validaciones
# =============================================================================

def tamano_remoto(url):
    """Consulta el tamaño de un archivo en el servidor sin descargarlo.

    Regresa el tamaño en bytes, o None si el archivo no existe o no hay conexión.
    """
    peticion = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(peticion, timeout=60) as respuesta:
            return int(respuesta.headers["Content-Length"])
    except Exception as error:
        print(f"   [no disponible] {url}\n   motivo: {error}")
        return None


# Primeros bytes ("número mágico") que identifican cada formato
FIRMAS = {
    "h5ad": b"\x89HDF\r\n\x1a\n",     # .h5ad es un archivo HDF5
    "bigwig": b"\x26\xfc\x8f\x88",    # 0x888FFC26 en little-endian
    "gzip": b"\x1f\x8b",
}


def validar_formato(ruta, tipo):
    """Comprueba que el archivo descargado tenga el formato esperado."""
    firma = FIRMAS[tipo]
    with open(ruta, "rb") as f:
        inicio = f.read(len(firma))
    if tipo == "bigwig":   # el número mágico puede venir en cualquier orden de bytes
        return inicio in (firma, firma[::-1])
    return inicio == firma


def descarga_completa(ruta, tamano_esperado):
    """True si el archivo ya existe localmente con el tamaño correcto."""
    return os.path.exists(ruta) and os.path.getsize(ruta) == tamano_esperado


# =============================================================================
# 3. Descarga
# =============================================================================

def descargar(url, destino):
    """Descarga un archivo. Usa wget -c (reanuda si se corta) y, si no está
    instalado, urllib.request.urlretrieve."""
    os.makedirs(os.path.dirname(destino), exist_ok=True)

    if shutil.which("wget"):
        subprocess.run(["wget", "-c", "-q", "--show-progress", "-O", destino, url],
                       check=True)
    else:
        print("   wget no está instalado; usando urllib (no reanuda descargas)")
        urllib.request.urlretrieve(url, destino)


def procesar_archivos(archivos):
    """Ciclo principal: valida existencia, descarga y verifica cada archivo.

    Regresa un DataFrame con el resumen de lo que pasó con cada archivo.
    """
    resumen = []

    for i, archivo in enumerate(archivos, start=1):
        nombre = os.path.basename(archivo["destino"])
        print(f"\n[{i}/{len(archivos)}] {nombre}")

        tamano = tamano_remoto(archivo["url"])
        if tamano is None:
            resumen.append({"archivo": nombre, "GB": None, "estado": "no existe en el servidor"})
            continue

        gb = round(tamano / 1e9, 2)
        print(f"   tamaño: {gb} GB")

        if descarga_completa(archivo["destino"], tamano):
            estado = "ya estaba descargado"
        elif SOLO_REVISAR:
            estado = "pendiente (modo revisión)"
        else:
            try:
                descargar(archivo["url"], archivo["destino"])
            except Exception as error:
                print(f"   error al descargar: {error}")
                resumen.append({"archivo": nombre, "GB": gb, "estado": "error en la descarga"})
                continue

            if not descarga_completa(archivo["destino"], tamano):
                estado = "incompleto: vuelve a correr el script para reanudar"
            elif not validar_formato(archivo["destino"], archivo["tipo"]):
                estado = f"formato inválido (no parece {archivo['tipo']})"
            else:
                estado = "descargado y validado"

        print(f"   estado: {estado}")
        resumen.append({"archivo": nombre, "GB": gb, "estado": estado})

    return pd.DataFrame(resumen)


# =============================================================================
# Programa principal
# =============================================================================

def main():
    print("=" * 70)
    print("Descarga de datos" + ("  [MODO REVISIÓN: no se descarga nada]" if SOLO_REVISAR else ""))
    print("=" * 70)

    archivos = construir_lista_archivos()
    resumen = procesar_archivos(archivos)

    print("\n" + "=" * 70)
    print("Resumen")
    print("=" * 70)
    print(resumen.to_string(index=False))
    print(f"\nTotal: {resumen['GB'].sum():.2f} GB")

    os.makedirs(DIR_SALIDA, exist_ok=True)
    ruta_resumen = os.path.join(DIR_SALIDA, "resumen_descargas.csv")
    resumen.to_csv(ruta_resumen, index=False)
    print(f"Resumen guardado en {ruta_resumen}")


if __name__ == "__main__":
    main()
