#!/usr/bin/env python3
"""
01_analisis_metadatos.py — Analisis de metadatos para el proyecto NeuroEpigenomics_AD
==========================================================================
 
Analiza los metadatos de snATAC-seq del giro temporal medio (MTG) de SEA-AD,
de donadores con Alzheimer (AD) y controles.
 
Fuente de SEA-AD: bucket publico de AWS S3 "sea-ad-single-cell-profiling"
(Gabitto et al. 2024, Nat Neurosci, sección "Data availability").
 
Entradas:
  1.- SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad
      Datos filtrados (solo nucleos finales), se abre en modo lectura (backed='r').
  2.- METADATA/SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv
      Metadatos de todos los nucleos, incluyendo los eliminados por la filtracion.
 
Pasos:
  1.- Carga el h5ad y el csv y verifica que los IDs de los nucleos coincidan.
  2.- Une ambos metadatos. Si una columna esta en los dos, se conserva la anotacion final del h5ad.
  3.- Define el grupo de severidad (Severity_Group) a partir de 'Overall AD neuropathological Change' y 'Neurotypical reference'.
  4.- Revisa si las metricas de calidad de ATAC están calculadas por nucleo o por biblioteca.
  5.- Resume la informacion clinica por donador y verifica que cada donador tenga un solo perfil clinico.
  6.- Extrae las coordenadas de los picos de ATAC en formato BED (cromosomas chr1-22, X, Y y M), para cruzarlas con la tabla rmsk.
 
Salidas (en results/01_metadata_output/):
  1.- metadata_nucleos_completa.csv        Metadatos unidos, una fila por nucleo.
  2.- metadata_donadores_extendida.csv     Variables clinicas por donador y numero de nucleos (n_nuclei), ordenado por Continuous Pseudo-progression Score.
  3.- atac_peaks.bed                       Coordenadas de los picos (chr, start, end, peak_id), ordenadas.
  4.- 00_qc_and_merge_log.txt              Registro de todo lo anterior.
 
Uso:
  El archivo se puede ejecutar directamente desde la linea de comandos, o bien
  importarse como modulo y llamar a la funcion main().
 
 
"""

#==================================================
#Liberias utilizadas
#==================================================
 
import scanpy as sc
import pandas as pd
from pathlib import Path
 
#==================================================
#Definición de variables
#==================================================
 
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
OUT_DIR = PROJECT_ROOT / "results" / "01_metadata_output"

h5ad_path = DATA_DIR / "SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad"
csv_path = DATA_DIR / "metadata" / "SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv"

nuclei_out = OUT_DIR / "metadata_nucleos_completa.csv"
qc_log_path = OUT_DIR / "00_qc_and_merge_log.txt"
 
log_lines = []
 
#==================================================
#Funciones
#==================================================
 
def log_print(msg):
    """Guarda lo impreso en pantalla en un log_lines para llevar un registro de corridas."""
    print(msg)
    log_lines.append(msg)
 
def cargar_datos(ruta_h5ad, ruta_csv):
    """Carga los datos desde el archivo h5ad y el csv de metadatos."""
    log_print(f"Cargando objeto (backed='r'): {ruta_h5ad.name}")
    adata = sc.read_h5ad(ruta_h5ad, backed="r")
    obs_h5ad = adata.obs.copy()
 
    log_print(f"Cargando csv de metadata: {ruta_csv.name}")
    df_meta = pd.read_csv(ruta_csv, index_col=0, low_memory=False)
 
    return adata, obs_h5ad, df_meta
 
def comparar_indices(obs_h5ad, df_meta):
    """Compara los indices de obs_h5ad y df_meta para verificar coincidencias."""
    n_h5ad = len(obs_h5ad)
    n_csv = len(df_meta)
    coincidencias = obs_h5ad.index.intersection(df_meta.index)
 
    log_print(f"Nucleos en h5ad: {n_h5ad}")
    log_print(f"Nucleos en csv: {n_csv}")
    log_print(f"Coincidencias exactas: {len(coincidencias)}")
 
    assert len(coincidencias) > 0, "Error: El indice del csv no coincide con el del h5ad."
    if len(coincidencias) < n_h5ad:
        log_print(f"Asvertencia: Faltan {n_h5ad - len(coincidencias)} nucleos del h5ad en el csv.")
 
def mezclar_metadatos(obs_h5ad, df_meta):
    """Integra los metadatos del CSV en obs_h5ad, priorizando la anotacion final del h5ad."""
    cols_nuevas = [c for c in df_meta.columns if c not in obs_h5ad.columns]
    obs = obs_h5ad.join(df_meta[cols_nuevas], how='left')
    return obs
 
def limpieza_y_definicion_severidad(obs):
    """Limpia y define los grupos de severidad en el DataFrame obs creado en mezclar_metadatos."""
    log_print("\nConstruyendo grupos de severidad...")
    obs['Severity_Group'] = obs['Overall AD neuropathological Change'].astype(str)
    mask_neurotypical = obs['Neurotypical reference'].isin(['Yes', 'True', True, 1, '1'])
    obs.loc[mask_neurotypical, 'Severity_Group'] = 'Not AD (Neurotypical)'
    obs['Severity_Group'] = obs['Severity_Group'].replace({'nan': 'Unknown', 'None': 'Unknown'})
 
    log_print("Conteos de Severidad:")
    log_print(obs['Severity_Group'].value_counts(dropna=False).to_string())
    return obs
 
def guardar_tabla_nivel_nucleos(obs, output_path):
    """Guarda la tabla de metadatos a nivel de nucleos en un archivo CSV."""
    obs.to_csv(output_path)
    log_print(f"-> Metadata a nivel nucleo exportada en: {output_path.name}")
 
def calidad_atac_seq(obs):
    """Analiza las metricas de calidad ATAC-seq y determina si son por nucleo o biblioteca."""
    log_print("\nAnalizando metricas de calidad ATAC-seq...")
    atac_qc_cols = [
        'ATAC_Median_high_quality_fragments_per_cell', 
        'ATAC_TSS_enrichment_score', 
        'ATAC_Fraction_of_high_quality_fragments_in_cells',
        'Doublet score'
    ]
    valid_qc_cols = [c for c in atac_qc_cols if c in obs.columns]
 
    lib_col = 'library_prep' if 'library_prep' in obs.columns else 'sample_id'
    for col in valid_qc_cols:
        max_unique_per_lib = obs.groupby(lib_col, observed=True)[col].nunique().max()
        tipo = "BIBLIOTECA" if max_unique_per_lib == 1 else "NÚCLEO"
        log_print(f" - {col}: Metrica calculada a nivel de {tipo}")
 
def exportar_metadatos_pdonador(obs, out_dir):
    """Extrae las variables clinicas a nivel de donador y verifica perfiles unicos."""
    log_print("\nGenerando tabla a nivel donador...")
    donador_cols_deseadas = [
        'Donor ID', 'Continuous Pseudo-progression Score', 'Overall AD neuropathological Change', 
        'Braak', 'Thal', 'CERAD score', 'Sex', 'Age at Death', 'APOE Genotype',
        'Severity_Group', 'Neurotypical reference'
    ]
    
    faltantes = [c for c in donador_cols_deseadas if c not in obs.columns]
    if faltantes:
        log_print(f"Aviso: Las siguientes columnas clinicas no se encontraron y se omitiran: {faltantes}")
 
    donador_cols = [c for c in donador_cols_deseadas if c in obs.columns]
    metadata_donador = obs[donador_cols].drop_duplicates()
 
    cps_counts = metadata_donador.groupby('Donor ID').size()
    if (cps_counts > 1).any():
        raise ValueError(f"Error: Donadores con multiples perfiles clinicos detectados:\n{cps_counts[cps_counts > 1]}")
 
    conteos_nuc = obs['Donor ID'].value_counts().rename("n_nuclei")
    metadata_donador = metadata_donador.merge(conteos_nuc, left_on='Donor ID', right_index=True)
    metadata_donador = metadata_donador.sort_values('Continuous Pseudo-progression Score').reset_index(drop=True)
 
    donor_out = out_dir / "metadata_donadores_extendida.csv"
    metadata_donador.to_csv(donor_out, index=False)
    log_print(f"-> Metadata de donadores exportada en: {donor_out.name}")
 
def extraccion_BED(adata, out_dir):
    """Extrae las coordenadas de los picos en formato BED validado."""
    log_print("\nProcesando BED de picos...")
    var_df = adata.var.copy().reset_index().rename(columns={'index': 'peak_id'})
 
    coord_split = var_df['peak_id'].str.extract(r'(chr[\w\d]+)[:\-](\d+)[:\-](\d+)')
    var_df['chr'] = coord_split[0]
    var_df['start'] = coord_split[1]
    var_df['end'] = coord_split[2]
 
    fallos = var_df['chr'].isna().sum()
    log_print(f"Total de variables (picos): {len(var_df)}")
    if fallos > 0:
        log_print(f"Advertencia: {fallos} picos no pudieron ser parseados y se omitiran.")
 
    bed_df = var_df.dropna(subset=['chr', 'start', 'end']).copy()
    bed_df = bed_df[bed_df['chr'].str.match(r'^chr([\d]+|[XYM])$')]
    bed_df = bed_df[['chr', 'start', 'end', 'peak_id']]
 
    bed_df['start'] = pd.to_numeric(bed_df['start'], errors='coerce')
    bed_df['end'] = pd.to_numeric(bed_df['end'], errors='coerce')
    bed_df = bed_df.dropna(subset=['start', 'end'])
    bed_df['start'] = bed_df['start'].astype(int)
    bed_df['end'] = bed_df['end'].astype(int)
 
    bed_df = bed_df.sort_values(['chr', 'start', 'end'])
 
    var_out = out_dir / "atac_peaks.bed"
    bed_df.to_csv(var_out, sep='\t', header=False, index=False)
    log_print(f"Archivo BED validado ({len(bed_df)} picos limpios) guardado en: {var_out.name}")
 
#==================================================
#Main
#==================================================
def main():
 
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    #Carga y cruce seguro de datos (h5ad vs CSV)
    adata, obs_h5ad, df_meta = cargar_datos(h5ad_path, csv_path)
 
    # Verificación de integridad en el cruce de IDs
    comparar_indices(obs_h5ad, df_meta)
 
    # Integrar datos priorizando la anotación final del h5ad
    obs = mezclar_metadatos(obs_h5ad, df_meta)
 
    # Limpieza y definición de severidad
    obs = limpieza_y_definicion_severidad(obs)
 
    # Guardar la tabla a nivel de nucleo
    guardar_tabla_nivel_nucleos(obs, nuclei_out)
 
    #Exploración de calidad de ATAC 
    calidad_atac_seq(obs)
 
    #Exportación de metadatas por donador
    exportar_metadatos_pdonador(obs,OUT_DIR)
    
    #Extracción de BED para cruce con rmsk
    extraccion_BED(adata, OUT_DIR)
 
    # Guardar logs
    with open(qc_log_path, 'w') as f:
        f.write("\n".join(log_lines))
 
 
if __name__ == "__main__":
    main()