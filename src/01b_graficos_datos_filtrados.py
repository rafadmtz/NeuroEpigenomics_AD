#!/usr/bin/env python3
"""
01b_graficos_datos_filtrados.py — Creacion de graficos de datos filtrados para el proyecto NeuroEpigenomics_AD
==========================================================================

Analiza los metadatos de snATAC-seq del giro temporal medio (MTG) de SEA-AD,
de donadores con Alzheimer (AD) y controles creando graficos.

Fuente de SEA-AD: bucket público de AWS S3 "sea-ad-single-cell-profiling"
(Gabitto et al. 2024, Nat Neurosci, seccion "Data availability").

Requiere haber corrido antes src/01_analisis_metadatos.py

Entrada:
  results/01_metadata_output/metadata_nucleos_completa.csv
      Metadatos unidos (h5ad + CSV all-nuclei), una fila por núcleo, con la columna Severity_Group ya definida

Salidas (en results/01b_graficas/):
  1.- 01b_composicion_subclases_absoluto.png
  2.- 01b_supertipos_composicion_[subclass].png (L2/3 IT, Sst, Vip)
  3.- 01b_graficos_log.txt
 
Uso:
  El archivo se puede ejecutar directamente desde la linea de comandos, o bien
  importarse como modulo y llamar a la funcion main().
 
 
"""

#==================================================
#Liberias utilizadas
#==================================================

import matplotlib
matplotlib.use('Agg') # Entorno headless para servidor

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path

#==================================================
#Definición de variables
#==================================================


PROJECT_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = PROJECT_ROOT / "results" / "01_metadata_output"
OUT_DIR = PROJECT_ROOT / "results" / "01b_graficas"
nucleos_csv = METADATA_DIR / "metadata_nucleos_completa.csv"

POBLACIONES = ['L2/3 IT', 'Sst', 'Vip']
UMBRAL_NUCLEOS = 50

COLORES_SEV = {
    'Not AD (Neurotypical)': '#067BC2',
    'Not AD': '#84BCDA',
    'Low': '#ECC30B',
    'Intermediate': '#F37748',
    'High': '#D56062',
    'Unknown': "#050505"
}

usecols = ['index', 'Donor ID', 'Severity_Group', 'Continuous Pseudo-progression Score', 'Subclass', 'Supertype']

log_lines = []

#==================================================
#Funciones
#==================================================

def log_print(msg):
    print(msg)
    log_lines.append(msg)

def graficar_subclass_vs_severidad(obs):
    """Genera un grafico de barras apiladas del numero absoluto de celulas por Subclass y Severidad."""
    log_print("\nGenerando conteo absoluto de celulas por Subclass apilado por severidad...")
    
    # Definir el orden correcto de severidad incluyendo referencias neurotípicas
    orden_severidad = ['Not AD (Neurotypical)', 'Not AD', 'Low', 'Intermediate', 'High']
    
    # Crear tabla de contingencia Subclass x Severidad
    tabla = pd.crosstab(obs['Subclass'], obs['Severity_Group'])
    
    # Filtrar solo las columnas de severidad que existan en los datos reales y ordenarlas
    columnas_existentes = [col for col in orden_severidad if col in tabla.columns]
    tabla = tabla[columnas_existentes]
    
    # Ordenar las subclases (filas) de mayor a menor cantidad total de células
    tabla = tabla.loc[tabla.sum(axis=1).sort_values(ascending=False).index]
    
    # Extraer los colores exactos del diccionario global para los grupos que existen
    colores_usar = [COLORES_SEV[c] for c in columnas_existentes]
    
    # Generar el gráfico
    fig, ax = plt.subplots(figsize=(14, 7))
    tabla.plot(kind='bar', stacked=True, color=colores_usar, width=0.85, ax=ax)
    
    # Formatear el eje Y para separar miles con comas (ej. 100,000)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: f'{int(y):,}'))
    
    plt.title('Número de células por Subclass, apilado por grupo de severidad', fontsize=14)
    plt.xlabel('Subclass', fontsize=12)
    plt.ylabel('Número de células', fontsize=12)
    plt.legend(title='Grupo de severidad', title_fontsize=11)
    plt.xticks(rotation=90)
    plt.tight_layout()
    
    out_path = OUT_DIR / "01b_composicion_subclases_absoluto.png"
    fig.savefig(out_path, dpi=300)
    plt.close()
    log_print(f"-> Guardada en: {out_path.name}")

def graficar_composicion_supertipos(obs):
    """Barras apiladas de supertipos con marcas de severidad y advertencias de baja N."""
    log_print("\nGenerando composicion de Supertipos...")
    
    # Excluir explícitamente donadores sin CPS y ordenar
    obs_valid = obs.dropna(subset=['Continuous Pseudo-progression Score'])
    cps_por_donador = obs_valid.groupby('Donor ID')['Continuous Pseudo-progression Score'].first().sort_values()
    donadores_ordenados = cps_por_donador.index
    sev_dict = obs_valid.set_index('Donor ID')['Severity_Group'].to_dict()
    
    for subc in POBLACIONES:
        df_sub = obs_valid[obs_valid['Subclass'] == subc].copy()
        if df_sub.empty: 
            log_print(f"Sin datos para la subclase {subc}. Omitiendo...")
            continue
            
        # Agrupar supertipos con menos del 3% en "Otros"
        st_conteos = df_sub['Supertype'].value_counts(normalize=True)
        supertipos_raros = st_conteos[st_conteos < 0.03].index
        df_sub['Supertype_Grouped'] = df_sub['Supertype'].replace(supertipos_raros, 'Otros')
        
        # Tabla de contingencia Donador x Supertipo (Normalizada por fila)
        tabla_st = pd.crosstab(df_sub['Donor ID'], df_sub['Supertype_Grouped'], normalize='index') * 100
        idx_ordenados = [d for d in donadores_ordenados if d in tabla_st.index]
        if not idx_ordenados: continue
        
        tabla_plot = tabla_st.reindex(idx_ordenados)
        cols = list(tabla_plot.columns)
        if 'Otros' in cols:
            cols.remove('Otros')
            cols.append('Otros')
            tabla_plot = tabla_plot[cols]
            
        fig, ax = plt.subplots(figsize=(14, 7))
        tabla_plot.plot(kind='bar', stacked=True, colormap='tab20', width=0.85, ax=ax)
        
        # Franja de color para severidad (cuadrados debajo del eje X)
        colores_eje = [COLORES_SEV.get(sev_dict.get(d, 'Unknown'), '#808080') for d in tabla_plot.index]
        ax.scatter(range(len(tabla_plot)), [-4] * len(tabla_plot), color=colores_eje, marker='s', s=80, clip_on=False)
        
        # Advertencias de Baja N (sobre cada barra)
        # Calculamos el número de núcleos de este donador y subclase al vuelo
        n_dict = df_sub['Donor ID'].value_counts().to_dict()
        for i, d in enumerate(tabla_plot.index):
            n = n_dict.get(d, 0)
            color_txt = 'red' if n < UMBRAL_NUCLEOS else 'black'
            peso_txt = 'bold' if n < UMBRAL_NUCLEOS else 'normal'
            ax.text(i, 102, f"n={n}", ha='center', va='bottom', color=color_txt, fontsize=7, rotation=90, fontweight=peso_txt)
        
        ax.set_ylim(0, 120) # Espacio extra vertical para el texto "n=..."
        ax.set_title(f'Composición de Supertipos en {subc}', fontsize=14)
        ax.set_xlabel('Donor ID (Ordenado por CPS)', fontsize=12, labelpad=15)
        ax.set_ylabel('Proporción (%)', fontsize=12)
        ax.legend(title='Supertype', bbox_to_anchor=(1.02, 1), loc='upper left')
        ax.set_xticks(range(len(tabla_plot)))
        ax.set_xticklabels(tabla_plot.index, rotation=90, fontsize=6)
        
        plt.tight_layout()
        out_path = OUT_DIR / f"01b_supertipos_composicion_{subc.replace('/', '_')}.png"
        fig.savefig(out_path, dpi=300)
        plt.close()
        log_print(f"-> Guardada en: {out_path.name}")

#==================================================
#Main
#==================================================

def main():
    if not nucleos_csv.exists():
        print(f"Error: No se encuentra {nucleos_csv}")
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("Cargando metadatos...")
    obs = pd.read_csv(nucleos_csv, index_col=0, usecols=usecols, low_memory=False)
    
    # Graficar el conteo absoluto por Subclass 
    graficar_subclass_vs_severidad(obs)
    
    # Graficar las proporciones de supertipos
    graficar_composicion_supertipos(obs)
    
    # Guardar logs
    with open(OUT_DIR / "01b_graficas_log.txt", 'w') as f:
        f.write("\n".join(log_lines))
    print("\nProceso finalizado")

if __name__ == "__main__":
    main()