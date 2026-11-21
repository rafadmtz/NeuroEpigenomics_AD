#!/usr/bin/env python3
"""
01b_filtered_data_graphs.py — Creation of filtered data graphs for the NeuroEpigenomics_AD project
==========================================================================

Analyzes snATAC-seq metadata from the middle temporal gyrus (MTG) of SEA-AD,
from Alzheimer's (AD) donors and controls by creating graphs.

SEA-AD Source: AWS S3 public bucket "sea-ad-single-cell-profiling"
(Gabitto et al. 2024, Nat Neurosci, "Data availability" section).

Requires having previously run src/01_metadata_analysis.py

Input:
  results/01_metadata_output/complete_nuclei_metadata.csv
      Merged metadata (h5ad + CSV all-nuclei), one row per nucleus, with the Severity_Group column already defined

Outputs (in results/01b_graphs/):
  1.- 01b_absolute_subclasses_composition.png
  2.- 01b_supertypes_composition_[subclass].png (L2/3 IT, Sst, Vip)
  3.- 01b_graphs_log.txt
 
Usage:
  The file can be executed directly from the command line, or
  imported as a module to call the main() function.


"""

#==================================================
# Libraries used
#==================================================

import matplotlib
matplotlib.use('Agg') # Headless environment for server

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from pathlib import Path

#==================================================
# Variable definitions
#==================================================


PROJECT_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = PROJECT_ROOT / "results" / "01_metadata_output"
OUT_DIR = PROJECT_ROOT / "results" / "01b_graphs"
nuclei_csv = METADATA_DIR / "complete_nuclei_metadata.csv"

POPULATIONS = ['L2/3 IT', 'Sst', 'Vip']
NUCLEI_THRESHOLD = 50

SEV_COLORS = {
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
# Functions
#==================================================

def log_print(msg):
    print(msg)
    log_lines.append(msg)

def graph_subclass_vs_severity(obs):
    """Generates a stacked bar chart of the absolute number of cells by Subclass and Severity."""
    log_print("\nGenerating absolute cell count by Subclass stacked by severity...")
    
    # Define the correct order of severity including neurotypical references
    severity_order = ['Not AD (Neurotypical)', 'Not AD', 'Low', 'Intermediate', 'High']
    
    # Create contingency table Subclass x Severity
    table = pd.crosstab(obs['Subclass'], obs['Severity_Group'])
    
    # Filter only the severity columns that exist in the real data and order them
    existing_columns = [col for col in severity_order if col in table.columns]
    table = table[existing_columns]
    
    # Order the subclasses (rows) from highest to lowest total cell count
    table = table.loc[table.sum(axis=1).sort_values(ascending=False).index]
    
    # Extract the exact colors from the global dictionary for the groups that exist
    colors_to_use = [SEV_COLORS[c] for c in existing_columns]
    
    # Generate the graph
    fig, ax = plt.subplots(figsize=(14, 7))
    table.plot(kind='bar', stacked=True, color=colors_to_use, width=0.85, ax=ax)
    
    # Format the Y axis to separate thousands with commas (e.g. 100,000)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: f'{int(y):,}'))
    
    plt.title('Number of cells by Subclass, stacked by severity group', fontsize=14)
    plt.xlabel('Subclass', fontsize=12)
    plt.ylabel('Number of cells', fontsize=12)
    plt.legend(title='Severity group', title_fontsize=11)
    plt.xticks(rotation=90)
    plt.tight_layout()
    
    out_path = OUT_DIR / "01b_absolute_subclasses_composition.png"
    fig.savefig(out_path, dpi=300)
    plt.close()
    log_print(f"-> Saved to: {out_path.name}")

def graph_supertypes_composition(obs):
    """Stacked bars of supertypes with severity marks and low N warnings."""
    log_print("\nGenerating Supertypes composition...")
    
    # Explicitly exclude donors without CPS and sort
    obs_valid = obs.dropna(subset=['Continuous Pseudo-progression Score'])
    cps_by_donor = obs_valid.groupby('Donor ID')['Continuous Pseudo-progression Score'].first().sort_values()
    sorted_donors = cps_by_donor.index
    sev_dict = obs_valid.set_index('Donor ID')['Severity_Group'].to_dict()
    
    for subc in POPULATIONS:
        df_sub = obs_valid[obs_valid['Subclass'] == subc].copy()
        if df_sub.empty: 
            log_print(f"No data for subclass {subc}. Skipping...")
            continue
            
        # Group supertypes with less than 3% in "Others"
        st_counts = df_sub['Supertype'].value_counts(normalize=True)
        rare_supertypes = st_counts[st_counts < 0.03].index
        df_sub['Supertype_Grouped'] = df_sub['Supertype'].replace(rare_supertypes, 'Others')
        
        # Contingency table Donor x Supertype (Row normalized)
        st_table = pd.crosstab(df_sub['Donor ID'], df_sub['Supertype_Grouped'], normalize='index') * 100
        sorted_idx = [d for d in sorted_donors if d in st_table.index]
        if not sorted_idx: continue
        
        plot_table = st_table.reindex(sorted_idx)
        cols = list(plot_table.columns)
        if 'Others' in cols:
            cols.remove('Others')
            cols.append('Others')
            plot_table = plot_table[cols]
            
        fig, ax = plt.subplots(figsize=(14, 7))
        plot_table.plot(kind='bar', stacked=True, colormap='tab20', width=0.85, ax=ax)
        
        # Color stripe for severity (squares below the X axis)[cite: 1]
        axis_colors = [SEV_COLORS.get(sev_dict.get(d, 'Unknown'), '#808080') for d in plot_table.index]
        ax.scatter(range(len(plot_table)), [-4] * len(plot_table), color=axis_colors, marker='s', s=80, clip_on=False)
        
        # Low N Warnings (above each bar)[cite: 1]
        # We calculate the number of nuclei for this donor and subclass on the fly
        n_dict = df_sub['Donor ID'].value_counts().to_dict()
        for i, d in enumerate(plot_table.index):
            n = n_dict.get(d, 0)
            txt_color = 'red' if n < NUCLEI_THRESHOLD else 'black'
            txt_weight = 'bold' if n < NUCLEI_THRESHOLD else 'normal'
            ax.text(i, 102, f"n={n}", ha='center', va='bottom', color=txt_color, fontsize=7, rotation=90, fontweight=txt_weight)
        
        ax.set_ylim(0, 120) # Extra vertical space for the text "n=..."
        ax.set_title(f'Supertypes Composition in {subc}', fontsize=14)
        ax.set_xlabel('Donor ID (Sorted by CPS)', fontsize=12, labelpad=15)
        ax.set_ylabel('Proportion (%)', fontsize=12)
        ax.legend(title='Supertype', bbox_to_anchor=(1.02, 1), loc='upper left')
        ax.set_xticks(range(len(plot_table)))
        ax.set_xticklabels(plot_table.index, rotation=90, fontsize=6)
        
        plt.tight_layout()
        out_path = OUT_DIR / f"01b_supertypes_composition_{subc.replace('/', '_')}.png"
        fig.savefig(out_path, dpi=300)
        plt.close()
        log_print(f"-> Saved to: {out_path.name}")

#==================================================
# Main
#==================================================

def main():
    if not nuclei_csv.exists():
        print(f"Error: {nuclei_csv} not found")
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("Loading metadata...")
    obs = pd.read_csv(nuclei_csv, index_col=0, usecols=usecols, low_memory=False)
    
    # Graph the absolute count by Subclass 
    graph_subclass_vs_severity(obs)
    
    # Graph the proportions of supertypes
    graph_supertypes_composition(obs)
    
    # Save logs
    with open(OUT_DIR / "01b_graphs_log.txt", 'w') as f:
        f.write("\n".join(log_lines))
    print("\nProcess finished")

if __name__ == "__main__":
    main()