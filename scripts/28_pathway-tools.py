#!/usr/bin/env python3
"""
Unified Comprehensive Microbiome Analysis Pipeline
Combines differential analysis, visualization, functional enrichment, and network analysis.
"""

import sys
import os
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import statsmodels.stats.multitest as smm
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from typing import Optional, Dict, List, Tuple

# Configuration
class Config:
    def __init__(self):
        self.PW_PATH = os.environ.get("PW_PATH", "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/picrust2_out/pathways_out/path_abun_unstrat.tsv")
        self.ASV_TABLE = os.environ.get("ASV_TABLE", "/Volumes/edKwamiBackUP/analysis_workflow/scripts/asv_table.tsv")
        self.TAXONOMY = os.environ.get("TAXONOMY", "/Volumes/edKwamiBackUP/analysis_workflow/scripts/asv_taxonomy.tsv")
        self.METADATA = os.environ.get("METADATA", "/Volumes/edKwamiBackUP/analysis_workflow/scripts/metadata.tsv")
        self.OUTDIR = os.environ.get("OUTDIR", "unified_outputs")
        self.TOPN = int(os.environ.get("TOPN", "25"))
        self.FDR_ALPHA = float(os.environ.get("FDR_ALPHA", "0.1"))

config = Config()

# Comprehensive pathway mapping
PATHWAY_NAMES = {
    'PWY-5367': 'GDP-mannose biosynthesis',
    'PWY-7874': 'superpathway of L-threonine metabolism',
    'PWY-7873': 'superpathway of L-serine and glycine biosynthesis II',
    'P101-PWY': 'xylose degradation',
    'PWY-6562': 'norspermidine biosynthesis',
    'GLYCOCAT-PWY': 'superpathway of glycol metabolism and degradation',
    'PWY-6293': 'superpathway of UDP-glucose-derived O-antigen building blocks biosynthesis',
    'PWY-5182': 'L-histidine degradation III',
    'REDCITCYC': 'TCA cycle VIII (helicobacter)',
    'PWY-801': 'superpathway of purine nucleotides salvage',
    'P461-PWY': 'hexitol fermentation to lactate, formate, ethanol and acetate',
    'PWY-6163': 'chorismate biosynthesis from 3-dehydroquinate',
    'PWY-6630': 'superpathway of L-aspartate and L-asparagine biosynthesis',
    'PWY-6901': 'superpathway of glucose and xylose degradation',
    'PWY0-1261': 'anhydromuropeptides recycling',
    'PWY-1296': 'purine ribonucleosides degradation',
    'PWY-5265': 'peptidoglycan biosynthesis I (meso-diaminopimelate containing)',
    'PWY-7117': 'C4 photosynthetic carbon assimilation cycle, NADP-ME type',
    'PWY-7222': 'guanosine deoxyribonucleotides de novo biosynthesis II',
    'PWY-7220': 'adenosine deoxyribonucleotides de novo biosynthesis II',
    'MET-SAM-PWY': 'superpathway of S-adenosyl-L-methionine biosynthesis',
    'VALSYN-PWY': 'L-valine biosynthesis',
    'COMPLETE-ARO-PWY': 'superpathway of aromatic amino acid biosynthesis',
    'GLYCOLYSIS': 'glycolysis I (from glucose 6-phosphate)',
    'TCA': 'TCA cycle I (prokaryotic)',
    'PENTOSE-P-PWY': 'pentose phosphate pathway'
}

# Functional categories
KEGG_CATEGORIES = {
    'Carbohydrate metabolism': ['glycolysis', 'gluconeogenesis', 'pentose', 'fructose', 'galactose', 'starch', 'glucose', 'mannose', 'xylose'],
    'Amino acid metabolism': ['alanine', 'arginine', 'asparagine', 'aspartate', 'cysteine', 'glutamate', 'glycine', 'histidine', 'isoleucine', 'leucine', 'lysine', 'methionine', 'phenylalanine', 'serine', 'threonine', 'tryptophan', 'tyrosine', 'valine'],
    'Nucleotide metabolism': ['purine', 'pyrimidine', 'adenine', 'guanine', 'nucleotide', 'nucleoside', 'adenosine', 'guanosine'],
    'Energy metabolism': ['TCA', 'citrate', 'oxidative phosphorylation', 'ATP synthesis'],
    'Cell wall and membrane': ['peptidoglycan', 'lipopolysaccharide', 'membrane', 'cell wall']
}

GO_CATEGORIES = {
    'Metabolic process': ['biosynthesis', 'degradation', 'metabolism', 'catabolism', 'anabolism'],
    'Transport': ['transport', 'uptake', 'efflux', 'import', 'export'],
    'Cell cycle': ['division', 'replication', 'cycle', 'growth'],
    'Response to stimulus': ['stress', 'response', 'adaptation', 'regulation']
}

# Core analysis functions
def load_metadata_groups(metadata_path: str) -> pd.Series:
    md = pd.read_csv(metadata_path, sep="\t", dtype=str)
    md = md.set_index("Sample")
    return md["Group"]

def normalize_relative(df: pd.DataFrame) -> pd.DataFrame:
    col_sums = df.sum().replace(0, 1.0)
    return (df / col_sums).fillna(0.0)

def differential_mwu(table: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    common = table.columns.intersection(groups.index)
    X = table[common]
    g = groups.loc[common]
    
    case_mask = (g == "Case")
    control_mask = (g == "Control")
    
    results = []
    for feat in X.index:
        case_vals = X.loc[feat, case_mask].values
        control_vals = X.loc[feat, control_mask].values
        
        try:
            _, p = mannwhitneyu(case_vals, control_vals, alternative="two-sided")
        except:
            p = 1.0
        
        results.append({
            "feature": feat,
            "mean_case": float(case_vals.mean()),
            "mean_control": float(control_vals.mean()),
            "diff_mean": float(case_vals.mean() - control_vals.mean()),
            "p_value": float(p)
        })
    
    df = pd.DataFrame(results).sort_values("p_value")
    df["fdr_bh"] = smm.multipletests(df["p_value"], method="fdr_bh")[1]
    df["significant"] = df["fdr_bh"] <= config.FDR_ALPHA
    return df

def wrap_text(text, width=20):
    words = text.split()
    lines = []
    current_line = []
    current_length = 0
    
    for word in words:
        if current_length + len(word) + 1 <= width:
            current_line.append(word)
            current_length += len(word) + 1
        else:
            if current_line:
                lines.append(' '.join(current_line))
            current_line = [word]
            current_length = len(word)
    
    if current_line:
        lines.append(' '.join(current_line))
    
    return '\n'.join(lines)

# Visualization functions
def plot_volcano(diff_df: pd.DataFrame, out_png: str, label_col: str) -> None:
    """Create publication-ready volcano plot."""
    df = diff_df.copy()
    df["neglog10_fdr"] = -np.log10(df["fdr_bh"].replace(0, np.finfo(float).eps))
    
    fig, ax = plt.subplots(figsize=(14, 10))
    plt.subplots_adjust(left=0.15, right=0.75, top=0.9, bottom=0.1)
    
    # Separate data
    ns = df[~df["significant"]]
    sig_case = df[df["significant"] & (df["diff_mean"] > 0)]
    sig_control = df[df["significant"] & (df["diff_mean"] < 0)]
    
    # Plot points
    ax.scatter(ns["diff_mean"], ns["neglog10_fdr"], c="lightgray", alpha=0.6, s=25, zorder=1)
    ax.scatter(sig_case["diff_mean"], sig_case["neglog10_fdr"], c="#1f77b4", alpha=0.8, s=40, zorder=2)
    ax.scatter(sig_control["diff_mean"], sig_control["neglog10_fdr"], c="#d62728", alpha=0.8, s=40, zorder=2)
    
    # Threshold line
    ax.axhline(-np.log10(config.FDR_ALPHA), color="black", linestyle="--", alpha=0.7, zorder=0)
    ax.axvline(0, color="black", linestyle="-", alpha=0.3, zorder=0)
    
    # Select top features and position labels
    top_case = sig_case.nlargest(3, "neglog10_fdr") if len(sig_case) > 0 else pd.DataFrame()
    top_control = sig_control.nlargest(2, "neglog10_fdr") if len(sig_control) > 0 else pd.DataFrame()
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_range = xlim[1] - xlim[0]
    y_range = ylim[1] - ylim[0]
    
    # Case labels
    case_x_positions = [xlim[1] - 0.5 * x_range, xlim[1] - 0.15 * x_range, xlim[1] - 0.3 * x_range]
    case_y_positions = [ylim[1] - 0.05 * y_range, ylim[1] - 0.35 * y_range, ylim[1] - 0.65 * y_range]
    
    for i, (_, row) in enumerate(top_case.iterrows()):
        if i < len(case_x_positions):
            wrapped_name = wrap_text(row[label_col], 20)
            ax.annotate(wrapped_name, xy=(row["diff_mean"], row["neglog10_fdr"]),
                       xytext=(case_x_positions[i], case_y_positions[i]), 
                       fontsize=8, ha='left', va='center',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.9),
                       arrowprops=dict(arrowstyle='->', color='blue', alpha=0.8))
    
    # Control labels
    control_x_positions = [xlim[0] + 0.2 * x_range, xlim[0] + 0.3 * x_range]
    control_y_positions = [ylim[1] - 0.1 * y_range, ylim[1] - 0.25 * y_range]
    
    for i, (_, row) in enumerate(top_control.iterrows()):
        if i < len(control_x_positions):
            wrapped_name = wrap_text(row[label_col], 20)
            ax.annotate(wrapped_name, xy=(row["diff_mean"], row["neglog10_fdr"]),
                       xytext=(control_x_positions[i], control_y_positions[i]),
                       fontsize=8, ha='right', va='center',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='lightcoral', alpha=0.9),
                       arrowprops=dict(arrowstyle='->', color='red', alpha=0.8))
    
    ax.set_xlabel("Effect Size (Case - Control)", fontsize=14, fontweight='bold')
    ax.set_ylabel("-log10(FDR-corrected p-value)", fontsize=14, fontweight='bold')
    ax.set_title("Differential Abundance Analysis", fontsize=16, fontweight='bold')
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#1f77b4', alpha=0.8, label=f'Enriched in Case (n={len(sig_case)})'),
        Patch(facecolor='#d62728', alpha=0.8, label=f'Enriched in Control (n={len(sig_control)})'),
        Patch(facecolor='lightgray', alpha=0.6, label=f'Not Significant (n={len(ns)})'),
        plt.Line2D([0], [0], color='black', linestyle='--', alpha=0.7, 
                  label=f'FDR threshold = {config.FDR_ALPHA}')
    ]
    
    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.02, 0.3), fontsize=11, frameon=True)
    
    # Directional indicators
    ax.text(xlim[1] * 0.85, ylim[1] * 0.95, 'Higher in Case', 
           ha='center', va='center', fontsize=11, style='italic',
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.6))
    ax.text(xlim[0] * 0.85, ylim[1] * 0.95, 'Higher in Control',
           ha='center', va='center', fontsize=11, style='italic', 
           bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.6))
    
    ax.grid(True, alpha=0.3)
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[Volcano] Saved: {out_png}")

def plot_heatmap(diff_df: pd.DataFrame, table_rel: pd.DataFrame, groups: pd.Series, 
                out_png: str, title: str, topn: int = 20) -> None:
    """Create heatmap of top differential features."""
    top_features = diff_df.nlargest(topn, "diff_mean").copy()
    
    common = table_rel.columns.intersection(groups.index)
    X = table_rel.loc[top_features["feature"], common]
    g = groups.loc[common]
    
    case_samples = g[g == "Case"].index
    control_samples = g[g == "Control"].index
    ordered_samples = list(case_samples) + list(control_samples)
    
    heatmap_data = X[ordered_samples]
    
    plt.figure(figsize=(12, 8))
    sns.heatmap(heatmap_data, cmap='RdYlBu_r', center=0, 
                cbar_kws={'label': 'Relative Abundance'})
    
    plt.axvline(x=len(case_samples), color='black', linewidth=2)
    plt.text(len(case_samples)/2, -1, 'Case', ha='center', fontweight='bold')
    plt.text(len(case_samples) + len(control_samples)/2, -1, 'Control', ha='center', fontweight='bold')
    
    plt.title(f'{title} - Top {topn} Features', fontsize=14, fontweight='bold')
    plt.xlabel('Samples', fontsize=12)
    plt.ylabel('Features', fontsize=12)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[Heatmap] Saved: {out_png}")

# Functional analysis
def kegg_pathway_enrichment(pathways_df: pd.DataFrame) -> pd.DataFrame:
    significant_pathways = pathways_df[pathways_df['significant'] == True]
    enrichment_results = []
    
    for category, keywords in KEGG_CATEGORIES.items():
        category_pathways = []
        for _, row in significant_pathways.iterrows():
            pathway_name = row.get('pathway_name', row['feature']).lower()
            if any(keyword in pathway_name for keyword in keywords):
                category_pathways.append(row['feature'])
        
        if len(category_pathways) > 0:
            total_significant = len(significant_pathways)
            category_count = len(category_pathways)
            enrichment_ratio = category_count / total_significant
            
            category_effects = significant_pathways[
                significant_pathways['feature'].isin(category_pathways)
            ]['diff_mean'].abs().mean()
            
            enrichment_results.append({
                'KEGG_Category': category,
                'Pathway_Count': category_count,
                'Enrichment_Ratio': enrichment_ratio,
                'Mean_Effect_Size': category_effects,
                'Pathways': '; '.join(category_pathways[:5])
            })
    
    return pd.DataFrame(enrichment_results).sort_values('Enrichment_Ratio', ascending=False)

def plot_kegg_enrichment(kegg_df: pd.DataFrame, output_path: str) -> None:
    if kegg_df.empty:
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    bars1 = ax1.barh(kegg_df['KEGG_Category'], kegg_df['Enrichment_Ratio'], 
                     color='steelblue', alpha=0.8)
    ax1.set_xlabel('Enrichment Ratio', fontweight='bold')
    ax1.set_title('KEGG Pathway Category Enrichment', fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)
    
    bars2 = ax2.barh(kegg_df['KEGG_Category'], kegg_df['Pathway_Count'], 
                     color='darkgreen', alpha=0.8)
    ax2.set_xlabel('Number of Pathways', fontweight='bold')
    ax2.set_title('KEGG Pathway Category Counts', fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[KEGG] Saved: {output_path}")

# Network analysis
def create_pathway_network(pathways_df: pd.DataFrame) -> nx.Graph:
    significant_pathways = pathways_df[pathways_df['significant'] == True]
    G = nx.Graph()
    
    for _, row in significant_pathways.iterrows():
        pathway_id = row['feature']
        pathway_name = row.get('pathway_name', pathway_id)
        effect_size = abs(row['diff_mean'])
        direction = 'Case' if row['diff_mean'] > 0 else 'Control'
        
        G.add_node(pathway_id, name=pathway_name, effect_size=effect_size,
                  direction=direction, fdr=row['fdr_bh'])
    
    # Add edges based on functional similarity
    pathways_list = list(significant_pathways['feature'])
    
    for i, pathway1 in enumerate(pathways_list):
        for j, pathway2 in enumerate(pathways_list[i+1:], i+1):
            name1 = significant_pathways[significant_pathways['feature'] == pathway1].iloc[0].get('pathway_name', pathway1).lower()
            name2 = significant_pathways[significant_pathways['feature'] == pathway2].iloc[0].get('pathway_name', pathway2).lower()
            
            words1 = set(name1.split()) - {'of', 'and', 'the', 'a', 'superpathway', 'biosynthesis'}
            words2 = set(name2.split()) - {'of', 'and', 'the', 'a', 'superpathway', 'biosynthesis'}
            
            if len(words1) > 0 and len(words2) > 0:
                similarity = len(words1.intersection(words2)) / len(words1.union(words2))
                if similarity > 0.2:
                    G.add_edge(pathway1, pathway2, weight=similarity)
    
    return G

def plot_pathway_network(G: nx.Graph, output_path: str) -> None:
    if len(G.nodes()) == 0:
        return
    
    plt.figure(figsize=(16, 12))
    pos = nx.spring_layout(G, k=3, iterations=50)
    
    node_sizes = [G.nodes[node]['effect_size'] * 5000 for node in G.nodes()]
    node_colors = ['#1f77b4' if G.nodes[node]['direction'] == 'Case' else '#d62728' 
                   for node in G.nodes()]
    edge_weights = [G[u][v]['weight'] * 5 for u, v in G.edges()]
    
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, 
                          alpha=0.8, edgecolors='black', linewidths=0.5)
    nx.draw_networkx_edges(G, pos, width=edge_weights, alpha=0.6, edge_color='gray')
    
    plt.title('Pathway Interaction Network', fontsize=14, fontweight='bold', pad=20)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[Network] Saved: {output_path}")

# Data loading and taxonomy
def load_taxonomy(tax_path: str, asv_index: pd.Index) -> Optional[pd.DataFrame]:
    if not os.path.isfile(tax_path):
        return None
    
    tax = pd.read_csv(tax_path, sep="\t", dtype=str)
    asv_set = set(asv_index.astype(str))
    best_col, best_matches = None, 0
    
    for col in tax.columns:
        matches = len(set(tax[col].astype(str)).intersection(asv_set))
        if matches > best_matches:
            best_matches = matches
            best_col = col
    
    genus_col = None
    for col in ["Genus", "genus", "taxon"]:
        if col in tax.columns:
            genus_col = col
            break
    
    if best_col and genus_col:
        return tax[[best_col, genus_col]].rename(columns={best_col: "ASV", genus_col: "Genus"})
    return None

def export_results(diff_df: pd.DataFrame, filename_base: str) -> None:
    diff_df.to_csv(f"{filename_base}.csv", index=False)
    diff_df.to_csv(f"{filename_base}.tsv", sep='\t', index=False)
    print(f"[Export] Saved: {filename_base}.[csv|tsv]")

def create_summary_stats(diff_df: pd.DataFrame, analysis_type: str) -> pd.DataFrame:
    stats = {
        'Analysis_Type': analysis_type,
        'Total_Features': len(diff_df),
        'Significant_Features': diff_df["significant"].sum(),
        'Enriched_in_Case': ((diff_df["significant"]) & (diff_df["diff_mean"] > 0)).sum(),
        'Enriched_in_Control': ((diff_df["significant"]) & (diff_df["diff_mean"] < 0)).sum(),
        'Mean_P_Value': diff_df["p_value"].mean(),
        'Mean_FDR': diff_df["fdr_bh"].mean(),
        'Max_Effect_Size': diff_df["diff_mean"].abs().max(),
        'Min_P_Value': diff_df["p_value"].min(),
        'Min_FDR': diff_df["fdr_bh"].min()
    }
    return pd.DataFrame([stats])

def main():
    """Unified comprehensive analysis pipeline."""
    print("[Info] Starting unified comprehensive analysis...")
    os.makedirs(config.OUTDIR, exist_ok=True)
    
    # Load metadata
    groups = load_metadata_groups(config.METADATA)
    
    # === PATHWAY ANALYSIS ===
    print("[Info] Pathway analysis...")
    pw_data = pd.read_csv(config.PW_PATH, sep="\t", index_col=0)
    pw_data = pw_data.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    
    # Check orientation
    col_overlap = len(set(pw_data.columns).intersection(set(groups.index)))
    idx_overlap = len(set(pw_data.index).intersection(set(groups.index)))
    if idx_overlap > col_overlap:
        pw_data = pw_data.T
    
    pw_rel = normalize_relative(pw_data)
    pw_diff = differential_mwu(pw_rel, groups)
    pw_diff["pathway_name"] = pw_diff["feature"].map(PATHWAY_NAMES).fillna(pw_diff["feature"])
    
    # Pathway visualizations
    plot_volcano(pw_diff, os.path.join(config.OUTDIR, "pathway_volcano.png"), "pathway_name")
    plot_heatmap(pw_diff, pw_rel, groups, os.path.join(config.OUTDIR, "pathway_heatmap.png"), "Pathway Analysis")
    
    # Functional enrichment
    kegg_enrichment = kegg_pathway_enrichment(pw_diff)
    kegg_enrichment.to_csv(os.path.join(config.OUTDIR, "kegg_enrichment.csv"), index=False)
    plot_kegg_enrichment(kegg_enrichment, os.path.join(config.OUTDIR, "kegg_enrichment.png"))
    
    # Network analysis
    pathway_network = create_pathway_network(pw_diff)
    plot_pathway_network(pathway_network, os.path.join(config.OUTDIR, "pathway_network.png"))
    
    # Export pathway results
    export_results(pw_diff, os.path.join(config.OUTDIR, "pathway_results"))
    pw_stats = create_summary_stats(pw_diff, "Pathways")
    
    # === GENUS ANALYSIS ===
    print("[Info] Genus analysis...")
    genus_stats = pd.DataFrame()
    
    if os.path.isfile(config.ASV_TABLE):
        asv_data = pd.read_csv(config.ASV_TABLE, sep="\t", index_col=0)
        asv_data = asv_data.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        
        col_overlap = len(set(asv_data.columns).intersection(set(groups.index)))
        idx_overlap = len(set(asv_data.index).intersection(set(groups.index)))
        if idx_overlap > col_overlap:
            asv_data = asv_data.T
        
        taxonomy = load_taxonomy(config.TAXONOMY, asv_data.index)
        
        if taxonomy is not None:
            merged = taxonomy.merge(asv_data, left_on="ASV", right_index=True, how="inner")
            genus_counts = merged.drop(columns=["ASV"]).groupby("Genus").sum()
            genus_rel = normalize_relative(genus_counts)
            genus_diff = differential_mwu(genus_rel, groups)
            genus_diff["genus"] = genus_diff["feature"]
            
            plot_volcano(genus_diff, os.path.join(config.OUTDIR, "genus_volcano.png"), "genus")
            plot_heatmap(genus_diff, genus_rel, groups, os.path.join(config.OUTDIR, "genus_heatmap.png"), "Genus Analysis")
            
            export_results(genus_diff, os.path.join(config.OUTDIR, "genus_results"))
            genus_stats = create_summary_stats(genus_diff, "Genera")
    
    # === SUMMARY ===
    if not genus_stats.empty:
        combined_stats = pd.concat([pw_stats, genus_stats], ignore_index=True)
    else:
        combined_stats = pw_stats
    
    combined_stats.to_csv(os.path.join(config.OUTDIR, "summary_statistics.csv"), index=False)
    
    print(f"\n[Complete] Unified analysis saved to {config.OUTDIR}")
    print("Generated files:")
    print("  - pathway_volcano.png - Publication-ready volcano plot")
    print("  - pathway_heatmap.png - Clustered heatmap")
    print("  - kegg_enrichment.png - Functional category enrichment")
    print("  - pathway_network.png - Pathway interaction network")
    print("  - genus_volcano.png - Genus differential abundance")
    print("  - genus_heatmap.png - Genus abundance heatmap")
    print("  - pathway_results.[csv|tsv] - Complete pathway statistics")
    print("  - genus_results.[csv|tsv] - Complete genus statistics")
    print("  - kegg_enrichment.csv - KEGG enrichment results")
    print("  - summary_statistics.csv - Analysis overview")

if __name__ == "__main__":
    main()
