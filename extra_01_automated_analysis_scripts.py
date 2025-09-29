# Standard library imports
import os
import subprocess
import tempfile
import json

# Third-party imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import cross_val_score, StratifiedKFold, cross_val_predict
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (confusion_matrix, classification_report, roc_curve, 
                           auc, precision_recall_curve, average_precision_score)
from sklearn.feature_selection import SelectFromModel
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import networkx as nx
from sklearn.manifold import MDS, TSNE
from umap import UMAP
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova
from statsmodels.stats.multitest import multipletests

# Pathway-Tools integration
try:
    import ptools
    from ptools import PathwayTools
    PATHWAY_TOOLS_AVAILABLE = True
    print("Pathway-Tools Python API detected and loaded successfully.")
except ImportError:
    PATHWAY_TOOLS_AVAILABLE = False
    print("Pathway-Tools Python API not found. Using command-line interface.")

def check_pathway_tools_installation():
    """Check if Pathway-Tools is properly installed and accessible"""
    try:
        # Try to run pathway-tools command
        result = subprocess.run(['pathway-tools', '--version'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            print(f"Pathway-Tools found: {result.stdout.strip()}")
            return True
        else:
            print("Pathway-Tools command failed")
            return False
    except (subprocess.TimeoutExpired, FileNotFoundError):
        print("Pathway-Tools command not found in PATH")
        return False

def create_metabolic_network_dataframe(df_case, df_ctrl, output_dir):
    """Create metabolic network data from taxonomic abundance data"""
    metabolic_dir = os.path.join(output_dir, 'metabolic_analysis')
    os.makedirs(metabolic_dir, exist_ok=True)
    
    # Define common metabolic pathways and their associated taxa
    # This is a simplified mapping - in practice, you'd use a comprehensive database
    metabolic_pathways = {
        'Glycolysis': ['Streptococcus', 'Lactobacillus', 'Bifidobacterium'],
        'TCA_Cycle': ['Pseudomonas', 'Escherichia', 'Salmonella'],
        'Amino_Acid_Metabolism': ['Bacteroides', 'Clostridium', 'Escherichia'],
        'Fatty_Acid_Metabolism': ['Bifidobacterium', 'Lactobacillus', 'Bacteroides'],
        'Vitamin_Synthesis': ['Bifidobacterium', 'Lactobacillus', 'Escherichia'],
        'Antibiotic_Resistance': ['Pseudomonas', 'Staphylococcus', 'Enterococcus'],
        'Short_Chain_Fatty_Acid_Production': ['Bacteroides', 'Clostridium', 'Bifidobacterium'],
        'Bile_Acid_Metabolism': ['Bacteroides', 'Clostridium', 'Lactobacillus'],
        'Neurotransmitter_Metabolism': ['Bifidobacterium', 'Lactobacillus', 'Escherichia'],
        'Inflammation_Modulation': ['Bacteroides', 'Bifidobacterium', 'Lactobacillus']
    }
    
    # Calculate pathway abundances based on taxonomic composition
    pathway_data = {}
    
    for pathway, associated_taxa in metabolic_pathways.items():
        case_pathway_abundance = []
        ctrl_pathway_abundance = []
        
        # Calculate pathway abundance for each sample
        for col in df_case.columns:
            sample_abundance = 0
            for taxon in associated_taxa:
                matching_taxa = df_case.index[df_case.index.str.contains(taxon, case=False, na=False)]
                if len(matching_taxa) > 0:
                    sample_abundance += df_case.loc[matching_taxa, col].sum()
            case_pathway_abundance.append(sample_abundance)
        
        for col in df_ctrl.columns:
            sample_abundance = 0
            for taxon in associated_taxa:
                matching_taxa = df_ctrl.index[df_ctrl.index.str.contains(taxon, case=False, na=False)]
                if len(matching_taxa) > 0:
                    sample_abundance += df_ctrl.loc[matching_taxa, col].sum()
            ctrl_pathway_abundance.append(sample_abundance)
        
        pathway_data[pathway] = {
            'case_abundance': case_pathway_abundance,
            'control_abundance': ctrl_pathway_abundance
        }
    
    # Create DataFrame for pathway analysis
    # First, determine the maximum length needed
    max_case_length = max(len(data['case_abundance']) for data in pathway_data.values())
    max_ctrl_length = max(len(data['control_abundance']) for data in pathway_data.values())
    max_length = max(max_case_length, max_ctrl_length)
    
    # Create DataFrame with proper index
    pathway_df = pd.DataFrame(index=range(max_length))
    
    for pathway, data in pathway_data.items():
        case_abundance = data['case_abundance']
        ctrl_abundance = data['control_abundance']
        
        # Pad case abundance if needed
        if len(case_abundance) < max_length:
            case_abundance = case_abundance + [0] * (max_length - len(case_abundance))
        
        # Pad control abundance if needed
        if len(ctrl_abundance) < max_length:
            ctrl_abundance = ctrl_abundance + [0] * (max_length - len(ctrl_abundance))
        
        pathway_df[f'{pathway}_Case'] = case_abundance
        pathway_df[f'{pathway}_Control'] = ctrl_abundance
    
    # Save pathway data
    pathway_df.to_csv(os.path.join(metabolic_dir, 'metabolic_pathway_abundances.csv'))
    
    return pathway_df, metabolic_pathways

def analyze_metabolic_pathways(pathway_df, output_dir):
    """Analyze metabolic pathway differences between groups"""
    metabolic_dir = os.path.join(output_dir, 'metabolic_analysis')
    os.makedirs(metabolic_dir, exist_ok=True)
    
    results = []
    
    # Analyze each pathway
    for col in pathway_df.columns:
        if col.endswith('_Case'):
            pathway_name = col.replace('_Case', '')
            ctrl_col = f'{pathway_name}_Control'
            
            if ctrl_col in pathway_df.columns:
                case_values = pathway_df[col]
                ctrl_values = pathway_df[ctrl_col]
                
                # Statistical test
                stat, pval = stats.mannwhitneyu(case_values, ctrl_values, alternative='two-sided')
                
                # Calculate effect size
                case_mean = case_values.mean()
                ctrl_mean = ctrl_values.mean()
                log2fc = np.log2((case_mean + 0.1) / (ctrl_mean + 0.1))
                
                results.append({
                    'Pathway': pathway_name,
                    'Case_mean': case_mean,
                    'Control_mean': ctrl_mean,
                    'Log2FoldChange': log2fc,
                    'P_value': pval,
                    'Effect_size': abs(log2fc)
                })
                
                # Create pathway-specific plots
                plt.figure(figsize=(8, 6))
                plt.boxplot([case_values, ctrl_values], labels=['Case', 'Control'])
                plt.title(f'{pathway_name} Pathway Abundance\n(p = {pval:.3f})')
                plt.ylabel('Relative Abundance')
                plt.tight_layout()
                plt.savefig(os.path.join(metabolic_dir, f'{pathway_name.lower()}_pathway.png'))
                plt.close()
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    results_df['FDR_P_value'] = fdrcorrection(results_df['P_value'])[1]
    results_df = results_df.sort_values('FDR_P_value')
    
    # Create volcano plot for pathways
    plt.figure(figsize=(10, 8))
    plt.scatter(results_df['Log2FoldChange'], 
               -np.log10(results_df['FDR_P_value']),
               alpha=0.7, s=100)
    
    # Add labels for significant pathways
    significant = results_df['FDR_P_value'] < 0.05
    for idx, row in results_df[significant].iterrows():
        plt.annotate(row['Pathway'], 
                    (row['Log2FoldChange'], -np.log10(row['FDR_P_value'])),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=10, fontweight='bold')
    
    plt.axhline(-np.log10(0.05), color='red', linestyle='--', label='FDR = 0.05')
    plt.axvline(0, color='gray', linestyle='-', alpha=0.5)
    plt.xlabel('Log2 Fold Change (Case/Control)')
    plt.ylabel('-log10(FDR adjusted p-value)')
    plt.title('Metabolic Pathway Volcano Plot')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(metabolic_dir, 'pathway_volcano_plot.png'))
    plt.close()
    
    # Create pathway heatmap
    pathway_heatmap_data = []
    pathway_names = []
    
    for _, row in results_df.iterrows():
        pathway_heatmap_data.append([row['Case_mean'], row['Control_mean']])
        pathway_names.append(row['Pathway'])
    
    if pathway_heatmap_data:
        plt.figure(figsize=(8, len(pathway_names) * 0.4))
        sns.heatmap(pathway_heatmap_data, 
                   xticklabels=['Case', 'Control'],
                   yticklabels=pathway_names,
                   cmap='YlOrRd', annot=True, fmt='.3f')
        plt.title('Metabolic Pathway Abundance Heatmap')
        plt.tight_layout()
        plt.savefig(os.path.join(metabolic_dir, 'pathway_heatmap.png'))
        plt.close()
    
    # Save results
    results_df.to_csv(os.path.join(metabolic_dir, 'metabolic_pathway_analysis.csv'), index=False)
    
    return results_df

def create_pathway_network_analysis(pathway_df, output_dir):
    """Create metabolic network analysis using Pathway-Tools concepts"""
    metabolic_dir = os.path.join(output_dir, 'metabolic_analysis')
    os.makedirs(metabolic_dir, exist_ok=True)
    
    # Calculate correlations between pathways
    pathway_corr = pathway_df.corr(method='spearman')
    
    # Create network
    G = nx.Graph()
    
    # Add edges for significant correlations
    correlation_threshold = 0.6
    for i in range(len(pathway_corr.columns)):
        for j in range(i+1, len(pathway_corr.columns)):
            if abs(pathway_corr.iloc[i,j]) > correlation_threshold:
                G.add_edge(pathway_corr.columns[i], 
                          pathway_corr.columns[j],
                          weight=abs(pathway_corr.iloc[i,j]))
    
    # Calculate network metrics
    network_metrics = {
        'Nodes': len(G.nodes()),
        'Edges': len(G.edges()),
        'Average_Degree': sum(dict(G.degree()).values()) / len(G.nodes()) if len(G.nodes()) > 0 else 0,
        'Density': nx.density(G),
        'Average_Clustering': nx.average_clustering(G) if len(G.nodes()) > 0 else 0,
        'Number_of_Components': nx.number_connected_components(G)
    }
    
    # Plot metabolic network
    if len(G.nodes()) > 0:
        plt.figure(figsize=(12, 10))
        pos = nx.spring_layout(G, k=1, iterations=50)
        
        # Draw nodes
        nx.draw_networkx_nodes(G, pos, node_color='lightgreen',
                             node_size=2000, alpha=0.7)
        
        # Draw edges with varying widths
        edge_widths = [G[u][v]['weight'] * 3 for u, v in G.edges()]
        nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.5, edge_color='gray')
        
        # Add labels
        nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold')
        
        plt.title(f'Metabolic Pathway Network\n(correlations > {correlation_threshold})')
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(os.path.join(metabolic_dir, 'metabolic_network.png'))
        plt.close()
    
    # Save network metrics
    pd.DataFrame([network_metrics]).to_csv(
        os.path.join(metabolic_dir, 'metabolic_network_metrics.csv'), index=False)
    
    return network_metrics

def generate_pathway_enrichment_analysis(df_case, df_ctrl, output_dir):
    """Generate pathway enrichment analysis using Pathway-Tools concepts"""
    metabolic_dir = os.path.join(output_dir, 'metabolic_analysis')
    os.makedirs(metabolic_dir, exist_ok=True)
    
    # Define metabolic subsystems and their associated functions
    metabolic_subsystems = {
        'Energy_Metabolism': {
            'taxa': ['Streptococcus', 'Lactobacillus', 'Bifidobacterium', 'Escherichia'],
            'functions': ['glycolysis', 'fermentation', 'respiration']
        },
        'Amino_Acid_Metabolism': {
            'taxa': ['Bacteroides', 'Clostridium', 'Escherichia', 'Pseudomonas'],
            'functions': ['biosynthesis', 'degradation', 'transport']
        },
        'Lipid_Metabolism': {
            'taxa': ['Bifidobacterium', 'Lactobacillus', 'Bacteroides'],
            'functions': ['fatty_acid_synthesis', 'bile_acid_metabolism']
        },
        'Vitamin_Metabolism': {
            'taxa': ['Bifidobacterium', 'Lactobacillus', 'Escherichia'],
            'functions': ['vitamin_b_synthesis', 'vitamin_k_synthesis']
        },
        'Neurotransmitter_Metabolism': {
            'taxa': ['Bifidobacterium', 'Lactobacillus', 'Escherichia'],
            'functions': ['serotonin_synthesis', 'gaba_metabolism']
        }
    }
    
    enrichment_results = []
    
    for subsystem, info in metabolic_subsystems.items():
        # Calculate subsystem abundance
        case_subsystem_abundance = []
        ctrl_subsystem_abundance = []
        
        for col in df_case.columns:
            abundance = 0
            for taxon in info['taxa']:
                matching_taxa = df_case.index[df_case.index.str.contains(taxon, case=False, na=False)]
                if len(matching_taxa) > 0:
                    abundance += df_case.loc[matching_taxa, col].sum()
            case_subsystem_abundance.append(abundance)
        
        for col in df_ctrl.columns:
            abundance = 0
            for taxon in info['taxa']:
                matching_taxa = df_ctrl.index[df_ctrl.index.str.contains(taxon, case=False, na=False)]
                if len(matching_taxa) > 0:
                    abundance += df_ctrl.loc[matching_taxa, col].sum()
            ctrl_subsystem_abundance.append(abundance)
        
        # Statistical test
        stat, pval = stats.mannwhitneyu(case_subsystem_abundance, ctrl_subsystem_abundance)
        
        # Calculate enrichment score
        case_mean = np.mean(case_subsystem_abundance)
        ctrl_mean = np.mean(ctrl_subsystem_abundance)
        enrichment_score = np.log2((case_mean + 0.1) / (ctrl_mean + 0.1))
        
        enrichment_results.append({
            'Subsystem': subsystem,
            'Case_mean': case_mean,
            'Control_mean': ctrl_mean,
            'Enrichment_score': enrichment_score,
            'P_value': pval,
            'Functions': ', '.join(info['functions'])
        })
    
    # Create enrichment results DataFrame
    enrichment_df = pd.DataFrame(enrichment_results)
    enrichment_df['FDR_P_value'] = fdrcorrection(enrichment_df['P_value'])[1]
    enrichment_df = enrichment_df.sort_values('FDR_P_value')
    
    # Create enrichment plot
    plt.figure(figsize=(12, 8))
    colors = ['red' if row['FDR_P_value'] < 0.05 else 'blue' for _, row in enrichment_df.iterrows()]
    plt.barh(enrichment_df['Subsystem'], enrichment_df['Enrichment_score'], color=colors, alpha=0.7)
    plt.axvline(0, color='black', linestyle='-', alpha=0.5)
    plt.xlabel('Log2 Enrichment Score (Case/Control)')
    plt.title('Metabolic Subsystem Enrichment Analysis')
    plt.tight_layout()
    plt.savefig(os.path.join(metabolic_dir, 'subsystem_enrichment.png'))
    plt.close()
    
    # Save results
    enrichment_df.to_csv(os.path.join(metabolic_dir, 'subsystem_enrichment_analysis.csv'), index=False)
    
    return enrichment_df

def create_pathway_tools_script(taxa_list, output_dir):
    """Create a Pathway-Tools script for metabolic analysis"""
    metabolic_dir = os.path.join(output_dir, 'metabolic_analysis')
    os.makedirs(metabolic_dir, exist_ok=True)
    
    script_content = f"""#!/usr/bin/env python
# Pathway-Tools Analysis Script
# Generated for microbiome analysis

import os
import sys

# Add Pathway-Tools to path if needed
# sys.path.append('/path/to/pathway-tools')

def analyze_metabolic_capabilities(taxa_list):
    \"\"\"
    Analyze metabolic capabilities of the given taxa using Pathway-Tools
    \"\"\"
    
    # This is a template for Pathway-Tools analysis
    # You would need to customize this based on your specific Pathway-Tools setup
    
    print("Analyzing metabolic capabilities for taxa:")
    for taxon in taxa_list:
        print(f"  - {{taxon}}")
    
    # Example Pathway-Tools commands (these would need to be adapted)
    # ptools = PathwayTools()
    # for taxon in taxa_list:
    #     pathways = ptools.get_pathways_for_organism(taxon)
    #     reactions = ptools.get_reactions_for_organism(taxon)
    #     compounds = ptools.get_compounds_for_organism(taxon)
    
    print("\\nMetabolic analysis complete!")

if __name__ == "__main__":
    taxa_list = {taxa_list}
    analyze_metabolic_capabilities(taxa_list)
"""
    
    script_path = os.path.join(metabolic_dir, 'pathway_tools_analysis.py')
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    # Make script executable
    os.chmod(script_path, 0o755)
    
    return script_path

def load_data(case_file, control_file):
    """Load and preprocess case and control data"""
    df_case = pd.read_csv(case_file)
    df_ctrl = pd.read_csv(control_file)
    
    # Set taxonomic names as index if they're in a column
    if 'Taxon' in df_case.columns:
        df_case = df_case.set_index('Taxon')
    if 'Taxon' in df_ctrl.columns:
        df_ctrl = df_ctrl.set_index('Taxon')
    
    # Convert to numeric and handle missing values
    for col in df_case.columns:
        df_case[col] = pd.to_numeric(df_case[col], errors='coerce')
    for col in df_ctrl.columns:
        df_ctrl[col] = pd.to_numeric(df_ctrl[col], errors='coerce')
    
    # Fill NaN with 0 and convert to relative abundance
    df_case = df_case.fillna(0)
    df_ctrl = df_ctrl.fillna(0)
    df_case = df_case.div(df_case.sum()) * 100
    df_ctrl = df_ctrl.div(df_ctrl.sum()) * 100
    
    return df_case, df_ctrl

def analyze_beta_diversity(df_case, df_ctrl, output_dir):
    """Perform comprehensive beta diversity analysis using multiple distance metrics"""
    beta_dir = os.path.join(output_dir, 'beta_diversity')
    os.makedirs(beta_dir, exist_ok=True)
    
    # Preprocess data for numerical stability
    scaler = StandardScaler()
    combined_data = pd.concat([df_case.T, df_ctrl.T])
    scaled_data = scaler.fit_transform(combined_data)
    scaled_data = np.clip(scaled_data, -10, 10)  # Clip extreme values
    scaled_df = pd.DataFrame(scaled_data, index=combined_data.index, columns=combined_data.columns)
    
    groups = ['Case'] * len(df_case.T) + ['Control'] * len(df_ctrl.T)
    
    # Set consistent color scheme
    group_colors = {
        'Case': '#E41A1C',    # Strong red
        'Control': '#377EB8'   # Strong blue
    }
    
    # Set color maps for different plot types
    distance_cmap = 'YlOrRd'  # Yellow-Orange-Red for distance matrices
    cluster_cmap = 'viridis'  # Viridis for clustering
    
    # Calculate multiple distance metrics
    metrics = {
        'braycurtis': 'Bray-Curtis',
        'jaccard': 'Jaccard',
        'euclidean': 'Euclidean',
        'correlation': 'Correlation',
        'cosine': 'Cosine',
        'canberra': 'Canberra',
        'chebyshev': 'Chebyshev',
        'cityblock': 'Manhattan'
    }
    
    results = []
    
    for metric_name, metric_label in metrics.items():
        try:
            # Calculate distance matrix
            dist_matrix = pdist(scaled_df, metric=metric_name)
            dist_matrix_square = squareform(dist_matrix)
            
            # Calculate within and between group distances
            case_indices = np.where(np.array(groups) == 'Case')[0]
            ctrl_indices = np.where(np.array(groups) == 'Control')[0]
            
            within_case = dist_matrix_square[np.ix_(case_indices, case_indices)]
            within_ctrl = dist_matrix_square[np.ix_(ctrl_indices, ctrl_indices)]
            between = dist_matrix_square[np.ix_(case_indices, ctrl_indices)]
            
            # Calculate summary statistics
            results.append({
                'Metric': metric_label,
                'Mean_Within_Case': np.mean(within_case[np.triu_indices_from(within_case, k=1)]),
                'SD_Within_Case': np.std(within_case[np.triu_indices_from(within_case, k=1)]),
                'Mean_Within_Control': np.mean(within_ctrl[np.triu_indices_from(within_ctrl, k=1)]),
                'SD_Within_Control': np.std(within_ctrl[np.triu_indices_from(within_ctrl, k=1)]),
                'Mean_Between': np.mean(between),
                'SD_Between': np.std(between)
            })
            
            # Plot distance matrix heatmap
            plt.figure(figsize=(10, 8))
            sns.heatmap(dist_matrix_square, cmap=distance_cmap)
            plt.title(f'{metric_label} Distance Matrix')
            plt.tight_layout()
            plt.savefig(os.path.join(beta_dir, f'distance_matrix_{metric_name}.png'))
            plt.close()
            
            # Perform hierarchical clustering
            linkage_matrix = linkage(dist_matrix, method='ward')
            plt.figure(figsize=(12, 8))
            dendrogram(linkage_matrix, labels=groups, leaf_rotation=90,
                      leaf_font_size=8, link_color_func=lambda k: 'black')
            plt.title(f'Hierarchical Clustering ({metric_label} Distance)')
            plt.tight_layout()
            plt.savefig(os.path.join(beta_dir, f'dendrogram_{metric_name}.png'))
            plt.close()
            
            # Perform PCoA
            pcoa_results = pcoa(dist_matrix_square)
            
            # Plot PCoA results
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')
            
            for group in ['Case', 'Control']:
                mask = np.array(groups) == group
                ax.scatter(pcoa_results.samples.values[mask, 0], 
                          pcoa_results.samples.values[mask, 1], 
                          pcoa_results.samples.values[mask, 2],
                          c=[group_colors[group]], label=group, alpha=0.7)
            
            ax.set_xlabel(f'PC1 ({pcoa_results.proportion_explained[0]:.1%})')
            ax.set_ylabel(f'PC2 ({pcoa_results.proportion_explained[1]:.1%})')
            ax.set_zlabel(f'PC3 ({pcoa_results.proportion_explained[2]:.1%})')
            ax.set_title(f'PCoA Plot ({metric_label} Distance)')
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(beta_dir, f'pcoa_{metric_name}.png'))
            plt.close()
            
            # Perform PERMANOVA
            f_stat, p_val = permanova(dist_matrix_square, groups, permutations=999)
            with open(os.path.join(beta_dir, f'permanova_{metric_name}.txt'), 'w') as f:
                f.write(f"PERMANOVA Results ({metric_label} Distance):\n")
                f.write(f"F-statistic: {f_stat:.3f}\n")
                f.write(f"p-value: {p_val:.3e}\n")
            
        except Exception as e:
            print(f"Warning: {metric_label} distance produced invalid values, skipping...")
            continue
    
    return pd.DataFrame(results)

def analyze_core_microbiome(df_case, df_ctrl, output_dir, prevalence_threshold=0.5, abundance_threshold=0.1):
    """Analyze core microbiome in each group"""
    core_dir = os.path.join(output_dir, 'core_microbiome')
    os.makedirs(core_dir, exist_ok=True)
    
    def get_core_taxa(df, prev_thresh, abund_thresh):
        # Calculate prevalence
        prevalence = (df > abund_thresh).mean(axis=1)
        # Calculate mean abundance
        mean_abundance = df.mean(axis=1)
        # Get core taxa
        core = df.index[prevalence >= prev_thresh]
        return core, prevalence, mean_abundance
    
    # Get core taxa for each group
    core_case, prev_case, abund_case = get_core_taxa(df_case, prevalence_threshold, abundance_threshold)
    core_ctrl, prev_ctrl, abund_ctrl = get_core_taxa(df_ctrl, prevalence_threshold, abundance_threshold)
    
    # Create Venn diagram data
    core_sets = {
        'Case': set(core_case),
        'Control': set(core_ctrl)
    }
    
    # Save core taxa lists
    pd.Series(list(core_sets['Case'])).to_csv(
        os.path.join(core_dir, 'core_taxa_case.csv'), index=False, header=['Taxon'])
    pd.Series(list(core_sets['Control'])).to_csv(
        os.path.join(core_dir, 'core_taxa_control.csv'), index=False, header=['Taxon'])
    
    # Create prevalence vs abundance plots
    for name, (prevalence, abundance) in [('Case', (prev_case, abund_case)), 
                                        ('Control', (prev_ctrl, abund_ctrl))]:
        plt.figure(figsize=(10, 8))
        plt.scatter(prevalence, abundance, alpha=0.6)
        plt.axvline(prevalence_threshold, color='red', linestyle='--', 
                   label=f'Prevalence threshold ({prevalence_threshold})')
        plt.axhline(abundance_threshold, color='blue', linestyle='--',
                   label=f'Abundance threshold ({abundance_threshold})')
        
        # Add labels for core taxa
        core_taxa = core_case if name == 'Case' else core_ctrl
        for taxon in core_taxa:
            plt.annotate(taxon, 
                        (prevalence[taxon], abundance[taxon]),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=8)
        
        plt.xlabel('Prevalence')
        plt.ylabel('Mean Relative Abundance')
        plt.title(f'{name} Group: Prevalence vs Abundance')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(core_dir, f'{name.lower()}_prevalence_abundance.png'))
        plt.close()
    
    return core_sets

def analyze_rare_taxa(df_case, df_ctrl, output_dir, rarity_threshold=0.1):
    """Analyze rare taxa in each group"""
    rare_dir = os.path.join(output_dir, 'rare_taxa')
    os.makedirs(rare_dir, exist_ok=True)
    
    def get_rare_taxa(df, threshold):
        # Calculate mean abundance
        mean_abundance = df.mean(axis=1)
        # Get rare taxa
        rare = df.index[mean_abundance < threshold]
        return rare, mean_abundance
    
    # Get rare taxa for each group
    rare_case, abund_case = get_rare_taxa(df_case, rarity_threshold)
    rare_ctrl, abund_ctrl = get_rare_taxa(df_ctrl, rarity_threshold)
    
    # Create summary statistics
    rare_stats = pd.DataFrame({
        'Group': ['Case', 'Control'],
        'Total_Taxa': [len(df_case.index), len(df_ctrl.index)],
        'Rare_Taxa': [len(rare_case), len(rare_ctrl)],
        'Rare_Percentage': [len(rare_case)/len(df_case.index)*100, 
                          len(rare_ctrl)/len(df_ctrl.index)*100]
    })
    
    rare_stats.to_csv(os.path.join(rare_dir, 'rare_taxa_statistics.csv'), index=False)
    
    # Plot abundance distributions
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Case group
    sns.histplot(data=abund_case, ax=ax1, bins=30)
    ax1.axvline(rarity_threshold, color='red', linestyle='--', 
                label=f'Rarity threshold ({rarity_threshold})')
    ax1.set_title('Case Group: Abundance Distribution')
    ax1.set_xlabel('Mean Relative Abundance')
    ax1.set_ylabel('Count')
    ax1.legend()
    
    # Control group
    sns.histplot(data=abund_ctrl, ax=ax2, bins=30)
    ax2.axvline(rarity_threshold, color='red', linestyle='--',
                label=f'Rarity threshold ({rarity_threshold})')
    ax2.set_title('Control Group: Abundance Distribution')
    ax2.set_xlabel('Mean Relative Abundance')
    ax2.set_ylabel('Count')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(rare_dir, 'abundance_distributions.png'))
    plt.close()
    
    return rare_stats

def analyze_taxonomic_ratios(df_case, df_ctrl):
    """Calculate taxonomic ratios between groups"""
    ratios = {}
    
    # Add pseudocount to avoid division by zero
    pseudocount = 1e-5
    df_case = df_case + pseudocount
    df_ctrl = df_ctrl + pseudocount
    
    # Calculate Firmicutes to Bacteroidetes ratio
    firmicutes_case = df_case.filter(like='Firmicutes').sum()
    bacteroidetes_case = df_case.filter(like='Bacteroidetes').sum()
    fb_ratio_case = np.log2(firmicutes_case / bacteroidetes_case)  # Log transform for better numerical stability
    
    firmicutes_ctrl = df_ctrl.filter(like='Firmicutes').sum()
    bacteroidetes_ctrl = df_ctrl.filter(like='Bacteroidetes').sum()
    fb_ratio_ctrl = np.log2(firmicutes_ctrl / bacteroidetes_ctrl)
    
    ratios['Firmicutes_Bacteroidetes'] = {
        'case_values': fb_ratio_case.values,
        'control_values': fb_ratio_ctrl.values
    }
    
    # Calculate Proteobacteria to Bacteroidetes ratio
    proteobacteria_case = df_case.filter(like='Proteobacteria').sum()
    pb_ratio_case = np.log2(proteobacteria_case / bacteroidetes_case)
    
    proteobacteria_ctrl = df_ctrl.filter(like='Proteobacteria').sum()
    pb_ratio_ctrl = np.log2(proteobacteria_ctrl / bacteroidetes_ctrl)
    
    ratios['Proteobacteria_Bacteroidetes'] = {
        'case_values': pb_ratio_case.values,
        'control_values': pb_ratio_ctrl.values
    }
    
    return ratios

def analyze_alpha_diversity(df_case, df_ctrl, output_dir):
    """Perform comprehensive alpha diversity analysis"""
    alpha_dir = os.path.join(output_dir, 'alpha_diversity')
    os.makedirs(alpha_dir, exist_ok=True)
    
    def calculate_shannon(x):
        p = x[x > 0] / x[x > 0].sum()
        return -(p * np.log(p)).sum() if len(p) > 0 else 0
    
    def calculate_simpson(x):
        p = x[x > 0] / x[x > 0].sum()
        return 1 - (p ** 2).sum() if len(p) > 0 else 0
    
    def calculate_evenness(x):
        shannon = calculate_shannon(x)
        return shannon / np.log(len(x[x > 0])) if len(x[x > 0]) > 0 else 0
    
    # Calculate diversity metrics
    metrics = {
        'Shannon': (calculate_shannon, 'Shannon Diversity'),
        'Simpson': (calculate_simpson, 'Simpson Diversity'),
        'Evenness': (calculate_evenness, 'Pielou Evenness'),
        'Richness': (lambda x: sum(x > 0), 'Species Richness')
    }
    
    results = pd.DataFrame()
    for metric, (func, title) in metrics.items():
        case_values = df_case.apply(func)
        ctrl_values = df_ctrl.apply(func)
        
        # Statistical test
        stat, pval = stats.mannwhitneyu(case_values, ctrl_values)
        
        # Plot distributions
        plt.figure(figsize=(8, 6))
        plt.boxplot([case_values, ctrl_values], tick_labels=['Case', 'Control'])
        plt.title(f'{title}\n(p = {pval:.3f})')
        plt.ylabel(metric)
        plt.tight_layout()
        plt.savefig(os.path.join(alpha_dir, f'{metric.lower()}_diversity.png'))
        plt.close()
        
        # Store results
        results[f'{metric}_Case'] = case_values
        results[f'{metric}_Control'] = ctrl_values
        results[f'{metric}_pvalue'] = pval
    
    # Save results
    results.to_csv(os.path.join(alpha_dir, 'alpha_diversity_metrics.csv'))
    return results

def analyze_machine_learning(df_case, df_ctrl, output_dir):
    """Perform enhanced machine learning analysis using Random Forest"""
    ml_dir = os.path.join(output_dir, 'machine_learning')
    os.makedirs(ml_dir, exist_ok=True)
    
    # Prepare data
    X = pd.concat([df_case.T, df_ctrl.T])
    y = np.array(['Case'] * len(df_case.T) + ['Control'] * len(df_ctrl.T))
    
    # Scale features and handle NaN values
    scaler = StandardScaler()
    X_filled = X.fillna(X.mean())
    X_scaled = scaler.fit_transform(X_filled)
    X_scaled = pd.DataFrame(X_scaled, index=X.index, columns=X.columns)
    
    # Initialize model and cross-validation
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    # Cross-validation predictions and probabilities
    y_pred = cross_val_predict(rf, X_scaled, y, cv=cv)
    y_prob = cross_val_predict(rf, X_scaled, y, cv=cv, method='predict_proba')
    
    # Calculate ROC curve and AUC
    fpr, tpr, _ = roc_curve(y == 'Case', y_prob[:, 1])
    roc_auc = auc(fpr, tpr)
    
    # Plot ROC curve
    plt.figure(figsize=(8, 8))
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(os.path.join(ml_dir, 'roc_curve.png'))
    plt.close()
    
    # Calculate and plot Precision-Recall curve
    precision, recall, _ = precision_recall_curve(y == 'Case', y_prob[:, 1])
    avg_precision = average_precision_score(y == 'Case', y_prob[:, 1])
    
    plt.figure(figsize=(8, 8))
    plt.plot(recall, precision, color='blue', lw=2, 
             label=f'Precision-Recall curve (AP = {avg_precision:.2f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.legend(loc="lower left")
    plt.tight_layout()
    plt.savefig(os.path.join(ml_dir, 'precision_recall_curve.png'))
    plt.close()
    
    # Plot confusion matrix
    cm = confusion_matrix(y, y_pred)
    plt.figure(figsize=(8, 8))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
                xticklabels=['Case', 'Control'],
                yticklabels=['Case', 'Control'])
    plt.title('Confusion Matrix')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.tight_layout()
    plt.savefig(os.path.join(ml_dir, 'confusion_matrix.png'))
    plt.close()
    
    # Feature importance analysis
    rf.fit(X_scaled, y)
    feature_importance = pd.DataFrame({
        'feature': X.columns,
        'importance': rf.feature_importances_
    }).sort_values('importance', ascending=False)
    
    # Plot feature importance (top 20)
    plt.figure(figsize=(12, 6))
    sns.barplot(data=feature_importance.head(20), x='importance', y='feature')
    plt.title('Top 20 Most Important Features')
    plt.xlabel('Feature Importance')
    plt.tight_layout()
    plt.savefig(os.path.join(ml_dir, 'feature_importance.png'))
    plt.close()
    
    # Feature importance distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(data=feature_importance, x='importance', bins=30)
    plt.title('Distribution of Feature Importance Scores')
    plt.xlabel('Feature Importance')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(os.path.join(ml_dir, 'feature_importance_distribution.png'))
    plt.close()
    
    # Dimensionality reduction and clustering visualization
    # t-SNE with perplexity adjustment for small datasets
    n_samples = len(X_scaled)
    perplexity = min(30, n_samples - 1)  # Adjust perplexity based on sample size
    
    tsne = TSNE(n_components=2, random_state=42, perplexity=perplexity)
    X_tsne = tsne.fit_transform(X_scaled)
    
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=[0 if label == 'Case' else 1 for label in y],
                         cmap='coolwarm', alpha=0.6)
    plt.colorbar(scatter)
    plt.title('t-SNE Visualization of Samples')
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.tight_layout()
    plt.savefig(os.path.join(ml_dir, 'tsne_visualization.png'))
    plt.close()
    
    # UMAP
    umap = UMAP(random_state=42, n_neighbors=min(15, n_samples - 1))
    X_umap = umap.fit_transform(X_scaled)
    
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(X_umap[:, 0], X_umap[:, 1], c=[0 if label == 'Case' else 1 for label in y],
                         cmap='coolwarm', alpha=0.6)
    plt.colorbar(scatter)
    plt.title('UMAP Visualization of Samples')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.tight_layout()
    plt.savefig(os.path.join(ml_dir, 'umap_visualization.png'))
    plt.close()
    
    # Cross-validation performance distribution
    cv_scores = cross_val_score(rf, X_scaled, y, cv=cv, scoring='accuracy')
    
    plt.figure(figsize=(8, 6))
    sns.boxplot(y=cv_scores)
    plt.title('Cross-validation Performance Distribution')
    plt.ylabel('Accuracy Score')
    plt.tight_layout()
    plt.savefig(os.path.join(ml_dir, 'cv_performance_distribution.png'))
    plt.close()
    
    # Save detailed classification report
    report = classification_report(y, y_pred)
    with open(os.path.join(ml_dir, 'classification_report.txt'), 'w') as f:
        f.write('Classification Report\n')
        f.write('=====================\n\n')
        f.write(report)
    
    # Save feature importance rankings
    feature_importance.to_csv(os.path.join(ml_dir, 'feature_importance.csv'), index=False)
    
    return {
        'cv_scores': cv_scores,
        'mean_cv_score': cv_scores.mean(),
        'std_cv_score': cv_scores.std(),
        'feature_importance': list(zip(feature_importance['feature'], feature_importance['importance'])),
        'roc_auc': roc_auc,
        'avg_precision': avg_precision,
        'confusion_matrix': cm
    }

def analyze_differential_abundance(df_case, df_ctrl, output_dir):
    """Perform differential abundance analysis"""
    da_dir = os.path.join(output_dir, 'differential_abundance')
    os.makedirs(da_dir, exist_ok=True)
    
    results = []
    for taxon in set(df_case.index) | set(df_ctrl.index):
        case_values = df_case.loc[taxon] if taxon in df_case.index else pd.Series([0] * len(df_case.columns))
        ctrl_values = df_ctrl.loc[taxon] if taxon in df_ctrl.index else pd.Series([0] * len(df_ctrl.columns))
        
        # Calculate statistics
        stat, pval = stats.mannwhitneyu(case_values, ctrl_values, alternative='two-sided')
        
        # Calculate effect size (log2 fold change)
        case_mean = case_values.mean()
        ctrl_mean = ctrl_values.mean()
        log2fc = np.log2((case_mean + 0.1) / (ctrl_mean + 0.1))
        
        results.append({
            'Taxon': taxon,
            'Case_mean': case_mean,
            'Control_mean': ctrl_mean,
            'Log2FoldChange': log2fc,
            'P_value': pval
        })
    
    # Create DataFrame and adjust p-values
    results_df = pd.DataFrame(results)
    results_df['FDR_P_value'] = fdrcorrection(results_df['P_value'])[1]
    results_df = results_df.sort_values('FDR_P_value')
    
    # Create volcano plot
    plt.figure(figsize=(10, 8))
    plt.scatter(results_df['Log2FoldChange'], 
               -np.log10(results_df['FDR_P_value']),
               alpha=0.5)
    
    # Add labels for significant taxa
    significant = results_df['FDR_P_value'] < 0.05
    for idx, row in results_df[significant].iterrows():
        plt.annotate(row['Taxon'], 
                    (row['Log2FoldChange'], -np.log10(row['FDR_P_value'])),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=8)
    
    plt.axhline(-np.log10(0.05), color='red', linestyle='--', label='FDR = 0.05')
    plt.xlabel('Log2 Fold Change (Case/Control)')
    plt.ylabel('-log10(FDR adjusted p-value)')
    plt.title('Volcano Plot of Differential Abundance')
    plt.tight_layout()
    plt.savefig(os.path.join(da_dir, 'volcano_plot.png'))
    plt.close()
    
    # Save results
    results_df.to_csv(os.path.join(da_dir, 'differential_abundance.csv'), index=False)
    return results_df

def analyze_taxonomic_profiles(df_case, df_ctrl, output_dir):
    """Analyze taxonomic profiles and generate visualizations"""
    tax_dir = os.path.join(output_dir, 'taxonomic_profiles')
    os.makedirs(tax_dir, exist_ok=True)
    
    # Calculate mean abundances across all samples
    combined_data = pd.concat([df_case, df_ctrl], axis=1)
    mean_abundances = combined_data.mean(axis=1).sort_values(ascending=False)
    
    # Get top 20 taxa
    top_taxa = list(zip(mean_abundances.index[:20], mean_abundances.values[:20]))
    
    # Plot top 20 taxa abundances
    plt.figure(figsize=(12, 6))
    taxa, abundances = zip(*top_taxa[:20])
    plt.bar(range(len(taxa)), abundances)
    plt.xticks(range(len(taxa)), taxa, rotation=45, ha='right')
    plt.title('Top 20 Most Abundant Taxa')
    plt.ylabel('Mean Relative Abundance')
    plt.tight_layout()
    plt.savefig(os.path.join(tax_dir, 'top_taxa_abundance.png'))
    plt.close()
    
    # Calculate group-specific abundances
    case_means = df_case.mean(axis=1)
    ctrl_means = df_ctrl.mean(axis=1)
    
    # Create heatmap of top taxa
    plt.figure(figsize=(10, 8))
    top_taxa_df = pd.DataFrame({
        'Case': case_means[mean_abundances.index[:20]],
        'Control': ctrl_means[mean_abundances.index[:20]]
    })
    sns.heatmap(top_taxa_df, cmap='YlOrRd', annot=True, fmt='.3f')
    plt.title('Top 20 Taxa Abundance by Group')
    plt.tight_layout()
    plt.savefig(os.path.join(tax_dir, 'top_taxa_heatmap.png'))
    plt.close()
    
    return {
        'top_taxa': top_taxa,
        'case_means': case_means,
        'ctrl_means': ctrl_means
    }

def perform_network_analysis(df_case, df_ctrl, output_dir, correlation_threshold=0.6):
    """Perform network analysis based on correlations"""
    network_dir = os.path.join(output_dir, 'network_analysis')
    os.makedirs(network_dir, exist_ok=True)
    
    def create_network(df, group_name):
        # Calculate correlations
        corr = df.T.corr(method='spearman')
        
        # Create network
        G = nx.Graph()
        
        # Add edges for significant correlations
        for i in range(len(corr.columns)):
            for j in range(i+1, len(corr.columns)):
                if abs(corr.iloc[i,j]) > correlation_threshold:
                    G.add_edge(corr.columns[i], 
                             corr.columns[j],
                             weight=abs(corr.iloc[i,j]))
        
        # Calculate network metrics
        metrics = {
            'Nodes': len(G.nodes()),
            'Edges': len(G.edges()),
            'Average_Degree': sum(dict(G.degree()).values()) / len(G.nodes()) if len(G.nodes()) > 0 else 0,
            'Density': nx.density(G),
            'Average_Clustering': nx.average_clustering(G) if len(G.nodes()) > 0 else 0
        }
        
        # Plot network
        if len(G.nodes()) > 0:
            plt.figure(figsize=(12, 12))
            pos = nx.spring_layout(G)
            
            # Draw nodes
            nx.draw_networkx_nodes(G, pos, node_color='lightblue',
                                 node_size=1000, alpha=0.6)
            
            # Draw edges with varying widths
            edge_widths = [G[u][v]['weight'] * 2 for u, v in G.edges()]
            nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.4)
            
            # Add labels
            nx.draw_networkx_labels(G, pos, font_size=8)
            
            plt.title(f'{group_name} Correlation Network\n(correlations > {correlation_threshold})')
            plt.axis('off')
            plt.tight_layout()
            plt.savefig(os.path.join(network_dir, f'{group_name.lower()}_network.png'))
            plt.close()
        
        return G, metrics
    
    # Create networks for both groups
    G_case, metrics_case = create_network(df_case, 'Case')
    G_ctrl, metrics_ctrl = create_network(df_ctrl, 'Control')
    
    # Save network metrics
    metrics_df = pd.DataFrame({
        'Metric': metrics_case.keys(),
        'Case': metrics_case.values(),
        'Control': metrics_ctrl.values()
    })
    metrics_df.to_csv(os.path.join(network_dir, 'network_metrics.csv'), index=False)
    
    return metrics_df

def generate_analysis_report(alpha_results, summary_data, ml_results, da_results, 
                        tax_results, network_results, core_sets, rare_stats, ratios, 
                        pathway_results, metabolic_network_metrics, enrichment_results, output_dir):
    """Generate a comprehensive analysis report in text format"""
    report_path = os.path.join(output_dir, 'complete_analysis_results.txt')
    
    with open(report_path, 'w') as f:
        f.write("COMPREHENSIVE MICROBIOME ANALYSIS REPORT\n")
        f.write("=====================================\n\n")
        
        # Alpha Diversity
        f.write("1. ALPHA DIVERSITY ANALYSIS\n")
        f.write("-------------------------\n")
        for metric in ['Shannon', 'Simpson', 'Evenness', 'Richness']:
            case_values = alpha_results[f'{metric}_Case']
            ctrl_values = alpha_results[f'{metric}_Control']
            f.write(f"\n{metric} Index:\n")
            f.write(f"  Case:    Mean = {np.mean(case_values):.3f} ± {np.std(case_values):.3f}\n")
            f.write(f"  Control: Mean = {np.mean(ctrl_values):.3f} ± {np.std(ctrl_values):.3f}\n")
            stat, pval = stats.mannwhitneyu(case_values, ctrl_values)
            f.write(f"  Mann-Whitney U test p-value: {pval:.3e}\n")
        
        # Beta Diversity
        f.write("\n2. BETA DIVERSITY ANALYSIS\n")
        f.write("------------------------\n")
        summary_df = pd.DataFrame(summary_data)
        for _, row in summary_df.iterrows():
            f.write(f"\n{row['Metric']} Distance:\n")
            f.write(f"  Within Case:    {row['Mean_Within_Case']:.3f} ± {row['SD_Within_Case']:.3f}\n")
            f.write(f"  Within Control: {row['Mean_Within_Control']:.3f} ± {row['SD_Within_Control']:.3f}\n")
            f.write(f"  Between Groups: {row['Mean_Between']:.3f} ± {row['SD_Between']:.3f}\n")
        
        # Machine Learning
        f.write("\n3. MACHINE LEARNING ANALYSIS\n")
        f.write("--------------------------\n")
        f.write(f"Model Performance:\n")
        f.write(f"  Cross-validation Score: {ml_results['mean_cv_score']:.3f} ± {ml_results['std_cv_score']:.3f}\n")
        f.write(f"  ROC AUC Score: {ml_results['roc_auc']:.3f}\n")
        f.write(f"  Average Precision Score: {ml_results['avg_precision']:.3f}\n")
        
        f.write("\nConfusion Matrix:\n")
        cm = ml_results['confusion_matrix']
        f.write(f"  True Case, Predicted Case: {cm[0,0]}\n")
        f.write(f"  True Case, Predicted Control: {cm[0,1]}\n")
        f.write(f"  True Control, Predicted Case: {cm[1,0]}\n")
        f.write(f"  True Control, Predicted Control: {cm[1,1]}\n")
        
        f.write("\nTop 10 Important Features:\n")
        for feature, importance in ml_results['feature_importance'][:10]:
            f.write(f"  {feature}: {importance:.3f}\n")
        
        # Differential Abundance
        f.write("\n4. DIFFERENTIAL ABUNDANCE ANALYSIS\n")
        f.write("--------------------------------\n")
        sig_taxa = da_results[da_results['FDR_P_value'] < 0.05]
        f.write(f"Number of significantly different taxa (FDR < 0.05): {len(sig_taxa)}\n")
        if len(sig_taxa) > 0:
            f.write("\nTop significant taxa:\n")
            for _, row in sig_taxa.head().iterrows():
                f.write(f"  {row['Taxon']}: log2FC = {row['Log2FoldChange']:.2f}, FDR = {row['FDR_P_value']:.2e}\n")
        
        # Taxonomic Profiles
        f.write("\n5. TAXONOMIC PROFILE ANALYSIS\n")
        f.write("----------------------------\n")
        f.write("Top 10 abundant taxa:\n")
        for taxon, abundance in tax_results['top_taxa'][:10]:
            f.write(f"  {taxon}: {abundance:.3%}\n")
        
        # Network Analysis
        f.write("\n6. NETWORK ANALYSIS\n")
        f.write("------------------\n")
        for metric in network_results['Metric'].unique():
            case_val = network_results[network_results['Metric'] == metric]['Case'].values[0]
            ctrl_val = network_results[network_results['Metric'] == metric]['Control'].values[0]
            f.write(f"{metric}:\n")
            f.write(f"  Case:    {case_val:.3f}\n")
            f.write(f"  Control: {ctrl_val:.3f}\n")
        
        # Core Microbiome
        f.write("\n7. CORE MICROBIOME ANALYSIS\n")
        f.write("-------------------------\n")
        f.write(f"Number of core taxa in Case group: {len(core_sets['Case'])}\n")
        f.write(f"Number of core taxa in Control group: {len(core_sets['Control'])}\n")
        f.write(f"Number of shared core taxa: {len(core_sets['Case'] & core_sets['Control'])}\n")
        
        # Rare Taxa
        f.write("\n8. RARE TAXA ANALYSIS\n")
        f.write("--------------------\n")
        f.write(rare_stats.to_string())
        
        # Taxonomic Ratios
        f.write("\n\n9. TAXONOMIC RATIO ANALYSIS\n")
        f.write("-------------------------\n")
        for ratio_name, ratio_data in ratios.items():
            f.write(f"\n{ratio_name}:\n")
            f.write(f"  Case:    Mean = {np.mean(ratio_data['case_values']):.3f} ± {np.std(ratio_data['case_values']):.3f}\n")
            f.write(f"  Control: Mean = {np.mean(ratio_data['control_values']):.3f} ± {np.std(ratio_data['control_values']):.3f}\n")
            stat, pval = stats.mannwhitneyu(ratio_data['case_values'], ratio_data['control_values'])
            f.write(f"  Mann-Whitney U test p-value: {pval:.3e}\n")
        
        # Metabolic Pathway Analysis
        f.write("\n\n10. METABOLIC PATHWAY ANALYSIS\n")
        f.write("-----------------------------\n")
        sig_pathways = pathway_results[pathway_results['FDR_P_value'] < 0.05]
        f.write(f"Number of significantly different pathways (FDR < 0.05): {len(sig_pathways)}\n")
        if len(sig_pathways) > 0:
            f.write("\nTop significant pathways:\n")
            for _, row in sig_pathways.head().iterrows():
                f.write(f"  {row['Pathway']}: log2FC = {row['Log2FoldChange']:.2f}, FDR = {row['FDR_P_value']:.2e}\n")
        
        f.write("\nAll pathway analysis results:\n")
        for _, row in pathway_results.iterrows():
            f.write(f"  {row['Pathway']}: log2FC = {row['Log2FoldChange']:.2f}, p = {row['P_value']:.3e}, FDR = {row['FDR_P_value']:.3e}\n")
        
        # Metabolic Network Analysis
        f.write("\n\n11. METABOLIC NETWORK ANALYSIS\n")
        f.write("-----------------------------\n")
        f.write("Metabolic pathway network metrics:\n")
        for metric, value in metabolic_network_metrics.items():
            f.write(f"  {metric}: {value:.3f}\n")
        
        # Pathway Enrichment Analysis
        f.write("\n\n12. PATHWAY ENRICHMENT ANALYSIS\n")
        f.write("------------------------------\n")
        sig_subsystems = enrichment_results[enrichment_results['FDR_P_value'] < 0.05]
        f.write(f"Number of significantly enriched subsystems (FDR < 0.05): {len(sig_subsystems)}\n")
        if len(sig_subsystems) > 0:
            f.write("\nSignificantly enriched subsystems:\n")
            for _, row in sig_subsystems.iterrows():
                f.write(f"  {row['Subsystem']}: Enrichment = {row['Enrichment_score']:.2f}, FDR = {row['FDR_P_value']:.2e}\n")
                f.write(f"    Functions: {row['Functions']}\n")
        
        f.write("\nAll subsystem enrichment results:\n")
        for _, row in enrichment_results.iterrows():
            f.write(f"  {row['Subsystem']}: Enrichment = {row['Enrichment_score']:.2f}, p = {row['P_value']:.3e}, FDR = {row['FDR_P_value']:.3e}\n")
        
        f.write("\nAnalysis complete! Results saved in visualization folder.\n")

def advanced_pathway_tools_analysis(taxa_list, output_dir):
    """Perform advanced Pathway-Tools analysis when available"""
    metabolic_dir = os.path.join(output_dir, 'metabolic_analysis')
    os.makedirs(metabolic_dir, exist_ok=True)
    
    if not PATHWAY_TOOLS_AVAILABLE:
        print("Pathway-Tools Python API not available. Skipping advanced analysis.")
        return None
    
    print("Starting advanced Pathway-Tools analysis...")
    
    try:
        # Initialize Pathway-Tools
        pt = PathwayTools()
        pt.start_pathway_tools()
        
        # Dictionary to store results
        advanced_results = {
            'pathways': {},
            'reactions': {},
            'compounds': {},
            'enzymes': {},
            'metabolic_networks': {}
        }
        
        # Analyze each taxon
        for taxon in taxa_list:
            print(f"Analyzing {taxon}...")
            
            try:
                # Get pathways for this organism
                pathways = pt.get_pathways_for_organism(taxon)
                advanced_results['pathways'][taxon] = pathways
                
                # Get reactions for this organism
                reactions = pt.get_reactions_for_organism(taxon)
                advanced_results['reactions'][taxon] = reactions
                
                # Get compounds for this organism
                compounds = pt.get_compounds_for_organism(taxon)
                advanced_results['compounds'][taxon] = compounds
                
                # Get enzymes for this organism
                enzymes = pt.get_enzymes_for_organism(taxon)
                advanced_results['enzymes'][taxon] = enzymes
                
            except Exception as e:
                print(f"Warning: Could not analyze {taxon}: {e}")
                continue
        
        # Save advanced results
        with open(os.path.join(metabolic_dir, 'advanced_pathway_tools_results.json'), 'w') as f:
            json.dump(advanced_results, f, indent=2, default=str)
        
        # Generate summary report
        generate_advanced_ptools_report(advanced_results, metabolic_dir)
        
        # Close Pathway-Tools
        pt.stop_pathway_tools()
        
        return advanced_results
        
    except Exception as e:
        print(f"Error in advanced Pathway-Tools analysis: {e}")
        return None

def generate_advanced_ptools_report(advanced_results, output_dir):
    """Generate a detailed report from advanced Pathway-Tools analysis"""
    report_path = os.path.join(output_dir, 'advanced_pathway_tools_report.txt')
    
    with open(report_path, 'w') as f:
        f.write("ADVANCED PATHWAY-TOOLS ANALYSIS REPORT\n")
        f.write("====================================\n\n")
        
        # Summary statistics
        f.write("SUMMARY STATISTICS\n")
        f.write("-----------------\n")
        f.write(f"Total taxa analyzed: {len(advanced_results['pathways'])}\n")
        
        # Pathway analysis
        f.write("\nPATHWAY ANALYSIS\n")
        f.write("----------------\n")
        for taxon, pathways in advanced_results['pathways'].items():
            f.write(f"\n{taxon}:\n")
            f.write(f"  Number of pathways: {len(pathways) if pathways else 0}\n")
            if pathways:
                f.write("  Top pathways:\n")
                for pathway in pathways[:5]:  # Show top 5
                    f.write(f"    - {pathway}\n")
        
        # Reaction analysis
        f.write("\nREACTION ANALYSIS\n")
        f.write("----------------\n")
        for taxon, reactions in advanced_results['reactions'].items():
            f.write(f"\n{taxon}:\n")
            f.write(f"  Number of reactions: {len(reactions) if reactions else 0}\n")
        
        # Compound analysis
        f.write("\nCOMPOUND ANALYSIS\n")
        f.write("-----------------\n")
        for taxon, compounds in advanced_results['compounds'].items():
            f.write(f"\n{taxon}:\n")
            f.write(f"  Number of compounds: {len(compounds) if compounds else 0}\n")
        
        # Enzyme analysis
        f.write("\nENZYME ANALYSIS\n")
        f.write("---------------\n")
        for taxon, enzymes in advanced_results['enzymes'].items():
            f.write(f"\n{taxon}:\n")
            f.write(f"  Number of enzymes: {len(enzymes) if enzymes else 0}\n")
        
        f.write("\nAdvanced Pathway-Tools analysis complete!\n")

def create_metabolic_network_from_ptools(advanced_results, output_dir):
    """Create metabolic network visualization from Pathway-Tools data"""
    if not advanced_results:
        return None
    
    metabolic_dir = os.path.join(output_dir, 'metabolic_analysis')
    
    # Create network from reactions
    G = nx.Graph()
    
    for taxon, reactions in advanced_results['reactions'].items():
        if reactions:
            for reaction in reactions:
                # Add reaction as node
                G.add_node(f"{taxon}_{reaction}", type='reaction', taxon=taxon)
                
                # Add connections between reactions (simplified)
                # In practice, you'd parse reaction equations to find shared compounds
    
    # Save network
    if len(G.nodes()) > 0:
        plt.figure(figsize=(15, 12))
        pos = nx.spring_layout(G, k=2, iterations=50)
        
        # Color nodes by type
        node_colors = ['lightblue' if G.nodes[node]['type'] == 'reaction' else 'lightgreen' 
                      for node in G.nodes()]
        
        nx.draw(G, pos, node_color=node_colors, node_size=500, 
               with_labels=False, alpha=0.7)
        
        plt.title('Metabolic Network from Pathway-Tools')
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(os.path.join(metabolic_dir, 'advanced_metabolic_network.png'))
        plt.close()
        
        # Save network data
        nx.write_gexf(G, os.path.join(metabolic_dir, 'advanced_metabolic_network.gexf'))
    
    return G

def analyze_pathway_enrichment_with_ptools(advanced_results, df_case, df_ctrl, output_dir):
    """Analyze pathway enrichment using actual Pathway-Tools data"""
    if not advanced_results:
        return None
    
    metabolic_dir = os.path.join(output_dir, 'metabolic_analysis')
    
    # Create pathway abundance matrix
    all_pathways = set()
    for pathways in advanced_results['pathways'].values():
        if pathways:
            all_pathways.update(pathways)
    
    pathway_abundance = pd.DataFrame(index=list(all_pathways), 
                                   columns=df_case.columns.tolist() + df_ctrl.columns.tolist())
    
    # Calculate pathway abundances based on taxonomic composition
    for pathway in all_pathways:
        for col in df_case.columns:
            abundance = 0
            for taxon in df_case.index:
                if taxon in advanced_results['pathways'] and pathway in advanced_results['pathways'][taxon]:
                    abundance += df_case.loc[taxon, col]
            pathway_abundance.loc[pathway, col] = abundance
        
        for col in df_ctrl.columns:
            abundance = 0
            for taxon in df_ctrl.index:
                if taxon in advanced_results['pathways'] and pathway in advanced_results['pathways'][taxon]:
                    abundance += df_ctrl.loc[taxon, col]
            pathway_abundance.loc[pathway, col] = abundance
    
    # Perform differential abundance analysis
    results = []
    for pathway in all_pathways:
        case_values = pathway_abundance.loc[pathway, df_case.columns]
        ctrl_values = pathway_abundance.loc[pathway, df_ctrl.columns]
        
        if case_values.sum() > 0 or ctrl_values.sum() > 0:
            stat, pval = stats.mannwhitneyu(case_values, ctrl_values, alternative='two-sided')
            
            case_mean = case_values.mean()
            ctrl_mean = ctrl_values.mean()
            log2fc = np.log2((case_mean + 0.1) / (ctrl_mean + 0.1))
            
            results.append({
                'Pathway': pathway,
                'Case_mean': case_mean,
                'Control_mean': ctrl_mean,
                'Log2FoldChange': log2fc,
                'P_value': pval
            })
    
    if results:
        results_df = pd.DataFrame(results)
        results_df['FDR_P_value'] = fdrcorrection(results_df['P_value'])[1]
        results_df = results_df.sort_values('FDR_P_value')
        
        # Save results
        results_df.to_csv(os.path.join(metabolic_dir, 'advanced_pathway_enrichment.csv'), index=False)
        
        # Create volcano plot
        plt.figure(figsize=(12, 8))
        plt.scatter(results_df['Log2FoldChange'], 
                   -np.log10(results_df['FDR_P_value']),
                   alpha=0.7, s=100)
        
        # Add labels for significant pathways
        significant = results_df['FDR_P_value'] < 0.05
        for idx, row in results_df[significant].iterrows():
            plt.annotate(row['Pathway'], 
                        (row['Log2FoldChange'], -np.log10(row['FDR_P_value'])),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=8, fontweight='bold')
        
        plt.axhline(-np.log10(0.05), color='red', linestyle='--', label='FDR = 0.05')
        plt.axvline(0, color='gray', linestyle='-', alpha=0.5)
        plt.xlabel('Log2 Fold Change (Case/Control)')
        plt.ylabel('-log10(FDR adjusted p-value)')
        plt.title('Advanced Pathway Enrichment Analysis')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(metabolic_dir, 'advanced_pathway_volcano.png'))
        plt.close()
        
        return results_df
    
    return None

def main():
    # Define paths
    case_file = '/Volumes/edKwamiBackUP/Manuscripts/In_Prepartion/AD.Oral/oral/species/case.csv'
    control_file = '/Volumes/edKwamiBackUP/Manuscripts/In_Prepartion/AD.Oral/oral/species/control.csv'
    output_dir = '/Volumes/edKwamiBackUP/Manuscripts/In_Prepartion/AD.Oral/bacteria_visualizations'
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Check Pathway-Tools installation
    print("Checking Pathway-Tools installation...")
    ptools_available = check_pathway_tools_installation()
    
    print("Loading data...")
    df_case, df_ctrl = load_data(case_file, control_file)
    
    print("Analyzing alpha diversity...")
    alpha_results = analyze_alpha_diversity(df_case, df_ctrl, output_dir)
    
    print("Analyzing beta diversity...")
    summary_data = analyze_beta_diversity(df_case, df_ctrl, output_dir)
    
    print("Performing machine learning analysis...")
    ml_results = analyze_machine_learning(df_case, df_ctrl, output_dir)
    
    print("Analyzing differential abundance...")
    da_results = analyze_differential_abundance(df_case, df_ctrl, output_dir)
    
    print("Analyzing taxonomic profiles...")
    tax_results = analyze_taxonomic_profiles(df_case, df_ctrl, output_dir)
    
    print("Performing network analysis...")
    network_results = perform_network_analysis(df_case, df_ctrl, output_dir)
    
    print("Analyzing core microbiome...")
    core_sets = analyze_core_microbiome(df_case, df_ctrl, output_dir)
    
    print("Analyzing rare taxa...")
    rare_stats = analyze_rare_taxa(df_case, df_ctrl, output_dir)
    
    print("Analyzing taxonomic ratios...")
    ratios = analyze_taxonomic_ratios(df_case, df_ctrl)
    
    # Pathway-Tools based metabolic analysis
    print("Performing metabolic pathway analysis...")
    pathway_df, metabolic_pathways = create_metabolic_network_dataframe(df_case, df_ctrl, output_dir)
    
    print("Analyzing metabolic pathway differences...")
    pathway_results = analyze_metabolic_pathways(pathway_df, output_dir)
    
    print("Creating metabolic network analysis...")
    metabolic_network_metrics = create_pathway_network_analysis(pathway_df, output_dir)
    
    print("Generating pathway enrichment analysis...")
    enrichment_results = generate_pathway_enrichment_analysis(df_case, df_ctrl, output_dir)
    
    # Create Pathway-Tools script for further analysis
    all_taxa = list(set(df_case.index) | set(df_ctrl.index))
    print("Creating Pathway-Tools analysis script...")
    ptools_script = create_pathway_tools_script(all_taxa, output_dir)
    
    # Advanced Pathway-Tools analysis (if available)
    print("Performing advanced Pathway-Tools analysis...")
    advanced_results = advanced_pathway_tools_analysis(all_taxa, output_dir)
    
    if advanced_results:
        print("Creating advanced metabolic network...")
        advanced_network = create_metabolic_network_from_ptools(advanced_results, output_dir)
        
        print("Analyzing advanced pathway enrichment...")
        advanced_enrichment = analyze_pathway_enrichment_with_ptools(advanced_results, df_case, df_ctrl, output_dir)
    
    print("Generating comprehensive analysis report...")
    generate_analysis_report(alpha_results, summary_data, ml_results, da_results, 
                           tax_results, network_results, core_sets, rare_stats, ratios, 
                           pathway_results, metabolic_network_metrics, enrichment_results, output_dir)
    
    print("\nAnalysis complete! Results saved in:", output_dir)
    
    # Print key findings
    print("\nKey findings:")
    print(f"1. Number of core taxa in Case group: {len(core_sets['Case'])}")
    print(f"2. Number of core taxa in Control group: {len(core_sets['Control'])}")
    print(f"3. Number of shared core taxa: {len(core_sets['Case'] & core_sets['Control'])}")
    print(f"\nMachine Learning Results:")
    print(f"Mean CV Score: {ml_results['mean_cv_score']:.3f} ± {ml_results['std_cv_score']:.3f}")
    print(f"\nDifferential Abundance Results:")
    print(f"Number of significant taxa (FDR < 0.05): {len(da_results[da_results['FDR_P_value'] < 0.05])}")
    print(f"\nNetwork Analysis Results:")
    print(network_results[['Metric', 'Case', 'Control']].to_string(index=False))
    print(f"\nRare taxa summary:")
    print(rare_stats.to_string(index=False))
    
    # Print metabolic analysis results
    print(f"\nMetabolic Pathway Analysis Results:")
    sig_pathways = pathway_results[pathway_results['FDR_P_value'] < 0.05]
    print(f"Number of significantly different pathways (FDR < 0.05): {len(sig_pathways)}")
    if len(sig_pathways) > 0:
        print("Top significant pathways:")
        for _, row in sig_pathways.head().iterrows():
            print(f"  {row['Pathway']}: log2FC = {row['Log2FoldChange']:.2f}, FDR = {row['FDR_P_value']:.2e}")
    
    print(f"\nPathway-Tools integration:")
    if ptools_available:
        print("✓ Pathway-Tools is available for advanced metabolic analysis")
        print(f"✓ Analysis script created: {ptools_script}")
    else:
        print("⚠ Pathway-Tools not found - using simplified metabolic analysis")
        print(f"✓ Analysis script template created: {ptools_script}")
    
    if advanced_results:
        print(f"\nAdvanced Pathway-Tools Analysis Results:")
        print(f"✓ Advanced analysis completed successfully")
        print(f"✓ Detailed pathway, reaction, and compound data saved")
        print(f"✓ Advanced metabolic network created")
        print(f"✓ Advanced pathway enrichment analysis completed")
    else:
        print(f"\nAdvanced Pathway-Tools Analysis:")
        print(f"⚠ Advanced analysis not performed (Pathway-Tools Python API not available)")
        print(f"  To enable advanced features, install Pathway-Tools and the Python API")

if __name__ == "__main__":
    main() 
