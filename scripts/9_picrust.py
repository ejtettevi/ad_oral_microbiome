import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import requests
import subprocess
import sys
import shutil
import re
from bs4 import BeautifulSoup
import time
from scipy.stats import mannwhitneyu, ttest_ind
from scipy.stats import fisher_exact
import matplotlib.patches as mpatches

# === Plot Annotation Helpers ===
def annotate_pvalue(ax, x, y, p, yoffset_frac=0.06, fontsize=11):
    """Place significance stars above a point/bar with automatic y-limit adjustment"""
    try:
        p_val = float(p)
        if p_val < 0.001:
            p_text = "***"
        elif p_val < 0.01:
            p_text = "**"
        elif p_val < 0.05:
            p_text = "*"
        else:
            p_text = "ns"
    except Exception:
        p_text = "ns"
    ylims = ax.get_ylim()
    span = ylims[1] - ylims[0] if (ylims[1] - ylims[0]) != 0 else max(1.0, y)
    y_text = y + yoffset_frac * span
    if y_text > ylims[1]*0.98:
        ax.set_ylim(ylims[0], y_text * 1.15)
    ax.text(x, y_text, p_text, ha='center', va='bottom', fontweight='bold', fontsize=fontsize)
    ax.margins(y=0.1)

# =============================================================================
# CONFIGURATION SECTION
# =============================================================================

print("="*80)
print("COMPREHENSIVE PICRUSt2 FUNCTIONAL ANALYSIS PIPELINE")
print("="*80)

# File paths - Update these according to your directory structure
base_dir = "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir"
input_tsv = os.path.join(base_dir, "dada2/asv_table.tsv")
output_dir = os.path.join(base_dir, "picrust2_out")
downstream_dir = os.path.join(base_dir, "picrust2_analysis")
fasta_file = os.path.join(base_dir, "dada2/rep_seqs.fasta")
metadata_tsv = os.path.join(base_dir, "metadata.tsv")

# Derived file paths
transposed_tsv = os.path.join(output_dir, "asv_table_transposed.tsv")
biom_table = os.path.join(output_dir, "asv_table_transposed.biom")
pathway_tsv = os.path.join(output_dir, "pathways_out/path_abun_unstrat.tsv")

# Analysis parameters
TOP_PATHWAYS_COUNT = 15  # Number of top pathways to analyze
STATISTICAL_THRESHOLD = 0.05  # P-value threshold for significance
FIGURE_DPI = 300  # Resolution for saved figures

# Ensure output directories exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(downstream_dir, exist_ok=True)

# =============================================================================
# STEP 1: DATA PREPROCESSING (from _6_picrust2.py)
# =============================================================================

print("\n1. DATA PREPROCESSING")
print("-" * 40)

# Check if transposed table already exists
if not os.path.exists(transposed_tsv):
    print("Transposing ASV table...")
    try:
        df_asv = pd.read_csv(input_tsv, sep='\t', index_col=0)
        df_t = df_asv.transpose()
        df_t.index.name = "#OTU ID"
        df_t.to_csv(transposed_tsv, sep='\t')
        print(f"✓ Transposed table saved to {transposed_tsv}")
    except Exception as e:
        print(f"✗ Error transposing ASV table: {e}")
        sys.exit(1)
else:
    print("✓ Transposed ASV table already exists")

# Check if BIOM table already exists
if not os.path.exists(biom_table):
    print("Converting transposed TSV to BIOM format...")
    convert_cmd = [
        "biom", "convert",
        "-i", transposed_tsv,
        "-o", biom_table,
        "--to-hdf5",
        "--table-type=OTU table"
    ]
    try:
        subprocess.run(convert_cmd, check=True)
        print(f"✓ Converted {transposed_tsv} to {biom_table}")
    except subprocess.CalledProcessError as e:
        print(f"✗ Error converting TSV to BIOM: {e}")
        sys.exit(1)
else:
    print("✓ BIOM table already exists")

# =============================================================================
# STEP 2: PICRUSt2 EXECUTION
# =============================================================================

print("\n2. RUNNING PICRUSt2 PIPELINE")
print("-" * 40)

if not os.path.exists(pathway_tsv):
    if os.path.exists(output_dir) and os.path.exists(os.path.join(output_dir, "pathways_out")):
        print(f"Removing existing PICRUSt2 output directory: {output_dir}")
        shutil.rmtree(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        
        # Re-create preprocessing files
        print("Re-creating preprocessing files...")
        df_asv = pd.read_csv(input_tsv, sep='\t', index_col=0)
        df_t = df_asv.transpose()
        df_t.index.name = "#OTU ID"
        df_t.to_csv(transposed_tsv, sep='\t')
        
        subprocess.run([
            "biom", "convert",
            "-i", transposed_tsv,
            "-o", biom_table,
            "--to-hdf5",
            "--table-type=OTU table"
        ], check=True)
    
    print("Running PICRUSt2 to generate pathway abundance table...")
    picrust2_cmd = [
        "picrust2_pipeline.py",
        "-s", fasta_file,
        "-i", biom_table,
        "-o", output_dir,
        "-p", "4"
    ]
    try:
        subprocess.run(picrust2_cmd, check=True)
        print("✓ PICRUSt2 run completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"✗ Error running PICRUSt2: {e}")
        sys.exit(1)
else:
    print("✓ Pathway abundance file already exists. Skipping PICRUSt2 run.")

# =============================================================================
# STEP 3: DATA LOADING AND PREPROCESSING
# =============================================================================

print("\n3. LOADING AND PREPROCESSING DATA")
print("-" * 40)

# Read pathway abundance data
try:
    df = pd.read_csv(pathway_tsv, sep='\t', index_col=0)
    print(f"✓ Loaded pathway abundance data: {df.shape}")
except Exception as e:
    print(f"✗ Error reading pathway abundance file: {e}")
    sys.exit(1)

# Read metadata (auto-detect delimiter)
try:
    meta = pd.read_csv(metadata_tsv, sep='\t')
    print(f"✓ Loaded metadata with tab delimiter")
except Exception:
    try:
        meta = pd.read_csv(metadata_tsv)
        print(f"✓ Loaded metadata with comma delimiter")
    except Exception as e:
        print(f"✗ Error reading metadata file: {e}")
        sys.exit(1)

# Set up metadata index
if 'Sample' in meta.columns:
    meta = meta.set_index('Sample')
elif 'sample' in meta.columns:
    meta = meta.set_index('sample')
else:
    raise ValueError("Metadata file must have a column named 'Sample' or 'sample'.")

# Ensure data consistency
df = df.loc[:, meta.index.intersection(df.columns)]
meta = meta.loc[df.columns]

print(f"✓ Data aligned: {df.shape} pathways × {len(df.columns)} samples")
print(f"✓ Group distribution: {meta['Group'].value_counts().to_dict()}")

# =============================================================================
# STEP 4: PATHWAY SELECTION AND ANNOTATION
# =============================================================================

print("\n4. PATHWAY SELECTION AND ANNOTATION")
print("-" * 40)

# Get top pathways by abundance
pathway_sums = df.sum(axis=1)
top_pathways = pathway_sums.sort_values(ascending=False).head(TOP_PATHWAYS_COUNT)
top_pathway_ids = top_pathways.index.tolist()
df_top = df.loc[top_pathway_ids]

print(f"✓ Selected top {TOP_PATHWAYS_COUNT} pathways by abundance")

# Pathway annotation functions
def normalize_pwy_id(pid):
    """Normalize pathway IDs to standard MetaCyc format"""
    pid = str(pid).strip()
    match = re.match(r'PWY0-(\d+)', pid)
    if match:
        return f"PWY-{int(match.group(1)):05d}"
    return pid

def get_metacyc_pathway_name(pwy_id, max_retries=3):
    """Query MetaCyc for pathway name with retry logic"""
    url = f"https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object={pwy_id}"
    headers = {"User-Agent": "Mozilla/5.0"}
    
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, headers=headers, timeout=15)
            if resp.status_code != 200:
                print(f"  Warning: Failed to fetch {pwy_id}: HTTP {resp.status_code}")
                return None
            
            soup = BeautifulSoup(resp.text, "html.parser")
            
            # Try different methods to extract pathway name
            title = soup.find("title")
            if title and "Pathway" in title.text:
                name = title.text.split(" Pathway")[0].strip()
                if name and name != pwy_id:
                    return name
            
            h1 = soup.find("h1")
            if h1 and h1.text.strip() != pwy_id:
                return h1.text.strip()
            
            h2 = soup.find("h2")
            if h2 and h2.text.strip() != pwy_id:
                return h2.text.strip()
            
            return None
            
        except Exception as e:
            if attempt < max_retries - 1:
                print(f"  Retry {attempt + 1} for {pwy_id} due to: {e}")
                time.sleep(2)
            else:
                print(f"  Error fetching {pwy_id} after {max_retries} attempts: {e}")
                return None
    
    return None

# Build pathway annotation mapping
print("Annotating pathways with MetaCyc database...")
pathway_label_map = {}
annotated_count = 0

for i, pid in enumerate(df_top.index, 1):
    print(f"  Processing {i}/{len(df_top.index)}: {pid}")
    
    pid_norm = normalize_pwy_id(pid)
    name = get_metacyc_pathway_name(pid_norm)
    
    if name:
        pathway_label_map[pid] = f"{pid_norm}: {name}"
        annotated_count += 1
    else:
        pathway_label_map[pid] = pid
    
    time.sleep(1)  # Be respectful to MetaCyc server

print(f"✓ Successfully annotated {annotated_count}/{len(df_top.index)} pathways")

# Apply annotations to dataframe
df_annotated = df_top.copy()
df_annotated.index = [pathway_label_map[pid] for pid in df_annotated.index]

# =============================================================================
# STEP 5: GROUP SEPARATION AND STATISTICAL ANALYSIS
# =============================================================================

print("\n5. STATISTICAL ANALYSIS")
print("-" * 40)

# Separate groups
case_samples = meta[meta['Group'] == 'Case'].index
control_samples = meta[meta['Group'] == 'Control'].index

print(f"✓ Case samples: {len(case_samples)}")
print(f"✓ Control samples: {len(control_samples)}")

if len(case_samples) == 0 or len(control_samples) == 0:
    print("✗ Error: Need samples from both groups for comparison")
    sys.exit(1)

# Perform statistical tests
print("Performing statistical tests...")
results = []

for pathway in df_annotated.index:
    original_pid = [k for k, v in pathway_label_map.items() if v == pathway][0]
    case_values = df_top.loc[original_pid, case_samples]
    control_values = df_top.loc[original_pid, control_samples]
    
    # Mann-Whitney U test (non-parametric)
    try:
        stat, pval = mannwhitneyu(case_values, control_values, alternative='two-sided')
    except:
        pval = 1.0
        stat = np.nan
    
    # Calculate fold change
    case_mean = case_values.mean()
    control_mean = control_values.mean()
    fold_change = case_mean / (control_mean + 1e-10)  # Avoid division by zero
    log2_fc = np.log2(fold_change + 1e-10)
    
    # Effect size (Cohen's d approximation)
    pooled_std = np.sqrt(((case_values.std()**2) + (control_values.std()**2)) / 2)
    cohens_d = (case_mean - control_mean) / (pooled_std + 1e-10)
    
    results.append({
        'Pathway': pathway,
        'Original_ID': original_pid,
        'Case_Mean': case_mean,
        'Control_Mean': control_mean,
        'Case_Std': case_values.std(),
        'Control_Std': control_values.std(),
        'Fold_Change': fold_change,
        'Log2_FC': log2_fc,
        'Cohens_D': cohens_d,
        'Mann_Whitney_U': stat,
        'P_value': pval,
        'Neg_Log10_P': -np.log10(pval + 1e-10),
        'Significant': pval < STATISTICAL_THRESHOLD
    })

results_df = pd.DataFrame(results)
results_df = results_df.sort_values('P_value')

# Save comprehensive statistical results
stats_file = os.path.join(downstream_dir, "comprehensive_statistical_results.tsv")
results_df.to_csv(stats_file, sep='\t', index=False)

n_significant = results_df['Significant'].sum()
n_higher_case = ((results_df['Significant']) & (results_df['Log2_FC'] > 0)).sum()
n_higher_control = ((results_df['Significant']) & (results_df['Log2_FC'] < 0)).sum()

print(f"✓ Statistical analysis complete:")
print(f"  - Significantly different pathways: {n_significant}/{len(results_df)}")
print(f"  - Higher in Case: {n_higher_case}")
print(f"  - Higher in Control: {n_higher_control}")

# =============================================================================
# STEP 6: COMPREHENSIVE VISUALIZATIONS
# =============================================================================

print("\n6. CREATING COMPREHENSIVE VISUALIZATIONS")
print("-" * 40)

# Set style for better plots
plt.style.use('default')
sns.set_palette("husl")

# --- Plot 1: Enhanced Group Comparison Boxplot ---
print("Creating enhanced group comparison boxplot...")
plt.figure(figsize=(18, 10))

# Prepare data for boxplot
plot_data = []
for pathway in df_annotated.index:
    original_pid = [k for k, v in pathway_label_map.items() if v == pathway][0]
    pathway_short = pathway[:60] + '...' if len(pathway) > 60 else pathway
    
    for sample in case_samples:
        plot_data.append({
            'Pathway': pathway_short,
            'Group': 'Case', 
            'Abundance': df_top.loc[original_pid, sample]
        })
    for sample in control_samples:
        plot_data.append({
            'Pathway': pathway_short,
            'Group': 'Control', 
            'Abundance': df_top.loc[original_pid, sample]
        })

plot_df = pd.DataFrame(plot_data)

ax = sns.boxplot(data=plot_df, x='Pathway', y='Abundance', hue='Group', 
                 palette=['lightcoral', 'skyblue'])
plt.title(f'Top {TOP_PATHWAYS_COUNT} Pathways: Case vs Control Comparison', 
          fontsize=18, fontweight='bold', pad=20)
plt.xlabel('Pathway', fontsize=14, fontweight='bold')
plt.ylabel('Abundance', fontsize=14, fontweight='bold')
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.legend(title='Group', fontsize=12, title_fontsize=12)

# Add group legend
group_legend = plt.legend(title='Group', fontsize=12, title_fontsize=12, loc='upper right')
plt.gca().add_artist(group_legend)

# Add significance legend in the middle
sig_legend_elements = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=0, label='*** p < 0.001'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=0, label='** p < 0.01'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=0, label='* p < 0.05'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=0, label='ns p ≥ 0.05')
]
plt.legend(handles=sig_legend_elements, title='Significance', loc='center right', bbox_to_anchor=(1.0, 0.5))

# Add significance annotations (exact p-values)
for i, pathway in enumerate(df_annotated.index):
    pval = results_df[results_df['Pathway'] == pathway]['P_value'].iloc[0]
    y_max = plot_df[plot_df['Pathway'] == plot_df['Pathway'].unique()[i]]['Abundance'].max()
    annotate_pvalue(ax, i, y_max, pval)

plt.tight_layout()
plt.savefig(os.path.join(downstream_dir, "enhanced_boxplot_comparison.png"), 
            dpi=FIGURE_DPI, bbox_inches='tight')
plt.close()

# --- Plot 2: Enhanced Volcano Plot ---
print("Creating volcano plot...")
plt.figure(figsize=(14, 10))

# Create color mapping based on significance and fold change
colors = []
for idx, row in results_df.iterrows():
    if row['P_value'] < 0.001:
        colors.append('darkred' if row['Log2_FC'] > 0 else 'darkblue')
    elif row['P_value'] < 0.01:
        colors.append('red' if row['Log2_FC'] > 0 else 'blue')
    elif row['P_value'] < 0.05:
        colors.append('lightcoral' if row['Log2_FC'] > 0 else 'lightblue')
    else:
        colors.append('lightgray')

plt.scatter(results_df['Log2_FC'], results_df['Neg_Log10_P'], 
           c=colors, alpha=0.7, s=60)

# Add significance threshold lines
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5, label='p=0.05')
plt.axhline(y=-np.log10(0.01), color='black', linestyle=':', alpha=0.5, label='p=0.01')
plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)

# Add pathway labels for most significant ones
top_significant = results_df.head(8)
for idx, row in top_significant.iterrows():
    pathway_short = row['Pathway'][:40] + '...' if len(row['Pathway']) > 40 else row['Pathway']
    plt.annotate(pathway_short, (row['Log2_FC'], row['Neg_Log10_P']), 
                xytext=(5, 5), textcoords='offset points', fontsize=9,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))

plt.xlabel('Log₂ Fold Change (Case/Control)', fontsize=14, fontweight='bold')
plt.ylabel('-Log₁₀ P-value', fontsize=14, fontweight='bold')
plt.title('Pathway Differential Abundance: Case vs Control', fontsize=16, fontweight='bold')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(downstream_dir, "enhanced_volcano_plot.png"), 
            dpi=FIGURE_DPI, bbox_inches='tight')
plt.close()

# --- Plot 3: Annotated Heatmap with Enhanced Clustering ---
print("Creating annotated heatmap...")
plt.figure(figsize=(16, 12))

# Create group color annotation
group_colors = meta['Group'].map({'Case': 'red', 'Control': 'blue'})

# Create significance annotation for pathways
pathway_sig_colors = []
for pathway in df_annotated.index:
    pval = results_df[results_df['Pathway'] == pathway]['P_value'].iloc[0]
    if pval < 0.001:
        pathway_sig_colors.append('darkgreen')
    elif pval < 0.01:
        pathway_sig_colors.append('green')
    elif pval < 0.05:
        pathway_sig_colors.append('lightgreen')
    else:
        pathway_sig_colors.append('white')

# Create the clustermap
g = sns.clustermap(df_annotated, 
                   col_colors=group_colors,
                   row_colors=pathway_sig_colors,
                   cmap='viridis', 
                   figsize=(16, 12),
                   cbar_kws={'label': 'Abundance'},
                   linewidths=0.1)

g.ax_col_dendrogram.set_title(f'Top {TOP_PATHWAYS_COUNT} Annotated Pathways Heatmap', 
                              fontsize=16, fontweight='bold', pad=20)

# Add Case-Control legend
case_patch = mpatches.Patch(color='red', label='Case')
control_patch = mpatches.Patch(color='blue', label='Control')
g.fig.legend(handles=[case_patch, control_patch], 
            title='Groups', loc='upper right', bbox_to_anchor=(0.98, 0.98))

# Add significance legend
sig1_patch = mpatches.Patch(color='darkgreen', label='p < 0.001')
sig2_patch = mpatches.Patch(color='green', label='p < 0.01')
sig3_patch = mpatches.Patch(color='lightgreen', label='p < 0.05')
ns_patch = mpatches.Patch(color='white', label='p ≥ 0.05')
g.fig.legend(handles=[sig1_patch, sig2_patch, sig3_patch, ns_patch], 
            title='Significance', loc='lower left', bbox_to_anchor=(0.02, 0.02))
plt.savefig(os.path.join(downstream_dir, "annotated_heatmap.png"), 
            dpi=FIGURE_DPI, bbox_inches='tight')
plt.close()

# --- Plot 4: Enhanced Mean Abundance Comparison ---
print("Creating mean abundance comparison...")
original_ids = [results_df[results_df['Pathway'] == p]['Original_ID'].iloc[0] for p in df_annotated.index]
case_means = df_top.loc[original_ids, case_samples].mean(axis=1)
control_means = df_top.loc[original_ids, control_samples].mean(axis=1)

comparison_df = pd.DataFrame({
    'Pathway': df_annotated.index,
    'Case_Mean': case_means.values,
    'Control_Mean': control_means.values,
    'P_value': [results_df[results_df['Pathway'] == p]['P_value'].iloc[0] for p in df_annotated.index]
})

comparison_df['Pathway_Short'] = [p[:50] + '...' if len(p) > 50 else p for p in comparison_df['Pathway']]

plt.figure(figsize=(16, 10))
x = np.arange(len(comparison_df))
width = 0.35

bars1 = plt.bar(x - width/2, comparison_df['Case_Mean'], width, 
                label='Case', color='lightcoral', alpha=0.8)
bars2 = plt.bar(x + width/2, comparison_df['Control_Mean'], width, 
                label='Control', color='skyblue', alpha=0.8)

# Add exact p-values above bars
for i, pval in enumerate(comparison_df['P_value']):
    y_max = max(comparison_df.iloc[i]['Case_Mean'], comparison_df.iloc[i]['Control_Mean'])
    annotate_pvalue(plt.gca(), i, y_max, pval)

plt.xlabel('Pathway', fontsize=14, fontweight='bold')
plt.ylabel('Mean Abundance', fontsize=14, fontweight='bold')
plt.title(f'Mean Pathway Abundance Comparison: Case vs Control', fontsize=16, fontweight='bold')
plt.xticks(x, comparison_df['Pathway_Short'], rotation=45, ha='right')

# Add group legend
group_legend = plt.legend(fontsize=12, loc='upper right')
plt.gca().add_artist(group_legend)

# Add significance legend
sig_legend_elements = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=0, label='*** p < 0.001'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=0, label='** p < 0.01'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=0, label='* p < 0.05'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=0, label='ns p ≥ 0.05')
]
plt.legend(handles=sig_legend_elements, title='Significance', loc='center right', bbox_to_anchor=(0.98, 0.5))

plt.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(downstream_dir, "enhanced_mean_comparison.png"), 
            dpi=FIGURE_DPI, bbox_inches='tight')
plt.close()

# --- Plot 5: Top Significant Pathways Individual Analysis ---
print("Creating individual pathway analysis...")
top_significant = results_df[results_df['Significant']].head(6)

if len(top_significant) > 0:
    fig, axes = plt.subplots(2, 3, figsize=(24, 16))
    axes = axes.ravel()
    
    for i, (idx, pathway_data) in enumerate(top_significant.iterrows()):
        if i >= 6:
            break
            
        pathway = pathway_data['Pathway']
        original_pid = pathway_data['Original_ID']
        pathway_short = pathway[:35] + '...' if len(pathway) > 35 else pathway
        
        case_values = df_top.loc[original_pid, case_samples]
        control_values = df_top.loc[original_pid, control_samples]
        
        # Create enhanced violin plot with individual points
        data_for_violin = [control_values.values, case_values.values]
        parts = axes[i].violinplot(data_for_violin, positions=[1, 2], widths=0.7)
        
        # Color the violin plots
        parts['bodies'][0].set_facecolor('skyblue')
        parts['bodies'][0].set_alpha(0.7)
        parts['bodies'][1].set_facecolor('lightcoral')
        parts['bodies'][1].set_alpha(0.7)
        
        # Add individual data points
        axes[i].scatter([1] * len(control_values), control_values, 
                       color='blue', alpha=0.6, s=30)
        axes[i].scatter([2] * len(case_values), case_values, 
                       color='red', alpha=0.6, s=30)
        
        axes[i].set_xticks([1, 2])
        axes[i].set_xticklabels(['Control', 'Case'], fontsize=14, fontweight='bold')
        axes[i].set_ylabel('Abundance', fontsize=14, fontweight='bold')
        
        # Enhanced title with more statistics
        title = f'{pathway_short}\np={pathway_data["P_value"]:.2e}, FC={pathway_data["Fold_Change"]:.2f}'
        axes[i].set_title(title, fontsize=12, fontweight='bold', pad=15)
        axes[i].tick_params(axis='y', labelsize=12)
        axes[i].grid(True, alpha=0.3)
        
        # Add mean lines
        axes[i].axhline(y=control_values.mean(), xmin=0.15, xmax=0.35, 
                       color='blue', linewidth=2, alpha=0.8)
        axes[i].axhline(y=case_values.mean(), xmin=0.65, xmax=0.85, 
                       color='red', linewidth=2, alpha=0.8)
    
    # Hide unused subplots
    for i in range(len(top_significant), 6):
        axes[i].set_visible(False)
    
    plt.suptitle('Top 6 Most Significantly Different Pathways (Individual Analysis)', 
                 fontsize=24, fontweight='bold', y=0.95)
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.savefig(os.path.join(downstream_dir, "top_significant_individual.png"), 
                dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
else:
    print("  No significantly different pathways found for individual analysis.")

# =============================================================================
# STEP 7: GENERATE COMPREHENSIVE SUMMARY REPORT
# =============================================================================

print("\n7. GENERATING COMPREHENSIVE REPORT")
print("-" * 40)

# Create detailed summary report
summary = f"""
COMPREHENSIVE PICRUSt2 FUNCTIONAL ANALYSIS REPORT
================================================

Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}
Analysis Directory: {downstream_dir}

DATASET INFORMATION
-------------------
• Total pathways in dataset: {len(df)}
• Top pathways analyzed: {TOP_PATHWAYS_COUNT}
• Total samples: {len(df.columns)}
• Case samples: {len(case_samples)}
• Control samples: {len(control_samples)}

PATHWAY ANNOTATION
------------------
• Successfully annotated pathways: {annotated_count}/{len(df_top.index)}
• Annotation success rate: {(annotated_count/len(df_top.index)*100):.1f}%

STATISTICAL RESULTS
-------------------
• Significance threshold: p < {STATISTICAL_THRESHOLD}
• Significantly different pathways: {n_significant}/{len(results_df)}
• Pathways higher in Case: {n_higher_case}
• Pathways higher in Control: {n_higher_control}

TOP 10 MOST SIGNIFICANT PATHWAYS
---------------------------------
"""

for i, (idx, row) in enumerate(results_df.head(10).iterrows(), 1):
    direction = "Higher in Case" if row['Log2_FC'] > 0 else "Higher in Control"
    effect_size = "Large" if abs(row['Cohens_D']) > 0.8 else "Medium" if abs(row['Cohens_D']) > 0.5 else "Small"
    
    summary += f"\n{i:2d}. {row['Pathway'][:80]}"
    if len(row['Pathway']) > 80:
        summary += "..."
    summary += f"\n    • P-value: {row['P_value']:.2e}"
    summary += f"\n    • Fold Change: {row['Fold_Change']:.2f}"
    summary += f"\n    • Effect Size (Cohen's d): {row['Cohens_D']:.2f} ({effect_size})"
    summary += f"\n    • Direction: {direction}"
    summary += f"\n    • Case Mean ± SD: {row['Case_Mean']:.2f} ± {row['Case_Std']:.2f}"
    summary += f"\n    • Control Mean ± SD: {row['Control_Mean']:.2f} ± {row['Control_Std']:.2f}\n"

summary += f"""

FILES GENERATED
---------------
• enhanced_boxplot_comparison.png - Group comparison boxplot with significance
• enhanced_volcano_plot.png - Volcano plot showing fold change vs significance
• annotated_heatmap.png - Clustered heatmap with group and significance annotations
• enhanced_mean_comparison.png - Mean abundance comparison with significance
• top_significant_individual.png - Individual analysis of most significant pathways
• comprehensive_statistical_results.tsv - Complete statistical analysis results
• analysis_summary_report.txt - This summary report

METHODS SUMMARY
---------------
• Data preprocessing: ASV table transposition and BIOM conversion
• Functional prediction: PICRUSt2 pipeline
• Statistical test: Mann-Whitney U test (non-parametric)
• Multiple testing: Not corrected (consider FDR correction for multiple comparisons)
• Effect size: Cohen's d
• Visualization: Enhanced plots with statistical annotations
• Pathway annotation: MetaCyc database integration

RECOMMENDATIONS
---------------
1. Consider applying multiple testing correction (FDR/Bonferroni) for more stringent analysis
2. Validate top findings with additional functional analysis methods
3. Consider pathway enrichment analysis for biological interpretation
4. Examine individual pathway contributors for mechanistic insights

Analysis completed successfully!
"""

# Save summary report
summary_file = os.path.join(downstream_dir, "analysis_summary_report.txt")
with open(summary_file, 'w') as f:
    f.write(summary)

print("✓ Comprehensive summary report generated")

# =============================================================================
# FINAL OUTPUT SUMMARY
# =============================================================================

print("\n" + "="*80)
print("ANALYSIS COMPLETED SUCCESSFULLY!")
print("="*80)
print(f"Output directory: {downstream_dir}")
print(f"Total files generated: 7")
print("\nGenerated files:")
print("  📊 enhanced_boxplot_comparison.png")
print("  🌋 enhanced_volcano_plot.png") 
print("  🔥 annotated_heatmap.png")
print("  📈 enhanced_mean_comparison.png")
print("  🎯 top_significant_individual.png")
print("  📋 comprehensive_statistical_results.tsv")
print("  📄 analysis_summary_report.txt")

print(f"\n🧬 FUNCTIONAL ANALYSIS SUMMARY:")
print(f"  • Analyzed {TOP_PATHWAYS_COUNT} top pathways")
print(f"  • Found {n_significant} significantly different pathways (p < {STATISTICAL_THRESHOLD})")

if __name__ == "__main__":
    print("\n🚀 Pipeline executed successfully!")
    print("Check the output directories for results and visualizations.")
