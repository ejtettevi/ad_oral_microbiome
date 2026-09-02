# alpha_beta_analysis.py

import os
import pandas as pd
import numpy as np
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.distance import permanova, permdisp
from skbio.stats.ordination import pcoa
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, kruskal, ranksums
from statsmodels.stats.multitest import multipletests
import warnings
import sys

warnings.filterwarnings("ignore")
sns.set_theme(style="whitegrid", font_scale=1.0)
plt.rcParams['figure.dpi'] = 700

output_dir = "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/diversity"
os.makedirs(output_dir, exist_ok=True)

metadata_file = "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/metadata.tsv"
asv_table_file = "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/dada2/asv_table.tsv"
asv_taxonomy_file = "/Volumes/edKwamiBackUP/analysis_workflow/oral_project/base_dir/dada2/asv_taxonomy.tsv"

np.random.seed(42)

def add_pval_text_boxplot(ax, pval, n_groups):
    y_max = ax.get_ylim()[1]
    x_pos = (n_groups - 1) / 2
    if pval < 0.001:
        pval_str = "p < 0.001 ***"
    elif pval < 0.01:
        pval_str = f"p = {pval:.3f} **"
    elif pval < 0.05:
        pval_str = f"p = {pval:.3f} *"
    else:
        pval_str = f"p = {pval:.3f}"
    ax.text(x_pos, y_max * 0.98, pval_str, ha='center', va='top', fontsize=16, color='black', fontweight='bold')

def clean_names(names):
    return names.astype(str).str.strip().str.lower()

def clr_transform(df):
    df_transformed = df + 1e-6
    gm = np.exp(np.log(df_transformed).mean(axis=0))
    clr_df = np.log(df_transformed.divide(gm, axis=1))
    return clr_df

def chao1(counts):
    counts = np.array(counts)
    S_obs = np.sum(counts > 0)
    singletons = np.sum(counts == 1)
    doubletons = np.sum(counts == 2)
    if doubletons == 0:
        chao1_val = S_obs + (singletons * (singletons - 1)) / 2
    else:
        chao1_val = S_obs + (singletons ** 2) / (2 * doubletons)
    return chao1_val

def fisher_alpha(counts):
    counts = np.array(counts)
    S = np.sum(counts > 0)
    N = np.sum(counts)
    if S == 0 or N == 0:
        return np.nan
    from scipy.optimize import root_scalar
    def func(alpha):
        return alpha * np.log(1 + N / alpha) - S
    try:
        sol = root_scalar(func, bracket=[1e-6, 1000], method='bisect')
        return sol.root if sol.converged else np.nan
    except Exception:
        return np.nan

# --- Main Analysis Workflow ---

# --- Load DADA2 ASV table and taxonomy ---
try:
    asv_table = pd.read_csv(asv_table_file, sep='\t')
    if 'SampleID' in asv_table.columns:
        asv_table = asv_table.set_index('SampleID').T
    # Remove any unnamed columns
    asv_table = asv_table.loc[:, ~asv_table.columns.str.contains('^Unnamed')]
    print(f"Loaded ASV table: {asv_table.shape}")
except Exception as e:
    print(f"Error loading ASV table: {e}")
    sys.exit(1)

try:
    asv_tax = pd.read_csv(asv_taxonomy_file, sep='\t')
    print(f"Loaded ASV taxonomy: {asv_tax.shape}")
except Exception as e:
    print(f"Error loading ASV taxonomy: {e}")
    sys.exit(1)

# --- Load metadata and align samples ---
try:
    metadata = pd.read_csv(metadata_file, sep='\t')
    asv_table.columns = clean_names(asv_table.columns)
    metadata['Sample'] = clean_names(metadata['Sample'])
    metadata = metadata.set_index('Sample')
    common_samples = asv_table.columns.intersection(metadata.index)
    missing_in_metadata = set(asv_table.columns) - set(metadata.index)
    missing_in_counts = set(metadata.index) - set(asv_table.columns)
    if missing_in_metadata:
        print(f"Samples in ASV table but not in metadata: {sorted(missing_in_metadata)}")
    if missing_in_counts:
        print(f"Samples in metadata but not in ASV table: {sorted(missing_in_counts)}")
    print(f"Number of samples in both: {len(common_samples)}")
    if len(common_samples) == 0:
        print("ERROR: No samples are shared between ASV table and metadata. Exiting.")
        print(f"ASV table samples: {list(asv_table.columns)}")
        print(f"Metadata samples: {list(metadata.index)}")
        sys.exit(1)
    asv_table = asv_table.loc[:, common_samples]
    metadata = metadata.loc[common_samples]
except Exception as e:
    print(f"Error loading metadata: {e}")
    sys.exit(1)

counts = asv_table
counts_clr = clr_transform(counts)
print("Performed CLR transformation.")

# --- Alpha Diversity Calculations ---
shannon = alpha_diversity('shannon', counts.values.T, ids=counts.columns)
simpson = alpha_diversity('simpson', counts.values.T, ids=counts.columns)
observed = (counts > 0).sum(axis=0)
chao1_vals = [chao1(counts.iloc[:, i]) for i in range(counts.shape[1])]
fisher_vals = [fisher_alpha(counts.iloc[:, i]) for i in range(counts.shape[1])]

alpha_div_df = pd.DataFrame({
    'Sample': counts.columns,
    'Group': metadata['Group'],
    'Observed': observed,
    'Chao1': chao1_vals,
    'Shannon': shannon,
    'Simpson': simpson,
    'Fisher': fisher_vals
})
alpha_div_df.to_csv(f"{output_dir}/alpha_diversity_case_control.csv", index=False)
print("Calculated and saved alpha diversity.")

# Colorblind palette for all groups
n_groups = len(alpha_div_df['Group'].unique())
palette = sns.color_palette("colorblind", n_colors=max(n_groups, 2))
group_colors = dict(zip(sorted(alpha_div_df['Group'].unique()), palette))

# --- Alpha Diversity Violin Plots (only significant, PNG only, colorblind) ---
for metric in ['Observed', 'Chao1', 'Shannon', 'Simpson', 'Fisher']:
    groups = alpha_div_df['Group'].unique()
    vals = [alpha_div_df.loc[alpha_div_df['Group'] == g, metric].dropna() for g in groups]
    if len(groups) == 2:
        stat, pval = mannwhitneyu(vals[0], vals[1], alternative='two-sided')
    elif len(groups) > 2:
        stat, pval = kruskal(*vals)
    else:
        pval = np.nan

    if pval < 0.05:
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Create violin plot with better styling
        violin_parts = ax.violinplot([vals[i] for i in range(len(groups))], 
                                   positions=range(len(groups)), 
                                   widths=0.6, showmeans=True, showmedians=True)
        
        # Color the violin plots
        for i, pc in enumerate(violin_parts['bodies']):
            pc.set_facecolor(list(group_colors.values())[i])
            pc.set_alpha(0.7)
            pc.set_edgecolor('black')
            pc.set_linewidth(1.5)
        
        # Style the means and medians
        violin_parts['cmeans'].set_color('red')
        violin_parts['cmeans'].set_linewidth(2)
        violin_parts['cmedians'].set_color('darkred')
        violin_parts['cmedians'].set_linewidth(2)
        
        # Add individual data points with jitter
        for i, group in enumerate(groups):
            group_data = alpha_div_df.loc[alpha_div_df['Group'] == group, metric].dropna()
            x_jitter = np.random.normal(i, 0.1, len(group_data))
            ax.scatter(x_jitter, group_data, alpha=0.6, s=30, 
                      color=list(group_colors.values())[i], edgecolors='black', linewidth=0.5)
        
        # Add significance annotation with better positioning
        y_max = max([max(v) for v in vals if len(v) > 0])
        y_min = min([min(v) for v in vals if len(v) > 0])
        y_range = y_max - y_min
        y_text = y_max + 0.05 * y_range
        
        if pval < 0.001:
            pval_str = "p < 0.001 ***"
        elif pval < 0.01:
            pval_str = f"p = {pval:.3f} **"
        elif pval < 0.05:
            pval_str = f"p = {pval:.3f} *"
        else:
            pval_str = f"p = {pval:.3f}"
        
        ax.text(len(groups)/2 - 0.5, y_text, pval_str, ha='center', va='bottom', 
                fontsize=14, fontweight='bold', bbox=dict(boxstyle='round,pad=0.3', 
                facecolor='white', edgecolor='black', alpha=0.8))
        
        # Set labels and title with better spacing
        ax.set_title(f"{metric} Diversity by Group", fontsize=18, fontweight='bold', pad=20)
        ax.set_xlabel("Group", fontsize=16, fontweight='bold')
        ax.set_ylabel(f"{metric} Diversity", fontsize=16, fontweight='bold')
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, fontsize=14)
        ax.tick_params(axis='y', labelsize=12)
        
        # Add grid for better readability
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_axisbelow(True)
        
        # Create legend with better positioning
        handles = [
            plt.Line2D([0], [0], color=group_colors[g], lw=6, label=g)
            for g in sorted(alpha_div_df['Group'].unique())
        ]
        ax.legend(handles=handles, title="Group", loc='upper right', 
                 frameon=True, fancybox=True, shadow=True, fontsize=12)
        
        # Adjust layout to prevent label overlap
        plt.tight_layout()
        plt.subplots_adjust(top=0.9)
        plt.savefig(f"{output_dir}/{metric.lower()}_diversity_violin_colorblind.png", 
                   dpi=700, bbox_inches='tight', facecolor='white')
        plt.close()

# --- Beta Diversity Analyses ---
bray_curtis_dist = beta_diversity('braycurtis', counts.values.T, ids=counts.columns)
jaccard_dist = beta_diversity('jaccard', (counts > 0).values.T, ids=counts.columns)

# Save distance matrices
pd.DataFrame(bray_curtis_dist.data, index=bray_curtis_dist.ids, columns=bray_curtis_dist.ids).to_csv(f"{output_dir}/bray_curtis_distance_matrix.csv")
pd.DataFrame(jaccard_dist.data, index=jaccard_dist.ids, columns=jaccard_dist.ids).to_csv(f"{output_dir}/jaccard_distance_matrix.csv")

# --- Betadisper (homogeneity of dispersion) ---
def run_betadisper(dist_matrix, grouping):
    try:
        res = permdisp(dist_matrix, grouping)
        print(f"Betadisper p-value: {res['p-value']:.4g}")
        return res['p-value']
    except Exception as e:
        print(f"Betadisper failed: {e}")
        return None

betadisper_bray = run_betadisper(bray_curtis_dist, metadata['Group'])
betadisper_jaccard = run_betadisper(jaccard_dist, metadata['Group'])

# --- Ordination Analyses ---
ordination_results = {}

# Bray-Curtis PCoA
pcoa_res = pcoa(bray_curtis_dist)
ordination_results['bray_curtis'] = pcoa_res
pcoa_points = pcoa_res.samples.iloc[:, :2]
pcoa_df = pcoa_points.copy()
pcoa_df['Group'] = metadata['Group']
pcoa_df.to_csv(f"{output_dir}/pcoa_points_case_control.csv")
print("Calculated and saved Bray-Curtis PCoA points.")

# Jaccard PCoA
jaccard_pcoa = pcoa(jaccard_dist)
ordination_results['jaccard'] = jaccard_pcoa
jaccard_points = jaccard_pcoa.samples.iloc[:, :2]
jaccard_df = jaccard_points.copy()
jaccard_df['Group'] = metadata['Group']
jaccard_df.to_csv(f"{output_dir}/jaccard_pcoa_points_case_control.csv")
print("Calculated and saved Jaccard PCoA points.")

# --- NMDS Ordination (using scikit-learn) ---
def run_nmds(dist_matrix, ids, random_state=42):
    if np.isnan(dist_matrix.data).any():
        raise ValueError("Distance matrix contains NaN values. Cannot perform NMDS.")
    mds = MDS(n_components=2, metric=False, dissimilarity="precomputed", random_state=random_state, max_iter=300, n_init=10)
    nmds_res = mds.fit_transform(dist_matrix.data)
    stress = mds.stress_
    nmds_df = pd.DataFrame(nmds_res, columns=['NMDS1', 'NMDS2'], index=ids)
    return nmds_df, stress

try:
    bray_nmds_df, bray_stress = run_nmds(bray_curtis_dist, bray_curtis_dist.ids)
    bray_nmds_df['Group'] = metadata['Group']
    bray_nmds_df.to_csv(f"{output_dir}/bray_curtis_nmds_points_case_control.csv")
    print(f"Bray-Curtis NMDS stress: {bray_stress:.4f}")
except ValueError as e:
    print(f"Skipping Bray-Curtis NMDS: {e}")
    bray_nmds_df = None
    bray_stress = None

try:
    jaccard_nmds_df, jaccard_stress = run_nmds(jaccard_dist, jaccard_dist.ids)
    jaccard_nmds_df['Group'] = metadata['Group']
    jaccard_nmds_df.to_csv(f"{output_dir}/jaccard_nmds_points_case_control.csv")
    print(f"Jaccard NMDS stress: {jaccard_stress:.4f}")
except ValueError as e:
    print(f"Skipping Jaccard NMDS: {e}")
    jaccard_nmds_df = None
    jaccard_stress = None

print("Calculated and saved NMDS points (if possible).")

# --- PERMANOVA for Beta Diversity ---
def run_permanova(dist_matrix, grouping):
    if dist_matrix is not None and len(set(grouping)) > 1 and pd.Series(grouping).value_counts().min() > 1:
        try:
            res = permanova(distance_matrix=dist_matrix, grouping=grouping)
            return res['p-value']
        except Exception:
            return None
    return None

permanova_results = {
    'bray_curtis': run_permanova(bray_curtis_dist, metadata['Group']),
    'jaccard': run_permanova(jaccard_dist, metadata['Group'])
}

# --- Ordination Plots (PCoA and NMDS, only significant, PNG only, colorblind) ---
ordination_dfs = {
    'bray_curtis': pcoa_df,
    'jaccard': jaccard_df
}
for name, df in ordination_dfs.items():
    pval = permanova_results[name]
    betadisper_p = betadisper_bray if name == 'bray_curtis' else betadisper_jaccard
    if df is not None and pval is not None and pval < 0.05:
        axes = [col for col in df.columns if col not in ['Group']]
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create scatter plot with better styling
        for group, color in group_colors.items():
            sub = df[df['Group'] == group]
            ax.scatter(sub[axes[0]], sub[axes[1]], label=group, color=color, 
                      s=120, edgecolor='black', alpha=0.8, linewidth=1.5)
        
        # Add confidence ellipses for each group
        from matplotlib.patches import Ellipse
        for group, color in group_colors.items():
            sub = df[df['Group'] == group]
            if len(sub) > 2:  # Need at least 3 points for ellipse
                # Calculate ellipse parameters
                x_mean, y_mean = sub[axes[0]].mean(), sub[axes[1]].mean()
                x_std, y_std = sub[axes[0]].std(), sub[axes[1]].std()
                
                # Create ellipse
                ellipse = Ellipse((x_mean, y_mean), 2*x_std, 2*y_std, 
                                alpha=0.2, facecolor=color, edgecolor=color, linewidth=2)
                ax.add_patch(ellipse)
        
        explained_var = ordination_results[name].proportion_explained[:2] * 100
        
        # Create better title with improved formatting
        title = f"PCoA ({name.replace('_', ' ').title()})\n"
        title += f"Axis 1: {explained_var[0]:.1f}% variance explained, "
        title += f"Axis 2: {explained_var[1]:.1f}% variance explained\n"
        title += f"PERMANOVA p={pval:.3g}, Betadisper p={betadisper_p:.3g}"
        
        ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
        ax.set_xlabel(f"{axes[0]} ({explained_var[0]:.1f}%)", fontsize=14, fontweight='bold')
        ax.set_ylabel(f"{axes[1]} ({explained_var[1]:.1f}%)", fontsize=14, fontweight='bold')
        ax.tick_params(axis='both', labelsize=12)
        
        # Add grid for better readability
        ax.grid(True, alpha=0.3)
        ax.set_axisbelow(True)
        
        # Create legend with better positioning and styling
        leg = ax.legend(title="Group", loc='center left', bbox_to_anchor=(1, 0.5), 
                      frameon=True, fancybox=True, shadow=True, markerscale=1.5)
        leg.set_title("Group", prop={'size': 14, 'weight': 'bold'})
        for text in leg.get_texts():
            text.set_fontsize(12)
        
        # Add sample count annotation
        sample_counts = df['Group'].value_counts()
        count_text = "Sample counts:\n" + "\n".join([f"{group}: {count}" for group, count in sample_counts.items()])
        ax.text(0.02, 0.98, count_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', 
                facecolor='white', edgecolor='black', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{name}_pcoa_plot_colorblind.png", 
                   dpi=700, bbox_inches='tight', facecolor='white')
        plt.close()

# NMDS plots
nmds_dfs = {
    'bray_curtis': (bray_nmds_df, bray_stress),
    'jaccard': (jaccard_nmds_df, jaccard_stress)
}
for name, (df, stress) in nmds_dfs.items():
    pval = permanova_results[name]
    betadisper_p = betadisper_bray if name == 'bray_curtis' else betadisper_jaccard
    if df is not None and pval is not None and pval < 0.05:
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create scatter plot with better styling
        for group, color in group_colors.items():
            sub = df[df['Group'] == group]
            ax.scatter(sub['NMDS1'], sub['NMDS2'], label=group, color=color, 
                      s=120, edgecolor='black', alpha=0.8, linewidth=1.5)
        
        # Add confidence ellipses for each group
        from matplotlib.patches import Ellipse
        for group, color in group_colors.items():
            sub = df[df['Group'] == group]
            if len(sub) > 2:  # Need at least 3 points for ellipse
                # Calculate ellipse parameters
                x_mean, y_mean = sub['NMDS1'].mean(), sub['NMDS2'].mean()
                x_std, y_std = sub['NMDS1'].std(), sub['NMDS2'].std()
                
                # Create ellipse
                ellipse = Ellipse((x_mean, y_mean), 2*x_std, 2*y_std, 
                                alpha=0.2, facecolor=color, edgecolor=color, linewidth=2)
                ax.add_patch(ellipse)
        
        # Create better title with improved formatting
        title = f"NMDS ({name.replace('_', ' ').title()})\n"
        title += f"Stress: {stress:.4f}\n"
        title += f"PERMANOVA p={pval:.3g}, Betadisper p={betadisper_p:.3g}"
        
        ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
        ax.set_xlabel("NMDS1", fontsize=14, fontweight='bold')
        ax.set_ylabel("NMDS2", fontsize=14, fontweight='bold')
        ax.tick_params(axis='both', labelsize=12)
        
        # Add grid for better readability
        ax.grid(True, alpha=0.3)
        ax.set_axisbelow(True)
        
        # Create legend with better positioning and styling
        leg = ax.legend(title="Group", loc='center left', bbox_to_anchor=(1, 0.5), 
                      frameon=True, fancybox=True, shadow=True, markerscale=1.5)
        leg.set_title("Group", prop={'size': 14, 'weight': 'bold'})
        for text in leg.get_texts():
            text.set_fontsize(12)
        
        # Add sample count annotation
        sample_counts = df['Group'].value_counts()
        count_text = "Sample counts:\n" + "\n".join([f"{group}: {count}" for group, count in sample_counts.items()])
        ax.text(0.02, 0.98, count_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', 
                facecolor='white', edgecolor='black', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{name}_nmds_plot_colorblind.png", 
                   dpi=700, bbox_inches='tight', facecolor='white')
        plt.close()

# --- Differential Abundance Analysis (Wilcoxon rank-sum, FDR correction, volcano plot) ---
# NOTE: For compositional DA, use ANCOM-BC/ALDEx2 results from R (see step 7)
if len(alpha_div_df['Group'].unique()) == 2:
    group1, group2 = sorted(alpha_div_df['Group'].unique())
    group1_samples = alpha_div_df.loc[alpha_div_df['Group'] == group1, 'Sample']
    group2_samples = alpha_div_df.loc[alpha_div_df['Group'] == group2, 'Sample']
    results = []
    for asv in counts.index:
        vals1 = counts.loc[asv, group1_samples]
        vals2 = counts.loc[asv, group2_samples]
        if (vals1.sum() > 0 or vals2.sum() > 0):
            stat, pval = ranksums(vals1, vals2)
            mean1 = vals1.mean()
            mean2 = vals2.mean()
            log2fc = np.log2((mean2 + 1e-6) / (mean1 + 1e-6))
            results.append({'ASV': asv, 'log2FC': log2fc, 'pval': pval})
    results_df = pd.DataFrame(results)
    results_df['padj'] = multipletests(results_df['pval'], method='fdr_bh')[1]
    results_df.to_csv(f"{output_dir}/differential_abundance_wilcoxon.csv", index=False)
    print("Differential abundance analysis complete (Wilcoxon, for sensitivity only).")

    # Enhanced Volcano plot (colorblind, robust to group count)
    sig = (results_df['padj'] < 0.05)
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Define colors for different significance levels
    colors = []
    sizes = []
    for idx, row in results_df.iterrows():
        if row['padj'] < 0.001:
            colors.append('#d62728')  # Dark red for highly significant
            sizes.append(80)
        elif row['padj'] < 0.01:
            colors.append('#ff7f0e')  # Orange for very significant
            sizes.append(60)
        elif row['padj'] < 0.05:
            colors.append('#2ca02c')  # Green for significant
            sizes.append(40)
        else:
            colors.append('#cccccc')  # Gray for non-significant
            sizes.append(20)
    
    # Create scatter plot with varying sizes and colors
    scatter = ax.scatter(results_df['log2FC'], -np.log10(results_df['padj']),
                        c=colors, s=sizes, alpha=0.7, edgecolors='black', linewidth=0.5)
    
    # Add significance threshold lines
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=2, alpha=0.7, label='FDR = 0.05')
    ax.axhline(-np.log10(0.01), color='black', linestyle=':', linewidth=1.5, alpha=0.7, label='FDR = 0.01')
    ax.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
    
    # Add labels for top significant points to avoid overlap
    top_significant = results_df[results_df['padj'] < 0.05].nsmallest(10, 'padj')
    for idx, row in top_significant.iterrows():
        ax.annotate(f"{row['ASV'][:8]}...", 
                   (row['log2FC'], -np.log10(row['padj'])),
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=8, alpha=0.8,
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white', 
                            edgecolor='black', alpha=0.7))
    
    # Enhanced labels and title
    ax.set_xlabel('log₂ Fold Change', fontsize=16, fontweight='bold')
    ax.set_ylabel('-log₁₀ Adjusted p-value', fontsize=16, fontweight='bold')
    ax.set_title('Differential Abundance Volcano Plot\n(Wilcoxon Rank-Sum Test with FDR Correction)', 
                fontsize=18, fontweight='bold', pad=20)
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3)
    ax.set_axisbelow(True)
    
    # Create enhanced legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#d62728', label='p < 0.001 (***)'),
        Patch(facecolor='#ff7f0e', label='p < 0.01 (**)'),
        Patch(facecolor='#2ca02c', label='p < 0.05 (*)'),
        Patch(facecolor='#cccccc', label='p ≥ 0.05 (ns)'),
        plt.Line2D([0], [0], color='black', linestyle='--', label='FDR = 0.05'),
        plt.Line2D([0], [0], color='black', linestyle=':', label='FDR = 0.01')
    ]
    
    ax.legend(handles=legend_elements, loc='upper right', frameon=True, 
             fancybox=True, shadow=True, fontsize=12)
    
    # Add statistics annotation
    n_sig = sig.sum()
    n_total = len(results_df)
    stats_text = f"Total ASVs: {n_total}\nSignificant: {n_sig} ({n_sig/n_total*100:.1f}%)"
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', 
            facecolor='white', edgecolor='black', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/differential_abundance_volcano_colorblind.png", 
               dpi=700, bbox_inches='tight', facecolor='white')
    plt.close()

# --- Save environment info for reproducibility ---
with open(f"{output_dir}/python_env_info.txt", "w") as f:
    import platform
    f.write("Python version: " + platform.python_version() + "\n")
    f.write("Platform: " + platform.platform() + "\n")
    f.write("Random seed: 42\n")
    f.write("Packages:\n")
    for pkg in ['skbio', 'scikit-learn', 'pandas', 'seaborn', 'matplotlib', 'statsmodels']:
        try:
            mod = __import__(pkg.replace('-', '_'))
            f.write(f"{pkg}: {mod.__version__}\n")
        except Exception:
            f.write(f"{pkg}: not found\n")

print("\nAnalysis complete. All results, colorblind violin plots, ordination plots, and volcano plot saved to output directory.")
