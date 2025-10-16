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
        fig, ax = plt.subplots(figsize=(7, 5))
        sns.violinplot(
            x='Group', y=metric, data=alpha_div_df,
            palette=group_colors, inner='box', linewidth=2, ax=ax
        )
        add_pval_text_boxplot(ax, pval, n_groups=len(groups))
        ax.set_title(f"{metric} Diversity by Group (Violin)", fontsize=16, fontweight='bold')
        ax.set_xlabel("Group", fontsize=16, fontweight='bold')
        ax.set_ylabel(f"{metric} Diversity", fontsize=16, fontweight='bold')
        ax.tick_params(axis='both', labelsize=14)
        handles = [
            plt.Line2D([0], [0], color=group_colors[g], lw=8, label=g)
            for g in sorted(alpha_div_df['Group'].unique())
        ]
        ax.legend(handles=handles, title="Group", loc='best', frameon=True)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{metric.lower()}_diversity_violin_colorblind.png", dpi=700, bbox_inches='tight')
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
        fig, ax = plt.subplots(figsize=(7, 6))
        for group, color in group_colors.items():
            sub = df[df['Group'] == group]
            ax.scatter(sub[axes[0]], sub[axes[1]], label=group, color=color, s=80, edgecolor='k', alpha=0.85)
        explained_var = ordination_results[name].proportion_explained[:2] * 100
        title = f"PCoA ({name.replace('_', ' ').title()})\nAxis 1: {explained_var[0]:.1f}%, Axis 2: {explained_var[1]:.1f}%"
        title += f"\nPERMANOVA p={pval:.3g}, Betadisper p={betadisper_p:.3g}"
        ax.set_title(title, fontsize=16, fontweight='bold')
        ax.set_xlabel(axes[0], fontsize=16, fontweight='bold')
        ax.set_ylabel(axes[1], fontsize=16, fontweight='bold')
        ax.tick_params(axis='both', labelsize=14)
        ax.grid(False)
        leg = ax.legend(title="Group", loc='center left', bbox_to_anchor=(1, 0.5), frameon=True, markerscale=1.2)
        leg.set_title("Group", prop={'size': 14, 'weight': 'bold'})
        for text in leg.get_texts():
            text.set_fontsize(14)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{name}_pcoa_plot_colorblind.png", dpi=700, bbox_inches='tight')
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
        fig, ax = plt.subplots(figsize=(7, 6))
        for group, color in group_colors.items():
            sub = df[df['Group'] == group]
            ax.scatter(sub['NMDS1'], sub['NMDS2'], label=group, color=color, s=80, edgecolor='k', alpha=0.85)
        title = f"NMDS ({name.replace('_', ' ').title()})\nStress: {stress:.4f}\nPERMANOVA p={pval:.3g}, Betadisper p={betadisper_p:.3g}"
        ax.set_title(title, fontsize=16, fontweight='bold')
        ax.set_xlabel("NMDS1", fontsize=16, fontweight='bold')
        ax.set_ylabel("NMDS2", fontsize=16, fontweight='bold')
        ax.tick_params(axis='both', labelsize=14)
        ax.grid(False)
        leg = ax.legend(title="Group", loc='center left', bbox_to_anchor=(1, 0.5), frameon=True, markerscale=1.2)
        leg.set_title("Group", prop={'size': 14, 'weight': 'bold'})
        for text in leg.get_texts():
            text.set_fontsize(14)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{name}_nmds_plot_colorblind.png", dpi=700, bbox_inches='tight')
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

    # Volcano plot (colorblind, robust to group count)
    sig = (results_df['padj'] < 0.05)
    fig, ax = plt.subplots(figsize=(8, 6))
    sig_color = palette[0]
    nonsig_color = palette[1] if len(palette) > 1 else "#cccccc"
    ax.scatter(
        results_df['log2FC'], -np.log10(results_df['padj']),
        c=[sig_color if s else nonsig_color for s in sig], s=30, alpha=0.8, edgecolor='none', label=None
    )
    ax.set_xlabel('log2 Fold Change', fontsize=16, fontweight='bold')
    ax.set_ylabel('-log10 Adjusted p-value', fontsize=16, fontweight='bold')
    ax.set_title('Differential Abundance Volcano Plot (Wilcoxon, FDR)', fontsize=18, fontweight='bold')
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1.5, label='FDR=0.05')
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=sig_color, edgecolor='none', label='Significant (FDR<0.05)'),
        Patch(facecolor=nonsig_color, edgecolor='none', label='Not Significant'),
        Patch(facecolor='white', edgecolor='black', linestyle='--', label='FDR=0.05')
    ]
    ax.legend(handles=legend_elements[:2], loc='best', frameon=True)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/differential_abundance_volcano_colorblind.png", dpi=700, bbox_inches='tight')
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
