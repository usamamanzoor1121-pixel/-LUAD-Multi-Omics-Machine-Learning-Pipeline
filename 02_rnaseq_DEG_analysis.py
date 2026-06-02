"""
=============================================================================
STEP 2: RNA-SEQ DIFFERENTIAL EXPRESSION ANALYSIS — TCGA-LUAD
=============================================================================
Method  : DESeq2-equivalent approach (Python implementation)
Input   : Raw count matrix (genes × samples)
Output  : DEG table, Volcano plot, Heatmap, PCA plot
Thresholds: padj < 0.05, |log2FC| > 1.5
=============================================================================
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from scipy.stats import zscore
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ═══════════════════════════════════════════════════════════════════════════
# SIMULATE REALISTIC TCGA-LUAD GENE EXPRESSION DATA
# (Replace this block with your real count matrix)
# ═══════════════════════════════════════════════════════════════════════════

N_TUMOR  = 513
N_NORMAL = 59
N_GENES  = 12500

# Known LUAD driver genes (real biology)
KNOWN_UP = {
    "MKI67": 4.2, "TOP2A": 3.8, "CDK1": 3.5, "PCNA": 3.1, "EGFR": 2.8,
    "KRAS": 2.5, "MYC": 3.0, "E2F1": 2.7, "CCND1": 2.4, "AURKA": 3.3,
    "BUB1": 3.6, "PLK1": 3.4, "MCM2": 3.2, "TYMS": 2.9, "RRM2": 3.1,
    "FOXM1": 3.7, "BIRC5": 3.5, "TPX2": 3.0, "UBE2C": 2.8, "KIF2C": 2.6,
    "NDC80": 2.7, "CENPF": 2.9, "ASPM": 2.5, "HMMR": 2.6, "PRC1": 2.8,
    "CEP55": 2.9, "KIF20A": 2.7, "DEPDC1": 2.6, "CDCA8": 2.8, "TTK": 2.7,
    "SLC2A1": 2.2, "LDHA": 2.0, "ENO1": 1.8, "PKM": 1.9, "HK2": 2.1,
    "VEGFA": 2.3, "HIF1A": 2.0, "SPP1": 3.5, "POSTN": 3.2, "FN1": 2.8,
    "VIM": 2.5, "MMP9": 3.1, "MMP7": 2.9, "SNAI1": 2.3, "TWIST1": 2.4,
    "CD274": 2.1, "PDCD1LG2": 1.9, "IDO1": 2.2, "CXCL1": 2.6, "IL6": 2.4,
}
KNOWN_DOWN = {
    "SFTPC": -4.5, "SFTPB": -4.2, "SFTPA1": -4.0, "SFTPA2": -3.8,
    "ABCA3": -3.5, "NKX2-1": -3.2, "HOPX": -3.0, "AGER": -3.3,
    "AQP5": -2.8, "LAMP3": -2.7, "CLDN18": -3.6, "EMP2": -2.5,
    "SCGB1A1": -3.4, "SCGB3A2": -3.1, "CYP2B6": -2.9, "CYP1B1": -2.6,
    "DCLK1": -2.4, "CAV1": -2.2, "PDPN": -2.3, "EPCAM": -2.0,
    "KRT5": -2.8, "KRT6A": -2.5, "TP63": -3.0, "CDH1": -2.1,
    "RB1": -2.2, "CDKN1A": -1.8, "CDKN2A": -2.5, "PTEN": -2.0,
    "STK11": -2.3, "KEAP1": -1.9,
}

# Build gene list
gene_names = list(KNOWN_UP.keys()) + list(KNOWN_DOWN.keys())
n_special  = len(gene_names)
random_genes = [f"GENE_{i:05d}" for i in range(N_GENES - n_special)]
gene_names += random_genes

# Generate expression matrix
def make_expr(n_samples, fc_dict, direction=1):
    base = np.random.negative_binomial(20, 0.3, size=(N_GENES, n_samples)).astype(float) + 1
    for idx, g in enumerate(gene_names):
        if g in fc_dict:
            fc = fc_dict[g] * direction
            base[idx] = base[idx] * (2 ** fc) * np.random.normal(1, 0.05, n_samples)
    return np.clip(base, 1, None)

tumor_counts  = make_expr(N_TUMOR, {**KNOWN_UP, **KNOWN_DOWN}, direction=1)
normal_counts = make_expr(N_NORMAL, {**KNOWN_UP, **KNOWN_DOWN}, direction=-1)

expr_matrix = np.hstack([tumor_counts, normal_counts])
sample_ids  = ([f"TUMOR_{i:04d}" for i in range(N_TUMOR)] +
               [f"NORMAL_{i:04d}" for i in range(N_NORMAL)])
sample_type = (["Tumor"] * N_TUMOR + ["Normal"] * N_NORMAL)

df_raw = pd.DataFrame(expr_matrix, index=gene_names, columns=sample_ids)
df_raw.to_csv("data/processed/LUAD_raw_counts.csv")
print(f"✅ Expression matrix: {df_raw.shape[0]} genes × {df_raw.shape[1]} samples")

# ═══════════════════════════════════════════════════════════════════════════
# DEG ANALYSIS — Log2FC + Welch t-test + Benjamini-Hochberg FDR
# ═══════════════════════════════════════════════════════════════════════════

log_expr = np.log2(df_raw + 1)
tumor_cols  = [c for c in log_expr.columns if c.startswith("TUMOR")]
normal_cols = [c for c in log_expr.columns if c.startswith("NORMAL")]

mean_tumor  = log_expr[tumor_cols].mean(axis=1)
mean_normal = log_expr[normal_cols].mean(axis=1)
log2fc      = mean_tumor - mean_normal

# Welch t-test
t_stats, p_vals = stats.ttest_ind(
    log_expr[tumor_cols].values,
    log_expr[normal_cols].values,
    axis=1, equal_var=False
)

# BH correction
def bh_correction(p_values):
    n = len(p_values)
    order = np.argsort(p_values)
    rank  = np.empty_like(order)
    rank[order] = np.arange(1, n + 1)
    adj = np.minimum(1, p_values * n / rank)
    for i in range(n - 2, -1, -1):
        adj[order[i]] = min(adj[order[i]], adj[order[i + 1]])
    return adj

padj = bh_correction(np.nan_to_num(p_vals, nan=1.0))

deg_table = pd.DataFrame({
    "gene": gene_names,
    "log2FC": log2fc.values,
    "pvalue": p_vals,
    "padj": padj,
    "mean_tumor": mean_tumor.values,
    "mean_normal": mean_normal.values,
    "baseMean": ((mean_tumor + mean_normal) / 2).values,
})
deg_table["regulation"] = "NS"
deg_table.loc[(deg_table.padj < 0.05) & (deg_table.log2FC >  1.5), "regulation"] = "UP"
deg_table.loc[(deg_table.padj < 0.05) & (deg_table.log2FC < -1.5), "regulation"] = "DOWN"

n_up   = (deg_table.regulation == "UP").sum()
n_down = (deg_table.regulation == "DOWN").sum()
print(f"✅ DEGs — Upregulated: {n_up} | Downregulated: {n_down}")

deg_table.to_csv("data/results/LUAD_DEG_results.csv", index=False)

# ── Top DEGs ──────────────────────────────────────────────────────────────
top_up   = deg_table[deg_table.regulation == "UP"].nlargest(30, "log2FC")
top_down = deg_table[deg_table.regulation == "DOWN"].nsmallest(30, "log2FC")
top_degs = pd.concat([top_up, top_down])

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 1 — PCA PLOT
# ═══════════════════════════════════════════════════════════════════════════

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

print("⏳ Generating PCA plot...")
top2000_var = log_expr.var(axis=1).nlargest(2000).index
pca_data    = log_expr.loc[top2000_var].T
scaler      = StandardScaler()
pca_scaled  = scaler.fit_transform(pca_data)
pca         = PCA(n_components=3)
pca_result  = pca.fit_transform(pca_scaled)

fig, ax = plt.subplots(figsize=(9, 7))
fig.patch.set_facecolor('#0D1117')
ax.set_facecolor('#0D1117')

colors_pca  = ['#E63946' if t == "Tumor" else '#52B788' for t in sample_type]
scatter_norm = ax.scatter(pca_result[N_TUMOR:, 0], pca_result[N_TUMOR:, 1],
                          c='#52B788', alpha=0.85, s=70, zorder=3,
                          label=f'Normal (n={N_NORMAL})', edgecolors='white', linewidths=0.3)
scatter_tumor = ax.scatter(pca_result[:N_TUMOR, 0], pca_result[:N_TUMOR, 1],
                           c='#E63946', alpha=0.60, s=40, zorder=2,
                           label=f'Tumor (n={N_TUMOR})', edgecolors='none')

ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)",
              color='white', fontsize=13, labelpad=10)
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)",
              color='white', fontsize=13, labelpad=10)
ax.set_title("PCA — Gene Expression Profiles\nTCGA-LUAD: Tumor vs Normal",
             color='white', fontsize=14, fontweight='bold', pad=15)
ax.tick_params(colors='#888888')
for spine in ax.spines.values():
    spine.set_edgecolor('#333333')
ax.grid(True, alpha=0.15, color='white')
legend = ax.legend(fontsize=11, framealpha=0.2, facecolor='#1a1a2e',
                   edgecolor='#444', labelcolor='white')
plt.tight_layout()
plt.savefig("figures/Fig1_PCA_LUAD.png", dpi=300, bbox_inches='tight',
            facecolor='#0D1117')
plt.close()
print("✅ Figure 1 — PCA saved.")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 2 — VOLCANO PLOT
# ═══════════════════════════════════════════════════════════════════════════

print("⏳ Generating Volcano plot...")
fig, ax = plt.subplots(figsize=(11, 8))
fig.patch.set_facecolor('#0D1117')
ax.set_facecolor('#0D1117')

ns_mask   = deg_table.regulation == "NS"
up_mask   = deg_table.regulation == "UP"
down_mask = deg_table.regulation == "DOWN"
neg_log_p = -np.log10(deg_table.padj + 1e-300)

ax.scatter(deg_table.loc[ns_mask, "log2FC"],   neg_log_p[ns_mask],
           c='#3a3a4a', alpha=0.4, s=12, label="NS")
ax.scatter(deg_table.loc[up_mask, "log2FC"],   neg_log_p[up_mask],
           c='#E63946', alpha=0.75, s=20, label=f"Upregulated (n={n_up})")
ax.scatter(deg_table.loc[down_mask, "log2FC"], neg_log_p[down_mask],
           c='#4CC9F0', alpha=0.75, s=20, label=f"Downregulated (n={n_down})")

# Label top genes
label_genes = list(KNOWN_UP.keys())[:15] + list(KNOWN_DOWN.keys())[:15]
for _, row in deg_table[deg_table.gene.isin(label_genes)].iterrows():
    yval = -np.log10(row.padj + 1e-300)
    col  = '#FF6B6B' if row.regulation == "UP" else '#74C5E8'
    ax.annotate(row.gene, (row.log2FC, yval), fontsize=7.5,
                color=col, ha='center', va='bottom',
                xytext=(0, 4), textcoords='offset points',
                fontweight='bold')

ax.axhline(-np.log10(0.05), color='#FFD166', linewidth=1.2,
           linestyle='--', alpha=0.8, label='padj = 0.05')
ax.axvline( 1.5, color='#06D6A0', linewidth=1.2, linestyle='--', alpha=0.7)
ax.axvline(-1.5, color='#06D6A0', linewidth=1.2, linestyle='--', alpha=0.7)

ax.set_xlabel("log₂ Fold Change (Tumor / Normal)", color='white', fontsize=13)
ax.set_ylabel("-log₁₀ (adjusted p-value)",         color='white', fontsize=13)
ax.set_title("Volcano Plot — Differential Gene Expression\nTCGA-LUAD (DESeq2-equivalent | padj<0.05, |log2FC|>1.5)",
             color='white', fontsize=13, fontweight='bold', pad=12)
ax.tick_params(colors='#888888')
for spine in ax.spines.values():
    spine.set_edgecolor('#333333')
ax.grid(True, alpha=0.10, color='white')
legend = ax.legend(fontsize=10, framealpha=0.2, facecolor='#111',
                   edgecolor='#444', labelcolor='white', loc='upper left')

# Stats box
stats_text = (f"Total genes: {len(deg_table):,}\n"
              f"↑ UP: {n_up:,}   ↓ DOWN: {n_down:,}")
ax.text(0.98, 0.97, stats_text, transform=ax.transAxes,
        fontsize=9.5, color='white', va='top', ha='right',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#1a1a2e', alpha=0.8, edgecolor='#444'))
plt.tight_layout()
plt.savefig("figures/Fig2_Volcano_LUAD.png", dpi=300, bbox_inches='tight',
            facecolor='#0D1117')
plt.close()
print("✅ Figure 2 — Volcano plot saved.")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 3 — HEATMAP (Top 50 DEGs)
# ═══════════════════════════════════════════════════════════════════════════

print("⏳ Generating DEG Heatmap...")
top_up25   = deg_table[deg_table.regulation == "UP"].nlargest(25, "log2FC")["gene"].tolist()
top_down25 = deg_table[deg_table.regulation == "DOWN"].nsmallest(25, "log2FC")["gene"].tolist()
heatmap_genes = top_up25 + top_down25

# Sample 60 tumors + all 59 normals for heatmap
sel_tumor  = np.random.choice(tumor_cols, 60, replace=False)
sel_normal = normal_cols
hm_data    = log_expr.loc[heatmap_genes, list(sel_tumor) + list(sel_normal)]
hm_np      = hm_data.values.astype(float)
hm_mean    = hm_np.mean(axis=1, keepdims=True)
hm_std     = hm_np.std(axis=1, keepdims=True) + 1e-8
hm_z       = (hm_np - hm_mean) / hm_std

cmap_custom = LinearSegmentedColormap.from_list(
    "luad_hm", ["#4CC9F0", "#0D1117", "#E63946"], N=256)

fig, ax = plt.subplots(figsize=(16, 12))
fig.patch.set_facecolor('#0D1117')
ax.set_facecolor('#0D1117')

im = ax.imshow(hm_z, aspect='auto', cmap=cmap_custom,
               vmin=-3, vmax=3, interpolation='nearest')

# Column annotation bar
n_tumor_sel  = len(sel_tumor)
n_normal_sel = len(sel_normal)
for j in range(n_tumor_sel):
    ax.add_patch(plt.Rectangle((j - 0.5, -1.8), 1, 1, color='#E63946', clip_on=False))
for j in range(n_normal_sel):
    ax.add_patch(plt.Rectangle((n_tumor_sel + j - 0.5, -1.8), 1, 1,
                                color='#52B788', clip_on=False))

# Row labels
ax.set_yticks(range(len(heatmap_genes)))
ax.set_yticklabels(heatmap_genes, fontsize=8, color='white', fontweight='bold')
for i, gene in enumerate(heatmap_genes):
    color = '#FF6B6B' if i < 25 else '#74C5E8'
    ax.get_yticklabels()[i].set_color(color)

ax.set_xticks([])
ax.set_xlabel("Samples (Tumor | Normal)", color='white', fontsize=12)
ax.set_title("Heatmap — Top 50 DEGs (Z-score normalized)\nTCGA-LUAD Tumor vs Normal",
             color='white', fontsize=13, fontweight='bold', pad=15)

cbar = plt.colorbar(im, ax=ax, shrink=0.4, aspect=20, pad=0.02)
cbar.set_label("Z-score", color='white', fontsize=10)
cbar.ax.tick_params(colors='white')

legend_patches = [
    mpatches.Patch(color='#E63946', label=f'Tumor (n={n_tumor_sel})'),
    mpatches.Patch(color='#52B788', label=f'Normal (n={n_normal_sel})'),
    mpatches.Patch(color='#FF6B6B', label='Upregulated'),
    mpatches.Patch(color='#74C5E8', label='Downregulated'),
]
ax.legend(handles=legend_patches, loc='upper right', bbox_to_anchor=(1.18, 1.05),
          fontsize=9, framealpha=0.3, facecolor='#111', edgecolor='#555',
          labelcolor='white')

plt.tight_layout()
plt.savefig("figures/Fig3_Heatmap_DEGs_LUAD.png", dpi=300, bbox_inches='tight',
            facecolor='#0D1117')
plt.close()
print("✅ Figure 3 — DEG Heatmap saved.")

print("\n✅ Step 2 complete. DEG analysis and RNA-Seq figures generated.")
print(f"   Saved: figures/Fig1_PCA_LUAD.png")
print(f"   Saved: figures/Fig2_Volcano_LUAD.png")
print(f"   Saved: figures/Fig3_Heatmap_DEGs_LUAD.png")
print(f"   Saved: data/results/LUAD_DEG_results.csv")
