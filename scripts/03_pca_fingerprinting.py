import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import os

BASE = '/mnt/e/IBD-metagenomics-project3'
MATRIX = os.path.join(BASE, 'results', 'species_abundance_matrix.tsv')
OUT_FIG = os.path.join(BASE, 'figures')
OUT_RES = os.path.join(BASE, 'results')

# ── Load ──────────────────────────────────────────────────────
df = pd.read_csv(MATRIX, sep='\t', index_col='sample_id')
diagnoses = df['diagnosis']
X = df.drop(columns=['diagnosis'])

# ── Preprocessing ─────────────────────────────────────────────
# Log transform: reduces dominance of highly abundant species
# Add small pseudocount to avoid log(0)
X_log = np.log1p(X)

# Standardize: each species has mean=0, std=1
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_log)

# ── PCA ───────────────────────────────────────────────────────
pca = PCA()
X_pca = pca.fit_transform(X_scaled)

explained = pca.explained_variance_ratio_ * 100
cumulative = np.cumsum(explained)

print("PCA Explained Variance:")
for i, (e, c) in enumerate(zip(explained, cumulative)):
    print(f"  PC{i+1}: {e:.1f}%  (cumulative: {c:.1f}%)")

# PC loadings: which species drive each PC
loadings = pd.DataFrame(
    pca.components_.T,
    index=X.columns,
    columns=[f'PC{i+1}' for i in range(len(explained))]
)

print(f"\nTop species driving PC1 (IBD axis):")
print(loadings['PC1'].abs().sort_values(ascending=False).head(10))
print(f"\nTop species driving PC2:")
print(loadings['PC2'].abs().sort_values(ascending=False).head(10))

loadings.to_csv(os.path.join(OUT_RES, 'pca_loadings.tsv'), sep='\t')

# PCA coordinates per sample
coords = pd.DataFrame(X_pca[:, :4],
                      index=X.index,
                      columns=['PC1','PC2','PC3','PC4'])
coords['diagnosis'] = diagnoses
coords.to_csv(os.path.join(OUT_RES, 'pca_coordinates.tsv'), sep='\t')

# ── Plot ──────────────────────────────────────────────────────
COLORS = {'nonIBD': '#2ecc71', 'UC': '#e74c3c', 'CD': '#e67e22'}
MARKERS = {'nonIBD': 'o', 'UC': 's', 'CD': '^'}

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle('PCA Microbiome Fingerprinting — IBD vs nonIBD', fontsize=14, fontweight='bold')

# Panel 1: PC1 vs PC2
ax = axes[0]
for sample in X.index:
    diag = diagnoses[sample]
    ax.scatter(coords.loc[sample, 'PC1'], coords.loc[sample, 'PC2'],
               c=COLORS[diag], marker=MARKERS[diag], s=200, zorder=5,
               edgecolors='black', linewidths=0.8)
    ax.annotate(sample, (coords.loc[sample, 'PC1'], coords.loc[sample, 'PC2']),
                textcoords='offset points', xytext=(8, 4), fontsize=8)

ax.set_xlabel(f'PC1 ({explained[0]:.1f}% variance)', fontsize=11)
ax.set_ylabel(f'PC2 ({explained[1]:.1f}% variance)', fontsize=11)
ax.set_title('PC1 vs PC2', fontsize=11)
ax.axhline(0, color='grey', linewidth=0.5, linestyle='--')
ax.axvline(0, color='grey', linewidth=0.5, linestyle='--')
ax.grid(True, alpha=0.3)
handles = [mpatches.Patch(color=COLORS[d], label=d) for d in COLORS]
ax.legend(handles=handles, fontsize=9)

# Panel 2: PC1 vs PC3
ax = axes[1]
for sample in X.index:
    diag = diagnoses[sample]
    ax.scatter(coords.loc[sample, 'PC1'], coords.loc[sample, 'PC3'],
               c=COLORS[diag], marker=MARKERS[diag], s=200, zorder=5,
               edgecolors='black', linewidths=0.8)
    ax.annotate(sample, (coords.loc[sample, 'PC1'], coords.loc[sample, 'PC3']),
                textcoords='offset points', xytext=(8, 4), fontsize=8)

ax.set_xlabel(f'PC1 ({explained[0]:.1f}% variance)', fontsize=11)
ax.set_ylabel(f'PC3 ({explained[2]:.1f}% variance)', fontsize=11)
ax.set_title('PC1 vs PC3', fontsize=11)
ax.axhline(0, color='grey', linewidth=0.5, linestyle='--')
ax.axvline(0, color='grey', linewidth=0.5, linestyle='--')
ax.grid(True, alpha=0.3)
ax.legend(handles=handles, fontsize=9)

# Panel 3: Scree plot (variance explained)
ax = axes[2]
pcs = range(1, len(explained) + 1)
ax.bar(pcs, explained, color='#3498db', alpha=0.7, label='Per-PC variance')
ax.plot(pcs, cumulative, 'o-', color='#e74c3c', linewidth=2, label='Cumulative variance')
ax.axhline(80, color='grey', linestyle='--', alpha=0.5, label='80% threshold')
ax.set_xlabel('Principal Component', fontsize=11)
ax.set_ylabel('Variance Explained (%)', fontsize=11)
ax.set_title('Scree Plot', fontsize=11)
ax.set_xticks(list(pcs))
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
out_path = os.path.join(OUT_FIG, 'pca_fingerprinting.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"\nFigure saved to {out_path}")
plt.close()
