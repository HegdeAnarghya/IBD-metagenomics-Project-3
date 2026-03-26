import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import networkx as nx
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.metrics import roc_curve, roc_auc_score
from scipy.stats import spearmanr
import os

BASE = '/mnt/e/IBD-metagenomics-project3'
OUT_FIG = os.path.join(BASE, 'figures')
os.makedirs(OUT_FIG, exist_ok=True)

COLORS = {'nonIBD': '#2ecc71', 'UC': '#e74c3c', 'CD': '#e67e22'}

# ── Load data ─────────────────────────────────────────────────
df5  = pd.read_csv(os.path.join(BASE, 'results', 'species_abundance_matrix.tsv'),
                   sep='\t', index_col='sample_id')
df50 = pd.read_csv(os.path.join(BASE, 'results', 'species_matrix_50samples.tsv'),
                   sep='\t', index_col='sample_id')
metrics = pd.read_csv(os.path.join(BASE, 'results', 'network_metrics.tsv'),
                      sep='\t', index_col='species')
importances = pd.read_csv(os.path.join(BASE, 'results', 'feature_importances_50.tsv'),
                          sep='\t', index_col=0, names=['species','importance'],
                          header=0).squeeze()

# ── Figure layout ─────────────────────────────────────────────
fig = plt.figure(figsize=(20, 14))
fig.patch.set_facecolor('#0f0f0f')
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.32)

ax1 = fig.add_subplot(gs[0, 0])  # Species heatmap (5 samples)
ax2 = fig.add_subplot(gs[0, 1])  # PCA fingerprint (5 samples)
ax3 = fig.add_subplot(gs[0, 2])  # Network hub species
ax4 = fig.add_subplot(gs[1, 0])  # PCA 50 samples
ax5 = fig.add_subplot(gs[1, 1])  # Feature importances
ax6 = fig.add_subplot(gs[1, 2])  # ROC + AUC comparison

for ax in [ax1,ax2,ax3,ax4,ax5,ax6]:
    ax.set_facecolor('#1a1a2e')
    for spine in ax.spines.values():
        spine.set_color('#444466')

title_kw  = dict(color='white', fontsize=11, fontweight='bold', pad=8)
label_kw  = dict(color='#cccccc', fontsize=9)
tick_kw   = dict(colors='#aaaaaa', labelsize=8)

# ── Panel 1: Species heatmap (top 20, 5 samples) ─────────────
diag5 = df5['diagnosis']
X5 = df5.drop(columns=['diagnosis'])
top20 = X5.mean().sort_values(ascending=False).head(20).index
heatmap_data = X5[top20].T

im = ax1.imshow(heatmap_data.values, aspect='auto', cmap='YlOrRd',
                interpolation='nearest')
ax1.set_xticks(range(5))
ax1.set_xticklabels([f"{s}\n({diag5[s]})" for s in X5.index],
                    fontsize=7, color='#cccccc', rotation=20)
ax1.set_yticks(range(20))
ax1.set_yticklabels([s.replace('_',' ') for s in top20],
                    fontsize=7, color='#cccccc', style='italic')
ax1.set_title('Top 20 Species Abundance\n(5 samples)', **title_kw)
plt.colorbar(im, ax=ax1, label='Relative Abundance (%)',
             fraction=0.03).ax.yaxis.set_tick_params(colors='#aaaaaa',
             labelsize=7)

# ── Panel 2: PCA fingerprint (5 samples) ─────────────────────
X5_log = np.log1p(X5)
X5_sc  = StandardScaler().fit_transform(X5_log)
pca5   = PCA(n_components=3)
coords5 = pca5.fit_transform(X5_sc)
ev5 = pca5.explained_variance_ratio_ * 100

for i, sample in enumerate(X5.index):
    diag = diag5[sample]
    ax2.scatter(coords5[i,0], coords5[i,1],
                c=COLORS[diag], s=220, zorder=5,
                edgecolors='white', linewidths=0.8)
    ax2.annotate(sample, (coords5[i,0], coords5[i,1]),
                 textcoords='offset points', xytext=(7,4),
                 fontsize=7, color='white')

ax2.set_xlabel(f'PC1 ({ev5[0]:.1f}%)', **label_kw)
ax2.set_ylabel(f'PC2 ({ev5[1]:.1f}%)', **label_kw)
ax2.set_title('PCA Microbiome Fingerprint\n(5 samples)', **title_kw)
ax2.tick_params(**tick_kw)
ax2.axhline(0, color='#444466', lw=0.5); ax2.axvline(0, color='#444466', lw=0.5)
handles2 = [mpatches.Patch(color=COLORS[d], label=d) for d in COLORS]
ax2.legend(handles=handles2, fontsize=8,
           facecolor='#1a1a2e', labelcolor='white', edgecolor='#444466')

# ── Panel 3: Network hub species ─────────────────────────────
top15_hubs = metrics.head(15)
ibd_dep = {'Faecalibacterium_prausnitzii','Akkermansia_muciniphila','Roseburia_intestinalis'}
ibd_enr = {'Bacteroides_vulgatus','Ruminococcus_gnavus','Bacteroides_uniformis'}
hub_colors = ['#e74c3c' if s in ibd_enr else '#2ecc71' if s in ibd_dep
              else '#7f8c8d' for s in top15_hubs.index]

bars3 = ax3.barh(range(len(top15_hubs)), top15_hubs['degree'],
                 color=hub_colors, edgecolor='#1a1a2e', height=0.7)
ax3.set_yticks(range(len(top15_hubs)))
ax3.set_yticklabels([s.replace('_',' ') for s in top15_hubs.index],
                    fontsize=7, color='#cccccc', style='italic')
ax3.invert_yaxis()
ax3.set_xlabel('Network Degree (connections)', **label_kw)
ax3.set_title('Co-occurrence Network Hubs\n(|r|≥0.9, 5 samples)', **title_kw)
ax3.tick_params(**tick_kw)
ax3.grid(True, alpha=0.2, axis='x', color='#444466')
h3 = [mpatches.Patch(color='#e74c3c', label='IBD-enriched'),
      mpatches.Patch(color='#2ecc71', label='IBD-depleted'),
      mpatches.Patch(color='#7f8c8d', label='Other')]
ax3.legend(handles=h3, fontsize=7, facecolor='#1a1a2e',
           labelcolor='white', edgecolor='#444466')

# ── Panel 4: PCA 50 samples ───────────────────────────────────
diag50 = df50['diagnosis']
X50 = df50.drop(columns=['diagnosis'])
X50 = X50.loc[:, (X50 > 0).sum() >= 5]
X50_log = np.log1p(X50)
X50_sc  = StandardScaler().fit_transform(X50_log)
pca50   = PCA(n_components=2)
coords50 = pca50.fit_transform(X50_sc)
ev50 = pca50.explained_variance_ratio_ * 100

for i, sample in enumerate(X50.index):
    diag = diag50[sample]
    ax4.scatter(coords50[i,0], coords50[i,1],
                c=COLORS[diag], s=80, alpha=0.85, zorder=5,
                edgecolors='white', linewidths=0.4)

ax4.set_xlabel(f'PC1 ({ev50[0]:.1f}%)', **label_kw)
ax4.set_ylabel(f'PC2 ({ev50[1]:.1f}%)', **label_kw)
ax4.set_title('PCA Microbiome Fingerprint\n(50 samples)', **title_kw)
ax4.tick_params(**tick_kw)
ax4.axhline(0, color='#444466', lw=0.5); ax4.axvline(0, color='#444466', lw=0.5)
ax4.legend(handles=handles2, fontsize=8,
           facecolor='#1a1a2e', labelcolor='white', edgecolor='#444466')

# add sample counts
for diag, color in COLORS.items():
    n = (diag50 == diag).sum()
    ax4.scatter([], [], c=color, s=80, label=f'{diag} (n={n})')

# ── Panel 5: Feature importances (50-sample model) ───────────
top15_imp = importances.head(15)
imp_colors = ['#2ecc71' if s in ibd_dep else '#e74c3c' if s in ibd_enr
              else '#3498db' for s in top15_imp.index]
ax5.barh(range(len(top15_imp)), top15_imp.values,
         color=imp_colors, edgecolor='#1a1a2e', height=0.7)
ax5.set_yticks(range(len(top15_imp)))
ax5.set_yticklabels([s.replace('_',' ') for s in top15_imp.index],
                    fontsize=7, color='#cccccc', style='italic')
ax5.invert_yaxis()
ax5.set_xlabel('Feature Importance', **label_kw)
ax5.set_title('Top 15 ML-Predictive Species\n(Random Forest, n=50)', **title_kw)
ax5.tick_params(**tick_kw)
ax5.grid(True, alpha=0.2, axis='x', color='#444466')
h5 = [mpatches.Patch(color='#2ecc71', label='IBD-depleted'),
      mpatches.Patch(color='#e74c3c', label='IBD-enriched'),
      mpatches.Patch(color='#3498db', label='Other')]
ax5.legend(handles=h5, fontsize=7, facecolor='#1a1a2e',
           labelcolor='white', edgecolor='#444466')

# ── Panel 6: ROC + AUC comparison ────────────────────────────
# Recompute ROC for 50-sample LOOCV
y_binary = (diag50 != 'nonIBD').astype(int)
rf = RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=-1)
loo = LeaveOneOut()
y_prob = cross_val_predict(rf, X50_sc, y_binary, cv=loo, method='predict_proba')[:,1]
auc_val = roc_auc_score(y_binary, y_prob)
fpr, tpr, _ = roc_curve(y_binary, y_prob)

ax6.plot(fpr, tpr, color='#3498db', linewidth=2.5,
         label=f'WGS n=50 LOOCV (AUC={auc_val:.3f})')
ax6.plot([0,1],[0,1], 'w--', alpha=0.3, linewidth=1, label='Random (AUC=0.500)')
ax6.axhline(0.894, color='#e74c3c', linestyle=':', linewidth=1.5,
            label='16S n=1233 (AUC=0.894)')
ax6.axhline(0.000, color='#e67e22', linestyle=':', linewidth=1,
            label='WGS n=5 (AUC=0.000)')

ax6.set_xlabel('False Positive Rate', **label_kw)
ax6.set_ylabel('True Positive Rate', **label_kw)
ax6.set_title('ROC Curve — IBD vs nonIBD\nAUC Progression Across Projects', **title_kw)
ax6.tick_params(**tick_kw)
ax6.legend(fontsize=8, facecolor='#1a1a2e', labelcolor='white', edgecolor='#444466')
ax6.grid(True, alpha=0.2, color='#444466')

# ── Main title ────────────────────────────────────────────────
fig.text(0.5, 0.98,
         'Microbiome Network Disruption and ML Classification in IBD — HMP2 WGS Metagenomics',
         ha='center', va='top', color='white', fontsize=13, fontweight='bold')
fig.text(0.5, 0.955,
         '5-sample network & PCA pilot  |  50-sample ML expansion  |  IBD vs nonIBD classification',
         ha='center', va='top', color='#aaaaaa', fontsize=9)

out_path = os.path.join(OUT_FIG, 'final_dashboard.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='#0f0f0f')
print(f"Dashboard saved to {out_path}")
plt.close()
