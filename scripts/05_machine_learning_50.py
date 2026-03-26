import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut, StratifiedKFold, cross_val_predict
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve, classification_report
from sklearn.decomposition import PCA
import os

BASE = '/mnt/e/IBD-metagenomics-project3'
MATRIX = os.path.join(BASE, 'results', 'species_matrix_50samples.tsv')
OUT_FIG = os.path.join(BASE, 'figures')
OUT_RES = os.path.join(BASE, 'results')
os.makedirs(OUT_FIG, exist_ok=True)
os.makedirs(OUT_RES, exist_ok=True)

# ── Load ──────────────────────────────────────────────────────
df = pd.read_csv(MATRIX, sep='\t', index_col='sample_id')
diagnoses = df['diagnosis']
X = df.drop(columns=['diagnosis'])

# ── Preprocessing ─────────────────────────────────────────────
# Remove species present in fewer than 5 samples (reduce noise)
X = X.loc[:, (X > 0).sum() >= 5]
print(f"Species after prevalence filter (≥5 samples): {X.shape[1]}")

# Log transform + standardize
X_log = np.log1p(X)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_log)
X_scaled = pd.DataFrame(X_scaled, index=X.index, columns=X.columns)

# ── Task 1: IBD vs nonIBD (binary) ───────────────────────────
print("\n── Task 1: IBD vs nonIBD (binary) ──")
y_binary = (diagnoses != 'nonIBD').astype(int)
print(f"IBD: {y_binary.sum()}  nonIBD: {(y_binary==0).sum()}")

rf = RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=-1)

# LOOCV
loo = LeaveOneOut()
y_prob_loo = cross_val_predict(rf, X_scaled, y_binary, cv=loo, method='predict_proba')[:, 1]
auc_loo = roc_auc_score(y_binary, y_prob_loo)
print(f"LOOCV AUC (IBD vs nonIBD): {auc_loo:.3f}")

# 5-fold CV (more stable estimate)
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
y_prob_5fold = cross_val_predict(rf, X_scaled, y_binary, cv=skf, method='predict_proba')[:, 1]
auc_5fold = roc_auc_score(y_binary, y_prob_5fold)
print(f"5-fold CV AUC (IBD vs nonIBD): {auc_5fold:.3f}")

# ── Task 2: CD vs UC vs nonIBD (multiclass) ──────────────────
print("\n── Task 2: CD vs UC vs nonIBD (3-class) ──")
le = LabelEncoder()
y_multi = le.fit_transform(diagnoses)
print(f"Classes: {le.classes_}")
y_prob_multi = cross_val_predict(rf, X_scaled, y_multi, cv=skf, method='predict_proba')
auc_multi = roc_auc_score(pd.get_dummies(diagnoses), y_prob_multi, multi_class='ovr', average='macro')
print(f"5-fold CV AUC (macro, 3-class): {auc_multi:.3f}")

report = classification_report(y_multi,
    cross_val_predict(rf, X_scaled, y_multi, cv=skf),
    target_names=le.classes_)
print(f"\nClassification Report:\n{report}")

# ── Feature importances ───────────────────────────────────────
rf.fit(X_scaled, y_binary)
importances = pd.Series(rf.feature_importances_, index=X_scaled.columns)
importances = importances.sort_values(ascending=False)
importances.to_csv(os.path.join(OUT_RES, 'feature_importances_50.tsv'), sep='\t', header=['importance'])
print(f"\nTop 15 predictive species:")
print(importances.head(15).to_string())

# ── Summary table ─────────────────────────────────────────────
summary = pd.DataFrame({
    'Model': ['Project 1 (16S, n=1233)', 'Project 2 (WGS, n=5)', 'Project 3 (WGS, n=50) LOOCV', 'Project 3 (WGS, n=50) 5-fold'],
    'AUC': [0.894, 0.000, auc_loo, auc_5fold],
    'n_samples': [1233, 5, 50, 50],
    'Method': ['Random Forest', 'Random Forest LOOCV', 'Random Forest LOOCV', 'Random Forest 5-fold CV']
})
summary.to_csv(os.path.join(OUT_RES, 'ml_comparison.tsv'), sep='\t', index=False)
print(f"\n── AUC Comparison ──")
print(summary.to_string(index=False))

# ── Plots ─────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle('ML Classification — IBD vs nonIBD (n=50 WGS samples)', fontsize=14, fontweight='bold')

# Panel 1: ROC curves
ax = axes[0]
fpr_loo, tpr_loo, _ = roc_curve(y_binary, y_prob_loo)
fpr_5f, tpr_5f, _   = roc_curve(y_binary, y_prob_5fold)
ax.plot(fpr_loo, tpr_loo, 'b-', linewidth=2, label=f'WGS LOOCV (AUC={auc_loo:.3f})')
ax.plot(fpr_5f,  tpr_5f,  'g-', linewidth=2, label=f'WGS 5-fold (AUC={auc_5fold:.3f})')
ax.plot([0,1],[0,1], 'k--', alpha=0.5, label='Random (AUC=0.500)')
ax.axhline(y=0.894, color='red', linestyle=':', alpha=0.7, linewidth=1.5, label='Project 1 16S (AUC=0.894)')
ax.set_xlabel('False Positive Rate', fontsize=11)
ax.set_ylabel('True Positive Rate', fontsize=11)
ax.set_title('ROC Curves', fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 2: Feature importances (top 20)
ax = axes[1]
top20 = importances.head(20)
colors = ['#e74c3c' if s in {'Faecalibacterium_prausnitzii','Akkermansia_muciniphila','Roseburia_intestinalis'}
          else '#3498db' for s in top20.index]
ax.barh(range(len(top20)), top20.values, color=colors, edgecolor='white')
ax.set_yticks(range(len(top20)))
ax.set_yticklabels([s.replace('_',' ') for s in top20.index], fontsize=8, style='italic')
ax.invert_yaxis()
ax.set_xlabel('Feature Importance', fontsize=11)
ax.set_title('Top 20 Predictive Species', fontsize=11)
ax.grid(True, alpha=0.3, axis='x')

# Panel 3: AUC comparison bar chart
ax = axes[2]
colors_bar = ['#95a5a6', '#e74c3c', '#2ecc71', '#27ae60']
bars = ax.bar(summary['Model'], summary['AUC'], color=colors_bar, edgecolor='black', linewidth=0.8)
ax.set_ylabel('AUC', fontsize=11)
ax.set_title('AUC Comparison Across Projects', fontsize=11)
ax.set_ylim(0, 1.05)
ax.axhline(0.5, color='black', linestyle='--', alpha=0.4, label='Random baseline')
for bar, auc in zip(bars, summary['AUC']):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
            f'{auc:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
ax.set_xticklabels(summary['Model'], rotation=15, ha='right', fontsize=8)
ax.grid(True, alpha=0.3, axis='y')
ax.legend(fontsize=9)

plt.tight_layout()
out_path = os.path.join(OUT_FIG, 'ml_classification_50.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"\nFigure saved to {out_path}")
plt.close()
