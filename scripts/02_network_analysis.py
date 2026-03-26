import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import spearmanr
import os

BASE = '/mnt/e/IBD-metagenomics-project3'
MATRIX = os.path.join(BASE, 'results', 'species_abundance_matrix.tsv')
OUT_FIG = os.path.join(BASE, 'figures')
OUT_RES = os.path.join(BASE, 'results')
os.makedirs(OUT_FIG, exist_ok=True)
os.makedirs(OUT_RES, exist_ok=True)

# ── Load matrix ───────────────────────────────────────────────
df = pd.read_csv(MATRIX, sep='\t', index_col='sample_id')
diagnoses = df['diagnosis']
X = df.drop(columns=['diagnosis'])

# Keep only species present in at least 2 samples (reduces noise)
X = X.loc[:, (X > 0).sum() >= 2]
print(f"Species retained (present in ≥2 samples): {X.shape[1]}")

# ── Compute Spearman correlations between all species pairs ───
species = X.columns.tolist()
n = len(species)
corr_matrix = pd.DataFrame(np.zeros((n, n)), index=species, columns=species)

for i in range(n):
    for j in range(i, n):
        if i == j:
            corr_matrix.iloc[i, j] = 1.0
        else:
            r, _ = spearmanr(X.iloc[:, i], X.iloc[:, j])
            corr_matrix.iloc[i, j] = r
            corr_matrix.iloc[j, i] = r

corr_matrix.to_csv(os.path.join(OUT_RES, 'species_correlation_matrix.tsv'), sep='\t')
print(f"Correlation matrix saved ({n}x{n})")

# ── Build network: edge if |r| >= 0.9 ────────────────────────
# With n=5 samples, use high threshold to keep only strong signals
THRESHOLD = 0.9
G = nx.Graph()
G.add_nodes_from(species)

pos_edges, neg_edges = 0, 0
for i in range(n):
    for j in range(i+1, n):
        r = corr_matrix.iloc[i, j]
        if abs(r) >= THRESHOLD:
            G.add_edge(species[i], species[j],
                       weight=r,
                       sign='positive' if r > 0 else 'negative')
            if r > 0:
                pos_edges += 1
            else:
                neg_edges += 1

# Remove isolated nodes (no edges)
isolated = list(nx.isolates(G))
G.remove_nodes_from(isolated)

print(f"\nNetwork (|r| ≥ {THRESHOLD}):")
print(f"  Nodes (species): {G.number_of_nodes()}")
print(f"  Edges total:     {G.number_of_edges()}")
print(f"  Positive edges:  {pos_edges}")
print(f"  Negative edges:  {neg_edges}")

# ── Network metrics ───────────────────────────────────────────
degree      = dict(G.degree())
betweenness = nx.betweenness_centrality(G)
closeness   = nx.closeness_centrality(G)

metrics = pd.DataFrame({
    'degree':      degree,
    'betweenness': betweenness,
    'closeness':   closeness,
}).sort_values('degree', ascending=False)

metrics.index.name = 'species'
metrics.to_csv(os.path.join(OUT_RES, 'network_metrics.tsv'), sep='\t')

print(f"\nTop 10 Hub Species (by degree — most connected):")
print(metrics.head(10).to_string())

print(f"\nTop 10 Keystone Species (by betweenness — bridges communities):")
print(metrics.sort_values('betweenness', ascending=False).head(10).to_string())

# ── Visualize ─────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(18, 8))
fig.suptitle('Gut Microbiome Co-occurrence Network (5 HMP2 Samples)', fontsize=14, fontweight='bold')

# Known IBD-relevant species for coloring
ibd_enriched  = {'Bacteroides_vulgatus', 'Bacteroides_uniformis',
                 'Prevotella_stercorea', 'Bacteroides_stercoris'}
ibd_depleted  = {'Faecalibacterium_prausnitzii', 'Akkermansia_muciniphila',
                 'Roseburia_intestinalis', 'Ruminococcus_bromii'}

def node_color(sp):
    if sp in ibd_enriched:  return '#e74c3c'   # red = enriched in IBD
    if sp in ibd_depleted:  return '#2ecc71'   # green = depleted in IBD
    return '#95a5a6'                            # grey = neutral

# Panel 1: Full network
ax = axes[0]
pos = nx.spring_layout(G, seed=42, k=0.5)
node_colors = [node_color(n) for n in G.nodes()]
node_sizes  = [300 + degree[n] * 80 for n in G.nodes()]

pos_edge_list = [(u, v) for u, v, d in G.edges(data=True) if d['sign'] == 'positive']
neg_edge_list = [(u, v) for u, v, d in G.edges(data=True) if d['sign'] == 'negative']

nx.draw_networkx_edges(G, pos, edgelist=pos_edge_list, edge_color='#3498db',
                       alpha=0.5, width=1.2, ax=ax)
nx.draw_networkx_edges(G, pos, edgelist=neg_edge_list, edge_color='#e67e22',
                       alpha=0.5, width=1.2, style='dashed', ax=ax)
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes,
                       alpha=0.9, ax=ax)

# Label only top hubs
top_hubs = set(metrics.head(8).index)
labels = {n: n.replace('_', '\n') for n in G.nodes() if n in top_hubs}
nx.draw_networkx_labels(G, pos, labels, font_size=6, ax=ax)

ax.set_title(f'Co-occurrence Network\n({G.number_of_nodes()} species, {G.number_of_edges()} edges, |r|≥{THRESHOLD})', fontsize=11)
ax.axis('off')

legend_handles = [
    mpatches.Patch(color='#e74c3c', label='IBD-enriched species'),
    mpatches.Patch(color='#2ecc71', label='IBD-depleted species'),
    mpatches.Patch(color='#95a5a6', label='Other species'),
    plt.Line2D([0],[0], color='#3498db', label='Positive co-occurrence'),
    plt.Line2D([0],[0], color='#e67e22', linestyle='dashed', label='Negative co-occurrence'),
]
ax.legend(handles=legend_handles, loc='lower left', fontsize=8)

# Panel 2: Hub species bar chart
ax2 = axes[1]
top15 = metrics.head(15)
colors = [node_color(sp) for sp in top15.index]
bars = ax2.barh(range(len(top15)), top15['degree'], color=colors, edgecolor='white')
ax2.set_yticks(range(len(top15)))
ax2.set_yticklabels([s.replace('_', ' ') for s in top15.index], fontsize=9, style='italic')
ax2.invert_yaxis()
ax2.set_xlabel('Degree (number of co-occurrence partners)', fontsize=10)
ax2.set_title('Top 15 Hub Species\n(node size = connectivity)', fontsize=11)
ax2.axvline(x=top15['degree'].mean(), color='black', linestyle='--', alpha=0.5, label='Mean degree')
ax2.legend(fontsize=9)

plt.tight_layout()
out_path = os.path.join(OUT_FIG, 'network_analysis.png')
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"\nFigure saved to {out_path}")
plt.close()

# ── Network summary stats ─────────────────────────────────────
print(f"\n── Network Summary ──")
if nx.is_connected(G):
    print(f"  Network is connected")
    print(f"  Average shortest path: {nx.average_shortest_path_length(G):.3f}")
else:
    components = list(nx.connected_components(G))
    print(f"  Network has {len(components)} connected components")
    print(f"  Largest component: {max(len(c) for c in components)} species")

print(f"  Density: {nx.density(G):.4f}")
print(f"  Average clustering coefficient: {nx.average_clustering(G):.4f}")
