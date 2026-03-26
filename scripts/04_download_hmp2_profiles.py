import pandas as pd
import os
import sys

BASE = '/mnt/e/IBD-metagenomics-project3'
TARGET_CSV = os.path.join(BASE, 'data', 'target_samples.csv')
MERGED_LOCAL = os.path.join(BASE, 'data', 'hmp2_metaphlan_merged.tsv')
OUT_DIR = os.path.join(BASE, 'data', 'hmp2_50samples')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Step 1: Check merged file exists ─────────────────────────
if not os.path.exists(MERGED_LOCAL):
    print("Merged file not found.")
    print(f"Expected at: {MERGED_LOCAL}")
    print("Download taxonomic_profiles_3.tsv.gz from ibdmdb.org → Merged Table")
    print("Then gunzip and rename to hmp2_metaphlan_merged.tsv")
    sys.exit(1)
else:
    print(f"Merged file found: {MERGED_LOCAL}")

# ── Step 2: Load merged file ──────────────────────────────────
print("\nLoading merged file...")
merged = pd.read_csv(MERGED_LOCAL, sep='\t', index_col=0, low_memory=False)
print(f"Merged file shape: {merged.shape[0]} clades x {merged.shape[1]} samples")
print(f"First 3 clade names: {list(merged.index[:3])}")
print(f"First 3 sample IDs (raw): {list(merged.columns[:3])}")

# Strip _profile suffix from column names to match our target IDs
merged.columns = merged.columns.str.replace('_profile$', '', regex=True)
print(f"First 3 sample IDs (after strip): {list(merged.columns[:3])}")

# ── Step 3: Filter to target samples ─────────────────────────
targets = pd.read_csv(TARGET_CSV)
target_ids = targets['External ID'].tolist()

available = [s for s in target_ids if s in merged.columns]
missing   = [s for s in target_ids if s not in merged.columns]

print(f"\nTarget samples found in merged file: {len(available)}/50")
if missing:
    print(f"Missing sample IDs ({len(missing)}): {missing}")

# ── Step 4: Extract species-level rows only ───────────────────
species_rows = merged[
    merged.index.str.contains('s__') &
    ~merged.index.str.contains('t__')
]

# Simplify names: keep only the s__ part
species_rows.index = species_rows.index.map(lambda x: x.split('|s__')[-1])

# Keep only available target samples
matrix = species_rows[available].T  # samples x species
matrix.index.name = 'sample_id'

# Add diagnosis
diag_map = dict(zip(targets['External ID'], targets['diagnosis']))
matrix.insert(0, 'diagnosis', matrix.index.map(diag_map))

print(f"\nFinal matrix: {matrix.shape[0]} samples x {matrix.shape[1]-1} species")
print(f"Diagnosis breakdown:\n{matrix['diagnosis'].value_counts()}")

# Save
os.makedirs(os.path.join(BASE, 'results'), exist_ok=True)
out_path = os.path.join(BASE, 'results', 'species_matrix_50samples.tsv')
matrix.to_csv(out_path, sep='\t')
print(f"\nSaved to {out_path}")
