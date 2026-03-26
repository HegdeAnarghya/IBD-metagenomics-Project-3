import pandas as pd
import os

# Sample ID to diagnosis mapping (from Project 2)
SAMPLES = {
    'data/project2_samples/sample1_taxonomic.tsv': ('MSM5LLHV',  'UC'),
    'data/project2_samples/sample2_taxonomic.tsv': ('HSM7CZ2A',  'nonIBD'),
    'data/project2_samples/sample3_taxonomic.tsv': ('HSM6XRQE',  'UC'),
    'data/project2_samples/sample4_taxonomic.tsv': ('CSM5FZ4C',  'CD'),
    'data/project2_samples/sample5_taxonomic.tsv': ('CSM9X1ZO',  'UC'),
}

BASE = '/mnt/e/IBD-metagenomics-project3'

def parse_metaphlan(filepath):
    """Parse a MetaPhlAn3 profile, return species-level abundances as a dict."""
    abundances = {}
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            clade, tax_id, abundance = parts[0], parts[1], float(parts[2])
            # Keep only species level (s__) but not strain level (t__)
            if '|s__' in clade and '|t__' not in clade:
                species = clade.split('|s__')[-1]
                abundances[species] = abundance
    return abundances

# Build matrix
records = {}
for rel_path, (sample_id, diagnosis) in SAMPLES.items():
    full_path = os.path.join(BASE, rel_path)
    abundances = parse_metaphlan(full_path)
    records[sample_id] = abundances
    print(f"  {sample_id} ({diagnosis}): {len(abundances)} species detected")

# Create DataFrame: rows = samples, columns = species
df = pd.DataFrame(records).T.fillna(0)
df.index.name = 'sample_id'

# Add diagnosis column
diagnoses = {sid: diag for _, (sid, diag) in SAMPLES.items()}
df.insert(0, 'diagnosis', df.index.map(diagnoses))

print(f"\nMatrix shape: {df.shape[0]} samples x {df.shape[1]-1} species")
print(f"\nDiagnosis breakdown:\n{df['diagnosis'].value_counts()}")
print(f"\nTop 10 most abundant species (mean across samples):")
species_cols = df.columns[1:]
print(df[species_cols].mean().sort_values(ascending=False).head(10))

# Save
out_path = os.path.join(BASE, 'results', 'species_abundance_matrix.tsv')
os.makedirs(os.path.dirname(out_path), exist_ok=True)
df.to_csv(out_path, sep='\t')
print(f"\nSaved to {out_path}")
