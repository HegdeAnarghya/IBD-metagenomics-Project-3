# Data Download Instructions

The raw data for this project comes from the NIH Human Microbiome Project 2 (HMP2), also known as the Integrative Human Microbiome Project (iHMP). Due to file size and data use policies, raw data files are not included in this repository.

---

## Data Source

**Portal:** NIH Human Microbiome Project 2 — Inflammatory Bowel Disease Multi-omics Database (ibdmdb.org)
**Study:** HMP2 / iHMP
**Data type:** WGS metagenomics — pre-computed MetaPhlAn 3.0 taxonomic profiles

---

## File Used in This Project

| File | Description | Size |
|------|-------------|------|
| `taxonomic_profiles_3.tsv.gz` | MetaPhlAn 3.0 merged taxonomic profiles — all HMP2 WGS samples | ~5 MB |
| `hmp2_metadata.csv` | Sample metadata with diagnosis, visit number, participant ID | < 1 MB |

---

## Download Instructions

### Step 1: Accept the Data Use Agreement

Go to the HMP2 portal (search "ibdmdb" or "iHMP portal") and accept the data use terms before downloading any files.

### Step 2: Navigate to the Merged Taxonomic Table

On the portal:

```
Results page → Available Studies → HMP2 / Metagenomics(MGX) / 2018.18 → products → Merged Table
```

Download `taxonomic_profiles_3.tsv.gz`.

### Step 3: Download the Metadata

On the same Results page, click:

```
Download HMP2 Metadata
```

Save as `hmp2_metadata.csv`.

### Step 4: Place and decompress files

```bash
# Move downloaded files (replace YOUR_USERNAME with your Windows username)
cp /mnt/c/Users/YOUR_USERNAME/Downloads/taxonomic_profiles_3.tsv.gz \
   /mnt/e/IBD-metagenomics-project3/data/

cp /mnt/c/Users/YOUR_USERNAME/Downloads/hmp2_metadata.csv \
   /mnt/e/IBD-metagenomics-project3/data/

# Decompress and rename
gunzip /mnt/e/IBD-metagenomics-project3/data/taxonomic_profiles_3.tsv.gz

mv /mnt/e/IBD-metagenomics-project3/data/taxonomic_profiles_3.tsv \
   /mnt/e/IBD-metagenomics-project3/data/hmp2_metaphlan_merged.tsv
```

### Step 5: Verify

```bash
wc -l /mnt/e/IBD-metagenomics-project3/data/hmp2_metaphlan_merged.tsv
# Expected: ~933 lines (932 clades + header)

head -1 /mnt/e/IBD-metagenomics-project3/data/hmp2_metaphlan_merged.tsv | cut -c1-100
# Expected: Feature\Sample   CSM5FZ3N_P_profile   CSM5FZ3R_P_profile ...
```

---

## Sample Selection

This project uses 50 samples selected from the full HMP2 WGS cohort (130 unique subjects). Selection criteria:

- One sample per subject (baseline visit — lowest `visit_num`)
- Balanced across diagnosis groups: 20 CD, 17 UC, 13 nonIBD
- Sample IDs listed in `data/target_samples.csv`

The selection script (`04_download_hmp2_profiles.py`) filters the merged table automatically to these 50 samples.

---

## Option: Download Raw FASTQ Files

If you want to rerun from raw reads (not required for this project — pre-computed profiles are sufficient):

### Requirements
- 50–200 GB storage per sample
- 16+ GB RAM for HUMAnN
- Several hours of compute per sample

### Download via SRA Toolkit

```bash
conda install -c bioconda sra-tools

# SRA accessions are listed in hmp2_metadata.csv
prefetch SRR_ACCESSION
fastq-dump --split-files SRR_ACCESSION
```

### Full pipeline from FASTQ

```bash
# 1. Quality control
trimmomatic PE -threads 8 \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    sample_R1_trimmed.fastq.gz sample_R1_unpaired.fastq.gz \
    sample_R2_trimmed.fastq.gz sample_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36

# 2. Host read removal
bowtie2 -x human_genome_index \
    -1 sample_R1_trimmed.fastq.gz \
    -2 sample_R2_trimmed.fastq.gz \
    --un-conc-gz sample_host_removed.fastq.gz \
    -S /dev/null

# 3. Taxonomic profiling
metaphlan sample_host_removed_1.fastq.gz,sample_host_removed_2.fastq.gz \
    --input_type fastq \
    --nproc 8 \
    -o sample_metaphlan_profile.txt

# 4. Merge profiles
merge_metaphlan_tables.py *_metaphlan_profile.txt > hmp2_metaphlan_merged.tsv
```

---

## File Size Reference

| File | Approximate Size |
|------|-----------------|
| `taxonomic_profiles_3.tsv.gz` (download) | ~5 MB |
| `hmp2_metaphlan_merged.tsv` (decompressed) | ~20 MB |
| `hmp2_metadata.csv` | < 1 MB |
| Raw FASTQ per sample (paired) | 2–10 GB |

---

## Citation

Lloyd-Price, J., Arze, C., Ananthakrishnan, A.N. et al. **Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases.** *Nature* 569, 655–662 (2019). DOI: 10.1038/s41586-019-1237-9
