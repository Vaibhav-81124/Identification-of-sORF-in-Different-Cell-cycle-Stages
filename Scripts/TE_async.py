import pandas as pd

# Base project directory
BASE = "Project Directory"

# ─── FILE PATHS ─────────────────────────────────────────────

# Ribo-seq counts
ribo1_path = f"{BASE}/data/reference/sorf_total_ribo_counts_rep1.txt"
ribo2_path = f"{BASE}/data/reference/sorf_total_ribo_counts_rep2.txt"

# RNA-seq counts
rna1_path = f"{BASE}/results/03_aligned_rnaseq/sorf_rna_counts_rep1_clean.txt"
rna2_path = f"{BASE}/results/03_aligned_rnaseq/sorf_rna_counts_rep2_clean.txt"

# ─── LOAD DATA ─────────────────────────────────────────────

ribo1 = pd.read_csv(
    ribo1_path,
    sep=r"\s+",
    header=None,
    names=["orf_id", "ribo_rep1"]
)

ribo2 = pd.read_csv(
    ribo2_path,
    sep=r"\s+",
    header=None,
    names=["orf_id", "ribo_rep2"]
)

rna1 = pd.read_csv(
    rna1_path,
    sep=r"\s+",
    header=None,
    names=["orf_id", "rna_rep1"]
)

rna2 = pd.read_csv(
    rna2_path,
    sep=r"\s+",
    header=None,
    names=["orf_id", "rna_rep2"]
)

# ─── MERGE TABLES ──────────────────────────────────────────

df = ribo1.merge(ribo2, on="orf_id", how="outer")
df = df.merge(rna1, on="orf_id", how="outer")
df = df.merge(rna2, on="orf_id", how="outer")

# Replace missing values
df = df.fillna(0)

# ─── TRANSLATION EFFICIENCY ─────────────────────────────────

df["TE_rep1"] = df["ribo_rep1"] / (df["rna_rep1"] + 1)
df["TE_rep2"] = df["ribo_rep2"] / (df["rna_rep2"] + 1)

df["TE_mean"] = (df["TE_rep1"] + df["TE_rep2"]) / 2

# ─── SAVE RESULT ───────────────────────────────────────────

output_file = f"{BASE}/results/async_translation_efficiency.csv"

df.to_csv(output_file, index=False)

print("Created:", output_file)
