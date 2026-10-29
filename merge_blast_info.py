import pandas as pd
# blastp \                                                                                                            >   -db ./databases/trembl/uniprot_trembl.fasta \
# >   -out ./blast_results/uniprot_trembl_eval_0_01_hsp_perc_90 \
# >   -query ./clusters_representatives_for_corr_merged_pairs.fasta \
# >   -num_threads 10 \
# >   -outfmt "6 qseqid sseqid pident ppos length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qseq sseq" \
# >   -evalue 0.01 \
# >   -qcov_hsp_perc 90

# blastp \
# >   -db ./databases/hmmer_bsh.fa \
# >   -out ./blast_results/hmmer_bsh_eval_0_01_hsp_perc_90 \
# >   -query ./clusters_representatives_for_corr_merged_pairs.fasta \
# >   -num_threads 10 \
# >   -outfmt "6 qseqid sseqid pident ppos length mismatch gapopen qstart qend sstart send evalue bitsco
# re qlen slen qcovs qseq sseq" \
# >   -evalue 0.01 \
# >   -qcov_hsp_perc 90

# blastp \
# >   -db ./databases/hmmer_pav.fa \
# >   -out ./blast_results/hmmer_pav_eval_0_01_hsp_perc_90 \
# >   -query ./clusters_representatives_for_corr_merged_pairs.fasta \
# >   -num_threads 10 \
# >   -outfmt "6 qseqid sseqid pident ppos length mismatch gapopen qstart qend sstart send evalue bitsco
# re qlen slen qcovs qseq sseq" \
# >   -evalue 0.01 \
# >   -qcov_hsp_perc 90

# qseqid – Query sequence ID.
# sseqid – Subject protein ID.
# pident – Percentage of identical amino acids.
# ppos – Percentage of positive-scoring amino acid substitutions.
# length – Alignment length.
# evalue – Statistical significance of the alignment.
# bitscore – Bit score representing the strength of the alignment.
# qlen – Length of the query sequence.
# slen – Length of the subject protein.
# qcovs – Percentage of the query sequence covered by the alignment.



columns = [
    "qseqid", "sseqid", "pident", "ppos", "length",
    "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore"
]

files = {
    "PAV": "./blast_results/hmmer_pav_eval_0_01_hsp_perc_90",
    "BSH": "./blast_results/hmmer_bsh_eval_0_01_hsp_perc_90",
    "uniprot": "./blast_results/uniprot_trembl_eval_0_01_hsp_perc_90",
}

# Wczytanie wyników
dfs = []

for dataset, filename in files.items():
    df = pd.read_csv(filename, sep="\t", names=columns)
    df["dataset"] = dataset
    dfs.append(df)

# Połączenie wszystkich zestawów
df = pd.concat(dfs, ignore_index=True)

# Najlepsze trafienie w każdym zestawie dla każdego query
best = (
    df.sort_values("bitscore", ascending=False)
      .groupby(["qseqid", "dataset"])
      .head(1)
)

# Porównanie zestawów
comparison = best.pivot(
    index="qseqid",
    columns="dataset",
    values=["sseqid", "pident", "evalue", "bitscore"]
)

comparison.columns = [
    f"{metric}_{dataset}"
    for metric, dataset in comparison.columns
]

comparison = comparison.reset_index()

# Wybór zestawu z najwyższym bitscore
score_cols = [f"bitscore_{name}" for name in files]

comparison["best_dataset"] = (
    comparison[score_cols]
    .idxmax(axis=1)
    .str.replace("bitscore_", "", regex=False)
)

comparison.to_csv(
    "./blast_results/comparison_ABC.tsv",
    sep="\t",
    index=False
)
