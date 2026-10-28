import pandas as pd

columns = [
    "qseqid", "sseqid", "pident", "length",
    "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore"
]

files = {
    "A": "./blast_results/set_A.tsv",
    "B": "./blast_results/set_B.tsv",
    "C": "./blast_results/set_C.tsv",
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
