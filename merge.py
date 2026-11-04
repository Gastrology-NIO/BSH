import pandas as pd
import numpy as np
import re

# ============================================================
# Pliki wejściowe
# ============================================================

files = {
    "PAV": "./blast_results/hmmer_pav_eval_0_01_hsp_perc_90",
    "BSH": "./blast_results/hmmer_bsh_eval_0_01_hsp_perc_90",
    "uniprot": "./blast_results/uniprot_trembl_eval_0_01_hsp_perc_90",
}

pav_file = "./pav.tsv"
bsh_file = "./bsh_hmmsearch_converted.tsv"

columns = [
    "qseqid", "sseqid", "pident", "ppos", "length",
    "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "qcovs", "qseq", "sseq"
]

# ============================================================
# 1. Wczytanie i klasyfikacja białek z HMMER
# ============================================================

pav_hmmer = pd.read_csv(pav_file, sep="\t")
bsh_hmmer = pd.read_csv(bsh_file, sep="\t")

for hmmer_df in [pav_hmmer, bsh_hmmer]:
    hmmer_df.columns = hmmer_df.columns.str.strip()

    hmmer_df["Target Accession"] = (
        hmmer_df["Target Accession"]
        .astype("string")
        .str.strip()
    )

    hmmer_df["Score"] = pd.to_numeric(
        hmmer_df["Score"], errors="coerce"
    )

    hmmer_df["E-value"] = pd.to_numeric(
        hmmer_df["E-value"], errors="coerce"
    )


def best_hmmer_hit(hmmer_df):
    """
    Dla każdego białka wybiera najwyższy HMMER Score.
    Przy identycznym Score wybiera niższe E-value.
    """
    return (
        hmmer_df
        .dropna(subset=["Target Accession", "Score"])
        .sort_values(
            ["Target Accession", "Score", "E-value"],
            ascending=[True, False, True]
        )
        .drop_duplicates("Target Accession")
        .set_index("Target Accession")[["Score", "E-value"]]
    )


pav_scores = best_hmmer_hit(pav_hmmer).rename(columns={
    "Score": "PAV_score",
    "E-value": "PAV_evalue"
})

bsh_scores = best_hmmer_hit(bsh_hmmer).rename(columns={
    "Score": "BSH_score",
    "E-value": "BSH_evalue"
})

# Połączenie wyników dla wszystkich białek
protein_families = pav_scores.join(bsh_scores, how="outer")

pav_present = protein_families["PAV_score"].notna()
bsh_present = protein_families["BSH_score"].notna()
both = pav_present & bsh_present

protein_families["family"] = "unassigned"

# Białko występuje tylko w jednym pliku
protein_families.loc[
    pav_present & ~bsh_present, "family"
] = "PAV"

protein_families.loc[
    bsh_present & ~pav_present, "family"
] = "BSH"

# Białko występuje w obu plikach: porównanie HMMER Score
protein_families.loc[
    both
    & (protein_families["PAV_score"] > protein_families["BSH_score"]),
    "family"
] = "PAV"

protein_families.loc[
    both
    & (protein_families["BSH_score"] > protein_families["PAV_score"]),
    "family"
] = "BSH"

# Identyczny Score: porównanie E-value
score_tie = (
    both
    & (protein_families["PAV_score"] == protein_families["BSH_score"])
)

protein_families.loc[
    score_tie
    & (protein_families["PAV_evalue"] < protein_families["BSH_evalue"]),
    "family"
] = "PAV"

protein_families.loc[
    score_tie
    & (protein_families["BSH_evalue"] < protein_families["PAV_evalue"]),
    "family"
] = "BSH"

# Identyczny Score i E-value
protein_families.loc[
    score_tie
    & (protein_families["PAV_evalue"] == protein_families["BSH_evalue"]),
    "family"
] = "ambiguous"

protein_families = protein_families.reset_index()
protein_families = protein_families.rename(
    columns={"Target Accession": "accession"}
)

# Zapis kontroli klasyfikacji HMMER
protein_families.to_csv(
    "./blast_results/protein_families_from_hmmer.tsv",
    sep="\t",
    index=False
)

print("Klasyfikacja białek na podstawie HMMER:")
print(protein_families["family"].value_counts(dropna=False))

# Słownik: accession -> wyniki HMMER i rodzina
hmmer_lookup = protein_families.set_index("accession")


# ============================================================
# 2. Funkcje do rozpoznawania accession i rodziny referencji
# ============================================================

def extract_accession(identifier):
    """
    Obsługuje m.in.:
      A0A3R7AW25_9FIRM
      tr|A0ABV1BW92|A0ABV1BW92_9FIRM
      sp|P12345|PROTEIN_NAME
    """
    if pd.isna(identifier):
        return None

    identifier = str(identifier).strip()

    # Identyfikatory UniProt w formacie sp|ACC|... lub tr|ACC|...
    match = re.match(r"^(?:sp|tr)\|([^|]+)\|", identifier)
    if match:
        return match.group(1)

    # Najpierw sprawdź cały identyfikator
    if identifier in hmmer_lookup.index:
        return identifier

    # Np. A0A3R7AW25_9FIRM -> A0A3R7AW25
    accession = identifier.split("_")[0]

    if accession in hmmer_lookup.index:
        return accession

    return accession


def get_hmmer_family(identifier):
    """
    Zwraca rodzinę białka według wyników HMMER.
    """
    accession = extract_accession(identifier)

    if accession is None or accession not in hmmer_lookup.index:
        return "unassigned"

    family = hmmer_lookup.loc[accession, "family"]

    if pd.isna(family):
        return "unassigned"

    return family


def get_hmmer_scores(identifier):
    """
    Zwraca wyniki PAV i BSH dla białka referencyjnego.
    """
    accession = extract_accession(identifier)

    if accession is None or accession not in hmmer_lookup.index:
        return np.nan, np.nan

    row = hmmer_lookup.loc[accession]

    return row["PAV_score"], row["BSH_score"]


# ============================================================
# 3. Wczytanie wyników BLASTP
# ============================================================

dfs = []

for dataset, filename in files.items():
    df_part = pd.read_csv(
        filename,
        sep="\t",
        names=columns,
        header=None
    )

    for col in [
        "pident", "ppos", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue",
        "bitscore", "qlen", "slen", "qcovs"
    ]:
        df_part[col] = pd.to_numeric(
            df_part[col], errors="coerce"
        )

    df_part["dataset"] = dataset

    print(f"{dataset}: {df_part.shape[0]} wierszy")

    dfs.append(df_part)

df = pd.concat(dfs, ignore_index=True)


# ============================================================
# 4. Najlepsze trafienie dla każdego query w każdej bazie
# ============================================================

best = (
    df.dropna(subset=["bitscore"])
      .sort_values(
          ["qseqid", "dataset", "bitscore", "evalue"],
          ascending=[True, True, False, True]
      )
      .groupby(["qseqid", "dataset"], sort=False)
      .head(1)
)

comparison = best.pivot(
    index="qseqid",
    columns="dataset",
    values=["sseqid", "pident", "evalue", "bitscore", "qseq", "sseq"]
)

comparison.columns = [
    f"{metric}_{dataset}"
    for metric, dataset in comparison.columns
]

comparison = comparison.reset_index()

# Zapewnienie obecności wszystkich oczekiwanych kolumn
for dataset in files:
    for metric in ["sseqid", "pident", "evalue", "bitscore", "qseq", "sseq"]:
        col = f"{metric}_{dataset}"

        if col not in comparison.columns:
            comparison[col] = np.nan

# Jawna konwersja do wartości liczbowych — zapobiega błędowi isfinite
for col in ["bitscore_PAV", "bitscore_BSH", "bitscore_uniprot"]:
    comparison[col] = pd.to_numeric(
        comparison[col], errors="coerce"
    )


# ============================================================
# 5. Klasyfikacja najlepszego trafienia UniProt
# ============================================================

uniprot_best = best[best["dataset"] == "uniprot"].copy()

uniprot_best["best_uniprot_accession"] = (
    uniprot_best["sseqid"].apply(extract_accession)
)

uniprot_best["best_uniprot_family"] = (
    uniprot_best["best_uniprot_accession"]
    .map(hmmer_lookup["family"])
    .fillna("unassigned")
)

uniprot_family_cols = uniprot_best[
    ["qseqid", "best_uniprot_accession", "best_uniprot_family"]
]

comparison = comparison.merge(
    uniprot_family_cols,
    on="qseqid",
    how="left"
)

comparison["best_uniprot_family"] = (
    comparison["best_uniprot_family"].fillna("unassigned")
)


# ============================================================
# 6. Rodzina białek referencyjnych z trafień PAV i BSH
# ============================================================

for dataset in ["PAV", "BSH"]:
    comparison[f"reference_accession_{dataset}"] = (
        comparison[f"sseqid_{dataset}"].apply(extract_accession)
    )

    comparison[f"reference_family_{dataset}"] = (
        comparison[f"sseqid_{dataset}"].apply(get_hmmer_family)
    )

    scores = comparison[f"sseqid_{dataset}"].apply(get_hmmer_scores)

    comparison[f"reference_PAV_score_{dataset}"] = scores.apply(
        lambda x: x[0]
    )

    comparison[f"reference_BSH_score_{dataset}"] = scores.apply(
        lambda x: x[1]
    )


# ============================================================
# 7. Wybór zestawu z najwyższym bitscore
# ============================================================

score_cols = [
    "bitscore_PAV",
    "bitscore_BSH",
    "bitscore_uniprot"
]

comparison["best_dataset"] = (
    comparison[score_cols]
    .idxmax(axis=1)
    .str.replace("bitscore_", "", regex=False)
)

comparison.loc[
    comparison[score_cols].isna().all(axis=1),
    "best_dataset"
] = pd.NA


# ============================================================
# 8. Rozstrzyganie remisu bitscore PAV vs BSH
# ============================================================

pav_bsh_tie = (
    comparison["bitscore_PAV"].notna()
    & comparison["bitscore_BSH"].notna()
    & np.isclose(
        comparison["bitscore_PAV"].astype(float),
        comparison["bitscore_BSH"].astype(float)
    )
)

# Remis PAV/BSH rozstrzygamy tylko wtedy, gdy UniProt
# nie ma wyższego bitscore.
uniprot_not_better = (
    comparison["bitscore_uniprot"].isna()
    | (
        comparison["bitscore_PAV"]
        >= comparison["bitscore_uniprot"]
    )
)

resolve_tie = pav_bsh_tie & uniprot_not_better


def resolve_reference_tie(row):
    """
    W remisie BLASTP sprawdza klasyfikację referencji z HMMER.

    Jeśli oba trafienia dotyczą tej samej referencji, jej rodzina
    rozstrzyga remis. Jeśli rodziny referencji są zgodne, wybiera
    tę rodzinę. Jeśli są sprzeczne, wynik pozostaje ambiguous.
    """
    bsh_accession = row["reference_accession_BSH"]
    pav_accession = row["reference_accession_PAV"]

    bsh_family = row["reference_family_BSH"]
    pav_family = row["reference_family_PAV"]

    # To samo białko w obu bazach — używamy jego klasyfikacji HMMER
    if (
        pd.notna(bsh_accession)
        and pd.notna(pav_accession)
        and bsh_accession == pav_accession
    ):
        if bsh_family in ["PAV", "BSH"]:
            return bsh_family
        return "ambiguous"

    # Różne referencje, ale zgodna klasyfikacja HMMER
    if bsh_family == pav_family and bsh_family in ["PAV", "BSH"]:
        return bsh_family

    # Jeśli tylko jedna referencja ma wiarygodną klasyfikację,
    # można użyć jej rodziny.
    if bsh_family in ["PAV", "BSH"] and pav_family not in ["PAV", "BSH"]:
        return bsh_family

    if pav_family in ["PAV", "BSH"] and bsh_family not in ["PAV", "BSH"]:
        return pav_family

    # Konflikt klasyfikacji albo brak wystarczających danych
    return "ambiguous"


comparison.loc[resolve_tie, "best_dataset"] = (
    comparison.loc[resolve_tie]
    .apply(resolve_reference_tie, axis=1)
)

# ============================================================
# 9. Końcowa klasyfikacja i ID przypisanego białka
# ============================================================

comparison["assigned_family"] = comparison["best_dataset"]
comparison["assigned_protein_id"] = pd.NA
for col in ["sseqid_BSH", "sseqid_PAV", "best_uniprot_accession"]:
    if col in comparison.columns:
        comparison[col] = comparison[col].str.split("_").str[0]
# Najlepsze trafienie PAV
is_pav = comparison["best_dataset"].eq("PAV")
comparison.loc[is_pav, "assigned_family"] = "PAV"
comparison.loc[is_pav, "assigned_protein_id"] = (
    comparison.loc[is_pav, "sseqid_PAV"]
)

# Najlepsze trafienie BSH
is_bsh = comparison["best_dataset"].eq("BSH")
comparison.loc[is_bsh, "assigned_family"] = "BSH"
comparison.loc[is_bsh, "assigned_protein_id"] = (
    comparison.loc[is_bsh, "sseqid_BSH"]
)

# Najlepsze trafienie UniProt — wykorzystujemy rodzinę przypisaną
# do białka UniProt na podstawie pav.tsv i bsh.tsv
is_uniprot = comparison["best_dataset"].eq("uniprot")
comparison.loc[is_uniprot, "assigned_family"] = (
    comparison.loc[is_uniprot, "best_uniprot_family"]
)
comparison.loc[is_uniprot, "assigned_protein_id"] = (
    comparison.loc[is_uniprot, "best_uniprot_accession"]
)


# ============================================================
# 10. Wybór sekwencji z właściwego najlepszego trafienia BLAST
# ============================================================

# Zachowujemy qseq i sseq z oryginalnego wiersza BLAST
alignment_columns = best[
    ["qseqid", "dataset", "qseq", "sseq"]
].copy()

# Jednoznaczne nazwy kolumn
alignment_columns = alignment_columns.rename(
    columns={"dataset": "alignment_dataset"}
)

# Usuwamy ewentualne stare kolumny sekwencji
comparison = comparison.drop(
    columns=[
        col for col in ["qseq", "sseq"]
        if col in comparison.columns
    ]
)

# Dołączamy przyrównanie z wybranego zbioru
comparison = comparison.merge(
    alignment_columns,
    left_on=["qseqid", "best_dataset"],
    right_on=["qseqid", "alignment_dataset"],
    how="left",
    validate="one_to_one"
)

comparison = comparison.drop(columns=["alignment_dataset"])



# Zapis wyników
output_file = "./blast_results/comparison_ABC.tsv"

comparison.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("\nKońcowa klasyfikacja:")
print(comparison["assigned_family"].value_counts(dropna=False))

print("\nLiczba przypisanych ID:")
print(comparison.groupby("assigned_family")["assigned_protein_id"].count())

print(f"\nZapisano wyniki do: {output_file}")