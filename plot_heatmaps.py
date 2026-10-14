import numpy as np
import matplotlib

matplotlib.use("TkAgg")

import matplotlib.pyplot as plt

cancer_type="rak_zoladka"
file = f"./headers/corr_{cancer_type}_merged_paires_rhea_down_up_headers.csv"


data = np.genfromtxt(
    file,
    delimiter=";",
    dtype=None,
    names=True,
    encoding="utf-8"
)
data_significant = data[data["test_BH_bool"]]
print(len(data_significant))
unique, counts = np.unique(
    data_significant["metabolit"],
    return_counts=True
)

for metabolit, liczba in zip(unique, counts):
    print(metabolit, liczba)

metabolity = np.unique(data_significant["metabolit"])
grupy = np.unique(data_significant["grupa_sekwencji"])

# teraz: wiersze = grupy, kolumny = metabolity
heatmap = np.full((len(grupy), len(metabolity)), np.nan)

for i, grupa in enumerate(grupy):
    for j, metabolit in enumerate(metabolity):

        mask = (
            (data_significant["metabolit"] == metabolit) &
            (data_significant["grupa_sekwencji"] == grupa)
        )

        values = data_significant["pearson_corr"][mask]

        if len(values) > 0:
            heatmap[i, j] = values[0]

plt.figure(figsize=(12, 10))
plt.imshow(
    heatmap,
    aspect="auto",
    cmap="RdBu_r",
    vmin=-1,
    vmax=1
)

plt.xticks(
    np.arange(len(metabolity)),
    metabolity,
    rotation=90
)
plt.colorbar(label="Pearson correlation")

plt.yticks(
    np.arange(len(grupy)),
    grupy
)

plt.xlabel("Metabolit")
plt.ylabel("Grupa sekwencji")
plt.tight_layout()
plt.savefig(
    f"heatmap_rhea_pearson_corr__all_{cancer_type}.png",
    dpi=300,
    bbox_inches="tight"
)


# usunięcie nan w pident_rhea
def plot_by_group(data_significant, name="rhea", pident_name="pident_rhea",
                  database_name="Rhea_Identyfikator_rhea", protein_name="Nazwa_białka_rhea",
                  organism_name="Nazwa_organizmy_rhea"):
    mask = ~np.isnan(data_significant[pident_name])
    data_rhea = data_significant[mask]
    print(data_rhea.dtype.names)
    rhea_id = np.char.replace(
        data_rhea[database_name].astype(str),
        "./blast_results/",
        ""
    )
    def unique_value(x):
        ids = str(x).split(",")
        ids = sorted(set(ids))
        return ",".join(ids)

    rhea_id = np.array([
        unique_value(x) for x in rhea_id
    ])
    Nazwa_bialka_rhea=np.array([
        unique_value(x) for x in data_rhea[protein_name].astype(str)
    ])
    Nazwa_organizmy_rhea=np.array([
        unique_value(x) for x in data_rhea[organism_name].astype(str)
    ])
    # stworzenie nowej nazwy klastra
    cluster_name = (
        rhea_id
        + "_"
        + Nazwa_bialka_rhea
        + "_"
        + Nazwa_organizmy_rhea
    )

    # unikalne wartości
    metabolity = np.unique(data_rhea["metabolit"])
    klastry = np.unique(cluster_name)

    # macierz: y = metabolity, x = klastry
    heatmap = np.full((len(metabolity), len(klastry)), np.nan)

    for i, metabolit in enumerate(metabolity):
        for j, klaster in enumerate(klastry):

            mask = (
                (data_rhea["metabolit"] == metabolit) &
                (cluster_name == klaster)
            )

            values = data_rhea["pearson_corr"][mask]

            if len(values) > 0:
                heatmap[i, j] = values[0]


    # rysowanie
    plt.figure(figsize=(15, 10))
    print(data_rhea.dtype.names)
    plt.imshow(
        heatmap,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-1,
        vmax=1
    )

    plt.colorbar(label="Pearson correlation")

    plt.xticks(
        np.arange(len(klastry)),
        klastry,
        rotation=90,
        fontsize=8
    )

    plt.yticks(
        np.arange(len(metabolity)),
        metabolity
    )

    plt.xlabel(f"{name} protein cluster")
    plt.ylabel("Metabolit")

    plt.tight_layout()

    plt.savefig(
        f"heatmap_{name}_pearson_corr_{cancer_type}.png",
        dpi=300,
        bbox_inches="tight"
    )
