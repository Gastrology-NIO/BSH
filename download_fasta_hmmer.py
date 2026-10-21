import requests
import os

# def download_fasta_by_hmmer_res(uniprot_acc, output_file):
#   accessions = list(dict.fromkeys(uniprot_acc))

#   with open(output_file, "w") as fasta_file:
#       for acc in accessions:
#           url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"

#           try:
#               response = requests.get(url, timeout=30)
#               response.raise_for_status()

#               fasta = response.text.strip()

#               if fasta.startswith(">"):
#                   fasta_file.write(fasta + "\n")
#                   print(f"Pobrano: {acc}")
#               else:
#                   print(f"Brak sekwencji dla: {acc}")

#           except requests.RequestException as e:
#               print(f"Błąd dla {acc}: {e}")

#   print(f"Zapisano sekwencje do pliku: {output_file}")

def align_new_file(fasta_file, mafft_file):
  command = f"mafft --auto {fasta_file} > {mafft_file}"
  os.system(command)
  
def identify_BSH(aligned_file, new_fasta_file):
  important_positions={
        2: "C",
        18: "R",
        21: "D",
        (81, 82): "N", # często F, Y lub W
        175: "N",
        228: "R",
    }
  reference_sequence=""
  
def mutations_pva():
  # https://www.uniprot.org/uniprotkb/Q6D291/entry
    return {30: ("C", "S"), 
            48: ("R", "H"),
            51: ("D", "E"),
            51: ("D", "N"),
            52: ("W", "F"),
            52: ("W", "I"),
            116: ("W", "N"),
            116: ("W", "Y"),
            244: ("R", "L")
           }

def mutations_BSH():
  return {2:("C", "A"), 
         3: ("T", "A")}
  
def read_uniprots_hmmer(file):
  text=open(file).read()
  return list(set([
      line.split("\t")[1]
      for line in text.strip().splitlines()[1:]
      if line.strip()
  ]))
def phy_to_msa(input_file, output_file):
    with open(input_file, "r") as f:
        lines = f.readlines()

    # Pierwsza linia zawiera liczbę sekwencji i długość wyrównania
    n_sequences, alignment_length = map(int, lines[0].split()[:2])

    sequences = {}

    for line in lines[1:]:
        parts = line.split()
        if len(parts) >= 2:
            name = parts[0]
            seq = "".join(parts[1:])
            sequences[name] = sequences.get(name, "") + seq

    with open(output_file, "w") as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n{seq}\n")


if __name__=="__main__":
  phy_to_msa("phylip.txt", "long_phylip.txt")
  # hmmer_tsv_file = "all_hmmer_results1.tsv"
  # fasta_files="all_hmmer_sequences.csv"
  # align_file="mafft_all_hmmer.aln"
  # # uniprots = read_uniprots_hmmer(hmmer_tsv_file)  
  # # fasta=download_fasta_by_hmmer_res(uniprots, fasta_files)

  # align_new_file(fasta_files, align_file)
  # # identify_BSH(aligned_file, new_fasta_file)
