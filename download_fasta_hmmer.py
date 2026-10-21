def download_fasta_by_hmmer_res(uniprot_acc, output_file):
  accessions = list(dict.fromkeys(uniprot_acc))

  with open(output_file, "w") as fasta_file:
      for acc in accessions:
          url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"

          try:
              response = requests.get(url, timeout=30)
              response.raise_for_status()

              fasta = response.text.strip()

              if fasta.startswith(">"):
                  fasta_file.write(fasta + "\n")
                  print(f"Pobrano: {acc}")
              else:
                  print(f"Brak sekwencji dla: {acc}")

          except requests.RequestException as e:
              print(f"Błąd dla {acc}: {e}")

  print(f"Zapisano sekwencje do pliku: {output_file}")

def align_new_file(file, new_file):
  pass
  
def identify_BSH(aligned_file, new_fasta_file):
  pass

def read_uniprots_hmmer(file):
    return list(set([
      line.split("\t")[1]
      for line in text.strip().splitlines()[1:]
      if line.strip()
  ]))
  
if __name__=="__main__":
  hmmer_tsv_file = "all_hmmer_results1.tsv"
  fasta_files="all_hmmer_sequences.csv"
  uniprots = read_uniprots_hmmer(hmmer_tsv_file)  
  fasta=download_fasta_by_hmmer_res(uniprots, fasta_files)

  # align_new_file(file, new_file)
  # identify_BSH(aligned_file, new_fasta_file)
