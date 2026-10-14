import pandas as pd
from Bio import SeqIO

def read_tab_hmmer(file):
  df = pd.read_csv(file)
  return df
a
def download_fasta(df):
  fasta_sequences={}
  with open("proteins_tmp.fasta", "w") as f:
    for acc in df["Target Accession"]:
      url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
      response = requests.get(url)
      if response.ok:
        fasta=response.text
        fasta_sequences[acc]={"fasta":fasta, "header":fasta.split("\n")[0], "seq":fasta.split("\n")[1:]}
        f.write(fasta)
  return fasta_sequences

def significant_bsh():
  return {2:"C", 18:"R", 21:"D", 82:"N", 228:"R"}

def select_proteins_with_positions(protein_positions, proteins):
  pass
  
def recalculate_positions_based_on_significant(significant_aas, canonical_acc, proteins):
  new_positions={}
  #todo
  return new_positions

def save_fasta_from_mafft(proteins):
  pass

def calc_mafft(input_file, output_file):

  with open(output_file, "w") as out:
    subprocess.run(
        ["mafft", "--auto", input_file],
        stdout=out,
        check=True
    )


def read_mafft(file):
  #todo
  alignment = {}
  for record in SeqIO.parse(file, "fasta"):
      acc = record.id
      seq = str(record.seq)
      alignment[acc] = seq
  return alignment

def run_hmmer():
  #todo
  pass

def save_fasta(mafft_results):
  pass
  
#assign proteins and count reads, compare them and calc corr


def run(proteins_files):
  #to jest do szczegółowej analizy BSH. 
  canonical_sequence_uniprotACC=""
  # weryfikacja hmmer czy wynik jest ok
  # 1) wczytanie wyników
  hmmer_proteins_uniprotacc=read_tab_hmmer(proteins_files)
  # 2) pobranie fasta
  proteins = download_fasta(hmmer_proteins_uniprotacc)
  # 3) mafft
  alignment_file = "alignment_before_clean.fasta"
  input_file = "proteins_tmp.fasta"
  calc_mafft(input_file, alignment_file)
  mafft_result = read_mafft(alignment_file)
  # 4) lokalizacja białek po significant_bsh
  significant_aas = significant_bsh()
  new_positions = recalculate_positions_based_on_significant(significant_aas, canonical_sequence_uniprotACC, mafft_result)
  # 5) ew. lokalizacja białek z opcjonalnymi aminokwasami
  new_fasta = select_proteins_with_positions(new_positions, mafft_result)
  # 6) policzenie maffta 
  input_file = "new_protein_file_tmp.fasta"
  save_fasta_from_mafft(new_fasta, new_protein_file)
  alignment_file = "alignment_after_clean.fasta"
  calc_mafft(input_file, alignment_file)
  mafft_result = read_mafft(alignment_file)
  # 7) zapisanie do fasta
  
