import os

import click
import numpy as np
from scipy.stats import pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests
import numpy as np
import requests

def read_stats_file(count_genes_file):
    file_names = list()
    stats = {}
    with open(count_genes_file) as f:
        for l in f:
            if l.strip():
                if "genes" in l:
                    file_names = l.strip().split(";")[1:]
                    stats = {i: {} for i in file_names}
                else:
                    line = l.strip().split(";")
                    gene = line[0].split("_")[1]
                    genes_no = line[1:]
                    for e, file in enumerate(file_names):
                        stats[file][gene] = float(genes_no[e])
    return stats

def read_cluser_mapping_to_uniprot(file_mapping):
    result={}
    proteins_ids={}
    with open(file_mapping) as f:
        for line in f:
            if "qseqid" not in line and line.strip():
                split_line=line.split("\t")
                cluster_name=split_line[0].split(";")[0]
                family=split_line[-4]
                assigned_protein_id=split_line[-3]
                # if len(split_line)==15:
                #     pident=split_line[2]
                # else:
                pident=max([float(i) for i in split_line[4:6] if i])
                print(len(split_line))
                if "|" in assigned_protein_id:
                    assigned_protein_id=assigned_protein_id.split("|")[1]
                if assigned_protein_id not in result:
                    result[assigned_protein_id]=set()
                result[assigned_protein_id].add(cluster_name)
                proteins_ids[assigned_protein_id]={"pident_BSH":pident, "family":family}
    return result, proteins_ids

def merge_clusters_based_on_mapping(genes_counts, clusters_mapping):
    new_gene_counts={}
    clusters_info={}
    for probe, genes_count in genes_counts.items():
        new_gene_counts[probe]={}
        for protein, clusters in clusters_mapping.items():
            print(protein, clusters)
            clusters=list(clusters)
            gene_no=0
            first_cluster=sorted(clusters)[0]
            if first_cluster not in clusters_info:
                clusters_info[first_cluster]={"proteins":set(), "clusters":set()}
            clusters_info[first_cluster]["proteins"].add(protein)
            clusters_info[first_cluster]["clusters"]=clusters
    
            new_gene_counts[probe][first_cluster]=0
            for cl in clusters:
                new_gene_counts[probe][first_cluster]+=genes_count[cl]

    return new_gene_counts, clusters_info
    
def save_new_stats(merged_stats, out_file):
    file_list = list(merged_stats.keys())
    genes = set()
    for file, genes_data in merged_stats.items():
        genes = genes.union(set(genes_data.keys()))
    with open(out_file, "w") as f:
        f.write("genes;")
        f.write(";".join(file_list))
        f.write("\n")
        for gene in list(genes):
            f.write(f"{str(gene)};")
            for file_name in file_list:
                f.write(f"{str(merged_stats[file_name].get(gene, 0))};")
            f.write("\n")


def read_metabolites(count_metabolites_files):
    res = {}
    with open(count_metabolites_files) as f:
        header = f.readline()
        header = header.strip().split(";")[2:]
        for l in f:
            if l.strip():
                line = l.strip().split(";")
                probe_name = line[0].strip()
                res[probe_name] = {}
                res[probe_name]["metabolites"] = {}
                for e, i in enumerate(header):
                    if line[e + 2]:
                        res[probe_name]["metabolites"][i] = line[e + 2].replace(",", '.')
                    else:
                        res[probe_name]["metabolites"][i] = None

    return res


def select_files_of_metabolites(genes_counts, metabolites):
    for file, gene_data in genes_counts.items():
        # mgshot_S4358Nr27.1.fastq.gz.shotgun.tsv
        #print(file)
        file_name = file.split("_", 1)[1].split(".")[0].strip()
        #file_id = file.split("_", 1)[1].split(".")[1].strip()
        #print(file_name, file_id, file_name in metabolites)
        if file_name in metabolites:
            metabolites[file_name]["genes"] = gene_data
    return metabolites


def count_pearson_corr(metabolites):
    all_metabolites = set()
    all_genes = set()
    all_files = set()

    for files, metabolites_data in metabolites.items():
        all_metabolites = all_metabolites.union(set(metabolites_data["metabolites"].keys()))
        #print(files, metabolites_data)
        all_genes = all_genes.union(set(metabolites_data["genes"].keys()))
        all_files.add(files)
    all_files = list(all_files)
    tests = {}  # (metabo, gene, no):{metabo_values:[], gene_values:[], pval:[], corr_value:[]}
    metabolites_sets = {m: [] for m in list(all_metabolites)}
    genes_1 = {m: [] for m in list(all_genes)}
    for file in all_files:
        for metabolite, metabolites_data in metabolites[file]["metabolites"].items():
            metabolites_sets[metabolite].append(metabolites_data)
        for gene_name, gene_no in metabolites[file]["genes"].items():
            genes_1[gene_name].append(gene_no)
    for metabolite, metabolites_sets in metabolites_sets.items():
        p_val_list=[]
        for gene_name, gene_sets in genes_1.items():
            metabolites_list=[]
            genes_no_list=[]
            for e, gene_no in enumerate(gene_sets):
                if metabolites_sets[e] is not None and gene_no is not None:
                    metabolites_list.append(float(metabolites_sets[e]))
                    genes_no_list.append(float(gene_no))
            # CPM
            
            cpm = [i* 1e6 for i in genes_no_list]
            
            # log2(CPM)
            cpm = np.array(cpm, dtype=float)
            log_cpm = np.log2(cpm+1)
            #print(metabolites_list, genes_no_list)
            metabolites = np.array(metabolites_list, dtype=float)
            print(metabolite, metabolites)
            print("Min:", np.nanmin(metabolites))
            print("Max:", np.nanmax(metabolites))
            print("NaN:", np.isnan(metabolites).sum())
            print("Inf:", np.isinf(metabolites).sum())
            print("Negative:", (metabolites < 0).sum())
            metabolites=np.log2(metabolites+1)
            correlation = spearmanr(metabolites, log_cpm)
            tests[(metabolite, gene_name)] = {"metabo_values": metabolites_list, "gene_values": genes_no_list,
                                                 "pval": correlation.pvalue, "corr_value": correlation.statistic}

            p_val_list.append((gene_name, correlation.pvalue))
        #calculate_q_val

        
        pvals = [i[1] for i in p_val_list]
        
        # usunięcie nanów
        valid = []
        
        for g, p in p_val_list:
            try:
                p = float(p)
                if not np.isnan(p):
                    valid.append((p, g))
            except (ValueError, TypeError):
                print("not valid")
                continue        
        pvals_valid = [x[0] for x in valid]
        
        reject, pvals_corrected, _, _ = multipletests(
            pvals_valid,
            alpha=0.05,
            method='fdr_bh'
        )
        
        for e, (p_val, gene_name) in enumerate(valid):
            tests[(metabolite, gene_name)]["qval_value"] = pvals_corrected[e]
            tests[(metabolite, gene_name)]["qval_significance"] = reject[e]

            

    return tests

import requests

def get_protein_info(proteins_ids, batch_size=10):

    uniprot_accs = list(proteins_ids.keys())
    
    url = "https://rest.uniprot.org/uniprotkb/search"

    for i in range(0, len(uniprot_accs), batch_size):

        batch = uniprot_accs[i:i + batch_size]
        batch=[i for i in batch if i!="BSH"]

        params = {
            "query": "accession:(" + " OR ".join(batch) + ")",
            "format": "tsv",
            "fields": "accession,protein_name,organism_name",
            "size": 500
        }

        response = requests.get(url, params=params, timeout=120)

        if not response.ok:
            print(f"Error in batch {i // batch_size + 1}")
            print(response.text)
            response.raise_for_status()

        lines = response.text.strip().splitlines()

        for line in lines[1:]:
            if "deleted" not in line:
                print(line)
                accession, protein_name, organism_name = line.split("\t")
    
                proteins_ids[accession]["protein_name"] = protein_name
                proteins_ids[accession]["organism_name"] = organism_name

    return proteins_ids

def save_corr_stats(corr, out_file, clusters_info, proteins_ids):
    with open(out_file, "w") as f:
        f.write("metabolite;gene_name;corr_value;pval;qval_value;qval_significance;representative protein UniprotAcc;Protein name; Organism name;pident_blast; pident_bsh;family;similar gene clusters\n")
        for datasets_info, datasets_data in corr.items():
            f.write(';'.join([str(i) for i in list(datasets_info)]))
            f.write(";")
            f.write(str(datasets_data['corr_value']))
            f.write(";")
            f.write(str(datasets_data['pval']))
            f.write(";")
            f.write(str(datasets_data.get('qval_value', "nan")))
            f.write(";")
            f.write(str(datasets_data.get('qval_significance', "nan")))
            f.write(";")
            f.write(",".join(clusters_info[datasets_info[1]]["proteins"]))
            proteins=clusters_info[datasets_info[1]]["proteins"]
            f.write(";")
            f.write(','.join([proteins_ids[i].get("protein_name", "deleted") for i in proteins]))
            f.write(";")
            f.write(','.join([proteins_ids[i].get("organism_name", "deleted") for i in proteins]))
            f.write(";")
            f.write(','.join([proteins_ids[i].get("pident_BSH", "deleted") for i in proteins]))
            f.write(";")
            f.write(','.join([proteins_ids[i].get("family", "deleted") for i in proteins]))
            f.write(";")
            f.write(str(clusters_info[datasets_info[1]].get("pident_blast", "del")))
            f.write(";")

            f.write(",".join(clusters_info[datasets_info[1]]["clusters"]))
            # f.write(";")
            # f.write(','.join([str(i) for i in datasets_data['metabo_values']]))
            # f.write(";")
            # f.write(','.join([str(i) for i in datasets_data['gene_values']]))
            f.write("\n")

def read_blast_stats(bsh_stats, clusters_info):
    with open(bsh_stats) as f:
        for line in f:
            if line.strip():
                line_split=line.split("\t")
                cluster_no=line_split[0].split(";")[0]
                uniprot=line_split[1].split("|")[1]
                pident=line_split[3]
                if cluster_no in clusters_info:
                    if uniprot in clusters_info[cluster_no]["proteins"]:
                        if "pident_blast" not in clusters_info[cluster_no]:
                            clusters_info[cluster_no]["pident_blast"]=float(pident)
                        elif clusters_info[cluster_no]["pident_blast"]<float(pident):
                            clusters_info[cluster_no]["pident_blast"]=float(pident) 
    return(clusters_info)
    
@click.command()
@click.option('--count_genes_file', default="./", help='File with count of genes.')
@click.option('--count_metabolites_files', default="./", help='File with count of genes.')
@click.option('--out_file', default="./", help='Out file with normalised genes statistics.')
@click.option('--file_mapping', default="./comparison_ABC.csv", help='Out file with normalised genes statistics.')
@click.option('--bsh_stats', default="./", help='File with count of genes.')

def count_corr(count_genes_file, out_file, count_metabolites_files, file_mapping, bsh_stats):
    genes_counts = read_stats_file(count_genes_file)
    clusters_mapping, protein_ids=read_cluser_mapping_to_uniprot(file_mapping)
    proteins_ids=get_protein_info(protein_ids)
    genes_counts, clusters_info=merge_clusters_based_on_mapping(genes_counts, clusters_mapping)
    metabolites = read_metabolites(count_metabolites_files)
    genes_of_metabolites = select_files_of_metabolites(genes_counts, metabolites)
    new_metabo={}
    for file, metabo_data in genes_of_metabolites.items():
        if len(metabo_data.keys())>1:
            #print(metabo_data.keys())
            new_metabo[file]=metabo_data
    genes_of_metabolites=new_metabo
    corr = count_pearson_corr(genes_of_metabolites)
    clusters_info=read_blast_stats(bsh_stats, clusters_info)
    
    #input()
    #genes_of_metabolites = select_files_of_metabolites(genes_counts, metabolites)
    #print(genes_of_metabolites.keys())
    #print(genes_of_metabolites["S4358Nr1"])
    #print(genes_counts)
    
    save_corr_stats(corr, out_file, clusters_info, proteins_ids)


if __name__ == "__main__":
    count_corr()
