import os

import click
import numpy as np
from scipy.stats import pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests
import numpy as np

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
                    gene = line[0]
                    genes_no = line[1:]
                    for e, file in enumerate(file_names):
                        stats[file][gene] = float(genes_no[e])
    return stats


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
            metabolites=np.log2(metabolites+1)
            correlation = spearmanr(metabolites, log_cpm)
            tests[(metabolite, gene_name)] = {"metabo_values": metabolites_list, "gene_values": genes_no_list,
                                                 "pval": correlation.pvalue, "corr_value": correlation.statistic}
            if not np.isnan(correlation.pvalue):
                p_val_list.append((gene_name, correlation.pvalue))
            else:
                print(correlation.pvalue)
            print(correlation.pvalue)
        #calculate_q_val
        pvals = [i[1] for i in p_val_list]
        
        # usunięcie nanów
        
        
        #pvals_valid = [x[1] for x in valid
        print(pvals)
        if pvals:
            reject, pvals_corrected, _, _ = multipletests(
                pvals,
                alpha=0.05,
                method='fdr_bh'
            )
        
            for e, (gene_name,p_val) in enumerate(p_val_list):
                tests[(metabolite, gene_name)]["qval_value"] = pvals_corrected[e]
                tests[(metabolite, gene_name)]["qval_significance"] = reject[e]

            

    return tests


def save_corr_stats(corr, out_file):
    with open(out_file, "w") as f:
        f.write("metabolite;gene_name;corr_value;pval;metabo_list_values;gene_list_percentages\n")
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
            f.write(','.join([str(i) for i in datasets_data['metabo_values']]))
            f.write(";")
            f.write(','.join([str(i) for i in datasets_data['gene_values']]))
            f.write("\n")


@click.command()
@click.option('--count_genes_file', default="./", help='File with count of genes.')
@click.option('--count_metabolites_files', default="./", help='File with count of genes.')
@click.option('--out_file', default="./", help='Out file with normalised genes statistics.')
def count_corr(count_genes_file, out_file, count_metabolites_files):
    genes_counts = read_stats_file(count_genes_file)
    metabolites = read_metabolites(count_metabolites_files)
    genes_of_metabolites = select_files_of_metabolites(genes_counts, metabolites)
    new_metabo={}
    for file, metabo_data in genes_of_metabolites.items():
        if len(metabo_data.keys())>1:
            #print(metabo_data.keys())
            new_metabo[file]=metabo_data
    genes_of_metabolites=new_metabo
    #input()
    #genes_of_metabolites = select_files_of_metabolites(genes_counts, metabolites)
    #print(genes_of_metabolites.keys())
    #print(genes_of_metabolites["S4358Nr1"])
    #print(genes_counts)
    corr = count_pearson_corr(genes_of_metabolites)
    save_corr_stats(corr, out_file)


if __name__ == "__main__":
    count_corr()
