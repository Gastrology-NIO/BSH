import os

import click
import numpy as np
from ete3 import NCBITaxa
import six


class CorrInfo:
    def __init__(self, line, metabolite, cluster, pval, corr, pair, metabo_values_list, gene_values_list):
        #             f.write(','.join([str(i) for i in datasets_data['metabo_values']]))
        #             f.write(";")
        #             f.write(','.join([str(i) for i in datasets_data['gene_values']]))
        self.metabolite = metabolite
        self.cluster = cluster
        self.line = line
        self.pval = pval
        self.corr = corr
        self.pair = pair
        self.metabo_values_list = metabo_values_list.split(',')
        self.gene_values_list = gene_values_list.split(',')
        self.ctrl = None
        self.fisher = None
        self.zou = None
        self.lines = []

    def count(self):
        fisher = independent_corr(self.corr, self.ctrl.corr, len(self.metabo_values_list),
                                  len(self.ctrl.metabo_values_list), method='fisher')
        zou = independent_corr(self.corr, self.ctrl.corr, len(self.metabo_values_list),
                               len(self.ctrl.metabo_values_list), method='zou')
        self.fisher = f"{fisher[0]};{fisher[1]}"
        self.zou = f"{zou[0]};{zou[1]}"


def read_corr(corr_file):
    result = {}
    with open(corr_file) as f:
        for l in f:
            line = l.strip()
            if line:
                if "pval" not in line:
                    data = line.split(";")
                    metabolite = data[0]
                    cluster_no = data[1].split("_")[-1]
                    if cluster_no not in result:
                        result[cluster_no] = {}
                    if metabolite not in result[cluster_no]:
                        result[cluster_no][metabolite] = set()
                    result[cluster_no][metabolite].add(
                        CorrInfo(line=line.strip(), metabolite=data[0], cluster=data[1], pval=float(data[3]),
                                 corr=float(data[2]), pair=data[2], metabo_values_list=data[6],
                                 gene_values_list=data[7]))
    return result


def compare_cntrl_corr(corr_file, cntrl_corr_info):
    for cluster_no, metabolit_relation in corr_file.items():
        for metabolite, metabolite_cl_data in metabolit_relation.items():
            for cntrl_data in cntrl_corr_info[cluster_no][metabolite]:
                for searched_data in metabolite_cl_data:
                    searched_data.ctrl = cntrl_data
                    searched_data.count()
    return corr_file


def get_database_sequences_info(fasta_database, protein_names=set(), result=dict()):
    with open(fasta_database) as f:
        for l in f:
            line = l.strip()
            if line.startswith(">"):
                uniprot_id = line.split("|")[1]
                if protein_names:
                    if uniprot_id in protein_names:
                        protein_name = line.split("OS=")[0].split(" ", 1)[-1]
                        organism_name = line.split("OS=", 1)[1].split("OX=")[0]
                        organism_taxid = line.split("OX=")[-1].split(" ")[0]
                        result[uniprot_id] = {"protein_name": protein_name, "organism_name": organism_name,
                                              "organism_taxid": organism_taxid}
                else:
                    protein_name = line.split("OS=")[0].split(" ", 1)[-1]
                    organism_name = line.split("OS=", 1)[1].split("OX=")[0]
                    organism_taxid = line.split("OX=")[-1].split(" ")[0]
                    result[uniprot_id] = {"protein_name": protein_name, "organism_name": organism_name,
                                          "organism_taxid": organism_taxid}
    return result


def read_proteins_in_analyse(blast_table, all_proteins):
    with open(blast_table) as f:
        for l in f:
            line = l.strip()
            if line:
                line = line.split()
                cl_no = line[0].split(";")[0]
                protein = line[1].split("_")[0]
                all_proteins.add(protein)
    return all_proteins


def read_blast_result(blast_table, description, database_fasta_info, ncbi, result=None, main_taxids=None, representation_main_taxids=None, enzyme_proteins=None):
    # 5857992;4123967_2       tr|E3ZSC8|E3ZSC8_LISSE  57.143  112     48      0       1       112     1       112     2.67e-45        145
    # 2712365;4750170_2       A0A2K8V8T9_9FLAO        34.688  369     226     5       34      391     6       370     2.22e-75        237

    if result is None:
        result = {}
    with open(blast_table) as f:
        for l in f:
            line = l.strip()
            if line:
                line = line.split()
                print("read_blast_result", line)
                cl_no = line[0].split(";")[0]
                protein = line[1].split("_")[0]
                pident = float(line[2])
                eval = float(line[10])
                if enzyme_proteins is None:
                    pass
                else:

                    if protein in enzyme_proteins:
                        pass
                    else:
                        continue

                # print(database_fasta_info[protein]['organism_taxid'])
                if main_taxids:
                    print("main_taxids")

                    lineage =ncbi.get_lineage(database_fasta_info[protein]['organism_taxid'])
                    print("protein info", protein, database_fasta_info[protein]['organism_taxid'], database_fasta_info[protein]['organism_taxid'] in main_taxids, lineage, set(lineage).intersection(set(main_taxids)))
                    
                    if set(lineage).intersection(set(main_taxids)):
                        print("jest", cl_no in result)
                        if cl_no in result:
                            print("jest2", result[cl_no], pident)
                            if result[cl_no]["pident"] < pident:
                                result[cl_no]["pident"] = pident
                                result[cl_no]["protein"] = [protein]
                                result[cl_no]["description"] = [description]
                                result[cl_no]["blast_table"] = [blast_table]
                                result[cl_no]["eval"] = [eval]
                                result[cl_no]["representation"] = [representation_main_taxids[list(set(lineage).intersection(set(main_taxids)))[0]]]

                            if result[cl_no]["pident"] == pident:
                                if protein not in result[cl_no]["protein"]:
                                    result[cl_no]["pident"] = pident
                                    result[cl_no]["protein"].append(protein)
                                    result[cl_no]["description"].append(description)
                                    result[cl_no]["blast_table"].append(blast_table)
                                    result[cl_no]["eval"].append(eval)
                                    result[cl_no]["representation"].append(representation_main_taxids[list(set(lineage).intersection(set(main_taxids)))[0]])

                        else:
                            result[cl_no] = {"pident": pident, "protein": [protein], "description": [description],
                                                                                         "blast_table": [blast_table], "eval": [eval], 
                                            "representation":[representation_main_taxids[list(set(lineage).intersection(set(main_taxids)))[0]]]}
                else:
                    print("without tax")
                    try:
                        lineage = ncbi.get_lineage(database_fasta_info[protein]['organism_taxid'])
                        if 2759 not in lineage:
                            if cl_no in result:
                                if result[cl_no]["pident"] < pident:
                                    result[cl_no]["pident"] = pident
                                    result[cl_no]["protein"] = [protein]
                                    result[cl_no]["description"] = [description]
                                    result[cl_no]["blast_table"] = [blast_table]
                                    result[cl_no]["eval"] = [eval]
                                if result[cl_no]["pident"] == pident:
                                    if protein not in result[cl_no]["protein"]:
                                        result[cl_no]["protein"].append(protein)
                                        result[cl_no]["description"].append(description)
                                        result[cl_no]["blast_table"].append(blast_table)
                                        result[cl_no]["eval"].append(eval)
                            else:
                                result[cl_no] = {"pident": pident, "protein": [protein], "description": [description],
                                                 "blast_table": [blast_table], "eval": [eval]}
                    except Exception as e:
                        print(e)
                        if cl_no in result:
                            if result[cl_no]["pident"] < pident:
                                result[cl_no]["pident"] = pident
                                result[cl_no]["protein"] = [protein]
                                result[cl_no]["description"] = [description]
                                result[cl_no]["blast_table"] = [blast_table]
                        else:
                            result[cl_no] = {"pident": pident, "protein": [protein], "description": [description],
                                             "blast_table": [blast_table]}
    return result


import numpy as np
from scipy.stats import t, norm
from math import atanh, pow
from numpy import tanh


def rz_ci(r, n, conf_level=0.95):
    zr_se = pow(1 / (n - 3), .5)
    moe = norm.ppf(1 - (1 - conf_level) / float(2)) * zr_se
    zu = atanh(r) + moe
    zl = atanh(r) - moe
    return tanh((zl, zu))


def independent_corr(xy, ab, n, n2=None, twotailed=True, conf_level=0.95, method='fisher'):
    """
    Calculates the statistic significance between two independent correlation coefficients
    @param xy: correlation coefficient between x and y
    @param xz: correlation coefficient between a and b
    @param n: number of elements in xy
    @param n2: number of elements in ab (if distinct from n)
    @param twotailed: whether to calculate a one or two tailed test, only works for 'fisher' method
    @param conf_level: confidence level, only works for 'zou' method
    @param method: defines the method uses, 'fisher' or 'zou'
    @return: z and p-val
    """

    if method == 'fisher':
        xy_z = 0.5 * np.log((1 + xy) / (1 - xy))
        ab_z = 0.5 * np.log((1 + ab) / (1 - ab))
        if n2 is None:
            n2 = n

        se_diff_r = np.sqrt(1 / (n - 3) + 1 / (n2 - 3))
        diff = xy_z - ab_z
        z = abs(diff / se_diff_r)
        p = (1 - norm.cdf(z))
        if twotailed:
            p *= 2

        return z, p
    elif method == 'zou':
        L1 = rz_ci(xy, n, conf_level=conf_level)[0]
        U1 = rz_ci(xy, n, conf_level=conf_level)[1]
        L2 = rz_ci(ab, n2, conf_level=conf_level)[0]
        U2 = rz_ci(ab, n2, conf_level=conf_level)[1]
        lower = xy - ab - pow((pow((xy - L1), 2) + pow((U2 - ab), 2)), 0.5)
        upper = xy - ab + pow((pow((U1 - xy), 2) + pow((ab - L2), 2)), 0.5)
        return lower, upper


def fix_corr(corr_info, blast_result, database_info, save_old_line=True):
    for cluster_no, metabolit_relation in corr_info.items():
        for metabolite, metabolite_cl_data in metabolit_relation.items():
            for corr_data in metabolite_cl_data:
                line = corr_data.line
                if not save_old_line:
                    line = f"{corr_data.metabolite};{corr_data.cluster};{corr_data.pair};{corr_data.pval};{corr_data.corr}"
                    corr_data.line = line
                if cluster_no in blast_result:
                    line = corr_data.line
                    print(blast_result[cluster_no]['protein'], blast_result[cluster_no])
                    line += f";{blast_result[cluster_no]['pident']}"
                    line += f";{','.join(blast_result[cluster_no]['protein'])}"
                    # line += f";https://www.uniprot.org/uniprotkb/{blast_result[cluster_no]['protein'][e]}/entry"
                    line += f";{','.join(blast_result[cluster_no]['description'])}"
                    line += f";{','.join(blast_result[cluster_no]['blast_table'])}"
                    line += f";{','.join([database_info[i]['protein_name'] for i in blast_result[cluster_no]['protein']])}"
                    line += f";{','.join([database_info[i]['organism_name'] for i in blast_result[cluster_no]['protein']])}"
                    line += f";{','.join([database_info[i]['organism_taxid'] for i in blast_result[cluster_no]['protein']])}"
                    if "representation" in blast_result[cluster_no]:
                        line+='"'+f";{'|'.join(list(set(blast_result[cluster_no]['representation'])))}"+'"'

                    if corr_data.ctrl is not None:
                        line += f";{corr_data.ctrl.corr}"
                        line += f";{corr_data.ctrl.pval}"
                        line += f";{corr_data.zou}"
                        line += f";{corr_data.fisher}"
                        line += f";{str(len(corr_data.metabo_values_list))}"
                        line += f";{str(len(corr_data.gene_values_list))}"
                    # corr_data.lines.append(line)
                    yield line
                    line = corr_data.line
                else:
                    line += f";"
                    line += f";"
                    line += f";"
                    line += f";"
                    line += f";"
                    line += f";"
                    line += f";"
                    if corr_data.ctrl is not None:
                        line += f";{corr_data.ctrl.corr}"
                        line += f";{corr_data.ctrl.pval}"
                        line += f";{corr_data.zou}"
                        line += f";{corr_data.fisher}"
                        line += f";{str(len(corr_data.metabo_values_list))}"
                        line += f";{str(len(corr_data.gene_values_list))}"
                    # corr_data.lines.append(line)
                    yield line
                    line = corr_data.line
    return corr_info


def save_corr(fixed_corr, out_file):
    with open(out_file, "w") as f:
        for line in fixed_corr:
            f.write(line)
            f.write("\n")
    # with open(out_file, "w") as f:
    #     for cluset_no, metabolit_relation in fixed_corr.items():
    #         for metabolite, metabolite_cl_data in metabolit_relation.items():
    #             for corr_data in metabolite_cl_data:
    #                 if corr_data.lines:
    #                     for line in corr_data.lines:
    #                         f.write(line)
    #                         f.write("\n")
    #                 else:
    #                     f.write(corr_data.line)
    #                     f.write("\n")


@click.command()
@click.option('--corr_file', default="./", help='Folder with BLAST files.')
@click.option('--corr_ctrl_file', default="", help='Folder with BLAST files.')
@click.option('--blast_files', default={}, help='')
@click.option('--fasta_databases', default="./", help='')
@click.option('--out_file', default="./", help='Out file with genes statistics.')
@click.option('--main_taxids', default="", help='Out file with genes statistics.')
@click.option('--selected_proteins', default="", help='Out file with genes statistics.')

@click.option('--alfa', default=0.05, help='Out file with genes statistics.')
def main(corr_file, corr_ctrl_file, blast_files, fasta_databases, out_file, main_taxids,selected_proteins,  alfa):
    print("main taxids")
    main_taxids_representation={}
    if main_taxids:
        
        tmp_ids = set()
        padj_no = 0
        if os.path.exists(main_taxids):
            with open(main_taxids) as f:
                for l in f:
                    if l.strip():
                        line = l.split(";")
                        
                        if line[0]:
                            print(line, padj_no, line[4].replace(',', '.'))
                            if "Species" in l:
                                for e, i in enumerate(line):
                                    if i == "padj":
                                        padj_no = e
                                    
                        #print(padj_no, line, line[padj_no])
                            elif float(line[4].replace(',', '.'))< alfa:
                                main_taxids_representation[int(line[1].replace("OTU", "").replace('"', ''))]=l.strip().replace(",", ".").replace(";", ",")
                            
                                tmp_ids.add(int(line[1].replace("OTU", "").replace('"', '')))
                                print(tmp_ids)
        main_taxids = list(tmp_ids)
    enzyme_proteins=set()
    if selected_proteins:
        with open(selected_proteins) as f:
            for l in f:
                if not "Entry" in l:

                    #line=l.strip().split("\t")

                    line=l.strip().split("\n")
                    uniprot_acc=line[0]
                    enzyme_proteins.add(uniprot_acc)
    print("main taxids", main_taxids)
    blast_files = eval(blast_files)
    print("read corr")

    corr_info = read_corr(corr_file)
    print("read corr_cntrl")

    if os.path.exists(corr_ctrl_file):
        cntrl_corr_info = read_corr(corr_ctrl_file)
        corr_info = compare_cntrl_corr(corr_info, cntrl_corr_info)
    print("read sequences info")
    blast_result = dict()
    database_fasta_info = dict()
    all_proteins = set()
    for blast_table, description in blast_files.items():
        if blast_table:
            all_proteins = read_proteins_in_analyse(blast_table, all_proteins)
    print(len(list(all_proteins)))
    for fasta_db in fasta_databases.split(","):
        database_fasta_info = get_database_sequences_info(fasta_db, all_proteins, database_fasta_info)
    ncbi = NCBITaxa()
    # ncbi.update_taxonomy_database()
    print("read blast results")
    print(blast_files)
    for blast_table, description in blast_files.items():
        if blast_table:
            blast_result = read_blast_result(blast_table, description, database_fasta_info, ncbi, blast_result,
                                             main_taxids, main_taxids_representation, enzyme_proteins)
    print(blast_result)
    print("fix corr")

    fixed_corr = fix_corr(corr_info, blast_result, database_fasta_info)
    save_corr(fixed_corr, out_file)


if __name__ == '__main__':
    main()
