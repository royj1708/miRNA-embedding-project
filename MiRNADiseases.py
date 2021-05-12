import csv
from csv import reader
import os
from collections import defaultdict
import pandas as pd
from scipy.stats import pearsonr


# def create_miRNA_diseases_dataset(miRNA_diseases_path: str, new_file_path: str, fasta_file_path: str):
#     miRNA_diseases_dict = defaultdict(lambda: [])
#     first_line = True
#     with open(miRNA_diseases_path, 'r') as miRNA_diseases_file:
#         for record in miRNA_diseases_file:
#             if first_line:
#                 first_line = False
#                 continue
#
#             record = record.split()
#             hairpin = record[1]
#             disease = record[2]
#             miRNA_diseases_dict[hairpin].append(disease)
#
#     return miRNA_diseases_dict
#
#
# def get_matures_list_for_hairpins(fasta_file_path: str):
#     matures_hairpins_dict = defaultdict(lambda: [])
#     mimat_to_mirnas_dict = create_mature_names_dictionary(fasta_file_path)


def get_diseases_for_matures(matures_diseases_file_path: str):
    mirs_dict = {}
    with open(matures_diseases_file_path, 'r') as matures_diseases_file:
        for record in matures_diseases_file:
            try:
                record = record.split(",")
                mir = f"hsa-{record[0]}"
                diseases = [disease.replace('"', "").strip() for disease in record[1:]]
                mirs_dict[mir] = diseases
            except:
                continue

    return mirs_dict


# def create_mature_names_dictionary(fasta_file_path: str):
#     """
#     Creates a dictionary that maps the mature accession (MIMAT) to its name (hsa).
#
#     :param fasta_file_path: path to the fasta file that contains all of the different mature miRNAs, with both their
#     'MIMAT' and 'hsa' names.
#     :return: a dictionary that maps the mature accession (MIMAT) to its name (hsa).
#     """
#     mature_dict = {}
#     with open(fasta_file_path, 'r') as file:
#         for record in SeqIO.parse(file, 'fasta'):
#             description = record.description.split(" ")
#             mature_MIMAT = description[1]
#             mature_hsa = description[0]
#             if "hsa" not in mature_hsa:
#                 continue
#             mature_dict[mature_MIMAT] = mature_hsa
#
#     return mature_dict


def intersection_size(list_1: list, list_2: list):
    return len(set(list_1) & set(list_2))


def compute_jaccard_indice(mirs_diseases_dict: dict):
    with open('/Users/royjudes/Desktop/miRNA embedding project/kidney_miRNA_pairs_pvalues_cossim.csv', 'r') as pairs_file:
        first_line = True
        with open('/Users/royjudes/Desktop/miRNA embedding project/kidney_miRNA_pairs_pvalues_cossim_jaccard.csv', 'a') as pairs_file_with_jaccard:
            writer = csv.writer(pairs_file_with_jaccard)
            writer.writerow(['miRNA_a', 'avg_rpm_a', 'miRNA_b', 'avg_rpm_b', 'cos_sim', 'pvalue', 'jaccard'])
            for record in pairs_file:
                if first_line:
                    first_line = False
                    continue

                record = record.split(",")
                miRNA_a = record[0]
                miRNA_b = record[2]

                record[5] = record[5].replace("\n", "")

                try:
                    mir_a_diseases = mirs_diseases_dict[miRNA_a]
                    mir_b_diseases = mirs_diseases_dict[miRNA_b]
                    intersection = intersection_size(mir_a_diseases, mir_b_diseases)
                    jaccard_indice = intersection / (len(mir_a_diseases) + len(mir_b_diseases) - intersection)
                    record.append(jaccard_indice)
                    writer.writerow(record)
                    
                except:
                    record.append('None')
                    writer.writerow(record)


miRNA_diseases_file_path = "/Users/royjudes/Desktop/miRNA embedding project/mature_miRNA_HMDD.csv"
# mirna_diseases_dict = get_diseases_for_matures(miRNA_diseases_file_path)
# compute_jaccard_indice(mirna_diseases_dict)
