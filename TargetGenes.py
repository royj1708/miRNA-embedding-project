import csv
from csv import reader
from scipy.stats import hypergeom
import os

AMOUNT_OF_GENES = 43884  # was checked separately


def create_miRNA_targets_dataset(miRNA_targets_path: str, prob_threshold: float, new_file_path: str):
    miRNA_targets_dict = {}
    first_line = True
    with open(miRNA_targets_path, 'r') as targets_file:
        for record in targets_file:
            if first_line:
                first_line = False
                continue

            record = record.split()
            miRNA = record[0]
            gene = record[1]
            prob = float(record[-1])

            if prob < prob_threshold:
                continue

            if miRNA not in miRNA_targets_dict:
                miRNA_targets_dict[miRNA] = set()
                miRNA_targets_dict[miRNA].add(gene)

            else:
                miRNA_targets_set = miRNA_targets_dict[miRNA]
                miRNA_targets_set.add(gene)

    with open(new_file_path, 'a') as new_file:
        writer = csv.writer(new_file)
        for miR in miRNA_targets_dict:
            writer.writerow([miR] + list(miRNA_targets_dict[miR]))


def intersection_size(list_1: list, list_2: list):
    return len(set(list_1) & set(list_2))


def hypergeomtric_statistical_test(miRNA_targets_csv_path: str):
    miRNA_targets_dict = {}
    with open(miRNA_targets_csv_path, 'r') as miRNA_targets_file:
        for record in miRNA_targets_file:
            record = record.split(",")
            miRNA = record[0]
            targets = record[1:]
            miRNA_targets_dict[miRNA] = targets

    with open('/Users/royjudes/Desktop/miRNA embedding project/similar_miRNAs_by_targets_4.csv', 'a', newline='') as similar_by_targets_file:
        writer = csv.writer(similar_by_targets_file)
        writer.writerow(['miRNA_a', 'miRNA_b', 'n', 'N', 'X', 'p_value'])
        already_computed = []
        counter = 0
        for miRNA_a in miRNA_targets_dict:
            already_computed.append(miRNA_a)
            miRNA_a_targets = miRNA_targets_dict[miRNA_a]
            miRNA_a_targets_size = len(miRNA_a_targets)

            for miRNA_b in miRNA_targets_dict:
                if miRNA_b in already_computed:
                    continue

                counter += 1
                if counter <= 3145726:
                    continue

                # if counter == 3145726:
                #     return

                miRNA_b_targets = miRNA_targets_dict[miRNA_b]
                overlap = intersection_size(miRNA_a_targets, miRNA_b_targets)
                miRNA_b_targets_size = len(miRNA_b_targets)

                p_val = hypergeom.sf(overlap-1, AMOUNT_OF_GENES, miRNA_a_targets_size, miRNA_b_targets_size)
                writer.writerow([miRNA_a, miRNA_b, miRNA_a_targets_size, miRNA_b_targets_size, overlap, p_val])


def merge_cossim_and_pvalue_datasets(pvals_directory_path: str):
    miRNAs_pvalues = {}
    for (root, dirs, files) in os.walk(pvals_directory_path, topdown=True):
        for file in files:
            first_line = True
            try:
                with open(os.path.join(root, file), 'r') as pvals_file:
                    for record in pvals_file:
                        if first_line:
                            first_line = False
                            continue

                        record = record.split(",")
                        miRNA_a = record[0]
                        miRNA_b = record[1]
                        pval = record[-1][0:-1]
                        miRNAs_pvalues[(miRNA_a, miRNA_b)] = pval

            except:
                continue

    with open("/Users/royjudes/Desktop/miRNA embedding project/mature_miRNAs_cossim_kidney.csv", 'r') as csv_input:
        with open("/Users/royjudes/Desktop/miRNA embedding project/kidney_miRNA_pairs_pvalues_cossim.csv", 'w') as csv_output:
            writer = csv.writer(csv_output)
            reader = csv.reader(csv_input)
            all = []
            row = next(reader)
            row.append("pvalue")
            all.append(row)
            for row in reader:
                miRNA_a = row[0]
                miRNA_b = row[2]
                try:
                    pval = miRNAs_pvalues[(miRNA_a, miRNA_b)] if (miRNA_a, miRNA_b) in miRNAs_pvalues else miRNAs_pvalues[(miRNA_b, miRNA_a)]

                except:
                    pval = 'None'

                row.append(pval)
                all.append(row)

            writer.writerows(all)


# def compute(miRNA_targets_path: str):
#     targets = set()
#     first_line = True
#     with open(miRNA_targets_path, 'r') as targets_file:
#         for record in targets_file:
#             if first_line:
#                 first_line = False
#                 continue
#
#             record = record.split()
#             gene = record[1]
#             targets.add(gene)
#
#     print(len(targets))


def create_mirna_genes_amount(path: str):
    miRNAs_genes = {}
    with open(path, 'r') as miRNA_targets_file:
        for record in miRNA_targets_file:
            record = record.split(",")
            miRNA = record[0]
            targets = record[1:]
            miRNAs_genes[miRNA] = len(targets)

    with open("/Users/royjudes/Desktop/miRNA embedding project/miRNAs_targets_amounts.csv", 'a') as amount_file:
        writer = csv.writer(amount_file)
        writer.writerow(["miRNA", "Targets"])
        for mir in miRNAs_genes:
            writer.writerow([mir, miRNAs_genes[mir]])



miRNA_targets_path = "/Users/royjudes/Desktop/miRNA embedding project/hsa_miRWalk_3UTR.txt"
new_file_path = "/Users/royjudes/Desktop/miRNA embedding project/matures_target_genes_csv.csv"
prob_threshold = 0.95

# create_miRNA_targets_dataset(miRNA_targets_path, prob_threshold, new_file_path)
# hypergeomtric_statistical_test(new_file_path)
merge_cossim_and_pvalue_datasets("/Users/royjudes/Desktop/miRNA embedding project/pvalues files")
create_mirna_genes_amount(new_file_path)
