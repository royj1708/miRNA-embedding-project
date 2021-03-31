import csv
from csv import reader


def create_miRNA_targets_dataset(miRNA_targets_path: str, prob_threshold: float):
    miRNA_targets_dict = {}
    with open(miRNA_targets_path, 'r') as targets_file:
        writer = csv.writer(targets_file)
        for record in targets_file:
            record = record.split()
            miRNA = record[0]
            gene = record[1]
            prob = float(record[-1])

            if prob < prob_threshold:
                continue

            if miRNA not in miRNA_targets_dict:
                miRNA_targets_dict[miRNA] = [miRNA, gene]
            else:
                miRNA_targets_list = miRNA_targets_dict[miRNA]
                if gene not in miRNA_targets_list:
                    miRNA_targets_list.append(gene)

        for miR in miRNA_targets_dict:
            writer.writerow(miRNA_targets_dict[miR])
