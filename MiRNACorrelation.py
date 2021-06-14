import csv
from csv import reader
import os
from collections import defaultdict
import pandas as pd
from scipy.stats import pearsonr


def compute_mirna_expression_correlation(mirna_expressions_file_path: str):
    """
    Computes the Pearson correlation between each two miRNAs according to their counts in the profiles.
    Creates a dataset with each miRNA pair and their correlation as a record.
    :param mirna_expressions_file_path: a path to the dataset of the miRNA expression profiles
    """
    with open('/Users/royjudes/Desktop/miRNA embedding project/kidney_miRNA_pairs_pvalues_cossim_jaccard.csv', 'r') as pairs_file:
        first_line = True
        with open('/Users/royjudes/Desktop/miRNA embedding project/kidney_miRNA_pairs_pvalues_cossim_jaccard_correlation.csv', 'a') as pairs_file_with_correlation:
            writer = csv.writer(pairs_file_with_correlation)
            writer.writerow(['miRNA_a', 'avg_rpm_a', 'miRNA_b', 'avg_rpm_b', 'cos_sim', 'pvalue', 'jaccard', 'correlation'])
            df = pd.read_csv(mirna_expressions_file_path)

            for record in pairs_file:
                if first_line:
                    first_line = False
                    continue

                try:
                    record = record.split(",")
                    miRNA_a = record[0]
                    miRNA_b = record[2]

                    record[6] = record[6].replace("\n", "")

                    correlation, _ = pearsonr(df[miRNA_a], df[miRNA_b])
                    record.append(correlation)
                    writer.writerow(record)

                except:
                    print(record)


def load_miRNA_embeddings_dictionary(path: str):
    """
    Reads a file with the miRNAs and their embedding vectors and stores it in a dictionary.
    :param path: the path to the embeddings file
    """
    with open(path, 'r', newline='') as embeddings_file:
        csv_reader = reader(embeddings_file)
        prev_row = []
        embeddings_dictionary = {}
        for row in csv_reader:
            if len(row) == 300:
                vector = []
                for value in row:
                    vector.append(float(value))
                embeddings_dictionary[prev_row[0]] = vector

            prev_row = row

    return embeddings_dictionary


# def compute_mirna_expression_correlation_embeddings(mirna_embeddings_file_path: str):
#     with open('/Users/royjudes/Desktop/miRNA embedding project/kidney_miRNA_pairs_pvalues_cossim_jaccard_correlation.csv', 'r') as pairs_file:
#         first_line = True
#         with open('/Users/royjudes/Desktop/miRNA embedding project/kidney_miRNA_pairs_pvalues_cossim_jaccard_correlation_emb.csv', 'a') as pairs_file_with_embedding_correlation:
#             writer = csv.writer(pairs_file_with_embedding_correlation)
#             writer.writerow(['miRNA_a', 'avg_rpm_a', 'miRNA_b', 'avg_rpm_b', 'cos_sim', 'pvalue', 'jaccard', 'correlation', 'correlation_emb'])
#             embeddings_dict = load_miRNA_embeddings_dictionary(mirna_embeddings_file_path)
#
#             for record in pairs_file:
#                 if first_line:
#                     first_line = False
#                     continue
#
#                 try:
#                     record = record.split(",")
#                     miRNA_a = record[0]
#                     miRNA_b = record[2]
#
#                     record[7] = record[7].replace("\n", "")
#
#                     correlation, _ = pearsonr(embeddings_dict[miRNA_a], embeddings_dict[miRNA_b])
#                     record.append(correlation)
#                     writer.writerow(record)
#
#                 except:
#                     print(record)


mirna_expressions_file_path = "/Users/royjudes/Desktop/miRNA embedding project/profiles_matures_dataset.csv"
# compute_mirna_expression_correlation(mirna_expressions_file_path)
mirna_embeddings_file_path = "/Users/royjudes/Desktop/miRNA embedding project/mature_miRNA_embeddings_kidney.csv"
compute_mirna_expression_correlation_embeddings(mirna_embeddings_file_path)
