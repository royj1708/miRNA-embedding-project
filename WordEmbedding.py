import os
from gensim.models.word2vec import Word2Vec
from multiprocessing import cpu_count
import numpy as np
import xlrd
from datetime import datetime
import pandas as pd
import csv
from DocumentsReduction import create_health_condition_dictionary


DOCUMENTS_DIRECTORY = "C:\\Users\\royj1\\Desktop\\University\\הנדסת מערכות מידע\\שנה ד\\פרויקט גמר\\דאטה\\Testing\\b"
VOCABULARY_FILE = "C:\\Users\\royj1\\Desktop\\University\\הנדסת מערכות מידע\\שנה ד\\פרויקט גמר\\דאטה\\all_mirnas.xls"
# LOG_FILE_PATH = "C:\\Users\\royj1\\Desktop\\University\\הנדסת מערכות מידע\\שנה ד\\פרויקט גמר\\דאטה\\Testing\\Log.txt"


def load_vocabulary():
    """
    Reads the miRNAs from an excel file and stores them in a list.
    :return: list of the miRNAs (the vocabulary).
    """
    vocabulary = []
    vocab_file = xlrd.open_workbook(VOCABULARY_FILE).sheet_by_index(0)
    for i in range(vocab_file.nrows):
        vocabulary.append(vocab_file.cell_value(i, 0))

    return vocabulary


def create_list_from_document(path: str):
    document_as_list = []
    with open(path, 'r') as document:
        document_as_list = document.read().split()
    # for line in document:
    #     document_as_list.extend(line.split())

    return document_as_list


def convert_documents_to_lists():
    list_of_all_documents = []
    for (root, dirs, files) in os.walk(DOCUMENTS_DIRECTORY, topdown=True):
        for file in files:
            try:
                if file[-4:] == ".txt":
                    list_of_all_documents.append(create_list_from_document(os.path.join(DOCUMENTS_DIRECTORY, file)))

            except:
                continue

    return list_of_all_documents


def convert_documents_csv_to_lists(path_to_csv):
    list_of_all_documents = []
    samples_df = pd.read_csv(path_to_csv)
    samples_df = samples_df.drop('sample_type', axis=1)
    miRNAs = list(samples_df.columns)

    for index, profile in samples_df.iterrows():
        profile_vector = []
        for miRNA in miRNAs:
            profile_vector.append(profile[miRNA])

        list_of_all_documents.append(profile_vector)

    return list_of_all_documents


def convert_documents_to_dictionary():
    dict_of_all_documents = {}
    for (root, dirs, files) in os.walk(DOCUMENTS_DIRECTORY, topdown=True):
        for file in files:
            try:
                if file[-4:] == ".txt":
                    dict_of_all_documents[file] = (create_list_from_document(os.path.join(DOCUMENTS_DIRECTORY, file)))

            except:
                continue

    return dict_of_all_documents


def compute_wmd_for_profiles(profiles_as_list, keys, health_condition_dictionary, model):
    with open('profiles_distances.csv', 'a', newline='') as distances_file:
        writer = csv.writer(distances_file)
        for i in range(0, len(profiles_as_list), 1):
            profile_a = profiles_as_list[i]
            profile_a_condition = 'N' if health_condition_dictionary[f"{keys[i][:-13]}.txt"] == 'Normal' else 'T'
            if profile_a_condition == 'N':
                continue
            print(keys[i])
            for j in range(i+1, len(profiles_as_list), 1):
                profile_b = profiles_as_list[j]
                distance = model.wmdistance(profile_a, profile_b)
                profile_b_condition = 'N' if health_condition_dictionary[f"{keys[j][:-13]}.txt"] == 'Normal' else 'T'
                record = [profile_a_condition, profile_b_condition, distance]
                writer.writerow(record)


# def write_into_log_file(size: int, window: int, min_count: int, model: str, results):
#     content = f"Time: {datetime.now().strftime('%d/%m/%Y %H:%M:%S')}\n" \
#               f"Results - {model}:\n" \
#               f"Embedding Vector Size: {size}, Context Window Size: {window}, Minimal amount of appearances " \
#               f"to be included: {min_count}, Accuracy: {results}\n\n"
#
#     with open(LOG_FILE_PATH, 'a') as log_file:
#         log_file.write(content)


def main():
    all_profiles = convert_documents_to_dictionary()
    health_condition_dictionary = create_health_condition_dictionary(
        '/Users/royjudes/Desktop/miRNA embedding project/gdc_sample_sheet.tsv')
    profiles_as_list = list(all_profiles.values())
    keys = list(all_profiles.keys())
    model = Word2Vec(profiles_as_list, min_count=0, size=300, window=5, workers=cpu_count())
    compute_wmd_for_profiles(profiles_as_list, keys, health_condition_dictionary, model)


    # profiles = convert_documents_csv_to_lists('/Users/royjudes/Desktop/miRNA embedding project/kidney_samples.csv')

    # profile_normal_a = all_profiles[11]
    # profile_normal_b = all_profiles[12]
    # profile_tumor_a = all_profiles[1]
    # profile_tumor_b = all_profiles[2]

    # two_normals_distance = model.wmdistance(profile_normal_a, profile_normal_b)
    # two_tumors_distance = model.wmdistance(profile_tumor_a, profile_tumor_b)
    # normal_tumor_distance = model.wmdistance(profile_normal_a, profile_tumor_a)
    # tumor_normal_distance = model.wmdistance(profile_tumor_b, profile_normal_b)

    # print(f"distance between two normals: {two_normals_distance}")
    # print(f"distance between two tumors: {two_tumors_distance}")
    # print(f"distance between normal and tumor a: {normal_tumor_distance}")
    # print(f"distance between normal and tumor b: {tumor_normal_distance}")


    # model = Word2Vec(profiles, min_count=0, size=300, window=5, workers=cpu_count())

    # samples_file = '/Users/royjudes/Desktop/miRNA embedding project/kidney_samples.csv'
    # samples_df = pd.read_csv(samples_file)
    # sample_types = samples_df["sample_type"]
    # samples_df = samples_df.drop('sample_type', axis=1)
    # miRNAs = list(samples_df.columns)

    # for index, profile in samples_df.iterrows():
    #     miRNAs_amount_in_profile = 0
    #     profile_vector = []
    #     for i in range(300):
    #         profile_vector.append(0)
    #     for miRNA in miRNAs:
    #         try:
    #             miRNA_vector = model[miRNA]  # the embedding vector of the miRNA
    #             expression = float(profile[miRNA])  # the miRNA expression in the profile
    #             miRNAs_amount_in_profile += expression  # sum of all miRNAs in the profile
    #             profile_vector = np.add(profile_vector, np.array(miRNA_vector)*expression)

    #         except:
    #             continue

    #     profile_vector /= miRNAs_amount_in_profile
    #     profile_vector = list(profile_vector)
    #     profile_vector.append(sample_types[index])

    #     with open('embedded_profiles_0_300_5.csv', 'a', newline='') as file:
    #         writer = csv.writer(file)
    #         writer.writerow(profile_vector)


    # distance = model.wmdistance(all_profiles[0], all_profiles[1])
    # print(distance)

main()
    # for size in [100, 200, 300, 400, 500]:
    #     for window in [3, 4, 5, 6, 7]:
    #         for min_count in [0, 5, 10, 15, 20]:
    #             word2vec_model = Word2Vec(all_profiles, min_count=min_count, window=window, size=size, workers=cpu_count())
                # print(word2vec_model.most_similar(positive=['hsa-let-7a-1', 'hsa-let-7a-2'], negative=['hsa-let-7a-3'], topn=1))

                # write_into_log_file(size, window, min_count, None)
                # we need to take all the miRNA vectors and try the classification models, then record the accuracy rate
                # with the values of the parameters


