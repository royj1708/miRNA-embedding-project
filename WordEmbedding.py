import os
from gensim.models.word2vec import Word2Vec
from multiprocessing import cpu_count
import numpy as np
import xlrd
from datetime import datetime
import pandas as pd
import csv
from csv import reader
from numpy import dot
from numpy.linalg import norm


DOCUMENTS_DIRECTORY = '/Users/royjudes/Desktop/miRNA embedding project/breast/documents'


def create_list_from_document(path: str):
    """
    Reads a document into a list
    """
    with open(path, 'r') as document:
        document_as_list = document.read().split()

    return document_as_list


def convert_documents_to_lists():
    """
    Reads all the documents in the directory and converts each one to a list.
    :return: list of all the documents as lists
    """
    list_of_all_documents = []
    for (root, dirs, files) in os.walk(DOCUMENTS_DIRECTORY, topdown=True):
        for file in files:
            try:
                if file[-4:] == ".txt":
                    list_of_all_documents.append(create_list_from_document(os.path.join(DOCUMENTS_DIRECTORY, file)))

            except:
                continue

    return list_of_all_documents


# def convert_documents_csv_to_lists(path_to_csv):
#     """
#
#     :param path_to_csv:
#     :return:
#     """
#     list_of_all_documents = []
#     samples_df = pd.read_csv(path_to_csv)
#     samples_df = samples_df.drop('sample_type', axis=1)
#     miRNAs = list(samples_df.columns)
#
#     for index, profile in samples_df.iterrows():
#         profile_vector = []
#         for miRNA in miRNAs:
#             profile_vector.append(profile[miRNA])
#
#         list_of_all_documents.append(profile_vector)
#
#     return list_of_all_documents


def convert_documents_to_dictionary():
    """
    Converts all the documents to lists and stores them in a dictionary. The key is the document's filename and the
    value is the list form of the document.
    """
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
    """
    Given a language model, computes Word Mover Distance of each profile pair and writes it into a csv.
    :param profiles_as_list: all the profiles as lists
    :param keys: the filename of the documents (the profiles)
    :param health_condition_dictionary: a dictionary that maps between a profile and its health condition
    :param model: the language model
    """
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


# CREATING PROFILES WITH EMBEDDING REPRESENTATION - created for the classification step,
# in order to test the classification of profiles after the embedding step
def create_embedded_profiles(embedded_profiles_path: str, samples_file: str, model):
    """
    Creates a csv dataset in which every record is a profile, computed with the miRNA embeddings.
    :param embedded_profiles_path: the path of the file in which the profiles are written
    :param samples_file: path to the dataset with the raw profiles and the miRNA counts
    :param model: the language model
    """
    samples_df = pd.read_csv(samples_file)
    sample_types = samples_df["sample_type"]
    samples_df = samples_df.drop('sample_type', axis=1)
    miRNAs = list(samples_df.columns)

    with open(embedded_profiles_path, 'a', newline='') as file:
        writer = csv.writer(file)
        for index, profile in samples_df.iterrows():
            miRNAs_amount_in_profile = 0
            profile_vector = []
            for i in range(300):
                profile_vector.append(0)
            for miRNA in miRNAs:
                try:
                    miRNA_vector = model[miRNA]  # the embedding vector of the miRNA
                    expression = float(profile[miRNA])  # the miRNA expression in the profile
                    miRNAs_amount_in_profile += expression  # sum of all miRNAs in the profile
                    profile_vector = np.add(profile_vector, np.array(miRNA_vector)*expression)

                except:
                    continue

            profile_vector /= miRNAs_amount_in_profile
            profile_vector = list(profile_vector)
            profile_vector.append(sample_types[index])

            writer.writerow(profile_vector)


# CREATING CSV WITH MIRNA EMBEDDINGS
def create_miRNA_embeddings_file(embeddings_file: str, samples_file: str, model):
    """
    Creates a csv file with the miRNA and their embedding vectors.
    :param embeddings_file: the file in which to write
    :param samples_file: path to the dataset with the raw profiles and the miRNA counts
    :param model: the language model
    """
    samples_df = pd.read_csv(samples_file)
    samples_df = samples_df.drop('sample_type', axis=1)
    miRNAs = list(samples_df.columns)

    with open(embeddings_file, 'a', newline='') as file:
        writer = csv.writer(file)
        for miRNA in miRNAs:
            try:
                writer.writerow([miRNA])
                writer.writerow(model[miRNA])

            except:
                continue


def load_miRNA_embeddings_dictionary(path: str):
    """
    Reads the miRNA embeddings dataset and stores it into a dictionary, where each key is the miRNA and the value is the
    embedding vector.
    :param path: a path to the embeddings file
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


def compute_cossim_for_miRNAs(embeddings_dictionary, counts_dict):
    """
    Creates a dataset in which every record is a miRNA pair and the cosine similarity score, computed on their embedding
    vectors.
    :param embeddings_dictionary: dictionary with miRNA as key and its embedding vector as value
    :param counts_dict: a dictionary with miRNA as key and its average count as value.
    """
    with open('mature_miRNAs_cossim_breast.csv', 'a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['miRNA_a', 'avg_rpm_a', 'miRNA_b', 'avg_rpm_b', 'cos_sim'])
        already_computed = []
        for miRNA_a in embeddings_dictionary:
            already_computed.append(miRNA_a)
            for miRNA_b in embeddings_dictionary:
                if miRNA_b in already_computed:
                    continue

                mir_a_vec = embeddings_dictionary[miRNA_a]
                mir_b_vec = embeddings_dictionary[miRNA_b]
                cos_sim = dot(mir_a_vec, mir_b_vec) / (norm(mir_a_vec) * norm(mir_b_vec))

                writer.writerow([miRNA_a, counts_dict[miRNA_a], miRNA_b,counts_dict[miRNA_b], cos_sim])


def load_avg_mirna_counts_dict(path: str):
    """
    Reads a file with each miRNA and its average count and stores it in a dictionary.
    :param path: the path to the file
    """
    df = pd.read_csv(path)
    df = df.transpose()
    counts_dict = {}
    for mir, row in df.iterrows():
        counts_dict[mir] = row[0]

    return counts_dict


# ------------------------------------------------------ Main ----------------------------------------------------------
def main():
    all_profiles = convert_documents_to_dictionary()
    # health_condition_dictionary = create_health_condition_dictionary(
    #     '/Users/royjudes/Desktop/miRNA embedding project/gdc_sample_sheet.tsv')
    profiles_as_list = list(all_profiles.values())
    # keys = list(all_profiles.keys())
    model = Word2Vec(profiles_as_list, min_count=0, size=300, window=5, workers=cpu_count())
    # compute_wmd_for_profiles(profiles_as_list, keys, health_condition_dictionary, model)

    samples_file = "/Users/royjudes/Desktop/miRNA embedding project/breast/profiles_matures_dataset.csv"
    # embedded_profiles_path = 'embedded_mature_profiles_0_300_5.csv'
    # create_embedded_profiles(embedded_profiles_path, samples_file, model)

    embeddings_file_name = 'mature_miRNA_embeddings_breast.csv'
    create_miRNA_embeddings_file(embeddings_file_name, samples_file, model)


