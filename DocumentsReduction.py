import csv
import os
import random
import numpy as np
from Bio import SeqIO


def read_profile_into_dictionary(profile_path: str):
    """
    OLD VERSION - relevant to raw profiles of hairpin miRNAs.
    Gets a path to a text file of a miRNA profile and creates a dictionary that maps each miRNA to its expression.
    :param profile_path: the name of the text file
    :return: dictionary that maps each miRNA to its expression
    """
    miRNA_counts = {}
    first_line = True
    profile_file = open(profile_path, 'r')
    for miRNA_expression in profile_file:
        if first_line:
            first_line = False
            continue

        miRNA_expression = miRNA_expression.split()

        if np.ceil(float(miRNA_expression[2])) - float(miRNA_expression[2]) > 0.5:
            miRNA_counts[miRNA_expression[0]] = int(np.floor(float(miRNA_expression[2])))
        else:
            miRNA_counts[miRNA_expression[0]] = int(np.ceil(float(miRNA_expression[2])))

        # for unnormalized counts
        # miRNA_counts[miRNA_expression[0]] = int(miRNA_expression[1])

    return miRNA_counts


def convert_all_profiles(profiles_directory: str, documents_directory: str):
    """
    OLD VERSION - relevant to raw profiles of hairpin miRNAs.
    Reads all the text files of the miRNA profiles, then creates a document for each one. The document consists of
    the name of each miRNA the amount of times it appeared in the profile.
    """
    for (root, dirs, files) in os.walk(profiles_directory, topdown=True):
        for file in files:
            try:
                if file[-4:] == ".txt":
                    miRNA_counts = read_profile_into_dictionary(os.path.join(root, file))
                    filename = os.path.split(file)[1]
                    convert_profile_to_document(miRNA_counts, filename[:-4], documents_directory)
            except:
                continue


def read_mature_profile_into_dictionary(profile_path: str, mature_dict: dict):
    """
    NEW VERSION - relevant to formatted profiles of mature miRNAs.
    Gets a path to a csv file of a mature miRNA profile and creates a dictionary that maps each miRNA to its expression.

    :param mature_dict: dictionary that maps each 'MIMAT' to 'hsa'.
    :param profile_path: the name of the csv file.
    :return: dictionary that maps each miRNA to its expression.
    """
    miRNA_counts = {}
    first_line = True
    with open(profile_path, 'r') as profile:
        for miRNA_expression in profile:
            if first_line:
                first_line = False
                continue

            miRNA_expression = miRNA_expression.split(",")
            mature_MIMAT = miRNA_expression[0]
            miRNA_rpm = float(miRNA_expression[-2])

            if np.ceil(miRNA_rpm) - miRNA_rpm > 0.5:
                miRNA_counts[mature_dict[mature_MIMAT]] = int(np.floor(miRNA_rpm))
            else:
                miRNA_counts[mature_dict[mature_MIMAT]] = int(np.ceil(miRNA_rpm))

    return miRNA_counts


def convert_all_mature_profiles(profiles_directory: str, documents_directory: str, mature_dict: dict):
    """
    NEW VERSION - relevant to formatted profiles of mature miRNAs.
    Reads all the text files of the fixed formatted mature miRNA profiles, then creates a document for each one.
    The document consists of the name of each miRNA the amount of times it appeared in the profile.
    """
    for (root, dirs, files) in os.walk(profiles_directory, topdown=True):
        for file in files:
            try:
                if file[-4:] == ".csv":
                    miRNA_counts = read_mature_profile_into_dictionary(os.path.join(root, file), mature_dict)
                    filename = os.path.split(file)[1]
                    convert_profile_to_document(miRNA_counts, filename[:-4], documents_directory)
            except:
                continue


def create_mature_names_dictionary(fasta_file_path: str):
    """
    Creates a dictionary that maps the mature accession (MIMAT) to its name (hsa).

    :param fasta_file_path: path to the fasta file that contains all of the different mature miRNAs, with both their
    'MIMAT' and 'hsa' names.
    :return: a dictionary that maps the mature accession (MIMAT) to its name (hsa).
    """
    mature_dict = {}
    with open(fasta_file_path, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            description = record.description.split(" ")
            mature_MIMAT = description[1]
            mature_hsa = description[0]
            if "hsa" not in mature_hsa:
                continue
            mature_dict[mature_MIMAT] = mature_hsa

    return mature_dict


def convert_profile_to_document(miRNA_counts: dict, filename: str, documents_directory: str):
    """
    Creates a document for a given dictionary that maps a miRNA to its expression in a certain profile.

    :param miRNA_counts: the dictionary that maps a miRNA to its expression.
    :param filename: the name of the file of the profile (without extension).
    """
    path_to_file = os.path.join(documents_directory, f"{filename}_document.txt")
    document = open(path_to_file, 'w')
    miRNA_expression_as_list = convert_to_list_of_miRNA_expression(miRNA_counts)
    random.shuffle(miRNA_expression_as_list)
    for miRNA in miRNA_expression_as_list:
        document.write(miRNA + " ")


def convert_to_list_of_miRNA_expression(miRNA_counts: dict):
    """
    Gets a dictionary that maps each miRNA to its expression in a profile and converts it to a list that contains the
    miRNAs, each one the amount of times that corresponds to its expression.
    :param miRNA_counts: the dictionary of the miRNA expression
    :return: the list with the miRNAs
    """
    miRNA_expression = []
    for miRNA in miRNA_counts:
        for i in range(int(miRNA_counts[miRNA])):
            miRNA_expression.append(miRNA)

    return miRNA_expression


def create_health_condition_dictionary(path: str):
    """
    Creates a dictionary that maps each profile (by its filename) to the health condition of the profile
    (Normal or Tumor).

    :param path: path to the tsv file that contains the metadata of the profiles.
    :return: the dictionary with the profiles' health conditions.
    """
    health_condition_dict = {}
    with open(path, 'r') as samples_file:
        samples_list = csv.reader(samples_file, delimiter="\t")
        next(samples_list, None)
        for sample in samples_list:
            file_name = sample[1]
            sample_type = 'Normal' if 'Normal' in sample[7] else 'Tumor'
            project = sample[4]
            health_condition_dict[file_name] = sample_type, project
    return health_condition_dict


def create_dataset(profile_path: str, samples_file_path: str, organ_name: str):
    """
    Creates a csv file with the profiles as rows and miRNAs and health condition as columns. The file
    contains the normalized counts of each miRNA in each profile.

    :param profile_path: path to the directory of the profiles
    :param samples_file_path: path to the tsv file that contains the metadata of the profiles.
    :param organ_name: the name of the organ.
    """
    header = False
    health_condition_dict = create_health_condition_dictionary(samples_file_path)
    for (root, dirs, files) in os.walk(profile_path, topdown=True):
        for file in files:
            try:
                if file[-4:] == ".txt" and file != "MANIFEST.txt":
                    miRNA_counts = read_profile_into_dictionary(os.path.join(root, file))
                    if not header:
                        features_names = list(miRNA_counts.keys())
                        features_names.append("sample_type")
                        features_names.append("project")
                        write_to_csv_file(features_names, organ_name)
                        header = True
                    filename = os.path.split(file)[1]
                    sample_type, project = health_condition_dict[filename]
                    miRNA_counts["sample_type"] = sample_type
                    miRNA_counts["project"] = project
                    new_record = list(miRNA_counts.values())
                    write_to_csv_file(new_record, organ_name)
            except:
                continue


def write_to_csv_file(record, organ: str):
    with open(f'{organ}_samples_counts_project_names_health_condition.csv', 'a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(record)


# ------------------------------------------------------ Main ----------------------------------------------------------
# profiles_directory = '/Users/royjudes/Desktop/miRNA embedding project/profiles'
# documents_directory = '/Users/royjudes/Desktop/miRNA embedding project/documents'

profiles_directory = '/Users/royjudes/Desktop/miRNA embedding project/bronchus/profiles_formatted'
documents_directory = '/Users/royjudes/Desktop/miRNA embedding project/bronchus/documents'
fasta_file_path = "/Users/royjudes/Desktop/miRNA embedding project/a/mature.fa"

mature_dict = create_mature_names_dictionary(fasta_file_path)
convert_all_mature_profiles(profiles_directory, documents_directory, mature_dict)


# samples_file = "/Users/royjudes/Desktop/miRNA embedding project/gdc_sample_sheet.tsv"
# create_dataset(profiles_directory, samples_file, "kidney")

# create_dataset(profiles_directory, samples_file, 'kidney')
# convert_all_profiles(profiles_directory, documents_directory)



# pipeline:
# putting the formatted profiles in the folder in {profiles_directory}, the path to the folder of the documents in
# {documents_directory} and keep the fasta file as is. then, run create_matures_names_dictionary function in order to
# get the dictionary which maps the mirs from MIMAT to miR versions, then run convert_all_mature_profiles function
# to get the documents
