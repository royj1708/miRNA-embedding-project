import csv
import os
import random
import numpy as np
from Bio import SeqIO


def convert_all_profiles(profiles_directory: str, target_path: str):
    """
    Reads all the text files of the mature miRNA profiles, then creates a fixed formatted file for each one. The document consists of
    the name of each mature miRNA, the hairpin miRNAs it is associated with, its read counts and normalized read counts and whether
    it is a guide or star.

    :param profiles_directory: path to the directory of the profiles.
    :param target_path: path to the directory where the new profiles will be created.
    """
    for (root, dirs, files) in os.walk(profiles_directory, topdown=True):
        for file in files:
            try:
                if file[-4:] == ".txt":
                    filename = os.path.split(file)[1]
                    filename = f"{filename[:-4]}.csv"
                    create_mature_profile_from_isoforms(os.path.join(root, file), target_path, filename)

                    # miRNA_counts = read_profile_into_dictionary(os.path.join(root, file))
                    # convert_profile_to_document(miRNA_counts, filename[:-4], target_path)
            except:
                continue


def create_mature_profile_from_isoforms(profile_path: str, target_path: str, filename: str):
    """
    Reads an isoform profile and converts it to a profile with mature miRNA counts. Each record in the new profile has a
    mature miRNA as key and the rest of the record is the read counts and normalized counts of the mature and a list of
    hairpin miRNAs that is associated with the mature.

    :param profile_path: the path to the isoform profile.
    """
    first_line = True

    # the data structure that is built while reading the file is a dictionary that contains mature name as key and
    # another dictionary as value, which contains hairpin miRNA as key and a list of its counts as value.
    mature_to_hairpins_dict = {}
    # dictionary for mapping hairpins to their matures, in order to determine which mature is the guide
    hairpins_to_matures_dict = {}

    profile_file = open(profile_path, 'r')
    for record in profile_file:
        if first_line:
            first_line = False
            continue

        record = record.split()
        # filters out precursors and unnecessary types
        mature_name = record[5]
        if "MIMAT" not in mature_name:
            continue
        else:
            mature_name = mature_name.replace("mature,", "")

        read_counts = int(record[2])
        read_per_million = float(record[3])
        hairpin = record[0]

        # checks if this is the first time this mature is encountered
        if mature_name in mature_to_hairpins_dict:
            # checks if this is the first time that this mature has been associated with this hairpin
            if hairpin in mature_to_hairpins_dict[mature_name]:
                # if not then we sum the reads with the existing amount
                mature_to_hairpins_dict[mature_name][hairpin][0] += read_counts
                mature_to_hairpins_dict[mature_name][hairpin][1] += read_per_million

            else:
                # if yes then we add an entry for this mature with the hairpin and its counts
                mature_to_hairpins_dict[mature_name][hairpin] = [read_counts, read_per_million]

        else:
            # if yes then we add a new entry with the mature and a dictionary with the hairpin and its counts
            mature_to_hairpins_dict[mature_name] = {hairpin: [read_counts, read_per_million]}

        # does the same as the upper if-else block, but for mapping the hairpins to the matures. it only counts the
        # read counts this time, since it is only to determine which mature has the highest total count
        if hairpin in hairpins_to_matures_dict:
            if mature_name in hairpins_to_matures_dict[hairpin]:
                hairpins_to_matures_dict[hairpin][mature_name] += read_counts
            else:
                hairpins_to_matures_dict[hairpin][mature_name] = read_counts

        else:
            hairpins_to_matures_dict[hairpin] = {mature_name: read_counts}

    # creates a csv file with matures as keys and their hairpins and counts. for each mature-hairpin combination, it is
    # also written if the mature is the guide or star
    with open(os.path.join(target_path, filename), 'a') as profile:
        writer = csv.writer(profile)
        writer.writerow(["mature_name", "miRNA_IDs", "read_count", "read_per_million", "guide/star"])
        for mature in mature_to_hairpins_dict:
            # for each mature, finds the hairpins it is associated with, and then finds its maximum count and in which
            # hairpin.
            hairpins = list(mature_to_hairpins_dict[mature].keys())
            max_count = max(count for count in list(mature_to_hairpins_dict[mature].values()))
            max_index = list(mature_to_hairpins_dict[mature].values()).index(max_count)
            max_hairpin = list(mature_to_hairpins_dict[mature].keys())[max_index]

            # finds the mature of the hairpin we found in order to determine if the mature we are now checking is its
            # guide or star
            max_count_for_hairpin = max(count for count in list(hairpins_to_matures_dict[max_hairpin].values()))
            max_index_hairpin = list(hairpins_to_matures_dict[max_hairpin].values()).index(max_count_for_hairpin)
            max_mature = list(hairpins_to_matures_dict[max_hairpin].keys())[max_index_hairpin]

            if mature == max_mature:
                guide = True
            else:
                guide = False

            record = [mature, hairpins, max_count[0], max_count[1], "GUIDE" if guide else "STAR"]
            writer.writerow(record)


def read_profile_into_dictionary(profile_path: str, mature_dict: dict):
    """
    Gets a path to a text file of a mature miRNA fixed formatted profile and creates a dictionary that maps each
    'hsa' mature miRNA to its expression.

    :param profile_path: the name of the text file
    :param mature_dict: dictionary that maps each 'MIMAT' to 'hsa'
    :return: dictionary that maps each miRNA to its expression
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
            miRNA_counts[mature_dict[mature_MIMAT]] = float(miRNA_expression[-2])

    return miRNA_counts


def create_health_condition_dictionary(samples_file_path: str):
    """
    Creates a dictionary that maps each profile (by its filename) to the health condition of the profile
    (Normal or Tumor).

    :param samples_file_path: path to the tsv file that contains the metadata of the profiles.
    :return: the dictionary with the profiles' health conditions.
    """
    health_condition_dict = {}
    with open(samples_file_path, 'r') as samples_file:
        samples_list = csv.reader(samples_file, delimiter="\t")
        next(samples_list, None)
        for sample in samples_list:
            file_name = sample[1]
            sample_type = 'Normal' if 'Normal' in sample[7] else 'Tumor'
            project = sample[4]
            health_condition_dict[file_name] = sample_type, project
    return health_condition_dict


def create_mature_dataset(fasta_file_path: str, profiles_path: str, samples_file: str):
    """
    Creates a csv file with the profiles as rows and 'hsa' mature miRNAs and health condition as columns. The file
    contains the normalized counts of each miRNA in each profile.

    :param fasta_file_path: path to the fasta file that contains all of the different mature miRNAs, with both their
    'MIMAT' and 'hsa' names.
    :param profiles_path: path to the directory of the profiles.
    :param samples_file: path to the tsv file that contains the metadata of the profiles.
    """
    # creates a dictionary that maps the mature accession (MIMAT) to its name (hsa)
    mature_dict = {}
    with open(fasta_file_path, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            description = record.description.split(" ")
            mature_MIMAT = description[1]
            mature_hsa = description[0]
            if "hsa" not in mature_hsa:
                continue
            mature_dict[mature_MIMAT] = mature_hsa

    # creates a dataset with the profiles as rows and the matures miRNAs as columns
    mature_hsa_list = list(mature_dict.values())
    feature_names = list(mature_dict.values())
    feature_names.append("sample_type")

    profiles_health_condition_dict = create_health_condition_dictionary(samples_file)
    with open("profiles_matures_dataset.csv", 'a') as dataset:
        writer = csv.writer(dataset)
        writer.writerow(feature_names)

        # iterates over the profiles and fills the dataset with the counts of each mature
        for (root, dirs, files) in os.walk(profiles_path, topdown=True):
            for file in files:
                try:
                    if file[-4:] == ".csv":
                        with open(os.path.join(root, file), 'r') as profile:
                            # checks the MIMAT mature and converts it to the hsa name
                            profile_dict = read_profile_into_dictionary(os.path.join(root, file), mature_dict)
                            profile_counts_list = []
                            for hsa in mature_hsa_list:
                                if hsa in profile_dict:
                                    profile_counts_list.append(profile_dict[hsa])
                                else:
                                    profile_counts_list.append(0)

                            filename = os.path.split(file)[1].replace(".csv", ".txt")
                            profile_counts_list.append(profiles_health_condition_dict[filename][0])
                            writer.writerow(profile_counts_list)

                except:
                    continue


profiles_path_txt = "/Users/royjudes/Desktop/miRNA embedding project/bronchus/profiles_raw"
profiles_path_csv = "/Users/royjudes/Desktop/miRNA embedding project/bronchus/profiles_formatted"
convert_all_profiles(profiles_path_txt, profiles_path_csv)

samples_file = "/Users/royjudes/Desktop/miRNA embedding project/bronchus/gdc_sample_sheet_bronchus.tsv"
fasta_file_path = "/Users/royjudes/Desktop/miRNA embedding project/a/mature.fa"
create_mature_dataset(fasta_file_path, profiles_path_csv, samples_file)


# pipeline:
# put the path to the raw profiles folder in {profiles_path_txt} and the target folder to the formatted profiles in
# {profiles_path_csv} and run convert_all_profiles function.
# then, put in {samples_file} the path to the metadata about the profiles. {fasta_file_path} should stay the same.
# then, run the create_matures_dataset function with the parameters, and you get the csv dataset of mirna as columns
# and profiles as rows.
