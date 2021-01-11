import csv
import os
import random


def read_profile_into_dictionary(profile_path: str):
    """
    Gets a path to a text file of a miRNA profile and creates a dictionary that maps each miRNA to its expression.
    :param profile_path: the name of the text file
    :return: dictionary that maps each miRNA to its expression
    """
    # need to join the directories between the main dir and the file to filename
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

    return miRNA_counts


def convert_all_profiles(profiles_directory: str, documents_directory: str):
    """
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
    health_condition_dict = {}
    with open(path, 'r') as samples_file:
        samples_list = csv.reader(samples_file, delimiter="\t")
        next(samples_list, None)
        for sample in samples_list:
            file_name = sample[1]
            sample_type = 'Normal' if 'Normal' in sample[7] else 'Tumor'
            health_condition_dict[file_name] = sample_type
    return health_condition_dict


def create_dataset(profile_path: str, samples_file_path: str, organ_name: str):
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
                        write_to_csv_file(features_names, organ_name)
                        header = True
                    filename = os.path.split(file)[1]
                    sample_type = health_condition_dict[filename]
                    miRNA_counts["sample_type"] = sample_type
                    new_record = list(miRNA_counts.values())
                    write_to_csv_file(new_record, organ_name)
            except:
                continue


def write_to_csv_file(record, organ: str):
    with open(f'{organ}_samples.csv', 'a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(record)


# ------------------------------------------------------ Main ----------------------------------------------------------
profiles_directory = '/Users/royjudes/Desktop/miRNA embedding project/profiles'
documents_directory = '/Users/royjudes/Desktop/miRNA embedding project/documents'
# convert_all_profiles(profiles_directory, documents_directory)
# print(create_health_condition_dictionary(SAMPLES_FILE))
# PROFILES_DIRECTORY = "C:\\Users\\Naor\\Google Drive\\שנה ד'\\פרויקט גמר\\profiles\\Bronchus_and_lung"
# profiles_directory = ""
samples_file = "/Users/royjudes/Desktop/miRNA embedding project/gdc_sample_sheet.tsv"
# create_dataset(PROFILES_DIRECTORY, samples_file, "bronchus_and_lung")

# create_dataset(profiles_directory, samples_file, 'kidney')
# convert_all_profiles(profiles_directory, documents_directory)