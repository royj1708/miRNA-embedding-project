import csv
from csv import writer


def merge_cossim_datasets():
    """
    Merges the datasets of miRNA pairs and the cosine similarity score in each organ.
    """
    bronchus_dict = {}
    kidney_dict = {}
    breast_dict = {}
    first_line = True
    with open("/Users/royjudes/Desktop/miRNA embedding project/bronchus/mature_miRNAs_cossim_bronchus.csv", 'r') as bronchus:
        for record in bronchus:
            if first_line:
                first_line = False
                continue

            record = record.split(",")
            bronchus_dict[(record[0], record[2])] = record[-1].replace("\n", "")

    first_line = True
    with open("/Users/royjudes/Desktop/miRNA embedding project/mature_miRNAs_cossim_kidney.csv", 'r') as kidney:
        for record in kidney:
            if first_line:
                first_line = False
                continue

            record = record.split(",")
            kidney_dict[(record[0], record[2])] = record[-1].replace("\n", "")

    first_line = True
    with open("/Users/royjudes/Desktop/miRNA embedding project/breast/mature_miRNAs_cossim_breast.csv", 'r') as breast:
        for record in breast:
            if first_line:
                first_line = False
                continue

            record = record.split(",")
            breast_dict[(record[0], record[2])] = record[-1].replace("\n", "")

    # with open('matures_miRNAs_cossim_kidney_bronchus_breast.csv', 'a') as merged:
    #     writer = csv.writer(merged)
    #     writer.writerow(['miRNA_a', 'miRNA_b', 'cos_sim_kidney', 'cos_sim_bronchus', 'cos_sim_breast'])
    #     for mir_tuple in kidney_dict:
    #         if mir_tuple in bronchus_dict and mir_tuple in breast_dict:
    #             writer.writerow([mir_tuple[0], mir_tuple[1], kidney_dict[mir_tuple],
    #                              bronchus_dict[mir_tuple], breast_dict[mir_tuple]])
    print("number of miRNAs:")
    print()
    print(f"kidney: {len(kidney_dict)}")
    print(f"bronchus: {len(bronchus_dict)}")
    print(f"breast: {len(breast_dict)}")
    overlap = len([pair for pair in kidney_dict if pair in bronchus_dict and pair in breast_dict])
    print(f"overlap: {overlap}")


