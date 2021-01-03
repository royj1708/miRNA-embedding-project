from gensim import corpora
from pprint import pprint
from gensim.utils import simple_preprocess
from smart_open import smart_open
import os
from gensim.models.word2vec import Word2Vec
from multiprocessing import cpu_count
import gensim.downloader as api
import numpy as np
import xlrd
from datetime import datetime

DOCUMENTS_DIRECTORY = "C:\\Users\\royj1\\Desktop\\University\\הנדסת מערכות מידע\\שנה ד\\פרויקט גמר\\דאטה\\Testing\\b"
VOCABULARY_FILE = "C:\\Users\\royj1\\Desktop\\University\\הנדסת מערכות מידע\\שנה ד\\פרויקט גמר\\דאטה\\all_mirnas.xls"
LOG_FILE_PATH = "C:\\Users\\royj1\\Desktop\\University\\הנדסת מערכות מידע\\שנה ד\\פרויקט גמר\\דאטה\\Testing\\Log.txt"


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


def get_classification_accuracy():
    pass


def write_into_log_file(size: int, window: int, min_count: int, model: str, results):
    content = f"Time: {datetime.now().strftime('%d/%m/%Y %H:%M:%S')}\n" \
              f"Results - {model}:\n" \
              f"Embedding Vector Size: {size}, Context Window Size: {window}, Minimal amount of appearances " \
              f"to be included: {min_count}, Accuracy: {results}\n\n"

    with open(LOG_FILE_PATH, 'a') as log_file:
        log_file.write(content)


def main():
    all_profiles = convert_documents_to_lists()
    model = Word2Vec(all_profiles, min_count=0, size=300, window=5, workers=cpu_count())
    distance = model.wmdistance(all_profiles[0], all_profiles[1])
    print(distance)


main()
    # for size in [100, 200, 300, 400, 500]:
    #     for window in [3, 4, 5, 6, 7]:
    #         for min_count in [0, 5, 10, 15, 20]:
    #             word2vec_model = Word2Vec(all_profiles, min_count=min_count, window=window, size=size, workers=cpu_count())
                # print(word2vec_model.most_similar(positive=['hsa-let-7a-1', 'hsa-let-7a-2'], negative=['hsa-let-7a-3'], topn=1))

                # write_into_log_file(size, window, min_count, None)
                # we need to take all the miRNA vectors and try the classification models, then record the accuracy rate
                # with the values of the parameters

# dictionary = corpora.Dictionary(simple_preprocess(line, deacc=True) for line in open('doc.txt', encoding='utf-8'))

# # Download dataset
# data = []
# dataset = api.load("text8")
# data.extend([d for d in dataset])
#
# # Split the data into 2 parts. Part 2 will be used later to update the model
# data_part1 = data[:1000]
# data_part2 = data[1000:]
#
# # Train Word2Vec model. Defaults result vector size = 100
# model = Word2Vec(data, min_count=0, workers=cpu_count())
# result = model.most_similar(positive=['woman', 'king'], negative=['man'], topn=3)
# print(result)
#
#
#
# model.build_vocab(data_part2, update=True)
# model.train(data_part2, total_examples=model.corpus_count, epochs=model.iter)

