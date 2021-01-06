from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import pandas as pd
import numpy as np


def oversampling(features):
    """
    fills the data set with more 'Normal' labelled profiles in order to balance it.
    """
    # Class count
    count_class_tumor, count_class_normal = features.sample_type.value_counts()
    # Divide by class
    features_normal = features[features['sample_type'] == "Normal"]
    features_tumor = features[features['sample_type'] == "Tumor"]

    features_normal_over = features_normal.sample(count_class_tumor, replace=True)
    return pd.concat([features_tumor, features_normal_over], axis=0)

def get_random_forest_predictions(train_features, test_features, train_labels, n_estimators):
    random_forest = RandomForestClassifier(n_estimators=n_estimators)
    random_forest.fit(train_features, train_labels)
    predictions = random_forest.predict(test_features)
    return predictions


def get_svc_predictions(train_features, test_features, train_labels):
    svc = SVC(kernel='linear')
    svc.fit(train_features, train_labels)
    predictions = svc.predict(test_features)
    return predictions


def get_naive_bayes_predictions(train_features, test_features, train_labels):
    naive_bayes = GaussianNB()
    naive_bayes.fit(train_features, train_labels)
    predictions = naive_bayes.predict(test_features)
    return predictions

def get_metrics_scores(test_labels, predictions):
    """
    gets the labels of the test set (y_true) and a classification model's predictions (y_pred) and computes different
    evaluation metrics accordingly.
    the return value is a dictionary with each metric and its score.
    """
    accuracy = metrics.accuracy_score(test_labels, predictions)

    results = {'recall_normal': metrics.recall_score(test_labels, predictions, pos_label='Normal'),
               'recall_tumor': metrics.recall_score(test_labels, predictions, pos_label='Tumor'), 'accuracy': accuracy}

    test_labels_auc = []
    predictions_auc = []
    for i in range(len(test_labels)):
        if test_labels[i] == 'Normal':
            test_labels_auc.append(0)
        else:
            test_labels_auc.append(1)

        if predictions[i] == 'Normal':
            predictions_auc.append(0)
        else:
            predictions_auc.append(1)

    results['auc'] = metrics.roc_auc_score(test_labels_auc, predictions_auc)
    return results

def main(samples_csv_path):
    features = pd.read_csv(samples_csv_path)
    #oversampling - optional
    features = oversampling(features)
    # Labels are the values we want to predict
    labels = np.array(features['sample_type'])
    # Remove the labels from the features
    # axis 1 refers to the columns
    features= features.drop('sample_type', axis=1)
    # Saving feature names for later use
    feature_list = list(features.columns)
    # Convert to numpy array
    features = np.array(features)

    # Split the data into training and testing sets
    train_features, test_features, train_labels, test_labels = train_test_split(features, labels, test_size = 0.3,stratify=labels)

    # print('Training Features Shape:', train_features.shape)
    # print('Training Labels Shape:', train_labels.shape)
    # print('Testing Features Shape:', test_features.shape)
    # print('Testing Labels Shape:', test_labels.shape)
    results_random_forest = get_metrics_scores(test_labels, get_random_forest_predictions(train_features, test_features, train_labels, 100))
    print("--Random Forest model--")
    print(f"accuracy: {results_random_forest['accuracy']}")
    print(f"recall normal: {results_random_forest['recall_normal']}")
    print(f"recall tumor: {results_random_forest['recall_tumor']}")
    print(f"auc: {results_random_forest['auc']}", "\n")

    results_svc = get_metrics_scores(test_labels, get_svc_predictions(train_features, test_features, train_labels))
    print("--SVC model--")
    print(f"accuracy: {results_svc['accuracy']}")
    print(f"recall normal: {results_svc['recall_normal']}")
    print(f"recall tumor: {results_svc['recall_tumor']}")
    print(f"auc: {results_svc['auc']}", "\n")

    results_naive_bayes = get_metrics_scores(test_labels, get_naive_bayes_predictions(train_features, test_features, train_labels))
    print("--Naive Bayes model--")
    print(f"accuracy: {results_naive_bayes['accuracy']}")
    print(f"recall normal: {results_naive_bayes['recall_normal']}")
    print(f"recall tumor: {results_naive_bayes['recall_tumor']}")
    print(f"auc: {results_naive_bayes['auc']}")


kidney_samples_csv_path = "C:\\Users\\Naor\\Google Drive\\שנה ד'\\פרויקט גמר\\profiles\\Bronchus_and_lung\\bronchus_and_lung_samples.csv"
kidney_embedded_samples_csv_path = "C:\\Users\\Naor\\Google Drive\\שנה ד'\\פרויקט גמר\\profiles\Kidney\\embedded_profiles.csv"

# kidney_embedded_samples_csv_path = 'embedded_profiles_0_300_5.csv'
# kidney_samples_csv_path = '/Users/royjudes/Desktop/miRNA embedding project/kidney_samples.csv'

print("<<<<<<<<<<< RESULTS WITHOUT EMBEDDING: >>>>>>>>>>>")
main(kidney_samples_csv_path)
print()
print()
print("<<<<<<<<<<< RESULTS WITH EMBEDDING: >>>>>>>>>>>")
main(kidney_embedded_samples_csv_path)