from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import pandas as pd
import numpy as np
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import NearMiss
from collections import Counter


def oversampling(features):
    """
    fills the data set with more 'Normal' labelled profiles in order to balance it.
    """
    # Class count
    count_class_tumor, count_class_normal = features.sample_type.value_counts()
    # Divide by class
    features_normal = features[featur['sample_type'] == "Normal"]
    features_tumor = features[features['sample_type'] == "Tumor"]

    features_normal_over = features_normal.sample(count_class_tumor, replace=True)
    return pd.concat([features_tumor, features_normal_over], axis=0)


def get_classifier_predictions(classifier, train_features, test_features, train_labels):
    classifier.fit(train_features, train_labels)
    predictions = classifier.predict(test_features)
    return predictions


def get_feature_importance_list(classifier):
    return list(classifier.feature_importances_)


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
    # oversampling - optional
    # features = oversampling(features)
    # Labels are the values we want to predict
    labels = np.array(features['sample_type'])
    print(Counter(labels))
    # Remove the labels from the features
    # axis 1 refers to the columns
    features = features.drop('sample_type', axis=1)
    # Saving feature names for later use
    feature_list = list(features.columns)
    # Convert to numpy array
    # features = np.array(features)

    # SMOTE
    # features_smote, labels_smote = smote(features, labels)
    # print(Counter(labels_smote))

    # Split the data into training and testing sets
    train_features, test_features, train_labels, test_labels = train_test_split(features, labels, test_size=0.3) #, stratify=labels)

    # SMOTE
    smt = SMOTE()
    train_features, train_labels = smt.fit_sample(train_features, train_labels)
    print(Counter(train_labels))


    # print('Training Features Shape:', train_features.shape)
    # print('Training Labels Shape:', train_labels.shape)
    # print('Testing Features Shape:', test_features.shape)
    # print('Testing Labels Shape:', test_labels.shape)

    # -------------------------------------------- Random Forest --------------------------------------------
    random_forest = RandomForestClassifier(n_estimators=100)
    random_forest_predictions = get_classifier_predictions(random_forest, train_features, test_features, train_labels)
    results_random_forest = get_metrics_scores(test_labels, random_forest_predictions)
    random_forest_feature_importance = get_feature_importance_list(random_forest)
    feature_importance_dict = {}
    for index, feature in enumerate(feature_list):
        feature_importance_dict[feature] = random_forest_feature_importance[index]

    print("--Random Forest model--")
    print(f"accuracy: {results_random_forest['accuracy']}")
    print(f"recall normal: {results_random_forest['recall_normal']}")
    print(f"recall tumor: {results_random_forest['recall_tumor']}")
    print(f"auc: {results_random_forest['auc']}", "\n")
    print("Feature Importance:\n", random_forest_feature_importance, "\n")
    for index, feature in enumerate(sorted(feature_importance_dict, key=feature_importance_dict.get, reverse=True)):
        print(f"MIR: {feature}, Importance: {feature_importance_dict[feature]}")
        if index == 100:
            break

    # -------------------------------------------- SVC --------------------------------------------
    svc = SVC(kernel='linear')
    svc_predictions = get_classifier_predictions(svc, train_features, test_features, train_labels)
    results_svc = get_metrics_scores(test_labels, svc_predictions)
    print("--SVC model--")
    print(f"accuracy: {results_svc['accuracy']}")
    print(f"recall normal: {results_svc['recall_normal']}")
    print(f"recall tumor: {results_svc['recall_tumor']}")
    print(f"auc: {results_svc['auc']}", "\n")

    # -------------------------------------------- Naive Bayes --------------------------------------------
    naive_bayes = GaussianNB()
    naive_bayes_predictions = get_classifier_predictions(naive_bayes, train_features, test_features, train_labels)
    results_naive_bayes = get_metrics_scores(test_labels, naive_bayes_predictions)
    print("--Naive Bayes model--")
    print(f"accuracy: {results_naive_bayes['accuracy']}")
    print(f"recall normal: {results_naive_bayes['recall_normal']}")
    print(f"recall tumor: {results_naive_bayes['recall_tumor']}")
    print(f"auc: {results_naive_bayes['auc']}")


# kidney_samples_csv_path = "C:\\Users\\Naor\\Google Drive\\שנה ד'\\פרויקט גמר\\profiles\\Bronchus_and_lung\\bronchus_and_lung_samples.csv"
# kidney_embedded_samples_csv_path = "C:\\Users\\Naor\\Google Drive\\שנה ד'\\פרויקט גמר\\profiles\Kidney\\embedded_profiles.csv"

# kidney_embedded_samples_csv_path = 'embedded_profiles_0_300_5.csv'
# kidney_embedded_samples_csv_path = '/Users/royjudes/Desktop/miRNA embedding project/embedding configurations/0_300_2/embedded_profiles_0_300_2.csv'
# kidney_embedded_samples_csv_path = '/Users/royjudes/Desktop/miRNA embedding project/embedding configurations/0_300_15/embedded_profiles_0_300_15.csv'
kidney_embedded_samples_csv_path = '/Users/royjudes/Desktop/miRNA embedding project/embedded_mature_profiles_0_300_5.csv'
kidney_samples_csv_path = '/Users/royjudes/Desktop/miRNA embedding project/profiles_matures_dataset.csv'

print("<<<<<<<<<<< RESULTS WITHOUT EMBEDDING: >>>>>>>>>>>")
main(kidney_samples_csv_path)
print()
print()
print("<<<<<<<<<<< RESULTS WITH EMBEDDING: >>>>>>>>>>>")
main(kidney_embedded_samples_csv_path)

