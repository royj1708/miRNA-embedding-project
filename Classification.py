import sklearn as sk
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import pandas as pd
import numpy as np


kidney_samples_csv_path = "C:\\Users\\Naor\\Google Drive\\שנה ד'\\פרויקט גמר\\profiles\\Bronchus_and_lung\\bronchus_and_lung_samples.csv"
embedded_path = "C:\\Users\\Naor\\Google Drive\\שנה ד'\\פרויקט גמר\\profiles\Kidney\\embedded_profiles.csv"
features = pd.read_csv(embedded_path)

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

print('Training Features Shape:', train_features.shape)
print('Training Labels Shape:', train_labels.shape)
print('Testing Features Shape:', test_features.shape)
print('Testing Labels Shape:', test_labels.shape)

random_forest = RandomForestClassifier(n_estimators=100)
random_forest.fit(train_features, train_labels)
predictions = random_forest.predict(test_features)
accuracy = metrics.accuracy_score(test_labels, predictions)
print(f"accuracy: {accuracy}")
recall = metrics.recall_score(test_labels, predictions, pos_label='Normal')
print(f"recall: {recall}")
recall = metrics.recall_score(test_labels, predictions, pos_label='Tumor')
print(f"recall: {recall}")

print(len(predictions))
counter = 0
for i in predictions:
    if i == 'Normal':
        counter += 1

print(counter)
