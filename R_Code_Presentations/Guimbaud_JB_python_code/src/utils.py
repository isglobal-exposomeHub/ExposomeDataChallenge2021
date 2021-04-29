# Author Jean-Baptiste Guimbaud
# Meersens

import pandas as pd
from enum import Enum
import numpy as np
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error

class Phenotypes(Enum):
    BIRTH_WEIGHT = "e3_bw"
    DIAGNOSED_ASTHMA = "hs_asthma"
    BODY_MASS_INDEX = "hs_zbmi_who"
    IQ = "hs_correct_raven"
    NEURO_BEHAVIOR = "hs_Gen_Tot"
    BODY_MASS_INDEX_CATEGORICAL = "hs_bmi_c_cat"

class Models(Enum):
    NAIVE = 0
    ELASTIC = 1
    RF = 2
    XGB = 3
    MLP = 4
    SVM = 5
    CART = 6

class DataType(Enum):
    POSTNATAL = 1
    PREGNANCY = 2
    BOTH = 3

def make_features(exposome, covariates, phenotype):
    print("joining dataframes on same ID")

    joined_ec = pd.merge(exposome, covariates, on="ID", how="inner")
    joined_ecp = pd.merge(joined_ec, phenotype, on="ID", how="inner")

    print("dataframe of:")
    print("  *", joined_ecp.shape[0], "rows total")
    print("  *", joined_ecp.shape[1], "features total")
    return joined_ecp
    
    
def standardize_data(train_features, train_labels, test_features=None, test_labels=None):
    features_mean = train_features.mean()
    features_std = train_features.std()
    train_features = (train_features - features_mean) / features_std
    test_features = (test_features - features_mean) / features_std


    labels_mean, labels_std = None, None
    if target_name != Phenotypes.BODY_MASS_INDEX_CATEGORICAL and target_name != Phenotypes.DIAGNOSED_ASTHMA:
        labels_mean = train_labels.mean()
        labels_std = train_labels.std()
        train_labels = (train_labels - labels_mean) / labels_std
        test_labels = (test_labels - labels_mean) / labels_std

    return {
        "train_features": train_features,
        "train_labels": train_labels,
        "test_features": test_features,
        "test_labels": test_labels,
        "features_mean": features_mean,
        "features_std": features_std,
        "labels_mean": labels_mean,
        "labels_std": labels_std
    }

# Removes all features from a dataframe that are correlated to a given predictor.
def find_features_correlated_to(feature_name, threshold, dataframe):
    col_corr = set()
    corr_matrix = dataframe.corr()
    for i in range(len(corr_matrix.columns)):
        if abs(corr_matrix[feature_name].iloc[i]) > threshold:
            colname = corr_matrix.columns[i]
            col_corr.add(colname)
    col_corr.remove(feature_name)
    return col_corr

def compute_residuals(features, labels, classification=False):
    feature_scaler, label_scaler = StandardScaler(), StandardScaler()
    features = feature_scaler.fit_transform(features)
    labels = label_scaler.fit_transform(labels)

    model = LinearRegression(normalize=False)
    if classification:
        model = LogisticRegression(penalty=None, class_weight='balanced', random_state=42)

    model.fit(features, labels)
    predictions = model.predict(features)
    
    predictions = label_scaler.inverse_transform(predictions)
    labels = label_scaler.inverse_transform(labels)

    labels = (labels - predictions)
    return labels