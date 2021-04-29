# Author Jean-Baptiste Guimbaud
# Meersens

# This is the main script. By default it computes mae score for predicting IQ with random forests.
# You can change the model to use and the target to predict at the begining of the sript (line 31).

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import shap
import graphviz
import os
import argparse
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2, f_regression, mutual_info_regression
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import export_graphviz

from src.analysis import cluster_corr, plot_corr
from src.utils import make_features, Phenotypes, Models, DataType, find_features_correlated_to
from src.utils import compute_residuals
from src.predictions import train_predict_and_test, cross_val, predict_naive
from src.insights import plot_features_importance, plot_predictions, shap_plots
from enum import Enum

# Here you can select the health outcome that you want to predict,
# the model, the data and wether you want to use feature selection
# technique or not.
TARGET = Phenotypes.IQ
MODEL = Models.RF
DATATYPE = DataType.BOTH
FEATURE_SELECTION = None # None, 'fwe', 'kbest', 'tree', 'corr

# Parse command line argument
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--cross",
                    action="store_true",
                    help="Compute score with cross validation and exit.")
parser.add_argument("-c", 
                    "--correlation",
                    action="store_true",
                    help="Print the correlation matrix")
parser.add_argument('-h',
                    '--help', 
                    action='help', 
                    default=argparse.SUPPRESS,
                    help='Show this help message and exit.')
args = parser.parse_args()
correlation = args.correlation
cross_validation = args.cross

# Get data
if DATATYPE == DataType.PREGNANCY:
    exposome = pd.read_csv("data/preprocessed/preg_exposome.csv")
elif DATATYPE == DataType.POSTNATAL:
    exposome = pd.read_csv("data/preprocessed/postnatal_exposome.csv")
else:
    exposome = pd.read_csv("data/preprocessed/exposome.csv")
covariates = pd.read_csv("data/preprocessed/covariates.csv")
phenotype = pd.read_csv("data/preprocessed/phenotype.csv")

# covariates.drop('e3_gac_None', axis=1, inplace=True)
if TARGET == Phenotypes.BODY_MASS_INDEX_CATEGORICAL or TARGET == Phenotypes.BODY_MASS_INDEX or TARGET == Phenotypes.BIRTH_WEIGHT:
    # if DATATYPE == DATATYPE.POSTNATAL or DATATYPE == DATATYPE.BOTH:
    covariates.drop('hs_c_weight_None', axis=1, inplace=True)
    covariates.drop('hs_c_height_None', axis=1, inplace=True)

# Drop IQ covariates
# Post natal
# if TARGET == Phenotypes.IQ:
#     # if DATATYPE == DataType.POSTNATAL or DATATYPE == DataType.BOTH:
#     covariates.drop('hs_child_age_None', axis=1, inplace=True)
#     covariates.drop('hs_c_height_None', axis=1, inplace=True)
#     covariates.drop('hs_c_weight_None', axis=1, inplace=True)

#     # pregnancy
#     # if DATATYPE == DataType.PREGNANCY or DATATYPE == DataType.BOTH:
#     covariates.drop('e3_yearbir_None', axis=1, inplace=True)

# Correlation matrix
if correlation:
    features = pd.merge(exposome, covariates, on="ID", how="inner")
    reindexed_features = cluster_corr(features)
    print(reindexed_features['ID'])
    print(features['ID'])
    reindexed_features = pd.merge(reindexed_features, phenotype, on="ID", how="inner")
    plot_corr(reindexed_features.corr())
    exit()

# Make features and target dataframes
features = make_features(exposome, covariates, phenotype)
target = features[TARGET.value]
features = features.drop(phenotype.columns, axis = 1)
features_columns = features.columns

# Removing specific features
# corr_features = find_features_correlated_to(feature_name="hs_child_age_None", threshold=0.9, dataframe=features)
# features = dataframe.drop(corr_features, axis=1)
# features = features.drop("hs_child_age_None", axis=1)
# features = features.drop("h_cohort", axis=1)

# predict residuals
# print(features['hs_child_age_None'].head())
# print(features["hs_child_age_None"].to_numpy().reshape(-1, 1))
# corr_features.add("hs_child_age_None")
# features_of_interest = features[corr_features]
# residuals = compute_residuals(features_of_interest, target.to_numpy().reshape(-1, 1))
# target = pd.Series(np.ravel(residuals))

# PCA reduction
# pca = PCA(n_components=40, svd_solver="full", random_state=42)
# pca.fit(features)
# print(features.shape)
# print(pca.n_components_)
# print(pca.explained_variance_ratio_)
# features = pca.transform(features)
# features = pd.DataFrame(features)


# Compute cross validation score
if cross_validation:
    scores = cross_val(MODEL, TARGET, features, target, features_selection=FEATURE_SELECTION)
    print(MODEL.name, "cross validation scores:")
    if TARGET == Phenotypes.DIAGNOSED_ASTHMA or TARGET == Phenotypes.BODY_MASS_INDEX_CATEGORICAL:
        print("  Mean Accuracy:", scores[0])
        print("  Mean Macro Roc Auc:", scores[1])
        print("  Mean Weighted (by prevalence) Roc Auc:", scores[2])
    else:
        print("  Mean MAE:", scores)
    exit()

# Compute prediction score on a test dataset
train_features, test_features, train_labels, test_labels = train_test_split(features,
                                                                            target,
                                                                            test_size = 0.25,
                                                                            random_state = 42)

print('Training Features Shape:', train_features.shape)
print('Training Labels Shape:', train_labels.shape)
print('Testing Features Shape:', test_features.shape)
print('Testing Labels Shape:', test_labels.shape)

# Train , predict and test
# naive model
predict_naive(train_features, train_labels, test_features, test_labels)

# model of interest
model, predictions = train_predict_and_test(MODEL,
                                            TARGET, 
                                            train_features, 
                                            train_labels, 
                                            test_features, 
                                            test_labels,
                                            FEATURE_SELECTION)

if MODEL == Models.CART:
    dot = export_graphviz(model, out_file=None, 
                         feature_names = features_columns,
                         class_names = True,
                         rounded = True, proportion = False, 
                         precision = 2, filled = True)
    graph = graphviz.Source(dot, format="png")
    graph.render('tree') 

# Features importance
if MODEL == Models.RF or MODEL == Models.XGB:
    plot_features_importance(model, train_features)


plot_predictions(TARGET.value, ground_truth=test_labels, predictions=predictions)
shap_plots(model, train_features, test_features, test_labels)