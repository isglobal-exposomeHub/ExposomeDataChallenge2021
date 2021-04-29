# Author Jean-Baptiste Guimbaud
# Meersens

# This script preprocess data (categorical data are one hot encoded, tertiles are replaced with the lower bound)

import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import xgboost as xgb
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor

LABELISE = False # Wether to labelise or one-hot encode categorical data

exposome = pd.read_csv("data/exposome/exposome.csv")
print("Loading exposome table:")
print(exposome.head())
print("dataframe of", exposome.shape[0], "rows total")


# Separate pregnancy and postnatal features
known_pregnancy_col = ['e3_alcpreg_yn_None', 'h_pamod_t3_None', 'h_pavig_t3_None', 'hs_tl_mdich_None',
                       'e3_asmokcigd_p_None', 'hs_cotinine_mcat_None', 'h_mbmi_None', 'hs_wgtgain_None',
                       'e3_gac_None', 'e3_sex_None', 'e3_yearbir_None', 'h_age_None', 'h_cohort',
                       'h_edumc_None', 'h_native_None', 'h_parity_None']
pregnancy_exposome, postnatal_exposome = pd.DataFrame(), pd.DataFrame()
pregnancy_exposome['ID'], postnatal_exposome['ID'] = exposome['ID'], exposome['ID']

# Separate prenatal and post natal data
for col in exposome.columns:
    if col.find("_preg_") != -1 or \
       col.find("_madj_") != -1 or \
       col.find("_m_") != -1 or \
       col in known_pregnancy_col:
        pregnancy_exposome[col] = exposome[col]
    else:
        postnatal_exposome[col] = exposome[col]

one_hot_columns = []
for col in exposome.columns:
    if exposome[col].dtypes == "object":
        if col.endswith("Ter"): # tertile data
            for index, value in exposome[col].items():
                lb = re.search('\((.*),', value)
                exposome.loc[index, col] = lb.group(1)
        else: # categorical data
            if LABELISE:
                label_encoder = preprocessing.LabelEncoder()
                exposome[col] = label_encoder.fit_transform(exposome[col])
            else: # one hot encoding
                one_hot_columns.append(col)
                dummies = pd.get_dummies(exposome[col], prefix=col)
                exposome = pd.concat((exposome, dummies), axis=1)
                exposome.drop([col], axis=1, inplace=True)
exposome.to_csv("data/preprocessed/exposome.csv", index=False)

for col in pregnancy_exposome.columns:
    if pregnancy_exposome[col].dtypes == "object":
        if col.endswith("Ter"): # tertile data
            for index, value in pregnancy_exposome[col].items():
                lb = re.search('\((.*),', value)
                pregnancy_exposome.loc[index, col] = lb.group(1)
        else: # categorical data
            if LABELISE:
                label_encoder = preprocessing.LabelEncoder()
                pregnancy_exposome[col] = label_encoder.fit_transform(pregnancy_exposome[col])
            else:
                dummies = pd.get_dummies(pregnancy_exposome[col], prefix=col)
                pregnancy_exposome = pd.concat((pregnancy_exposome, dummies), axis=1)
                pregnancy_exposome.drop([col], axis=1, inplace=True)
pregnancy_exposome.to_csv("data/preprocessed/preg_exposome.csv", index=False)

for col in postnatal_exposome.columns:
    if postnatal_exposome[col].dtypes == "object":
        if col.endswith("Ter"): # tertile data
            for index, value in postnatal_exposome[col].items():
                lb = re.search('\((.*),', value)
                postnatal_exposome.loc[index, col] = lb.group(1)
        else: # categorical data
            if LABELISE:
                label_encoder = preprocessing.LabelEncoder()
                postnatal_exposome[col] = label_encoder.fit_transform(postnatal_exposome[col])
            else:
                dummies = pd.get_dummies(postnatal_exposome[col], prefix=col)
                postnatal_exposome = pd.concat((postnatal_exposome, dummies), axis=1)
                postnatal_exposome.drop([col], axis=1, inplace=True)
postnatal_exposome.to_csv("data/preprocessed/postnatal_exposome.csv", index=False)


print("Loading covariates:")
covariates = pd.read_csv("data/exposome/covariates.csv")
print(covariates.head())
print("dataframe of", covariates.shape[0], "rows total")

pregnancy_covariates, postnatal_covariates = pd.DataFrame(), pd.DataFrame()
pregnancy_covariates['ID'], postnatal_covariates['ID'] = covariates['ID'], covariates['ID']

for col in covariates.columns:
    if col.find("_preg_") != -1 or \
       col.find("_madj_") != -1 or \
       col.find("_m_") != -1 or \
       col in known_pregnancy_col:
        pregnancy_covariates[col] = covariates[col]
    else:
        postnatal_covariates[col] = covariates[col]

for col in covariates.columns:
    if covariates[col].dtypes == "object": # Categorical variables
        if LABELISE:
            label_encoder = preprocessing.LabelEncoder()
            covariates[col] = label_encoder.fit_transform(covariates[col])
        else:
            dummies = pd.get_dummies(covariates[col], prefix=col)
            covariates = pd.concat((covariates, dummies), axis=1)
            covariates.drop([col], axis=1, inplace=True)
covariates.to_csv("data/preprocessed/covariates.csv", index=False)

for col in pregnancy_covariates.columns:
    if pregnancy_covariates[col].dtypes == "object": # Categorical variables
        if LABELISE:
            label_encoder = preprocessing.LabelEncoder()
            pregnancy_covariates[col] = label_encoder.fit_transform(pregnancy_covariates[col])
        else:
            dummies = pd.get_dummies(pregnancy_covariates[col], prefix=col)
            pregnancy_covariates = pd.concat((pregnancy_covariates, dummies), axis=1)
            pregnancy_covariates.drop([col], axis=1, inplace=True)
pregnancy_covariates.to_csv("data/preprocessed/pregnancy_covariates.csv", index=False)

for col in postnatal_covariates.columns:
    if postnatal_covariates[col].dtypes == "object": # Categorical variables
        if LABELISE:
            label_encoder = preprocessing.LabelEncoder()
            postnatal_covariates[col] = label_encoder.fit_transform(postnatal_covariates[col])
        else:
            dummies = pd.get_dummies(postnatal_covariates[col], prefix=col)
            postnatal_covariates = pd.concat((postnatal_covariates, dummies), axis=1)
            postnatal_covariates.drop([col], axis=1, inplace=True)
postnatal_covariates.to_csv("data/preprocessed/postnatal_covariates.csv", index=False)

print("Loading phenotypes:")
phenotype = pd.read_csv("data/exposome/phenotype.csv")
label_encoder = preprocessing.LabelEncoder()
label_encoder.fit(phenotype['hs_bmi_c_cat'])
phenotype['hs_bmi_c_cat'] = label_encoder.transform(phenotype['hs_bmi_c_cat'])
pheno_columns = phenotype.columns
print(phenotype.head())
print("dataframe of", phenotype.shape[0], "rows total")
phenotype.to_csv("data/preprocessed/phenotype.csv", index=False)
