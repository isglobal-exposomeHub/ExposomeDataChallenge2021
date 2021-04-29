# Author Jean-Baptiste Guimbaud
# Meersens

import os
import numpy as np
import matplotlib.pyplot as plt
import configparser
import pprint
import pandas as pd
from src.utils import Models

from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectFwe, SelectFromModel, SelectKBest, VarianceThreshold
from sklearn.feature_selection import f_regression, f_classif
from sklearn.linear_model import ElasticNet, LogisticRegression
from sklearn.neural_network import MLPRegressor, MLPClassifier
from sklearn.svm import SVC, SVR
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor, plot_tree
import xgboost as xgb
from xgboost.sklearn import XGBClassifier, XGBRegressor, XGBRFRegressor
import imblearn
import torch
import tensorflow as tf
from src.utils import Phenotypes
from src.insights import print_classification_metrics, print_regression_metrics
from src.models.mlp import MLP
from src.models.wrapper import WrappedClassifier, optimize_thresholds



def cross_val(model_name, target_name, features, labels, features_selection=None):
# def score_model(model, params, sampling_strategy, X_train, y_train):
    classification = (target_name == Phenotypes.DIAGNOSED_ASTHMA or
                      target_name == Phenotypes.BODY_MASS_INDEX_CATEGORICAL)
    cv = KFold(n_splits=5, shuffle=True, random_state=42)

    scores = []
    for train_fold_index, val_fold_index in cv.split(features, labels):

        # Get the training data
        X_train_fold, y_train_fold = features.iloc[train_fold_index], labels[train_fold_index]
        # Get the validation data
        X_val_fold, y_val_fold = features.iloc[val_fold_index], labels[val_fold_index]

        # Fit the model on the upsampled training data
        model, predictions = train_predict_and_test(model_name,
                                                    target_name,
                                                    X_train_fold, 
                                                    y_train_fold, 
                                                    X_val_fold, 
                                                    y_val_fold,
                                                    feature_selection=features_selection)
        # model_obj = model(**params).fit(X_train_fold_upsample, y_train_fold_upsample)
        
        if classification:
            acc, macro_roc_auc_ovr, weighted_roc_auc_ovr = print_classification_metrics(
                ground_truth=y_val_fold,
                predictions=predictions,
                num_classes=y_val_fold.nunique(),
                verbose=False)
            scores.append([acc, macro_roc_auc_ovr, weighted_roc_auc_ovr])
        else:
            mae = print_regression_metrics(ground_truth=y_val_fold, predictions=predictions, verbose=False)
            scores.append(mae)

    print(scores)
    if classification:
        return np.mean(scores, axis=0)        
    return np.mean(scores)


def train_predict_and_test(model,
                           target_name, 
                           train_features, 
                           train_labels, 
                           test_features, 
                           test_labels,
                           feature_selection=None):
    classification = (target_name == Phenotypes.DIAGNOSED_ASTHMA or
                      target_name == Phenotypes.BODY_MASS_INDEX_CATEGORICAL)
    
    # Standardize data
    standardized = False
    if model == Models.MLP or model == Models.SVM:
        print("Standardizing data..")
        standardized = True
        features_mean = train_features.mean()
        features_std = train_features.std()
        train_features = (train_features - features_mean) / features_std
        test_features = (test_features - features_mean) / features_std

        if not classification:
            labels_mean = train_labels.mean()
            labels_std = train_labels.std()
            train_labels = (train_labels - labels_mean) / labels_std
            test_labels = (test_labels - labels_mean) / labels_std
    
    # Load optimized params
    params = load_optimized_params(model, target_name)

    # Features selection
    feature_selector = VarianceThreshold(threshold=0).fit(train_features) # Removing features with 0 variance
    train_col, test_col = train_features.columns, test_features.columns
    train_features = pd.DataFrame(feature_selector.transform(train_features), columns=train_col)
    test_features = pd.DataFrame(feature_selector.transform(test_features), columns=test_col)
    if feature_selection == "fwe":
        print("Selecting features according to Familly Wise Error")
        # alpha = 5e-2
        alpha = 0.3
        if params is not None:
            try:
                alpha = params['transformer_alpha']
            except KeyError:
                print("Cannot find parameter alpha for FWE feature selector. Using default value")

        features_selector = SelectFwe(f_regression, alpha=alpha).fit(train_features, train_labels)
        train_features = features_selector.transform(train_features)
        test_features = features_selector.transform(test_features)
    elif feature_selection == "kbest":
        k = 150
        if params is not None:
            try:
                k = params['k']
            except KeyError:
                print("Cannot find parameter k for k-best feature selector. Using default value: k=", k)
        print("Selecting k-best features:", k)
        score_func = f_regression
        if classification:
            score_func = f_classif
        features_selector = SelectKBest(score_func=score_func, k=k)
        features_selector = features_selector.fit(train_features,train_labels)
        train_features = features_selector.transform(train_features)
        test_features = features_selector.transform(test_features)
    elif feature_selection == "tree":
        print("Selecting features from RF feature importance")
        clf = RandomForestRegressor(n_estimators=100).fit(train_features, train_labels)
        if classification:
            clf = RandomForestClassifier(n_estimators=100).fit(train_features, train_labels)
        features_selector = SelectFromModel(clf, prefit=True)
        train_features = features_selector.transform(train_features)
        test_features = features_selector.transform(test_features)
    elif feature_selection == "corr":
        threshold = 0.9 # Recommended default value
        col_corr = set()
        corr_matrix = train_features.corr()
        for i in range(len(corr_matrix.columns)):
            for j in range(i):
                if abs(corr_matrix.iloc[i, j]) > threshold:
                    colname = corr_matrix.columns[i]
                    col_corr.add(colname)
        train_features = train_features.drop(col_corr, axis = 1)
        test_features = test_features.drop(col_corr, axis = 1)

    # Oversampling
    if classification and model != Models.SVM and model != Models.CART and model != Models.ELASTIC:
        print("Oversampling features..")
        if target_name == Phenotypes.DIAGNOSED_ASTHMA:
            sampling_strat = 0.5
        else:
            sampling_strat = {0: np.max(np.bincount(train_labels))//4,
                              1: np.max(np.bincount(train_labels)),
                              2: np.max(np.bincount(train_labels)),
                              3: np.max(np.bincount(train_labels))//2}
        oversampler = imblearn.over_sampling.RandomOverSampler(sampling_strategy=sampling_strat,
                                                                random_state=42)
        # oversampler = imblearn.over_sampling.SMOTE(sampling_strategy=1.0,
        #                                          k_neighbors=5,
        #                                          random_state=42)
        train_features, train_labels = oversampler.fit_resample(train_features, train_labels)

    if model == Models.RF:
        if target_name == Phenotypes.BODY_MASS_INDEX_CATEGORICAL:
            # Create validation set for threshold optimization
            val_features, test_features, val_labels, test_labels = train_test_split(test_features, test_labels,
                                                        test_size=0.5,
                                                        random_state=42)
            model, predictions = _predict_rf(target_name, train_features, train_labels, val_features, val_labels)
        else:
            model, predictions = _predict_rf(target_name,
                                             train_features, 
                                             train_labels, 
                                             test_features, 
                                             test_labels, 
                                             params=params)
    elif model == Models.ELASTIC:
        model, predictions = predict_elastic_net(target_name,
                                                 train_features,
                                                 train_labels,
                                                 test_features,
                                                 test_labels)
    elif model == Models.XGB:
        model, predictions = _predict_xgb(target_name, 
                                          train_features, 
                                          train_labels, 
                                          test_features, 
                                          test_labels,
                                          params=params)
    elif model == Models.MLP:
        model, predictions = _predict_mlp(target_name, 
                                          train_features, 
                                          train_labels, 
                                          test_features, 
                                          test_labels,
                                          params=params)
    elif model == Models.SVM:
        model, predictions = _predict_svm(target_name, train_features, train_labels, test_features, test_labels)
    elif model == Models.CART:
        model, predictions = _predict_cart(target_name, train_features, train_labels, test_features, test_labels)
    elif model == Models.NAIVE:
        if not(classification):
            predictions = predict_naive(train_features, train_labels, test_features, test_labels)
        else:
            raise SystemExit("Cannot use naive model on classification task")
    else:
        raise SystemExit("Unkwown model:", model)

    # Destandardize results
    if standardized and not(classification):
        print("destandardize data..")
        predictions = (predictions * labels_std) + labels_mean
        test_labels = (test_labels * labels_std) + labels_mean

    # Print results
    if classification:
        print_classification_metrics(ground_truth=test_labels,
                                    predictions=predictions,
                                    num_classes=test_labels.nunique())
    else:
        print_regression_metrics(ground_truth=test_labels, predictions=predictions)

    return model, predictions


def load_optimized_params(model_name, target_name):
    if model_name == Models.MLP:
        params = load_params_mlp(section_name=target_name.name)
    elif model_name == Models.XGB:
        params = load_params_xgb(section_name=target_name.name)
    elif model_name == Models.RF:
        params = load_params_rf(section_name=target_name.name)
    else:
        print("No config file for model:", model_name.name)
        params = None
    return params


def load_params_xgb(section_name: str):
    config = configparser.ConfigParser()
    config.read('config/config_xgb.ini')
    section = config[section_name]
    params = {
        'learning_rate' : section.getfloat('learning_rate'),
        'n_estimators' : section.getint('n_estimators'),
        'max_depth': section.getint('max_depth'),
        'min_child_weight': section.getint('min_child_weight'),
        'gamma': section.getfloat('gamma'),
        'subsample': section.getfloat('subsample'),
        'colsample_bytree': section.getfloat('colsample_bytree'),
        'objective': section.get('objective'),
        'eval_metric': section.get('eval_metric'),
        'booster': section.get('booster'),
        'use_label_encoder': section.getboolean('use_label_encoder'),
        'seed': section.getint('seed')
    }
    if section_name == "BODY_MASS_INDEX_CATEGORICAL":
        params["num_class"] = section.getint("num_class")

    print("loaded hyperparameters: ")
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(params)
    return params

def load_params_mlp(section_name: str):
    config = configparser.ConfigParser()
    config.read('config/config_mlp.ini')
    section = config[section_name]
    params = {
        'n_layers' : section.getint('n_layers'),
        'n_units' : eval(section.get('n_units')),
        'learning_rate': section.getfloat('learning_rate'),
        'momentum': section.getfloat('momentum'),
        'dropout': section.getfloat('dropout'),
        'regularizer': section.get('regularizer'),
        'feature_selector': section.get('feature_selector'),
        'k': section.getint('k'),
        'transformer_alpha' : section.getfloat('transformer_alpha'),
        'optimizer' : section.get('optimizer'),
        'threshold' : section.getfloat('threshold')
    }

    print("loaded hyperparameters: ")
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(params)
    return params


def load_params_rf(section_name: str):
    config = configparser.ConfigParser()
    config.read('config/config_rf.ini')
    section = config[section_name]
    params = {
        'n_estimators' : section.getint('n_estimators'),
        'max_depth' : section.getint('max_depth'),
        'min_samples_split': section.getint('min_samples_split'),
        'min_samples_leaf': section.getint('min_samples_leaf')
    }
    if params['max_depth'] == -1:
        params['max_depth'] = None

    print("loaded hyperparameters: ")
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(params)
    return params
        

def predict_naive(train_features, train_labels, test_features, test_labels):
    naive_predictions = np.repeat(np.mean(test_labels), test_labels.shape[0])
    error = abs(naive_predictions - test_labels)
    print('naive mae :', round(np.mean(error), 2))
    return naive_predictions

def predict_elastic_net(target_name, train_features, train_labels, test_features, test_labels):
    if target_name == Phenotypes.DIAGNOSED_ASTHMA:
        regr = LogisticRegression(penalty='elasticnet',
                                 multi_class='ovr',
                                 class_weight='balanced',
                                 solver='saga',
                                 l1_ratio=0.5,
                                 random_state=42)
    elif target_name == Phenotypes.BODY_MASS_INDEX_CATEGORICAL:
        regr = LogisticRegression(penalty='elasticnet',
                                  multi_class='multinomial',
                                  class_weight='balanced',
                                  solver='saga',
                                  l1_ratio=0.5,
                                  random_state=42)
    else: # Classical regression task
        regr = ElasticNet(random_state=42)

    regr = regr.fit(train_features, train_labels)
    predictions = regr.predict(test_features)
    return regr, predictions


def _predict_rf(target_name, train_features, train_labels, test_features, test_labels, params=None):
    classification = (target_name == Phenotypes.DIAGNOSED_ASTHMA or
                    target_name == Phenotypes.BODY_MASS_INDEX_CATEGORICAL)
    if params is None:
        params = load_params_rf(target_name.name)
    if classification:
        rf = RandomForestClassifier(n_estimators=params['n_estimators'],
                                    max_depth=params['max_depth'],
                                    min_samples_split=params['min_samples_split'],
                                    min_samples_leaf=params['min_samples_leaf'],
                                    random_state = 42)
        rf = WrappedClassifier(rf)
    else:
        rf = RandomForestRegressor(n_estimators = params['n_estimators'],
                                   max_depth=params['max_depth'],
                                   min_samples_split=params['min_samples_split'],
                                   min_samples_leaf=params['min_samples_leaf'],
                                   random_state = 42)

    rf.fit(train_features, train_labels)
    if classification:
        thresholds = optimize_thresholds(proxy_model=rf,
                                         output_class_num=test_labels.nunique(),
                                         X_test=test_features,
                                         Y_test=test_labels)
        predictions = rf.predict(test_features, threshold_list=thresholds)
    else:
        predictions = rf.predict(test_features)
    # predictions = rf.predict(test_features)
    return rf, predictions


def _predict_xgb(target_name,
                 train_features, 
                 train_labels, 
                 test_features, 
                 test_labels,
                 params = None):
    if params is None:
        params = load_params_xgb(section_name=target_name.name)
    classification = (target_name == Phenotypes.DIAGNOSED_ASTHMA or
                      target_name == Phenotypes.BODY_MASS_INDEX_CATEGORICAL)
    if classification:
        model = XGBClassifier(**params)
        model = WrappedClassifier(model)
        eval_metric = 'mlogloss'
    else:
        model = XGBRegressor(**params)
        eval_metric = 'mae'

    eval_set = [(test_features, test_labels)]
    model.fit(train_features, train_labels, 
              eval_metric=eval_metric, 
              eval_set=eval_set,
              early_stopping_rounds=20, 
              verbose=True
              )

    if classification:
        thresholds = optimize_thresholds(proxy_model=model,
                                         output_class_num=test_labels.nunique(),
                                         X_test=test_features,
                                         Y_test=test_labels)
        predictions = model.predict(test_features, threshold_list=thresholds)
        # predictions = model.predict(test_features)
    else:
        predictions = model.predict(test_features)
    return model, predictions


# def _predict_mlp(target_name,
#                  train_features,
#                  train_labels,
#                  test_features,
#                  test_labels,
#                  params = None):
#     if params is None:
#         params = load_params_mlp(section_name=target_name.name)
#     classification = (target_name == Phenotypes.DIAGNOSED_ASTHMA or
#                     target_name == Phenotypes.BODY_MASS_INDEX_CATEGORICAL)

#     # Instanciate model
#     if classification:
#         model = MLPClassifier(solver='adam',
#                             alpha=params['alpha'],
#                             hidden_layer_sizes=params['n_units'],
#                             learning_rate_init=params['learning_rate'],
#                             activation='relu',
#                             batch_size=64,
#                             max_iter=1000,
#                             early_stopping=True,
#                             validation_fraction=0.2,
#                             n_iter_no_change=20,
#                             random_state=42)
#     else: # Classical regression task
#         model = MLPRegressor(solver='adam',
#                              alpha=params['alpha'],
#                              hidden_layer_sizes=params['n_units'],
#                              learning_rate_init=params['learning_rate'],
#                              activation='relu',
#                              batch_size=64,
#                              max_iter=1000,
#                              early_stopping=True,
#                              validation_fraction=0.2,
#                              n_iter_no_change=20,
#                              random_state=42)

#     model = model.fit(train_features, train_labels)
#     predictions = model.predict(test_features)

#     return model, predictions

def _predict_mlp(target_name,
                 train_features,
                 train_labels,
                 test_features,
                 test_labels,
                 params = None):
    if params is None:
        params = load_params_mlp(section_name=target_name.name)
    classification = (target_name == Phenotypes.DIAGNOSED_ASTHMA or
                    target_name == Phenotypes.BODY_MASS_INDEX_CATEGORICAL)

    regularizer = params['regularizer']
    if regularizer == 'None':
        regularizer = None
    model = MLP(layers_size=params['n_units'], 
                num_features=train_features.shape[1], 
                dropout=params['dropout'],
                regularizer=regularizer,
                activation='elu',
                classifier=classification)

    optimizer_type = params['optimizer']
    if optimizer_type == "Adam":
        optimizer = tf.keras.optimizers.Adam(learning_rate=params['learning_rate'], name='Adam')
    else:
        optimizer = tf.keras.optimizers.RMSprop(learning_rate=params['learning_rate'],
                                                momentum=params['momentum'],
                                                name='RMSprop')
    # Instanciate model
    if classification:
        loss_object = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=False)
        metric_object, metric_name, mode = "accuracy", 'accuracy', 'max'
    else: # Classical regression task
        loss_object = tf.losses.MeanSquaredError()
        metric_object, metric_name, mode = tf.metrics.MeanAbsoluteError(), 'mean_absolute_error', 'min'
    model.compile(loss=loss_object, optimizer=optimizer, metrics=[metric_object])
    callbacks = [tf.keras.callbacks.EarlyStopping(monitor=metric_name,
                                                  patience=10,
                                                  min_delta=0.0002,
                                                  mode=mode,
                                                  verbose=1, 
                                                  restore_best_weights=True),
                #  tf.keras.callbacks.ReduceLROnPlateau(monitor=metric_name,
                #                                       factor=0.5,
                #                                       patience=2,
                #                                       verbose=1, 
                #                                       mode='auto', 
                #                                       min_delta=0.0001, 
                #                                       cooldown=0,
                #                                       min_lr=1e-6)
                ]

    model.fit(train_features, train_labels,
              epochs=1000, 
              batch_size=64, 
              verbose=1, 
              validation_split=0.2,
              callbacks=callbacks
              )
    print(model.summary())
    predictions = model.predict(test_features)
    if classification:
        predictions = np.argmax(predictions, axis=1)

    return model, predictions


def _predict_svm(target_name,
                 train_features,
                 train_labels,
                 test_features,
                 test_labels):
    classification = (target_name == Phenotypes.DIAGNOSED_ASTHMA or
                      target_name == Phenotypes.BODY_MASS_INDEX_CATEGORICAL)

    # Instanciate model
    if classification:
        model = SVC(C=1,
                    gamma='scale',
                    class_weight="balanced",
                    random_state=42)
    else: # Classical regression task
        model = SVR(C=0.4,
                    gamma='scale')

    model.fit(train_features, train_labels)
    predictions = model.predict(test_features)

    return model, predictions


def _predict_cart(target_name,
                  train_features,
                  train_labels,
                  test_features,
                  test_labels):

    classification = (target_name == Phenotypes.DIAGNOSED_ASTHMA or
                      target_name == Phenotypes.BODY_MASS_INDEX_CATEGORICAL)

    # Instanciate model
    if classification:
        model = DecisionTreeClassifier(max_depth=10,
                                       max_features=0.9,
                                       class_weight="balanced",
                                       random_state=42)
    else: # Classical regression task
        model = DecisionTreeRegressor(max_depth=4,
                                      max_features=0.7,
                                      random_state=42)

    model.fit(train_features, train_labels)
    predictions = model.predict(test_features)

    # plot_tree(model)
    # plt.plot()

    return model, predictions


def _compute_accuracy_from_class_probability(y_true, y_pred):
    y_pred_softmax = torch.log_softmax(y_pred, dim = 1)
    _, y_pred_tags = torch.max(y_pred_softmax, dim = 1)    
    
    correct_pred = (y_pred_tags == y_true).float()
    acc = correct_pred.sum() / len(correct_pred)
    
    acc = torch.round(acc * 100)
    if acc > 100:
        print(acc)
    return acc