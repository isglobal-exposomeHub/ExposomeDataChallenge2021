# Author Jean-Baptiste Guimbaud
# Meersens

import numpy as np
import matplotlib.pyplot as plt
import xgboost as xgb
import shap

from src.utils import Models
from src.models.mlp import MLP
import xgboost as xgb
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import label_binarize
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.neural_network import MLPRegressor, MLPClassifier
from sklearn.linear_model import ElasticNet, LogisticRegression

from sklearn.metrics import roc_curve, auc, roc_auc_score
from scipy import interp
from itertools import cycle


def plot_features_importance(model, train_features):
    if isinstance(model, RandomForestRegressor) or isinstance(model, RandomForestClassifier):
        print("Ploting features importance")
        return _plot_features_importance_rf(model, train_features)
    elif isinstance(model, xgb.sklearn.XGBRegressor) or isinstance(model, xgb.sklearn.XGBClassifier):
        print("Ploting features importance")
        return _plot_features_importance_xgb(model)


def _plot_features_importance_rf(rf, train_features):
    # Feature importance plot
    importances = rf.feature_importances_
    std = np.std([tree.feature_importances_ for tree in rf.estimators_],
                axis=0)
    indices = np.argsort(importances)[::-1]

    # Print the feature ranking
    print("Feature ranking:")

    for f in range(50):
        print("%d. feature %s (%f)" % (f + 1, train_features.columns[indices[f]], importances[indices[f]]))
    print("...")

    # Plot the impurity-based feature importances of the forest
    plt.figure()
    plt.title("Feature importances (50 most)")
    plt.bar(range(50), importances[indices[0:50]],
            color="r", yerr=std[indices[0:50]], align="center")
    plt.xticks(range(50), train_features.columns[indices[0:50]])
    plt.xlim([-1, 50])
    plt.show()


def _plot_features_importance_xgb(xg):
    xgb.plot_importance(xg)
    plt.show()

def plot_predictions(label, ground_truth, predictions):
    plt.figure(label)
    plt.subplot(311)
    plt.plot(ground_truth.to_numpy(), 'bx', predictions, 'gx')
    plt.legend(["ground truth", "predictions"], fancybox=True, shadow=True)
    plt.ylabel("Value")
    plt.xlabel("Data points")
    plt.xticks(np.arange(0, ground_truth.shape[0], 5.0))

    plt.subplot(312)
    plt.scatter(ground_truth.to_numpy(), predictions)
    plt.plot([ground_truth.min(), ground_truth.max()], [ground_truth.min(), ground_truth.max()], 'k--', lw=4)
    plt.xlabel('Measured')
    plt.ylabel('Predicted')

    plt.subplot(313)
    plt.hist(abs(ground_truth.to_numpy() - predictions), bins=10, label='Reversed emp.')
    plt.ylabel("Occurence")
    plt.xlabel("Absolute error")
    plt.show()

def print_regression_metrics(ground_truth, predictions, verbose=True):
    absolute_error = abs(ground_truth - predictions)
    mae = round(np.mean(absolute_error), 2)
    # Print out the mean absolute error (mae)
    if verbose:
        print('Mean Absolute Error:', mae)

    # Calculate mean absolute percentage error (MAPE)
    # mape = np.mean(100 * (error / test_labels.replace(0, 0.01)))
    # print("Mean Absolute Percent Error:", mape)

    # Calculate and display accuracy
    # accuracy = 100 - mape
    # print('Accuracy:', round(accuracy, 2), '%.')
    return mae

def print_classification_metrics(ground_truth, predictions, num_classes=2, verbose=True):
    correct_pred = (predictions == ground_truth)
    acc = correct_pred.sum() / len(correct_pred)
    acc = np.round(acc * 100)

    if verbose:
        print("Accuracy:", acc)
        print("Confusion matrix: ")
        print(confusion_matrix(ground_truth, predictions))
        print("Value counts in labels:")
        print(ground_truth.value_counts())

    # Binarize predictions and labels
    if num_classes > 2:
        ground_truth = label_binarize(ground_truth, classes=range(num_classes))
        predictions = label_binarize(predictions, classes=range(num_classes))

    # Compute OvR ROC AUC scores
    macro_roc_auc_ovr = roc_auc_score(ground_truth, predictions, multi_class="ovr",
                                  average="macro")
    weighted_roc_auc_ovr = roc_auc_score(ground_truth, predictions, multi_class="ovr",
                                     average="weighted")

    if verbose:
        print("One-vs-Rest ROC AUC scores:\n{:.6f} (macro),\n{:.6f} "
        "(weighted by prevalence)"
        .format(macro_roc_auc_ovr, weighted_roc_auc_ovr))

        # Compute ROC curve and ROC area for each class
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        for i in range(num_classes):
            fpr[i], tpr[i], _ = roc_curve(ground_truth[:, i], predictions[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])

        # Compute micro-average ROC curve and ROC area
        fpr["micro"], tpr["micro"], _ = roc_curve(ground_truth.ravel(), predictions.ravel())
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

        # First aggregate all false positive rates
        all_fpr = np.unique(np.concatenate([fpr[i] for i in range(num_classes)]))

        # Then interpolate all ROC curves at this points
        mean_tpr = np.zeros_like(all_fpr)
        for i in range(num_classes):
            mean_tpr += interp(all_fpr, fpr[i], tpr[i])

        # Finally average it and compute AUC
        mean_tpr /= num_classes

        fpr["macro"] = all_fpr
        tpr["macro"] = mean_tpr
        roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

        # Plot all ROC curves
        plt.figure()
        plt.plot(fpr["micro"], tpr["micro"],
                label='micro-average ROC curve (area = {0:0.2f})'
                    ''.format(roc_auc["micro"]),
                color='deeppink', linestyle=':', linewidth=4)

        plt.plot(fpr["macro"], tpr["macro"],
                label='macro-average ROC curve (area = {0:0.2f})'
                    ''.format(roc_auc["macro"]),
                color='navy', linestyle=':', linewidth=4)

        lw = 2
        colors = cycle(['aqua', 'darkorange', 'cornflowerblue', 'green'])
        for i, color in zip(range(num_classes), colors):
            plt.plot(fpr[i], tpr[i], color=color, lw=lw,
                    label='ROC curve of class {0} (area = {1:0.2f})'
                    ''.format(i, roc_auc[i]))

        plt.plot([0, 1], [0, 1], 'k--', lw=lw)
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC')
        plt.legend(loc="lower right")
        plt.show()

    return acc, macro_roc_auc_ovr, weighted_roc_auc_ovr

def shap_plots(model, train_features, test_features, test_labels):
    print("Computing shapley values..")
    # compute SHAP values
    if isinstance(model, (MLP, MLPRegressor, MLPClassifier, ElasticNet, LogisticRegression)):
        train_sample = shap.sample(train_features, 10)
        explainer = shap.Explainer(model.predict, train_sample)
    elif isinstance(model, (RandomForestRegressor, RandomForestClassifier)):
        explainer = shap.TreeExplainer(model, train_features)
    else:
        explainer = shap.Explainer(model, train_features)
    
    shap_values = explainer(test_features)
    shap.plots.bar(shap_values, max_display=10)
    # shap.plots.bar(shap_values[0]) # Local 


    # beeswarm plot
    shap.plots.beeswarm(shap_values)

    # Decision plot
    expected_value = explainer.expected_value
    select = range(20)
    features_sample = test_features.iloc[select]
    shap.decision_plot(expected_value, explainer.shap_values(features_sample), features_sample)

    # Heatmap 
    shap.plots.heatmap(shap_values, max_display=10)

    # Scatter
    shap.plots.scatter(shap_values[:, "hs_child_age_None"],
                       color=shap_values,
                       alpha=0.8)

    # Feature clustering (redondant feature detection)
    clustering = shap.utils.hclust(test_features, test_labels) # by default this trains (X.shape[1] choose 2) 2-feature XGBoost models
    shap.plots.bar(shap_values, clustering=clustering, clustering_cutoff=0.5)