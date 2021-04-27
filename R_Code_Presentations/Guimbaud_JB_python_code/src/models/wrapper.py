import numpy as np
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import f1_score
from scipy import optimize
from sklearn.preprocessing import label_binarize

# Wrapper for classification threshold optimization
class WrappedClassifier():
    def __init__(self, origin_model):
        self.origin_model = origin_model

    def predict_proba(self, x, threshold_list=None):
        # get origin probability
        ori_proba = self.origin_model.predict_proba(x)

        # set default threshold
        if threshold_list is None:
            threshold_list = np.full(ori_proba[0].shape, 1)

        # get the output shape of threshold_list
        output_shape = np.array(threshold_list).shape

        print(threshold_list)
        # element-wise divide by the threshold of each classes
        new_proba = np.divide(ori_proba, threshold_list)

        # calculate the norm (sum of new probability of each classes)
        norm = np.linalg.norm(new_proba, ord=1, axis=1)

        # reshape the norm
        norm = np.broadcast_to(np.array([norm]).T, (norm.shape[0],output_shape[0]))

        # renormalize the new probability
        new_proba = np.divide(new_proba, norm)

        return new_proba

    def predict(self, x, threshold_list=None):
        return np.argmax(self.predict_proba(x, threshold_list), axis=1)
    
    def fit(self, X, y,              
            eval_metric=None, 
            eval_set=None,
            early_stopping_rounds=None,
            verbose=True):
        if eval_metric is None or eval_set is None or early_stopping_rounds is None:
            self.origin_model.fit(X, y)
        else:  # early stopping parameters for xgboost
            self.origin_model.fit(X, y,
                                eval_metric=eval_metric,
                                eval_set=eval_set,
                                early_stopping_rounds=early_stopping_rounds,
                                verbose=verbose)

def scoreFunc(model, X, y_true, threshold_list):
    y_pred = model.predict(X, threshold_list=threshold_list)
    y_pred_proba = model.predict_proba(X, threshold_list=threshold_list)

    # Binarize predictions and labels
    values, _ = np.unique(y_true, return_counts=True)
    num_classes = len(values)
    if num_classes > 2:
        y_true = label_binarize(y_true, classes=range(num_classes))
        y_pred = label_binarize(y_pred, classes=range(num_classes))

    # score = accuracy_score(y_true, y_pred)
    # score = roc_auc_score(y_true, y_pred_proba, average='macro', multi_class="ovr")
    # score = average_precision_score(y_true, y_pred_proba, average='macro')
    score = f1_score(y_true, y_pred, average='macro')
    return score


def weighted_score_with_threshold(threshold, model, X_test, Y_test, delta=5e-5):
    # if the sum of thresholds were not between 1+delta and 1-delta, 
    # return infinity (just for reduce the search space of the minimizaiton algorithm, 
    # because the sum of thresholds should be as close to 1 as possible).
    threshold_sum = np.sum(threshold)

    if threshold_sum > 1+delta:
        return np.inf

    if threshold_sum < 1-delta:
        return np.inf

    # to avoid objective function jump into nan solution
    if np.isnan(threshold_sum):
        print("threshold_sum is nan")
        return np.inf

    # renormalize: the sum of threshold should be 1
    normalized_threshold = threshold/threshold_sum

    # calculate scores based on thresholds
    # suppose it'll return 4 scores in a tuple: (accuracy, roc_auc, pr_auc, f1)
    score = scoreFunc(model, X_test, Y_test, threshold_list=normalized_threshold)    

    return -1 * score

def optimize_thresholds(proxy_model, output_class_num, X_test, Y_test):
    bounds = optimize.Bounds([1e-5]*output_class_num,[1]*output_class_num)

    result = optimize.differential_evolution(
        weighted_score_with_threshold,
        bounds,
        args=(proxy_model, X_test, Y_test),
        popsize=100,
        maxiter=5000,
        polish=True, 
        seed=42)

    # calculate threshold
    threshold = result.x/np.sum(result.x)

    # print the optimized score
    print(scoreFunc(proxy_model, X_test, Y_test, threshold_list=threshold))
    return threshold