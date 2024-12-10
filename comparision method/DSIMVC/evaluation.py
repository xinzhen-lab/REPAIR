from sklearn.metrics import v_measure_score, accuracy_score
from scipy.optimize import linear_sum_assignment
import numpy as np


def cluster_acc(y_true, y_pred):
    y_true = y_true.astype(np.int64)
    assert y_pred.size == y_true.size
    D = max(y_pred.max(), y_true.max()) + 1
    w = np.zeros((D, D), dtype=np.int64)
    for i in range(y_pred.size):
        w[y_pred[i], y_true[i]] += 1
    u = linear_sum_assignment(w.max() - w)
    ind = np.concatenate([u[0].reshape(u[0].shape[0], 1), u[1].reshape([u[0].shape[0], 1])], axis=1)
    acc=sum([w[i, j] for i, j in ind]) * 1.0 / y_pred.size
    TN=w[0,0]
    TP=w[1,1]
    FN=w[0,1]
    FP=w[1,0]
    sen = float(TP) / (float(TP + FN) + 0.000000001)
    spe = float(TN) / (float(TN + FP) + 0.000000001)
    bacc=(spe+sen)/2

    return acc,bacc,sen,spe


def purity(y_true, y_pred):
    y_voted_labels = np.zeros(y_true.shape)
    labels = np.unique(y_true)
    ordered_labels = np.arange(labels.shape[0])
    for k in range(labels.shape[0]):
        y_true[y_true == labels[k]] = ordered_labels[k]
    labels = np.unique(y_true)
    bins = np.concatenate((labels, [np.max(labels)+1]), axis=0)

    for cluster in np.unique(y_pred):
        hist, _ = np.histogram(y_true[y_pred == cluster], bins=bins)
        winner = np.argmax(hist)
        y_voted_labels[y_pred == cluster] = winner

    return accuracy_score(y_true, y_voted_labels)

def cal_sen_spe(y_actual, y_pred):
    TP = 0.0
    FP = 0.0
    TN = 0.0
    FN = 0.0

    for i in range(len(y_pred)):
        if y_actual[i] == y_pred[i] == 1:
            TP += 1
    for i in range(len(y_pred)):
        if y_actual[i] == 0 and y_actual[i] != y_pred[i]:
            FP += 1
    for i in range(len(y_pred)):
        if y_actual[i] == y_pred[i] == 0:
            TN += 1
    for i in range(len(y_pred)):
        if y_actual[i] == 1 and y_actual[i] != y_pred[i]:
            FN += 1

    if ((y_actual == y_pred).all() and np.sum(y_actual) == len(y_actual)) or ((y_actual == y_pred).all() and np.sum(y_actual) == 0):  # all 1s or 0s
        Sensitivity = 1.0
        Specificity = 1.0
    elif (y_actual == y_pred).any() == False:
        Sensitivity = 0.0
        Specificity = 0.0
    else:
        Sensitivity = float(TP) / (float(TP + FN) + 0.000000001)
        Specificity = float(TN) / (float(TN + FP) + 0.000000001)

    # Overall accuracy
    ACC = (TP + TN) / (TP + FP + FN + TN)
    bacc=(Sensitivity+Specificity)/2

    return Sensitivity, Specificity, bacc
def evaluate(label, pred):
    nmi = v_measure_score(label, pred)
    acc,bacc,Sensitivity, Specificity = cluster_acc(label, pred)
    pur = purity(label, pred)
    # Sensitivity, Specificity, bacc=cal_sen_spe(label, pred)

    return acc, nmi, pur,Sensitivity, Specificity,bacc

