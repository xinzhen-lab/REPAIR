# -*- coding: utf-8 -*-
#
# Copyright © dawnranger.
#
# 2018-05-08 10:15 <dawnranger123@gmail.com>
#
# Distributed under terms of the MIT license.
from __future__ import division, print_function
import numpy as np
import torch
from torch.utils.data import Dataset
from sklearn.metrics import normalized_mutual_info_score,f1_score, adjusted_rand_score, cluster,accuracy_score,precision_score,recall_score
from munkres import Munkres
pre = precision_score
rec = recall_score
Fscore = f1_score

def load_mnist(path='./data/mnist.npz'):
    f = np.load(path)

    x_train, y_train, x_test, y_test = f['x_train'], f['y_train'], f[
        'x_test'], f['y_test']
    f.close()
    x = np.concatenate((x_train, x_test))
    y = np.concatenate((y_train, y_test)).astype(np.int32)
    x = x.reshape((x.shape[0], -1)).astype(np.float32)
    x = np.divide(x, 255.)
    print('MNIST samples', x.shape)
    return x, y


class MnistDataset(Dataset):

    def __init__(self):
        self.x, self.y = load_mnist()

    def __len__(self):
        return self.x.shape[0]

    def __getitem__(self, idx):
        return torch.from_numpy(np.array(self.x[idx])), torch.from_numpy(
            np.array(self.y[idx])), torch.from_numpy(np.array(idx))


#######################################################
# Evaluate Critiron
#######################################################


def cluster_acc(y_true, y_pred):
    """
    Calculate clustering accuracy. Require scikit-learn installed

    # Arguments
        y: true labels, numpy.array with shape `(n_samples,)`
        y_pred: predicted labels, numpy.array with shape `(n_samples,)`

    # Return
        accuracy, in [0,1]
    """
    y_true = y_true.astype(np.int64)
    assert y_pred.size == y_true.size
    D = max(y_pred.max(), y_true.max()) + 1
    w = np.zeros((D, D), dtype=np.int64)
    for i in range(y_pred.size):
        w[y_pred[i], y_true[i]] += 1
    from scipy.optimize import linear_sum_assignment as linear_assignment
    # from sklearn.utils.linear_assignment_ import linear_assignment
    ind = linear_assignment(w.max() - w)
    ind1 = np.vstack(ind)
    ind2 = ind1.T
    return sum([w[i, j] for i, j in ind2]) * 1.0 / y_pred.size

def best_map(L1,L2):
    #L1 should be the groundtruth labels and L2 should be the clustering labels we got
    Label1 = np.unique(L1)
    nClass1 = len(Label1)
    Label2 = np.unique(L2)
    nClass2 = len(Label2)
    nClass = np.maximum(nClass1,nClass2)
    G = np.zeros((nClass,nClass))
    for i in range(nClass1):
        ind_cla1 = L1 == Label1[i]
        ind_cla1 = ind_cla1.astype(float)
        for j in range(nClass2):
            ind_cla2 = L2 == Label2[j]
            ind_cla2 = ind_cla2.astype(float)
            G[i,j] = np.sum(ind_cla2 * ind_cla1)
    m = Munkres()
    index = m.compute(-G.T)
    index = np.array(index)
    c = index[:,1]
    newL2 = np.zeros(L2.shape)
    for i in range(nClass2):
        newL2[L2 == Label2[i]] = Label1[c[i]]
    return newL2   

def acc_rate(gt_s, s):
    c_x = best_map(gt_s,s)
    err_x = np.sum(gt_s[:] == c_x[:])
    accrate = err_x.astype(float) / (gt_s.shape[0])

    return accrate 

def purity_score(y_true, y_pred):
    # compute contingency matrix (also called confusion matrix)
    contingency_matrix = cluster.contingency_matrix(y_true, y_pred)
    purity=np.sum(np.amax(contingency_matrix, axis=0)) / np.sum(contingency_matrix)
    # return purity
    return purity

def fscore_score(y_true, y_pred):
    # compute contingency matrix (also called confusion matrix)

    result_fscore = f1_score(y_true, y_pred, labels=None, pos_label=1, average='binary', sample_weight=None)
    # return purity
    return result_fscore

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

    return Sensitivity, Specificity