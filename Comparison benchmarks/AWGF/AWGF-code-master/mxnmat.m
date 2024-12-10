clear;
clc
load('F:\incomplete modle\1_Adaptive Weighted Graph Fusion Incomplete\multi-view-datasets-master\random remove\GBM6Flair_mxn1.mat')
load('F:\incomplete modle\1_Adaptive Weighted Graph Fusion Incomplete\multi-view-datasets-master\random remove\GBM6Flair_mxn_percent0.4.mat')
data=X;
miss40=folds;
save('GBM6Flair_mxn.mat','data','miss40','truth');