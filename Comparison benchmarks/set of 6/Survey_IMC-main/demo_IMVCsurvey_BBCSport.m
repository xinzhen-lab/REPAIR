% If you use the code, please cite the following references:
% [1] Jie Wen, Zheng Zhang, Lunke Fei, Bob Zhang, Yong Xu, Zhao Zhang, Jinxing Li, A Survey on Incomplete Multi-view Clustering, IEEE TRANSACTIONS ON SYSTEzhong, MAN, AND CYBERNETICS: SYSTEzhong, 2022.
% [2] Wen J, Xu Y, Liu H. Incomplete multiview spectral clustering with adaptive graph learning[J]. IEEE Transactions on Cybernetics, 2020, 50(4): 1418-1429.
% The following released codes are written  by Jie Wen. For any problezhong, contact jiewen_pr@126.com
% Sugguested Matlab version 2015a. For different versions, the results will be different. 
% 
clear all
clc
Dataname = 'LGG 1p19g1_mxn';
addpath('dataset');
addpath('clustering_metrics');
addpath('utils');
for percentDel = 0
    fprintf(['nfold=',num2str(percentDel),'\n'])
f = 1;
% -------- BSV ------ %
addpath('BSV_Concat')
CLU_BSVresult = BSV_clustering(Dataname,percentDel,f)

% -------- Concat ------ %
addpath('BSV_Concat')
CLU_Concatresult = Concat_clustering(Dataname,percentDel,f)

% -------- DAIMC -------- %
addpath('DAIMC')
CLU_DAIMCresult = DAIMC_clustering(Dataname,percentDel,f)
% 
% % -------- OPIMC -------- %
addpath('OPIMC')
CLU_OPIMCresult = OPIMC_clustering(Dataname,percentDel,f)
% 
% 
% % -------- IMSC_AGL -------- %
addpath('IMSC_AGL')
CLU_IMSCAGLresult = IMSC_AGL_clustering(Dataname,percentDel,f)

% -------- PIC --------- %
addpath('PIC')
CLU_PICresult = PIC_clustering(Dataname,percentDel,f)

end