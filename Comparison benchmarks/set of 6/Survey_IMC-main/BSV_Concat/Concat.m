function [Clu_result] = Concat(X,ind_folds,ground_label,numClust)
% If you use the code, please cite the following references:
% [1] Jie Wen, Zheng Zhang, Lunke Fei, Bob Zhang, Yong Xu, Zhao Zhang, Jinxing Li, A Survey on Incomplete Multi-view Clustering, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS: SYSTEMS, 2022.
% [2] Wen J, Xu Y, Liu H. Incomplete multiview spectral clustering with adaptive graph learning[J]. IEEE Transactions on Cybernetics, 2020, 50(4): 1418-1429.
% The following released codes are written  by Jie Wen. For any problems, contact jiewen_pr@126.com
 
new_X = [];
for iv = 1:length(X)
    X1 = X{iv}; 
    X1 = NormalizeFea(X1,0);
    ind_1 = find(ind_folds(:,iv) == 1);
    linshi_X = X1(:,ind_1);
    mean_linshi_X = mean(linshi_X,2);      

    ind_0 = find(ind_folds(:,iv) == 0);
    X1(:,ind_0) = repmat(mean_linshi_X,1,length(ind_0));
    new_X = [new_X;X1];
end
% csvwrite('F:\shiyan\interpretability\Concat_GBM7_X.csv',new_X);
rand('seed',7878)
pre_labels    = kmeans(new_X',numClust);
result_LatLRR = ClusteringMeasure(ground_label, pre_labels);   
% Clu_result.AUC = result_LatLRR(4);
Clu_result.AC    = result_LatLRR(1);
Clu_result.SPE = result_LatLRR(2);
Clu_result.SEN = result_LatLRR(3); 
Clu_result.NMI = result_LatLRR(5);
Clu_result.Purity= result_LatLRR(4);
Clu_result.Fscore= result_LatLRR(6); 
Clu_result.Bacc=(Clu_result.SPE+Clu_result.SEN)/2;
end