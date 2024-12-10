function results = spcclust(Z, k, Y)

cluster_iter = 10;
for i = 1:cluster_iter
    D   = diag(sum(Z,2));
    D2 = diag(1./sqrt(diag(D)+eps));
    L   = D2 - Z;
    Lnor = D2*L*D2;
    Lnor = (Lnor + Lnor')/2;
    %L = eye(116)-D2*Z*D2;

    options.tol = 1e-8;
    options.maxit = 30000;
    [V, val] = eigs(Lnor, k, 'sa', options);
    %V = D2*V;
    
    Ypred     = kmeans(V, k, 'replicates',100,'display','off');
    metric    = clusteringMeasure(Y, Ypred);
   f1_score(Y, Ypred);
    
    ACC(i)    = metric(1);
    NMI(i)    = metric(2);
    Fscore(i) = metric(3);
    Sensitivity(i)=metric(4);
    Specificity(i)=metric(5);
    purity(i)=metric(6);
    
end
meanACC     = mean(ACC);
meanNMI     = mean(NMI);
meanFscore  = mean(Fscore);
meanSensitivity = mean(Sensitivity);
meanSpecificity=mean(Specificity);
meanpurity=mean(purity);


results = [meanACC, meanSpecificity, meanSensitivity, meanpurity, meanNMI, meanFscore];

end
