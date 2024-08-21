function [ FinalData,reconData ] = PCA( dataSet, k )
    [m,n] = size(dataSet);

   %% Subtract the average
    dataSetMean = mean(dataSet);
    dataSetAdjust = zeros(m,n);
    for i = 1 : m
        dataSetAdjust(i , :) = dataSet(i , :) - dataSetMean;
    end

    %% Calculate the covariance matrix
    dataCov = cov(dataSetAdjust);

    %% The eigenvalues and eigenvectors of the covariance matrix are calculated
    [V, D] = eig(dataCov);
    
    d = zeros(1, n);
    for i = 1:n
        d(1,i) = D(i,i);
    end
    
    %% Sort the feature values
    [maxD, index] = sort(d);
    
    %% Select the first k largest feature values
    index_k = index(1, (n-k+1):n);
    V_k = zeros(n,k);
    for i = 1:k
        V_k(:,i) = V(:,index_k(1,i));
    end
    
    %% Convert to the new space
    FinalData = dataSetAdjust*V_k;  
    reconData = FinalData * V_k';
    for i = 1 : m
        reconData(i , :) = reconData(i , :) + dataSetMean;
    end
end
