function [pcaU] = generate_u_pca_matrix(U,lable)
U=U';
[n,m]=size(U);
pcaU=zeros(n,m);
for i=1:m
    new_U=U;
    new_U(:,i) = [];
    [pcaU(:,i),~]=PCA(new_U,1);
end 
end

