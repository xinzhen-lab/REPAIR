function [loss,U_g] = redundancy(U1,pcaU)
loss=0;
U=U1';
[n,m]=size(U);
U_g=zeros(n,m);
for i=1:m
    pcau=generate_u_relation_matrix(pcaU(:,i));
    [l,U_g(:,i)]=relevancy_f(U(:,i),pcau);
    loss=loss+l;
end
U_g=U_g';
end

