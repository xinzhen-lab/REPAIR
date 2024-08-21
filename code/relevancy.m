function [loss,U_g] = relevancy(U1,C)
loss=0;
U=U1';
[n,m]=size(U);
U_g=zeros(n,m);
zhang=[];
for i=1:m
    [l,U_g(:,i)]=relevancy_f(U(:,i),C);
    loss=loss+l;
    zhang=[zhang,loss];
end
loss=-loss;
U_g=-U_g';
end

