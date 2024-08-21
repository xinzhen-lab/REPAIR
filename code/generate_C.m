function [C]=generate_C(y)
y=y';
C=generate_u_relation_matrix(y);
[n1,n2]=size(C);
for i=1:n1
    for j=1:n2
        if C(i,j)~=1
            C(i,j)=0;
        end
    end
end        
end
