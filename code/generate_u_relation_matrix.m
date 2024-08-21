function [juzhen] = generate_u_relation_matrix(u)
n=size(u,1);
one=ones(1,n);
u_one=u*one;
juzhen=exp(-abs(u_one-u_one'));
end

