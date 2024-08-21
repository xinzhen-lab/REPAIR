function [functionvalue,gradient] = relevancy_f(f,C)
n=size(f,1);
E=eye(n);
one_1=zeros(n,n);
for i=1:n
    one_1(i,1)=1;
end
e=one_1;
i=ones(1,n);
E1=E./n;
C1=(C*e)*e'.*E;
for j=1:size(C1,1)
    C1(j,j)=1/C1(j,j);
end
T_0=f*i;
T_1=T_0-T_0';
T_2=exp(-abs(T_1));
t_3=log(2);
T_4=e*e';
T_5=T_0'-T_0;
T_6=exp(-abs(T_5));
T_7=T_4*T_6.*E1;
T_8=(T_2*e)*e'.*E1;
t_9=det(T_8);
t_10=n*t_3;
T_11=E.*C1;
T_12=(T_4*(T_6.*C')).*T_11;
T_13=sign(T_1);
T_14=((T_2.*C)*e)*e'.*T_11;
t_15=det(T_14);
T_16=sign(T_5);
a1=((det(T_7)* n)* t_3);
a2=((t_9 .* n) .* t_3);
a3=(t_10 .* t_15);
a4=(t_10 .* det(T_12));

gradient = ((((1 /a1) .* (((((det(T_7) .* inv(T_7)) .* E1))*(e))*(e') .* (T_2 .* T_13))*(i')) - ((1 / a2) .* (((T_4)*(((det(T_8) .* inv(T_8)) .* E1)) .* (T_6 .* T_16)))*(i'))) - (((1 / a4) .* (((((((det(T_12) .*inv(T_12)) .* E) .* C1))*(e))*(e') .* ((C .* T_2) .* T_13)))*(i')) - ((1 / a3) .* (((T_4)*((((det(T_14) .*inv(T_14)) .* E) .* C1)) .* ((C' .* T_6) .* T_16)))*(i'))));
functionvalue=log(t_15/t_9)/t_10;
end   


