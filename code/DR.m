function [M,M1] = DR(quanshu,uind,tr_n,tr_Y)
M=1/quanshu.*uind';
M=repmat(M,1,tr_n);   
M1=zeros(tr_n,tr_n);    
m1ind=uind.*tr_Y;    
m0ind=uind-m1ind;    
Y1ind=find(tr_Y==1);   
Y0ind=find(tr_Y==0);    
m11ind=find(m1ind==1);    
m00ind=find(m0ind==1);    
M1(m11ind,Y1ind)=1/sum(m1ind);    
M1(m00ind,Y0ind)=1/sum(m0ind);
end

