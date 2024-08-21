function [W,T,T1,T2,T3,T4] = qnnW(qushu,motai,m,mshu,X,tr_n,tr_Y,qudistan)
quekong=0;%X Number of missing columns
indkong=[ ];
T=zeros(m,tr_n);%The vacant position corresponding to X is 1
T1=ones(m,tr_n);
biaozhi=zeros(motai,tr_n);
for g=1:motai
    [hang1,hang2] = motaishu(g,mshu);
   for i=1:tr_n 
       if X(hang1:hang2,i)==0
           quekong=quekong+1;
           indkong(quekong)=i;     
           biaozhi(g,i)=1;            
           T(hang1:hang2,i)=1;
           T1(hang1:hang2,i)=0;          
       end  
   end
end    
    indkong=unique(sort(indkong));
    kongsize=size(indkong,2);
    kongtai=zeros(kongsize,motai);
    W=zeros(tr_n,tr_n);
    T2=zeros(tr_n,tr_n);%缺的对应Q的列
    T3=zeros(tr_n,tr_n);%不计入权重
    T4=zeros(tr_n,tr_n);%计入权重的,同W分布
for i=1:kongsize
        Xou=X;
        yangben=indkong(i);%第几个样本，即第几列
        T2(1:end,yangben)=1;
        T3(1:end,yangben)=T2(1:end,yangben);
        for g=1:motai
            [hang1,hang2] = motaishu(g,mshu);

            if Xou(hang1:hang2,yangben)==0
                Xou(hang1:hang2,1:end)=0;
                kongtai(i,g)=g;
            end
        end
        for g=1:motai
            if kongtai(i,g)~= 0
                for panind=1:tr_n 
                    if biaozhi(kongtai(i,g),panind)==1       
                        Xou(1:end,panind)=Xou(1:end,yangben);         
                    end              
                end      
            end      
        end
        
        distan=zeros(1,tr_n);  
        if strcmp(qudistan,'gaosi')
            for g=1:tr_n 
                distan(1,g)=sqrt(sum((Xou(1:end,g)-Xou(1:end,yangben)).^2));     % 
            end 
            qunan=find(distan==0);   
            distan(1,qunan)=nan;
            distannan=distan;
            distannan(1,:)=nan;
            len_Y=length(unique(tr_Y));
            for leibieshu=0:len_Y %The number of class
                if tr_Y(1,indkong(i))==leibieshu
                    qunan1=find(tr_Y==leibieshu);
                    distannan(1,qunan1)=distan(1,qunan1);
                end               
            end
            
            distan1=sort(distannan);
            qushu1=distan1(1,qushu+1);
            quind=find(distan<qushu1);
            quanzhong=exp(-distan(1,quind)); 
            
        elseif strcmp(qudistan,'oushi')
            fprintf('oushi.\n')
            for g=1:tr_n     
                distan(1,g)=1./sqrt(sum((Xou(1:end,g)-Xou(1:end,yangben)).^2));    
            end 
            distan(isinf (distan))=0;
            distan1=sort(distan,'descend');    
            qushu1=distan1(1,qushu+1);   
            quind=find(distan>qushu1);   
            quanzhong=distan(1,quind);                     
        end
        quanzhonghe=sum(quanzhong);
        W(quind,yangben)=(quanzhong./quanzhonghe)';
        if isnan(W(quind,yangben))
            W(quind,yangben)=0;            
        end
        
        
        T3(quind,yangben)=0;
        T4(quind,yangben)=1;
end
end
