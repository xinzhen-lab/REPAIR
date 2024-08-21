function [L,quanshu,uind] = generate_L(X,tr_Y,Zqushu,motai,tr_n,mshu)
L2=zeros(tr_n,tr_n);
    for g =1:motai
        shu=0;
        uind=zeros(1,tr_n);
        [hang1,hang2] = motaishu(g,mshu);
        for i=1:tr_n           
            if X(hang1:hang2,i)==0
                shu=shu+1;
            end
            if shu==0
                uind(1,i)=1;
            end
            shu=0;
        end        
        quanshu=sum(uind);
        quuind=find(uind==1); 
        Yuind=tr_Y(1,quuind);
        Xu=X(hang1:hang2,quuind);
        
        Zbing=zeros(tr_n,quanshu);
        Zbing1=zeros(tr_n,quanshu);
        Zju=Zbing1;
        deta=5;
        for ix=1:tr_n
            if find(X(hang1:hang2,ix)~=0)                
                for ju=1:quanshu      
                    Zju(ix,ju)=norm(X(hang1:hang2,ix)-Xu(1:end,ju),2);
                    Zbing1(ix,ju)=exp(-(Zju(ix,ju))^2/(deta^2));
                end  
                Zdistan=Zju(ix,1:end);
                Zdistan1=sort(Zdistan);
                Zqushu1=Zdistan1(1,Zqushu);
                Zquind=find(Zdistan>Zqushu1);
                Zbing1(ix,Zquind)=0; 
                if find(Zbing1(ix,1:end)~=0) 
                    Zbing(ix,1:end)=Zbing1(ix,1:end)./sum(Zbing1(ix,1:end),2);
                end
            end            
        end
        a=ones(tr_n,1);
        DJ=diag((Zbing'*a));  
        S=Zbing*(DJ^-1)*Zbing';
        D=diag(sum(S,2));
        L1=D-S;
        L2=L1+L2;        
    end
    L2=L2/motai;
    Stxt=zeros(tr_n,tr_n);
    for ix=1:tr_n
        for iy=1:tr_n
            if tr_Y(1,ix)==tr_Y(1,iy)
                Stxt(ix,iy)=1;
            end
        end
    end
    Dtxt=diag(sum(Stxt,2));
    Ltxt=Dtxt-Stxt;
    Ltxt=Ltxt/600;
    L=L2+Ltxt;
end

