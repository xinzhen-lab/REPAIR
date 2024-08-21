function [Vttend] = Vtest(motai,mshu,tt_X,Vquan,k,Q,T,P)
[m,tt_n]=size(tt_X);  
if strcmp(Vquan,'Vtyou')     
    Vbing=zeros(k*motai,tt_n);                                      
    Vttsum=zeros(k,tt_n);                               
    for g=1:motai        
        [hang1,hang2] = motaishu(g,mshu);       
        Vtt = P(1:end,hang1:hang2)*tt_X(hang1:hang2,1:end);                                               
        Vbing(1+(g-1)*k:g*k,1:end) =Vtt;                                              
        Vttsum=Vttsum+Vtt;                                     
    end   
    Vttend=zeros(k,tt_n);
    for g=1:tt_n                                            
        youtai=motai;                                                  
        for i=1:motai                                                                
            if Vbing(1+(i-1)*k:i*k,g)==0                                                                            
                youtai=youtai-1;                                                              
            end            
        end        
        Vttend(1:end,g)=Vttsum(1:end,g).*(motai/youtai);
    end    
    
elseif strcmp(Vquan,'Vjun')
    Vbing=zeros(k*motai,tt_n);          
    Vbingyuan=zeros(k*motai,tr_n);    
    for g=1:motai      
        [hang1,hang2] = motaishu(g,mshu);        
        Xend=tr_X1*Q.*T+tr_X1;        
        Vtr = P(1:end,hang1:hang2)*Xend(hang1:hang2,1:end);           
        Vtt = P(1:end,hang1:hang2)*tt_X(hang1:hang2,1:end);                                               
        Vbing(1+(g-1)*k:g*k,1:end) =Vtt;          
        Vbingyuan(1+(g-1)*k:g*k,1:end) =Vtr;        
    end            
    VV=sum(Vbingyuan,2)/tr_n;    
    for g=1:tt_n                                                                                          
        for i=1:motai                                                               
            if Vbing(1+(i-1)*k:i*k,g)==0                
                Vbing(1+(i-1)*k:i*k,g)=VV(1+(i-1)*k:i*k,1);                
            end            
        end        
    end    
    Vttend=zeros(k,tt_n);    
    for g=1:motai    
        Vttend=Vttend+Vbing(1+(g-1)*k:g*k,1:end);        
    end   
end
end
