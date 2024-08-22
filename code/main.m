clear;clc;

%% Parameter Settings
iter=20;% Iterations
Wqushu=5;% The nearest neighbors of Q
Zqushu=10; % The number of sample-to-pivot pairs :N_g
xQ=0.01;%  Study step
xV=0.1;
xH=0.1;
xP=0.1;
qudistan='gaosi';   % Gaussian distance
Vquan='Vtyou';   % Select the testing data imputation mode

mshu{1}=109;% The features number of modality
mshu{2}=109;
mshu{3}=109;
mshu{4}=109;
motai=size(mshu,2); %The number of modalities

k=30; %Dimension of latent representation
nfold=5; %5-fold cross-validation

qlambda1=1e-1;%The range is [1e-2,1e-1,1e-0]
qlambda2=1e+1;%The range is [1e+0,1e+1,1e+2]
qlambda3=1e+1;%The range is [1e+0,1e+1,1e+2]
qlambda4=1e+1;%The range is [1e+0,1e+1,1e+2]

lambda0=1e-1;%The range is [1e-2,1e-1,1e+0]
lambda1=1e+0;%The range is [1e-1,1e+0,1e+1]
lambda2=1e+1;%The range is [1e+0,1e+1,1e+2]
lambda3=1e+1;%The range is [1e+0,1e+1,1e+2]
lambda4=1e+2;%The range is [1e+1,1e+2,1e+3]
lambda5=1e+0;%The range is [1e-1,1e+0,1e+1]
lambda6=1e+1;%The range is [1e+0,1e+1,1e+2]
lambda7=1e+2;%The range is [1e+1,1e+2,1e+3]
lambda8=1e+1;%The range is [1e+0,1e+1,1e+2]

%% Procedure£º 5-fold cross-validation
for j=1:nfold          
    fprintf(['nfold=',num2str(j),'\n']);
    
% %  Input: Enter datasets that contain missing modalities and enter labels
%     tr_X: Training data,size: m x n, composite multimodal feature matrix,with feature dimensionality of m and patient sample size of n                                 
%     tr_Y: The label of training data,size:1 x n                                                                  
%     tt_X: Testing data                              
%     tt_Y: The label of testing data

% %  Initialization
%     Q: Size:m x n, Data imputation matrix, with sample size of n
%     V: Size:n x n, Latent representation matrix, with dimensionality of k and sample size of n
%     P: Size:k x n, Projection matrix, with latent shared space dimensionality of k and original feature space dimensionality of m
%     H: Size:m x k, Reconstruction matrix, with original feature space dimensionality of m and latent shared space dimensionality of k
    X=tr_X;                                      
    [m,tr_n] = size(X); 
    
% Multimodal Structural Calibration
    [L,quanshu,uind]=generate_L(X,tr_Y,Zqushu,motai,tr_n,mshu);
    
% Discriminative Regularization    
    [M,M1] = DR(quanshu,uind,tr_n,tr_Y);
    
% Prior weighting matrix W
    [W,T,T1,T2,T3,T4] = qnnW(Wqushu,motai,m,mshu,X,tr_n,tr_Y,qudistan);
 
    
    yanben1=sum(T,2);
    n1=1./yanben1;
    yanben2=tr_n-yanben1;
    n2=1./yanben2;
    a=ones(tr_n,1);
    
    Y = generate_C(tr_Y);
%% Iterative update
    for i=1:iter 
        % The objective function
        f1=qlambda1*norm((((X*Q).*T)*a).*n1-(X*a).*n2,2)^2;
        f2=qlambda2*norm(X*Q.*T1-X,'fro')^2;
        f3=qlambda3*norm(Q.*T2-W,'fro')^2;   
        f4=qlambda6*norm((Q.*T2)'*a-a,'fro')^2;
  
        f5=norm(V-P*(X*Q.*T+X),'fro')^2;
        f6=lambda0*norm(X*Q.*T+X-H*V,'fro')^2; 
        f7=lambda1*norm(V,'fro');
        f8=lambda2*sum(sqrt(sum(P.*P)));
        f9=lambda3*norm(H,'fro')^2;
        f10=lambda4*trace(V*L*V');
        f11=lambda5*norm((V-V*M1)*(V-V*M1)','fro')^2; 
        f12=-lambda6*norm((V*M-V*M1)*(V*M-V*M1)','fro')^2;
       
        [pcaV] = generate_u_pca_matrix(V,tr_Y);
        [f13,dV7] = relevancy(V,Y);% feature-task relevancy
        [f14,dV8] = redundancy(V,pcaV);% intra-feature redundancy 

        f13=lambda7*f13;
        f14=lambda8*f14;  
        ff=f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14;

        % Updating Q
        fprintf(['Updating Q','\n']);  
        dQ1=2*X'*(((((((X*Q).*T)*a).*n1-(X*a).*n2).*n1)*a').*T);
        dQ2=2*X'*(((X*Q).*T1-X).*T1);
        dQ3=2*(Q.*T2-W).*T2;
        dQ4=2*(a*(a'*(Q.*T2)-a')).*T2;       
        dQ5=-2*X'*((P'*(V-P*(X+(X*Q).*T))).*T);
        dQ6=2*X'*((X+(X*Q).*T-H*V).*T);
        dQ=qlambda1*dQ1+qlambda2*dQ2+qlambda3*dQ3+qlambda4*dQ4+dQ5+lambda0*dQ6;
        Q=Q-xQ*dQ;    
        
        % Updating V
        fprintf(['Updating V','\n']);    
        dV1=2*(V-P*(X+(X*Q).*T));
        dV2=lambda0*(-2)*H'*(X+(X*Q).*T-H*V);
        dV3=lambda1*V;
        dV4=2*V*L;
        dV4=lambda4*dV4;
        vm0=V-V*M1;
        vm1=vm0*(vm0')*vm0;
        dV5=(4*vm1-4*vm1*M1');
        dV5=lambda5*dV5;
        vm2=V*M-V*M1;
        vm3=vm2*(vm2')*vm2;
        dV6=(4*vm3*(M'-M1'));
        dV6=-lambda6*dV6;
        dV7=lambda7*dV7;
        dV8=lambda8*dV8;
        dV=dV1+dV2+dV3+dV4+dV5+dV6+dV7+dV8;       
        V=V-xV*dV;
        
        % Updating P
        fprintf(['Updating P','\n']);
        Xend=X*Q.*T+X;
        for g=1:size(P,2)           
            xi        =    Xend(g,:);
            cc         =    V-P*(X*Q.*T+X)+P(:,g)*xi;
            P(:,g)    =    (cc*xi')/(xi*xi')*max(1-lambda2/(2*norm(cc*xi',2)),0);
        end
        
        % Updating H
        fprintf(['Updating H','\n']);        
        dH1=lambda0*(-2)*(X+(X*Q).*T-H*V)*V';
        dH2=lambda3*H*2;
        dH=dH1+dH2;
        H=H-xH*dH;                  
    end    
    
%% Find V of the testing data       
    [Vttend] = Vtest(motai,mshu,tt_X,Vquan,k,Q,T,P);
    V_train = [tr_Y;V];         
    V_test = [tt_Y;Vttend];
    
%% Save V files     
    csvwrite(['.\percent0.2_Vtraintest\',num2str(Vquan),num2str(Vqushu),'_',num2str(iter),num2str(fanshi),num2str(Wqushu),'\V1_train4_',num2str(j),'_',num2str(k),'_',num2str(lambda0),'_',num2str(lambda1),'_',num2str(lambda2),'_',num2str(lambda3),'_',num2str(lambda4),'_',num2str(lambda5),'_',num2str(lambda6),'_',num2str(lambda7),'_',num2str(lambda8),'_',num2str(lambda9),'.csv'],V_train);                                       
    csvwrite(['.\percent0.2_Vtraintest\',num2str(Vquan),num2str(Vqushu),'_',num2str(iter),num2str(fanshi),num2str(Wqushu),'\V1_test4_',num2str(j),'_',num2str(k),'_',num2str(lambda0),'_',num2str(lambda1),'_',num2str(lambda2),'_',num2str(lambda3),'_',num2str(lambda4),'_',num2str(lambda5),'_',num2str(lambda6),'_',num2str(lambda7),'_',num2str(lambda8),'_',num2str(lambda9),'.csv'],V_test);
        
end     
  
disp('end')


