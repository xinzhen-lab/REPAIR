clear
clc
warning off;

path = './';
addpath(genpath(path));
dataName = 'LGG 1p19g1_nxn'; %%% flower17; flower102; CCV; caltech101_numofbasekernel_10
%% %% washington; wisconsin; texas; cornell
load([path,'datasets/',dataName,'_Kmatrix'],'KH','Y');
kk1=KH(:,:,1);
% KH1=KH(:,:,1);
% kH2=KH(:,:,3);
Y=abs(Y-2);
% load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numclass = length(unique(Y));%多少类
numker = size(KH,3);%多少模态？
num = size(KH,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KH = kcenter(KH);
KH = knorm(KH);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qnorm = 2;
% [H_normalized,gamma,obj] = mkkmeans_train(KH,numclass,qnorm);
% res_gnd = myNMIACC(H_normalized,Y,numclass);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsionset = [0.75];
% ie =1:length(epsionset)
for ie=1
    for iter = 1:10
        load([path,'missingratio/',dataName,'_missingRatio_',num2str(epsionset(ie)),...
            '_missingIndex_iter_',num2str(iter),'.mat'],'S');
        
% %         %%%%%%%%%%%--Zero-Filling--%%%%%%%%%%%%%%%%%%%%%%%%%%
%         tic
%         KH1 = algorithm2(KH,S);
%         [H_normalized1,gamma1,obj1] = mkkmeans_train(KH1,numclass,qnorm);
%         timingcost(1) = toc;
%         res(:,1) = myNMIACC(H_normalized1,Y,numclass);
%         
% %         %%%%%%%%%%%--mean-Filling--%%%%%%%%%%%%%%%%%%%%%%%%%%
%         tic;
%         KH2 = algorithm3(KH,S);
%         [H_normalized2,gamma2,obj2] = mkkmeans_train(KH2,numclass,qnorm);
%         timingcost(2) = toc;
%         res(:,2) = myNMIACC(H_normalized2,Y,numclass);
%         
%         
%         %%%%%%%%%%--knn-Filling--%%%%%%%%%%%%%%%%%%%%%%%%%%
%         tic;
%         KH3 = algorithm0(KH,S,7);
%         [H_normalized3,gamma3,obj3] = mkkmeans_train(KH3,numclass,qnorm);
%         timingcost(3) = toc;
%         res(:,3) = myNMIACC(H_normalized3,Y,numclass);
%         
%         % %         %%%%%%%%%%%---EM-filling---%%%%%%%%%%%%%%%%%%%%%%%
%         %         KH4 = algorithm6(KH,S);
%         %         [H_normalized4,gamma4,obj4] = mkkmeans_train(KH4,numclass,qnorm);
%         %         res(:,4) = myNMIACC(H_normalized4,Y,numclass);
%         %%%%%%%%%--Laplacian-filling----%%%%%%%%%%%%%%%%%%%%%%
%         tic;
%         alpha04 = 1e-3;
%         KH4 = algorithm4(KH,S,numclass,alpha04);
%         [H_normalized4,gamma4,obj4] = mkkmeans_train(KH4,numclass,qnorm);
%         timingcost(4) = toc;
%         res(:,4) = myNMIACC(H_normalized4,Y,numclass);
% %         
%         %%%%%%%%---Average---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        algorithm_choose1 = 'algorithm2';%zero
        [H_normalized1,gamma1,obj1,KH1] = myabsentmultikernelclustering(KH,S,numclass,qnorm,algorithm_choose1);
        kk=KH1(:,:,1);
        timingcost(1) = toc;
        res(:,1) = myNMIACC(H_normalized1,Y,numclass);
        
        tic;
        algorithm_choose2 = 'algorithm3';%mean
        [H_normalized2,gamma2,obj2,KH2] = myabsentmultikernelclustering(KH,S,numclass,qnorm,algorithm_choose2);
        timingcost(2) = toc;
        res(:,2) = myNMIACC(H_normalized2,Y,numclass);
        
        tic;
        algorithm_choose3 = 'algorithm0';%KNN
        [H_normalized3,gamma3,obj3,KH3] = myabsentmultikernelclustering(KH,S,numclass,qnorm,algorithm_choose3);
        timingcost(3) = toc;
        res(:,3) = myNMIACC(H_normalized3,Y,numclass);
        
% %         %%%%%%%%---AAAI18---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%         tic;
%         algorithm_choose1 = 'algorithm0';
%         lambdaset3 = 2.^[-15:2:15];
%         tauset3 = [0.2:0.1:0.8];
%         accval3 = zeros(length(tauset3),length(lambdaset3));
%         nmival3 = zeros(length(tauset3),length(lambdaset3));
%         purval3 = zeros(length(tauset3),length(lambdaset3));
%         for it =1:length(tauset3)
%             numSel = round(tauset3(it)*num);
%             A3 = genarateNeighborhood(avgKer,numSel);
%             HE3 = calHessian(KH,A3);
%             for il =1:length(lambdaset3)
%                 [H_normalized4,gamma4,obj4,KH4] = myabsentlocalizedmultikernelclustering(KH,HE3,A3,...
%                     numclass,lambdaset3(il),numSel,S,qnorm,algorithm_choose1);
%                 res3 = myNMIACC(H_normalized3,Y,numclass);
%                 accval3(it,il) = res3(1);
%                 nmival3(it,il) = res3(2);
%                 purval3(it,il) = res3(3);
%             end
%         end            
%         
%         [H_normalized5,gamma5,obj5,KH5] = myabsentmultikernelclustering(KH,S,numclass,qnorm,algorithm_choose1);
%         timingcost(5) = toc;
%         res(:,8) = myNMIACC(H_normalized5,Y,numclass);
%                 
        
% % %         %%%%%%%%---IJCAI2017-----%%%%%%%%%%%%%%%%
%         algorithm_choose8 = 'algorithm3';
%         lambdaset8 = 2.^[-15:2:11];
%         accval8 = zeros(length(lambdaset8),1);
%         nmival8 = zeros(length(lambdaset8),1);
%         purval8 = zeros(length(lambdaset8),1);
%         algval8 = zeros(length(lambdaset8),1);
%         for il =1:length(lambdaset8)
%             tic;
%             [H_normalized8,gamma8,obj8,KH8] = myamkcwithlambda(KH,S,numclass,qnorm,algorithm_choose8,lambdaset8(il));
%             timingcost(8) = toc;
%             res8 = myNMIACC(H_normalized8,Y,numclass);
%             accval8(il) = res8(1);
%             nmival8(il) = res8(2);
%             purval8(il) = res8(3);
%             algval8(il) = calKernelAlignment(KH,KH8)'*gamma8;
%         end
%         res(:,8) = [max(accval8); max(nmival8); max(purval8)];
% %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        alignment(1) = calKernelAlignment(KH,KH1)'*gamma1; 
        alignment(2) = calKernelAlignment(KH,KH2)'*gamma2; 
        alignment(3) = calKernelAlignment(KH,KH3)'*gamma3;
%         alignment(4) = calKernelAlignment(KH,KH4)'*gamma4;
%         alignment(1) = calKernelAlignment(KH,KH5)'*gamma5;
%         alignment(2) = calKernelAlignment(KH,KH6)'*gamma6;
%         alignment(3) = calKernelAlignment(KH,KH7)'*gamma7;
%         alignment(8) = max(algval8);
%         
        save([path,'myRes/',dataName,'_missingRatio_',num2str(epsionset(ie)),'_norm_',num2str(qnorm),...
            '_clustering_iter_',num2str(iter),'.mat'],'res','timingcost');
    end
end