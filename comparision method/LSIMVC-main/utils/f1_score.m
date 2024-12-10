function [SEN, SPE] = f1_score(label, predict)
   M = confusionmat(label, predict);
   
   SPE = M(2,2) / (M(2,1) + M(2,2)); 
   SEN = M(1,1) / (M(1,1) + M(1,2)); %SP: TN/(TN+FP)
end
