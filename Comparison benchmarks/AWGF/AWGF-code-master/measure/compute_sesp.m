function [Sensitivity, Specificity] = compute_sesp(y_actual, y_pred)
    TP = 0;
    FP = 0;
    TN = 0;
    FN = 0;
    biaoqian=unique(y_actual);
    if sum(biaoqian)~=1
        y_actual=-y_actual+2;
        y_pred=-y_pred+2;
    end
    
    for i =1:size(y_pred,2)
        if y_actual(1,i) ==1 &&y_pred(1,i) == 1
            TP =TP+ 1;
        end
    end
    for i =1:size(y_pred,2)
        if y_actual(1,i) == 0 &&  y_pred(1,i)==1
            FP =FP+ 1;
        end
    end
    
    for i =1:size(y_pred,2)
        if y_actual(1,i)  ==0&& y_pred(1,i)  == 0
            TN =TN + 1;
        end
    end
    for i  =1:size(y_pred,2)
        if y_actual(1,i) == 1 &&  y_pred(1,i)== 0
            FN = FN +1;
        end
    end    
        Sensitivity = TP / (TP + FN);
        Specificity = TN / (TN + FP);

    

end