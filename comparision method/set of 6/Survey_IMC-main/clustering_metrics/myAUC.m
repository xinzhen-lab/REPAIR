function [AUC] = myAUC(predict,label)
predict=abs(predict-2);
label=abs(label-2);
%计算AUC,predic为预测值，label为真实值，均为行向量形式
num=size(label,2);%计算样本数目
positive=label(1,1);%将第一个标记为阳性，
positive_num=sum(label==positive);%计算阳性个数
negative_num=num-positive_num;%计算阴性个数
AUC=0;
X=[];
Y=[];
for i=1:num
    temp=(predict>=predict(1,i));
    %在所有实际为阳性的样本中，被正确地判断为阳性的样本比率。
    ture_positive_rate=sum(temp==positive&label==positive)/positive_num;
    %在所有实际为阴性的样本中，被错误地判断为阳性的样本比率。
    false_positive_rate=sum(temp==1&label==0)/negative_num;
    X=[X false_positive_rate];
    Y=[Y ture_positive_rate];
end
Z=sortrows([X;Y]');
num=size(Z,1);
where=[];
for i=1:num-1
    if Z(i+1,1)==Z(i,1)
        where=[where i+1];
    end
end
Z(where,:)=[];
AUC=1-trapz(Z(:,1),Z(:,2));
plot(Z(:,1),Z(:,2))
string = {['AUC=' num2str(AUC)]};
title(string)
