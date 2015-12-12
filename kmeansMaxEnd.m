function [indx C D] = kmeansMaxEnd(Dic,numClus1)

[indx C dummy D] = kmeans(Dic',numClus1,'distance','cosine','emptyaction','singleton');
s = sparse(1:length(indx),indx,1,length(indx),size(D,2),length(indx));
clusteredD = D.*s + 1-s;
[dummy in] = min(clusteredD);
C = Dic(:,in)';

for i = 1:numClus1
    D(:,i) = distan(Dic,repmat(C(i,:)',[1 size(Dic,2)]));
end

[dummy indx] = min(D,[],2);

