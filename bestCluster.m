function [indx centroid bound] = bestCluster(Dic, clusterSize)

indx = zeros(size(Dic,2),1);

DisMatrix = zeros(size(Dic,2));

for i = 1:size(DisMatrix,1)
    DisMatrix(i,i:size(DisMatrix,1)) = distan(repmat(Dic(:,i),[1 size(DisMatrix,1)-i+1]),Dic(:,i:size(DisMatrix,1)));
end
DisMatrix = DisMatrix + DisMatrix';

[sortDisMatrix sortI] = sort(DisMatrix);

[dummy minCluster] = min(sum(sortDisMatrix(1:clusterSize,:)));

in = sortI(1:clusterSize,minCluster);

centroid = mean(Dic(:,in),2);
bound = max(distan(Dic(:,in),repmat(centroid(:,end),[1 length(in)])));
indx(in) = 1;