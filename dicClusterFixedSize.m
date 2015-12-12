function [indx centroid bound] = dicClusterFixedSize(Dic, numOfClusters, lastIndx)

bound = [];
centroid = [];
indx = zeros(size(Dic,2),1);

clusterSize = ceil(size(Dic,2)/numOfClusters);

[I C D] = kmeansMaxEnd(Dic,numOfClusters);

num = hist(I, 1:numOfClusters);

bigClusters = find(num>=clusterSize);

for i = 1:length(bigClusters)
    [I2 centroid2 bound2] = bestCluster(Dic(:,I==bigClusters(i)), clusterSize);
    centroid = [centroid centroid2];
    bound = [bound; bound2];
    lastIndx = lastIndx+1;
    indx(I==bigClusters(i)) = I2*lastIndx;
end

if sum(indx==0) > clusterSize
    [indx(indx==0) centroid2 bound2] = dicClusterFixedSize(Dic(:,indx==0), numOfClusters-length(bound), lastIndx);
    centroid = [centroid centroid2];
    bound = [bound; bound2];
elseif sum(indx==0)~=0
    centroid(:,end+1) = mean(Dic(:,indx==0),2);
    bound(end+1,1) = max(distan(Dic(:,indx==0),repmat(centroid(:,end),[1 sum(indx==0)])));
    lastIndx = lastIndx+1;
    indx(indx==0) = lastIndx;
end