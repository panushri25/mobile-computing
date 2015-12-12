function [indx centroid bound] = bestClusterAprox(Dic, centroid, clusterSize)

indx = zeros(size(Dic,2),1);

DisMatrix = (Dic'*centroid').^2;


[sortDisMatrix sortI] = sort(DisMatrix);


in = sortI(1:clusterSize);

centroid = mean(Dic(:,in),2);
bound = max(distan(Dic(:,in),repmat(centroid(:,end),[1 length(in)])));
indx(in) = 1;