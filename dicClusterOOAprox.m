function tree = dicClusterOOAprox(Dic,numOfClusters)

[indx, centroid, bound] = dicClusterFixedSizeAprox(Dic, numOfClusters(1),0);

tree = Node(centroid, bound, indx, numOfClusters(1));
if length(numOfClusters)>1
    for i = 1:numOfClusters(1)
        newBranch = dicClusterAprox(Dic(:,indx==i),numOfClusters(2:end));
        tree = addBranch(tree,newBranch);
    end
end