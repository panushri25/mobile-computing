function tree = dicClusterDevOO(tree,Dic,numOfClusters)
    
if length(numOfClusters)>1
    if tree.numBranches == numOfClusters(1)
        for i = 1:numOfClusters(1)
            i
            tree.branches(i) = {dicClusterDevOO(tree.branches{i,1},Dic(:,tree.dicIndx==i),numOfClusters(2:end))};
        end
    else
        for i = 1:numOfClusters(1)
            i
            newBranch = dicClusterOO(Dic(:,tree.dicIndx==i),numOfClusters(2:end));
            tree = addBranch(tree,newBranch);
        end
    end
end