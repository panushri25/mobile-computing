function [coef numProd] = bbmatchOO(Dic,tree,numClus, X, a, numProd)

testToCentroid = zeros(size(X,2),numClus(1));
for i = 1:numClus(1)
    testToCentroid(:,i) = distan(X,repmat(tree.centroids(:,i),[1 size(X,2)]))';
end

numProd = numProd + numClus(1)*size(X,2);

b = repmat(tree.bounds(1:numClus(1))',[size(X,2) 1]);
upperBound = testToCentroid.*sqrt(1-b.^2) + b.*sqrt(1-testToCentroid.^2); 
upperBound(asin(testToCentroid)+asin(b)> (pi/2)) = 1;
lowerBound = max(testToCentroid.*sqrt(1-b.^2) - b.*sqrt(1-testToCentroid.^2),0);

selectedX = repmat(min(upperBound,[],2),[1 numClus(1)])>=lowerBound;


if a~=1
    for k = 1:size(X,2)
        selectedClustersK = find(selectedX(k,:)==1);
        if length(selectedClustersK)>floor(a*numClus(1))
            [val, s] = sort(testToCentroid(k,selectedClustersK));
            selectedX(k,selectedClustersK(s(floor(a*numClus(1))+1:end))) = 0; 
        end
    end
end

coefTemp = zeros(size(Dic,2),size(X,2));

for i = 1:numClus(1)
    selectedXI = selectedX(:,i)==1;

    if length(numClus)>1
        if sum(selectedXI)~=0
            [coefTemp(tree.dicIndx==i,selectedXI) numProd] = bbmatchOO(Dic(:,tree.dicIndx==i),tree.branches{i,1},numClus(2:end), X(:,selectedXI), a, numProd);
        end
    else
        if sum(selectedXI)~=0
            coefTemp(tree.dicIndx==i,selectedXI) = myOMP(Dic(:,tree.dicIndx==i),X(:,selectedXI),1);
            numProd = numProd + sum(tree.dicIndx==i)*sum(selectedXI);
        end
    end 
end


[val ind] = max(abs(coefTemp),[],1);

coefI = sparse(ind,1:size(X,2),ones(size(val)),size(Dic,2),size(X,2),size(X,2));
coef = coefI.*coefTemp;