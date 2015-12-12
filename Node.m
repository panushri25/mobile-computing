classdef Node
    properties
       numBranches = 0;
       dicIndx = [];
       centroids = []
       disMat = [];
       bounds = [];
       branches = {};
    end
    methods
        function newNode = Node(newCentroids, newBounds, newIndx, newNumBranches)
            newNode.centroids = newCentroids;
            newNode.bounds = newBounds;
            newNode.disMat = zeros(size(newCentroids,2),size(newCentroids,2));
            for i = 1:size(newCentroids,2)-1
                newNode.disMat(i+1:end,i) = distan(newCentroids(:,i+1:end),repmat(newCentroids(:,i),[1 size(newCentroids,2)-i]));
            end
            newNode.disMat = newNode.disMat + newNode.disMat';
            if nargin >2
                newNode.dicIndx = newIndx;
                if nargin==5
                    newNode.branches = cell(newNumBranches,1);
                end
            end
        end
        function obj = addBranch(obj,newBranch)
            obj.branches(obj.numBranches+1,1) = {newBranch};
            obj.numBranches = obj.numBranches + 1;
        end
    end
end