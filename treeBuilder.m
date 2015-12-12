clc
%%close all
clear

loadPath = 'Dictionaries/Image_16x16_40k';
savePath = 'Dictionaries/tree';
savePath2 = 'Dictionaries/100x10x10';

load(loadPath);

numOfClusters = [100 10 10]; 
%[sldict] = normc(sldict);

tic
[treeClus streeClus]  = dicClusterOO(dict,numOfClusters);
%[treeClus streeClus] = dicClusterDevOO(treeClus,sldict,numOfClusters);
treeTime = toc;
%s = struct(treeClus);
save(savePath,'treeClus','treeTime','numOfClusters');
save(savePath2,'streeClus');