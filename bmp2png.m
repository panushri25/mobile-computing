clc
close all
clear


path = '../ImagesB/';

pathDir = ls(path);


for i = 3:size(pathDir,1)
    [dummy,name,fmt] = fileparts(pathDir(i,:));
    if strcmp(fmt(1:4),'.jpg')
        I = imread([path pathDir(i,:)]);
        imwrite(I,[path name '.png'],'png');
    end
end