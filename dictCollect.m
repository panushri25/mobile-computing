clc
close all
clear

dictPath = '..\lightFieldDicForAshok\';

names = ls(dictPath);

dict = [];
for i= 3:size(names,1);
    load([dictPath names(i,:)]);
    dict = [dict D];
end

save('Dictionaries/LF_16x16x5x5_146k')