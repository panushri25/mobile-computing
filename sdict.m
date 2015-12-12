loadPath = 'Dictionaries/Image_16x16_40k';
savePath = 'Dictionaries/dict10k'  ;
savePath1 = 'Dictionaries/ainput75';

load(loadPath);

[dict] = normc(dict);

[sldict] = dict(:,1:10000);
[input] = dict(:, 3001:3075);

save(savePath,'sldict');
save(savePath1,'input');