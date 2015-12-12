clc
close all
clear

ImgPath = 'C:\Users\Ali\Dropbox\FastSparseRepresentation\Latex\Paper\Figures\SourceImages\Video_SubSamp_GPSR.png';
zoomImgPath = 'C:\Users\Ali\Dropbox\FastSparseRepresentation\Latex\Paper\Figures\SourceImages\Video_SubSamp_GPSR_zoom.png';

patchLoc = [170 280];

%% Video Denoising: patchLoc = [20 50];

Img = imread(ImgPath);
zoomImg = Img(patchLoc(1)+(1:75),patchLoc(2)+(1:120),:);

imshow(zoomImg)
imwrite(zoomImg,zoomImgPath);