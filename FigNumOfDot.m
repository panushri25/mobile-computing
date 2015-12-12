clc
close all
clear

a = [0.11 0.125 0.15 0.175 0.2];

lvls = 10;

N(1) = 100000;

c = 100;

num = zeros(lvls, length(a));

for i = 1:lvls
    for j = 1:length(a)
        num(i,j) = numOfDot(a(j),c,N(i));
    end
    c = [c 10];
    N(i+1) = N(i)*10;
end

num = [num N(1:end-1)'];

figure(1)
plot(repmat(N(1:end-1)',[1 length(a)]),num(:,1:end-1),'LineWidth',2)
% hold on
% plot(N(1:5),num(1:5,end))
% hold off
title('Number of products sublinearity');
xlabel('Number of Dictionary atoms');
ylabel('Number of products');
legend('\alpha = 0.11','\alpha = 0.125','\alpha = 0.15','\alpha = 0.175','\alpha = 0.2','Location','NorthWest');


load('C:\Users\Ali\Documents\My Version\ObjOrt\Results\LF\Denoising\LF01\results_BB.mat')
figure(2)
[AX,H1,H2] = plotyy(a,timeBB(:,end)/60,a,snrBB(:,end));
xlabel('\alpha');
title('Time and Performance vs. \alpha');
set(H1,'LineWidth',2);
set(H2,'LineWidth',2);
set(get(AX(1),'Ylabel'),'String','time(min)') 
set(get(AX(2),'Ylabel'),'String','SNR(dB)') 
text(0.01,45,'Time for OMP = 1 day and 12 hrs','FontSize',12);