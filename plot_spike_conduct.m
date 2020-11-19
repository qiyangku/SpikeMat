%%
clear
clc
close all

load data_diffapprox/para_search_conduct.mat
%indx = (Rmax-Rmin)<1; %stable uniform solution
%R2=R;
%R2(~indx)=NaN;

subplot(1,2,1)
imagesc(ae,ai,r0'*1e3,[0 5])
xlabel('a_e');
ylabel('a_i');
set(gca,'ydir','normal');
colorbar

load data_spike/combined_dataset02.mat
subplot(1,2,2)

imagesc(Ae,Ai,meanR',[0 1]);
colorbar
xlabel('a_e');
ylabel('a_i');

set(gca,'ydir','normal');
set(gca,'color','w')

