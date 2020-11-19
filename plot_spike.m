%%
%load combined_dataset1.mat
clear
close all

dataset=1;

%L=load(['./dataset' num2str(dataset) '/data_trial1.mat']);
L1=load('data_spike/combined_dataset01.mat');
L2=load('data_diffapprox/para_search_current');

D=[1 7]; %subplot dimension %BUG for D=[1 8+]???

m=length(L1.Ae)*length(L1.Ai);

figure
[pos,H,W]=subplotmod(D,1);
close all
figure
a = 120; %conversion factor from natural units to pixel
set(gcf,'position',[100 100 a*W a*H],'color','w')


for trial=1:m
    
    [i,j]=ind2sub([length(L1.Ae) length(L1.Ai)],trial);        
    
    if L1.meanR(i,j)<60 && L1.meanR(i,j)>0.01%only plot if mean firing rate is reasonable
    clc

    load(['./data_spike/dataset' num2str(dataset,'%02u') '/data_trial' num2str(trial,'%03u') '.mat'],'Spk');
    subplot(D(1),D(2),1)
    imagesc(L1.Ai,L1.Ae,L1.meanR);
    hold on
    plot(L1.Ai(j),L1.Ae(i),'.g','markersize',20)
    hold off
    xlabel('a_i');ylabel('a_e');
    subplotmod(D,1);
    
    
    subplot(D(1),D(2),2)
    imagesc(L2.ai,L2.ae,L2.maxR,[0 40]);
    %axis([Ai(1) Ai(end) Ae(1) Ae(end)])
    hold on
    plot(L2.ai(j),L2.ae(i),'.g','markersize',20)
    hold off
    set(gca,'ytick',[]);
    subplotmod(D,2);
               
    subplot(D(1),D(2),3)
    plot(Spk(:,1),Spk(:,2),'.',[1 1],[-pi pi],'r','linewidth',2)
    xlim([0 5]) %limit to T<5 s otherwise too many points - can't see!
    ylim([-pi  pi])
    title(['trial = ' num2str(trial) ', meanR = ' num2str(L1.meanR(i,j)) ' Hz, CV=' ...
        num2str(L1.CV(i,j)) ', ae = ' num2str(L1.Ae(i)) ', ai = ' num2str(L1.Ai(j))])
    [pos,H,W]=subplotmod(D,3);
    set(gca,'position',axiscrd([pos(1)+0.2 pos(2)+0.2 0.6 5.8],H,W)); %subplot across multiple fields
    xlabel('t (s)')
    ylabel('x');
    
    pause
    
    %subplot(1,2,2)
    title('Loading... Please wait...');%set(gca,'position',[h/w 0.25 (w-1.2*h)/w 0.5])
    drawnow
    
    end
    
end


