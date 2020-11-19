%do parameter search of fixedpt.m
%identify (WE,WI) range where firing rate who
%doesn't blow up or equal to trivial solution
clear
close all
m=21; %number of parameter values


dataset=3;


switch dataset
    case 1
p=parameters([],'Vlb',-150,'nV',2001);
WE = linspace(0.02,0.04,m);
WI = 0.3;
    case 2
        p=parameters([],'model_type','conduct');
WE = linspace(0.02,0.024,m);
WI = 0.2;
    case 3
        p=parameters([],'model_type','conduct');
WE = linspace(0.0211,0.0212,m);
WI = 0.2;
    case 4 %replicate Yifan's simulation
        p=parameters([],'model_type','conduct','separateEIpop',true);
wee = 0.5*p.aee;%0.042;%linspace(0.01,0.06,m);
wie = p.aie;%0.05;
%wei = 0.1;
g = linspace(0.1,3.1,m);
WEI = p.aei*g;
wii = p.aii;
    case 5 %only 2d no separateEIpop
        p=parameters([],'model_type','conduct');
        WE = linspace(0.0262,0.0267,m);
        WI=  0.2;
    case 31
        p=parameters('yifan');
        wee = 0.4775*p.aee;%0.042;%linspace(0.01,0.06,m);
        wie = p.aie;%0.05;
        %wei = 0.1;
        g = linspace(0.5,1.5,m);
        WEI = p.aei*g;
        wii = p.aii;
end

%%
T=1e3;%1e4; %ms; simulation time

%NI=NE;
%I0=500; %[strength rate (Hz)]constant baseline current %fluctuation driven needs this < 15

switch p.dim
    case 1
R = zeros(p.nx,m);
    case 2
R = zeros(p.nx*p.nx,m);        
end

FP =R;
%RT = zeros(length(WE),1e4);
%U = R;
%S = R;
%FP = R;
maxdR = zeros(1,m);
maxR = maxdR;
minR = maxR;
meanR = maxR;
WHM=maxdR; %width at half maximum
totR = maxdR; %total firing rate
UI = R; %mean total synaptic current
VI = UI; %var total synaptic current
meanV = maxR;
stdV = maxR;
skV = maxR;

tic
for i=1:m
    %for j=1:length(WI)
        if p.separateEIpop
                %NB [aee aei aie aii]
                %[r0,opt]=fixedpt_ode23([WE(i) wie],[wei wii],T,[],p);
                [r0,opt]=fixedpt_ode23([wee wie],[WEI(i) wii],T,[],p);
                n = length(opt.r)/2;
                opt.r = opt.r(1:n); %only look at pop E for now
                opt.fp = opt.fp(1:n);
        else
                [r0,opt]=fixedpt_ode23(WE(i),WI,T,[],p);
                
        end
                
        maxdR(i)=max(abs(opt.r-opt.fp));
        maxr=max(abs(opt.r));
        maxR(i)=maxr;
        minr=min(abs(opt.r));
        minR(i)=minr;
        meanR(i)=mean(opt.r);
        %if maxr>1e4
        %    break
        %end
        
        v = opt.prb.V;
        dV=opt.prb.V(2)-opt.prb.V(1);
        p0 = mean(opt.prb.p); %mean of histograms
        meanV(i) = sum(v.*p0*dV); %population mean voltage
        stdV(i) = sqrt(sum( (v-meanV(i)).^2.*p0*dV)); %std of popualtion V
        skV(i) = sum(  ((v-meanV(i))./stdV(i) ).^3.*p0*dV); %skewness        
               
        
        
        indx = find(opt.r>((maxr-minr)/2+minr));
        if isempty(indx)
            WHM(i)=0;
        else
        switch p.dim
            case 1
        WHM(i)=opt.x(indx(end))-opt.x(indx(1));
            case 2
                [indxi,indxj]=ind2sub(indx,[p.nx p.nx]);
                %indxj = indxj(indxi == floor(p.nx/2));
                WHM(i)=max(opt.x(indxj))-min(opt.x(indxj));
        end
        end
        
%         dx = opt.x(2)-opt.x(1);
%         switch p.dim
%             case 1
%         totR(i)=sum(opt.r)*dx;
%             case 2
%         totR(i)=sum(opt.r)*dx;
%         end
        
        %size(r)
        R(:,i)=opt.r;
        FP(:,i)=opt.fp;
        %UI(:,i)=opt.prb.utot;
        %VI(:,i)=opt.prb.vtot;
        
        clc
        toc
        progress = i
    %end
end

return

%% Compute uniform fixed point solution

if p.separateEIpop
    r0e=nan(3,m);    
    r0i=r0e;
    %test_conv=r0; %test convergence at Vlb
    for i=1:m
        %[rr,prb]=fixedpt_uniform([WE(i) wie], [wei wii],p);
        [C,prb]=fixedpt_uniform([wee wie], [WEI(i) wii],p);
        [ii,jj]=size(C); %ii=number of fixed poitns; jj =2
        ii = min(ii,3);
        r0e(1:ii,i)=C(1:ii,1);
        r0i(1:ii,i)=C(1:ii,2);
        %test_conv(i)=prb.p(1)/max(prb.p);
    end
else    
    r0=zeros(size(WE));
    test_conv=r0; %test convergence at Vlb
    for i=1:length(WE)
        [r0(i),prb]=fixedpt_uniform(WE(i),WI,p);
        test_conv(i)=prb.p(1)/max(prb.p);
    end
end

mybeep
return
%%
close all
subplot(3,2,1)


subplot(3,2,2)
plot(WE,maxdR*1e3,'.');
xlabel('a_e');ylabel('Error (Hz)')
axis([WE(1) WE(end) 0 0.01])

subplot(3,2,3)
%indx=maxR<100;
%plot(WE,r0*1e3,WE,minR*1e3,'.',WE,maxR*1e3,'.',WE,meanR*1e3);
plot(WE,r0*1e3,WE,minR*1e3,'.',WE,maxR*1e3,'.');
xlabel('a_e');ylabel('Firing rate (Hz)')
axis([WE(1) WE(end) 0 10])
legend('Uniform','Trough','Peak')

subplot(3,2,4)
%xx=linspace(0.024,0.027,51);
%yy=pi-exp(3.322)*sqrt(xx-0.024);
yy=linspace(1,pi,31);
xx=(pi-yy)*exp(-3.35);xx=xx.*xx+0.0241;
plot(WE,WHM,'.',xx,yy);
hold on
xlabel('a_e');ylabel('Bump width');
axis([WE(1) WE(end) 0 pi])


subplot(3,2,5)
indx=45;%25;

plot(opt.x,R(:,indx)*1e3,opt.x(1:20:end),FP(1:20:end,indx)*1e3,'.')
xlim([-pi pi])
xlabel('x');ylabel('Firing rate (Hz)')
title(['Bump at ae=' num2str(WE(indx))])

subplot(3,2,6)
%plot(opt.x,UI(:,indx),opt.x,UI(:,indx)+sqrt(VI(:,indx)));legend('Mean','Mean+Std')
plot(opt.x,UI(:,indx),opt.x,sqrt(VI(:,indx)));legend('Mean','Std')
%plot(opt.x,opt.prb.utot,opt.x,sqrt(opt.prb.vtot))

hold on
plot([-pi pi],(p.Vth-p.El)./p.tau*[1 1],'--k'); %firing threshold (mv)
hold off
xlim([-pi pi])
%ylim([-2 12])
xlabel('x');ylabel('Total synaptic current')
legend('Mean','Std')

set(gcf,'position',[731  79   560   802],'color','w')

%%
%save(['FDB_we_cond' num2str(floor(now)) '_' num2str(rem(now,1)) '.mat'],'FP','R','I0','K','maxdR','maxR','minR','N','WE','WI','opt','totR','WHM','T');