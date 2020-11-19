%do parameter search of fixedpt.m
%identify (WE,WI) range where firing rate who
%doesn't blow up or equal to trivial solution
clear
close all



dataset=2;


switch dataset
    case 1
p=parameters([],'Vlb',-150,'nV',2001);
ae = linspace(0.01,0.05,21);
ai = linspace(0.1,0.5,21);

    case 2
        p=parameters([],'model_type','conduct');
ae = linspace(0,0.05,21);
ai = linspace(0.1,0.3,21);
end



T=1e3;

%R = zeros(length(WE),length(WI),N);
maxdR = zeros(length(ae),length(ai));
maxR = maxdR;
minR = maxR;
meanR=maxR;
test_conv=maxR; %record convergence of probability
%U = R;
%S = R;
tic
for i=1:length(ae)
    for j=1:length(ai)
        %[r,flags]=fixedpt(WE(i),WI(j),N,K,I0);
        %%>>>> ERROR DETECTED 27/2/2017>>>>
        %[rt,opt]=fixedpt_ode23(WE(i),WI(i),N,K,I0,T,[]);
        %<<<< indext should be WI(j) <<<<
        [R,opt]=fixedpt_ode23(ae(i),ai(j),T,[],p);
        maxdR(i,j)=max(abs(opt.r-opt.fp));
        maxR(i,j)=max(opt.r);
        minR(i,j)=min(opt.r);
        meanR(i,j)=mean(opt.r);
        test_conv(i,j)=max(opt.prb.p(:,1)); %record maximum P(Vlb) among all neurons
        %size(r)
    
        %R(i,j,:)=opt.re;
        %U(i,j,:)=opt.ue;
        %S(i,j,:)=opt.se;
        clc
        toc
        progress = [i j]
    end 
end

%save(['./output/para_search.mat'])
