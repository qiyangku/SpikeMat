function [C,prb]=fixedpt_uniform(ae,ai,p)
%[R,opt]=fixedpt_uniform(ae,ai,T,ini_cond,p)
%Find the fixed point of uniform stationary solution
%INPUT: ae,ai = excitatory/inhibitory synaptic strength
%p: parameters
%OUTPUT: C = uniform fixed points. (#points x E/I pop)
if isempty(p) %load default parameters if no input parameters
    disp('No input parameters!')
    p=parameters([]);
end

%METHOD 1: fsolve. Sometimes may converge to negative solution.
%Negative firing rate corresponds to swapping E/I.
%test speed: 0.08 sec
%options = optimoptions('fsolve','FunctionTolerance',1e-9);
%[X,FVAL,EXITFLAG]=fsolve(@(x) -x+fokkerplanck(p.K*x+p.I0,p.K*x,ae,ai,p),0.01);

%METHOD 2: lsqnonlin. Least square with constraint r >=0.
%Not suitable for complex valued function.
%test speed: 0.04 sec
options = optimoptions('lsqnonlin','FunctionTolerance',1e-16,'OptimalityTolerance',1e-16,'Display', 'off');





if p.separateEIpop
    
        
    %number of differential initial conditions for root finding
    init=logspace(-5,-1,p.ninit); %10^-5 ~ 10^-1 kHz => 0.01~100 Hz
    r=zeros(p.ninit,2);
    dr=r;
    %prb=struct([]);
    %prb(p.ninit).p=NaN; %preallocate memory
    
    for i=1:p.ninit
        %[X,FVAL,RES,EXITFLAG] = lsqnonlin(@(x) myrhsEI_uniform(x,ae,ai,p),[0.6e-3 4e-3],[0.4e-3; 2e-3],[1e-3; 4e-3],options);%init(i)*[1;1]
        [X,FVAL,RES,EXITFLAG] = lsqnonlin(@(x) myrhsEI_uniform(x,ae,ai,p),init(i)*[1;1],[0; 0],[inf; inf],options);%
        r(i,:)=X;        
        dr(i,:)=abs(FVAL); %error of root finding
    end
    
    C=subclust(r,1e-5);  %substractive clustering
    
    %order C by magnitude
    C=sortrows(C,1); %sort fixed pts based on firing rate of excitatory neuron
    
    %evaluate again for probability density (for the largest firing rate)
    %only record distribution of excitatory neurons
    %maxr=max(r);
    maxr=C(end,:);
    [r0,xe] = fokkerplanck(p.KEE*maxr(1),p.KEI*maxr(2),p.Ie,[ae(1) p.see],[ai(1) p.sei],p.jee,p);
    [r0,xi] = fokkerplanck(p.KIE*maxr(1),p.KII*maxr(2),p.Ii,[ae(2) p.sie],[ai(2) p.sii],p.jie,p);
    %dri = fokkerplanck(p.KIE*maxr(1)+p.Ii,p.KEI*maxr(2),ae(2),ai(2),p);
    
    prb.xe=xe;
    prb.xi=xi;
    
    prb.xe.Re=p.KEE*maxr(1);
    prb.xe.Ri=p.KEI*maxr(2);
    prb.xe.Rext = p.Ie;
    prb.xe.aext = p.jee;
    prb.xi.Re=p.KIE*maxr(1);
    prb.xi.Ri=p.KII*maxr(2);
    prb.xi.Rext = p.Ii;
    prb.xi.aext = p.jie;
    
    %all fixed points found
    prb.r0=C;
    prb.dr=dr;
    %prb.p=prb0.p;
    %prb.dp=prb0.dp;
    
    
    prb.exitflag = EXITFLAG;
    
else
    
    %number of differential initial conditions for root finding
    init=logspace(-5,1,p.ninit); %10^-5 ~ 10^-1 kHz => 0.01~100 Hz
    X=zeros(p.ninit,1);
    dr=X;
    
    for i=1:p.ninit
        [X(i),FVAL,RES,EXITFLAG] = lsqnonlin(@(x) -x+fokkerplanck(p.K*x,p.K*x,p.I0,ae,ai,ae,p),init(i),0,inf,options);
        dr(i)=abs(FVAL); %error of root finding
    end
    
    
    %evaluate again for probability density (for the largest firing rate)
    maxr=max(X);
    [r0,prb]=fokkerplanck(p.K*maxr,p.K*maxr,p.I0,ae,ai,ae,p);
    
    C=subclust(X,1e-5);  %substractive clustering
    %order C by magnitude
    C=sort(C); %sort fixed pts based on firing rate of excitatory neuron
    
    prb.Re=p.K*maxr;
    prb.Ri=p.K*maxr;
    prb.exitflag = EXITFLAG;
    prb.Rext = p.I0;
    prb.aext = ae;
    %all fixed points found
    prb.r0=X;
    prb.dr=dr;
    %prb.p=prb0.p;
    %prb.dp=prb0.dp;
    
end


while false %debug only! Verify susceptibility
    %compute dr for different static dR.
    %should be consistent with low frequency limit of linear response
    dR = p.K*maxr*linspace(0.001,0.005,20);  %fraction of perturbation
    dr = zeros(size(dR));
    for i=1:length(dR)
    [r0,prbe]=fokkerplanck(p.K*maxr,p.K*maxr+dR(i),p.I0,ae,ai,ae,p);
    dr(i)=max(r0)-maxr;
    end
    prb.dr = dr;
    prb.dR = dR;
    break
end


function dr = myrhsEI_uniform(r,ae,ai,p)
%ae = [aee aie]; ai = [aei aii];
re=r(1);ri=r(2);
dre = -re+fokkerplanck(p.KEE*re,p.KEI*ri,p.Ie,[ae(1) p.see],[ai(1) p.sei],p.jee,p);
dri = -ri+fokkerplanck(p.KIE*re,p.KII*ri,p.Ii,[ae(2) p.sie],[ai(2) p.sii],p.jie,p);
dr = [dre;dri];
%dr = log(dre.*dre+dri.*dri);


