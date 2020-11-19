function [x,R,dr]=fixedpt_rootfinding(~)
%fine fixed point using optimization algorithm
p=parameters([]);
ae=0.015;
ai=0.3;


x=linspace(-pi,pi,101)';
dx=x(2)-x(1);

ke = p.K*vonMises(x,p.de)*dx;
ki = p.K*vonMises(x,p.di)*dx;

%initial condition
r0=50*exp(-x.*x/2/0.2/0.2); %very skinny initial condition
r0(abs(x)>0.5)=0;


%LB=zeros(size(r0));    %firing rates must be positive
LB=0.1*r0;             %ignore uniform trivial solution
UB=100*ones(size(r0)); %exclude firing rates that are too high


options = optimset('MaxIter',1e3*length(x),'MaxFunEvals',1e5*length(x));


tic
%[R,dr]=fminsearch(@(r) objfun(r,ke,ki,ae,ai,p),r0);
[R,dr]=fminsearchbnd(@(r) objfun(r,ke,ki,ae,ai,p),r0,LB,UB,options);
%[R,dr]=fsolve(@(r) objfun(r,ke,ki,ae,ai,p),r0);
toc

plot(x,R)
dr


function dr=objfun(r,ke,ki,ae,ai,p)
Re=convnfft_c(ke,r);
Ri=convnfft_c(ki,r);
[fp,prb]=fokkerplanck(Re,Ri,ae,ai,p);
dr = r-fp;
dr = dr'*dr;