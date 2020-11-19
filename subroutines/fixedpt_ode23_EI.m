function [R,opt]=fixedpt_ode23_EI(ae,ai,T,ini_cond)
%Same as fixedpt_ode23.m except E/I populations are separate.
%INPUT: ae,ai = excitatory/inhibitory synaptic strength
%T: simulation time
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%Documentation
%Start with cond. based model
%C*dV/dt = -gL*(V-El) - ge*(V-Ee)*Se - gi*(V-Ei)*Si (1)
%where S(t) is instantaneous firing rate (spike count per unit time)
%Divide C through Eq.1
%dV/dt = -(V-El)/tau - ae*(V-Ee)*Se - ai*(V-Ei)*Si (2)

%Case 1: Full spiking neuron
%S(t) = Poisson spike train. spike count = poissonrnd(R*dt);
%Careful scaling property delta(a*t)=delta(t)/a;

%Case 2: diffusion approx.
%S(t) = Re + sqrt(Re)*eta;

%Reudction to current-based model
%Choose a holding potential Ei<v0<Ee and replace:
%we=-ae*(v0-Ee)>0; wi=-ai*(v0-Ei)<0;
%dV/dt = -(V-El)/tau + we*Se + wi*Si;

%Case 1: Full spiking neuron as before;
%Case 2: Diffusion approximation. Merge two noise source into one.
%dV/dt = -(V-El)/tau + u + sigma*eta; (Ito sense)
%where u = we*Re*+wi*Ri; var = sigma^2 = we^2*Re + wi^2*Ri;
%FPE is (Ito sense)
%dP/dt = D[(V-El-tau*u)/tau*P] +  1/2*var*D2[P]
%Or dP/dt = -dJ/dV (*)
% -J = (V-El-tau*u)/tau*P +  1/2*var*dP/dV (**)
%J(V) = probability 'flux'. No. of particles crossing V per unit time
%For stationary solution, we have:
%(**)=>  -dP/dV =  [(V-El-tau*u)/tau*P + J]/var*2
%(*)=> -dJ/dV = 0; except at boundary:
% -dJ/dV = r*delta for V-->Vth
% -dJ/dV = -r*delta for V-->Vreset
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

p=parameters;

%we = -ae*(p.v0-p.Ee);
%wi = -ai*(p.v0-p.Ei);
%W = [ae ai;ae ai]; %[wee wei;wie wii]

if mod(p.nfp,2)==0
    disp('Number of neuron changed to odd so convnfft_c works correctly. (Fix it?)')
    p.nfp=p.nfp+1;
end

x=linspace(-pi,pi,p.nfp)'; %column vector
dx=x(2)-x(1);


%Exponential weight function
KE=p.K;KI=p.K;
ke = KE*vonMises(x,p.de)*dx;
ki = KI*vonMises(x,p.di)*dx;


if isempty(ini_cond)
r=100*cos(x/2).^2;
r=[r;r/2];
else
    r=ini_cond;r=r(:); %make sure its column vector
end

TSPAN=linspace(0,T,1e2); %ms

[TSPAN,R] = ode23(@(t,r) myrhs(r,[we wi],ke,ki,p),TSPAN,r);

R=R'; %MATLAB is horrible; require column vector for input, output row vector
re=R(1:p.nfp,end);ri=R(p.nfp+1:end,end);

Re = convnfft_c(ke,re) + p.I0; %NB: ke = kie=kee; excitatory input poisson rate
Ri = convnfft_c(ki,ri);      %inhibitory input poisson rate

% WEE=W(1,1);WEI=W(1,2);WIE=W(2,1);WII=W(2,2);%[wee wei;wie wii]
% 
% uee = WEE*ue;uei = WEI*ui;  %currents
% uie = WIE*ue;uii = WII*ui;
% 
% ue = uee + uei; %mean current
% ui = uie + uii;
% vare = WEE*uee+WEI*uei; %var current
% vari = WIE*uie+WII*uii;
% 
% se=sqrt(vare);
% si=sqrt(vari);

fpe=fokkerplanck(Re,Ri,ae,ai,p);
fpi=fpe;
%fpe=fokkerplanck(ue,vare,101,p);
%fpi=fokkerplanck(ui,vari,101,p);

opt.re=re;opt.ri=ri; %final state
% opt.ue=ue;opt.ui=ui;%mean current
% opt.se=se;opt.si=si; %std current
% opt.fpe=fpe;opt.fpi=fpi; %fp transfer function output; check convergence with r
% opt.x=x;
% opt.t=TSPAN;
% opt.pe=ke; %connection probablity; sum(pe) should be equal to KE
% opt.pi=ki;
% opt.wee=WEE;opt.wii=WII;
% opt.wei=WEI;opt.wie=WIE;
% opt.iext=I0;
% opt.dx=2*pi/p.nfp;
% opt.w=WEE*ke-WEI*ki;
% opt.de=dE;
% opt.di=dI;


function phi=myrhs(r,W,ke,ki,p)
%INPUT: p = structure array of parameters
WEE=W(1,1);WEI=W(1,2);
%WIE=W(2,1);WII=W(2,2);%[wee wei;wie wii]

re = r(1:length(r)/2);
ri = r((length(r)/2+1):end);

Re = convnfft_c(ke,re) + p.I0; %NB: ke = kie=kee; total excitatory rate
Ri = convnfft_c(ki,ri);

fpe=fokkerplanck(Re,Ri,WEE,WEI,p);
fpi=fpe;
%fpi=fokkerplanck(Re,Ri,WIE,WII);
fp = [fpe;fpi];
phi = (-r+fp)/p.tau;