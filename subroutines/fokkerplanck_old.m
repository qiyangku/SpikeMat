function [r,prb]=fokkerplanck0(Re,Ri,ae,ai,p)
%INPUT: mean current, variance of current, discretization
%OUTPUT: mean firint rate
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




switch p.model_type
    case 'current'
        we = -ae.*(p.v0-p.Ee); wi= -ai.*(p.v0-p.Ei);
        u0 = we.*Re + wi.*Ri;
        var = we.*we.*Re + wi.*wi.*Ri;
    case 'conduct'    
        be=1-exp(-ae);bi=1-exp(-ai);
end

%n=101; %discretization of membrane potential
dV = (p.Vth-p.Vlb)/(p.nV-1);

kref = (p.Vref-p.Vlb)/dV+1; %index for refractory potential %must be integer

if mod(kref,1)~=0
    disp('kref is not integer! Rounding to nearest integer...')
    kref=round(kref);
end


V=linspace(p.Vlb,p.Vth,p.nV);
P=zeros(size(Re,1),size(Re,2),p.nV);
%J=zeros(1,n);


%initial condition
pp=0;
j=1;
v=V(end);
P(:,:,end)=pp;
%J(end)=j;
for k=(p.nV-1):-1:1
    
    switch p.model_type
        case 'current'            
            u= (v-p.El)/p.tau - u0;
        case 'conduct'
            %var current (Eq.64)>>
            var = Re.*(be.*(v-p.Ee)).^2 + Ri.*(bi.*(v-p.Ei)).^2;
            %mean current (stratonovich) Eq.68
            u = (v-p.El)/p.tau + Re.*(be+0.5*be.*be).*(v-p.Ee) + Ri.*(bi+0.5*bi.*bi).*(v-p.Ei);          
    end
    G=dV*u./var*2;
    A = exp(G);
    B=(A-1)./G;
    B(G==0)=1;
    B=dV*B./var*2;
    
    j = j - (k==(kref+1));
    pp = pp.*A + j.*B;
    v = V(k);
    %size(p)
    P(:,:,k)=pp;
    %J(i)=j;
end

%firing rate
r = 1./(sum(P,3)*dV);
prb.V=V;
try
    P=r.*P;  %reuires MATLAB2016+    
    prb.p=P;    
    prb.dp=r.*dp;
catch
end
r = r*1e3; %convert unit from ms-1 to Hz

if ~p.stats
    return
end

disp('Computing additional statistics...')
%>> additional statistics
%>> Voltage stats
uV = sum(reshape(V,[1 1 n]).*P*dV,3); %expected value of E[V]
vVe = sum( reshape((V-p.Ee).^2,[1 1 n]).*P*dV,3);
vVi = sum( reshape((V-p.Ei).^2,[1 1 n]).*P*dV,3);
vVei = sum( reshape((V-p.Ei).*(V-p.Ee),[1 1 n]).*P*dV,3); %covariance term

switch p.model_type
    case 'current'
        ue = we.*Re;     ui = wi.*Ri;
        ve = we.*we.*Re; vi= wi.*wi.*Ri;
        vtot=var;
    case 'conduct'
        %THESE FORMULAE NEEDS CORRECTION!
        ue = -(be+0.5*be.*be)*Re.*(uV-p.Ee); %mean excitatory current
        ui = -(be+0.5*be.*be)*Ri.*(uV-p.Ei); %mean inhibitory current
        %ul = -(uV-p.El); %mean leak current
        utot = ue+ui; %total synaptic current
        %vl = sum( reshape((V-p.El).^2,[1 1 n]).*P*dV,3);%leak var
        ve = be*be*Re.*(1+Re).*vVe - ue.*ue; %double check this formula is correct
        vi = bi*bi*Ri.*(1+Ri).*vVi - ui.*ue;
        vtot = be*be*Re.*(1+Re).*vVe + bi*bi*Ri.*(1+Ri).*vVi; %var part
        vtot = vtot + 2*be*bi*vVei.*Re.*Ri; %plus the covariance part
        vtot = vtot -utot.*utot; %total variance of current
        %NB these variances are not additive as V in each term are not independent
end

prb.ue=ue;prb.ve=ve;
prb.ui=ui;prb.vi=vi;
%prb.ul=ul;prb.vl=vl;
prb.vtot=vtot;
prb.uV=uV;
disp('Done!')