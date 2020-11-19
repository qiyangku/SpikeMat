function [R,opt]=fixedpt_ode23(ae,ai,T,ini_cond,p)
%Same as fixedpt_ode23 except E/I populations are the same
%INPUT: ae,ai = excitatory/inhibitory synaptic strength
%T: simulation time
%%%%%%%%%% CHANGE LOG %%%%%%%%%%%%
%ver0.5:3/3/2017; FPE verified. Consistent with Langevin & diff approx
%Changed expression comparable with current-based
%combined E/I population.
%ver0.3;23/02/2017; separated E/I population. 
%Connection probability only depend on pre-synaptic neuron (ke=kie=kee)
%Weight also depend only on pre-synaptic neuron
%But only E population receives I0 (the only asymmetry)
%ver0.4;26/02/2017; add variance of I0. mean(I0) = strength*rate;
%var(I0)=strength^2*rate; 
%assume every E neuron receives idd Poisson input spikes
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%Documentation
%Start with cond. based model
%C*dV/dt = -gL*(V-El) - ge*(V-Ee)*Se - gi*(V-Ei)*Si (1)
%where S(t) is instantaneous firing rate (spike count per unit time)
%Divide C through Eq.1
%dV/dt = -(V-El)/tau - ae*(V-Ee)*Se - ai*(V-Ei)*Si (2)
%
%Case 1: Full spiking neuron
%S(t) = Poisson spike train. spike count = poissonrnd(R*dt);
%Careful scaling property delta(a*t)=delta(t)/a;
%
%Case 2: diffusion approx.
%S(t) = Re + sqrt(Re)*eta;
%
%Reudction to current-based model
%Choose a holding potential Ei<v0<Ee and replace:
%we=-ae*(v0-Ee)>0; wi=-ai*(v0-Ei)<0;
%dV/dt = -(V-El)/tau + we*Se + wi*Si;
%
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

%p=parameters;

%W = [we wi;we wi]; %[wee wei;wie wii]

if mod(p.nx,2)==0
    disp('Number of neuron changed to odd so convnfft_c works correctly. (Fix it?)')
    p.nx=p.nx+1;
end

x=linspace(-pi,pi,p.nx)'; %column vector
dx=x(2)-x(1);


%initial condition
if isempty(ini_cond)
%r=50e-3*cos(x/2).^2; %kHz
r=50e-3*exp(0.5*cos(x));
%r=[r;r/2];
else
    r=ini_cond;r=r(:); %make sure its column vector
end


%weight function

switch p.dim
    case 1
        KE=p.K;KI=p.K;
        ke = KE*coupling_fun(x,p.de,p)*dx;
        ki = KI*coupling_fun(x,p.di,p)*dx;
    case 2
        
        if p.separateEIpop
            %nn=length(r)/2;
            %r0e = r(1:nn);r0i=r(nn+1:end);
            r0e = 0.1*r;
            r0i = r;
            r0e = sqrt(r0e*r0e'); r0e=r0e(:);%outer product so r(x,y)=sqrt(r(x)r(y))
            r0i = sqrt(r0i*r0i'); r0i=r0i(:);%outer product so r(x,y)=sqrt(r(x)r(y))
            k.ee = coupling_fun(x,p.dee,p)*dx*dx; %ke = ke*ke';
            k.ei = coupling_fun(x,p.dei,p)*dx*dx; %ki = ki*ki';
            k.ie = coupling_fun(x,p.die,p)*dx*dx; %ke = ke*ke';
            k.ii = coupling_fun(x,p.dii,p)*dx*dx; %ki = ki*ki';
            a.ee = [ae(1) p.see]; a.ie=[ae(2) p.sie];
            a.ei = [ai(1) p.sei]; a.ii=[ai(2) p.sii];
        else
            r = sqrt(r*r'); r=r(:);%outer product so r(x,y)=sqrt(r(x)r(y))
        end
        
        %kee = vonMises(x,p.dee)*dx;
        %kei = vonMises(x,p.dei)*dx;
        %kie = vonMises(x,p.die)*dx;
        %kii = vonMises(x,p.dii)*dx;
        
        %kee = p.KEE*(kee*kee');
        %kei = p.KEI*(kei*kei');
        %kie = p.KIE*(kie*kie');
        %kii = p.KII*(kii*kii');
        
end




TSPAN=linspace(0,T,1e2); %ms

options = odeset('Events',@myEventsFcn);

%tic
disp('Starting relaxation process...')
switch p.dim
    case 1
        [TSPAN,R] = ode23(@(t,r) myrhs(r,[ae ai],ke,ki,p),TSPAN,r,options);
    case 2
%         if length(ae)~=2
%             disp('WARNING: not enough parameters!')
%         end
        if p.separateEIpop
            %ae = [aee aei]; ai=[aie aii]
            [TSPAN,R] = ode23(@(t,r) myrhs2d_ei(r,a,k,p),TSPAN,[r0e;r0i],options);
        else
            [TSPAN,R] = ode23(@(t,r) myrhs2d(r,[ae ai],ke,ki,p),TSPAN,r,options);
        end
        
end
%toc

disp('Done! Outputing results...')

R=R'; %MATLAB is horrible; require column vector for input, output row vector
r=R(:,end);%ri=R(N+1:end,end);


switch p.dim
    case 1
        Re = convnfft_c(ke,r) + p.I0; %NB: ke = kie=kee; excitatory input poisson rate
        Ri = convnfft_c(ki,r);      %inhibitory input poisson rate
        [fp,prb]=fokkerplanck(Re,Ri,ae,ai,p);
    case 2
        if p.separateEIpop
            n=length(r)/2;re = r(1:n);ri = r(n+1:end);
            re = reshape(re,[p.nx p.nx]); %convert to 2d array
            ri = reshape(ri,[p.nx p.nx]);
            
            %Re = ; %spatial coupling depend only on the type of presynaptic neuron
            %Ri = convnfft_c(ki,ri);
            
            Ree = p.KEE*convnfft_c(k.ee,re);
            Rei = p.KEI*convnfft_c(k.ei,ri);
            Rie = p.KIE*convnfft_c(k.ie,re); %NB: ke = kie=kee; total excitatory rate
            Rii = p.KII*convnfft_c(k.ii,ri);
            Ree=Ree(:);Rei=Rei(:); %r=r(:);%convert back to 1d array
            Rie=Rie(:);Rii=Rii(:);
            %W = [aee aei aie aii]
            [fpe,prb]=fokkerplanck(Ree,Rei,p.Ie,a.ee,a.ei,p.jee,p);
            [fpi,prb]=fokkerplanck(Rie,Rii,p.Ii,a.ie,a.ii,p.jie,p);
            fp = [fpe(:);fpi(:)];
        else
            r = reshape(r,[p.nx p.nx]); %convert to 2d array
            Re = p.K*convnfft_c(ke,r) + p.I0; %NB: ke = kie=kee; total excitatory rate
            Ri = p.K*convnfft_c(ki,r);
            Re=Re(:);Ri=Ri(:); r=r(:);%convert back to 1d array
            [fp,prb]=fokkerplanck(Re,Ri,p.I0,ae,ai,p.jee,p);
        end
        
end



opt.r=r;
opt.R.ee=Ree;opt.R.ie=Rie;
opt.R.ei=Rei;opt.R.ii=Rii;
opt.fp=fp;
opt.x=x;
opt.t=TSPAN;
opt.k=k; %connection probablity; sum(pe) should be equal to KE
opt.dx=2*pi/p.nx;
opt.prb=prb;


function [value,isterminal,direction] = myEventsFcn(t,y)
%Terminate solver if:
%maximum firing rate is too high
value = (max(y(:))-1);
isterminal = 1;  % Halt integration 
direction = 0;

function phi=myrhs(r,W,ke,ki,p)
Re = convnfft_c(ke,r) + p.I0; %NB: ke = kie=kee; total excitatory rate
Ri = convnfft_c(ki,r);
[fp,prb]=fokkerplanck(Re,Ri,W(1),W(2),p);
phi = (-r+fp)/p.tau;

function phi=myrhs2d(r,W,ke,ki,p)
%two-dim version of myrhs
r = reshape(r,[p.nx p.nx]); %convert to 2d array
Re = p.K*convnfft_c(ke,r) + p.I0; %NB: ke = kie=kee; total excitatory rate
Ri = p.K*convnfft_c(ki,r);
Re=Re(:);Ri=Ri(:); r=r(:);%convert back to 1d array

[fp,prb]=fokkerplanck(Re,Ri,W(1),W(2),p);
phi = (-r+fp)/p.tau;

function phi=myrhs2d_ei(r,a,k,p)
%W = [aee aie aei aii]
%retrive E/I firing rates
n=length(r)/2;re = r(1:n);ri = r(n+1:end);
re = reshape(re,[p.nx p.nx]); %convert to 2d array
ri = reshape(ri,[p.nx p.nx]);

Ree = p.KEE*convnfft_c(k.ee,re);
Rei = p.KEI*convnfft_c(k.ei,ri);
Rie = p.KIE*convnfft_c(k.ie,re); %NB: ke = kie=kee; total excitatory rate
Rii = p.KII*convnfft_c(k.ii,ri);

Ree=Ree(:);Rei=Rei(:); %r=r(:);%convert back to 1d array
Rie=Rie(:);Rii=Rii(:);

[fpe,prb]=fokkerplanck(Ree,Rei,p.Ie,a.ee,a.ei,p.jee,p);
[fpi,prb]=fokkerplanck(Rie,Rii,p.Ii,a.ie,a.ii,p.jie,p);
fp = [fpe(:);fpi(:)];

phi = (-r+fp)/p.tau;



