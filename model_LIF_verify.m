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

clear
close all

m=6; %number of different parameters
dataset=4;


T=1; %simulation time (s)
dt = 0.1; %ms, time step
numit = floor(T*1e3/dt);
switch dataset
    case 1
p=parameters([],'Vlb',-150,'nV',2001);
re=linspace(1,1.5,m); ri=0.5;%linspace(0.5,1,m);%(kHz = ms^-1)
ae=0.025;ai=0.3;

zmax=15;
    case 2
        p=parameters([],'model_type','conduct');
 re=linspace(10,15,m); ri=linspace(4,6,m);%(kHz=ms^-1)
 ae=0.004;ai=0.026;
 zmax=7;
    case 3
        p=parameters([],'model_type','conduct');
        %re=p.K*linspace(1.5e-3,2.5e-3,m)+p.I0; ri=p.K*2e-3;%linspace(4,6,m);%(kHz=ms^-1)
        re=linspace(10,15,m); ri=4;%linspace(4,6,m);
        %re=linspace(5,10,m); ri=1;%linspace(4,6,m);
        ae=0.004;ai=0.026;
        %ae=0.006;ai=0.026;
        %zmax=7;
    case 4
        p=parameters([],'model_type','conduct');
        m=10;
        re=3.85*ones(1,m); ri=2.31;%linspace(4,6,m);
        ae=0.0235;ai=0.2;
        
        T=10; numit = floor(T*1e3/dt);
        
        V=zeros(m,numit);Z=V;
        IV=V;IZ=V;
end





[Re,Ri]=ndgrid(re,ri); %consant total E/I synaptic Poisson rate



switch p.model_type
    case 'current'
        we = -ae.*(p.v0-p.Ee); wi= -ai.*(p.v0-p.Ei);
        u = we.*Re + wi.*Ri;
        var = we.*we.*Re + wi.*wi.*Ri;
        sigma = sqrt(var);
    case 'conduct'
        be = 1-exp(-ae); %ae in non-dimensional unit relative to C
        bi = 1-exp(-ai);

end


%>>Set initial condition
v=p.Vref+3*randn(size(Re)); %initial membrane potential
z=v; %Lengevin euqation using diffusion approximation
fire_old=false(size(v));
firez_old=false(size(v));

%NB apparently double sparse is the fastest
%<<



%V=zeros(N,numit);
Spk=false(length(re),length(ri),numit); %spikes
Spkz=Spk;

%V = zeros(numit,m); dont save all these data too many
%Z = V;

bin=linspace(p.Vlb,p.Vth,201);
dV = bin(2)-bin(1);
%bin = bin+ dV/2; %centre of bins; make sure Vref is a boundary
hV = zeros(length(bin),m);
hZ = hV;


flag = 1; %1=Ito or 0=Stratonovich %no difference for current-based
disp('Start simulation of neuron...')
tic
for i=1:numit
    
    %Ne = poissrnd(Re*dt);
    %Ni = poissrnd(Ri*dt);
    Ne = rand(size(Re))<(Re*dt);  %valid when Re*dt<<1
    Ni = rand(size(Ri))<(Ri*dt);
    
    
    
    switch p.model_type
        case 'current'
            dv = -(v-p.El)/p.tau*dt + we.*Ne + wi.*Ni;
            %<<< this is the correct integration! the units are correct!
            %Ito integral:
            dz = -(z-p.El)/p.tau*dt + u*dt + sigma.*randn(size(z))*sqrt(dt);
        case 'conduct'
            Iv = -be*(v-p.Ee).*Ne - bi*(v-p.Ei).*Ni;
            dv = -(v-p.El)/p.tau*dt +Iv;
            
             
            %<<< this is the correct integration! the units are correct!
            %>>>integrate using Milstein method
            dWe = sqrt(dt)*randn(size(z));
            dWi = sqrt(dt)*randn(size(z));
            f = -(z-p.El)/p.tau -be*(z-p.Ee).*Re -bi*(z-p.Ei).*Ri;
            ge= -be*(v-p.Ee).*sqrt(Re).*dWe + 0.5*be*be*(v-p.Ee).*Re.*(dWe.*dWe-flag*dt);
            gi= -bi*(v-p.Ei).*sqrt(Ri).*dWi + 0.5*bi*bi*(v-p.Ei).*Ri.*(dWi.*dWi-flag*dt);
            dz = f*dt + ge + gi;
            
            Iz = -be*(z-p.Ee).*(Re+sqrt(Re).*dWe/sqrt(dt)) - bi*(z-p.Ei).*(Ri+sqrt(Ri).*dWi/sqrt(dt));
            %Iz = sqrt(be*be*(z-p.Ee).*(z-p.Ee).*(Re+sqrt(Re).*dWe/sqrt(dt)) - bi*(z-p.Ei).*(Ri+sqrt(Ri).*dWi/sqrt(dt));
            %<<<
    end
    v=v+dv;
    z=z+dz;
    %vi=vi+dvi;
    
    fire = v>p.Vth;
    firez = z>p.Vth;
        
    %record spikes
    Spk(:,:,i)= fire;
    Spkz(:,:,i)= firez;
    
    %>>> reset previously firing neuron to reset potential
    v(fire_old)=p.Vref;
    z(firez_old)=p.Vref;
    fire_old = fire;
    firez_old = firez;
    %<<<
        
    %record voltage after reset previous spikes
    %switch dataset
    %    case 4
    V(:,i)=v;IV(:,i)=Iv;
    Z(:,i)=z;IZ(:,i)=Iz;
    %end
    
    for j=1:m
    indx = v(j)<bin+dV/2 & v(j)>bin-dV/2 ;
    hV(indx,j) =hV(indx,j)+1; 
    indx = z(j)<bin+dV/2 & z(j)>bin-dV/2;
    hZ(indx,j) = hZ(indx,j)+1;    
    end
    
    
end
toc
return

%Firing rate
R=sum(Spk,3)/T;
Rz=sum(Spkz,3)/T;
hV = hV/(numit*dV);
hZ = hZ/(numit*dV);

%%
% %empirical voltage distribution
% bin=linspace(p.Vlb,p.Vth,101);
% dV = bin(2)-bin(1);
% bin = bin + dV/2; %centre of bins; make sure Vref is a boundary
% hV = zeros(length(bin),m);
% hZ = hV;
% for i=1:m
%     hV(:,i)=hist(V(:,i),bin)/(numit*dV);
%     hZ(:,i)=hist(Z(:,i),bin)/(numit*dV);
% end

%calculate empircial current J(Vth)
%J(Vth) = -D*dP - F*P; This should be equal to r0


%u=u*tau;  %make spiking simulation consistent with diffusion approx
%v = v*tau*tau;
%%
tic
disp('Start solving FPE...')
[r,prb]=fokkerplanck(Re(:),Ri(:),ae,ai,p);%parameters([],'model_type','conduct','nV',10001));
%[r0_euler,prb]=fokkerplanck0_euler(Re(:),Ri(:),ae,ai,parameters([],'model_type','conduct','nV',10001));
toc
r=reshape(r,size(Re));
%rr=tmp(u,var,p.n);

%gamma = (Re*ae*ae+Ri*ai*ai)./(tau+Re*ae+Ri*ai)/2; %diffusion approx is good for gamma<<1
return
%% PLOT 1D curve (dataset 3)
close all
%plot(re,r*1e3,re,R,re,Rz)
plot(re,Rz,'.',re,R,'.')
hold on
errorbar(re,r*1e3,sqrt(r*1e3/T)); %Theoretical firing rate & error
hold off
%% PLOT distribution
ii=3;
plot(bin,hZ(:,ii),'.',bin,hV(:,ii),'.',prb.V,prb.p(ii,:),[-60 -60],[0 max(prb.p(ii,:))*1.3],'--')
%plot(bin,hZ(:,ii),'.',prb.V,prb.p(ii,:),[-60 -60],[0 max(prb.p(ii,:))*1.3],'--')
%% PLOT V(t)
close;plot(Z(1:1e5,end));
%% PLOT 2D mesh

close all
subplot(1,3,1)
mesh(ri,re,R);title('LIF neuron')
ylabel('R_e (kHz)');xlabel('R_i (kHz)');zlabel('r_{out} (Hz)')
zlim([0 zmax])
subplot(1,3,3)
mesh(ri,re,r*1e3);title('PFE')
zlim([0 zmax])
subplot(1,3,2)
mesh(ri,re,Rz);title('Diffusion approx.')
zlim([0 zmax])
set(gcf,'position',[100 100 800 250],'color','w')

%
%subplot(2,2,4)
%mesh(rr);title('PFE v1')
%plot(sum(Spk))

% %calculate CV of ISI
% isi=[];
% for i=1:N
%     t=find(Spk(i,:));
%     if length(t)>1
%         isi=[isi diff(t)];
%     end
% end
%
% CV= std(isi)/mean(isi);
% disp(['CV of ISI is ' num2str(CV)])
%
% %mean firing rate of all neuron
% rr=sum(Spk(:))/T/N;
% disp(['Mean firing rate is ' num2str(rr) ' Hz'])
%
% %mean spatial firing rate profile
% nbin = 10; %number of bins
% binsize = floor(N/nbin);
% R = zeros(1,nbin);
% for i=1:nbin
% ss = Spk(i:(i+binsize-1),:);
% R(i)=sum(ss(:))/T/binsize;
% end
%
% %output spike data as points
% [xx,tt]=find(Spk);
% tt=tt*dt;
% x=x(xx);
% Spk=[tt(:) x(:)];
%
%



%subplot(2,2,1)
%imagesc(V,[Vref Vth])
%colorbar
%subplot(2,2,2)
%plot(tt,xx,'.');axis([0 dt*numit 0 2*pi])
%imagesc(Spk,[0 1])
%subplot(2,2,3)
%plot(x,sum(Spk,2)/T,'.')
%subplot(2,2,4)
%plot(V(2500,:))

