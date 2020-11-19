% For verifying Asem's theory of super-diffusion approximation
%Simulate Leaky-Integrate-and-Fire neuron.
%INPUT: ae,ai = ge/C >0,gi/C >0, synaptic weights; %T = simulation time (s)
%p = parameters (generated using the function parameters);
%OUTPUT: CV = CV of inter-spike-interval; Spk = [spike time x spike position];
%R=mean firing rate spatial profile; meanr = total mean firing rate (Hz)
clear
p=parameters('yifan'); %parameters

T = 10; %simulation time (s)

Ttrans=1; %duration of transient external input (s)
dt = 0.1; %ms, time step

ae = 0.01625;
ai = 0.175;


%model_type='current';


%switch p.model_type
%    case 'current'
        wee= -ae*(p.v0-p.Ee); %current-based model weight
        wei= -ai*(p.v0-p.Ei);
        wie=wee;
        wii=wei;
%     case 'conduct'
%         be=1-exp(-ae);
%         bi=1-exp(-ai);
%         jee=be;jie=be; %conductance-based model weight
%         jei=bi;jii=bi;
% end

p.I0=0.8; %kHz %p.I0*1e-3; % convert unit to ms^-1
p.N = 1e3; %number of neurons

x=linspace(-pi,pi,p.N+1);
x=x(1:p.N)';

%K=900;
%tic
disp('Generating coupling matrix...')
% Cee=coupling(p.K,p.de,p.N); %binary connectivity matrix
% Cie=coupling(p.K,p.de,p.N);
% Cei=coupling(p.K,p.di,p.N);
% Cii=coupling(p.K,p.di,p.N);

pp = 0.01;
Cee = (rand(p.N)<pp).*(wee + sqrt(wee)*randn(p.N));
Cei = (rand(p.N)<pp).*(wei + sqrt(wei)*randn(p.N));
Cie = (rand(p.N)<pp).*(wie + sqrt(wie)*randn(p.N));
Cii = (rand(p.N)<pp).*(wii + sqrt(wii)*randn(p.N));

Cee = sparse(Cee); Cei = sparse(Cei);
Cie = sparse(Cie); Cii = sparse(Cii);

%disp('Done.')
%toc
disp('Starting simulation...')


ve=p.El+(p.Vth-p.El)*rand(p.N,1); %initial membrane potential
vi=p.El+(p.Vth-p.El)*rand(p.N,1);
firee = sparse(double(ve>p.Vth));    firei = sparse(double(vi>p.Vth));

numit = floor(T*1e3/dt);

%V=zeros(N,numit);
Spk=false(p.N,numit); %spikes

%external current rate (ms^-1)
%Iext = 1 + 0.2*exp(-x.*x/2/1/1);
Iext=p.I0;

mark = p.print_interval; %how often print message

for i=1:numit
    
    if i>floor(Ttrans*1e3/dt)
        Iext=p.I0*ones(size(ve));
    end
    
    Nee = sparse(poissrnd(Iext*dt)); %external input spike to E neurons
    Nie = sparse(poissrnd(Iext*dt)); %external input spike to I neurons
    
    Ree = Cee*firee  + wee*Nee; %total excitatory Poisson rate to E neurons
    Rei = Cei*firei;        %total inhibitory Poisson rate to I neurons
    Rie = Cie*firee  + wie*Nie;
    Rii = Cii*firei;
    
%     switch p.model_type
%         case 'current'
%     ue =  wee*Ree + wei*Rei;
%     ui =  wie*Rie + wii*Rii;  
%         case 'conduct'
%             
%     ue =  -jee*(ve-p.Ee).*Ree - jei*(ve-p.Ei).*Rei;
%     ui =  -jie*(vi-p.Ee).*Rie - jii*(vi-p.Ei).*Rii;    
%     end
    
    dve = -(ve-p.El)/p.tau*dt + Ree + Rei;
    dvi = -(vi-p.El)/p.tau*dt + Rie + Rii;
    ve=ve+dve;
    vi=vi+dvi;
    
    firee = ve>p.Vth;    firei = vi>p.Vth;
    ve(firee)=p.Vref;    vi(firei)=p.Vref;
    Spk(:,i)= firee;      
    
    
        if i>floor(2*Ttrans*1e3/dt) && sum(firee)> p.N*dt %If mean istantaneous rate exceeds 1 kHz
        disp('Activity blowed up!')
        Spk=inf;R=inf;CV=NaN;meanr=inf;
        return
        end
    
    
    
    firee=double(sparse(firee));
    firei=double(sparse(firei));
    
    
    progress=(i-1)/numit*100; %percentage progress
    if progress > mark
        disp(['Progress: ' num2str(floor(progress)) '%...'])
        mark = mark + 5;
        %toc
    end

end
disp('Simulation completed!')
%toc
return


%calculate CV of ISI
isi=[];
for i=1:p.N
    t=find(Spk(i,:));
    if length(t)>1
        isi=[isi diff(t)];
    end
end

CV= std(isi)/mean(isi);
disp(['CV of ISI is ' num2str(CV)])

%mean firing rate of all neuron
meanr=sum(Spk(:))/T/p.N;
disp(['Mean firing rate is ' num2str(meanr) ' Hz'])

%mean spatial firing rate profile
nbin = 10; %number of bins
binsize = floor(p.N/nbin);
R = zeros(1,nbin);
for i=1:nbin
ss = Spk(i:(i+binsize-1),:);
R(i)=sum(ss(:))/T/binsize;
end

%output spike data as points
[xx,tt]=find(Spk);
tt=tt*dt*1e-3;
x=x(xx);
Spk=[tt(:) x(:)];

return
plot(Spk(:,1),Spk(:,2),'.',Ttrans*[1 1],[-pi pi],'--')
xlabel('t (ms)')
ylabel('x')
ylim([-pi pi])

