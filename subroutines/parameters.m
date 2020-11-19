function p=parameters(preset,varargin)
%INPUT: p = old parameter values; varargin = new parameter name/value pairs
%OUTPUT: p = new parameter values

if isempty(preset)
    %disp('Loading default parameters...')
    preset='qi1d';
end
%preset='qi2d';
%preset='yifan';




switch preset
    case 'qi1d'
        disp('Loading parameter preset: qi1d')
        p.dim=1; %spatial dimensions of network
        %network parameters
        p.K=500;  %number of synapses
        p.N=5000; %number of neurons (does not matter too much)
        p.I0=0.5; %kHz
        
        p.de=1;       %Excitatory synaptic range
        p.di=2*p.de;       %Inhibitory synaptic range
        p.coupling_type='vonMises';
        
        p.separateEIpop = false;
    case 'qi2d'
        disp('Loading parameter preset: qi2d')
        p.coupling_type='expo';
        p.dim=2; %spatial dimensions of network
        p.aee=0.02;%0.01625; %Excitatory synaptic weight; ge/C
        p.aie=0.02;%0.175;   %Inhibitory synaptic weight; gi/C
        p.aei=0.2;
        p.aii=0.2;
        
        p.dee=1;%0.8107;
        p.dei=2;
        p.die=1;
        p.dii=2;
        
        p.KEE = 500;
        p.KEI = 500;
        p.KIE = 500;
        p.KII = 500;
        %external input parameters
        p.Ie = 0.5; %rate to E pop (kHz)
        p.Ii = 0.5; %rate to I pop
        p.jee = p.aee; %strength to E pop
        p.jie = p.jee; %strength to I pop
        p.separateEIpop = true;%false;
    case 'yifan'
        disp('Loading parameter preset: yifan')
        p.coupling_type='expo';
        
        p.dim=2; %spatial dimensions of network
        p.aee=0.052;%0.01625; %Excitatory synaptic weight; ge/C
        p.aie=0.065;%0.175;   %Inhibitory synaptic weight; gi/C
        p.aei=0.105;
        p.aii=0.195;
        
        p.dee=8/31*pi;%0.8107;
        p.dei=20/31*pi;
        p.die=10/31*pi;
        p.dii=p.dei;
        p.KEE = 640;
        p.KEI = 200;
        p.KIE = 800;
        p.KII = 400;
        %external input parameters
        p.Ie = 0.85; %rate to E pop (kHz)
        p.Ii = 1.0; %rate to I pop
        p.jee = 0.5*p.aee; %strength to E pop
        p.jie = p.jee; %strength to I pop
        p.separateEIpop = true;%false;
        
        p.see = 1.9/4; %CV of strength (std/mean)
        p.sei = 0.25;
        p.sie = 0;
        p.sii = 0;
        
end




%neural parameters
%p.model_type='current';
p.model_type='conduct';
p.Vlb = -100; %lower bound (mV) for solving FPE
p.Vth = -50; %firing threshold (mV)
p.Vref = -60; %refractory potential (mV)
p.Ee =   0; %excitatory reversal potential
p.Ei = -80; %inhibitory reversal potential
p.El = -70; %leak reversal potential
p.tau=15;   %membrane time constant (ms) = C/gL
p.v0=-65; %holding potential for current-based model



%find fixed point parameters
%if ~isfield(p,'dim')
%    disp('checkpoint')

%end

switch p.dim
    case 1
        p.nx = 1001; %spatial discretisation ('number of neurons')
    case 2
        p.nx = 51; %spatial discretisation ('number of neurons')
end

p.ninit = 6; %number of initial conditions for root finding

%FPE solver parameters
p.nV=1001; %discretization in V
%(current-based model is not sensitive to p.n as numerical scheme exact)
p.stats = false; %compute additional stats from FPE?

%display parameters
p.print_interval=10; %interval of printing progress (unit percentage)

p.preset=preset;

if ~isempty(varargin)
    disp('Updating values...')
    for i=1:2:length(varargin)
        p.(varargin{i})=varargin{i+1};
    end
    disp('Done. Calling this function again will reset to default parameters!')
else
    disp('Done!')
end




