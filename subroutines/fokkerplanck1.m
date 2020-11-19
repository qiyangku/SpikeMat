function [r1,prb] = fokkerplanck1(w,x,flag,ae,ai,p)
%[r1,prb] = fokkerplanck1(w,x,p)
%INPUTS: x = stationary probabilty density/its derivative/current/parameters
%       flag = 1 indicates perturb. in Re; 0 indicates perturb. in Ri
%       w = angular frequency
%       p = other parameters
%OUTPUT: r1 = susceptibility function of firing rate

%P1a=zeros(p.nV,1);
%P1b=P1a;
w=1i*w;  %

Re=x.Re;Ri=x.Ri;Rext=x.Rext;  %kHz
Re1=flag; Ri1=1-flag;
%Re1=0.01*Re1*Re; Ri1=Ri1*Ri;%amplitude of perturbation



switch p.model_type
    case 'current'
        we = -ae.*(p.v0-p.Ee); wi= -ai.*(p.v0-p.Ei);
        u0 = we.*Re + wi.*Ri;
        F1= -we*Re1 - wi*Ri1;
        D0 = 0.5*(we.*we.*Re + wi.*wi.*Ri);
        D1 = 0.5*(we.*we.*Re1 + wi.*wi.*Ri1);
    case 'conduct'    
        be=1-exp(-ae(1));bi=1-exp(-ai(1));bext=1-exp(-x.aext);
        be2 = be*be*(1+ae(2)*ae(2));
        bi2 = bi*bi*(1+ai(2)*ai(2));
end

%n=101; %discretization of membrane potential
dV = (p.Vth-p.Vlb)/(p.nV-1);
kref = (p.Vref-p.Vlb)/dV+1; %index for refractory potential %must be integer



%initial conditions; two pairs of variables
p1a=0;
p1b=0;
j1a=1;
j1b=0;


P1a=zeros(length(w),p.nV);
P1b=P1a;
J1a=P1a;J1b=P1a;

for k=p.nV:-1:1
    
    v = x.V(k);
    dp0 = x.dp(k); %derivative of stationary distribution
    p0  = x.p(k);  %stationary distribution
    
    
    %P(:,k)=pp;
    
    %>>>Calculate drift & diffusion terms for current v
    switch p.model_type
        case 'current'            
            F0= (v-p.El)/p.tau - u0;
        case 'conduct'            
            
            we = (be+be2).*(v-p.Ee); wi =(bi+bi2).*(v-p.Ei);
            wext = (bext+bext.*bext).*(v-p.Ee);
            ue = be2.*(v-p.Ee).^2;     ui = bi2.*(v-p.Ei).^2;
            uext = (bext.*(v-p.Ee)).^2;
            
            D0 = 0.5*( Re.*ue + Ri.*ui + Rext.*uext);
            F0 = (v-p.El)/p.tau + Re.*we + Ri.*wi +Rext.*wext;
            F1= Re1.*we + Ri1.*wi;
            D1 = 0.5*(Re1.*ue + Ri1.*ui);
    end
    %<<<
    
    %store value of current step    
    
    
    %>>integration (p1a,j1a)
    G=dV*F0./D0;
    A = exp(G);
    B=(A-1)./G;
    B(G==0)=1;
    B=dV*B./D0;    
    
    j1a_new = j1a + dV.*(w.*p1a) - (k==(kref+1));  %Euler scheme
    %j1 and p1 are coupled (unlike j0 and p0). Exponetial integration
    %scheme needs to be updated accordingly? (Richardson2007 used Euler)
    p1a_new = p1a.*A + j1a.*B;
        
    j1a = j1a_new;
    p1a = p1a_new;
    
    %<<
    
    %>>integration (p1b,j1b)
    G=dV*F0./D0;
    A = exp(G);
    B=(A-1)./G;
    B(G==0)=1;
    B=dV*B./D0;
    
    j1b_new = j1b + dV.*(w.*p1b) ;  %Euler scheme
    p1b_new = p1b.*A + (j1b + D1.*dp0+F1.*p0).*B;
           
    j1b = j1b_new;
    p1b = p1b_new;
    
    
    
    
    P1a(:,k)=p1a;
    P1b(:,k)=p1b;
    J1a(:,k)=j1a;
    J1b(:,k)=j1b;
    
end

r1 = -j1b./j1a; %in unit of kHz???
%change in output rate = r1 * change in input rate
prb.p1=r1.*P1a+P1b;
prb.p1a=P1a;
prb.p1b=P1b;
prb.j1a=J1a;
prb.j1b=J1b;


%return

%asymptote for high frequency limite w=0+i*frequency
if length(ae)>1 
    warning('Only support calculating asymptote for scalar ae & ai.')
    return
end
    
we = (be+be.*be).*(x.V-p.Ee); wi =(bi+bi.*bi).*(x.V-p.Ei);
%wext = (bext+bext.*bext).*(x.V-p.Ee);
ue = (be.*(x.V-p.Ee)).^2;     ui = (bi.*(x.V-p.Ei)).^2;
uext = (bext.*(x.V-p.Ee)).^2;

D0 = 0.5*( Re.*ue + Ri.*ui + Rext.*uext);
%F0 = (v-p.El)/p.tau + Re.*we + Ri.*wi +Rext.*wext;
F1= Re1.*we + Ri1.*wi;
D1 = 0.5*(Re1.*ue + Ri1.*ui);

Fa=D1.*x.dp+F1.*x.p;
%abs(Fa(end))
dFa = (0.5*Fa(end-2) -2*Fa(end-1) + 1.5*Fa(end))/dV; %backward finite difference


prb.limit = exp(1i*4*pi)*sqrt(D0(end)./(-1i*w))*dFa - Fa(end);


%r1 = r1*1e3; %convert to Hz

