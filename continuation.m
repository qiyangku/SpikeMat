function dat=continuation(varargin)
%Numeric continuation

p=parameters('yifan');
options = optimoptions('lsqnonlin','FunctionTolerance',1e-8,'OptimalityTolerance',1e-4,'Display', 'off');
%options = optimoptions('fsolve','FunctionTolerance',1e-8,'OptimalityTolerance',1e-8,'Display', 'off');

if isempty(varargin)
    disp('Calculating h for g=1')
    dat.h= lsqnonlin( @(x) objfun(1,x,p) , 0.5 , 0 , 5,options);
    return
end

switch varargin{1}
    
    case 3   %Search critical points in the (wei,Rext) plane
        n=30;
        p.aee = 0.4775*p.aee;
        
        h=linspace(0.8,1.2,2*n+1);
        g=zeros(size(h));
        resnorm=g;
        
        %Initialisation; calculate crit point for h(n+1)
        disp('Initialising...')
        
        tic
        [g(n+1),resnorm(n+1)]= lsqnonlin( @(x) objfun3(x,h(n+1),p) , 1 , 0 , 5,options);
        %[g(n+1),resnorm(n+1)]= fsolve( @(x) objfun(x,h(n+1),p) , 1 ,options);
        toc
        disp('Start numerical continuation...')
        tic
        for i=1:n
            
            %calculate crit point for n-i
            [g(n+1-i),resnorm(n+1-i)]= lsqnonlin( @(x) objfun3(x,h(n+1-i),p) , g(n+2-i) , 0 , 5,options);
            %[g(n+1-i),resnorm(n+1+i)]= fsolve( @(x) objfun(x,h(n+1-i),p) , g(n+2-i) ,options);
            %calculate crit point for n+i
            [g(n+1+i),resnorm(n+1+i)]= lsqnonlin( @(x) objfun3(x,h(n+1+i),p) , g(n+i) , 0 , 5,options);
            %[g(n+1+i),resnorm(n+1+i)]= fsolve( @(x) objfun(x,h(n+1+i),p) , g(n+i) ,options);
            clc
            toc
            progress=i/n
        end
        
        dat.g=g;
        dat.h=h;
        dat.resnorm=resnorm;
        disp('Numerical continuation complete!')
        
    case 1   %Search critical points in the (wei,wee) plane
        n=30;
        h=linspace(0.3,0.6,2*n+1);
        g=zeros(size(h));
        resnorm=g;
        
        %Initialisation; calculate crit point for h(n+1)
        disp('Initialising...')
        
        tic
        [g(n+1),resnorm(n+1)]= lsqnonlin( @(x) objfun(x,h(n+1),p) , 1 , 0 , 5,options);
        %[g(n+1),resnorm(n+1)]= fsolve( @(x) objfun(x,h(n+1),p) , 1 ,options);
        toc
        disp('Start numerical continuation...')
        tic
        for i=1:n
            
            %calculate crit point for n-i
            [g(n+1-i),resnorm(n+1-i)]= lsqnonlin( @(x) objfun(x,h(n+1-i),p) , g(n+2-i) , 0 , 5,options);
            %[g(n+1-i),resnorm(n+1+i)]= fsolve( @(x) objfun(x,h(n+1-i),p) , g(n+2-i) ,options);
            %calculate crit point for n+i
            [g(n+1+i),resnorm(n+1+i)]= lsqnonlin( @(x) objfun(x,h(n+1+i),p) , g(n+i) , 0 , 5,options);
            %[g(n+1+i),resnorm(n+1+i)]= fsolve( @(x) objfun(x,h(n+1+i),p) , g(n+i) ,options);
            clc
            toc
            progress=i/n
        end
        
        dat.g=g;
        dat.h=h;
        dat.resnorm=resnorm;
        disp('Numerical continuation complete!')
        
        
    case 4   %Search critical points in the (wei,wee) plane; for plane wave
        n=30;
        h=linspace(0.3,0.6,2*n+1);
        g=zeros(size(h));
        resnorm=g;
        
        %Initialisation; calculate crit point for h(n+1)
        disp('Initialising...')
        
        tic
        [g(n+1),resnorm(n+1)]= lsqnonlin( @(x) objfun_pw(x,h(n+1),p) , 1 , 0 , 5,options);
        %[g(n+1),resnorm(n+1)]= fsolve( @(x) objfun(x,h(n+1),p) , 1 ,options);
        toc
        disp('Start numerical continuation...')
        tic
        for i=1:n
            
            %calculate crit point for n-i
            [g(n+1-i),resnorm(n+1-i)]= lsqnonlin( @(x) objfun_pw(x,h(n+1-i),p) , g(n+2-i) , 0 , 5,options);
            %[g(n+1-i),resnorm(n+1+i)]= fsolve( @(x) objfun(x,h(n+1-i),p) , g(n+2-i) ,options);
            %calculate crit point for n+i
            [g(n+1+i),resnorm(n+1+i)]= lsqnonlin( @(x) objfun_pw(x,h(n+1+i),p) , g(n+i) , 0 , 5,options);
            %[g(n+1+i),resnorm(n+1+i)]= fsolve( @(x) objfun(x,h(n+1+i),p) , g(n+i) ,options);
            clc
            toc
            progress=i/n
        end
        
        dat.g=g;
        dat.h=h;
        dat.resnorm=resnorm;
        disp('Numerical continuation complete!')
        
    case 2  %under construction
        %search through (h,g) plane, where h,g scales ALL exct/inhi strength respectively
        %the goal is to find a point that is: 1. is a critical point 2. same
        %firing rate as the original model.
        %the underlying assumption is that, a common scaling factor gives
        %exact phenomenological match between analysis & original network
        
        n=20;
        h=linspace(0.5,1.5,2*n+1);
        g=zeros(size(h));
        resnorm=g;
        
        %Initialisation; calculate crit point for h(n+1)
        disp('Initialising...')
        
        tic
        [g(n+1),resnorm(n+1)]= lsqnonlin( @(x) objfun(x,h(n+1),p) , 1 , 0 , 5,options);
        %[g(n+1),resnorm(n+1)]= fsolve( @(x) objfun(x,h(n+1),p) , 1 ,options);
        toc
        disp('Start numerical continuation...')
        tic
        for i=1:n
            
            %calculate crit point for n-i
            [g(n+1-i),resnorm(n+1-i)]= lsqnonlin( @(x) objfun(x,h(n+1-i),p) , g(n+2-i) , 0 , 5,options);
            %[g(n+1-i),resnorm(n+1+i)]= fsolve( @(x) objfun(x,h(n+1-i),p) , g(n+2-i) ,options);
            %calculate crit point for n+i
            [g(n+1+i),resnorm(n+1+i)]= lsqnonlin( @(x) objfun(x,h(n+1+i),p) , g(n+i) , 0 , 5,options);
            %[g(n+1+i),resnorm(n+1+i)]= fsolve( @(x) objfun(x,h(n+1+i),p) , g(n+i) ,options);
            clc
            toc
            progress=i/n
        end
        
        dat.g=g;
        dat.h=h;
        dat.resnorm=resnorm;
        disp('Numerical continuation complete!')
        
end




function l=objfun3(g,h,p)
p.Ie= h*p.Ie; %scale external firing rate to E pop
if p.separateEIpop
    %[C,L,prb]=linear_response_uniform([WEE(i) wie],[wei wii],p);
    [C,L,prb]=linear_response_uniform([p.aee p.aie],[g*p.aei p.aii],p);
    %r0(i,:)=C(end,:); %C(end,:) = fixed pt with largest firing rate
else
    [C,L,prb]=linear_response_uniform(ae(i),ai(i),p);
    %r0(i)=max(C);
end
L(prb.exitflag<=0)=-inf; %remove data that fail to converge
L(prb.resnorm>1e-12)=-inf; %remove data that are likely to be local minimum

[ll,indx]=max(real(L(:)));
%Lmax.real(i)=ll;
%Lmax.imag(i)=imag(L(indx));
%[indxi,indxj]=ind2sub(size(L),indx);
%Lmax.k(i) = indxj;
%Lmax.exitflag(i)=prb.exitflag(indx);
%Lmax.resnorm(i)=prb.resnorm(indx);

l = abs(L(indx));

function l=objfun(g,h,p)
if p.separateEIpop
    %[C,L,prb]=linear_response_uniform([WEE(i) wie],[wei wii],p);
    [C,L,prb]=linear_response_uniform([h*p.aee p.aie],[g*p.aei p.aii],p);
    %r0(i,:)=C(end,:); %C(end,:) = fixed pt with largest firing rate
else
    [C,L,prb]=linear_response_uniform(ae(i),ai(i),p);
    %r0(i)=max(C);
end
L(prb.exitflag<=0)=-inf; %remove data that fail to converge
L(prb.resnorm>1e-12)=-inf; %remove data that are likely to be local minimum

[ll,indx]=max(real(L(:)));
%Lmax.real(i)=ll;
%Lmax.imag(i)=imag(L(indx));
%[indxi,indxj]=ind2sub(size(L),indx);
%Lmax.k(i) = indxj;
%Lmax.exitflag(i)=prb.exitflag(indx);
%Lmax.resnorm(i)=prb.resnorm(indx);

l = abs(L(indx));

function l=objfun_pw(g,h,p)
% same as objfun but now work for plane wave
if p.separateEIpop
    %[C,L,prb]=linear_response_uniform([WEE(i) wie],[wei wii],p);
    [C,L,prb]=linear_response_uniform_plane_wave([h*p.aee p.aie],[g*p.aei p.aii],p);
    %r0(i,:)=C(end,:); %C(end,:) = fixed pt with largest firing rate
else
    [C,L,prb]=linear_response_uniform_plane_wave(ae(i),ai(i),p);
    %r0(i)=max(C);
end
L(prb.exitflag<=0)=-inf; %remove data that fail to converge
L(prb.resnorm>1e-12)=-inf; %remove data that are likely to be local minimum

[ll,indx]=max(real(L(:)));
%Lmax.real(i)=ll;
%Lmax.imag(i)=imag(L(indx));
%[indxi,indxj]=ind2sub(size(L),indx);
%Lmax.k(i) = indxj;
%Lmax.exitflag(i)=prb.exitflag(indx);
%Lmax.resnorm(i)=prb.resnorm(indx);

l = abs(L(indx));