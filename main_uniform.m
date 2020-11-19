%Find stationary uniform solution and its stability based on linear response
%Generate bifurcation diagram
%if ~exist('i','var') %if simulation hasn't started
clear
close all
clc


m=51; %number of parameter values

dataset=27;
switch dataset
    case 1 %current-based, Turing bifurcation
        p=parameters([],'Vlb',-150,'nV',2001);
        ae = linspace(0.005,0.02,m);
        ai = 0.026*ones(size(ae));
    case 2 %conductance-based;
        p=parameters('qi1d');
        ae = linspace(0.02,0.024,m);
        ai = 0.2*ones(size(ae));
    case 3 %current-based, uniform solution is Turing stable. Exist bump - not calculated here
        p=parameters([],'Vlb',-150,'nV',2001);
        ae = linspace(0.02,0.04,m);
        ai = 0.3*ones(size(ae));
    case 4 %combine 1&3 and sweep through both ae and ai
        p=parameters([],'Vlb',-150,'nV',2001);
        ae=linspace(0.005,0.04,m);
        ai=linspace(0.02,0.3,m);
        [ae,ai]=ndgrid(ae,ai);
        ae=ae(:);ai=ai(:);
    case 5 %conductance-based; 2d survey
        p=parameters('qi1d');
        ae = linspace(0.01,0.04,m);
        ai=linspace(0,0.3,16);ai=ai(2:end);
        [ae,ai]=ndgrid(ae,ai);
        ae=ae(:);ai=ai(:);
        m=length(ae);
        
    case 6 %same as 2;
        p=parameters([],'model_type','conduct');
        ae = linspace(0.01,0.025,m);
        ai = 0.03*ones(size(ae));
    case 11 %separate EI pop %under construction
        
        p=parameters([],'model_type','conduct','separateEIpop',true);
        WE = linspace(0.01,0.06,m);
        wie = p.aie;
        wei = 0.1;
        wii = 0.2;
    case 12 %no separate EI pop; 2dim network %under construction
        p=parameters([],'model_type','conduct');
        ae = linspace(0.01,0.04,m);
        ai=linspace(0.02,0.3,m);
        [ae,ai]=ndgrid(ae,ai);
        ae=ae(:);ai=ai(:);
    case 21 %for yifan's paper
        p=parameters('yifan');
        wee = 0.4775*p.aee;%0.042;%linspace(0.01,0.06,m);
        wie = p.aie;%0.05;
        %wei = 0.1;
        g = linspace(0.5,1.5,m);
        WEI = p.aei*g;
        wii = p.aii;
        
        
    case 22
        p=parameters([],'model_type','conduct','separateEIpop',true);
        WEE = linspace(0.01,0.04,m);
        wie = p.aie;%0.05;
        wei = 0.075;
        wii = p.aii;

    case 23
        p=parameters('yifan');
        h = 0.1:0.05:1.1;
        g = linspace(0.5,1.5,m);
        
        WEE = p.aee*h;%0.042;%linspace(0.01,0.06,m);
        wie = p.aie;%0.05;
        WEI = p.aei*g;
        wii = p.aii;
        
        [WEE,WEI]=ndgrid(WEE,WEI);
        WEE=WEE(:);WEI=WEI(:);
        m=length(WEE);
    case 24
        p=parameters('yifan');
        
        g = linspace(0.5,1.5,m);
        h = g;
        
        [H,G]=ndgrid(h,g);
        H=H(:);G=G(:);
        
        WEE = p.aee*H;%0.042;%linspace(0.01,0.06,m);
        WIE = p.aie*H;%0.05;
        WEI = p.aei*G;
        WII = p.aii*G;       
        
        jee0=p.jee;
        jie0=p.jie;

        %[WEE,WEI]=ndgrid(WEE,WEI);
        %WEE=WEE(:);WEI=WEI(:);
        m=length(WEE);
        case 25
        p=parameters('yifan');
        
        H = linspace(0.1,0.51,m);
        %h = g;
        
        %[H,G]=ndgrid(h,g);
        %H=H(:);G=G(:);
        
        WEE = p.aee*H;%0.042;%linspace(0.01,0.06,m);
        WIE = p.aie*H;%0.05;
        WEI = p.aei*H;
        WII = p.aii*H;       
        
        jee0=p.jee;
        jie0=p.jie;

        %[WEE,WEI]=ndgrid(WEE,WEI);
        %WEE=WEE(:);WEI=WEI(:);
        %m=length(WEE);
        
    case 26 %repeat 21, for Yifan's paper
        
        p=parameters('yifan');
        
        g = linspace(0.5,1.5,m);
        h = 0.4775;
        
        [H,G]=ndgrid(h,g);
        H=H(:);G=G(:);
        
        WEE = p.aee*H;%0.042;%linspace(0.01,0.06,m);
        WIE = p.aie*ones(size(H));%0.05;
        WEI = p.aei*G;
        WII = p.aii*ones(size(G));
        
        %jee0=p.jee;
        %jie0=p.jie;
        
        %[WEE,WEI]=ndgrid(WEE,WEI);
        %WEE=WEE(:);WEI=WEI(:);
        m=length(WEE);
    case 27 %repeat 26, for Yifan's paper, but with h=1
        
        p=parameters('yifan');
        
        g = linspace(0.5,1.5,m);
        h = 1;
        
        [H,G]=ndgrid(h,g);
        H=H(:);G=G(:);
        
        WEE = p.aee*H;%0.042;%linspace(0.01,0.06,m);
        WIE = p.aie*ones(size(H));%0.05;
        WEI = p.aei*G;
        WII = p.aii*ones(size(G));
        
        %jee0=p.jee;
        %jie0=p.jie;
        
        %[WEE,WEI]=ndgrid(WEE,WEI);
        %WEE=WEE(:);WEI=WEI(:);
        m=length(WEE);
        
        
        
end


r0=zeros(m,p.dim);
test_conv=r0; %test convergence at Vlb
Lmax.real = zeros(m,2); %temporal eigenvalue with largest real part
Lmax.imag = zeros(m,2); %imaginary part of max L
Lmax.k = zeros(m,2);    %spatial wave number corresponding to maxl
Lmax.exitflag = zeros(m,2); %convergence of the root finding

L2 = cell(m,2);

P0=zeros(p.nV,m);  %voltage distribution
Re0=r0;     %excitatory rate (kHz)
Ri0=r0;     %inhibitory rate (kHz)

%init=1;
tic
%else %if simulation is interrupted by error and need to resume
%    init=i+1;
%end


for i=1:m

if p.separateEIpop
%[C,L,prb]=linear_response_uniform([WEE(i) wie],[wei wii],p);
%[C,L,prb]=linear_response_uniform([WEE(i) wie],[WEI(i) wii],p);
%[C,L,prb]=linear_response_uniform([wee wie],[WEI(i) wii],p);

%p.jee=jee0*H(i);p.jie=jie0*H(i); for 25

[C,L,prb]=linear_response_uniform([WEE(i) WIE(i)],[WEI(i) WII(i)],p);
[C_pw,L_pw,prb_pw]=linear_response_uniform_plane_wave([WEE(i) WIE(i)],[WEI(i) WII(i)],p);

r0(i,:)=C(end,:); %C(end,:) = fixed pt with largest firing rate

L2{i,1} = L; %bump
L2{i,2} = L_pw; %plane wave

else
[C,L,prb]=linear_response_uniform(ae(i),ai(i),p);
r0(i)=max(C);
end
%test_conv(i)=prb.p(1)/max(prb.p);

%P0(:,i) = prb.p;
%Re0(i) = prb.Re;
%Ri0(i) = prb.Ri;

%indx = find(real(L)==max(real(L)));

%if isempty(indx)
%    k(i)=NaN;
%    Lmax(i)=NaN;
%else
%k(i)=indx-1;
%Lmax(i) = L(indx);
%end

%>>>>>>>>>>> bump instability
L(prb.exitflag<=0)=-inf; %remove data that fail to converge
L(prb.resnorm>1e-12)=-inf; %remove data that are likely to be local minimum
%L(:,1)=-inf; %remove data for 0-th fourier mode (uniform perturbation)

[ll,indx]=max(real(L(:)));
Lmax.real(i,1)=ll;
Lmax.imag(i,1)=imag(L(indx));
[indxi,indxj]=ind2sub(size(L),indx);
Lmax.k(i,1) = indxj;
Lmax.exitflag(i,1)=prb.exitflag(indx);
Lmax.resnorm(i,1)=prb.resnorm(indx);
%>>>>>>>>>>> plane wave instability
L_pw(prb_pw.exitflag<=0)=-inf; %remove data that fail to converge
L_pw(prb_pw.resnorm>1e-12)=-inf; %remove data that are likely to be local minimum

[ll,indx]=max(real(L_pw(:)));
Lmax.real(i,2)=ll;
Lmax.imag(i,2)=imag(L_pw(indx));
[indxi,indxj]=ind2sub(size(L_pw),indx);
Lmax.k(i,2) = indxj;
Lmax.exitflag(i,2)=prb_pw.exitflag(indx);
Lmax.resnorm(i,2)=prb_pw.resnorm(indx);
%<<<<<<<<<<


clc
toc
progress = i/m
end

return

if dataset==4 || dataset==5
    r0 = reshape(r0,[m m]);
    test_conv = reshape(test_conv,[m m]);
    k = reshape(k,[m m]);
    Lmax = reshape(Lmax,[m m]);
    
    ae=reshape(ae,[m m]);ai=reshape(ai,[m m]);
    ae = ae(:,1);ai=ai(1,:);
    
end
clear i
return

%% plot image
close all
%subplot(1,2,1)
imagesc(ai,ae,r0*1e3,[0 10]);colorbar
%imagesc(ai,ae,log(r0),[-6 -4]);
xlabel('a_i');ylabel('a_e')
set(gcf,'color','w')
%subplot(1,2,2)
figure
imagesc(ai,ae,real(Lmax))
xlabel('a_i');ylabel('a_e')
title('Turing instability')
colorbar
colormap(b2r( -0.02, 0.02 ) )
set(gcf,'color','w')
%% plot curve
close all
subplot(2,1,1)
plot(ae,r0*1e3,'.')
ylabel('r_0 (Hz)')
ylim([0 20])
subplot(2,1,2)
plot(ae,real(Lmax))
ylim([-0.1 0.1])
ylabel('Max eigenvalue')
xlabel('a_e')