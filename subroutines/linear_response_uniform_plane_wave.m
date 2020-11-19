function [r0,L,prb]=linear_response_uniform_plane_wave(ae,ai,p)
%linear_response_uniform(ae,ai,p)
% Same as linear_response_uniform, but now eigenvalues are calculated for
% plane wave perturbation.
%INPUT: ae,ai = E/I synaptic strength; p = parameters
%OUTPUT: r0 = uniform firing rate; lambda = temporal eigenvalues;
%           prb = stationary distribution
%p=parameters([]);
%ae=0.03;ai=0.3;
%ae=0.01;ai=0.026; working parameters

%STEP 1: find uniform solution
[r0,prb]=fixedpt_uniform(ae,ai,p);

% %% STEP1a: inspect stationary distribution
% subplot(2,1,1)
% plot(prb.V,prb.p,p.Vref*[1 1],[0 0.14],'--')
% title(['r0 = ' num2str(R) ' Hz'])
% ylabel('P_0(v)')
% subplot(2,1,2)
% plot(prb.V,prb.dp,p.Vref*[1 1],[-0.05 0.05],'--')
% set(gcf,'color','w')
% xlabel('v (mV)')
% ylabel('P_0''(v)')
%
% %R
% %return
%
% %%


%STEP 2: spatial Fourier modes & impulse response
%x=linspace(0,pi,ceil(p.nx/2))';%x=x(1:end-1); %column vector
x=linspace(0,pi,51)';
x=[x;-x(end-1:-1:2)];
dx=x(2)-x(1);

switch p.dim
    case 1
        
        KE=p.K;KI=p.K;
        
        ke = KE*coupling_fun(x,p.de,p);%Synaptic coupling density (integral(ke)=KE)
        ki = KI*coupling_fun(x,p.di,p);
        
        we = fft(ke)*dx; %complex Fourier coefficient * 2pi
        wi = fft(ki)*dx;
        
        %plot(abs(we))
        %we = fft(ke)*(dx/pi)*dx; %we(1)=we(1)/2;
        %wi = fft(ki)*(dx/pi)*dx; %wi(1)=wi(1)/2;
        
        prb.imp = impulse_response(prb,ae(1),ai(1),p);
        
    case 2
        kee = p.KEE*coupling_fun(x,p.dee,p); %ke = ke*ke';
        kie = p.KIE*coupling_fun(x,p.die,p);
        kei = p.KEI*coupling_fun(x,p.dei,p);
        kii = p.KII*coupling_fun(x,p.dii,p);
        
        wee = fft2(kee)*dx*dx; %F[]=int(f*exp(-ikx)*dx)
        wei = fft2(kei)*dx*dx;
        wie = fft2(kie)*dx*dx;%W/numel(kie)*pi*pi
        wii = fft2(kii)*dx*dx;
        %wee = fft2(kee)*(dx*dx/pi/pi);wee(1,:)=wee(1,:)/2;wee(:,1)=wee(:,1)/2;
        %wei = fft2(kei)*(dx*dx/pi/pi);wei(1,:)=wei(1,:)/2;wei(:,1)=wei(:,1)/2;
        %wie = fft2(kie)*(dx*dx/pi/pi);wie(1,:)=wie(1,:)/2;wie(:,1)=wie(:,1)/2;
        %wii = fft2(kii)*(dx*dx/pi/pi);wii(1,:)=wii(1,:)/2;wii(:,1)=wii(:,1)/2;
        
        prb.impe = impulse_response(prb.xe,ae(1),ai(1),p);
        prb.impi = impulse_response(prb.xi,ae(2),ai(2),p);
end

%plot(1:10,diag(wee(1:10,1:10)),1:10,diag(wei(1:10,1:10)),'.')

% imagesc(real(wee));colorbar
% return
while false %for debug only; set to true for debug
% % STEP2b: (debug only) verify by reconstructing ke
% switch p.dim
%     case 1
%         ke2= 0;%real(we(1));
%         for k=1:floor(length(we)/2)
%             %ke2= ke2 + we(k)*exp(1i*k*x);
%             ke2=ke2+real(we(k))*cos((k-1)*x);
%         end
%         plot(x,ke,'o',x,real(ke2),'.')
%     case 2
%         [Y,X]=ndgrid(x,x);
%         kee2= 0;%real(we(1));
%         for m=1:size(wee,1)
%             for n=1:size(wee,2)
%                 %kee2=kee2+real(wee(m,n))*cos((m-1)*Y).*cos((n-1)*X);
%                 kee2=kee2+wee(m,n)*exp(1i*(m-1)*Y+1i*(n-1)*X);kee2=real(kee2);
%
%             end
%         end
%         plot(x,kee(10,:),'*',x,kee2(10,:),'.')
%         %plot(x,kee(1,:)-kee2(1,:),'.')
% end
% return


% % %% STEP3a: inspect fokkerplanck1 (debug only)
% close all
disp('Plotting linear response')
%X=5*linspace(0,1,51)'; %real part of temporal eigenvalue ms^-1
%Y=2.5*linspace(-1,1,51)'; %imaginary part of temporal eigenvalue
X=1*linspace(-1,1,101)'; %real part of temporal eigenvalue ms^-1
Y=1*linspace(-1,1,101)'; %imaginary part of temporal eigenvalue


[Y,X]=ndgrid(Y,X);
lambda=X+1i*Y;
% %
% % %dV=prb.V(2)-prb.V(1);
% %
% % indx=sub2ind([51,51],26,26-6);
% %
% % if p.separateEIpop
% % figure
% % [r1e,prb1] = fokkerplanck1(2*pi*lambda(:),prb.xe,1,ae(1),ai(1),p);
% % subplot(2,2,1)
% % plot(prb.xe.V,prb.xe.p,prb.xe.V,prb.xe.p+0.1*prb1.p1(indx,:))
% % axis([-90 -50 0 0.2])
% % title('P_0(V) under \delta{R_e}')
% % [r1i,prb1] = fokkerplanck1(2*pi*lambda(:),prb.xi,0,ae(2),ai(2),p);
% % subplot(2,2,2)
% % plot(prb.xi.V,prb.xi.p,prb.xi.V,prb.xi.p+0.1*prb1.p1(indx,:))
% % axis([-90 -50 0 0.2])
% % title('P_0(V) under \delta{R_i}')
% % xlabel('v (mV)')
% % set(gcf,'color','w')
% % else
% % figure
% % [r1e,prb1] = fokkerplanck1(2*pi*lambda(:),prb,1,ae,ai,p);
% % subplot(2,2,1)
% % plot(prb.V,prb.p,prb.V,prb.p+0.1*prb1.p1(indx,:))
% % axis([-90 -50 0 0.2])
% % title('P_0(V) under \delta{R_e}')
% % [r1e,prb1] = fokkerplanck1(2*pi*lambda(:),prb,0,ae,ai,p);
% % subplot(2,2,2)
% % plot(prb.V,prb.p,prb.V,prb.p+0.1*prb1.p1(indx,:))
% % axis([-90 -50 0 0.2])
% % title('P_0(V) under \delta{R_i}')
% % xlabel('v (mV)')
% % set(gcf,'color','w')
% % end
% %
% %
% % imagesc(abs(reshape(r1e,size(X))))
% % colorbar
% % %,[0 0.01])
% %
% % %
% STEP3b: inspect evan_fun (debug only)
disp('Plotting Evan function for a particular wave number')
k=10;

if p.separateEIpop
    prb.F=[wee(k,k) wei(k,k);wie(k,k) wii(k,k)];
    E=evan_fun([real(lambda(:)) imag(lambda(:))],prb,ae,ai,p);
    %plot(lambda,real(E),'.');
    E = E(:,1)+1i*E(:,2);
    figure
    %subplot(2,2,4) %abs(E) consistently>=1?? no zero unless E<--E-1???
    imagesc(X(1,:),Y(:,1),abs(reshape(E,size(X))),[0 1]);colorbar
    set(gcf,'color','w');set(gca,'ydir','normal')
    xlabel('Re(\lambda) (kHz)')
    ylabel('Im(\lambda) (kHz)')
    title(['Evans function for k = ' num2str(k)])
    
else
    prb.F=[we(k) wi(k)];
    E=evan_fun([real(lambda(:)) imag(lambda(:))],prb,ae,ai,p);
    %plot(lambda,real(E),'.');
    E = E(:,1)+1i*E(:,2);
    
    %subplot(2,2,4)
    imagesc(X(1,:),Y(:,1),abs(reshape(E,size(X))),[0 1]);
    set(gcf,'color','w');set(gca,'ydir','normal')
    xlabel('Re(\lambda) (kHz)')
    ylabel('Im(\lambda) (kHz)')
    title(['Evans function for k = ' num2str(k)])
    %colorbar
end

break

end
% 

%STEP 3: for every wave number, find temporal frequency
%disp('Finding roots of Evan function for each spatial wave number')

%options = optimoptions('lsqnonlin','FunctionTolerance',1e-16,'OptimalityTolerance',1e-16);
options = optimoptions('lsqnonlin','FunctionTolerance',1e-16,'OptimalityTolerance',1e-16,'Display', 'off');
%init = [-0.01 0.01];
ninit = 16;
L=zeros(ninit,10);%try first 10 wave numbers (for k>10, spectral power is small)
EXITFLAG=L;
RESNORM=L;
init = [linspace(-0.5,0.5,ninit) ; 0.5*ones(1,ninit)]';
%options = optimset('Display', 'off');
%tic
if p.separateEIpop
    for k=1:size(L,2) 
        %Evan function
        % Plane wave case (wave number in y direction = 0)
        prb.F=[wee(1,k) wei(1,k); wie(1,k) wii(1,k)]; % 2 dimensional
        
        for l=1:ninit
            %[ll,FVAL,EXITFLAG(l,k)]= fsolve( @(x) evan_fun(x,prb,ae,ai,p) , init(l,:) ,options);
            [ll,RESNORM(l,k),RESIDUAL,EXITFLAG(l,k)]= lsqnonlin( @(x) evan_fun(x,prb,p) , init(l,:),[-1 -1],[1 1],options);
            L(l,k)=ll(1)+1i*ll(2);
        end
        %clc
        %toc
        %k
        
    end
else
    for k=1:size(L,2) %try first 10 wave numbers
        %Evan function
        prb.F=[we(k) wi(k);];
        for l=1:ninit
        [ll,RESNORM(l,k),RESIDUAL,EXITFLAG(l,k)]= lsqnonlin( @(x) evan_fun(x,prb,p) , init(l,:),[-1 -1],[1 1],options);
        %[ll,FVAL,EXITFLAG(k)]= fsolve( @(x) evan_fun(x,prb,ae,ai,p) , init ,options);
        L(l,k)=ll(1)+1i*ll(2);
        end
    end
end
%subplot(2,2,3)
%plot(0:length(L)-1,real(L)*1e3,'o-')
%L(EXITFLAG~=1)=NaN;
prb.exitflag=EXITFLAG;
prb.resnorm=RESNORM;
% plot(real(L)*1e3,imag(L)*1e3,'o',[0 0],[-1 1],'--')
%
% %set(gcf,'color','w')
% xlabel('k')
% ylabel('real(\lambda) (Hz)')
return