function imp = impulse_response(prb,ae,ai,p)
%IMPULSE_RESPONSE: 
%INPUT: prb = stationary distribution etc; 
%       ae,ai = E/I synaptic weight (scalar);p = other parameters
%OUTPUT: impulse response

T = 50; %required time length (ms)
dt = 0.1; %required temporal resolution (ms)
w=linspace(0,1/dt*pi,floor(T/dt)+1)';
w=[w; -w(end-1:-1:2)];
dw = w(2)-w(1);

w(1)=dw*1e-2;%w(end)=w(1); %approximate susceptibility at w=0 with finite limit


[r1e,prbe] = fokkerplanck1(w,prb,1,[ae p.see],[ai p.sei],p);
[r1i,prbi] = fokkerplanck1(w,prb,0,[ae p.sie],[ai p.sii],p);

r1e(1)=real(r1e(1)); %set imag part at w=0 to 0
r1i(1)=real(r1i(1));

%remove imag part and double real part
%this way, we obtain two-sided impulse response
%this eliminates discountinuity due to imag part of susceptibility at w=pi
r1e = real(r1e)*2;
r1i = real(r1i)*2;


t = 0:dt:dt*(length(w)-1);t=t(:);
impe=ifft(r1e,'symmetric')/dt;
impi=ifft(r1i,'symmetric')/dt;


%extract one-sided impulse response
indx=1:floor(length(t)/2);
t=t(indx);
impe=impe(indx);
impi=impi(indx);

%delta spike near t=0, as spectrum doesn't decay to zero due to finite dt
impe(1)=impe(2);
impi(1)=impi(2);


%save data
imp.e=impe;
imp.i=impi;
imp.t=t;
imp.dt=dt;

%debug only
while false %test ifft against analytical result
    tau = 5; %ms
    F = 1./(1/tau+1i*w);  %analytical result for F[exp(-t)]
    F = real(F)*2;    %which is also F(w) + F(-w)
    f=ifft(F,'symmetric')/dt;
    f0=exp(-t/tau);
    figure
    plot(t,real(f),t,f0);
    
    %test laplace transform
    %indx = 1:floor(length(t)/2);
    %lambda = 10+1i*10;
    %L = sum(f(indx).*exp(-lambda*t(indx)'))*dt    
    %compare to:
    %1/(1+lambda)
    
    return
end


    
end

