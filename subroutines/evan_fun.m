function E = evan_fun(lambda,prb,p)
%E = evan_fun(f,k)
%INPUT: lambda = temporal exponent; lambda = [real_part imag_part]
%F=spatial Fourier coefficient
%OUTPUT: [real_part imag_part]

lambda = lambda(:,1) + 1i*lambda(:,2); %convert inputs to complex number
lambda = lambda(:)';

if p.separateEIpop

    
impe = prb.impe;%impulse_response(prb.xe,ae(1),ai(1),p);
impi = prb.impi;%impulse_response(prb.xi,ae(2),ai(2),p);

r1ee = sum(impe.e.*exp(-lambda.*impe.t)).*impe.dt;
r1ei = sum(impe.i.*exp(-lambda.*impe.t)).*impe.dt;
r1ie = sum(impi.e.*exp(-lambda.*impi.t)).*impi.dt;
r1ii = sum(impi.i.*exp(-lambda.*impi.t)).*impi.dt;

%[r1ee,prb1] = fokkerplanck1(lambda,prb.xe,1,ae(1),ai(1),p); 
%[r1ei,prb1] = fokkerplanck1(lambda,prb.xe,0,ae(1),ai(1),p); 

%[r1ie,prb1] = fokkerplanck1(lambda,prb.xi,1,ae(2),ai(2),p); 
%[r1ii,prb1] = fokkerplanck1(lambda,prb.xi,0,ae(2),ai(2),p); 


%M = [r1ee*prb.F(1,1) r1ei*prb.F(1,2); r1ie*prb.F(2,1) r1ii*prb.F(2,2)];
%z = det(M-eye(2));
z = (r1ee*prb.F(1,1)-1).*(r1ii*prb.F(2,2)-1) - r1ei.*r1ie*prb.F(2,1)*prb.F(1,2);

else
%[r1e,prb1] = fokkerplanck1(lambda,prb,1,ae,ai,p);
%[r1i,prb1] = fokkerplanck1(lambda,prb,0,ae,ai,p);
    
%imp = impulse_response(prb,ae(1),ai(1),p);
imp = prb.imp;
r1e = sum(imp.e.*exp(-lambda*imp.t))*imp.dt;
r1i = sum(imp.i.*exp(-lambda*imp.t))*imp.dt;
z = r1e*prb.F(1) + r1i*prb.F(2) - 1;
end



E=[real(z)' imag(z)'];
