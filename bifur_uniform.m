%Fine the fixed point. of uniform stationary activity
%Find Fourier transform of weight function (analytical results)?
%For every wave number k, find lambda (complex).

clear
close all

m=51;
        p=parameters([],'model_type','conduct');
        ae = linspace(0.01,0.04,m);
        ai=linspace(0.02,0.3,m);
        %[ae,ai]=ndgrid(ae,ai);
        %ae=ae(:);ai=ai(:);
        
R=zeros(m,m);
%ae=0.01;ai=0.026; working parameters

%ae=0.015;ai=0.3;
for i=1:m%length(ae)
for j=1:m
    
    [R(i,j),prb]=fixedpt_uniform(ae(i),ai(j),p);
end
end


%plot(ae,R)


mybeep;