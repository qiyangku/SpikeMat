function y = coupling_fun(x,d,p)
%c = coupling_fun(x,d,p)
%INPUT: x = spatial coordinates;d=range;p=parameters
%OUTPUT: y = coupling probability distribution

switch p.coupling_type
    case 'vonMises'
        y=vonMises(x,d);
    case 'expo'
        y=expo(x,d,p);
end

%make sure discrete distribution is normalised to 1
dx=x(2)-x(1);
y = y/(sum(y(:))*dx.^p.dim);


function y=expo(x,d,p)
switch p.dim
    case 1
        y = exp(-abs(x)/d)/d/2;
       
    case 2
        [Y,X]=ndgrid(x,x);
        r=sqrt(Y.*Y+X.*X);    
        y0=2*pi*d*d; %normalization factor
        y = exp(-r/d)/y0;   

end



function y=wrapped_gaussian(d,N)
x=linspace(0,1,N+1);
y=zeros(1,N+1);
for n=-10:10
    z=(x+n)/d;
    y=y + exp(-z.*z)/sqrt(2*pi)/d;
end
y=y(1:end-1)';

function y=vonMises(x,d)
%von Mises distribution
%INPUT: d = width of distribution; x= spatial coordinate
k=1/d/d;
y = exp(cos(x)*k)/2/pi/besseli(0,k);

