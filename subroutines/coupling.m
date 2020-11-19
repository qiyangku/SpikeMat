function C=coupling(k,d,N)
%INPUT: k=avg in-degree; d=range of connection; N=number of neurons
%OUTPUT: binary sparse coupling matrix

x=linspace(0,2*pi,N+1); %N=number of neurons
x=x(1:N)';
dx=x(2)-x(1);

p=vonMises(x,d);
p=k*p*dx;

C=zeros(N,N);

for i=1:N    
    C(i,:) = rand(N,1)<p;
    p = circshift(p,1);
end

C=sparse(C);