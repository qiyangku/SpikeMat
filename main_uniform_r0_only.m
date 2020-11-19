clear
close all
clc

m=101;
n=101;
p=parameters('yifan');

h = linspace(0.2,1.2,n); 
g = linspace(0.5,1.5,m); %in this window, h>0.72 is zero 

r0e=zeros(m,n,3);
r0i=r0e;

tic

%don't need to compute h>0.72 as we know they are unstable ~0 Hz
indx = find(h>0.72);
indx = indx(1);


for i=1:indx
    
    for j=1:n
        
        
        
        [C,prb]=fixedpt_uniform([p.aee*h(i) p.aie], [p.aei*g(j) p.aii],p);
        [ii,jj]=size(C); %ii=number of fixed poitns; jj =2
        ii = min(ii,3);
        r0e(i,j,1:ii)=C(1:ii,1);
        r0i(i,j,1:ii)=C(1:ii,2);
        
    end
    clc
    toc
    progress=i/m
end

save data_diffapprox/uniform_r0only31c.mat
return
%%
close all
for i=1:m
    plot(g,r0i(:,i,1),'.',g,r0i(:,i,1),'.')
    plot(g,r0i(:,i,2),'.',g,r0i(:,i,2),'.')
    hold on
    plot(g,r0e(:,i,1),'*',g,r0e(:,i,1),'*')
    plot(g,r0e(:,i,2),'*',g,r0e(:,i,2),'*')
    title(num2str(h(i)))
    hold off
    ylim([0 100e-3])
    pause
end