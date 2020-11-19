function [pos,H,W]=subplotmod(D,indx)
%subplotmod(D,indx)
%INPUT: D = dimension of subplots; 
%       indx = horizontal-wise index of the subplot
%OUTPUT: pos = position of subplot 
%        H,W = height,width of figure; in natural coordinates.


%NB all length are in natural units & each subplot is a unit square!
t=0.25;
l=0.25;
r=0.25;
b=0.25;
L=0.6; %extra left margin
B=0.5; %extra bottom margin


[i,j]=ind2sub(D',indx); %change subplot index to position
%NB. matlab subplot indx uses horizontal convention.

%height and width of figure in natural unit (default unit square)
pos(3)=1; 
pos(4)=1;

%>>>>>> Possible to define different layouts
%height/width of figure
H = D(1) + t +B +(D(1)-1)*b;
W = D(2) + (D(2)-1)*l +r+L;
%position of subplot depends on [i,j] and margin sizes
pos(1) = t + (i-1)*(1+b);
pos(2) = L + (j-1)*(1+l);
%<<<<<

set(gca,'position',axiscrd(pos,H,W))