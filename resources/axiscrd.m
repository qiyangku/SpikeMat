function pos=axiscrd(pos,H,W)
%pos=axiscrd(pos,H,W)
%Convert coordinate of subplot from 'natural' coordinate to MATLAB default
%realtive coordinates.
%INPUT: pos = [y x height width]. y x are relative to topleft corner of subplot
%       H,W defines the size of the figure.


%We split this function from custom subplot function...
%because sometimes it's desirable to manually adjust some subplo

x= pos(2)/W;
y= 1 - (pos(1)+pos(3))/H;

w=pos(4)/W;
h=pos(3)/H;

pos=[x y w h]; %[x y width height] relative to figure







