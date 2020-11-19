function mybeep(~)
fs = 6000; %sampling rate
T = 2/pi; % 2 seconds duration
t = 0:(1/fs):T;
f = 440; %sound frequency
a = 1.1;  %amplitude
y = a*sin(2*pi*f*t);
sound(y, fs);
