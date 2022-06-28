function [scaleFactor,ditherPattern] = modulateAmplitude(power,freq,sz)

if power== 0
  ditherPattern = zeros(sz);
  scaleFactor = 0;
  return
end
sinFunc2Inv = @(a,b,c,d,y) 1*((asin(sqrt(y)/sqrt(a))+c) +d)/b;

%Experimental Values
% a = 1.0023; 
% b=1.5515;   
% c=-0.0198;    
% d=0.0025;
%values 15/11
% a = 1.0258;
% b = 1.4193;
% c = 0.0434;
% d = -0.0042;

%values 01-02-2022
a = 0.9946;
b = 1.4491;
c = 0.0434;
d = -0.0088;

scaleFactor = sinFunc2Inv(a,b,c,d,power);

xres = linspace(-1,1, sz(1));
xrescos=cos(xres.*(freq));
yrescos=cos((xres.').*freq);
ditherPattern = scaleFactor*(pi/2).*sign(xrescos.*yrescos);

end