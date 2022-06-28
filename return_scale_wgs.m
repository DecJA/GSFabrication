function [scale] = return_scale_wgs(dist)

% a = 2.136e-7;
% a = 2.075e-7; % 09/05/22
a = 2.914e-7; %20/06/22 
%Assume no offset
b = 0*4.548e-5; %not sure if this is needed yet
scale = (dist-b)./a;

end