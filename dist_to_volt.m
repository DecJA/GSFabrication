function [send_voltage] = dist_to_volt(d_dist)

%Converts requested chnage in distance to DAQ acquisition card voltage
%change. Moves in relative distance only, not absolute. Voltage offset does
%not change with manual DC offset on controller

offset = 0;
%scale_factor = 7.0904e-06;
% scale_factor = 9.753e-06;
% 
% d_voltage = d_dist/scale_factor;
% send_voltage = d_voltage;

%%300um travel stage
send_voltage = d_dist/3e-5 + offset;

end